#' Proximal two-stage-least-squares with survival outcome using competing risks as negative control outcomes
#'
#' This function computes the counterfactual marginal cumulative cause-specific hazard functions and cumulative incidence functions for the proximal 
#' two-stage-least-squares with survival outcome, using competing risks as the negative control outcomes
#'
#' @param times an n-vector of observed time to event outcomes
#' @param cause an n-vector of indicators for causes of events; in each entry, 0 indicates primary event, 1 indicates competing risks, and -1 indicates censoring
#' @param A an n-vector of primary exposure
#' @param a a unidimensional scalar value to fix A=a for calculating the counterfactual survival curve
#' @param X an n * nX matrix of adjusted covariates , can be empty
#' @param Z an n * nZ matrix of negative control exposure (NCE)
#' @param surv_correct logical, whether the survival function will be corrected post hoc to be non-increasing and
#' fall between 0 and 1, default to be TRUE.
#' @param tmax The maximal time to calculate the counterfactual marginal survival functions. Default to be the 
#' maximal event time multiplied by 1.05;
#' @param nt Number of time points to evaluate the survival functions, evenly distributed between 0 and tmax. 
#' @returns A data frame with two columns: t: time points where the counterfactual marginal survival functions are evaluated; 
#' survfunc: value of the counterfactual marginal survival function
#' @examples
#' set.seed(2050)
#' N <- 2000
#' U <- runif(N); X <- runif(N)
#' expit  <-  function(x) exp(x)/(1 + exp(x))
#' A <- rbinom(N, 1, expit(-3 + 5 * U + 1 * X))
#' Y1 <- rexp(N, 0.5 * U + 0.2 * X + 0.2 * A)  ## event times for primary event of interest
#' Y2 <- rexp(N, 0.3 * U + 0.6 * X)  ## event times for the competing risk
#' C <- runif(N, 3, 5)  ## censoring times
#' times <- pmin(Y1, Y2, C)
#' cause <- ifelse(Y1 <= C, ifelse(Y1 <= Y2, 0, 1), ifelse(Y2 <= C, 1, -1))
#' Z <- rnorm(N, 2 * U + 0.5 * X)
#' ## Obtain the counterfactual marginal survival curves under a = 0, 1
#'
#' p2sls_rslt_a1 <- p2sls.cprisk(times = times, cause = cause, A = A, a = 1, X = X, Z = Z)
#' p2sls_rslt_a0 <- p2sls.cprisk(times = times, cause = cause, A = A, a = 0, X = X, Z = Z)
#' 
#' 
#' ## Plot the counterfactual marginal cause-specific cumulative hazard function for the primary event 
#' plot(cumhaz0 ~ t, data = p2sls_rslt_a1, type = "l", lwd = 2,
#'     xlab = "Time", ylab = "Cumulative hazard function", col = "red")
#' lines(cumhaz0 ~ t, data = p2sls_rslt_a0, lwd = 2) 
#' legend("topleft", legend = c("a = 0", "a = 1"), lwd = 2,
#'       col = c("black", "red"))
#'       
#' ## Plot the counterfactual marginal cumulative incidence functions for the primary event 
#' plot(cif0 ~ t, data = p2sls_rslt_a1, type = "l", lwd = 2,
#'     xlab = "Time", ylab = "Cumulative incidence", col = "red")
#' lines(cif0 ~ t, data = p2sls_rslt_a0, lwd = 2) 
#' legend("topleft", legend = c("a = 0", "a = 1"), lwd = 2,
#'       col = c("black", "red"))     
#'       
#' @export
#' 
p2sls.cprisk <- function(times, cause, A, a, X, Z, 
                         surv_correct = T,
                         tmax = NULL, nt = 1000) {
  
  if (is.null(tmax)) {
    tmax <- max(times) * 1.05
  }
  
  tseq <- seq(0, tmax, length.out = nt)
  
  nn <- length(times)
  # clean data type
  Y <- W <- times
  D <- as.numeric(cause == 0)
  Dc <- as.numeric(cause == 1)
  
  
  ## first stage model
  mod1 <- lin_ah(time = W, event = Dc, covariates = cbind(A, X, Z))
  param1 <- mod1$ESTIMATE 
  mu1 <- as.numeric(cbind(A, X, Z) %*% mod1$ESTIMATE)

  
  cumhaz1_func <- stepfun(x = unique(sort(times)),
                          y = c(0, mod1$CUMHAZ_K))
  cumhaz1_tseq <- cumhaz1_func(tseq)
  
  ## second stage model
  mod2s <- lin_ah(time = Y, event = D, covariates = cbind(A, X, mu1))
  
  beta_a <- mod2s$ESTIMATE[1]  
  
  ## marginal ah model for the primary event
  mod0 <- lin_ah(time = Y, event = D, covariates = cbind(A, X, Z))
  param0 <- mod0$ESTIMATE
  mu0 <- as.numeric(cbind(A, X, Z) %*% param0)  ## Linear predictors in the first stage model
 
  cumhaz0_func <- stepfun(x = unique(sort(times)),
                          y = c(0, mod0$CUMHAZ_K))
  cumhaz0_tseq <- cumhaz0_func(tseq) 
  
 
  ## obtain the estimated conditional survival functions
  mod_axz <- lin_ah(time = times,
                    event = as.numeric(cause >= 0),
                    covariates = cbind(A, X, Z))
  
  cumhaz_func <- stepfun(x = unique(sort(times)),
                         y = c(0, mod_axz$CUMHAZ_K))
  cumhaz_tseq <- cumhaz_func(tseq)
  
  mu_axz <- c(cbind(A, X, Z) %*% mod_axz$ESTIMATE)
  
  ## calculate the counterfactual marginal survival functions separately
  tint <- tmax / (nt - 1) ## length of the interval
  
  surva_tseq <- cumhaz1_a_tseq <- cumhaz0_a_tseq <- 
    haz1_a_num_tseq <- haz0_a_num_tseq <-
    cif1_a_tseq <- cif0_a_tseq <- rep(NA, nt)
  
  for (tt in seq_along(tseq)) {
    surv_tt <- exp(-cumhaz_tseq[tt] - mu_axz * tseq[tt])
    surva_tseq[tt] <- mean(exp(-beta_a * a * tseq[tt] + beta_a * A * tseq[tt]) * surv_tt)
    
    if (surv_correct == T & tt > 1) {
      surva_tseq[tt] <- min(surva_tseq[tt - 1], surva_tseq[tt])
    }
    
    haz1_a_num_tseq[tt] <- mean(mu1 * exp(-beta_a * a * tseq[tt] + beta_a * A * tseq[tt]) * surv_tt)
    haz0_a_num_tseq[tt] <- mean(mu0 * exp(-beta_a * a * tseq[tt] + beta_a * A * tseq[tt]) * surv_tt)
  }
  
  haz1_a_tseq <- haz1_a_num_tseq / surva_tseq
  haz0_a_tseq <- haz0_a_num_tseq / surva_tseq
  
  ## obtain the counterfactual marginal cumulative cause-specific hazard functions
  for (tt in seq_along(tseq)) {
    if (tt == 1) {
      cumhaz0_a_tseq[tt] <- 0; cumhaz1_a_tseq[tt] <- 0
    } else {
      cumhaz0_a_tseq[tt] <- cumhaz0_a_tseq[tt - 1] + tint * haz0_a_tseq[tt]
      cumhaz1_a_tseq[tt] <- cumhaz0_a_tseq[tt - 1] + tint * haz1_a_tseq[tt]
    }
  }
  
  cumhaz0_a <- cumhaz0_a_tseq + cumhaz0_tseq
  cumhaz1_a <- cumhaz1_a_tseq + cumhaz1_tseq

  ## obtain the counterfactual marginal cumulative incidence functions
  
  for (tt in seq_along(tseq)) {
    if (tt == 1) {
      cif0_a_tseq[tt] <- 0; cif1_a_tseq[tt] <- 0
    } else {
      cif0_a_tseq[tt] <- cif0_a_tseq[tt - 1] + surva_tseq[tt] * (cumhaz0_tseq[tt] - cumhaz0_tseq[tt - 1]) + tint * haz0_a_num_tseq[tt]
      cif1_a_tseq[tt] <- cif1_a_tseq[tt - 1] + surva_tseq[tt] * (cumhaz1_tseq[tt] - cumhaz1_tseq[tt - 1]) + tint * haz0_a_num_tseq[tt]
    }
  }
  

  out <- data.frame(t = tseq, 
                    survfunc = surva_tseq,
                    cumhaz0 = cumhaz0_a,
                    cumhaz1 = cumhaz1_a,
                    cif0 = cif0_a_tseq, 
                    cif1 = cif1_a_tseq)
  
  return(out)
}
