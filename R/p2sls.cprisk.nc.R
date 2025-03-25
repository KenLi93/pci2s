#' Counterfactual marginal CIF for proximal two-stage-least-squares with survival outcome using competing risks as negative control outcomes, but the exposure has an effect to the competing risk during the initial period after the treatment
#'
#' This function computes the counterfactual marginal cumulative cause-specific hazard functions and cumulative incidence functions for the proximal 
#' two-stage-least-squares with survival outcome, using competing risks as the negative control outcomes
#'
#' @param times an n-vector of observed time to event outcomes
#' @param cause an n-vector of indicators for causes of events; in each entry, 0 indicates primary event, 1 indicates competing risks, and -1 indicates censoring
#' @param A an n-vector of primary exposure
#' @param X an n * nX matrix of adjusted covariates , can be empty
#' @param Z an n * nZ matrix of negative control exposure (NCE)
#' @param nc_time time after which the exposure has no effect on the competing risk
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
#' p2sls_rslt <- p2sls.cprisk.nc(times = times, cause = cause, nc_time = 2, A = A, X = X, Z = Z, bootstrap = T, nboot = 100)
#' @export
#' 
p2sls.cprisk.nc <- function(times, cause, A, X1, X2=NULL, Z, nc_time, 
                            bootstrap = F, nboot = 1000, conf.level = 0.95) {
  if (!is.null(X1) && is.null(X2)) {
        X2  <- X1
  }
  
  est <- p2sls.cprisk.nc.est(times, cause, A, X1,X2, Z, nc_time)
  
  if (bootstrap == T) {
    boot_ci <- p2sls.cprisk.nc.boot(times, cause, A, X1,X2, Z, nc_time, conf.level, nboot)
    out <- c(est, boot_ci)
  } else {
    out <- est
  }

  return(out)
}



## point estimates
p2sls.cprisk.nc.est <- function(times, cause, A, X1,X2, Z, nc_time) {
  nn <- length(times)
  # clean data type
  Y <- W <- times
  D <- as.numeric(cause == 0)
  Dc <- as.numeric(cause == 1)
  
  
  ## first stage model
  time_post <- times[times > nc_time]
  cov_post <- cbind(A, X1, Z)[times > nc_time,]
  D_post <- D[times > nc_time]
  Dc_post <- Dc[times > nc_time]
  
  mod1_nc <- lin_ah(time = time_post, event = Dc_post, covariates = cov_post)
  param1_nc <- mod1_nc$ESTIMATE 
  mu1_nc <- as.numeric(as.matrix(cbind(A, X1, Z)) %*% param1_nc)
  
  ## second stage model
  mod2s <- lin_ah(time = Y, event = D, covariates = cbind(A, X2, mu1_nc))
  
  beta_a <- mod2s$ESTIMATE[1]  
  
  param_2s <- mod2s$ESTIMATE
  
  return(list(beta_a = beta_a, param1 = param1_nc, param2 = param_2s))
}


## bootstrap inference
p2sls.cprisk.nc.boot <- function(times, cause, A, X1, X2, Z, nc_time, conf.level = 0.95, nboot = 1000) {
  
  nn <- length(times)
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  Z <- as.matrix(Z)
  # clean data type
  
  boot_result <- vector("list", length = nboot)
  
  for (bb in 1:nboot) {
    boot_id <- sample(nn, replace = T)
    times_boot <- times[boot_id]
    cause_boot <- cause[boot_id]
    A_boot <- A[boot_id]
    X1_boot <- X1[boot_id,, drop = F]
    X2_boot <- X2[boot_id,, drop = F]
    Z_boot <- Z[boot_id,, drop = F] 
    
    boot_result[[bb]] <- p2sls.cprisk.nc.est(times_boot, cause_boot, A_boot, X1_boot, X2_boot, Z_boot, nc_time)
  }
  alp <- 1 - conf.level
  all_beta_a <- sapply(boot_result, function(x) x$beta_a)
  beta_a_ci <- quantile(all_beta_a, c(alp / 2, 1 - alp / 2))
  
  all_param1 <- sapply(boot_result, function(x) x$param1)
  param1_ci <- apply(all_param1, 1, function(x) quantile(x, c(alp / 2, 1 - alp / 2)))
  
  all_param2 <- sapply(boot_result, function(x) x$param2)
  param2_ci <- apply(all_param2, 1, function(x) quantile(x, c(alp / 2, 1 - alp / 2)))
  
  return(list(beta_a_ci = beta_a_ci, param1_ci = param1_ci, param2_ci = param2_ci))
}



#' Counterfactual marginal CIF for proximal two-stage-least-squares with survival outcome using competing risks as negative control outcomes, but the exposure has an effect to the competing risk during the initial period after the treatment
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
#' @param nc_time time after which the exposure has no effect on the competing risk
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
#' p2sls_rslt_a1 <- p2sls.cprisk.nc.cif(times = times, cause = cause, nc_time = 2, A = A, a = 1, X = X, Z = Z)
#' p2sls_rslt_a0 <- p2sls.cprisk.nc.cif(times = times, cause = cause, nc_time = 2, A = A, a = 0, X = X, Z = Z)
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
p2sls.cprisk.nc.cif <- function(times, cause, A, a, X, Z, nc_time, 
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
  time_post <- times[times > nc_time]
  cov_post <- cbind(A, X, Z)[times > nc_time,]
  D_post <- D[times > nc_time]
  Dc_post <- Dc[times > nc_time]
  
  mod1_nc <- lin_ah(time = time_post, event = Dc_post, covariates = cov_post)
  param1_nc <- mod1_nc$ESTIMATE 
  mu1_nc <- as.numeric(as.matrix(cbind(A, X, Z)) %*% param1_nc)
  
  ## second stage model
  mod2s <- lin_ah(time = Y, event = D, covariates = cbind(A, X, mu1_nc))
  
  beta_a <- mod2s$ESTIMATE[1]  
  
  ## marginal ah model for the primary event
  mod0 <- lin_ah(time = Y, event = D, covariates = cbind(A, X, Z))
  param0 <- mod0$ESTIMATE
  mu0 <- as.numeric(as.matrix(cbind(A, X, Z)) %*% param0)  ## Linear predictors in the first stage model
  
  cumhaz0_func <- stepfun(x = unique(sort(times)),
                          y = c(0, mod0$CUMHAZ_K))
  cumhaz0_tseq <- cumhaz0_func(tseq) 
  
  
  ## obtain the estimated conditional survival functions
  mod_axz_pre <- lin_ah(time = pmin(times, nc_time),
                        event = as.numeric(cause >= 0 & times <= nc_time),
                        covariates = cbind(A, X, Z))
  
  haz_axz_pre <- data.frame(time = unique(sort(pmin(times, nc_time))),
                            haz = mod_axz_pre$HAZ)
  
  mu_axz_pre <- c(as.matrix(cbind(A, X, Z)) %*% mod_axz_pre$ESTIMATE)
  
  mod_axz_post <- lin_ah(time = time_post,
                         event = D_post,
                         covariates = cov_post)
  
  haz_axz_post <- data.frame(time = unique(sort(time_post)),
                             haz = mod_axz_post$HAZ)
  
  mu_axz_post <- c(as.matrix(cbind(A, X, Z)) %*% mod_axz_post$ESTIMATE)
  
  haz_mat <- rbind(haz_axz_pre, haz_axz_post)
  haz_mat$cumhaz <- cumsum(haz_mat$haz)
  
  cumhaz_func <- stepfun(haz_mat$time,
                         c(0, haz_mat$cumhaz))
  
  cumhaz_tseq <- cumhaz_func(tseq)
  
  ## calculate the counterfactual marginal survival functions separately
  tint <- tmax / (nt - 1) ## length of the interval
  
  surva_tseq <- cumhaz0_a_tseq <- 
    haz0_a_num_tseq <- cif0_a_tseq <- rep(NA, nt)
  
  for (tt in seq_along(tseq)) {
    t1 <- tseq[tt]
    if (t1 <= nc_time) {
      surv_tt <- exp(-cumhaz_tseq[tt] - mu_axz_pre * t1)
    } else {
      surv_tt <- exp(-cumhaz_tseq[tt] - mu_axz_pre * nc_time - mu_axz_post * (t1 - nc_time))
    }
    surva_tseq[tt] <- mean(exp(-beta_a * a * t1 + beta_a * A * t1) * surv_tt)
    
    if (surv_correct == T & tt > 1) {
      surva_tseq[tt] <- min(surva_tseq[tt - 1], surva_tseq[tt])
    }
    
    haz0_a_num_tseq[tt] <- mean((mu0 + beta_a * a - beta_a * A) * exp(-beta_a * a * t1 + beta_a * A * t1) * surv_tt)
  }
  
  
  haz0_a_tseq <- haz0_a_num_tseq / surva_tseq
  
  ## obtain the counterfactual marginal cumulative cause-specific hazard functions
  for (tt in seq_along(tseq)) {
    if (tt == 1) {
      cumhaz0_a_tseq[tt] <- 0; 
    } else {
      cumhaz0_a_tseq[tt] <- cumhaz0_a_tseq[tt - 1] + tint * haz0_a_tseq[tt]
    }
  }
  
  cumhaz0_a <- cumhaz0_a_tseq + cumhaz0_tseq
  
  
  ## obtain the counterfactual marginal cumulative incidence functions
  
  for (tt in seq_along(tseq)) {
    if (tt == 1) {
      cif0_a_tseq[tt] <- 0; 
    } else {
      cif0_a_tseq[tt] <- cif0_a_tseq[tt - 1] + surva_tseq[tt] * (cumhaz0_tseq[tt] - cumhaz0_tseq[tt - 1]) + tint * haz0_a_num_tseq[tt]
    }
  }
  
  
  out <- data.frame(t = tseq, 
                    survfunc = surva_tseq,
                    cumhaz0 = cumhaz0_a,
                    cif0 = cif0_a_tseq)
  
  return(out)
}
