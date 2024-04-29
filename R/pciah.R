#' Proximal causal inference with survival outcome using two-stage-least-square
#'
#' 'pciah2s()' computes the adjusted hazard difference using additive hazards model by
#' Lin \& Ying 1994, using a pair of negative control variable to control for unmeasured
#' confounding
#'
#' @param Y an n-vector of observed time to event outcomes
#' @param D an n-vector of event indicators for Y
#' @param A an n-vector of primary exposure
#' @param X an n * nX matrix of adjusted covariates , can be empty
#' @param W an n * nW matrix of negative control outcome (NCO)
#' @param Z an n * nZ matrix of negative control exposure (NCE). Alternatively,
#' the NCE can be included in the design matrix of the first-stage model by specifying
#' Xw
#' @param Xw a matrix with n rows or a list of length nW, each element is a design matrix of
#' the corresponding NCO. Default is a matrix with columns of A, Z, X, and interactions between
#' A, Z and A, X
#' @param Xy a matrix with n rows, covariates to be adjusted for in the second-stage
#' model besides A and predictors for W. This may include interaction between A and columns of X.
#' Default to be X
#' @param nco_type An nW-vector for the types of models used in the first-stage, including
#' "linear" (linear model), "loglin" (log-linear model), "poisson" (Poisson regression model),
#' "negbin" (negative-binomial regression model), and "ah" (Lin \& Ying's additive hazards model)
#' @param nco_args A list of length nW containing additional arguments for the first-stage models.
#' Each element in the list is a sublist and should contain a vector named "offset" (default 0).
#' If nco_type == "ah", the sublist needs to include another vector named "event" as the event indicator
#' (default 1). If nco_type == "negbin", the sublist needs to include a vector "init" as the initial values,
#' which can be NA.
#' @param se_method Method to compute the standard error, can be one of "all", "analytic", "exponential multiplier",
#' "gaussian multiplier" or "none".
#' @returns A list with three elements: ESTIMATE includes the parameter estimates from the second-stage
#' model; SE produces their standard error, and PARAM includes all parameter estimates
#' @examples
#' N <- 2000
#' U <- runif(N); X <- runif(N)
#' expit  <-  function(x) exp(x)/(1 + exp(x))
#' A <- rbinom(N, 1, expit(-3 + 5 * U + 1 * X))
#' Y <- rexp(N, 0.2 + 1.5 * U + 0.2 * X + 0.2 * A)
#' D <- as.numeric(Y < 5)
#' Y[D == 0] <- 5
#' Z <- rnorm(N, 2 * U + 0.5 * X)
#' W2 <- rexp(N, 0.1 + 1 * U + 0.1 * X)
#' D2 <- as.numeric(W2 < 5)
#' W2[D2 == 0] <- 5
#' W <- cbind(rnbinom(N, size = 25,
#'                    mu = exp(2.5 * U + 0.2 * X)),
#'            W2)
#' pciah2s(Y = Y, D = D, A = A, X = X,
#'        W = W, Z = Z, nboot = 100,
#'        nco_type = c("negbin", "ah"),
#'        nco_args = list(list(offset = rep(0, N)),
#'                        list(offset = rep(0, N),
#'                             event = D2)))
#' @export
pciah2s <- function(Y, D, A, X = NULL, W, Z = NULL, Xw = NULL,
                    Xy = NULL,
                    nboot = 2000,
                    nco_type = NULL,
                    nco_args = NULL,
                    se_method = "all", verbose =F) {



  nn <- length(Y)
  # clean data type
  Y0 <- as.numeric(Y)
  D0 <- as.numeric(D)
  W0 <- matrix(as.matrix(W), nrow = length(Y0))
  nW <- ncol(W0)
  if (!is.null(X)) {
    X0 <- as.matrix(X)
  } else {
    X0 <- X
  }

  A0 <- as.numeric(A)

  if (is.null(nco_type)) {
    nco_type <- rep("linear", nW)
  }

  if (length(nco_type) != nW) {
    stop("length of nco_type should equal the number of columns of W.")
  }

  #
  if (!all(nco_type %in% c("linear", "poisson", "negbin", "loglin", "ah"))) {
    stop("'nco_type' should be left empty or one of 'linear', 'poisson', 'negbin',
         'loglin', and 'ah'.")
  }


  if (is.null(Xw)) {
    if (is.null(Z)) {
      stop("One of Z and Xw needs to be specified")
    } else {
      ## Default
      Z0 <- as.matrix(Z)
      if (is.null(X0)) {
        Xw0 <- lapply(1:nW, function(i) {  ## if no covariates, design matrix is A, Z and their interaction
          if (nco_type[i] != "ah") {
            cbind(1, A0, Z0, Z0 * A0)
          } else {
            cbind(A0, Z0, Z0 * A0)
          }
        })
      } else {
        Xw0 <- lapply(1:nW, function(i) { ## with covariates, design matrix is A, Z, X and AZ, AX interactions
          if (nco_type[i] != "ah") {
            cbind(1, A0, X0, Z0, X0 * A0, Z0 * A0)
          } else {
            cbind(A0, X0, Z0, X0 * A0, Z0 * A0)
          }
        })
      }
    }
  } else {
    if (is.data.frame(Xw) | is.matrix(Xw)) {  ## same design matrix for every NCO
      Xw0 <- lapply(1:nW, function(i) as.matrix(Xw))
    } else {
      if (is.list(Xw)) {
        Xw0 <- lapply(Xw, as.matrix)
      }
    }
  }

  if (is.null(Xy)) {
    Xy0 <- matrix(as.matrix(X), nrow = nn)
  } else {
    Xy0 <- matrix(as.matrix(Xy), nrow = nn)
  }


  if (is.null(nco_args)) {
    nco_args <- lapply(1:nW,
                       function(i) {
                         if (nco_type[i] == "ah") {
                           list(offset = rep(0, nn),
                                event = rep(1, nn))
                         } else if (nco_type[i] == "negbin") {
                           list(offset = rep(0, nn),
                                init = NA)
                         } else {
                           list(offset = rep(0, nn))
                         }
                       })
  } else {
    if (!is.list(nco_args)) {
      stop("nco_args should either be left empty or a list of length ncol(W)")
    } else if (length(nco_args) != nW) {
      stop("nco_args should be a list of length ncol(W)")
    }
  }

  ## sort the data based on time to event
  Y_order <- order(Y0)
  t1 <- Y0[Y_order]
  d1 <- D0[Y_order]




  if (!is.null(X0)) X1 <- X0[Y_order, ]

  t2 <- unique(t1)
  ntime <- length(t2) # length of unique time points
  dtime <- c(t2[1], diff(t2)) # increments of unique time points

  # obtain the order of time
  o1 <- dplyr::left_join(data.frame(t1 = t1),
                  data.frame(t1 = t2, o1 = 1:length(unique(t1))),
                  by = "t1")$o1
  tmin <- sapply(1:ntime, function(k) min(which(o1 == k)))


  W1 <- matrix(W0[Y_order,], nrow = nn)
  Xw1 <- lapply(Xw0, function(xx) matrix(xx[Y_order,], nrow = nn))
  Xy1 <- matrix(Xy0[Y_order,], nrow = nn)
  A1 <- A0[Y_order]
  # fitting the 2SLS model

  ## first-stage model: for jth entry of W, fit the corresponding model
  ## record the parameters, estimating equations, Jacobian, and number of
  ## regression coefficients and nuisance parameters (for negative binomial regression)
  param_1s <- U1j <- J1j <- vector("list", length = nW)

  nparam1_main <- nparam1_nuisance <- rep(NA, nW)
  ## predictor of W
  W_hat <- matrix(nrow = nn, ncol = nW)
  ## in the custom functions, make sure the nuisance parameters are before the regression
  ## coefficients
  for (j in 1:nW) {
    Wj <- W1[, j]
    Xwj <- Xw1[[j]]
    offset_j <- nco_args[[j]]$offset
    if (is.null(offset_j)) {
      offset_j <- rep(0, nn)
    }
    event_j <- nco_args[[j]]$event
    if (is.null(event_j)) {
      if (nco_type[j] == "ah") {
        warning("Event indicator for the NCO not specified -- assume no censoring.")
      }
      event_j <- rep(1, nn)
    }
    init_j <- nco_args[[j]]$init
    if (is.null(init_j)) {
      init_j <- NA
    }
    if (nco_type[j] == "linear") {
      W_model <- linear_fit(y = Wj, x = Xwj, offset = offset_j,
                            variance = T)
      ## no nuisance parameter
      param_1s[[j]] <- W_model$ESTIMATE
      U1j[[j]] <- W_model$EST_FUNC
      J1j[[j]] <- W_model$JACOBIAN
      nparam1_main[j] <- length(W_model$ESTIMATE)
      nparam1_nuisance[j] <- 0
      W_hat[, j] <- c(Xwj %*% param_1s[[j]])
    } else if (nco_type[j] == "loglin") {
      W_model <- loglin_fit(y = Wj, x = Xwj, offset = offset_j,
                            variance = T)

      ## no nuisance parameter
      param_1s[[j]] <- W_model$ESTIMATE
      U1j[[j]] <- W_model$EST_FUNC
      J1j[[j]] <- W_model$JACOBIAN
      nparam1_main[j] <- length(W_model$ESTIMATE)
      nparam1_nuisance[j] <- 0
      W_hat[, j] <- c(Xwj %*% param_1s[[j]])
    } else if (nco_type[j] == "poisson") {
      W_model <- poisson_fit(y = Wj, x = Xwj, offset = offset_j,
                             variance = T)

      ## no nuisance parameter
      param_1s[[j]] <- W_model$ESTIMATE
      U1j[[j]] <- W_model$EST_FUNC
      J1j[[j]] <- W_model$JACOBIAN
      nparam1_main[j] <- length(W_model$ESTIMATE)
      nparam1_nuisance[j] <- 0
      W_hat[, j] <- c(Xwj %*% param_1s[[j]])
    } else if (nco_type[j] == "ah") {
      W_model <- lin_ah(time = Wj, event = event_j,
                        covariates = Xwj, offset = offset_j)

      ## no nuisance parameter
      param_1s[[j]] <- W_model$ESTIMATE
      U1j[[j]] <- W_model$EST_FUNC
      J1j[[j]] <- W_model$JACOBIAN
      nparam1_main[j] <- length(W_model$ESTIMATE)
      nparam1_nuisance[j] <- 0
      W_hat[, j] <- c(Xwj %*% param_1s[[j]])
    } else if (nco_type[j] == "negbin") {
      W_model <- negbin_fit(y = Wj, x = Xwj, offset = offset_j,
                            variance = T, init = init_j)

      ## one nuisance parameter
      param_1s[[j]] <- W_model$ESTIMATE
      U1j[[j]] <- W_model$EST_FUNC
      J1j[[j]] <- W_model$JACOBIAN
      nparam1_main[j] <- length(W_model$ESTIMATE) - 1
      nparam1_nuisance[j] <- 1
      W_hat[, j] <- c(Xwj %*% param_1s[[j]][-1])
    }
  }


  ## second-stage model: fit the additive hazards model with the outcome against A, X and
  ## predictors of W
  ah_result <- lin_ah(time = t1, event = d1,
                      covariates = cbind(A1, Xy1, W_hat))
  param_2s <- ah_result$ESTIMATE

  params <- as.numeric(c(unlist(param_1s), param_2s))

  if (se_method %in% c("all", "analytic", "gaussian multiplier")) {

    np_s2 <- 1 + nW + ncol(Xy1) ## number of parameters in the second stage model
    par_W <- tail(params, nW)  ## parameters associated with W_hat

    U1 <- suppressMessages(as.matrix(dplyr::bind_cols(U1j)))
    J11 <- as.matrix(Matrix::bdiag(J1j))



    U2 <- ah_result$EST_FUNC# make predictors of W


    U <- cbind(U1, U2)


    ## Jacobian matrix
    ## make predictors of W
    S2 <- as.matrix(cbind(A1, Xy1, W_hat))

    J12 <- matrix(0, nrow = sum(nparam1_main) + sum(nparam1_nuisance), ncol = np_s2)

    J21_i <- array(0, dim = c(np_s2, sum(nparam1_main), nn))

    ## get incremental cumulative hazard and hazard function at each time point
    cumhaz <- ah_result$CUMHAZ_K
    haz <- ah_result$HAZ

    ## the last nW rows of Jacobian of M against gamma for each individual
    dW_i <- array(0, dim = c(nW, sum(nparam1_main), nn))

    for (i in 1:nn) {
      S1i <- lapply(Xw1, function(x) t(x[i,]))
      dW_i[, , i] <- as.matrix(Matrix::bdiag(S1i))
    }

    ## the last nW rows of Jacobian of \bar S2(T_{(k)}) against gamma for
    ## each individual at each unique event time
    dbW_Tk <- dsW_Tk <- array(0, dim = c(nW, sum(nparam1_main), ntime))


    for (k in ntime:1) {
      if (k == ntime) {
        dsW_Tk[, , k] <- apply(dW_i[, ,tmin[k]:nn, drop = F], c(1, 2), sum)
        dbW_Tk[, , k] <- dsW_Tk[, , k] / sum(o1 == ntime)
      } else {
        dsW_Tk[, , k] <- dsW_Tk[, , k + 1] +
          apply(dW_i[, ,tmin[k]:(tmin[k + 1] - 1), drop = F], c(1, 2), sum)
        dbW_Tk[, , k] <- dsW_Tk[, , k] / sum(o1 >= k)
      }
    }

    # integral of dbM_Tk by time
    int_dbW_Tk <- dbW_Tk

    for (k in 1:ntime) {
      int_dbW_Tk[, , k] <- int_dbW_Tk[, , k] * dtime[k]

      if (k >= 2) {
        int_dbW_Tk[, , k] <- int_dbW_Tk[, , k] + int_dbW_Tk[, , k - 1]
      }
    }

    ## Jacobian of the incremental cumulative hazard against gamma at each
    ## event time
    dhaz <- array(0, dim = c(1, sum(nparam1_main), ntime))

    for (k in 1:ntime) {
      dhaz[, , k] <- - par_W %*% dbW_Tk[, , k] * dtime[k]
    }

    ## average of S2 among the at-risk and its cumulative
    Zb_Tk <- matrix(NA, nrow = ntime, ncol = ncol(S2))
    for (k in 1:ntime) {
      Zb_Tk[k, ] <- colSums(matrix(S2[tmin[k]:nn,], ncol = ncol(S2))) /
        sum(o1 >= k)
    }

    int_Zb_Tk <- apply(Zb_Tk * dtime, 2, cumsum)

    ## Jacobian of the cumulative hazard against gamma at each
    ## event time

    dcumhaz <- dhaz
    for (k in 2:ntime) {
      dcumhaz[, , k] <- dcumhaz[, , k] + dcumhaz[, , k - 1]
    }

    ## integral of Zb_Tk with respect to dhaz and intergral of
    ## dbM_Tk with respect to haz
    int_dW_haz <- array(0, dim = c(nW, sum(nparam1_main), ntime))

    int_Zb_dhaz <- array(0, dim = c(ncol(S2), sum(nparam1_main), ntime))

    for (k in 1:ntime) {
      int_Zb_dhaz[, , k] <- Zb_Tk[k, ] %*% t(dhaz[, , k])
      int_dW_haz[, , k] <- dbW_Tk[, , k] * haz[k]

      if (k >= 2) {
        int_Zb_dhaz[, , k] <- int_Zb_dhaz[, , k] + int_Zb_dhaz[, , k - 1]
        int_dW_haz[, , k] <- int_dW_haz[, , k] + int_dW_haz[, , k - 1]
      }
    }




    ## Jacobian for every individual
    for (i in 1:nn) {
      dS2_i <- rbind(matrix(0, nrow = 1 + ncol(Xy1), ncol = sum(nparam1_main)),
                     dW_i[, , i])
      dbS2_Ti <- rbind(matrix(0, nrow = 1 + ncol(Xy1), ncol = sum(nparam1_main)),
                       dbW_Tk[, , o1[i]])

      J21_i1 <- dS2_i * (d1[i] - cumhaz[o1[i]] - c(param_2s %*% S2[i, ]) * t1[i])

      J21_i2 <- S2[i, ] %*% (-dcumhaz[, , o1[i]] - par_W %*% dW_i[, , i] * t1[i])

      J21_i3 <- d1[i] * dbS2_Ti

      J21_i4 <- rbind(matrix(0, nrow = 1 + ncol(Xy1), ncol = sum(nparam1_main)),
                      int_dW_haz[, , o1[i]])

      J21_i5 <- int_Zb_dhaz[, , o1[i]]
      J21_i6 <- int_Zb_Tk[o1[i], ] %*% t(par_W) %*% dW_i[, , i]
      J21_i7 <- c(param_2s %*% S2[i, ]) *
        rbind(matrix(0, nrow = 1 + ncol(Xy1), ncol = sum(nparam1_main)),
              int_dbW_Tk[, , o1[i]])

      J21_i[, , i] <- J21_i1 + J21_i2 - J21_i3 + J21_i4 + J21_i5 + J21_i6 +
        J21_i7
    }

    J21_main <- apply(J21_i, c(1, 2), sum)

    ## add columns of zero for nuisance parameter
    J21 <- NULL; nmain <- 0;

    for (w in 1:nW) {
      if (nparam1_nuisance[w] > 0) {
        J21 <- cbind(J21, matrix(0, nrow = np_s2, ncol = nparam1_nuisance[w]))
      }
      J21 <- cbind(J21, J21_main[, (nmain + 1):(nmain + nparam1_main[w])])
      nmain <- nmain + nparam1_main[w]
    }

    ## average derivative matrix among the at-risk at each time point


    J22 <- ah_result$JACOBIAN

    JJ <- rbind(cbind(J11, J12), cbind(J21, J22))

    if (verbose) {
        dout  <- svd(JJ)$d
        print(sprintf('Condition number of the Fisher information matrix is %1.1e', dout[1]/dout[length(dout)]))
    }
    Jinv <- solve(JJ)
  }



  ## parameter of interests are regression coefficients in the second stage
  ## model
  if (se_method %in% c("all", "analytic")) {
    DD <- t(U) %*% U

    VAR_analytic <- Jinv %*% DD %*% t(Jinv)
    SE_analytic <- tail(sqrt(diag(VAR_analytic)), np_s2)
  }

  if (se_method %in% c("all", "gaussian multiplier")) {
    IF <- nn * (U %*% t(Jinv))

    est_boot <- matrix(nrow = nboot,
                       ncol = sum(nparam1_main) + sum(nparam1_nuisance) + np_s2)
    for (b in 1:nboot) {
      est_boot[b, ] <- colMeans(IF * rnorm(nn))
    }

    SE_gaussian_multiplier <- tail(apply(est_boot, 2, sd), np_s2)
  }

  if (se_method %in% c("all", "exponential multiplier")) {
    est_boot <- matrix(nrow = nboot, ncol = np_s2)

    for (b in 1:nboot) {
      est_boot[b, ] <- lin_ah(time = Y, event = D,
                              covariates = S2,
                              weights = rexp(nn), variance = FALSE)$ESTIMATE
    }

    SE_exponential_multiplier <- apply(est_boot, 2, sd) ## 0.5864
  }

  if (se_method == "all") {
    SE <- rbind(SE_analytic, SE_gaussian_multiplier, SE_exponential_multiplier)
    rownames(SE) <- c("analytic", "gaussian multiplier",
                      "exponential multiplier")
  } else if (se_method == "analytic") {
    SE <- SE_analytic
  } else if (se_method == "gaussian multiplier") {
    SE <- SE_gaussian_multiplier
  } else if (se_method == "exponential multiplier") {
    SE <- SE_exponential_multiplier
  } else if (se_method == "none") {
    SE <- NULL
  }
  names(param_2s)  <-  c( 'A', sprintf('X%d', 1:ncol(Xy0)), sprintf('W%d', 1:ncol(W_hat)))

  return(list(ESTIMATE = param_2s,
              SE = SE,
              PARAM = unlist(c(param_1s, param_2s))))
}



