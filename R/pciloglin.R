#' Proximal Two-stage-least-squares with count outcomes
#'
#' This function computes the adjusted rate ratio using log-linear regression model,
#' using a pair of negative control variable
#'
#' @param Y an n-vector of count outcomes
#' @param offset an n-vector of offset for Y
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
#' @param variance whether to return the variance and components for the calculation
#' @returns A list with three elements: ESTIMATE includes the parameter estimates from the second-stage
#' model; SE produces their standard error, and PARAM includes all parameter estimates
#' @examples
#' N <- 2000
#' expit  <-  function(x) exp(x)/(1 + exp(x))
#' U <- runif(N); X <- runif(N)
#' A <- rbinom(N, 1, expit(-3 + 5 * U + 1 * X))
#' Y <- rpois(N, exp(0.2 + 1.5 * U + 0.2 * X + 0.2 * A))
#' Z <- rnorm(N, 2 * U + 0.5 * X)
#' W2 <- rexp(N, 0.1 + 1 * U + 0.1 * X)
#' D2 <- as.numeric(W2 < 5)
#' W2[D2 == 0] <- 5
#' W <- cbind(W1 = rnbinom(N, size = 25,
#'                    mu = exp(2.5 * U + 0.2 * X)),
#'            W2)
#' p2sls_result <- p2sls.loglin(Y = Y, A = A, X = X,
#'        W = W, Z = Z,
#'        nco_type = c("negbin", "ah"),
#'        nco_args = list(list(offset = rep(0, N)),
#'                        list(offset = rep(0, N),
#'                             event = D2)))
#' p2sls_result$summary_first_stage
#' p2sls_result$summary_second_stage
#' @export
p2sls.loglin <- function(Y, offset = rep(0, length(Y)),
                       A, X = NULL, W, Z = NULL, Xw = NULL,
                       Xy = NULL,
                       nco_type = NULL,
                       nco_args = NULL,
                       variance = TRUE, verbose = F) {



  nn <- length(Y)
  # clean data type
  Y0 <- as.numeric(Y)
  eta0 <- as.numeric(offset)
  W0 <- matrix(as.matrix(W), nrow = length(Y0))
  nW <- ncol(W0)

  ## add column names to W if needed
  if (nW == 1) {
    colnames(W0) <- "W"
  } else {
    if (is.null(colnames(W))) {
      colnames(W0) <- paste0("W", 1:nW)
    } else {
      colnames(W0) <- colnames(as.data.frame(W))
    }
  }

  if (!is.null(X)) {
    X0 <- as.matrix(X)
  } else {
    X0 <- X
  }

  nX <- ncol(X0)
  ## add column names to X if needed
  if (!is.null(X0)) {
    if (nX == 1) {
      colnames(X0) <- "X"
    } else if (nX > 0) {
      if (is.null(colnames(X))) {
        colnames(X0) <- paste0("X", 1:nX)
      } else {
        colnames(X0) <- colnames(as.data.frame(X))
      }
    }
  }

  A0 <- matrix(as.numeric(A), nrow = nn)
  if (ncol(A0) == 1){
    colnames(A0) <- "A"
  } else {
    colnames(A0) <- paste0("A", 1:ncol(A0))
  }

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

      nZ <- ncol(Z0)
      ## add column names to Z if needed
      if (!is.null(Z0)) {
        if (nZ == 1) {
          colnames(Z0) <- "Z"
        } else if (nZ > 0) {
          if (is.null(colnames(Z))) {
            colnames(Z0) <- paste0("Z", 1:nZ)
          } else {
            colnames(Z0) <- colnames(as.data.frame(Z))
          }
        }
      }

      paste_int <- function(a1, a2) paste(a1, a2, sep = ":")
      int_names <- function(names1, names2) c(outer(names1, names2, paste_int))
      int_mat <- function(a1, a2) {
        dplyr::bind_cols(apply(a2, 2, function(aa) a1 * aa, simplify = F))
      }

      if (is.null(X0)) {
        Xw0 <- lapply(1:nW, function(i) {  ## if no covariates, design matrix is A, Z and their interaction
          if (nco_type[i] != "ah") {
            xw <- as.matrix(cbind(1, A0, Z0, int_mat(Z0, A0)))
            colnames(xw) <- c("(Intercept)", colnames(A0), colnames(Z0),
                              int_names(colnames(Z0), colnames(A0)))
          } else {
            xw <- as.matrix(cbind(A0, Z0, int_mat(Z0, A0)))
            colnames(xw) <- c(colnames(A0), colnames(Z0), int_names(colnames(Z0), colnames(A0)))
          }
          return(xw)
        })
      } else {
        Xw0 <- lapply(1:nW, function(i) { ## with covariates, design matrix is A, Z, X and AZ, AX interactions
          if (nco_type[i] != "ah") {
            xw <- as.matrix(cbind(1, A0, X0, Z0, int_mat(X0, A0), int_mat(Z0, A0)))
            colnames(xw) <- c("(Intercept)", colnames(A0), colnames(X0), colnames(Z0),
                              int_names(colnames(X0), colnames(A0)),
                              int_names(colnames(Z0), colnames(A0)))
          } else {
            xw <- as.matrix(cbind(A0, X0, Z0, int_mat(X0, A0), int_mat(Z0, A0)))
            colnames(xw) <- c(colnames(A0), colnames(X0), colnames(Z0),
                              int_names(colnames(X0), colnames(A0)),
                              int_names(colnames(Z0), colnames(A0)))
          }
          return(xw)
        })
      }
    }
  } else {
    if (is.data.frame(Xw) | is.matrix(Xw)) {  ## same design matrix for every NCO
      Xw0 <- lapply(1:nW, function(i) {
        xw <- as.matrix(Xw)
        colnames(xw) <- colnames(Xw)
        return(xw)
      })
    } else {
      if (is.list(Xw)) {
        Xw0 <- lapply(Xw, function(xwi) {
          xw <- as.matrix(xwi)
          colnames(xw) <- colnames(xwi)
          return(xw)
        })
      }
    }
  }

  if (is.null(Xy)) {
    Xy0 <- X0
  } else {
    Xy0 <- matrix(as.matrix(Xy), nrow = nn)
    if (ncol(Xy0) == 1)  {
      colnames(Xy0) <- "X"
    } else {
      colnames(Xy0) <- paste0("X", 1:ncol(Xy0))
    }
  }

  if (!is.null(Z)) {
    Z0 <- as.matrix(Z)
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

  # fitting the 2SLS model

  ## first-stage model: for jth entry of W, fit the corresponding model
  ## record the parameters, estimating equations, Jacobian, and number of
  ## regression coefficients and nuisance parameters (for negative binomial regression)
  W_model <- param_1s <- U1j <- J1j <- vector("list", length = nW)

  nparam1_main <- nparam1_nuisance <- rep(NA, nW)
  ## predictor of W
  W_hat <- matrix(nrow = nn, ncol = nW)
  colnames(W_hat) <- colnames(W0)
  ## in the custom functions, make sure the nuisance parameters are before the regression
  ## coefficients
  for (j in 1:nW) {
    Wj <- W0[, j]
    Xwj <- Xw0[[j]]
    offset_j <- nco_args[[j]]$offset
    if (is.null(offset_j)) {
      offset_j <- rep(0, nn)
    }
    event_j <- nco_args[[j]]$event
    if (nco_type[j] == "ah" & is.null(event_j)) {
      warning("Event indicator for the NCO not specified -- assume no censoring.")
      event_j <- rep(1, nn)
    }
    init_j <- nco_args[[j]]$init
    if (is.null(init_j)) {
      init_j <- NA
    }

    if (nco_type[j] == "linear") {
      W_model[[j]] <- linear_fit(y = Wj, x = Xwj, offset = offset_j,
                                 variance = T)
      ## no nuisance parameter
      param_1s[[j]] <- W_model[[j]]$ESTIMATE
      U1j[[j]] <- W_model[[j]]$EST_FUNC
      J1j[[j]] <- W_model[[j]]$JACOBIAN
      nparam1_main[j] <- length(W_model[[j]]$ESTIMATE)
      nparam1_nuisance[j] <- 0
      W_hat[, j] <- c(Xwj %*% param_1s[[j]])
    } else if (nco_type[j] == "loglin") {
      W_model[[j]] <- loglin_fit(y = Wj, x = Xwj, offset = offset_j,
                                 variance = T)

      ## no nuisance parameter
      param_1s[[j]] <- W_model[[j]]$ESTIMATE
      U1j[[j]] <- W_model[[j]]$EST_FUNC
      J1j[[j]] <- W_model[[j]]$JACOBIAN
      nparam1_main[j] <- length(W_model[[j]]$ESTIMATE)
      nparam1_nuisance[j] <- 0
      W_hat[, j] <- c(Xwj %*% param_1s[[j]])
    } else if (nco_type[j] == "poisson") {
      W_model[[j]] <- poisson_fit(y = Wj, x = Xwj, offset = offset_j,
                                  variance = T)

      ## no nuisance parameter
      param_1s[[j]] <- W_model[[j]]$ESTIMATE
      U1j[[j]] <- W_model[[j]]$EST_FUNC
      J1j[[j]] <- W_model[[j]]$JACOBIAN
      nparam1_main[j] <- length(W_model[[j]]$ESTIMATE)
      nparam1_nuisance[j] <- 0
      W_hat[, j] <- c(Xwj %*% param_1s[[j]])
    } else if (nco_type[j] == "ah") {
      W_model[[j]] <- lin_ah(time = Wj, event = event_j,
                             covariates = Xwj, offset = offset_j)

      ## no nuisance parameter
      param_1s[[j]] <- W_model[[j]]$ESTIMATE
      U1j[[j]] <- W_model[[j]]$EST_FUNC
      J1j[[j]] <- W_model[[j]]$JACOBIAN
      nparam1_main[j] <- length(W_model[[j]]$ESTIMATE)
      nparam1_nuisance[j] <- 0
      W_hat[, j] <- c(Xwj %*% param_1s[[j]])
    } else if (nco_type[j] == "negbin") {
      W_model[[j]] <- negbin_fit(y = Wj, x = Xwj, offset = offset_j,
                                 variance = T, init = init_j)

      ## one nuisance parameter
      param_1s[[j]] <- W_model[[j]]$ESTIMATE
      U1j[[j]] <- W_model[[j]]$EST_FUNC
      J1j[[j]] <- W_model[[j]]$JACOBIAN
      nparam1_main[j] <- length(W_model[[j]]$ESTIMATE) - 1
      nparam1_nuisance[j] <- 1
      W_hat[, j] <- c(Xwj %*% param_1s[[j]][-1])
    }
  }


  ## second-stage model: fit the additive hazards model with the outcome against A, X and
  ## predictors of W
  loglin_result <- loglin_fit(y = Y0, x = cbind(1, A0, Xy0, W_hat), offset = eta0)
  param_2s <- loglin_result$ESTIMATE
  params <- as.numeric(c(unlist(param_1s), param_2s))


  if (verbose) {
    # The number of stage 1 parameters is only correct for when stage 1 is negative binomial and additive hazard models
    negbin.stage1.nparam  <-  (1 + 1 + ncol(A0) + ncol(X0) + ncol(Z0) + ncol(A0) * ncol(X0) + ncol(A0) * ncol(Z0))*sum(nco_type == 'negbin') # Each negbin has 1 nuisance + intercept + A + X + Z + AX + AZ
    ah.stage1.nparam  <- (ncol(A0) + ncol(X0) + ncol(Z0) +  ncol(A0) * ncol(X0) + ncol(A0) * ncol(Z0)) * sum(nco_type == 'ah')  # Each ah has  A + X + Z + AX + AZ (no intercept)
    stage2.nparam  <- (1 + ncol(A0) + ncol(Xy0) + ncol(W_hat)) # nuisance + intercept + A + X + W_hat
    print(sprintf('Number of params: Negbin stage 1: %d, AH stage 1: %d, stage2: %d, total: %d, check: %d', negbin.stage1.nparam, ah.stage1.nparam, stage2.nparam, negbin.stage1.nparam+ah.stage1.nparam+stage2.nparam, length(params)))
  }


  if (variance) {

    np_s2 <- length(param_2s) ## number of parameters in the second stage model
    par_W <- tail(params, nW)  ## parameters associated with W_hat

    U1 <- suppressMessages(as.matrix(dplyr::bind_cols(U1j)))
    J11 <- as.matrix(Matrix::bdiag(J1j))



    U2 <- loglin_result$EST_FUNC # make predictors of W


    U <- cbind(U1, U2)


    ## Jacobian matrix
    ## make predictors of W
    np_s2 <- length(param_2s)
    S2 <- as.matrix(cbind(1, A0, Xy0, W_hat))

    J12 <- matrix(0, nrow = sum(nparam1_main) + sum(nparam1_nuisance), ncol = np_s2)

    J21_i <- array(0, dim = c(np_s2, sum(nparam1_main), nn))

    ## the last nW rows of Jacobian of M against gamma for each individual
    dW_i <- array(0, dim = c(nW, sum(nparam1_main), nn))

    for (i in 1:nn) {
      S1i <- lapply(Xw0, function(x) t(x[i,]))
      dW_i[, , i] <- as.matrix(Matrix::bdiag(S1i))
    }
    ## Jacobian for every individual
    for (i in 1:nn) {
      dS2_i <- rbind(matrix(0, nrow = ifelse(is.null(Xy0), 1 + ncol(A0),
                                             1 + ncol(A0) + ncol(Xy0)),
                            ncol = sum(nparam1_main)),
                     dW_i[, , i])
      mu_i <- exp(eta0[i] + c(S2[i, ] %*% param_2s))

      J21_i[, , i] <- dS2_i * (Y0[i] - mu_i) - c(S2[i, ]) %*% t(param_2s) %*% dS2_i * mu_i
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


    J22 <- loglin_result$JACOBIAN

    JJ <- rbind(cbind(J11, J12), cbind(J21, J22))
    if (verbose) {
      dout  <- svd(JJ)$d
      print(sprintf('Condition number of the Fisher information matrix is %1.1e', dout[1]/dout[length(dout)]))
    }

    Jinv <- solve(JJ)

    DD <- t(U) %*% U

    VAR <- Jinv %*% DD %*% t(Jinv)
    all_se <- sqrt(diag(VAR))
    se_2s <- tail(sqrt(diag(VAR)), np_s2)

    summ_main <- matrix(NA, nrow = length(param_2s), ncol = 4)
    summ_main[, 1] <- param_2s
    summ_main[, 2] <- se_2s
    summ_main[, 3] <- param_2s / se_2s
    summ_main[, 4] <- pchisq((param_2s / se_2s) ^ 2, df = 1, lower.tail = F)


    colnames(summ_main) <- c("Estimate", "Std. Error", "z value",
                             "Pr(>|z|)")
    rownames(summ_main) <- c("(Intercept)", colnames(A0), colnames(Xy0), colnames(W_hat))
    summ_nuisance <- lapply(W_model, function(x) x$summary)
    if (nW > 1) names(summ_nuisance) <- paste0("W", 1:nW)

    return(list(ESTIMATE = param_2s,
                SE = se_2s,
                summary_first_stage = summ_nuisance,
                summary_second_stage = summ_main))
  }


  return(list(ESTIMATE = param_2s,
              PARAM = unlist(c(param_1s, param_2s))))
}
