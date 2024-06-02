#' Proximal causal inference with categorical outcomes using two-stage-least-square
#'
#' This function computes the adjusted odds ratio using multiple logistic regression model,
#' using a pair of negative control variable
#'
#' @param Y an n-vector of categorical outcomes with nY distinct values
#' @param offset an n-vector of offset for Y
#' @param A an n-vector of primary exposure
#' @param X an n * nX matrix of adjusted covariates , can be empty
#' @param W an n-vector of categorical NCO with nW distinct values
#' @param Z an n * nZ matrix of negative control exposure (NCE). Alternatively,
#' the NCE can be included in the design matrix of the first-stage model by specifying
#' Xw
#' @param Xw a matrix with n rows or a list of length nW, each element is a design matrix of
#' the corresponding category of W. Default is a matrix with columns of A, Z, X, and Y
#' @param Xw_Y0 a matrix with n rows or a list of length nW, each element is a design matrix of
#' the corresponding category of W, with the value of Y fixed at zero. Default is a matrix with columns of A, Z, X,
#' and Y replaced by zero.
#' @param Xy a matrix with n rows or a list of length nY, each element is a design matrix of
#' the corresponding category of Y besides A, W and linear predictors of W. Default to be X
#' @param nco_offset An n-vector of offset for W.
#' @param variance whether to return the variance and components for the calculation
#' @returns A list with three elements: ESTIMATE includes the parameter estimates from the second-stage
#' model; SE produces their standard error, and PARAM includes all parameter estimates
#' @examples
#' set.seed(2002)
#' N <- 2000
#' expit  <-  function(x) exp(x)/(1 + exp(x))
#' U <- runif(N, -1, 1); X <- runif(N, -1, 1)
#' A <- rbinom(N, 1, expit(-3 + 5 * U + 1 * X))
#' Y <- rep(NA, N)
#' for (i in 1:N) {
#'   Yi <- rmultinom(1, 1, c(1,
#'                           exp(0.2 + 1.5 * U[i] + 0.2 * X[i] + 0.2 * A[i]),
#'                           exp(0.1 + 3 * U[i] + 0.1 * X[i] + 0.4 * A[i])))
#'   Y[i] <- which(Yi == 1) - 1
#' }
#'
#' W <- rep(NA, N)
#' for (i in 1:N) {
#'   W[i] <- rbinom(1, 1, expit(0.1 + 1 * U[i] + 0.1 * X[i]))
#' }
#' Z <- rnorm(N, 2 * U + 0.5 * X)
#'
#' pci_result <- pci.logitreg(Y = Y, A = A, X = X,
#'        W = W, Z = Z)
#' pci_result$summary_first_stage
#' pci_result$summary_second_stage
#' @export
pci.logitreg <- function(Y, offset = rep(0, length(Y)),
                       A, X = NULL, W, Z = NULL, Xw = NULL,
                       Xw_Y0 = NULL,
                       Xy = NULL,
                       nco_offset = rep(0, length(Y)),
                       variance = TRUE) {



  nn <- length(Y)

  # clean data type
  Y_levels <- levels(as.factor(Y))
  nY <- length(Y_levels)
  Y0 <- matrix(NA, nrow = nn, ncol = nY - 1)

  for (i in 1:(nY - 1)) {
    Y0[, i] <- as.numeric(Y == Y_levels[i + 1])
  }

  Y00 <- Y0 * 0
  ## add column names to Y if needed
  if (nY == 2) {
    colnames(Y0) <- "Y"
  } else {
    colnames(Y0) <- paste0("Y_", Y_levels[-1])
  }

  eta0 <- as.numeric(offset)

  W_levels <- levels(as.factor(W))
  nW <- length(W_levels)
  W0 <- matrix(NA, nrow = nn, ncol = nW - 1)

  for (i in 1:(nW - 1)) {
    W0[, i] <- as.numeric(W == W_levels[i + 1])
  }

  etaW <- as.numeric(nco_offset)

  ## add column names to W if needed
  if (nW == 2) {
    colnames(W0) <- "W"
  } else {
    colnames(W0) <- paste0("W_", W_levels[-1])
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

      if (is.null(X0)) {
        Xw0 <- lapply(1:(nW - 1), function(i) {  ## if no covariates, design matrix is A, Z and their interaction

          xw <- as.matrix(cbind(1, A0, Z0, Y0))
          colnames(xw) <- c("(Intercept)", colnames(A0), colnames(Z0), colnames(Y0))
          return(xw)
        })

        Xw00 <- lapply(1:(nW - 1), function(i) {  ## if no covariates, design matrix is A, Z and their interaction

          xw <- as.matrix(cbind(1, A0, Z0, Y00))
          colnames(xw) <- c("(Intercept)", colnames(A0), colnames(Z0), colnames(Y0))
          return(xw)
        })
      } else {
        Xw0 <- lapply(1:(nW - 1), function(i) { ## with covariates, design matrix is A, Z, X and AZ, AX interactions
          xw <- as.matrix(cbind(1, A0, X0, Z0, Y0))
          colnames(xw) <- c("(Intercept)", colnames(A0), colnames(X0), colnames(Z0), colnames(Y0))
          return(xw)
        })

        Xw00 <- lapply(1:(nW - 1), function(i) { ## with covariates, design matrix is A, Z, X and AZ, AX interactions
          xw <- as.matrix(cbind(1, A0, X0, Z0, Y00))
          colnames(xw) <- c("(Intercept)", colnames(A0), colnames(X0), colnames(Z0), colnames(Y0))
          return(xw)
        })
      }
    }
  } else {
    if (is.null(Xw_Y0)) {
      stop("If Xw is given, then Xw_Y0 also need to be specified.")
    }

    if (is.data.frame(Xw) | is.matrix(Xw)) {  ## same design matrix for every NCO
      Xw0 <- lapply(1:(nW - 1), function(i) {
        xw <- as.matrix(Xw)
        colnames(xw) <- colnames(Xw)
        return(xw)
      })
    } else {
      if (is.list(Xw)) {
        Xw0 <- lapply(Xw, function(xwi) {
          xw <- as.matrix(xwi)
          colnames(xw) <- colnames(xwi)
        })
      }
    }

    if (is.data.frame(Xw_Y0) | is.matrix(Xw_Y0)) {  ## same design matrix for every NCO
      Xw00 <- lapply(1:(nW - 1), function(i) {
        xw <- as.matrix(Xw_Y0)
        colnames(xw) <- colnames(Xw)
        return(xw)
      })
    } else {
      if (is.list(Xw_Y0)) {
        Xw00 <- lapply(Xw_Y0, function(xwi) {
          xw <- as.matrix(xwi)
          colnames(xw) <- colnames(xwi)
        })
      }
    }
  }


  if (is.null(Xy)) {
    Xy0 <- lapply(1:(nY - 1), function(i) {  ## if no covariates, the remaining part in the design matrix is W
      xy <- X0
      return(xy)
    })
  } else {
    if (is.data.frame(Xy) | is.matrix(Xy)) {  ## same design matrix for every NCO
      Xy0 <- lapply(1:(nY - 1), function(i) {
        xy <- as.matrix(Xy)
        colnames(xy) <- colnames(Xy)
        return(xy)
      })
    } else {
      if (is.list(Xy)) {
        Xy0 <- lapply(Xy, function(xyi) {
          xy <- as.matrix(xwyi)
          colnames(xy) <- colnames(xyi)
          return(xy)
        })
      }
    }
  }

  # fitting the 2SLS model

  ## first-stage model: fit a multiple logistic regression model for W
  W_model <- param_1s <- U1j <- J1j <- vector("list", length = nW - 1)

  nparam1_main <- rep(NA, nW - 1)
  ## linear predictor of W
  W_hat <- matrix(nrow = nn, ncol = nW - 1)
  colnames(W_hat) <- paste0(colnames(W0), "_hat")
  ## in the custom functions, make sure the nuisance parameters are before the regression
  ## coefficients
  for (j in 1:(nW - 1)) {
    Wj <- W0[, j]
    Xwj <- Xw0[[j]]
    Xw00j <- Xw00[[j]]
    offset_j <- etaW
    if (is.null(offset_j)) {
      offset_j <- rep(0, nn)
    }


    W_model[[j]] <- logit_reg(y = Wj, x = Xwj, offset = offset_j,
                              variance = T)

    ## no nuisance parameter
    param_1s[[j]] <- W_model[[j]]$ESTIMATE
    U1j[[j]] <- W_model[[j]]$EST_FUNC
    J1j[[j]] <- W_model[[j]]$JACOBIAN
    nparam1_main[j] <- length(W_model[[j]]$ESTIMATE)
    W_hat[, j] <- c(Xw00j %*% param_1s[[j]])

  }


  ## second-stage model: for each non-reference category of Y, fit the
  ## logistic regression model with the outcome against A, X, W, and
  ## predictors of W
  if (variance) {
    summary_second_stage <- vector("list", length = nY - 1)
  }

  for (k in 1:(nY - 1)) {
    Yk <- Y0[, k]
    logit_reg_result <- logit_reg(y = Yk, x = cbind(1, A0, Xy0[[k]], W0, W_hat), offset = eta0)
    param_2s <- logit_reg_result$ESTIMATE
    params <- as.numeric(c(unlist(param_1s), param_2s))


    if (variance) {

      np_s2 <- length(param_2s) ## number of parameters in the second stage model
      par_W <- tail(params, nW - 1)  ## parameters associated with W_hat

      U1 <- suppressMessages(as.matrix(dplyr::bind_cols(U1j)))
      J11 <- as.matrix(Matrix::bdiag(J1j))



      U2 <- logit_reg_result$EST_FUNC # make predictors of W


      U <- cbind(U1, U2)


      ## Jacobian matrix
      ## make predictors of W
      np_s2 <- length(param_2s)
      S2 <- as.matrix(cbind(1, A0, Xy0[[k]], W0, W_hat))

      J12 <- matrix(0, nrow = sum(nparam1_main), ncol = np_s2)

      J21_i <- array(0, dim = c(np_s2, sum(nparam1_main), nn))

      ## the last nW rows of Jacobian of M against gamma for each individual
      dW_i <- array(0, dim = c(nW - 1, sum(nparam1_main), nn))

      for (i in 1:nn) {
        S1i <- lapply(Xw00, function(x) t(x[i,]))
        dW_i[, , i] <- as.matrix(Matrix::bdiag(S1i))
      }
      ## Jacobian for every individual
      for (i in 1:nn) {
        dS2_i <- rbind(matrix(0, nrow = 1 + ncol(A0) + ncol(Xy0[[k]]) + ncol(W0), ncol = sum(nparam1_main)),
                       dW_i[, , i])
        mu_i <- expit(eta0[i] + c(S2[i, ] %*% param_2s))

        J21_i[, , i] <- dS2_i * (Y0[i] - mu_i) -
          c(S2[i, ]) %*% t(param_2s) %*% dS2_i *
          exp(eta0[i] + c(S2[i, ] %*% param_2s)) / (1 + exp(eta0[i] + c(S2[i, ] %*% param_2s))) ^ 2
      }

      J21 <- apply(J21_i, c(1, 2), sum)


      ## average derivative matrix among the at-risk at each time point


      J22 <- logit_reg_result$JACOBIAN

      JJ <- rbind(cbind(J11, J12), cbind(J21, J22))

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
      rownames(summ_main) <- c("(Intercept)", colnames(A0), colnames(Xy0[[k]]),
                               colnames(W0), colnames(W_hat))
      summary_second_stage[[k]] <- summ_main
    }
  }

  if (variance) {
    summ_nuisance <- lapply(W_model, function(x) x$summary)
    names(summ_nuisance) <- colnames(W0)
    names(summary_second_stage) <- colnames(Y0)

    return(list(ESTIMATE = param_2s,
                SE = se_2s,
                summary_first_stage = summ_nuisance,
                summary_second_stage = summary_second_stage))
  }


  return(list(ESTIMATE = as.numeric(param_2s),
              PARAM = as.numeric(unlist(c(param_1s, param_2s)))))
}
