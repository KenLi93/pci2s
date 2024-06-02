#' Proximal Two-stage-least-squares with survival outcome
#'
#' This function computes the adjusted hazard difference using additive hazards model by
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
#' @param variance whether to return the variance and components for the calculation
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
#' W <- cbind(W1 = rnbinom(N, size = 25, mu = exp(2.5 * U + 0.2 * X)),
#'            W2)
#' p2sls_result <- p2sls.ah(Y = Y, D = D, A = A, X = X,
#'        W = W, Z = Z, variance = TRUE,
#'        nco_type = c("negbin", "ah"),
#'        nco_args = list(list(offset = rep(0, N)),
#'                        list(offset = rep(0, N),
#'                             event = D2)))
#' p2sls_result$summary_first_stage
#' p2sls_result$summary_second_stage
#'
#' @export
p2sls.ah <- function(Y, D, A, X = NULL, W, Z = NULL, Xw = NULL,
                    Xy = NULL,
                    nco_type = NULL,
                    nco_args = NULL,
                    variance = TRUE, verbose = FALSE) {



  nn <- length(Y)
  # clean data type
  Y0 <- as.numeric(Y)
  D0 <- as.numeric(D)
  W0 <- matrix(as.matrix(as.data.frame(W)), nrow = length(Y0))
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


  W1 <- W0[Y_order, , drop = F]
  Xw1 <- lapply(Xw0, function(xx) xx[Y_order, , drop = F])
  Xy1 <- Xy0[Y_order, , drop = F]
  A1 <- A0[Y_order, , drop = F]
  # fitting the 2SLS model

  ## first-stage model: for jth entry of W, fit the corresponding model
  ## record the parameters, estimating equations, Jacobian, and number of
  ## regression coefficients and nuisance parameters (for negative binomial regression)
  W_model <- param_1s <- U1j <- J1j <- vector("list", length = nW)

  nparam1_main <- nparam1_nuisance <- rep(NA, nW)
  ## predictor of W
  W_hat <- matrix(nrow = nn, ncol = nW)
  colnames(W_hat) <- colnames(W1)
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
  ah_result <- lin_ah(time = t1, event = d1,
                      covariates = cbind(A1, Xy1, W_hat))
  param_2s <- ah_result$ESTIMATE

  params <- as.numeric(c(unlist(param_1s), param_2s))

  if (variance) {

    np_s2 <- ncol(A1) + nW + ncol(Xy1) ## number of parameters in the second stage model
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
      dS2_i <- rbind(matrix(0, nrow = ifelse(is.null(Xy1), ncol(A1), ncol(A1) + ncol(Xy1)), ncol = sum(nparam1_main)),
                     dW_i[, , i])
      dbS2_Ti <- rbind(matrix(0, nrow = ifelse(is.null(Xy1), ncol(A1), ncol(A1) + ncol(Xy1)), ncol = sum(nparam1_main)),
                       dbW_Tk[, , o1[i]])

      J21_i1 <- dS2_i * (d1[i] - cumhaz[o1[i]] - c(param_2s %*% S2[i, ]) * t1[i])

      J21_i2 <- S2[i, ] %*% (-dcumhaz[, , o1[i]] - par_W %*% dW_i[, , i] * t1[i])

      J21_i3 <- d1[i] * dbS2_Ti

      J21_i4 <- rbind(matrix(0, nrow = ifelse(is.null(Xy1), ncol(A1), ncol(A1) + ncol(Xy1)), ncol = sum(nparam1_main)),
                      int_dW_haz[, , o1[i]])

      J21_i5 <- int_Zb_dhaz[, , o1[i]]
      J21_i6 <- int_Zb_Tk[o1[i], ] %*% t(par_W) %*% dW_i[, , i]
      J21_i7 <- c(param_2s %*% S2[i, ]) *
        rbind(matrix(0, nrow = ifelse(is.null(Xy1), ncol(A1), ncol(A1) + ncol(Xy1)), ncol = sum(nparam1_main)),
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
      print(sprintf('Condition number of the Fisher information matrix is %1.1e', dout[1] / dout[length(dout)]))
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
    rownames(summ_main) <- c(colnames(A1), colnames(Xy1), colnames(W_hat))
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



