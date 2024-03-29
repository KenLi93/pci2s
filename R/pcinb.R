#' Proximal causal inference with count outcomes using two-stage-least-square
#'
#' 'pcinb2s()' computes the adjusted rate ratio using negative-binomial regression model,
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
#' @param se_method Method to compute the standard error, can be "analytic", "gaussian multiplier" or "none".
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
#' W <- cbind(rnbinom(N, size = 25,
#'                    mu = exp(2.5 * U + 0.2 * X)),
#'            W2)
#' pcinb2s(Y = Y, A = A, X = X,
#'        W = W, Z = Z, nboot = 100,
#'        nco_type = c("negbin", "ah"),
#'        nco_args = list(list(offset = rep(0, N)),
#'                        list(offset = rep(0, N),
#'                             event = D2)))
#' @export
pcinb2s <- function(Y, offset = rep(0, length(Y)),
                    A, X = NULL, W, Z = NULL, Xw = NULL,
                    Xy = NULL,
                    nboot = 2000,
                    nco_type = NULL,
                    nco_args = NULL,
                    nb_init = NA,
                    se_method = "analytic", verbose=F) {



  nn <- length(Y)
  # clean data type
  Y0 <- as.numeric(Y)
  eta0 <- as.numeric(offset)
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
    Xy0 <- matrix(cbind(1, as.matrix(X)), nrow = nn)
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
  negbin_result <- negbin_fit(y = Y0, x = cbind(A0, Xy0, W_hat),
                              init = nb_init, offset = eta0)
  param_2s <- negbin_result$ESTIMATE
  param_2_theta <- param_2s[1] ## size parameter of the negative binomial regression
  param_2_beta <- param_2s[-1] ## regression coefficients
  params <- as.numeric(c(unlist(param_1s), param_2s))


  if (verbose) {
      # The number of stage 1 parameters is only correct for when stage 1 is negative binomial and additive hazard models
      negbin.stage1.nparam  <-  (1+1+ 1+ncol(X0)+ ncol(Z0)+ncol(X0)+ ncol(Z0))*sum(nco_type=='negbin') # Each negbin has 1 nuisance + intercept + A + X + Z + AX + AZ
      ah.stage1.nparam  <- (1+ncol(X0)+ ncol(Z0)+ncol(X0)+ ncol(Z0))*sum(nco_type=='ah')  # Each ah has  A + X + Z + AX + AZ (no intercept)
      stage2.nparam  <-  (1+ 1 + ncol(Xy0) + ncol(W_hat)) # nuisance + intercept + A + X + W_hat
      print(sprintf('Number of params: Negbin stage 1: %d, AH stage 1: %d, stage2: %d, total: %d, check: %d', negbin.stage1.nparam, ah.stage1.nparam, stage2.nparam, negbin.stage1.nparam+ah.stage1.nparam+stage2.nparam, length(params)))
  }


  if (se_method %in% c("analytic", "gaussian multiplier")) {

    np_s2 <- length(param_2s) ## number of parameters in the second stage model
    par_W <- tail(params, nW)  ## parameters associated with W_hat

    U1 <- suppressMessages(as.matrix(dplyr::bind_cols(U1j)))
    J11 <- as.matrix(Matrix::bdiag(J1j))



    U2 <- negbin_result$EST_FUNC # make predictors of W


    U <- cbind(U1, U2)


    ## Jacobian matrix
    ## make predictors of W
    np_s2 <- length(param_2s)
    S2 <- as.matrix(cbind(A0, Xy0, W_hat))

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
      dS2_i <- rbind(matrix(0, nrow = 1 + ncol(Xy0), ncol = sum(nparam1_main)),
                     dW_i[, , i])
      mu_i <- exp(eta0[i] + c(S2[i, ] %*% param_2_beta))
        # if (mu_i > 100){
        #       print(sprintf('%d, %.0e',i, mu_i))
    # }
      J21_i1 <- (param_2_theta + mu_i) / param_2_theta * mu_i *
        t(param_2_beta) %*% dS2_i +
        (param_2_theta + Y0[i]) / (param_2_theta + mu_i) ^ 2 * mu_i *
        t(param_2_beta) %*% dS2_i

      J21_i2 <- (Y0[i] - mu_i) / (param_2_theta + mu_i) * param_2_theta * dS2_i -
        param_2_theta * mu_i / (param_2_theta + mu_i) *
        c(S2[i, ]) %*% t(param_2_beta) %*% dS2_i -
        (Y0[i] - mu_i) / (param_2_theta + mu_i) ^ 2 * param_2_theta * mu_i *
        c(S2[i, ]) %*% t(param_2_beta) %*% dS2_i

      J21_i[, , i] <- rbind(J21_i1, J21_i2)
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


    J22 <- negbin_result$JACOBIAN

    JJ <- rbind(cbind(J11, J12), cbind(J21, J22))
    if (verbose) {
        dout  <- svd(JJ)$d
        print(sprintf('Condition number of the Fisher information matrix is %1.1e', dout[1]/dout[length(dout)]))
    }

    Jinv <- solve(JJ)
  }



  names(param_2s)  <-  c('size', 'A', sprintf('X%d', 1:ncol(Xy0)), sprintf('W%d', 1:ncol(W_hat)))
  ## parameter of interests are regression coefficients in the second stage
  ## model
  if (se_method == "analytic") {
    DD <- t(U) %*% U

    VAR_analytic <- Jinv %*% DD %*% t(Jinv)
    SE <- tail(sqrt(diag(VAR_analytic)), np_s2)
    names(SE)  <-  names(param_2s)
  }

  if (se_method %in% c("gaussian multiplier")) {
    IF <- nn * (U %*% t(Jinv))

    est_boot <- matrix(nrow = nboot,
                       ncol = sum(nparam1_main) + sum(nparam1_nuisance) + np_s2)
    for (b in 1:nboot) {
      est_boot[b, ] <- colMeans(IF * rnorm(nn))
    }

    SE <- tail(apply(est_boot, 2, sd), np_s2)
    names(SE)  <-  names(param_2s)
  }

  if (se_method == "none") {
    SE <- NULL
  }


  return(list(ESTIMATE = param_2s,
              SE = SE,
              PARAM = unlist(c(param_1s, param_2s))))
}
