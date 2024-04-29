#' Semiparametric log-linear model
#'
#' Fit semiparametric log-linear model E[Y | X] = exp(offset + beta * x)
#' @param y n-vector of outcome
#' @param x design matrix
#' @param offset n-vector of offset
#' @param variance whether to return the variance and components for the calculation
#' @export
loglin_fit <- function(y, x, offset = rep(0, length(y)),
                       variance = TRUE, interval = c(-100, 100)) {

  x <- matrix(x, nrow = length(y))

  ## estimate the regression coefficients
  U_loglin <- function(par) {
    uu <- x * c(y - exp(offset + x %*% par))
    return(uu)
  }

  est_func <- function(par) {
    sum(colMeans(U_loglin(par)) ^ 2)
  }


  if (ncol(x) == 1) {
    param <- optim(par = rep(0, ncol(x)), fn = est_func,
                   lower = interval[1], upper = interval[2],
                   method = "Brent")$par
  } else {
    param <- optim(par = rep(0, ncol(x)), fn = est_func)$par
  }

  PREDICT <- c(x %*% param)
  if (variance) {
    ## Jacobian matrix of the estimating function
    J <- - t(x) %*% (x * exp(c(offset) + c(x %*% param)))
    U <- U_loglin(param)
    Umat <- t(U) %*% U

    Jinv <- solve(J)

    VAR <- Jinv %*% Umat %*% t(Jinv)
    SE <- sqrt(diag(VAR))
    return(list(ESTIMATE = param,
                PREDICT = PREDICT,
                SE = SE,
                JACOBIAN = J,
                EST_FUNC = U,
                Umat = Umat))
  } else {
    return(list(ESTIMATE = param, PREDICT = PREDICT))
  }
}

#' Poisson regression model
#'
#' Fit Poisson regression model Y | X ~ Poisson(exp(offset + beta * x))
#' @param y n-vector of outcome
#' @param x design matrix
#' @param offset n-vector of offset
#' @param variance whether to return the variance and components for the calculation
#' @export
poisson_fit <- function(y, x, offset = rep(0, length(y)),
                        variance = TRUE, interval = c(-100, 100)) {

  x <- matrix(x, nrow = length(y))

  ## objective function: negative log-likelihood (omitting constants)
  negloglik <- function(par) {
    - sum(y * c(x %*% par) - exp(c(offset) + c(x %*% par)))
  }

  ## score function: gradient of the negative log likelihood
  negscore <- function(par) {
    - colSums(x * c(y - exp(c(offset) + c(x %*% par))))
  }

  if (ncol(x) == 1) {
    param <- optim(par = rep(0, ncol(x)), fn = negloglik,
                   lower = interval[1], upper = interval[2],
                   method = "Brent")$par
  } else {
    param <- optim(par = rep(0, ncol(x)), fn = negloglik,
                   gr = negscore)$par
  }

  PREDICT <- c(x %*% param)
  if (variance) {
    ## Jacobian matrix of the estimating function
    J <- - t(x) %*% (x * exp(c(offset) + c(x %*% param)))
    U <- x * c(y - exp(c(offset) + c(x %*% param)))
    Umat <- t(U) %*% U

    Jinv <- solve(J)

    VAR <- Jinv %*% Umat %*% t(Jinv)
    SE <- sqrt(diag(VAR))
    return(list(ESTIMATE = param,
                PREDICT = PREDICT,
                SE = SE,
                JACOBIAN = J,
                EST_FUNC = U,
                Umat = Umat))
  } else {
    return(list(ESTIMATE = param, PREDICT = PREDICT))
  }
}




#' Negative-binomial regression model
#'
#' Fit negative-binomial regression model Y | X ~ NegativeBinomial(theta, exp(offset + beta * x))
#' @param y n-vector of outcome
#' @param x design matrix
#' @param offset n-vector of offset
#' @param variance whether to return the variance and components for the calculation
#' @param init a vector of length ncol(x)+1, initial values for the negative-binomial regression.
#' Default is crude estimates of the result.
#' @export
negbin_fit <- function(y, x, offset = rep(0, length(y)),
                       variance = TRUE, interval = c(-100, 100),
                       init = NA) {
  nn <- length(y)
  x <- matrix(x, nrow = nn)

  ## objective function: negative log likelihood
  negloglik <- function(par) {
    thet <- par[1]; bb <- par[-1]
    mu <- exp(c(offset) + c(x %*% bb))

    loglik <- rep(NA, nn)

    for (i in 1:nn) {
      loglik[i] <- (y[i] > 0) * ifelse(y[i] == 0, 0, sum(log(0:(y[i] - 1) + thet))) -
        (y[i] > 0) * ifelse(y[i] == 0, 0, sum(log(1:y[i]))) +
        thet * log(thet) +
        y[i] * log(mu[i]) -
        (thet + y[i]) * log(thet + mu[i])
    }

    return(-sum(loglik))
  }

  ## negative log likelihood for theta, with fixed beta
  negloglik_thet <- function(thet) {
    negloglik(c(thet, beta_init))
  }

   ## score equation
  U_loglin <- function(par) {
    thet <- par[1]; bb <- par[-1]
    mu <- exp(c(offset) + c(x %*% bb))
    u_thet <- rep(NA, nn)

    for (i in 1:nn) {
      u_thet[i] <- log(thet / (thet + mu[i])) +
        ifelse(y[i] == 0, 0, sum(1 / (0:(y[i] - 1) + thet))) +
        (mu[i] - y[i]) / (thet + mu[i])
    }

    u_b <- x * c((y - mu) *(thet / (thet + mu)))

    uu <- cbind(u_thet, u_b)
    return(uu)
  }

  negscore <- function(par) {
    U <- U_loglin(par)
    return(-colSums(U))
  }


  est_func <- function(par) {
    sum(colMeans(U_loglin(par)) ^ 2)
  }
  if (is.na(init)) {
    ## initial value: by log-linear model
    beta_init <- loglin_fit(y = y, x = x, offset = offset, variance = FALSE)$ESTIMATE
    thet_init <- optimize(negloglik_thet, c(0.001, 100))$minimum

    param_init <- optim(par = c(thet_init, beta_init),
                        fn = est_func,
                        lower = c(0.001, rep(-Inf, ncol(x))),
                        method = "L-BFGS-B")$par
  } else {
    init <- init
  }
  param <- optim(par = param_init,
                 fn = negloglik, gr = negscore,
                 lower = c(0.001, rep(-Inf, ncol(x))),
                 method = "L-BFGS-B")$par


  thet <- param[1]; bb <- param[-1]; mu <- exp(c(offset) + c(x %*% bb))
  PREDICT <- c(x %*% bb)
  if (variance) {
    ## Jacobian matrix of the estimating function
    J <- matrix(NA, nrow = 1 + ncol(x), ncol = 1 + ncol(x))

    J11_i <- rep(NA, nn);

    for (i in 1:nn) {
      J11_i[i] <- ifelse(y[i] == 0, 0, - sum(1 / (0:(y[i] - 1) + thet) ^ 2)) +
        (mu[i] ^ 2 + thet * y[i]) / (thet * (thet + mu[i]) ^ 2)
    }

    J11 <- sum(J11_i)
    J21 <- t(x) %*% c(mu * (y - mu) / (thet + mu) ^ 2)
    J12 <- t(J21)
    J22 <- - t(x) %*% (x * c(thet * mu * (thet + y) / (thet + mu) ^ 2))
    J <- rbind(cbind(J11, J12), cbind(J21, J22))
    U <- U_loglin(param)
    Umat <- t(U) %*% U

    Jinv <- solve(J)

    VAR <- Jinv %*% Umat %*% t(Jinv)
    SE <- sqrt(diag(VAR))
    return(list(ESTIMATE = param,
                PREDICT = PREDICT,
                SE = SE,
                JACOBIAN = J,
                EST_FUNC = U,
                Umat = Umat))
  } else {
    return(list(ESTIMATE = param, PREDICT = PREDICT))
  }
}
