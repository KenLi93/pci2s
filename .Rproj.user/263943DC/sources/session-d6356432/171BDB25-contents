#' Semiparametric log-linear model
#'
#' Fit semiparametric log-linear model E[Y | X] = exp(offset + beta * x).
#' This is a utility function. For formal log-linear regression models,
#' the glm function is recommended
#' @param y n-vector of outcome
#' @param x design matrix
#' @param offset n-vector of offset
#' @param variance whether to return the variance and components for the calculation
#' @examples
#' x <- cbind(1, rep(1:3, each = 10))
#' tt <- rep(1:2, 15)
#' y <- rpois(30, lambda = tt * exp(x[, 1] + 0.5 * x[, 2]))
#'
#' loglin_result <- loglin_fit(y, x, offset = log(tt))
#' loglin_result$summary
#' @export
loglin_fit <- function(y, x, offset = rep(0, length(y)),
                       variance = TRUE, interval = c(-100, 100)) {

  if (class(x)[[1]] %in% c("matrix", "data.frame", "array")) {
    x <- as.matrix(x)
  } else {
    x <- matrix(x, nrow = nn)
  }

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

    summ <- matrix(NA, nrow = length(param), ncol = 4)
    summ[, 1] <- param
    summ[, 2] <- SE
    summ[, 3] <- param / SE
    summ[, 4] <- pchisq((param / SE) ^ 2, df = 1, lower.tail = F)
    if (is.null(colnames(x))) {
      if (ncol(x) == 1) {
        rownames(summ) <- "x"
      } else {
        rownames(summ) <- paste0("x", 1:ncol(x))
      }
    } else {
      rownames(summ) <- colnames(x)
    }

    colnames(summ) <- c("Estimate", "Std. Error", "z value",
                        "Pr(>|z|)")
    return(list(ESTIMATE = param,
                PREDICT = PREDICT,
                SE = SE,
                JACOBIAN = J,
                EST_FUNC = U,
                Umat = Umat,
                summary = summ))
  } else {
    return(list(ESTIMATE = param, PREDICT = PREDICT))
  }
}

#' Poisson regression model
#'
#' Fit Poisson regression model Y | X ~ Poisson(exp(offset + beta * x))
#' This is a utility function. For formal Poisson regression models,
#' the glm function is recommended
#' @param y n-vector of outcome
#' @param x design matrix
#' @param offset n-vector of offset
#' @param variance whether to return the variance and components for the calculation
#' @examples
#' x <- cbind(1, rep(1:3, each = 10))
#' tt <- rep(1:2, 15)
#' y <- rpois(30, lambda = tt * exp(x[, 1] + 0.5 * x[, 2]))
#'
#' poisreg_result <- poisson_fit(y, x, offset = log(tt))
#' poisreg_result$summary
#' @export
poisson_fit <- function(y, x, offset = rep(0, length(y)),
                        variance = TRUE, interval = c(-100, 100)) {

  if (class(x)[[1]] %in% c("matrix", "data.frame", "array")) {
    x <- as.matrix(x)
  } else {
    x <- matrix(x, nrow = nn)
  }

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

    summ <- matrix(NA, nrow = length(param), ncol = 4)
    summ[, 1] <- param
    summ[, 2] <- SE
    summ[, 3] <- param / SE
    summ[, 4] <- pchisq((param / SE) ^ 2, df = 1, lower.tail = F)
    if (is.null(colnames(x))) {
      if (ncol(x) == 1) {
        rownames(summ) <- "x"
      } else {
        rownames(summ) <- paste0("x", 1:ncol(x))
      }
    } else {
      rownames(summ) <- colnames(x)
    }

    colnames(summ) <- c("Estimate", "Std. Error", "z value",
                        "Pr(>|z|)")

    return(list(ESTIMATE = param,
                PREDICT = PREDICT,
                SE = SE,
                JACOBIAN = J,
                EST_FUNC = U,
                Umat = Umat,
                summary = summ))
  } else {
    return(list(ESTIMATE = param, PREDICT = PREDICT))
  }
}




#' Negative-binomial regression model
#'
#' Fit negative-binomial regression model Y | X ~ NegativeBinomial(theta, exp(offset + beta * x)).
#' This is a utility function. For formal negative-binomial regression models,
#' MASS::glm.nb is recommended
#'
#' @param y n-vector of outcome
#' @param x design matrix
#' @param offset n-vector of offset
#' @param variance whether to return the variance and components for the calculation
#' @param init a vector of length ncol(x)+1, initial values for the negative-binomial regression.
#' Default is crude estimates of the result.
#' @examples
#' x <- cbind(1, rep(1:3, each = 20))
#' tt <- rep(1:2, 30)
#' y <- rnbinom(60, size = 2, mu = tt * exp(x[, 1] + 0.5 * x[, 2]))
#'
#' negbin_result <- negbin_fit(y, x, offset = log(tt))
#' negbin_result$summary
#' @export
negbin_fit <- function(y, x, offset = rep(0, length(y)),
                       variance = TRUE, interval = c(-100, 100),
                       init = NA) {
  nn <- length(y)
  if (class(x)[[1]] %in% c("matrix", "data.frame", "array")) {
    x <- as.matrix(x)
  } else {
    x <- matrix(x, nrow = nn)
  }

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
  if (length(init) == 1 && is.na(init)) {
    ## initial value: by log-linear model
    ## Use the MASS::glm.nb function
    nb_model <- MASS::glm.nb(y ~ 0 + x + stats::offset(offset))
    beta_init <- nb_model$coefficients[1:ncol(x)]
    thet_init <- nb_model$theta

    param_init <- optim(par = c(thet_init, beta_init),
                        fn = est_func,
                        lower = c(0.001, rep(-Inf, ncol(x))),
                        method = "L-BFGS-B")$par
  } else {
      if (length(init) != ncol(x) + 1) {
        stop(sprintf("length(init) must be %d", ncol(x) + 1))
      }
    param_init <- init
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

    summ <- matrix(NA, nrow = length(param), ncol = 4)
    summ[, 1] <- param
    summ[, 2] <- SE
    summ[, 3] <- param / SE
    summ[, 4] <- pchisq((param / SE) ^ 2, df = 1, lower.tail = F)
    if (is.null(colnames(x))) {
      if (ncol(x) == 1) {
        rownames(summ) <- "x"
      } else {
        rownames(summ) <- paste0("x", 1:ncol(x))
      }
    } else {
      rownames(summ) <- c("size", colnames(x))
    }

    colnames(summ) <- c("Estimate", "Std. Error", "z value",
                        "Pr(>|z|)")

    return(list(ESTIMATE = param,
                PREDICT = PREDICT,
                SE = SE,
                JACOBIAN = J,
                EST_FUNC = U,
                Umat = Umat,
                summary = summ))
  } else {
    return(list(ESTIMATE = param, PREDICT = PREDICT))
  }
}
