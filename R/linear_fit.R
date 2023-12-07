#' Linear model
#'
#' Fit linear model E[Y | X] = offset + beta * x
#' @param y n-vector of outcome
#' @param x design matrix
#' @param offset n-vector of offset
#' @param variance whether to return the variance and components for the calculation
#' @export
linear_fit <- function(y, x, offset = rep(0, length(y)),
                        variance = TRUE) {

  x <- matrix(x, nrow = length(y))

  ## estimate the regression coefficients

  param <- c(solve(t(x) %*% x) %*% t(x) %*% (y - offset))

  ## linear predictor
  PREDICT <- c(x %*% param)
  if (variance) {
    ## Jacobian matrix of the estimating function
    J <- - t(x) %*% x
    U <- x * c(y - offset - x %*% param)
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
