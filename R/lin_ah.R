################################################################################
#=========== semiparametric additive hazards model =============================
################################################################################
# Semiparametric additive hazards model in Lin & Ying 1997 with time-independent
# covariates
# time: vector of time-to-event
# event: vector of event indicator; 0 = censored
# covariates: design matrix of covariates, intercept should not be included.

lin_ah <- function(time, event, covariates, weights = NULL,
                   offset = rep(0, length(time)), variance = TRUE) {
  t <- as.numeric(time)
  d <- as.numeric(event)
  x <- as.matrix(covariates)
  eta <- as.numeric(offset)
  
  nn <- length(t) # sample size
  
  # assign weights
  if (is.null(weights)) weights <- rep(1, nn)
  
  w <- as.numeric(weights)
  
  # sort the observations by time in increasing order
  dat1 <- cbind(t, d, w, eta, x)[order(t),]
  
  t1 <- dat1[, 1]
  d1 <- dat1[, 2]
  w1 <- dat1[, 3]
  eta1 <- dat1[, 4]
  x1 <- matrix(dat1[, -(1:4)], nrow = nrow(dat1))
  
  t2 <- unique(t1)
  ntime <- length(t2) # length of unique time points
  
  # obtain the order of time 
  o1 <- left_join(data.frame(t1 = t1), 
                  data.frame(t1 = t2, o1 = 1:length(unique(t1))),
                  by = "t1")$o1
  
  # the minimum index of those who have t1 equals the ith ordered t1
  tmin <- sapply(1:ntime, function(tt) min(which(o1 == tt)))
  dtime <- c(t2[1], diff(t2))
  Zb_Tk <- matrix(NA, nrow = ntime, ncol = ncol(x1))
  for (k in 1:ntime) {
    Zb_Tk[k, ] <- colSums(matrix(x1[tmin[k]:nn,], ncol = ncol(x1))) /
      sum(o1 >= k)
  }
  
  Zb_Ti <- Zb_Tk[o1, ]
  
  # integral of Zb_Tk over time
  int_Zb_Tk <- apply(Zb_Tk * dtime, 2, cumsum)
  int_Zb_Ti <- int_Zb_Tk[o1, ]
  ## calculative the Jacobian matrix

  
  A1 <- t(x1) %*% (x1 * t1 * w1)
  
  
  A2 <- t(x1) %*% (int_Zb_Ti * w1)
  
  A3 <- colSums((x1 - Zb_Ti) * d1 * w1 - x1 * eta1 * w1 * t1 + int_Zb_Ti * eta1 * w1)
  
  ESTIMATE <- as.numeric(solve(A1 - A2) %*% A3)
  
  PREDICT <- as.numeric(x1 %*% ESTIMATE)
  ## obtain the estimating equation
  # U01 <- (x1 - Zb_Ti) * d1 
  # 
  # U02 <- x1 * c(x1 %*% ESTIMATE) * t1
  # 
  # U03 <- A2_i * c(x1 %*% ESTIMATE)
  # 
  # 
  # U0 <- (U01 - U02 + U03) * w1
  
  ## Jacobian matrix of the estimating equation
  if (variance == TRUE) {
    J <- (- A1 + A2) 
    
    
    
    # baseline cumulative hazard function
    haz <- 1:ntime * NA
    
    for (k in 1:ntime) {
      haz[k] <- sum(d1 * (o1 == k) - (eta1 + x1 %*% ESTIMATE) * (o1 >= k) * dtime[k]) /
        sum(o1 >= k)
    }
    cumhaz <- cumsum(haz)
    cumhaz_i <- cumhaz[o1]
    
    # Estimating functions
    U1 <- x1 * c(d1 - cumhaz_i - c(eta1 + x1 %*% ESTIMATE) * t1)
    U2 <- Zb_Ti * d1
    U3_k <- apply(Zb_Tk * haz, 2, cumsum)
    U3 <- U3_k[o1,]
    U4_k <- apply(Zb_Tk * dtime, 2, cumsum)
    U4 <- U4_k[o1,] * c(eta1 + x1 %*% ESTIMATE)
    
    U <- U1 - U2 + U3 + U4  
    Umat <- t(U) %*% (U * w1)

    SE <- sqrt(diag(solve(J) %*% Umat %*% t(solve(J))))
    return(list(ESTIMATE = ESTIMATE, ## point estimate
                PREDICT = PREDICT, ## linear predictor
                SE = SE, ## model based standard error
                EST_FUNC = U, ## estimating functions
                Umat = Umat,  ## (weighted) "meat" matrix 
                HAZ = haz, ## hazard function at each interval
                CUMHAZ_K = cumhaz, ## cumulative hazard function at each event time
                JACOBIAN = J))  ## Jacobian matrix; "bread"
  } else {
    return(list(ESTIMATE = ESTIMATE, PREDICT = PREDICT))
  }
}
