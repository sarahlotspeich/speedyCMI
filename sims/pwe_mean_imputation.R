survInt_pwe <- function(centime, breaks, blhaz, X, beta) {
  X     <- as.matrix(X)
  lp    <- (X %*% beta)[, 1]
  J     <- length(breaks) - 1
  start <- breaks[1:J]
  end   <- breaks[2:(J+1)]
  diff  <- end - start
  
  fullInt  <- numeric(nrow(X))   ## \int_{0}^{\infty} S(t) dt
  lowerInt <- numeric(nrow(X))   ## \int_{0}^{c} S(t) dt
  cumHaz   <- numeric(nrow(X))
  
  for ( j in 1:J ) {
    
    ## Compute hazard for interval j
    haz_j   <- exp( log(blhaz[j]) + lp )
    
    ## Update full integral
    timeatrisk_j <- end[j] - start[j]
    fullInt      <- fullInt + ( 1 - exp( -haz_j * timeatrisk_j ) ) / haz_j * exp( -cumHaz )
    
    ## Compute lower integral
    indx_atrisk   <- centime > start[j]
    if ( sum(indx_atrisk) > 0 ) {
      timeatrisk_j  <- pmin( centime[indx_atrisk], end[j] ) - start[j]
      lowerInt[indx_atrisk] <- (
        lowerInt[indx_atrisk] 
        + ( 1 - exp( -haz_j[indx_atrisk] * timeatrisk_j ) ) / haz_j[indx_atrisk] * exp( -cumHaz[indx_atrisk] )
      )
    }
    
    ## Compute cumulative hazard for interval j
    cumHaz <- cumHaz + haz_j * (end[j] - start[j])
  }
  return(fullInt - lowerInt)
}

Spweph <- function(x, breaks, hazards) {
  J     <- length(breaks) - 1
  start <- breaks[1:J]
  end   <- breaks[2:(J+1)]
  diff  <- end - start
  
  ## Cumulative hazard computation
  n      <- length(x)
  cumHaz <- numeric(n)
  for ( j in 1:J ) {
    indx_surv           <- x >= end[j]
    cumHaz[indx_surv]   <- cumHaz[indx_surv] + hazards[j] * diff[j]
    
    indx_action         <- (x >= start[j]) & (x < end[j])
    cumHaz[indx_action] <- cumHaz[indx_action] + hazards[j] * (x[indx_action] - start[j])
  }
  exp(-cumHaz)
}

pwe_mean_imputation <- function(centime, breaks, blhaz, X, beta) {
  X     <- as.matrix(X)
  lp    <- (X %*% beta)[, 1]
  J     <- length(breaks) - 1
  start <- breaks[1:J]
  end   <- breaks[2:(J+1)]
  diff  <- end - start
  
  fullInt  <- numeric(nrow(X))   ## \int_{0}^{\infty} S(t) dt
  lowerInt <- numeric(nrow(X))   ## \int_{0}^{c} S(t) dt
  cumHaz   <- numeric(nrow(X))
  negLogSurvtime <- numeric(nrow(X))  ## -log S(c)
  
  for ( j in 1:J ) {
    
    ## Compute hazard for interval j
    haz_j   <- exp( log(blhaz[j]) + lp )
    
    ## Update full integral
    timeatrisk_j <- end[j] - start[j]
    fullInt      <- fullInt + ( 1 - exp( -haz_j * timeatrisk_j ) ) / haz_j * exp( -cumHaz )
    
    ## Compute lower integral and update hazard
    indx_atrisk   <- centime > start[j]
    if ( sum(indx_atrisk) > 0 ) {
      timeatrisk_j                <- pmin( centime[indx_atrisk], end[j] ) - start[j]
      cumhazcontrib_j             <- haz_j[indx_atrisk] * timeatrisk_j
      lowerInt[indx_atrisk] <- (
        lowerInt[indx_atrisk] 
        + ( 1 - exp( -cumhazcontrib_j ) ) / haz_j[indx_atrisk] * exp( -cumHaz[indx_atrisk] )
      )
      
      negLogSurvtime[indx_atrisk] <- negLogSurvtime[indx_atrisk] + cumhazcontrib_j
      
    }
    
    ## Compute cumulative hazard for interval j
    cumHaz <- cumHaz + haz_j * (end[j] - start[j])
  }
  Int <- fullInt - lowerInt
  
  ## Return L + [\int_{L}^{\infty} S(t) dt] / S(L)
  survTime <- exp(-negLogSurvtime)
  res      <- centime + Int / survTime
  return(res)
}