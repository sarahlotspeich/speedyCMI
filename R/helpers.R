eq15_integrand = function(t, Wi, lpi, dist, fit, use_cumulative_hazard) {
  if (use_cumulative_hazard) {
    ## H{Wi + (1 - t) / t|Zi}
    num = -log(1 - psurvreg(q = Wi + (1 - t) / t,
                            mean = lpi,
                            scale = fit$scale,
                            dist = dist))
    return(exp(- (num + 2 * log(t))))
  } else {
    ## S{Wi + (1 - t) / t|Zi}
    num = 1 - psurvreg(q = Wi + (1 - t) / t,
                       mean = lpi,
                       scale = fit$scale,
                       dist = dist)
    denom = t ^ 2
    return(num / denom)
  }
}

cmi_pwe_setup = function(imputation_formula, data) {
  Terms = terms(imputation_formula, data = data)
  if ( attr(Terms, 'response') == 0 )
    stop('formula must have a Surv response')
  mf = model.frame(imputation_formula, data = data)
  Y  = model.extract(mf, 'response')
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object")
  ## Extract names of variables
  Wname = all.vars(imputation_formula)[1]
  Deltaname = all.vars(imputation_formula)[2]   ## assumes right-censoring
  ## Return list
  list('Wname' = Wname,
       'Deltaname' = Deltaname)
}

loglik_pweph <- function(y, event, X, breaks, logblhaz, beta) {

  ## Compute hazard for each interval
  lp     <- (as.matrix(X) %*% beta)[, 1]

  ## Get interval start and end points
  J     <- length(breaks) - 1
  start <- breaks[1:J]
  end   <- breaks[2:(J+1)]

  ## Initialize log likelihood contribution
  loglik <- numeric(length(y))
  for ( j in 1:J ) {
    loghaz_j <- logblhaz[j] + lp
    atrisk_j <- ( y >= start[j] )
    action_j <- atrisk_j & ( y < end[j] ) & ( event == 1 )

    ## Compute survival likelihood contribution
    logrisktime_j    <- log( pmin(y[atrisk_j], end[j]) - start[j] )
    loglik[atrisk_j] <- loglik[atrisk_j] - exp(logrisktime_j + loghaz_j[atrisk_j] )

    ## Compute hazard likelihood contribution
    loglik[action_j] <- loglik[action_j] + loghaz_j[action_j]
  }

  return(sum(loglik))
}
