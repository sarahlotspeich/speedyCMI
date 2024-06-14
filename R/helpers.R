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
