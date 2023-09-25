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
