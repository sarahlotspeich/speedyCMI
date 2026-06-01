#' Fit analysis model to each imputed dataset and pool the results using Rubin's rules
#'
#' Fit analysis model to each imputed dataset and pool the results using Rubin's rules
#'
#' @param analysis_model imputation model formula (or coercible to formula) passed through to \code{survreg}, a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See \code{survreg} documentation for more details.
#' @param mult_imp list (or list of lists) of imputed datasets returned from one of the \code{cmi_fp} functions
#' @param rubins_rules logical, if \code{rubins_rules = TRUE} (the default) Rubin's Rules are used to calculate the pooled variance. If \code{rubins_rules = FALSE}, the approximation from Bernhardt et al. (2014) is used instead.
#'
#' @return A dataframe containing the pooled coefficients and standard errors for \code{analysis_model}

#'
#' @export
cmi_fp_pool_fit = function(analysis_model, mult_imp, rubins_rules = TRUE) {
  # Define the number of imputations
  B = length(mult_imp)

  # Fit analysis model to each imputed dataset
  mult_fit = do.call(what = rbind,
                     args = lapply(X = mult_imp,
                                   FUN = function(l) {
                                     summary(lm(formula = as.formula(analysis_model),
                                                data = l$imputed_data))$coefficients
                                   }
                     )
  )

  # Define the dimension of the coefficient vector
  p = nrow(summary(lm(formula = as.formula(analysis_model),
                      data = mult_imp[[1]]$imputed_data))$coefficients)

  # Separate coefficient and variances estimates
  beta_b = matrix(data = mult_fit[, 1],
                  nrow = B,
                  byrow = TRUE)
  vbeta_b = matrix(data = mult_fit[, 2]^2,
                   nrow = B,
                   byrow = TRUE)

  # Pool coefficient and variance estimates
  beta_pooled = colMeans(beta_b)
  beta_pooled_rep = matrix(data = beta_pooled,
                           nrow = B,
                           ncol = p,
                           byrow = TRUE)

  ## Pooled variance(s)
  vbeta_within = colMeans(vbeta_b)
  vbeta_between = colSums((beta_b - beta_pooled_rep) ^ 2) / (B - 1)
  if (rubins_rules) {
    vbeta_pooled = vbeta_within + (1 + 1 / B) * vbeta_between
  } else {
    ## Bernhardt's approximation to the correct formula (eq. 9)
    ## Fit complete-case model to get "initial" estimates and information
    cc_data = data[data[, Delta] == 1, ]
    cc_fit = lm(formula = analysis_model,
                data = cc_data)
    vbeta_initial = diag(vcov(cc_fit))
    vbeta_pooled = vbeta_within + (1 + 1 / B) * vbeta_between  +
      (vbeta_between / vbeta_within) * vbeta_initial * (vbeta_between / vbeta_within)
  }

  # Create table of estimates
  res = data.frame(Est = beta_pooled,
                   SE = sqrt(vbeta_pooled))
  res$t = res$Est / res$SE ## t-statistic
  res$p = 2 * pnorm(q = abs(res$t), mean = 0, sd = 1, lower.tail = FALSE) ## two-sided p-value
  rownames(res) = rownames(summary(lm(formula = as.formula(analysis_model),
                                      data = mult_imp[[1]]$imputed_data))$coefficients) ## take row names from one of the fits

  # Return it
  return(res)
}
