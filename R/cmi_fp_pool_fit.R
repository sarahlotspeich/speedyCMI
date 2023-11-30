#' Fit and pool models fit to multiply imputed datasets
#'
#' Fully parametric conditional mean imputation for a right-censored covariate using an accelerated failure-time model to estimate the conditional survival function and then uses an analytical solution to compute conditional means.
#'
#' @param formula analysis model formula (or coercible to formula), a formula expression as for other regression models.
#' @param family (optional) analysis model family, to be supplied to \code{glm()} to fit the model. Default is \code{family = gaussian}, but see \code{?family} for more options and details.
#' @param imp_data A list of lists returned from one of the following functions: \code{cmi_fp_original()}, \code{cmi_fp_stabilized()}, or \code{cmi_fp_analytical()}.
#'
#' @return A list (or, in the case of multiple imputation, a list of lists) containing:
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @export

cmi_fp_pool_fit = function(formula, family = gaussian, imp_data) {
  # Fit analysis model to each imputed dataset
  mult_fit = do.call(what = rbind,
                     args = lapply(X = imp_data,
                                   FUN = function(l) {
                                     summary(glm(formula = as.formula(formula),
                                                 data = l$imputed_data,
                                                 family = family))$coefficients
                                   }
                     )
  )

  # Separate coefficient and standard error estimates
  imp_params = matrix(data = mult_fit[, 1],
                      ncol = length(unique(rownames(mult_fit))),
                      byrow = TRUE)

  imp_vars = matrix(data = mult_fit[, 2] ^ 2,
                    ncol = length(unique(rownames(mult_fit))),
                    byrow = TRUE)

  # Build dataframe of pooled estimates
  return_coeff = data.frame(Estimate = colMeans(imp_params),
                            Standard.Error = sapply(X = 1:ncol(imp_vars),
                                                    FUN = function(c) sqrt(mean(imp_vars[, c]) + (B + 1) * mean((imp_params[, c] - mean(imp_params[, c])) ^ 2))))
  rownames(return_coeff) = unique(rownames(mult_fit))
  return(return_coeff)
}
