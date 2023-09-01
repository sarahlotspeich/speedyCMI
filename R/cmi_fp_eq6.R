#' Single, fully parametric conditional mean imputation for a right-censored covariate (Weibull distribution) with conditional means following Equation (6)
#'
#' Single, fully parametric conditional mean imputation for a right-censored covariate using an accelerated failure-time model with a Weibull distribution to estimate the conditional survival function and then uses an analytic solution to compute conditional means, as in Equation (6) of the manuscript.
#'
#' @param imputation_formula imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param W character, column name for observed values of the censored covariate
#' @param Delta character, column name for censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param max_iter (optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{100}.
#'
#' @return
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @export

cmi_fp_eq6 = function(imputation_formula, W, Delta, data, infinite_integral = TRUE, maxiter = 100) {
  # Initialize imputed values
  data$imp = data[, W] ## start with imp = W

  # Fit AFT imputation model for X ~ Z
  fit = survreg(formula = imputation_formula,
                data = data,
                dist = "weibull",
                maxiter = maxiter)

  # Calculate linear predictor for AFT imputation model
  lp = fit$linear.predictors ## linear predictors

  # Transform parameters to agree with R's weibull parameterization
  weib_shape = 1 / fit$scale
  weib_scale = exp(lp)

  # Create an indicator variable for being uncensored
  uncens = data[, Delta] == 1

  # Use closed-form to compute the conditional mean
  ## Transform parameters to agree with paper's parameterization
  alpha = weib_shape
  lambda = weib_scale ^ (- weib_shape)

  # Use formula from Appendix for right-censored
  ## Save quantities for use in formula
  inside_exp = lambda[which(!uncens)] * data[which(!uncens), W] ^ alpha ## inside exp() for Weibull survival function
  gamma_surv = pgamma(q = inside_exp,
                      shape = 1 / alpha,
                      scale = 1,
                      lower.tail = FALSE) ## survival function of a gamma
  data[which(!uncens), "imp"] = data[which(!uncens), W] * exp(- inside_exp) +
    gamma(1 / alpha) / (alpha * lambda[which(!uncens)] ^ (1 / alpha)) * gamma_surv ## start with numerator
  data[which(!uncens), "imp"] = data[which(!uncens), "imp"] / exp(- inside_exp) ## divide by denominator

  # Return input dataset with appended column imp containing imputed values
  return(list(imputed_data = data,
              code = !any(is.na(data$imp)))
         )
}
