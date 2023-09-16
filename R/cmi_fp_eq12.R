#' Single, fully parametric imputation for a right-censored covariate with conditional means following Equation (12)
#'
#' Single, fully parametric conditional mean imputation for a right-censored covariate using an accelerated failure-time model to estimate the conditional survival function and then integrates over the estimated survival function from \code{0} to \code{W} to compute conditional means, as in Equation (12) of the manuscript.
#'
#' @param imputation_formula imputation model formula (or coercible to formula) passed through to \code{survreg}, a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See \code{survreg} documentation for more details.
#' @param L character, column name for the lower limit of the censoring interval for the censored covariate
#' @param U character, column name for the upper limit of the censoring interval for the censored covariate
#' @param Delta character, column name for censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param max_iter (optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{100}.
#'
#' @return
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @export

cmi_fp_eq12 = function(imputation_formula, L, U, Delta, data, maxiter = 100) {
  # Fit AFT imputation model for X ~ Z
  fit = survreg(formula = imputation_formula,
                data = data,
                dist = "lognormal",
                maxiter = maxiter)

  # Initialize imputed values
  data$imp = data[, W] ## start with imp = W

  # Calculate linear predictor for AFT imputation model
  mu = fit$linear.predictors ## linear predictors from the survreg fit = location
  scale = fit$scale ## scale from the survreg fit = scale

  # Create an indicator variable for being uncensored
  uncens = data[, Delta] == 1

  # Compute standardized U and L values (to go into standard normal CDF)
  stand_U = (log(data[, U]) - mu) / scale
  stand_L = (log(data[, L]) - mu) / scale

  # Calculate conditional means
  ## Numerator
  num = pnorm(q = (stand_U - scale), mean = 0, sd = 1, lower.tail = TRUE) -
    pnorm(q = (stand_L - scale), mean = 0, sd = 1, lower.tail = TRUE)
  ## Denominator
  denom = pnorm(q = stand_U, mean = 0, sd = 1, lower.tail = TRUE) -
    pnorm(q = stand_L, mean = 0, sd = 1, lower.tail = TRUE)
  ## Multiplier
  mult = exp(mu + (scale ^ 2) / 2)
  ## Put them together
  condl_means = mult * num / denom

  # Calculate imputed value E(X|X>W,Z) = W + MRL(W)
  data[which(!uncens), "imp"] = condl_means[which(!uncens)]

  # Return input dataset with appended column imp containing imputed values
  return_list = list(imputed_data = data,
                     code = !any(is.na(data$imp)))
  return(return_list)
}
