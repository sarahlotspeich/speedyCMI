#' Single, fully parametric imputation for an interval-censored covariate with conditional means following Equation (2)
#'
#' Single, fully parametric conditional mean imputation for an interval-censored covariate using an accelerated failure-time model to estimate the conditional survival function and then integrates over the estimated survival function from \code{W} to \code{Inf} to compute conditional means, as in Equation (1) of the manuscript.
#'
#' @param imputation_formula imputation model formula (or coercible to formula) passed through to \code{survreg}, a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See \code{survreg} documentation for more details.
#' @param dist imputation model distribution passed through to \code{survreg}. See \code{survreg} documentation for more details.
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

cmi_fp_eq2 = function(imputation_formula, dist, L, U, Delta, data, maxiter = 100) {
  # Initialize imputed values
  data$imp = data[, L] ## start with imp = L

  # Fit AFT imputation model for X ~ Z
  fit = survreg(formula = imputation_formula,
                data = data,
                dist = dist,
                maxiter = maxiter)

  # Calculate linear predictor for AFT imputation model
  lp = fit$linear.predictors ## linear predictors from the survreg fit

  # Create an indicator variable for being uncensored
  uncens = data[, Delta] == 1

  # Use adaptive quadrature to estimate integral from X = Li to X = Ui
  int_surv_L_to_U = sapply(
    X = which(!uncens),
    FUN = function(i) {
      integrate(f = function(x, mean, scale, distribution) {
        1 - psurvreg(q = x, mean = mean, scale = scale, distribution = distribution)
        },
                lower = data[i, L],
                upper = data[i, W],
                mean = lp[i],
                scale = fit$scale,
                distribution = dist)$value
    }
  )

  # Calculate S(L|Z)
  surv_L = 1 - psurvreg(q = data[which(!uncens), L],
                        mean = lp[which(!uncens)],
                        scale = fit$scale,
                        distribution = dist)

  # Calculate S(U|Z)
  surv_U = 1 - psurvreg(q = data[which(!uncens), U],
                        mean = lp[which(!uncens)],
                        scale = fit$scale,
                        distribution = dist)

  # Calculate imputed value E(X|L < X <= U, Z) = imp_num / imp_denom
  ## Numerator: L x S(L|Z) - U x S(U|Z) + integral of S(x|Z) from x = L to x = U
  imp_num = data[which(!uncens), L] * surv_L + data[which(!uncens), U] * surv_U + int_surv_L_to_U

  ## Denominator: S(L|Z) - S(U|Z)
  imp_denom = (surv_L - surv_U)

  ## Put them together
  data[which(!uncens), "imp"] = imp_num / imp_denom

  # Return input dataset with appended column imp containing imputed values
  return_list = list(imputed_data = data,
                     code = !any(is.na(data$imp)))
  return(return_list)
}
