#' Single, fully parametric imputation for a right-censored covariate with conditional means following Equation (1)
#'
#' Single, fully parametric conditional mean imputation for a right-censored covariate using an accelerated failure-time model to estimate the conditional survival function and then integrates over the estimated survival function from \code{W} to \code{Inf} to compute conditional means, as in Equation (1) of the manuscript.
#'
#' @param imputation_formula imputation model formula (or coercible to formula) passed through to \code{survreg}, a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See \code{survreg} documentation for more details.
#' @param dist imputation model distribution passed through to \code{survreg}. See \code{survreg} documentation for more details.
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

cmi_fp_eq1 = function(imputation_formula, dist, W, Delta, data, maxiter = 100) {
  # Initialize imputed values
  data$imp = data[, W] ## start with imp = W

  # Fit AFT imputation model for X ~ Z
  fit = survreg(formula = imputation_formula,
                data = data,
                dist = dist,
                maxiter = maxiter)

  # Calculate linear predictor for AFT imputation model
  lp = fit$linear.predictors ## linear predictors from the survreg fit

  # Create an indicator variable for being uncensored
  uncens = data[, Delta] == 1

  # Use adaptive quadrature to estimate
  ## integral from X = 0 to X = Wi
  int_surv_W_to_Inf = sapply(
    X = which(!uncens),
    FUN = function(i) {
      integrate(f = function(x, mean, scale, distribution) {
        1 - psurvreg(q = x, mean = mean, scale = scale, distribution = distribution)
        },
                lower = data[i, W],
                upper = Inf,
                mean = lp[i],
                scale = fit$scale,
                distribution = dist)$value
    }
  )

  # Calculate S(W|Z)
  surv = 1 - psurvreg(q = data[which(!uncens), W],
                      mean = lp[which(!uncens)],
                      scale = fit$scale,
                      distribution = dist)

  # Calculate MRL(W) = int_surv / S(W|Z)
  est_mrl = int_surv_W_to_Inf / surv

  # Calculate imputed value E(X|X>W,Z) = W + MRL(W)
  data[which(!uncens), "imp"] = data[which(!uncens), W] + est_mrl

  # Return input dataset with appended column imp containing imputed values
  return(
    list(imputed_data = data,
         code = !any(is.na(data$imp)))
    )
}
