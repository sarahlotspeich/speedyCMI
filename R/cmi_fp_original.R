#' Fully parametric imputation for a right-censored covariate with conditional means using the original integral
#'
#' Fully parametric conditional mean imputation for a right-censored covariate using an accelerated failure-time model to estimate the conditional survival function and then integrates over the estimated survival function from \code{W} to \code{Inf} to compute conditional means, as in Equation (1) of the manuscript.
#'
#' @param imputation_formula imputation model formula (or coercible to formula) passed through to \code{survreg}, a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See \code{survreg} documentation for more details.
#' @param dist imputation model distribution passed through to \code{survreg}. See \code{survreg} documentation for more details.
#' @param W character, column name for observed values of the censored covariate
#' @param Delta character, column name for censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and any other variables in \code{imputation_formula}.
#' @param max_iter (optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{100}.
#' @param boots (optional) numeric, for multiple imputation supply the desired number of imputations (obtained via bootstrapping) to \code{boots}. Default is \code{0}, which is single imputation.
#' @param seed (optional) numeric, for multiple imputation set the random seed for the bootstrapping with \code{seed}. Default is \code{NULL}, which does not reset \code{seed}.
#'
#' @return A list (or, in the case of multiple imputation, a list of lists) containing:
#' \item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
#' \item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
#'
#' @export
#' @importFrom survival survreg
#' @importFrom survival Surv
#' @importFrom survival psurvreg

cmi_fp_original = function(imputation_formula, dist, W, Delta, data, maxiter = 100, boots = 0, seed = NULL) {
  if (boots == 0) {
    return_list = cmi_fp_eq1_single(imputation_formula = imputation_formula,
                                    dist = dist,
                                    W = W,
                                    Delta = Delta,
                                    data = data,
                                    maxiter = maxiter)
  } else {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    return_list = list()
    for (b in 1:boots) {
      re_data = data[sample(x = 1:nrow(data), size = nrow(data), replace = TRUE), ]
      return_list[[b]] = cmi_fp_eq1_single(imputation_formula = imputation_formula,
                                           dist = dist,
                                           W = W,
                                           Delta = Delta,
                                           data = re_data,
                                           maxiter = maxiter)
    }
  }
  return(return_list)
}
