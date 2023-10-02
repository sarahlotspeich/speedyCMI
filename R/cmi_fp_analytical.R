#' Fully parametric conditional mean imputation for a right-censored covariate with conditional means using the analytical solution
#'
#' Fully parametric conditional mean imputation for a right-censored covariate using an accelerated failure-time model to estimate the conditional survival function and then uses an analytical solution to compute conditional means.
#'
#' @param imputation_formula imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param dist imputation model distribution passed through to \code{survreg}. See \code{survreg} documentation for more details.
#' @param W character, column name for observed values of the censored covariate.
#' @param Delta character, column name for censoring indicators. Note that values of zero in \code{Delta} are interpreted as censored observations.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and any other variables in \code{imputation_formula}.
#' @param max_iter (optional) numeric, maximum iterations allowed in call to \code{survival::survreg()}. Default is \code{max_iter = 100}.
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

cmi_fp_analytical = function(imputation_formula, dist, W, Delta, data, maxiter = 100, boots = 0, seed = NULL) {
  # Single imputation
  if (boots == 0) {
    if (toupper(dist) == "WEIBULL") {
      return_list = cmi_fp_eq5_single(imputation_formula = imputation_formula,
                                      W = W,
                                      Delta = Delta,
                                      data = data,
                                      maxiter = maxiter)
    }
  } else { # Multiple imputation
    if (!is.null(seed)) {
      set.seed(seed)
    }
    return_list = list()
    for (b in 1:boots) {
      b_rows = sample(x = 1:nrow(data),
                      size = nrow(data),
                      replace = TRUE)
      b_data = data[b_rows, ]
      if (toupper(dist) == "WEIBULL") {
        return_list[[b]] = cmi_fp_eq5_single(imputation_formula = imputation_formula,
                                             W = W,
                                             Delta = Delta,
                                             data = b_data,
                                             maxiter = maxiter)
      }
    }
  }
  return(return_list)
}
