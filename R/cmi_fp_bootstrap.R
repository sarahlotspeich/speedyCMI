#' Multiple, parametric conditional mean imputation for a censored covariate
#'
#' Multiple, parametric conditional mean imputation for a censored covariate using parametric survival regression to estimate the conditional survival function.
#'
#' @param imputation_model Imputation model formula (or coercible to formula), a formula expression as for other regression models. The response is usually a survival object as returned by the \code{Surv} function. See the documentation for \code{Surv} for details.
#' @param dist Imputation model distribution passed through to \code{survreg}. See \code{survreg} documentation for more details.
#' @param analysis_model Analysis model formula (or coercible to formula), a formula expression as for other regression models. The response should be a continuous outcome for normal linear regression.
#' @param data Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.
#' @param type String name of the type of imputation procedure to be used. Options include \code{"original"} (the default), \code{"analytical"}, \code{"stabilized (with mean)"}, and \code{"stabilized (without mean)"}.
#' @param B numeric, number of imputations. Default is \code{10}.
#'
#' @return A dataframe containing the pooled coefficients and standard errors for \code{analysis_model}

#'
#' @export

cmi_fp_bootstrap = function(imputation_model, dist, analysis_model, data, type = "original", B = 10) {
  # Size of resample -----------------------------------------------------------
  n = nrow(data)

  # Make data a data frame (in case it's a tibble) -----------------------------
  data = data.frame(data)

  # Extract variable names from imputation_model -------------------------------
  W = all.vars(imputation_model)[1] ## censored covariate
  Delta = all.vars(imputation_model)[2] ## corresponding event indicator
  Z = all.vars(imputation_model)[-c(1:2)] ## additional covariates

  # Take names and dimension from naive fit
  data[, "imp"] = data[, W] ## start imputed value with observed
  naive_fit = lm(formula = analysis_model, data = data) ## fit naive model
  p = length(naive_fit$coefficients) ## dimension of beta vector

  # Initialize empty dataframe to hold results
  mult_fit = data.frame()

  # Loop through replicates
  for (b in 1:B) {
    # Sample with replacement from the original data ---------------------------
    re_rows = ceiling(runif(n = n, min = 0, max = 1) * nrow(data))
    re_data = data[re_rows, ]

    # Calculate proportion of censored observations in resampled data ----------
    re_prop_cens = 1 - mean(re_data[, Delta])

    # Check for censoring in resampled data
    if (0 < re_prop_cens) {
      # Impute censored x in re_data -------------------------------------------
      if (type == "analytical") {
        re_data_imp = cmi_fp_analytical(imputation_model = imputation_model,
                                        dist = dist,
                                        data = re_data,
                                        boots = 0)
      } else if (type == "stabilized (with mean)") {
        re_data_imp = cmi_fp_stabilized(imputation_model = imputation_model,
                                        dist = dist,
                                        data = re_data,
                                        with_mean = TRUE,
                                        boots = 0)
      } else if (type == "stabilized (without mean)") {
        re_data_imp = cmi_fp_stabilized(imputation_model = imputation_model,
                                        dist = dist,
                                        data = re_data,
                                        with_mean = FALSE,
                                        use_cumulative_hazard = TRUE,
                                        boots = 0)
      } else {
        re_data_imp = cmi_fp_original(imputation_model = imputation_model,
                                      dist = dist,
                                      data = re_data,
                                      boots = 0)
      }

      # After imputing, fit the analysis model ---------------------------------
      re_fit = lm(formula = analysis_model,
                  data = re_data_imp$imputed_data)
    } else {
      # If no censored, just fit the usual model -------------------------------
      re_data$imp = re_data[, W]
      re_fit = lm(formula = analysis_model,
                  data = re_data)
    }

    # Save coefficients to results matrix --------------------------------------
    mult_fit = rbind(mult_fit,
                     summary(re_fit)$coefficients)
  }

  ## Pool estimates
  pooled_est = rowsum(x = mult_fit[, 1],
                      group = rep(x = 1:p, times = B)) / B

  ## Calculate within-imputation variance
  within_var = rowsum(x = mult_fit[, 2] ^ 2,
                      group = rep(x = 1:p, times = B)) / B

  ## Calculate between-imputation variance
  between_var = rowsum(x = (mult_fit[, 1] - rep(pooled_est, times = B)) ^ 2,
                       group = rep(x = 1:p, times = B)) / (B - 1)

  ## Pool variance
  pooled_var = within_var + between_var + (between_var / B)

  # Return table of pooled estimates
  tab = data.frame(Coefficient = names(naive_fit$coefficients),
                   Est = pooled_est,
                   SE = sqrt(pooled_var))
  return(tab)
}
