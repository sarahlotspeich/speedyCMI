# Install package
# Run once:
# devtools::install_github(repo = "sarahlotspeich/SpeedyCMI")
# devtools::install_github(repo = "sarahlotspeich/imputeCensRd")

# Load packages
library(imputeCensRd) ## to impute (semiparametric)
library(speedyCMI) ## to impute (parametric)
library(ggplot2) ## to plot results

# /////////////////////////////////////////////////////////////////////////
# Data generation function for all simulations ////////////////////////////
# /////////////////////////////////////////////////////////////////////////
generate_data = function(n, censoring = "light") {
  z = rbinom(n = n, size = 1, prob = 0.5) # Uncensored covariate
  x = rlnorm(n = n, meanlog = 0 + 0.05 * z, sdlog = 0.5) # To-be-censored covariate
  e = rnorm(n = n, mean = 0, sd = 1) # Random errors
  y = 1 + 0.5 * x + 0.25 * z + e # Continuous outcome
  q = ifelse(test = censoring == "light",
             yes = 0.2, # ~ 20%
             no = ifelse(test = censoring == "heavy",
                         yes = 0.7, # ~50%
                         no = 1.67) # ~79%
  ) # Rate parameter for censoring
  c = rexp(n = n, rate = q) # Random censoring mechanism
  w = pmin(x, c) # Observed covariate value
  d = as.numeric(x <= c) # "Event" indicator
  dat = data.frame(x, z, w, y, d) # Construct data set
  return(dat)
}
sims = 1000 ## number of simulated replicates

# /////////////////////////////////////////////////////////////////////////
# Simulations using each CMI approach /////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////
## Parametric (Original Integral)
sett_old = expand.grid(sim = 1:sims,
                       n = 1000,
                       censoring = "heavy",
                       alpha_old = NA,
                       beta_old = NA,
                       gamma_old = NA,
                       time_old = NA,
                       num_imps = c(0, 5, 10, 20, 40))
for (s in 1:nrow(sett_old)) {
  # Generate data
  set.seed(sett_old$sim[s]) ## set seed = sim ID
  dat = generate_data(n = sett_old$n[s], ## Sample size
                      censoring = sett_old$censoring[s]) ## Censoring setting

  # Number of imputations
  B = sett_old$num_imps[s]

  ## Multiply impute censored covariates
  time_imp = system.time(
    mult_imp <- cmi_fp_original(imputation_model = Surv(time = w, event = d) ~ z,
                                dist = "lognormal",
                                data = dat,
                                boots = B
    )
  )

  ## Save computing time (for the imputations)
  sett_old[s, 7] = time_imp[3]

  if(sett_old$num_imps[s] > 0) {
    ## Fit model to each imputed dataset and pool results
    fit_pooled = cmi_fp_pool_fit(analysis_model = y ~ imp + z,
                                 mu = mult_imp)

    ## Save parameter estimates
    sett_old[s, c(4:6)] = fit_pooled$Est

    ## Save standard error estimates
    sett_old[s, c(8:10)] = fit_pooled$SE
  } else {
    ## Fit model to imputed data
    fit = lm(y ~ imp + z,
             data = mult_imp$imputed_data)

    ## Save parameter estimates
    sett_old[s, c(4:6)] = fit$coefficients

    ## Save standard error estimates
    sett_old[s, c(8:10)] = sqrt(diag(vcov(fit)))
  }
}

## Parametric (Stabilized Integral With Mean)
sett_new = expand.grid(sim = 1:sims,
                       n = 1000,
                       censoring = "heavy",
                       alpha_new = NA,
                       beta_new = NA,
                       gamma_new = NA,
                       time_new = NA,
                       num_imps = c(0, 5, 10, 20, 40))
for (s in 1:nrow(sett_new)) {
  # Generate data
  set.seed(sett_new$sim[s]) ## set seed = sim ID
  dat = generate_data(n = sett_new$n[s], ## Sample size
                      censoring = sett_new$censoring[s]) ## Censoring setting

  ## Multiply impute censored covariates
  time_imp = system.time(
    mult_imp <- cmi_fp_stabilized(imputation_model = Surv(time = w, event = d) ~ z,
                                  dist = "lognormal",
                                  data = dat,
                                  with_mean = TRUE,
                                  boots = sett_new$num_imps[s]
    )
  )

  ## Save computing time (for the imputations)
  sett_new[s, 7] = time_imp[3]

  if (sett_new$num_imps[s] > 0) {
    ## Fit model to each imputed dataset and pool results
    fit_pooled = cmi_fp_pool_fit(analysis_model = y ~ imp + z,
                                 mult_imp = mult_imp)

    ## Save parameter estimates
    sett_new[s, c(4:6)] = fit_pooled$Est

    ## Save standard error estimates
    sett_new[s, c(8:10)] = fit_pooled$SE
  } else {
    ## Fit model to imputed data
    fit = lm(y ~ imp + z, data = mult_imp$imputed_data)

    ## Save parameter estimates
    sett_new[s, c(4:6)] = fit$coefficients

    ## Save standard error estimates
    sett_new[s, c(8:10)] = sqrt(diag(vcov(fit)))
  }
}

## Parametric (Stabilized Integral Without Mean)
sett_new2 = expand.grid(sim = 1:sims,
                        n = 1000,
                        censoring = "heavy",
                        alpha_new2 = NA,
                        beta_new2 = NA,
                        gamma_new2 = NA,
                        time_new2 = NA,
                        num_imps = c(0, 5, 10, 20, 40))
for (s in 1:nrow(sett_new2)) {
  # Generate data
  set.seed(sett_new2$sim[s]) ## set seed = sim ID
  dat = generate_data(n = sett_new2$n[s], ## Sample size
                      censoring = sett_new2$censoring[s]) ## Censoring setting

  ## Multiply impute censored covariates
  time_imp = system.time(
    mult_imp <- cmi_fp_stabilized(imputation_model = Surv(time = w, event = d) ~ z,
                                  dist = "lognormal",
                                  data = dat,
                                  with_mean = FALSE,
                                  use_cumulative_hazard = TRUE,
                                  boots = sett_new2$num_imps[s]
    )
  )

  ## Save computing time (for the imputations)
  sett_new2[s, 7] = time_imp[3]

  if (sett_new2$num_imps[s] > 0) {
    ## Fit model to each imputed dataset and pool results
    fit_pooled = cmi_fp_pool_fit(analysis_model = y ~ imp + z,
                                 mu = mult_imp)

    ## Save parameter estimates
    sett_new2[s, c(4:6)] = fit_pooled$Est

    ## Save standard error estimates
    sett_new2[s, c(8:10)] = fit_pooled$SE
  } else {
    ## Fit model to imputed data
    fit = lm(y ~ imp + z, data = mult_imp$imputed_data)

    ## Save parameter estimates
    sett_new2[s, c(4:6)] = fit$coefficients

    ## Save standard error estimates
    sett_new2[s, c(8:10)] = sqrt(diag(vcov(fit)))
  }
}

## Parametric (Analytical Solution)
sett_analytical = expand.grid(sim = 1:sims,
                              n = 1000,
                              censoring = "heavy",
                              alpha_analytical = NA,
                              beta_analytical = NA,
                              gamma_analytical = NA,
                              time_analytical = NA,
                              num_imps = c(0, 5, 10, 20, 40))
for (s in 1:nrow(sett_analytical)) {
  # Generate data
  set.seed(sett_analytical$sim[s]) ## set seed = sim ID
  dat = generate_data(n = sett_analytical$n[s], ## Sample size
                      censoring = sett_analytical$censoring[s]) ## Censoring setting

  ## Multiply impute censored covariates
  time_imp = system.time(
    mult_imp <- cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                                  dist = "lognormal",
                                  data = dat,
                                  boots = sett_analytical$num_imps[s]
    )
  )

  ## Save computing time (for the imputations)
  sett_analytical[s, 7] = time_imp[3]

  if (sett_analytical$num_imps[s] > 0) {
    ## Fit model to each imputed dataset and pool results
    fit_pooled = cmi_fp_pool_fit(analysis_model = y ~ imp + z,
                                 mu = mult_imp)

    ## Save parameter estimates
    sett_analytical[s, c(4:6)] = fit_pooled$Est

    ## Save standard error estimates
    sett_analytical[s, c(8:10)] = fit_pooled$SE
  } else {
    ## Fit model to imputed data
    fit = lm(y ~ imp + z, data = mult_imp$imputed_data)

    ## Save parameter estimates
    sett_analytical[s, c(4:6)] = fit$coefficients

    ## Save standard error estimates
    sett_analytical[s, c(8:10)] = sqrt(diag(vcov(fit)))
  }
}

## Semiparametric
sett_sp = expand.grid(sim = 1:sims,
                      n = 1000,
                      censoring = "heavy",
                      alpha_sp = NA,
                      beta_sp = NA,
                      gamma_sp = NA,
                      time_sp = NA,
                      num_imps = c(0, 5, 10, 20, 40))
for (s in 1:nrow(sett_sp)) {
  # Generate data
  set.seed(sett_sp$sim[s]) ## set seed = sim ID
  dat = generate_data(n = sett_sp$n[s], ## Sample size
                      censoring = sett_sp$censoring[s]) ## Censoring setting

  if (sett_sp$num_imps[s] > 0) {
    ## Multiply impute censored covariates
    time_imp = system.time(
      imp_fit <- cmi_sp_bootstrap(imputation_model = Surv(time = w, event = d) ~ z,
                                  analysis_model = y ~ imp + z,
                                  data = dat,
                                  integral = "TR",
                                  surv_between = "cf",
                                  surv_beyond = "d",
                                  B = sett_sp$num_imps[s])
    )

    ## Save computing time
    sett_sp[s, 7] = time_imp[3]

    ## Save parameter estimates
    sett_sp[s, c(4:6)] = imp_fit$Est

    ## Save standard error estimates
    sett_sp[s, c(8:10)] = imp_fit$SE
  } else {
    ## Impute censored covariates
    time_imp = system.time(
      imp_dat <- cmi_sp(imputation_model = Surv(time = w, event = d) ~ z,
                        data = dat,
                        integral = "TR",
                        surv_between = "cf",
                        surv_beyond = "d")
    )

    ## Save computing time
    sett_sp[s, 7] = time_imp[3]

    ## Fit model to imputed data
    fit = lm(y ~ imp + z, data = imp_dat$imputed_data)

    ## Save parameter estimates
    sett_sp[s, c(4:6)] = fit$coefficients

    ## Save standard error estimates
    sett_sp[s, c(8:10)] = sqrt(diag(vcov(fit)))
  }
}

## Full cohort
sett_fc = expand.grid(sim = 1:sims,
                      n = c(500, 1000, 2500, 5000),
                      censoring = c("light", "heavy"),
                      alpha_fc = NA,
                      beta_fc = NA,
                      gamma_fc = NA,
                      time_fc = NA,
                      se_alpha_fc = NA,
                      se_beta_fc = NA,
                      se_gamma_fc = NA)
for (s in 1:nrow(sett_fc)) {
  # Generate data
  set.seed(sett_fc$sim[s]) ## set seed = sim ID
  dat = generate_data(n = sett_fc$n[s], ## Sample size
                      censoring = sett_fc$censoring[s]) ## Censoring setting

  ## Fit model to full cohort data
  fit = lm(y ~ x + z, data = dat)

  ## Save parameter estimates
  sett_fc[s, c(4:6)] = fit$coefficients

  ## Save standard error estimates
  sett_fc[s, c(8:10)] = sqrt(diag(vcov(fit)))
}

# /////////////////////////////////////////////////////////////////////////
# Combine and save simulation results from all methods ////////////////////
# /////////////////////////////////////////////////////////////////////////
sett = sett_old |>
  dplyr::left_join(sett_new,
                   by = dplyr::join_by(sim == sim,
                                       n == n,
                                       censoring == censoring,
                                       num_imps == num_imps)) |>
  dplyr::left_join(sett_new2,
                   by = dplyr::join_by(sim == sim,
                                       n == n,
                                       censoring == censoring,
                                       num_imps == num_imps)) |>
  dplyr::left_join(sett_analytical,
                   by = dplyr::join_by(sim == sim,
                                       n == n,
                                       censoring == censoring,
                                       num_imps == num_imps)) |>
  dplyr::left_join(sett_sp,
                   by = dplyr::join_by(sim == sim,
                                       n == n,
                                       censoring == censoring,
                                       num_imps == num_imps)) |>
  dplyr::left_join(sett_fc,
                   by = dplyr::join_by(sim == sim,
                                       n == n,
                                       censoring == censoring))
sett |>
  write.csv("speedyCMI/sims/multiple-imputation-sims.csv",
            row.names = F)
