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
  x = rweibull(n = n, shape = 0.75, scale = 0.25 + 0.25 * z)  # To-be-censored covariate
  e = rnorm(n = n, mean = 0, sd = 1) # Random errors
  y = 1 + 0.5 * x + 0.25 * z + e # Continuous outcome
  q = ifelse(test = censoring == "light",
             yes = 0.5, ## ~ 12%
             no = ifelse(test = censoring == "heavy",
                         yes = 2.9, ## ~ 41%
                         no = 20) ## ~ 78%
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
                       n = c(500, 1000, 2500, 5000),
                       censoring = c("light", "heavy"),
                       alpha_old = NA,
                       beta_old = NA,
                       gamma_old = NA,
                       time_old = NA,
                       se_alpha_old = NA,
                       se_beta_old = NA,
                       se_gamma_old = NA
)
for (s in 1:nrow(sett_old)) {
  # Generate data
  set.seed(sett_old$sim[s]) ## set seed = sim ID
  dat = generate_data(n = sett_old$n[s], ## Sample size
                      censoring = sett_old$censoring[s]) ## Censoring setting

  ## Impute censored covariates
  time_imp = system.time(
    imp_dat <- cmi_fp_original(imputation_formula = Surv(time = w, event = d) ~ z,
                               dist = "weibull",
                               W = "w",
                               Delta = "d",
                               data = dat)
  )

  ## Save computing time
  sett_old[s, 7] = time_imp[3]

  ## Fit model to imputed data
  fit = lm(y ~ imp + z, data = imp_dat$imputed_data)

  ## Save parameter estimates
  sett_old[s, c(4:6)] = fit$coefficients

  ## Save standard error estimates
  sett_old[s, c(8:10)] = sqrt(diag(vcov(fit)))
}

## Parametric (Stabilized Integral With Mean)
sett_new = expand.grid(sim = 1:sims,
                       n = c(500, 1000, 2500, 5000),
                       censoring = c("light", "heavy"),
                       alpha_new = NA,
                       beta_new = NA,
                       gamma_new = NA,
                       time_new = NA,
                       se_alpha_new = NA,
                       se_beta_new = NA,
                       se_gamma_new = NA)
for (s in 1:nrow(sett_new)) {
  # Generate data
  set.seed(sett_new$sim[s]) ## set seed = sim ID
  dat = generate_data(n = sett_new$n[s], ## Sample size
                      censoring = sett_new$censoring[s]) ## Censoring setting

  ## Impute censored covariates
  time_imp = system.time(
    imp_dat <- cmi_fp_stabilized(imputation_formula = Surv(time = w, event = d) ~ z,
                                 dist = "weibull",
                                 W = "w",
                                 Delta = "d",
                                 data = dat,
                                 with_mean = TRUE)
  )

  ## Save computing time
  sett_new[s, 7] = time_imp[3]

  ## Fit model to imputed data
  fit = lm(y ~ imp + z, data = imp_dat$imputed_data)

  ## Save parameter estimates
  sett_new[s, c(4:6)] = fit$coefficients

  ## Save standard error estimates
  sett_new[s, c(8:10)] = sqrt(diag(vcov(fit)))
}

## Parametric (Stabilized Integral Without Mean)
sett_new2 = expand.grid(sim = 1:sims,
                        n = c(500, 1000, 2500, 5000),
                        censoring = c("light", "heavy"),
                        alpha_new2 = NA,
                        beta_new2 = NA,
                        gamma_new2 = NA,
                        time_new2 = NA,
                        se_alpha_new2 = NA,
                        se_beta_new2 = NA,
                        se_gamma_new2 = NA)
for (s in 1:nrow(sett_new2)) {
  # Generate data
  set.seed(sett_new2$sim[s]) ## set seed = sim ID
  dat = generate_data(n = sett_new2$n[s], ## Sample size
                      censoring = sett_new2$censoring[s]) ## Censoring setting

  ## Impute censored covariates
  time_imp = system.time(
    imp_dat <- cmi_fp_stabilized(imputation_formula = Surv(time = w, event = d) ~ z,
                                 dist = "weibull",
                                 W = "w",
                                 Delta = "d",
                                 data = dat,
                                 with_mean = FALSE,
                                 use_cumulative_hazard = TRUE)
  )

  ## Save computing time
  sett_new2[s, 7] = time_imp[3]

  ## Fit model to imputed data
  fit = lm(y ~ imp + z, data = imp_dat$imputed_data)

  ## Save parameter estimates
  sett_new2[s, c(4:6)] = fit$coefficients

  ## Save standard error estimates
  sett_new2[s, c(8:10)] = sqrt(diag(vcov(fit)))
}

## Parametric (Analytical Solution)
sett_analytical = expand.grid(sim = 1:sims,
                              n = c(500, 1000, 2500, 5000),
                              censoring = c("light", "heavy"),
                              alpha_analytical = NA,
                              beta_analytical = NA,
                              gamma_analytical = NA,
                              time_analytical = NA,
                              se_alpha_analytical = NA,
                              se_beta_analytical = NA,
                              se_gamma_analytical = NA)
for (s in 1:nrow(sett_analytical)) {
  # Generate data
  set.seed(sett_analytical$sim[s]) ## set seed = sim ID
  dat = generate_data(n = sett_analytical$n[s], ## Sample size
                      censoring = sett_analytical$censoring[s]) ## Censoring setting

  ## Impute censored covariates
  time_imp = system.time(
    imp_dat <- cmi_fp_analytical(imputation_formula = Surv(time = w, event = d) ~ z,
                                 dist = "weibull",
                                 W = "w",
                                 Delta = "d",
                                 data = dat)
  )

  ## Save computing time
  sett_analytical[s, 7] = time_imp[3]

  ## Fit model to imputed data
  fit = lm(y ~ imp + z, data = imp_dat$imputed_data)

  ## Save parameter estimates
  sett_analytical[s, c(4:6)] = fit$coefficients

  ## Save standard error estimates
  sett_analytical[s, c(8:10)] = sqrt(diag(vcov(fit)))
}

## Semiparametric
sett_sp = expand.grid(sim = 1:sims,
                      n = c(500, 1000, 2500, 5000),
                      censoring = c("light", "heavy"),
                      alpha_sp = NA,
                      beta_sp = NA,
                      gamma_sp = NA,
                      time_sp = NA,
                      se_alpha_sp = NA,
                      se_beta_sp = NA,
                      se_gamma_sp = NA)
for (s in 1:nrow(sett_sp)) {
  # Generate data
  set.seed(sett_sp$sim[s]) ## set seed = sim ID
  dat = generate_data(n = sett_sp$n[s], ## Sample size
                      censoring = sett_sp$censoring[s]) ## Censoring setting

  ## Impute censored covariates
  time_imp = system.time(
    imp_dat <- cmi_sp(imputation_formula = Surv(time = w, event = d) ~ z,
                      W = "w",
                      Delta = "d",
                      Z = "z",
                      data = dat,
                      trapezoidal_rule = TRUE,
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

# /////////////////////////////////////////////////////////////////////////
# Combine and save simulation results from all methods ////////////////////
# /////////////////////////////////////////////////////////////////////////
sett = sett_old |>
  dplyr::left_join(sett_new,
                   by = dplyr::join_by(sim == sim,
                                       n == n,
                                       censoring == censoring)) |>
  dplyr::left_join(sett_new2,
                   by = dplyr::join_by(sim == sim,
                                       n == n,
                                       censoring == censoring)) |>
  dplyr::left_join(sett_analytical,
                   by = dplyr::join_by(sim == sim,
                                       n == n,
                                       censoring == censoring)) |>
  dplyr::left_join(sett_sp,
                   by = dplyr::join_by(sim == sim,
                                       n == n,
                                       censoring == censoring))
sett |>
  write.csv("speedyCMI/sims/single-imputation-sims.csv",
            row.names = F)
