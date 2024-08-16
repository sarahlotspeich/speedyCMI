# Install package
# Run once:
# devtools::install_github(repo = "sarahlotspeich/SpeedyCMI")

# Load packages
library(speedyCMI) ## to impute (parametric)

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
## Log-Normal
sett_lognorm = expand.grid(sim = 1:sims,
                           alpha_lognorm = NA,
                           beta_lognorm = NA,
                           gamma_lognorm = NA,
                           se_alpha_lognorm = NA,
                           se_beta_lognorm = NA,
                           se_gamma_lognorm = NA,
                           aic_lognorm = NA,
                           bic_lognorm = NA)
for (s in 1:nrow(sett_lognorm)) {
  # Generate data
  set.seed(sett_lognorm$sim[s]) ## set seed = sim ID
  dat = generate_data(n = 1000) ## Sample size

  # Singly impute censored covariates
  imp_dat = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                              dist = "lognormal",
                              data = dat)

  # Fit model to imputed data
  fit = lm(y ~ imp + z, data = imp_dat$imputed_data)

  # Save parameter estimates
  sett_lognorm[s, c(2:4)] = fit$coefficients

  # Save standard error estimates
  sett_lognorm[s, c(5:7)] = sqrt(diag(vcov(fit)))

  # Save model diagnostics AIC / BIC
  sett_lognorm[s, c(8:9)] = with(imp_dat, c(aic, bic))
}

## Log-Normal
sett_loglog = expand.grid(sim = 1:sims,
                          alpha_loglog = NA,
                          beta_loglog = NA,
                          gamma_loglog = NA,
                          se_alpha_loglog = NA,
                          se_beta_loglog = NA,
                          se_gamma_loglog = NA,
                          aic_loglog = NA,
                          bic_loglog = NA)
for (s in 1:nrow(sett_loglog)) {
  # Generate data
  set.seed(sett_loglog$sim[s]) ## set seed = sim ID
  dat = generate_data(n = 1000) ## Sample size

  # Singly impute censored covariates
  imp_dat = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                              dist = "loglogistic",
                              data = dat)

  # Fit model to imputed data
  fit = lm(y ~ imp + z, data = imp_dat$imputed_data)

  # Save parameter estimates
  sett_loglog[s, c(2:4)] = fit$coefficients

  # Save standard error estimates
  sett_loglog[s, c(5:7)] = sqrt(diag(vcov(fit)))

  # Save model diagnostics AIC / BIC
  sett_loglog[s, c(8:9)] = with(imp_dat, c(aic, bic))
}

## Exponential
sett_expo = expand.grid(sim = 1:sims,
                           alpha_expo = NA,
                           beta_expo = NA,
                           gamma_expo = NA,
                           se_alpha_expo = NA,
                           se_beta_expo = NA,
                           se_gamma_expo = NA,
                           aic_expo = NA,
                           bic_expo = NA)
for (s in 1:nrow(sett_expo)) {
  # Generate data
  set.seed(sett_expo$sim[s]) ## set seed = sim ID
  dat = generate_data(n = 1000) ## Sample size

  # Singly impute censored covariates
  imp_dat = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                              dist = "exponential",
                              data = dat)

  # Fit model to imputed data
  fit = lm(y ~ imp + z, data = imp_dat$imputed_data)

  # Save parameter estimates
  sett_expo[s, c(2:4)] = fit$coefficients

  # Save standard error estimates
  sett_expo[s, c(5:7)] = sqrt(diag(vcov(fit)))

  # Save model diagnostics AIC / BIC
  sett_expo[s, c(8:9)] = with(imp_dat, c(aic, bic))
}

## Weibull
sett_weibull = expand.grid(sim = 1:sims,
                           alpha_weibull = NA,
                           beta_weibull = NA,
                           gamma_weibull = NA,
                           se_alpha_weibull = NA,
                           se_beta_weibull = NA,
                           se_gamma_weibull = NA,
                           aic_weibull = NA,
                           bic_weibull = NA)
for (s in 1:nrow(sett_weibull)) {
  # Generate data
  set.seed(sett_weibull$sim[s]) ## set seed = sim ID
  dat = generate_data(n = 1000) ## Sample size

  # Singly impute censored covariates
  imp_dat = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                              dist = "weibull",
                              data = dat)

  # Fit model to imputed data
  fit = lm(y ~ imp + z, data = imp_dat$imputed_data)

  # Save parameter estimates
  sett_weibull[s, c(2:4)] = fit$coefficients

  # Save standard error estimates
  sett_weibull[s, c(5:7)] = sqrt(diag(vcov(fit)))

  # Save model diagnostics AIC / BIC
  sett_weibull[s, c(8:9)] = with(imp_dat, c(aic, bic))
}

## Piecewise exponential
sett_pwe = expand.grid(sim = 1:sims,
                           alpha_pwe = NA,
                           beta_pwe = NA,
                           gamma_pwe = NA,
                           se_alpha_pwe = NA,
                           se_beta_pwe = NA,
                           se_gamma_pwe = NA,
                           aic_pwe = NA,
                           bic_pwe = NA)
for (s in 1:nrow(sett_pwe)) {
  # Generate data
  set.seed(sett_pwe$sim[s]) ## set seed = sim ID
  dat = generate_data(n = 1000) ## Sample size

  # Singly impute censored covariates
  imp_dat = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                              dist = "pwe",
                              nintervals = 10,
                              data = dat)

  # Fit model to imputed data
  fit = lm(y ~ imp + z, data = imp_dat$imputed_data)

  # Save parameter estimates
  sett_pwe[s, c(2:4)] = fit$coefficients

  # Save standard error estimates
  sett_pwe[s, c(5:7)] = sqrt(diag(vcov(fit)))

  # Save model diagnostics AIC / BIC
  sett_pwe[s, c(8:9)] = with(imp_dat, c(aic, bic))
}

# /////////////////////////////////////////////////////////////////////////
# Combine and save simulation results from all distributions //////////////
# /////////////////////////////////////////////////////////////////////////
sett = sett_lognorm |>
  dplyr::left_join(sett_loglog,
                   by = dplyr::join_by(sim == sim)) |>
  dplyr::left_join(sett_expo,
                   by = dplyr::join_by(sim == sim)) |>
  dplyr::left_join(sett_weibull,
                   by = dplyr::join_by(sim == sim)) |>
  dplyr::left_join(sett_pwe,
                   by = dplyr::join_by(sim == sim))

sett |>
  write.csv("speedyCMI/sims/misspecification-sims-si.csv",
            row.names = F)
