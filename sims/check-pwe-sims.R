# Install package
# Run once:
# devtools::install_github(repo = "sarahlotspeich/SpeedyCMI")
# devtools::install_github(repo = "sarahlotspeich/imputeCensRd")

# Load packages
library(eha) ## to simulate data (piecewise exponential)
library(imputeCensRd) ## to impute (semiparametric)
library(speedyCMI) ## to impute (parametric)
library(ggplot2) ## to plot results

source("~/Documents/speedyCMI/sims/pwe_mean_imputation.R")

# /////////////////////////////////////////////////////////////////////////
# Data generation function for all simulations ////////////////////////////
# /////////////////////////////////////////////////////////////////////////
breaks = c(1, 2, 3)
beta = -0.5
haz0 = seq(1, length(breaks) + 1) / 10
haz1 = haz0 * exp(beta)
generate_data = function(n, censoring = "light") {
  z = rbinom(n = n, size = 1, prob = 0.5) # Uncensored covariate
  x0 = rpch(n = n, cuts = breaks, levels = haz0) # To-be-censored covariate (given Z = 0)
  x1 = rpch(n = n, cuts = breaks, levels = haz1) # To-be-censored covariate (given Z = 1)
  x = (1 - z) * x0 + z * x1 # To-be-censored covariate
  e = rnorm(n = n, mean = 0, sd = 1) # Random errors
  y = 1 + 0.5 * x + 0.25 * z + e # Continuous outcome
  q = ifelse(test = censoring == "light",
             yes = 0.03, ## ~ 12%
             no = ifelse(test = censoring == "heavy",
                         yes = 0.15, ## ~ 41%
                         no = 0.7) ## ~ 78%
  ) # Rate parameter for censoring
  c = rexp(n = n, rate = q) # Random censoring mechanism
  w = pmin(x, c) # Observed covariate value
  d = as.numeric(x <= c) # "Event" indicator
  dat = data.frame(x, z, w, y, d) # Construct data set
  return(dat)
}
sims = 100 ## number of simulated replicates

# /////////////////////////////////////////////////////////////////////////
# Simulations using each CMI approach /////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////
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

  imp_dat <- cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                               dist = "pwe",
                               nintervals = 10,
                               data = dat)
  # ## Impute censored covariates
  # ## Get data-driven breaks
  # J           <- 10
  # data.breaks <- dat %>%
  #   filter(d == 1) %>%
  #   select(x) %>%
  #   unlist %>%
  #   quantile(probs = seq(1, J-1) / J )
  # data.breaks <- round( c(0, data.breaks, Inf), 3 )
  #
  # fit.mle <- pchreg(formula = Surv(w, d) ~ z,
  #                   cuts = data.breaks,
  #                   data = dat)
  #
  # mle                <- c(fit.mle$hazards, fit.mle$coefficients)
  # names(mle)[1:J]    <- paste0('blhaz', 1:J)
  # names(mle)[-(1:J)] <- paste0('beta', seq(1, length(mle) - J))
  #
  # ## Vectorized test
  # wcen <- dat$w[dat$d == 0]
  # zcen <- dat$z[dat$d == 0]
  # impX = pwe_mean_imputation(centime = wcen,
  #                            breaks = data.breaks,
  #                            blhaz = fit.mle$hazards,
  #                          X = zcen,
  #                            beta = fit.mle$coefficients)
  # imp_dat = dat
  # imp_dat$imp = imp_dat$x
  # imp_dat$imp[imp_dat$d == 0] = impX

  ## Fit model to imputed data
  fit = lm(y ~ imp + z, data = imp_dat$imputed_data)

  ## Save parameter estimates
  sett_analytical[s, c(4:6)] = fit$coefficients

  ## Save standard error estimates
  sett_analytical[s, c(8:10)] = sqrt(diag(vcov(fit)))

  ## Save progress
  sett_analytical |>
    write.csv("speedyCMI/sims/single-imputation-pwe-sims.csv",
              row.names = F)
}
