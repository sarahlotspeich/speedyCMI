Statistically and computationally efficient conditional mean imputation
for censored covariates
================

This repository contains R code and simulation data to reproduce results
from the [manuscript]() by Sarah C. Lotspeich and Ethan M. Alt (2023+).

## Package Installation

Installation of the `speedyCMI` package from GitHub requires the
[`devtools`](https://www.r-project.org/nosvn/pandoc/devtools.html)
package and can be done in the following way.

``` r
# Install the package
devtools::install_github(repo = "sarahlotspeich/speedyCMI", 
                         ref = "main")
```

Once installed, the package can then be loaded in the usual way.

``` r
# Load the package
library(speedyCMI)
```

## Functions

The `speedyCMI` package contains separate functions to run conditional
mean imputation (CMI) for censored covariates using the different
formulas discussed in the corresponding paper. Using the following
simulated dataset, we provide a brief overview of these functions and
their use below.

``` r
set.seed(918) # For reproducibility
n = 1000 # Sample size
z = rbinom(n = n, size = 1, prob = 0.5) # Uncensored covariate
x = rweibull(n = n, shape = 0.75, scale = 0.25 + 0.25 * z)  # To-be-censored covariate
e = rnorm(n = n, mean = 0, sd = 1) # Random errors
y = 1 + 0.5 * x + 0.25 * z + e # Continuous outcome
q = 0.5 # Rate parameter for censoring
c = rexp(n = n, rate = q) # Random censoring mechanism
w = pmin(x, c) # Observed covariate value
d = as.numeric(x <= c) # "Event" indicator
dat = data.frame(x, z, w, y, d) # Construct data set
head(dat) # Inspect data set
```

    ##             x z           w         y d
    ## 1 0.539817268 1 0.539817268 1.5004468 1
    ## 2 0.712698933 1 0.712698933 2.7059625 1
    ## 3 0.138822892 1 0.138822892 1.4781689 1
    ## 4 3.174653770 0 2.664720071 2.9043014 0
    ## 5 0.009582614 1 0.009582614 1.8282313 1
    ## 6 0.095243557 0 0.095243557 0.9791424 1

The functions are first introduced for single imputation, before
providing an example with multiple imputation.

### 1. Parametric CMI (Original Integral)

``` r
# Single parametric CMI using the original integral formula
imp_dat = cmi_fp_original(imputation_formula = Surv(time = w, event = d) ~ z,
                          dist = "weibull",
                          W = "w",
                          Delta = "d",
                          data = dat)

## Inspect the return object (list)
head(imp_dat$imputed_data) ### completed dataset with new column "imp" 
```

    ##             x z           w         y d         imp
    ## 1 0.539817268 1 0.539817268 1.5004468 1 0.539817268
    ## 2 0.712698933 1 0.712698933 2.7059625 1 0.712698933
    ## 3 0.138822892 1 0.138822892 1.4781689 1 0.138822892
    ## 4 3.174653770 0 2.664720071 2.9043014 0 3.249785438
    ## 5 0.009582614 1 0.009582614 1.8282313 1 0.009582614
    ## 6 0.095243557 0 0.095243557 0.9791424 1 0.095243557

``` r
imp_dat$code ### logical indicator for whether or not imputation model converged
```

    ## [1] TRUE

``` r
## Fit the analysis model to singly imputed data
lm(y ~ imp + z, 
   data = imp_dat$imputed_data)
```

    ## 
    ## Call:
    ## lm(formula = y ~ imp + z, data = imp_dat$imputed_data)
    ## 
    ## Coefficients:
    ## (Intercept)          imp            z  
    ##      1.0720       0.3867       0.2462

### 2. Parametric CMI (Stabilized Integral Without Mean)

``` r
# Single parametric CMI using the stabilized integral formula without mean
imp_dat = cmi_fp_stabilized(imputation_formula = Surv(time = w, event = d) ~ z,
                            dist = "weibull",
                            data = dat, 
                            with_mean = FALSE)

## Fit the analysis model to singly imputed data
lm(y ~ imp + z, 
   data = imp_dat$imputed_data)
```

    ## 
    ## Call:
    ## lm(formula = y ~ imp + z, data = imp_dat$imputed_data)
    ## 
    ## Coefficients:
    ## (Intercept)          imp            z  
    ##      1.0720       0.3867       0.2462

### 3. Parametric CMI (Stabilized Integral With Mean)

``` r
# Single parametric CMI using the stabilized integral formula with mean
imp_dat = cmi_fp_stabilized(imputation_formula = Surv(time = w, event = d) ~ z,
                            dist = "weibull",
                            data = dat, 
                            with_mean = TRUE)

## Fit the analysis model to singly imputed data
lm(y ~ imp + z, 
   data = imp_dat$imputed_data)
```

    ## 
    ## Call:
    ## lm(formula = y ~ imp + z, data = imp_dat$imputed_data)
    ## 
    ## Coefficients:
    ## (Intercept)          imp            z  
    ##      1.0720       0.3867       0.2462

### 4. Parametric CMI (Analytical Solution)

``` r
# Single parametric CMI using the analytical solution
imp_dat = cmi_fp_analytical(imputation_formula = Surv(time = w, event = d) ~ z,
                            dist = "weibull",
                            data = dat)

## Fit the analysis model to singly imputed data
lm(y ~ imp + z, 
   data = imp_dat$imputed_data)
```

    ## 
    ## Call:
    ## lm(formula = y ~ imp + z, data = imp_dat$imputed_data)
    ## 
    ## Coefficients:
    ## (Intercept)          imp            z  
    ##      1.0720       0.3867       0.2462

### 6. Parametric CMI (Multiple Imputation)

Using the analytical solution as our example, the following code snippet
outlines how the use of the `cmi_fp_analytical()` function can be
modified for multiple imputation via bootstrap resampling. Notice that
the first modification is in supplying the desired number of imputations
to the `boots` argument.

``` r
# Multiple parametric CMI using the analytical solution
B = 20 ## Desired number of imputations 
mult_imp = cmi_fp_analytical(imputation_formula = Surv(time = w, event = d) ~ z, 
                             dist = "weibull",
                             data = dat, 
                             boots = B) 

## Inspect the return object (list of lists)
head(mult_imp[[1]]$imputed_data) ### completed dataset with new column "imp" from imputation 1
```

    ##             x z         w          y d       imp
    ## 605 0.3963907 0 0.3963907 0.06337772 1 0.3963907
    ## 670 0.2185118 0 0.2185118 0.11257288 1 0.2185118
    ## 592 1.1698449 1 0.1521263 0.61717024 0 0.7488695
    ## 10  0.4303014 1 0.4303014 1.35695814 1 0.4303014
    ## 626 1.6291951 1 1.6291951 1.77973530 1 1.6291951
    ## 641 0.1534294 0 0.1534294 0.87258973 1 0.1534294

``` r
mult_imp[[1]]$code ### logical indicator for whether or not imputation model converged from imputation 1
```

    ## [1] TRUE

``` r
## Fit the analysis model to multiplt imputed data
cmi_fp_pool_fit(formula = y ~ imp + z, 
                imp_data = mult_imp)
```

    ##              Estimate Standard.Error
    ## (Intercept) 1.0835455      0.2254262
    ## imp         0.3923314      0.2706760
    ## z           0.2194103      0.3151813

## Figures

**Figure 1.** Average computing runtime per replicate for single imputation simulations (in seconds). The solid and dashed lines connect the mean and median per-replicate computing times, respectively, across 1000 replicates.

- [Script (Run Simulations)](sims/single-imputation-sims.R)
- [Script (Make Figure)](figures/fig1-average-computing-time-weibull-single-imp.R)
- [Data (Simulation Results)](sims/single-imputation-sims.csv)

**Figure 2.** Average computing runtime per-replicate for imputation simulations (in seconds) with an increasing number of imputations $B$. The solid and dashed lines connect the mean and median per-replicate computing times, respectively. Censoring was heavy, and $n = 1000$ subjects were simulated per replicate.

- [Script (Run Simulations)](sims/multiple-imputation-sims.R)
- [Script (Make Figure)](figures/fig2-average-computing-time-multiple-imp.R)
- [Data (Simulation Results)](sims/multiple-imputation-sims.csv)

**Figure 3.** Estimates of $\beta_1$, the parameter on the censored covariate $X$ in the linear regression analysis model, resulting from the multiple ($B = 10$) and single parametric CMI approaches assuming different distributions for the imputation model of $X$ given $Z$. The true distribution for $X$ given $Z$ was log-normal, and the horizontal dashed line denotes the true value of $\beta_1 = 0.5$. 

- [Script (Run Simulations)](sims/misspecification-sims-mi.R)
- [Script (Make Figure)](figures/fig3-betas-misspecification.R)
- [Data (Simulation Results)](sims/misspecification-sims-mi.csv)

**Figure 4.** Model diagnostics Akaike information criterion (AIC) and Bayesian information criterion (BIC) for the parametric imputation model of $X$ given $Z$. The true distribution was log-normal.

- [Script (Run Simulations)](sims/misspecification-sims-mi.R)
- [Script (Make Figure)](figures/fig4-diagnostics-misspecification.R)
- [Data (Simulation Results)](sims/misspecification-sims-mi.csv)

**Figure 5.** Model diagnostics Akaike information criterion (AIC) and Bayesian information criterion (BIC) for the parametric imputation model of $X$ given $Z$. The true distribution was log-normal.

- [Script (Fit Models and Make Figure)](framingham/framingham.Rmd)
- [Data (Framingham Teaching data)](https://rdrr.io/cran/riskCommunicator/man/framingham.html)

**Figure S1.** Total computing runtime across 1000 replicates for single imputation simulations (in seconds).

- [Script (Run Simulations)](sims/single-imputation-sims.R)
- [Script (Make Figure)](figures/figS1-total-computing-time-single-imp.R)
- [Data (Simulation Results)](sims/single-imputation-sims.csv)

**Figure S2.** Estimates of $\beta_1$, the parameter on the censored covariate $X$ in the linear regression analysis model, resulting from each single imputation approach. The horizontal dashed line denotes the true value of $\beta_1 = 0.5$.

- [Script (Run Simulations)](sims/single-imputation-sims.R)
- [Script (Make Figure)](figures/figS2-betas-single-imp.R)
- [Data (Simulation Results)](sims/single-imputation-sims.csv)

**Figure S3.** Total computing runtime across 1000 replication for imputation simulations (in seconds) with an increasing number of imputations $B$. Censoring was heavy, and $n = 1000$ subjects were simulated per replicate.

- [Script (Run Simulations)](sims/multiple-imputation-sims.R)
- [Script (Make Figure)](figures/figS3-total-computing-time-multiple-imp.R)
- [Data (Simulation Results)](sims/multiple-imputation-sims.csv)

**Figure S4.** Estimates of $\beta_1$, the parameter on the censored covariate $X$ in the linear regression analysis model, resulting from each imputation approach with an increasing number of imputations $B$. The horizontal dashed line denotes the true value of $\beta_1 = 0.5$. Censoring was heavy, and $n = 1000$ subjects were simulated per replicate.

- [Script (Run Simulations)](sims/multiple-imputation-sims.R)
- [Script (Make Figure)](figures/figS4-betas-multiple-imp.R)
- [Data (Simulation Results)](sims/multiple-imputation-sims.csv)

**Figure S5.** Empirical density of time from first visit to hypertension diagnosis (observed or singly imputed) in the Framingham teaching dataset using various distributions for the imputation model. The vertical dashed line denotes $TIME = 25$ years to hypertension diagnosis, which was the end of follow-up.

- [Script (Fit Models and Make Figure)](framingham/framingham.Rmd)
- [Data (Framingham Teaching data)](https://rdrr.io/cran/riskCommunicator/man/framingham.html)
