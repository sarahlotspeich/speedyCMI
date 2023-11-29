Statistically and computationally efficient conditional mean imputation
for censored covariates
================

This repository contains R code and simulation data to reproduce results
from the [manuscript]() by Sarah C. Lotspeich and Ethan M. Alt (2023+).

## Package Installation

Installation of the `speedyCMI` from GitHub requires the
[`devtools`](https://www.r-project.org/nosvn/pandoc/devtools.html)
package and can be done in the following way.

``` r
# Install the package
devtools::install_github(repo = "sarahlotspeich/speedyCMI")
```

``` r
# Load the package
library(speedyCMI)
```

    ## Loading required package: survival

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

### 1. Parametric (Original Integral)

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
lm(y ~ imp + z, data = imp_dat$imputed_data)
```

    ## 
    ## Call:
    ## lm(formula = y ~ imp + z, data = imp_dat$imputed_data)
    ## 
    ## Coefficients:
    ## (Intercept)          imp            z  
    ##      1.0720       0.3867       0.2462

### 2. Parametric (Stabilized Integral Without Mean)

### 3. Parametric (Stabilized Integral With Mean)

### 4. Parametric (Analytical Solution)

### Multiple Imputation

## Figures

**Figure 1.** Average computing runtime per-replication for single
imputation simulations (in seconds). The solid and dashed lines connect
the mean and median per-replicate computing times, respectively.

<figure>
<img src="figures/fig1-average-computing-time-weibull-single-imp.png"
alt="alt text" />
<figcaption aria-hidden="true">alt text</figcaption>
</figure>

- [Script (Run Simulations)](sims/single-imputation-sims.R)
- [Script (Make
  Figure)](figures/fig1-average-computing-time-weibull-single-imp.R)
- [Data (Simulation Results)](sims/single-imputation-sims.csv)

**Figure 2.** Estimates of $\beta_1$, the parameter on the censored
covariate $X$ in the linear regression analysis model, resulting from
each single imputation approach. The horizontal dashed line denotes the
true value of $\beta_1 = 0.5$.

<figure>
<img src="figures/fig2-betas-weibull-single-imp.png" alt="alt text" />
<figcaption aria-hidden="true">alt text</figcaption>
</figure>

- [Script (Run Simulations)](sims/single-imputation-sims.R)
- [Script (Make Figure)](figures/fig2-betas-weibull-single-imp.R)
- [Data (Simulation Results)](sims/single-imputation-sims.csv)

**Figure 3.** Average computing runtime per-replication for imputation
simulations (in seconds) with an increasing number of imputations $B$.
The solid and dashed lines connect the mean and median per-replicate
computing times, respectively.

<figure>
<img src="figures/fig3-average-computing-time-weibull-multiple-imp.png"
alt="alt text" />
<figcaption aria-hidden="true">alt text</figcaption>
</figure>

- [Script (Run Simulations)](sims/multiple-imputation-sims.R)
- [Script (Make
  Figure)](figures/fig3-average-computing-time-weibull-multiple-imp.R)
- [Data (Simulation Results)](sims/multiple-imputation-sims.csv)

**Figure 4.** Estimates of $\beta_1$, the parameter on the censored
covariate $X$ in the linear regression analysis model, resulting from
each imputation approach with an increasing number of imputations $B$.
The horizontal dashed line denotes the true value of $\beta_1 = 0.5$.

<figure>
<img src="figures/fig4-betas-weibull-multiple-imp.png" alt="alt text" />
<figcaption aria-hidden="true">alt text</figcaption>
</figure>

- [Script (Run Simulations)](sims/multiple-imputation-sims.R)
- [Script (Make Figure)](figures/fig4-betas-weibull-multiple-imp.R)
- [Data (Simulation Results)](sims/multiple-imputation-sims.csv)

**Figure S1.** Total computing runtime across 1000 replicates for single
imputation simulations (in seconds).

<figure>
<img src="figures/figS1-total-computing-time-weibull-single-imp.png"
alt="alt text" />
<figcaption aria-hidden="true">alt text</figcaption>
</figure>

- [Script (Run Simulations)](sims/single-imputation-sims.R)
- [Script (Make
  Figure)](figures/figS1-total-computing-time-weibull-single-imp.R)
- [Data (Simulation Results)](sims/single-imputation-sims.csv)

**Figure S2.** Total computing runtime across 1000 replicates for
imputation simulations (in seconds) with an increasing number of
imputations $B$.

<figure>
<img src="figures/figS2-total-computing-time-weibull-multiple-imp.png"
alt="alt text" />
<figcaption aria-hidden="true">alt text</figcaption>
</figure>

- [Script (Run Simulations)](sims/multiple-imputation-sims.R)
- [Script (Make
  Figure)](figures/figS2-total-computing-time-weibull-multiple-imp.R)
- [Data (Simulation Results)](sims/multiple-imputation-sims.csv)
