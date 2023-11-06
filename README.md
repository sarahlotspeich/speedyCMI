# Statistically and computationally efficient conditional mean imputation for censored covariates

This repository contains R code and simulation data to reproduce results from the [manuscript]() by Sarah C. Lotspeich and Ethan M. Alt (2023+). 

## Package Installation

Installation of the `speedyCMI` from GitHub requires the
[`devtools`](https://www.r-project.org/nosvn/pandoc/devtools.html)
package and can be done in the following way.

``` r
# Install the package
devtools::install_github(repo = "sarahlotspeich/speedyCMI", 
                         ref = "main")
```

## Figures 

**Figure 1.** Average computing runtime per-replication for single imputation simulations (in seconds). The solid and dashed lines connect the mean and median per-replicate computing times, respectively.

  - [Script (Run Simulations)](sims/)
  - [Script (Make Figure)](figures/fig1-average-computing-time-weibull-single-imp.R)
  - [Data (Simulation Results)](sims/si_single_imp.csv)
