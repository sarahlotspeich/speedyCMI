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

![alt text](figures/fig1-average-computing-time-weibull-single-imp.png)

  - [Script (Run Simulations)](sims/single-imputation-sims.R)
  - [Script (Make Figure)](figures/fig1-average-computing-time-weibull-single-imp.R)
  - [Data (Simulation Results)](sims/single-imputation-sims.csv)

**Figure 2.** Estimates of $\beta_1$, the parameter on the censored covariate $X$ in the linear regression analysis model, resulting from each single imputation approach. The horizontal dashed line denotes the true value of $\beta_1 = 0.5$. 

![alt text](figures/fig2-betas-weibull-single-imp.png)

  - [Script (Run Simulations)](sims/single-imputation-sims.R)
  - [Script (Make Figure)](figures/fig2-betas-weibull-single-imp.R)
  - [Data (Simulation Results)](sims/single-imputation-sims.csv)

**Figure 3.** Average computing runtime per-replication for imputation simulations (in seconds) with an increasing number of imputations $B$. The solid and dashed lines connect the mean and median per-replicate computing times, respectively.

![alt text](figures/fig3-average-computing-time-weibull-multiple-imp.png)

  - [Script (Run Simulations)](sims/multiple-imputation-sims.R)
  - [Script (Make Figure)](figures/fig3-average-computing-time-weibull-multiple-imp.R)
  - [Data (Simulation Results)](sims/multiple-imputation-sims.csv)

**Figure 4.** Estimates of $\beta_1$, the parameter on the censored covariate $X$ in the linear regression analysis model, resulting from each imputation approach with an increasing number of imputations $B$. The horizontal dashed line denotes the true value of $\beta_1 = 0.5$.

![alt text](figures/fig4-betas-weibull-multiple-imp.png)

  - [Script (Run Simulations)](sims/multiple-imputation-sims.R)
  - [Script (Make Figure)](figures/fig4-betas-weibull-multiple-imp.R)
  - [Data (Simulation Results)](sims/multiple-imputation-sims.csv)
