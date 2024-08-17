#############################################################################################
# Create FIGURE 3 from the manuscript, with caption "Estimates of $\beta_1$, the parameter ##
# on the censored covariate $X$ in the linear regression analysis model, resulting from #####
# Parametric CMI using different distributions. The horizontal dashed line denotes the true #
# value of $\beta_1 = 0.5$. #################################################################
#############################################################################################

# Load packages
library(ggplot2) ## for plots

# Load simulation results from GitHub
## Multiple imputation
mi_sett = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/speedyCMI/master/sims/misspecification-sims-mi.csv")
## Single imputation
si_sett = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/speedyCMI/master/sims/misspecification-sims-si.csv")
## Combine them
sett = mi_sett |>
  dplyr::mutate(Imputation = "Multiple Imputation") |>
  dplyr::bind_rows(
    si_sett |>
      dplyr::mutate(Imputation = "Single Imputation")
  )

# Create plot
sett |>
  dplyr::select(sim, Imputation, dplyr::starts_with("beta")) |>
  tidyr::pivot_longer(cols = dplyr::starts_with("beta"),
                      names_to = "dist",
                      values_to = "beta") |>
  dplyr::mutate(dist = sub(pattern = "beta_",
                           replacement = "",
                           x = dist),
                dist = factor(x = dist,
                              levels = c("expo", "weibull", "pwe", "lognorm", "loglog"),
                              labels = c("Exponential", "Weibull", "Piecewise\nExponential", "Log-Normal", "Log-Logistic"))
                ) |>
  ggplot(aes(x = dist,
             y = beta,
             col = dist,
             fill = dist)) +
  geom_boxplot(alpha = 0.6) +
  geom_hline(yintercept = 0.5,
             linetype = 2) +
  theme_minimal(base_size = 14) +
  xlab("Distribution for the Imputation Model") +
  ylab("Parameter Estimate") +
  facet_grid(cols = vars(Imputation)) +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white")) +
  ggthemes::scale_fill_colorblind() +
  ggthemes::scale_color_colorblind() +
  guides(fill = "none",
         color = "none")
ggsave(filename = "speedyCMI/figures/fig3-betas-misspecification.png",
       device = "png", width = 10, height = 6, units = "in")
