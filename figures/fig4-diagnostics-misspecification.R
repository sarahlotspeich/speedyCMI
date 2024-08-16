#############################################################################################
# Create FIGURE 3 from the manuscript, with caption "Estimates of $\beta_1$, the parameter ##
# on the censored covariate $X$ in the linear regression analysis model, resulting from #####
# Parametric CMI using different distributions. The horizontal dashed line denotes the true #
# value of $\beta_1 = 0.5$. #################################################################
#############################################################################################

# Load packages
library(ggplot2) ## for plots

# Load simulation results from GitHub
sett = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/speedyCMI/master/sims/misspecification-sims-mi.csv")

# Numeric summary
aic_mat = sett |>
  dplyr::select(sim, dplyr::starts_with("aic"))
which_min_aic = apply(X = aic_mat[, -1], MARGIN = 1, FUN = which.min)
table(which_min_aic)
## AIC was lowest for log-normal (correctly specified) in 982 / 1000 replications

bic_mat = sett |>
  dplyr::select(sim, dplyr::starts_with("bic"))
which_min_bic = apply(X = bic_mat[, -1], MARGIN = 1, FUN = which.min)
table(which_min_bic)
## Same for BIC: lowest for log-normal (correctly specified) in 982 / 1000 replications

# Create plot
sett |>
  dplyr::select(sim, dplyr::starts_with(c("aic", "bic"))) |>
  tidyr::pivot_longer(cols = dplyr::starts_with(c("aic", "bic")),
                      names_to = "dist",
                      values_to = "val") |>
  dplyr::mutate(diag = toupper(sub(pattern = "_.*", replacement = "", x = dist)),
                dist = sub(pattern = "aic_",
                           replacement = "",
                           x = dist),
                dist = sub(pattern = "bic_",
                           replacement = "",
                           x = dist),
                dist = factor(x = dist,
                              levels = c("expo", "weibull", "pwe", "lognorm", "loglog"),
                              labels = c("Exponential", "Weibull", "Piecewise\nExponential", "Log-Normal", "Log-Logistic"))
                ) |>
  ggplot(aes(x = dist,
             y = val,
             col = dist,
             fill = dist)) +
  geom_boxplot(alpha = 0.6) +
  facet_wrap(~diag) +
  theme_minimal(base_size = 14) +
  xlab("Diagnostic Value") +
  ylab("Parameter Estimate") +
  theme(axis.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white")) +
  ggthemes::scale_fill_colorblind() +
  ggthemes::scale_color_colorblind() +
  guides(fill = "none",
         color = "none")
ggsave(filename = "~/Documents/speedyCMI/figures/fig4-diagnostics-misspecification.png",
       device = "png", width = 10, height = 6, units = "in")
