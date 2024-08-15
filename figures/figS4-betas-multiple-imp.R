#############################################################################################
# Create FIGURE S4 from the manuscript, with caption "Estimates of $\beta_1$, the parameter #
# on the censored covariate $X$ in the linear regression analysis model, resulting from #####
# each imputation approach with an increasing number of imputations $B$. The horizontal #####
# dashed line denotes the true value of $\beta_1 = 0.5$. ####################################
#############################################################################################

# Load packages
library(ggplot2) ## for plots

# Load simulation results from GitHub
sett = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/speedyCMI/master/sims/multiple-imputation-sims.csv")

# Create plot
sett |>
  dplyr::select(sim, num_imps, censoring, dplyr::starts_with("beta")) |>
  tidyr::gather("integral", "beta", -c(1:3)) |>
  dplyr::mutate(num_imps = factor(x = num_imps,
                                  levels = c(0, 5, 10, 20, 40),
                                  labels = c("B = 1", "B = 5", "B = 10", "B = 20", "B = 40")),
                integral = factor(x = integral,
                                  levels = c("beta_sp", "beta_old", "beta_new2", "beta_new", "beta_analytical"),
                                  labels = c("Semiparametric", "Parametric\n(Original Integral)", "Parametric\n(Stabilized Integral\nWithout Mean)", "Parametric\n(Stabilized Integral\n With Mean)", "Parametric\n(Analytical Solution)"))
  ) |>
  ggplot(aes(x = num_imps,
             y = beta,
             col = integral,
             fill = integral)) +
  geom_boxplot(alpha = 0.6) +
  geom_hline(yintercept = 0.5,
             linetype = 2) +
  theme_minimal(base_size = 14) +
  xlab("Number of Imputations") +
  ylab("Parameter Estimate") +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold")) +
  ggthemes::scale_fill_colorblind(name = "Imputation\nApproach:") +
  ggthemes::scale_color_colorblind(name = "Imputation\nApproach:") +
  guides(fill = guide_legend(ncol=1, byrow=TRUE),
         color = guide_legend(ncol=1, byrow=TRUE))
ggsave(filename = "~/Documents/speedyCMI/figures/figS4-betas-multiple-imp.png",
       device = "png", width = 10, height = 6, units = "in")
