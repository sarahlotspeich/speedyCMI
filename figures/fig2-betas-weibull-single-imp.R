#############################################################################################
# Create FIGURE 2 from the manuscript, with caption "Estimates of $\beta_1$, the parameter on
## the censored covariate $X$ in the linear regression analysis model, resulting from each
## single imputation approach. The horizontal dashed line denotes the true value of $\beta_1 = 0.5$.
#############################################################################################

# Load packages
library(ggplot2) ## for plots

# Load simulation results from GitHub
sett = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/speedyCMI/master/sims/single-imputation-sims.csv")

# Create plot
sett |>
  dplyr::select(sim, n, censoring, dplyr::starts_with("beta")) |>
  tidyr::gather("integral", "beta", -c(1:3)) |>
  dplyr::mutate(n = factor(x = n,
                           levels = c(500, 1000, 2500, 5000),
                           labels = c("n = 500", "n = 1000", "n = 2500", "n = 5000")),
                censoring = factor(x = censoring,
                                   levels = c("light", "heavy", "extra heavy"),
                                   labels = c("Light Censoring", "Heavy Censoring", "Extra Heavy Censoring")),
                integral = factor(x = integral,
                                  levels = c("beta_sp", "beta_old", "beta_new2", "beta_new", "beta_analytical"),
                                  labels = c("Semiparametric", "Parametric\n(Original Integral)", "Parametric\n(Stabilized Integral\nWithout Mean)", "Parametric\n(Stabilized Integral\n With Mean)", "Parametric\n(Analytical Solution)"))
  ) |>
  ggplot(aes(x = n,
             y = beta,
             col = integral,
             fill = integral)) +
  geom_boxplot(alpha = 0.6) +
  geom_hline(yintercept = 0.5,
             linetype = 2) +
  theme_minimal(base_size = 14) +
  facet_grid(cols = vars(censoring)) +
  xlab("Sample Size") +
  ylab("Parameter Estimate") +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white"),
        axis.title = element_text(face = "bold")) +
  ggthemes::scale_fill_colorblind(name = "Imputation\nApproach:") +
  ggthemes::scale_color_colorblind(name = "Imputation\nApproach:") +
  guides(fill = guide_legend(ncol=1, byrow=TRUE),
         color = guide_legend(ncol=1, byrow=TRUE))
ggsave(filename = "speedyCMI/figures/fig2-betas-weibull-single-imp.png",
       device = "png", width = 10, height = 6, units = "in")
