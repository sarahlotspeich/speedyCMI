#############################################################################################
# Create FIGURE 1 from the manuscript, with caption beginning ###############################
## "Average computing runtime per-replication for single imputation simulations..." #########
#############################################################################################

# Load packages
library(ggplot2) ## for plots

# Load simulation results from GitHub
sett = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/speedyCMI/master/sims/single-imputation-sims.csv")

# Create plot
sett |>
  dplyr::select(sim, n, censoring, dplyr::starts_with("time")) |>
  tidyr::gather("integral", "time", -c(1:3)) |>
  dplyr::mutate(n = factor(x = n,
                           levels = c(500, 1000, 2500, 5000),
                           labels = c("n = 500", "n = 1000", "n = 2500", "n = 5000")),
                censoring = factor(x = censoring,
                                   levels = c("light", "heavy", "extra heavy"),
                                   labels = c("Light Censoring", "Heavy Censoring", "Extra Heavy Censoring")),
                integral = factor(x = integral,
                                  levels = c("time_sp", "time_old", "time_new2", "time_new", "time_analytical"),
                                  labels = c("Semiparametric", "Parametric\n(Original Integral)", "Parametric\n(Stabilized Integral\nWithout Mean)", "Parametric\n(Stabilized Integral\n With Mean)", "Parametric\n(Analytical Solution)"))
  ) |>
  dplyr::group_by(n, censoring, integral) |>
  dplyr::summarize(avg_time = mean(time),
                   med_time = median(time)) |>
  ggplot(aes(x = n,
             y = avg_time,
             col = integral,
             group = integral)) +
  geom_point() +
  geom_line() +
  geom_point(aes(y = med_time)) +
  geom_line(aes(y = med_time), linetype = 2) +
  theme_minimal(base_size = 14) +
  facet_grid(cols = vars(censoring)) +
  xlab(label = "Sample Size") +
  ylab("Average Per-Replicate Computing Time (in Seconds)") +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white"),
        axis.title = element_text(face = "bold")) +
  ggthemes::scale_color_colorblind(name = "Imputation\nApproach:") +
  guides(color = guide_legend(ncol=1, byrow=TRUE))
ggsave(filename = "speedyCMI/figures/fig1-average-computing-time-single-imp.png",
       device = "png", width = 10, height = 6, units = "in")
