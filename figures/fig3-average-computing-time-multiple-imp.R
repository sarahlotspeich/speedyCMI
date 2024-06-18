#############################################################################################
# Create FIGURE 3 from the manuscript, with caption beginning ###############################
## "Average computing runtime per-replication for imputation simulations..." ################
#############################################################################################

# Load packages
library(ggplot2) ## for plots

# Load simulation results from GitHub
sett = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/speedyCMI/master/sims/multiple-imputation-sims.csv")

# Create plot
sett |>
  dplyr::select(sim, num_imps, censoring, dplyr::starts_with("time")) |>
  tidyr::gather("integral", "time", -c(1:3)) |>
  dplyr::mutate(num_imps = factor(x = num_imps,
                                  levels = c(0, 5, 10, 20, 40),
                                  labels = c("B = 1", "B = 5", "B = 10", "B = 20", "B = 40")),
                integral = factor(x = integral,
                                  levels = c("time_sp", "time_old", "time_new2", "time_new", "time_analytical"),
                                  labels = c("Semiparametric", "Parametric\n(Original Integral)", "Parametric\n(Stabilized Integral\nWithout Mean)", "Parametric\n(Stabilized Integral\n With Mean)", "Parametric\n(Analytical Solution)"))
  ) |>
  dplyr::group_by(num_imps, integral) |>
  dplyr::summarize(avg_time = mean(time),
                   med_time = median(time)) |>
  ggplot(aes(x = num_imps,
             y = avg_time,
             col = integral,
             group = integral)) +
  geom_point() +
  geom_line() +
  geom_point(aes(y = med_time)) +
  geom_line(aes(y = med_time), linetype = 2) +
  theme_minimal(base_size = 14) +
  xlab("Number of Imputations") +
  ylab("Average Per-Replicate Computing Time (in Seconds)") +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold")) +
  ggthemes::scale_color_colorblind(name = "Imputation\nApproach:") +
  guides(color = guide_legend(ncol=1, byrow=TRUE))
ggsave(filename = "speedyCMI/figures/fig3-average-computing-time-multiple-imp.png",
       device = "png", width = 10, height = 6, units = "in")
