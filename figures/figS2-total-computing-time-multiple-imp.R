#############################################################################################
# Create FIGURE S2 from the manuscript, with caption beginning ##############################
## "Total computing runtime per-replication for imputation simulations ..." #################
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
                                  labels = c("B = 1 Imputation", "B = 5 Imputations", "B = 10 Imputations", "B = 20 Imputations", "B = 40 Imputations")),
                censoring = factor(x = censoring,
                                   levels = c("light", "heavy", "extra heavy"),
                                   labels = c("Light Censoring", "Heavy Censoring", "Extra Heavy Censoring")),
                integral = factor(x = integral,
                                  levels = c("time_sp", "time_old", "time_new2", "time_new", "time_analytical"),
                                  labels = c("Semiparametric", "Parametric\n(Original Integral)", "Parametric\n(Stabilized Integral\nWithout Mean)", "Parametric\n(Stabilized Integral\n With Mean)", "Parametric\n(Analytical Solution)"))
  ) |>
  dplyr::group_by(num_imps, censoring, integral) |>
  dplyr::summarize(total_time = sum(time)) |>
  ggplot(aes(x = num_imps,
             y = total_time,
             col = integral,
             group = integral)) +
  geom_point() +
  geom_line() +
  theme_bw(base_size = 14) +
  xlab("Number of Imputations") +
  ylab("Total Computing Time for 1000 Replicates (in Seconds)") +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white")) +
  ggthemes::scale_color_colorblind(name = "Imputation\nApproach:") +
  guides(color = guide_legend(ncol=1, byrow=TRUE))

ggsave(filename = "speedyCMI/figures/figS2-total-computing-time-multiple-imp.png",
       device = "png", width = 10, height = 6, units = "in")
