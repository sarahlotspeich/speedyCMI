# ///////////////////////////////////////////////////////////////////////
# Replicate Table 3 /////////////////////////////////////////////////////
# Caption begins "Estimates of $\beta_1$, the parameter on the censored /
# covariate $X$ in the linear regression analysis model, resulting from /
# the full cohort analysis (i.e., no censored values) and parametric ////
# and semiparametric multiple imputation." //////////////////////////////
# ///////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) # To wrangle data
library(tidyr) # To transform data
library(kableExtra) # To format pretty tables

# Read in simulation results
## Full cohort and semiparametric multiple imputation
oth_res = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/speedyCMI/master/sims/multiple-imputation-sims.csv") |>
  filter(num_imps == 10) |> ## subset to B = 10 imputations
  select(sim, ends_with(c("fc", "sp")), -starts_with("time")) ## keep only relevant columns
oth_res_summ = oth_res |>
  dplyr::summarize(bias_fc = mean(beta_fc - 0.5),
                   ese_fc = sd(beta_fc),
                   bias_sp = mean(beta_sp - 0.5),
                   ese_sp = sd(beta_sp),
                   ase_sp = mean(se_beta_sp),
                   cp_sp = mean((beta_sp - 1.96 * se_beta_sp) <= 0.5 &
                                  0.5 <= (beta_sp + 1.96 * se_beta_sp),
                                na.rm = TRUE)) |>
  mutate(perc_bias_fc = paste0("($", format(round(bias_fc / 0.5 * 100, 2), nsmall = 2), "$)"),
         perc_bias_sp = paste0("($", format(round(bias_sp / 0.5 * 100, 2), nsmall = 2), "$)"),
         re_sp = ese_fc ^ 2 / ese_sp ^ 2,
         mid1 = "", mid2 = "", mid3 = "")
## Parametric multiple imputation
res = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/speedyCMI/master/sims/misspecification-sims-mi.csv") |>
  select(sim, starts_with(c("beta", "se_beta"))) |>
  pivot_longer(cols = starts_with(c("beta", "se_beta")),
               names_to = "distribution",
               values_to = "value") |>
  mutate(quantity = case_when(grepl(pattern = "se", x = distribution) ~ "se_beta_analytical",
                              .default = "beta_analytical"),
         distribution = toupper(sub(pattern = ".*_",
                                    replacement = "",
                                    x = distribution)),
         distribution = factor(x = distribution,
                               levels = c("EXPO", "WEIBULL", "LOGNORM", "LOGLOG", "PWE"))) |>
  pivot_wider(names_from = quantity,
              values_from = value)
res_summ = res |>
  group_by(distribution) |>
  dplyr::summarize(bias_analytical = mean(beta_analytical - 0.5),
            ese_analytical = sd(beta_analytical),
            ase_analytical = mean(se_beta_analytical),
            cp_analytical = mean((beta_analytical - 1.96 * se_beta_analytical) <= 0.5 &
                                   0.5 <= (beta_analytical + 1.96 * se_beta_analytical),
                                 na.rm = TRUE),
  ) |>
  left_join(data.frame(distribution = unique(res$distribution),
                       ese_fc = oth_res_summ$ese_fc)) |>
  mutate(perc_bias_analytical = paste0("($", format(round(bias_analytical / 0.5 * 100, 2), nsmall = 2), "$)"),
         re_analytical = ese_fc ^ 2 / ese_analytical ^ 2,
         mid1 = "", mid2 = "", mid3 = "")

# //////////////////////////////////////////////////////////////////////
# Summarize simulation results by setting //////////////////////////////
# //////////////////////////////////////////////////////////////////////
res_summ = res_summ |>
  left_join(cbind(oth_res_summ, distribution = unique(res$distribution))) |>
  select(distribution, starts_with(c("bias", "perc_bias", "ese", "ase", "cp", "re")), everything()) |>
  select(distribution, mid1,
         ends_with("_fc"), mid2,
         ends_with("_analytical"), mid3,
         ends_with("_sp"))

# //////////////////////////////////////////////////////////////////////
# Format table for export to LaTex /////////////////////////////////////
# //////////////////////////////////////////////////////////////////////

# Write function to add "padded" zeros and wrap with $$ for consistency
format_num = function(num) {
  paste0("$", format(round(num, 3), nsmall = 3), "$")
}

# Format res_summ_wide for LaTex
res_summ |>
  mutate_at(.vars = c("bias_fc", "ese_fc",
                      "bias_analytical", "ese_analytical", "ase_analytical", "cp_analytical", "re_analytical",
                      "bias_sp", "ese_sp", "ase_sp", "cp_sp", "re_sp"),
            .funs = format_num) |>
  kable(format = "latex", booktabs = TRUE, escape = FALSE,
        align = "lrcccccccrccccc") |>
  kable_styling()
## Note: For visual reasons, the \addlinespace were manually deleted in LaTex
## And a \multicolumn used to separate the three parameters
