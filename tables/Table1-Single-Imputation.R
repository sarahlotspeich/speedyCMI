# //////////////////////////////////////////////////////////////////////
# Replicate Table 1 ////////////////////////////////////////////////////
# Caption begins "Simulation results for $\pmb{\beta}$ (the coefficient/
# on censored Weibull $X$) from the full cohort analysis and ///////////
# conditional mean imputation (CMI) approaches." ///////////////////////
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) # To wrangle data
library(kableExtra) # To format pretty tables

# Read in simulation results 
res = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/speedyCMI/master/sims/single-imputation-sims.csv")

# //////////////////////////////////////////////////////////////////////
# Summarize simulation results by setting //////////////////////////////
# //////////////////////////////////////////////////////////////////////
res_summ = res |> 
  group_by(censoring, n) |> 
  summarize(bias_fc = mean(beta_fc - 0.5), 
            ese_fc = sd(beta_fc),
            bias_analytical = mean(beta_analytical - 0.5), 
            ese_analytical = sd(beta_analytical), 
            ase_analytical = mean(se_beta_analytical), 
            cp_analytical = mean((beta_analytical - 1.96 * se_beta_analytical) <= 0.5 & 
                               0.5 <= (beta_analytical + 1.96 * se_beta_analytical), 
                             na.rm = TRUE),
            bias_sp = mean(beta_sp - 0.5), 
            ese_sp = sd(beta_sp), 
            ase_sp = mean(se_beta_sp), 
            cp_sp = mean((beta_sp - 1.96 * se_beta_sp) <= 0.5 & 
                               0.5 <= (beta_sp + 1.96 * se_beta_sp), 
                             na.rm = TRUE)
            ) |> 
  mutate(perc_bias_fc = paste0("($", format(round(bias_fc / 0.5 * 100, 2), nsmall = 2), "$)"),
         perc_bias_analytical = paste0("($", format(round(bias_analytical / 0.5 * 100, 2), nsmall = 2), "$)"),
         perc_bias_sp = paste0("($", format(round(bias_sp / 0.5 * 100, 2), nsmall = 2), "$)"),
         censoring = factor(x = censoring, 
                            levels = c("light", "heavy", "extra_heavy"), 
                            labels = c("Light", "Heavy", "Extra Heavy")),
         re_analytical = ese_fc ^ 2 / ese_analytical ^ 2,
         re_sp = ese_fc ^ 2 / ese_sp ^ 2,
         mid1 = "", mid2 = "", mid3 = "") |> 
  select(censoring, n, starts_with(c("bias", "perc_bias", "ese", "ase", "cp", "re")), everything()) |> 
  arrange(censoring, n) |> 
  select(censoring, n, mid1,
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
        align = "llrcccccccrccccc") |> 
  kable_styling() 
## Note: For visual reasons, the \addlinespace were manually deleted in LaTex
## And a \multicolumn used to separate the three parameters
