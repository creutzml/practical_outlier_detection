#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### April 1st, 2023                                               ###
### Project: New Methods for Functional  Outlier Detection        ###
###                                                               ###
###   Assessing the simulation results for the procedure (which   ###
### was implemented through CU alpine).                           ###
#####################################################################

## Libraries
library(tidyverse)
library(ggplot2)
library(fdaoutlier) 
library(fda)
library(progress)
library(corrplot)
library(glmnet)
library(magrittr)
library(kableExtra)
library(ggthemes)

### First step: investigating best way to set number of intervals 
###  to use for each iteration
#####################################################################
## Load the resulting data frames in together and combine
# Directory paths
dir_path <- file.path(here::here(), "R Code")
data_path <- file.path(here::here(), "Data")

## The following block of code allows for the user to compile all new
## simulations into one data file from the saved files the simulation
## creates. E.g., on line 311 in run_outlier_sim.R, a folder is 
## created that saves the simulation results. Copy the name of the 
## folder and place below where it says *INSERT FOLDER NAME*. 

# # File path to the folder
# fold_path <- file.path(data_path, "*INSERT FOLDER NAME*")
# 
# # List the files in the folder
# fold_files <- list.files(fold_path, full.names = T)

# # Read in first data frame
# load(file = file.path(fold_files[1]))
# sim_results_all <- sim_results_sum
# 
# # Loop through all files, binding the rows to the first data frame
# p_bar <- progress::progress_bar$new(total = length(fold_files) - 1)
# for (f in 2:length(fold_files)) {
#   # Load the next data file in
#   load(file = file.path(fold_files[f]))
#   
#   # Combine into one data frame
#   sim_results_all <- dplyr::bind_rows(sim_results_all, 
#                                       sim_results_sum)
#   
#   # Remove the single file
#   rm(sim_results_sum)
#   
#   # Update progress
#   p_bar$tick()
# }

# # Save that^
# save(sim_results_all,
#      file = file.path(fold_path, "sim_results_all.RData"))

# load(file.path(fold_path, "sim_results_all.RData"))

## Load the simulation results for the type simulation
load(file.path(data_path, "sim_results_all.RData"))
#####################################################################


## Some summary info on the results
#####################################################################
# Average over everything except the methods
sim_results_summary_all <- sim_results_all %>%
  mutate(mcc = ((n_tp*n_tn - n_fp*n_fn)/
                  (sqrt((n_tp + n_fp)*(n_tp + n_fn)*
                          (n_tn + n_fp)*(n_tn + n_fn))))) %>%
  dplyr::mutate(mcc = ifelse(is.finite(mcc), mcc, NA)) %>%
  group_by(out_meth) %>%
  summarize(across(.cols = c(run_times, sensitivity:precision, mcc), 
                   list(mean = mean, sd = sd), 
                   na.rm = T, 
                   .names = "{.col}.{.fn}"), 
            n_prop_used = sum(!is.na(mcc))/n()) %>%
  dplyr::mutate(across(where(is.numeric), round, 3)) %>%
  dplyr::mutate(across(where(is.numeric), format, nsmall = 3)) %>%
  dplyr::filter(!(out_meth %in% c("easy_outliers_adj_tuk",
                                  "easy_outliers_adj_mix",
                                  "easy_outliers_c1.5_mix",
                                  "easy_outliers_c3_mix",
                                  "easy_outliers_c3_tuk",
                                  "easy_outliers_new_mix",
                                  "ff_outliers",
                                  ""))) %>%
  dplyr::mutate(out_meth = factor(
    out_meth, 
    labels = c("easy_outliers_c1.5_tuk" = "POD (Tuk)", 
               "easy_outliers_new_tuk" = "POD (user)",
               "ms_outliers" = "MS-Plot", 
               "muod_outliers_box" = "MUOD (box)",
               "muod_outliers_tan" = "MUOD (tan)",
               "tvd_outliers" = "TVD")
  ))

# Make a publishable form
sim_results_summary_all_pretty <- sim_results_summary_all %>%
  mutate(Method = out_meth,
         `Run Times (s)` = paste0(run_times.mean, " (", 
                                  run_times.sd, ") "), 
         `Sensitivity` = paste0(sensitivity.mean, " (", 
                                  sensitivity.sd, ") "), 
         `Specificity` = paste0(specificity.mean, " (", 
                                  specificity.sd, ") "), 
         `Accuracy` = paste0(accuracy.mean, " (", 
                                  accuracy.sd, ") "), 
         `Precision` = paste0(precision.mean, " (", 
                                  precision.sd, ") "), 
         `MCC` = paste0(mcc.mean, " (", 
                                  mcc.sd, ") ")) %>%
  dplyr::select(c(Method, `Run Times (s)`, Sensitivity, Specificity, 
                  Accuracy, Precision, MCC))

## Reproduces Table 1 in the manuscript in Section 3.2
# Produce latex code
kableExtra::kbl(sim_results_summary_all_pretty, 
                booktabs = TRUE, 
                format = "latex") %>%
  kableExtra::kable_styling(latex_options = c("scale_down", 
                                              "hold_position"))

# # Save that^
# write.csv(x = sim_results_summary_all,
#           file = paste0("/Users/creutzml/Library/Mobile Documents",
#                         "/com~apple~CloudDocs/Documents",
#                         "/Dissertation/functional_data_analysis",
#                         "/sim_results/outlier_methods",
#                         "/sim_results_avg_over_all_3_24_23.csv"),
#           row.names = F)

## Average with respect to n, T, and method for a plot
sim_results_summary_nT_filt <- sim_results_all %>%
  mutate(mcc = ((n_tp*n_tn - n_fp*n_fn)/
                  (sqrt((n_tp + n_fp)*(n_tp + n_fn)*
                          (n_tn + n_fp)*(n_tn + n_fn))))) %>%
  dplyr::mutate(mcc = ifelse(is.finite(mcc), mcc, NA)) %>%
  group_by(n_obs, n_smpl_pts, out_meth) %>%
  dplyr::filter(!(out_meth %in% c("easy_outliers_adj_tuk",
                                  "easy_outliers_adj_mix",
                                  "easy_outliers_c1.5_mix",
                                  "easy_outliers_c3_mix",
                                  "easy_outliers_c3_tuk",
                                  "easy_outliers_new_mix", 
                                  "", 
                                  "muod_outliers_tan",
                                  "ff_outliers"))) %>%
  summarize(across(.cols = c(run_times, 
                             sensitivity:precision, 
                             mcc), 
                   list(mean = mean), 
                   na.rm = T, 
                   .names = "{.col}.{.fn}"))

sim_results_summary_n_filt %>%
  dplyr::group_by(n_obs, out_meth) %>%
  dplyr::summarize(mcc.mean = mean(mcc.mean)) %>%
  View()
  
sim_results_summary_nT_filt_long <- sim_results_summary_nT_filt %>%
  pivot_longer(cols = c(run_times.mean:mcc.mean), 
               names_to = "Metric", 
               values_to = "Estimate") %>%
  na.omit() %>%
  dplyr::filter(Metric == "mcc.mean") %>%
  dplyr::mutate(out_meth = factor(
    out_meth, 
    labels = c("easy_outliers_c1.5_tuk" = "POD (Tuk)", 
               "easy_outliers_new_tuk" = "POD (user)",
               "ms_outliers" = "MS-Plot", 
               "muod_outliers_box" = "MUOD (box)",
               "tvd_outliers" = "TVD")
  )) %>%
  dplyr::rename("n" = "n_obs")


# Colorblind pallette
cb_pallette <- c("#999999", "#E69F00", "#56B4E9", 
                 "#009E73", "#0072B2")

## Reproduces Figure 2 in Section 3.2 of the manuscript
# Plot, fcaeted by sample size, with diagnostic trending over T
ggplot() +
  geom_line(aes(x = n_smpl_pts, 
                y = Estimate,
                color = out_meth, 
                linetype = out_meth),
            size = 1.5,
            data = sim_results_summary_nT_filt_long) + 
  geom_point(aes(x = n_smpl_pts, 
                 y = Estimate, 
                 color = out_meth, 
                 shape = out_meth), 
             size = 4,
             data = sim_results_summary_nT_filt_long) + 
  facet_grid(cols = vars(n), 
             # rows = vars(Metric),
             scales = "free", 
             labeller = label_both) +
  scale_color_manual(values = cb_pallette) +
  labs(x = "T",
       y = "Matthews Correlation Coefficient (MCC)",
       title = "Comparison of Outlier Detection Methods",
       color = "Method", 
       shape = "Method", 
       linetype = "Method") +
  theme_bw(base_size = 18)


## Average with respect to alpha, T, and method for a plot
sim_results_summary_alphaT_filt <- sim_results_all %>%
  mutate(mcc = ((n_tp*n_tn - n_fp*n_fn)/
                  (sqrt((n_tp + n_fp)*(n_tp + n_fn)*
                          (n_tn + n_fp)*(n_tn + n_fn))))) %>%
  group_by(cov_alpha, n_smpl_pts, out_meth) %>%
  dplyr::filter(!(out_meth %in% c("easy_outliers_adj_tuk",
                                  "easy_outliers_adj_mix",
                                  "easy_outliers_c1.5_mix",
                                  "easy_outliers_c3_mix",
                                  "easy_outliers_c3_tuk",
                                  "easy_outliers_new_mix", 
                                  "", 
                                  "muod_outliers_tan",
                                  "ff_outliers"))) %>%
  summarize(across(.cols = c(run_times, 
                             sensitivity:precision, 
                             mcc), 
                   list(mean = mean), 
                   na.rm = T, 
                   .names = "{.col}.{.fn}"))

# sim_results_summary_alphaT_filt %>%
#   dplyr::group_by(cov_alpha, out_meth) %>%
#   dplyr::filter(n_smpl_pts >= 100) %>%
#   dplyr::summarize(mcc.mean = mean(mcc.mean)) %>%
#   View()

sim_results_summary_alphaT_filt_long <- 
  sim_results_summary_alphaT_filt %>%
  pivot_longer(cols = c(run_times.mean:mcc.mean), 
               names_to = "Metric", 
               values_to = "Estimate") %>%
  na.omit() %>%
  dplyr::filter(Metric == "mcc.mean") %>%
  dplyr::mutate(out_meth = factor(
    out_meth, 
    labels = c("easy_outliers_c1.5_tuk" = "POD (Tuk)", 
               "easy_outliers_new_tuk" = "POD (user)",
               "ms_outliers" = "MS-Plot", 
               "muod_outliers_box" = "MUOD (box)",
               "tvd_outliers" = "TVD")
  )) %>%
  dplyr::rename("alpha" = "cov_alpha")

## Reproduces Figure 3 in Section 3.3
# Plot, fcaeted by sample size, with diagnostic trending over T
ggplot() +
  geom_line(aes(x = n_smpl_pts, 
                y = Estimate,
                color = out_meth, 
                linetype = out_meth),
            size = 1.5,
            data = sim_results_summary_alphaT_filt_long) + 
  geom_point(aes(x = n_smpl_pts, 
                 y = Estimate, 
                 color = out_meth, 
                 shape = out_meth), 
             size = 4,
             data = sim_results_summary_alphaT_filt_long) + 
  facet_grid(cols = vars(alpha), 
             # rows = vars(Metric),
             scales = "free", 
             labeller = label_both) +
  scale_color_manual(values = cb_pallette) +
  labs(x = "T",
       y = "Matthews Correlation Coefficient (MCC)",
       title = "Comparison of Outlier Detection Methods",
       color = "Method", 
       shape = "Method", 
       linetype = "Method") +
  theme_bw(base_size = 21)
  # theme(legend.position = c(0.1, 0.15), 
  #       legend.background = element_rect(fill = "white", 
  #                                        color = "black"))
# coord_cartesian(ylim = c(0.75, 1))


## Average with respect to r, T, and method for a plot
sim_results_summary_rT_filt_long <- sim_results_all %>%
  mutate(mcc = ((n_tp*n_tn - n_fp*n_fn)/
                  (sqrt((n_tp + n_fp)*(n_tp + n_fn)*
                          (n_tn + n_fp)*(n_tn + n_fn))))) %>%
  group_by(outlier_rate, n_smpl_pts, out_meth) %>%
  dplyr::filter(!(out_meth %in% c("easy_outliers_adj_tuk",
                                  "easy_outliers_adj_mix",
                                  "easy_outliers_c1.5_mix",
                                  "easy_outliers_c3_mix",
                                  "easy_outliers_c3_tuk",
                                  "easy_outliers_new_mix", 
                                  "", 
                                  "muod_outliers_tan",
                                  "ff_outliers"))) %>%
  summarize(across(.cols = c(run_times, 
                             sensitivity:precision, 
                             mcc), 
                   list(mean = mean), 
                   na.rm = T, 
                   .names = "{.col}.{.fn}")) %>%
  pivot_longer(cols = c(run_times.mean:mcc.mean), 
               names_to = "Metric", 
               values_to = "Estimate") %>%
  na.omit() %>%
  dplyr::filter(Metric == "mcc.mean") %>%
  dplyr::mutate(out_meth = factor(
    out_meth, 
    labels = c("easy_outliers_c1.5_tuk" = "POD (Tuk)", 
               "easy_outliers_new_tuk" = "POD (user)",
               "ms_outliers" = "MS-Plot", 
               "muod_outliers_box" = "MUOD (box)",
               "tvd_outliers" = "TVD")
  ))

## Produces a plot showing how the diagnostics change with respect to
## the proportion of outliers present (last statement of Section 3.2
## is a summary from this plot).
# Plot, faceted by sample size, with diagnostic trending over T
ggplot() +
  geom_line(aes(x = n_smpl_pts, 
                y = Estimate,
                color = out_meth, 
                linetype = out_meth),
            # size = 0.25,
            data = sim_results_summary_rT_filt_long) + 
  geom_point(aes(x = n_smpl_pts, 
                 y = Estimate, 
                 color = out_meth, 
                 shape = out_meth), 
             # size = 0.5,
             data = sim_results_summary_rT_filt_long) + 
  facet_grid(cols = vars(outlier_rate), 
             # rows = vars(Metric),
             scales = "free") +
  labs(x = "T",
       y = "Matthews Correlation Coefficient",
       title = "Comparison of Outlier Detection Methods",
       color = "Method", 
       shape = "Method", 
       linetype = "Method") +
  theme_bw(base_size = 14)
# coord_cartesian(ylim = c(0.75, 1))

## Average with respect to r, T, and method for a plot
sim_results_summary_modT_filt_long <- sim_results_all %>%
  mutate(mcc = ((n_tp*n_tn - n_fp*n_fn)/
                  (sqrt((n_tp + n_fp)*(n_tp + n_fn)*
                          (n_tn + n_fp)*(n_tn + n_fn))))) %>%
  group_by(sim_model, n_smpl_pts, out_meth) %>%
  dplyr::filter(!(out_meth %in% c("easy_outliers_adj_tuk",
                                  "easy_outliers_adj_mix",
                                  "easy_outliers_c1.5_mix",
                                  "easy_outliers_c3_mix",
                                  "easy_outliers_c3_tuk",
                                  "easy_outliers_new_mix", 
                                  "", 
                                  "muod_outliers_tan",
                                  "ff_outliers"))) %>%
  summarize(across(.cols = c(run_times, 
                             sensitivity:precision, 
                             mcc), 
                   list(mean = mean), 
                   na.rm = T, 
                   .names = "{.col}.{.fn}")) %>%
  pivot_longer(cols = c(run_times.mean:mcc.mean), 
               names_to = "Metric", 
               values_to = "Estimate") %>%
  na.omit() %>%
  dplyr::filter(Metric == "mcc.mean") %>%
  dplyr::mutate(out_meth = factor(
    out_meth, 
    labels = c("easy_outliers_c1.5_tuk" = "POD (Tuk)", 
               "easy_outliers_new_tuk" = "POD (user)",
               "ms_outliers" = "MS-Plot", 
               "muod_outliers_box" = "MUOD (box)",
               "tvd_outliers" = "TVD")
  ))


## Shows the results for each model used in the simulation
# Plot, fcaeted by sample size, with diagnostic trending over T
ggplot() +
  geom_line(aes(x = n_smpl_pts, 
                y = Estimate,
                color = out_meth, 
                linetype = out_meth),
            # size = 0.25,
            data = sim_results_summary_modT_filt_long) + 
  geom_point(aes(x = n_smpl_pts, 
                 y = Estimate, 
                 color = out_meth, 
                 shape = out_meth), 
             # size = 0.5,
             data = sim_results_summary_modT_filt_long) + 
  facet_grid(cols = vars(sim_model), 
             # rows = vars(Metric),
             scales = "free") +
  labs(x = "T",
       y = "Matthews Correlation Coefficient",
       title = "Comparison of Outlier Detection Methods",
       color = "Method", 
       shape = "Method", 
       linetype = "Method") +
  theme_bw(base_size = 14)
# coord_cartesian(ylim = c(0.75, 1))
#####################################################################

