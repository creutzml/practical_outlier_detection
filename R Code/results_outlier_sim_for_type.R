#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### April 10th, 2023                                              ###
### Project: New Methods for Functional  Outlier Detection        ###
###                                                               ###
###   Assessing the simulation results for the procedure (which   ###
### was implemented through CU alpine). How well does each method ###
### predict the type of functional outlier?                       ###
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



### First step: combine all results into one file
#####################################################################
## Load the resulting data frames in together and combine
# Directory paths
dir_path <- file.path(here::here(), "R Code")
data_path <- file.path(here::here(), "Data")

## The following block of code allows for the user to compile all new
## simulations into one data file from the saved files the simulation
## creates. E.g., on line 517, a folder is created that saves the
## simulation results. Copy the name of the folder and place below 
## where it says *INSERT FOLDER NAME*. 

# # File path to the folder
# fold_path <- file.path(data_path, "*INSERT FOLDER NAME*")
# 
# # List the files in the folder
# fold_files <- list.files(fold_path, full.names = T)
# 
# # Read in first data frame
# load(file = file.path(fold_files[1]))
# sim_results_all <- sim_results_sum
# 
# # Loop through all files, binding the rows to the first data frame
# p_bar <- progress::progress_bar$new(total = length(fold_files))
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

## Load the simulation results for the type simulation
load(file.path(data_path, "sim_results_for_type.RData"))
#####################################################################



### Some summary info on the results
### Averaged over method
#####################################################################
## Create the MCC metric, then get some summary data frames:
# Summarized everything per method
sim_results_summary_meth <- sim_results_all %>%
  mutate(mcc = ((n_tp*n_tn - n_fp*n_fn)/
                  (sqrt((n_tp + n_fp)*(n_tp + n_fn)*
                          (n_tn + n_fp)*(n_tn + n_fn)))), 
         mcc_shift = ((n_tp_shift*n_tn_shift - 
                         n_fp_shift*n_fn_shift)/
                        (sqrt((n_tp_shift + n_fp_shift)*
                                (n_tp_shift + n_fn_shift)*
                                (n_tn_shift + n_fp_shift)*
                                (n_tn_shift + n_fn_shift)))), 
         mcc_shape = ((n_tp_shape*n_tn_shape - 
                         n_fp_shape*n_fn_shape)/
                        (sqrt((n_tp_shape + n_fp_shape)*
                                (n_tp_shape + n_fn_shape)*
                                (n_tn_shape + n_fp_shape)*
                                (n_tn_shape + n_fn_shape))))) %>%
  dplyr::mutate(mcc = ifelse(is.finite(mcc), mcc, NA), 
                mcc_shift = ifelse(is.finite(mcc_shift), 
                                   mcc_shift, NA), 
                mcc_shape = ifelse(is.finite(mcc_shape), 
                                   mcc_shape, NA)) %>%
  group_by(out_meth) %>%
  summarize(across(.cols = c(sensitivity:mcc_shape), 
                   list(mean = mean, sd = sd), 
                   na.rm = T, 
                   .names = "{.col}.{.fn}")) 

# Split based on shift vs shape
sim_results_summary_meth_shift <- sim_results_summary_meth %>%
  dplyr::select(out_meth, 
                sensitivity_shift.mean:precision_shift.sd,
                mcc_shift.mean, mcc_shift.sd) %>%
  dplyr::mutate(across(where(is.numeric), round, 3)) %>%
  dplyr::mutate(across(where(is.numeric), format, nsmall = 3)) %>%
  dplyr::mutate(Method = dplyr::case_when(
    out_meth == "pod_outliers_c1.5" ~ "POD (Tuk)", 
    out_meth == "pod_outliers_user" ~ "POD (user)", 
    out_meth == "muod_outliers_box" ~ "MUOD (box)",
    out_meth == "muod_outliers_tan" ~ "MUOD (tan)",
    out_meth == "tvd_outliers" ~ "TVD"),
    `Sensitivity` = paste0(sensitivity_shift.mean, " (", 
                           sensitivity_shift.sd, ") "), 
    `Specificity` = paste0(specificity_shift.mean, " (", 
                           specificity_shift.sd, ") "), 
    `Accuracy` = paste0(accuracy_shift.mean, " (", 
                        accuracy_shift.sd, ") "), 
    `Precision` = paste0(precision_shift.mean, " (", 
                         precision_shift.sd, ") "), 
    `MCC` = paste0(mcc_shift.mean, " (", 
                   mcc_shift.sd, ") ")) %>%
  dplyr::select(c(Method, Sensitivity, Specificity, 
                  Accuracy, Precision, MCC))

## Reproduce "Magnitude" table of Table 2 in manuscript
# Produce latex code
kableExtra::kbl(sim_results_summary_meth_shift, 
                booktabs = TRUE, 
                format = "latex") %>%
  kableExtra::kable_styling(latex_options = c("scale_down", 
                                              "hold_position"))


# Split based on shift vs shape
sim_results_summary_meth_shape <- sim_results_summary_meth %>%
  dplyr::select(out_meth, 
                sensitivity_shape.mean:precision_shape.sd,
                mcc_shape.mean, mcc_shape.sd) %>%
  dplyr::mutate(across(where(is.numeric), round, 3)) %>%
  dplyr::mutate(across(where(is.numeric), format, nsmall = 3)) %>%
  dplyr::mutate(Method = dplyr::case_when(
    out_meth == "pod_outliers_c1.5" ~ "POD (Tuk)", 
    out_meth == "pod_outliers_user" ~ "POD (user)", 
    out_meth == "muod_outliers_box" ~ "MUOD (box)",
    out_meth == "muod_outliers_tan" ~ "MUOD (tan)",
    out_meth == "tvd_outliers" ~ "TVD"),
    `Sensitivity` = paste0(sensitivity_shape.mean, " (", 
                           sensitivity_shape.sd, ") "), 
    `Specificity` = paste0(specificity_shape.mean, " (", 
                           specificity_shape.sd, ") "), 
    `Accuracy` = paste0(accuracy_shape.mean, " (", 
                        accuracy_shape.sd, ") "), 
    `Precision` = paste0(precision_shape.mean, " (", 
                         precision_shape.sd, ") "), 
    `MCC` = paste0(mcc_shape.mean, " (", 
                   mcc_shape.sd, ") ")) %>%
  dplyr::select(c(Method, Sensitivity, Specificity, 
                  Accuracy, Precision, MCC))

## Reproduce "Shape" table of Table 2 in manuscript
# Produce latex code
kableExtra::kbl(sim_results_summary_meth_shape, 
                booktabs = TRUE, 
                format = "latex") %>%
  kableExtra::kable_styling(latex_options = c("scale_down", 
                                              "hold_position"))
#####################################################################



### Averaged over n, T, and method:
#####################################################################
## Create the MCC metric, then get some summary data frames:
# Summarized over n, T, and method
sim_results_summary_nTmeth <- sim_results_all %>%
  mutate(mcc = ((n_tp*n_tn - n_fp*n_fn)/
                  (sqrt((n_tp + n_fp)*(n_tp + n_fn)*
                          (n_tn + n_fp)*(n_tn + n_fn)))), 
         mcc_shift = ((n_tp_shift*n_tn_shift - 
                         n_fp_shift*n_fn_shift)/
                        (sqrt((n_tp_shift + n_fp_shift)*
                                (n_tp_shift + n_fn_shift)*
                                (n_tn_shift + n_fp_shift)*
                                (n_tn_shift + n_fn_shift)))), 
         mcc_shape = ((n_tp_shape*n_tn_shape - 
                         n_fp_shape*n_fn_shape)/
                        (sqrt((n_tp_shape + n_fp_shape)*
                                (n_tp_shape + n_fn_shape)*
                                (n_tn_shape + n_fp_shape)*
                                (n_tn_shape + n_fn_shape))))) %>%
  dplyr::mutate(mcc = ifelse(is.finite(mcc), mcc, NA), 
                mcc_shift = ifelse(is.finite(mcc_shift), 
                                   mcc_shift, NA), 
                mcc_shape = ifelse(is.finite(mcc_shape), 
                                   mcc_shape, NA)) %>%
  group_by(n_obs, n_smpl_pts, out_meth) %>%
  summarize(across(.cols = c(sensitivity:mcc_shape), 
                   list(mean = mean), 
                   na.rm = T, 
                   .names = "{.col}.{.fn}")) 


# Summarized over n, T, and method; pivoted for just shape outliers
sim_results_summary_nTmeth_shape <- sim_results_summary_nTmeth %>%
  dplyr::select(-c(sensitivity.mean:precision_shift.mean, 
                   mcc.mean, mcc_shift.mean)) %>%
  pivot_longer(cols = c(sensitivity_shape.mean:precision_shape.mean,
                        mcc_shape.mean), 
               names_to = "Metric", 
               values_to = "Estimate") %>%
  dplyr::mutate(Method = dplyr::case_when(
    out_meth == "pod_outliers_c1.5" ~ "POD (Tuk)", 
    out_meth == "pod_outliers_user" ~ "POD (user)", 
    out_meth == "muod_outliers_box" ~ "MUOD (box)",
    out_meth == "muod_outliers_tan" ~ "MUOD (tan)",
    out_meth == "tvd_outliers" ~ "TVD")) %>%
  dplyr::filter(Metric == "mcc_shape.mean") %>%
  dplyr::rename("n" = "n_obs") %>%
  dplyr::mutate(Out_Type = "Shape")

# sim_results_summary_nTmeth_shape %>%
#   group_by(n, n_smpl_pts, Method) %>%
#   dplyr::filter(Metric == "mcc_shape.mean") %>%
#   dplyr::summarise(avg_mcc = mean(Estimate, na.rm = T)) %>%
#   View()

# change the order of factor levels for method to match rest of paper
sim_results_summary_nTmeth_shape$Method <- factor(
  sim_results_summary_nTmeth_shape$Method,
  levels = c("POD (Tuk)", 
             "POD (user)", 
             "MUOD (box)",
             "MUOD (tan)",
             "TVD")
)

cb_pallette <- c("#999999", "#E69F00", "#56B4E9", 
                 "#009E73", "#0072B2")

# ## Reproduces the "Compare Classification of Shape Outliers" plot
# ## in Figure 4 of the manuscript
# # Plot, fcaeted by sample size, with diagnostic trending over T
# shape_plot <- ggplot() +
#   geom_line(aes(x = n_smpl_pts, 
#                 y = Estimate,
#                 color = Method, 
#                 linetype = Method),
#             linewidth = 1.5,
#             data = sim_results_summary_nTmeth_shape) + 
#   geom_point(aes(x = n_smpl_pts, 
#                  y = Estimate, 
#                  color = Method, 
#                  shape = Method), 
#              size = 4,
#              data = sim_results_summary_nTmeth_shape) + 
#   facet_grid(cols = vars(n),
#              scales = "free",
#              labeller = purrr::partial(label_both, sep = " = ")) +
#   scale_color_manual(values = cb_pallette) +
#   labs(y = "Matthews Correlation Coefficient (MCC)",
#        x = "T", 
#        title = "Comparing Classification of Shape Outliers") +
#   # theme(strip.text = element_text(size = 20))
#   theme_bw(base_size = 18)


# Summarized over n, T, and method; pivoted for just shift outliers
sim_results_summary_nTmeth_shift <- sim_results_summary_nTmeth %>%
  dplyr::select(-c(sensitivity.mean:precision.mean,
                   sensitivity_shape.mean:precision_shape.mean,
                   mcc.mean, mcc_shape.mean)) %>%
  pivot_longer(cols = c(sensitivity_shift.mean:precision_shift.mean,
                        mcc_shift.mean), 
               names_to = "Metric", 
               values_to = "Estimate") %>%
  dplyr::mutate(Method = dplyr::case_when(
    out_meth == "pod_outliers_c1.5" ~ "POD (Tuk)", 
    out_meth == "pod_outliers_user" ~ "POD (user)", 
    out_meth == "muod_outliers_box" ~ "MUOD (box)",
    out_meth == "muod_outliers_tan" ~ "MUOD (tan)",
    out_meth == "tvd_outliers" ~ "TVD")) %>%
  dplyr::filter(Metric == "mcc_shift.mean") %>%
  dplyr::rename("n" = "n_obs") %>%
  dplyr::mutate(Out_Type = "Magnitude")

# sim_results_summary_nTmeth_shift %>%
#   group_by(n, n_smpl_pts, Method) %>%
#   dplyr::filter(Metric == "mcc_shift.mean") %>%
#   dplyr::summarise(avg_mcc = mean(Estimate, na.rm = T)) %>%
#   View()

sim_results_summary_nTmeth_shift$Method <- factor(
  sim_results_summary_nTmeth_shift$Method,
  levels = c("POD (Tuk)", 
             "POD (user)", 
             "MUOD (box)",
             "MUOD (tan)",
             "TVD")
)

# ## Reproduces the "Compare Classification of Magnitude Outliers" plot
# ## in Figure 4 of the manuscript
# # Plot, fcaeted by sample size, with diagnostic trending over T
# mag_plot <- ggplot() +
#   geom_line(aes(x = n_smpl_pts, 
#                 y = Estimate,
#                 color = Method, 
#                 linetype = Method),
#             size = 1.5,
#             data = sim_results_summary_nTmeth_shift) + 
#   geom_point(aes(x = n_smpl_pts, 
#                  y = Estimate, 
#                  color = Method, 
#                  shape = Method), 
#              size = 4,
#              data = sim_results_summary_nTmeth_shift) + 
#   facet_grid(cols = vars(n), 
#              # rows = vars(Metric),
#              scales = "free",
#              labeller = purrr::partial(label_both, sep = " = ")) +
#   scale_color_manual(values = cb_pallette) +
#   labs(y = "Matthews Correlation Coefficient (MCC)",
#        x = "T", 
#        title = "Comparing Classification of Magnitude Outliers") +
#   theme_bw(base_size = 18)


## Make it as a single plot instead
sim_results_summary_nTmeth_all <- dplyr::bind_rows(
  sim_results_summary_nTmeth_shift, sim_results_summary_nTmeth_shape
)


## Make Figure 4
ggplot() +
  geom_line(aes(x = n_smpl_pts, 
                y = Estimate,
                color = Method, 
                linetype = Method),
            size = 1,
            data = sim_results_summary_nTmeth_all) + 
  geom_point(aes(x = n_smpl_pts, 
                 y = Estimate, 
                 color = Method, 
                 shape = Method), 
             size = 2.5,
             data = sim_results_summary_nTmeth_all) + 
  facet_grid(cols = vars(n), 
             rows = vars(Out_Type),
             scales = "free",
             labeller = label_bquote(
               cols = "n = "*.(n), 
               rows = .(Out_Type)*" Outlier")) +
  scale_color_manual(values = cb_pallette) +
  labs(y = "Matthews Correlation Coefficient (MCC)",
       x = "T") +
  theme_bw(base_size = 16) +
  theme(panel.spacing = unit(1, "lines"), 
        strip.text = element_text(size = 16), 
        legend.position = "bottom", 
        legend.key.size = unit(1, "cm"))
#####################################################################




## The follow block of code is used to assess the affect of \beta on 
## the performance of each method in the simulation (summary of the
## results is the last paragraph of Section 3.3)
### Averaged over beta, T, and method:
#####################################################################
## Create the MCC metric, then get some summary data frames:
# Summarized over beta, T, and method
sim_results_summary_betaTmeth <- sim_results_all %>%
  mutate(mcc = ((n_tp*n_tn - n_fp*n_fn)/
                  (sqrt((n_tp + n_fp)*(n_tp + n_fn)*
                          (n_tn + n_fp)*(n_tn + n_fn)))), 
         mcc_shift = ((n_tp_shift*n_tn_shift - 
                         n_fp_shift*n_fn_shift)/
                        (sqrt((n_tp_shift + n_fp_shift)*
                                (n_tp_shift + n_fn_shift)*
                                (n_tn_shift + n_fp_shift)*
                                (n_tn_shift + n_fn_shift)))), 
         mcc_shape = ((n_tp_shape*n_tn_shape - 
                         n_fp_shape*n_fn_shape)/
                        (sqrt((n_tp_shape + n_fp_shape)*
                                (n_tp_shape + n_fn_shape)*
                                (n_tn_shape + n_fp_shape)*
                                (n_tn_shape + n_fn_shape))))) %>%
  group_by(cov_beta, n_smpl_pts, out_meth) %>%
  summarize(across(.cols = c(sensitivity:mcc_shape), 
                   list(mean = mean), 
                   na.rm = T, 
                   .names = "{.col}.{.fn}")) 


# Summarized over beta, T, and method; pivoted for just shape outliers
sim_results_summary_betaTmeth_shape <- sim_results_summary_betaTmeth %>%
  dplyr::select(-c(sensitivity.mean:precision_shift.mean, 
                   mcc.mean, mcc_shift.mean)) %>%
  pivot_longer(cols = c(sensitivity_shape.mean:precision_shape.mean,
                        mcc_shape.mean), 
               names_to = "Metric", 
               values_to = "Estimate") %>%
  dplyr::mutate(Method = dplyr::case_when(
    out_meth == "pod_outliers_c1.5" ~ "POD (Tuk)", 
    out_meth == "pod_outliers_user" ~ "POD (user)", 
    out_meth == "muod_outliers_box" ~ "MUOD (box)",
    out_meth == "muod_outliers_tan" ~ "MUOD (tan)",
    out_meth == "tvd_outliers" ~ "TVD")) %>%
  dplyr::filter(Metric == "mcc_shape.mean") %>%
  dplyr::rename("beta" = "cov_beta")

# Plot, fcaeted by beta, with diagnostic trending over T
ggplot() +
  geom_line(aes(x = n_smpl_pts, 
                y = Estimate,
                color = Method, 
                linetype = Method),
            size = 0.5,
            data = sim_results_summary_betaTmeth_shape) + 
  geom_point(aes(x = n_smpl_pts, 
                 y = Estimate, 
                 color = Method, 
                 shape = Method), 
             size = 0.5,
             data = sim_results_summary_betaTmeth_shape) + 
  facet_grid(cols = vars(beta), 
             # rows = vars(Metric),
             scales = "free",
             labeller = purrr::partial(label_both, sep = " = ")) +
  labs(y = "Matthews Correlation Coefficient (MCC)",
       x = "T", 
       title = "Comparing Classification of Shape Outliers") +
  theme_bw(base_size = 18)


# Summarized over beta, T, and method; pivoted for just shift outliers
sim_results_summary_betaTmeth_shift <- sim_results_summary_betaTmeth %>%
  dplyr::select(-c(sensitivity.mean:precision.mean,
                   sensitivity_shape.mean:precision_shape.mean,
                   mcc.mean, mcc_shape.mean)) %>%
  pivot_longer(cols = c(sensitivity_shift.mean:precision_shift.mean,
                        mcc_shift.mean), 
               names_to = "Metric", 
               values_to = "Estimate") %>%
  dplyr::mutate(Method = dplyr::case_when(
    out_meth == "pod_outliers_c1.5" ~ "POD (Tuk)", 
    out_meth == "pod_outliers_user" ~ "POD (user)", 
    out_meth == "muod_outliers_box" ~ "MUOD (box)",
    out_meth == "muod_outliers_tan" ~ "MUOD (tan)",
    out_meth == "tvd_outliers" ~ "TVD")) %>%
  dplyr::filter(Metric == "mcc_shift.mean") %>%
  dplyr::rename("beta" = "cov_beta")

# Plot, fcaeted by beta, with diagnostic trending over T
ggplot() +
  geom_line(aes(x = n_smpl_pts, 
                y = Estimate,
                color = Method, 
                linetype = Method),
            size = 0.5,
            data = sim_results_summary_betaTmeth_shift) + 
  geom_point(aes(x = n_smpl_pts, 
                 y = Estimate, 
                 color = Method, 
                 shape = Method), 
             size = 0.5,
             data = sim_results_summary_betaTmeth_shift) + 
  facet_grid(cols = vars(beta), 
             # rows = vars(Metric),
             scales = "free",
             labeller = purrr::partial(label_both, sep = " = ")) +
  labs(y = "Matthews Correlation Coefficient (MCC)",
       x = "T", 
       title = "Comparing Classification of Magnitude Outliers") +
  theme_bw(base_size = 18)
#####################################################################



### Averaged over out rate, T, and method:
#####################################################################
## Create the MCC metric, then get some summary data frames:
# Summarized over r, T, and method
sim_results_summary_rTmeth <- sim_results_all %>%
  mutate(mcc = ((n_tp*n_tn - n_fp*n_fn)/
                  (sqrt((n_tp + n_fp)*(n_tp + n_fn)*
                          (n_tn + n_fp)*(n_tn + n_fn)))), 
         mcc_shift = ((n_tp_shift*n_tn_shift - 
                         n_fp_shift*n_fn_shift)/
                        (sqrt((n_tp_shift + n_fp_shift)*
                                (n_tp_shift + n_fn_shift)*
                                (n_tn_shift + n_fp_shift)*
                                (n_tn_shift + n_fn_shift)))), 
         mcc_shape = ((n_tp_shape*n_tn_shape - 
                         n_fp_shape*n_fn_shape)/
                        (sqrt((n_tp_shape + n_fp_shape)*
                                (n_tp_shape + n_fn_shape)*
                                (n_tn_shape + n_fp_shape)*
                                (n_tn_shape + n_fn_shape))))) %>%
  group_by(outlier_rate, n_smpl_pts, out_meth) %>%
  dplyr::mutate(mcc = ifelse(is.finite(mcc), mcc, NA), 
                mcc_shift = ifelse(is.finite(mcc_shift), 
                                   mcc_shift, NA), 
                mcc_shape = ifelse(is.finite(mcc_shape), 
                                   mcc_shape, NA)) %>%
  summarize(across(.cols = c(sensitivity:mcc_shape), 
                   list(mean = mean), 
                   na.rm = T, 
                   .names = "{.col}.{.fn}")) 


# Summarized over beta, T, and method; pivoted for just shape outliers
sim_results_summary_rTmeth_shape <- sim_results_summary_rTmeth %>%
  dplyr::select(-c(sensitivity.mean:precision_shift.mean, 
                   mcc.mean, mcc_shift.mean)) %>%
  pivot_longer(cols = c(sensitivity_shape.mean:precision_shape.mean,
                        mcc_shape.mean), 
               names_to = "Metric", 
               values_to = "Estimate") %>%
  dplyr::mutate(Method = dplyr::case_when(
    out_meth == "pod_outliers_c1.5" ~ "POD (Tuk)", 
    out_meth == "pod_outliers_user" ~ "POD (user)", 
    out_meth == "muod_outliers_box" ~ "MUOD (box)",
    out_meth == "muod_outliers_tan" ~ "MUOD (tan)",
    out_meth == "tvd_outliers" ~ "TVD")) %>%
  dplyr::filter(Metric == "mcc_shape.mean") %>%
  dplyr::rename("r" = "outlier_rate") %>%
  dplyr::mutate(Out_Type = "Shape")

# sim_results_summary_rTmeth_shape %>%
#   group_by(r, n_smpl_pts, Method) %>%
#   dplyr::filter(Metric == "mcc_shape.mean") %>%
#   dplyr::summarise(avg_mcc = mean(Estimate, na.rm = T)) %>%
#   View()

sim_results_summary_rTmeth_shape$Method <- factor(
  sim_results_summary_rTmeth_shape$Method,
  levels = c("POD (Tuk)", 
             "POD (user)", 
             "MUOD (box)",
             "MUOD (tan)",
             "TVD")
)

# ## Recreates the "Comparing Classification of Shape Outliers" plot in
# ## Figure 8 found in Appendix A of the manuscript
# # Plot, faceted by r, with diagnostic trending over T
# ggplot() +
#   geom_line(aes(x = n_smpl_pts, 
#                 y = Estimate,
#                 color = Method, 
#                 linetype = Method),
#             size = 1.5,
#             data = sim_results_summary_rTmeth_shape) + 
#   geom_point(aes(x = n_smpl_pts, 
#                  y = Estimate, 
#                  color = Method, 
#                  shape = Method), 
#              size = 4,
#              data = sim_results_summary_rTmeth_shape) + 
#   facet_grid(cols = vars(r), 
#              # rows = vars(Metric),
#              scales = "free",
#              labeller = purrr::partial(label_both, sep = " = ")) +
#   scale_color_manual(values = cb_pallette) +
#   labs(y = "Matthews Correlation Coefficient (MCC)",
#        x = "T", 
#        title = "Comparing Classification of Shape Outliers") +
#   theme_bw(base_size = 18)


# Summarized over beta, T, and method; pivoted for just shift outliers
sim_results_summary_rTmeth_shift <- sim_results_summary_rTmeth %>%
  dplyr::select(-c(sensitivity.mean:precision.mean,
                   sensitivity_shape.mean:precision_shape.mean,
                   mcc.mean, mcc_shape.mean)) %>%
  pivot_longer(cols = c(sensitivity_shift.mean:precision_shift.mean,
                        mcc_shift.mean), 
               names_to = "Metric", 
               values_to = "Estimate") %>%
  dplyr::mutate(Method = dplyr::case_when(
    out_meth == "pod_outliers_c1.5" ~ "POD (Tuk)", 
    out_meth == "pod_outliers_user" ~ "POD (user)", 
    out_meth == "muod_outliers_box" ~ "MUOD (box)",
    out_meth == "muod_outliers_tan" ~ "MUOD (tan)",
    out_meth == "tvd_outliers" ~ "TVD")) %>%
  dplyr::filter(Metric == "mcc_shift.mean") %>%
  dplyr::rename("r" = "outlier_rate") %>%
  dplyr::mutate(Out_Type = "Magnitude")

# sim_results_summary_rTmeth_shift %>%
#   group_by(r, n_smpl_pts, Method) %>%
#   dplyr::filter(Metric == "mcc_shift.mean") %>%
#   dplyr::summarise(avg_mcc = mean(Estimate, na.rm = T)) %>%
#   View()

sim_results_summary_rTmeth_shift$Method <- factor(
  sim_results_summary_rTmeth_shift$Method,
  levels = c("POD (Tuk)", 
             "POD (user)", 
             "MUOD (box)",
             "MUOD (tan)",
             "TVD")
)

# ## Recreates the "Comparing Classification of Magnitude Outliers" plot
# ## in Figure 8 found in Appendix A of the manuscript
# # Plot, fcaeted by beta, with diagnostic trending over T
# ggplot() +
#   geom_line(aes(x = n_smpl_pts, 
#                 y = Estimate,
#                 color = Method, 
#                 linetype = Method),
#             size = 1.5,
#             data = sim_results_summary_rTmeth_shift) + 
#   geom_point(aes(x = n_smpl_pts, 
#                  y = Estimate, 
#                  color = Method, 
#                  shape = Method), 
#              size = 4,
#              data = sim_results_summary_rTmeth_shift) + 
#   facet_grid(cols = vars(r), 
#              # rows = vars(Metric),
#              scales = "free",
#              labeller = purrr::partial(label_both, sep = " = ")) +
#   scale_color_manual(values = cb_pallette) +
#   labs(y = "Matthews Correlation Coefficient (MCC)",
#        x = "T", 
#        title = "Comparing Classification of Magnitude Outliers") +
#   theme_bw(base_size = 18)


# Combine into single data frame for plotting
sim_results_summary_rTmeth_all <- dplyr::bind_rows(
  sim_results_summary_rTmeth_shift, sim_results_summary_rTmeth_shape
)

# Figure 8 in Appendix A:
ggplot() +
  geom_line(aes(x = n_smpl_pts, 
                y = Estimate,
                color = Method, 
                linetype = Method),
            size = 1.5,
            data = sim_results_summary_rTmeth_all) + 
  geom_point(aes(x = n_smpl_pts, 
                 y = Estimate, 
                 color = Method, 
                 shape = Method), 
             size = 4,
             data = sim_results_summary_rTmeth_all) + 
  facet_grid(cols = vars(r), 
             rows = vars(Out_Type),
             scales = "free",
             labeller = label_bquote(
               cols = "r = "*.(r), 
               rows = .(Out_Type)*" Outlier")) +
  scale_color_manual(values = cb_pallette) +
  labs(y = "Matthews Correlation Coefficient (MCC)",
       x = "T", 
       title = "Comparing Classification of Magnitude Outliers") +
  theme_bw(base_size = 18)
#####################################################################