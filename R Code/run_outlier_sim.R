#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### March 23rd, 2023                                              ###
###                                                               ###
###   Testing the accuracy of practical functional outlier        ###
### detection through a simulation, using data generated from the ###
### simulation models (1-9) provided in package `fdaoutlier`.     ###
### This script is setup to be implemented using HPC alpine from  ###
### CU Boulder.                                                   ###
#####################################################################

## Libraries
library(tidyverse)
library(MASS)
library(fda)
library(fdaoutlier)
library(pracma)
library(R.utils)
options(dplyr.summarise.inform = FALSE)
library(progress)


## Backend setup
#####################################################################
## Code that was necessary to run the simulation on an HPC with the 
## use of bash scripting
# # Read in command arguments from bash files
# args <- commandArgs(trailingOnly=TRUE)
# if (length(args) != 6) {
#   stop("Looping parameters failed to initialize in 'commandArgs'.")
# }
# 
# # Make new objects from input
# sim_mod <- as.numeric(args[1])
# n_obs <- as.numeric(args[2])
# n_smpl_pts <- as.numeric(args[3])
# out_rate <- as.numeric(args[4])
# alpha_cov <- as.numeric(args[5])


# Simulation parameters to use
{sim_mod = 4; n_obs = 100; n_smpl_pts = 100; out_rate = 0.05; 
alpha_cov = 1}


# File name for saved results later
file_name <- paste0(
  "Simulation_Model", sim_mod, "_n_obs=", n_obs, "_n_smpl_pts=", 
  n_smpl_pts, "_out_rate=", out_rate, "_cov_alpha=", alpha_cov,
  ".RData"
)

base::cat(file_name, "\n")

# Directories
dir_path <- file.path(here::here(), "R Code")
data_path <- file.path(here::here(), "Data")
#####################################################################



## Helper function and main method
#####################################################################
# Source the functions
source(file.path(dir_path, "practical_outlier.R"))
source(file.path(dir_path, "tvdmss_mc.R"))

# Specify the simulation function for getting new data
sim_mod_fct <- match.fun(paste0("simulation_model", sim_mod))
#####################################################################



## Start the simulation
#####################################################################
# Number of simulation iterations
n_iters <- 1

# Number of methods being used:
n_meths <- 6

# Number of total observations being recorded
n_length <- n_iters*n_meths

# Empty data frame to store results
sim_results <- data.frame(
  iter_n = rep(1:n_iters, each = n_meths),
  sim_model = rep(sim_mod, n_length),
  n_obs = rep(n_obs, n_length),
  n_smpl_pts = rep(n_smpl_pts, n_length),
  outlier_rate = rep(out_rate, n_length),
  cov_alpha = rep(alpha_cov, n_length),
  out_meth = rep("", n_length),
  n_outliers = rep(0, n_length),
  n_found = rep(0, n_length),
  n_tp = rep(0, n_length),
  n_fp = rep(0, n_length),
  n_fn = rep(0, n_length),
  run_times = rep(0, n_length),
  stringsAsFactors = F
)

# sss_time <- Sys.time()

# Progress bar to keep track of simulation progress
p_bar <- progress_bar$new(total = n_iters)

# Code leftover from bash script
# # Add a stopper for the max wall time that slurm allows, in order to
# # save whatever iterations completed by that time
# withTimeout({
  
# Loop over iterations
for (iter in 1:n_iters) {
  
  
  ## Simulate the data used for testing
  ###################################################################
  # Use the function to simulate data
  test_data_obj <- sim_mod_fct(
    n = n_obs, p = n_smpl_pts, outlier_rate = out_rate, 
    cov_alpha = alpha_cov, plot = F
  )
  
  # Obtain just the matrix of data
  test_data <- test_data_obj$data
  
  # Obtain the list of true outliers and a count
  true_outs <- test_data_obj$true_outliers
  n_outliers <- length(true_outs)
  
  # Format the outlier names
  true_outs <- paste0("Obs", str_pad(true_outs, 3, "left", "0"))
  ###################################################################
  
  
  
  ## Test the method using the simulated data
  ###################################################################
  
  ### Need to run the five methods, which results in 13 unique lists
  ### of outliers
  
  ## Run easy outlier method (new rule) with user specified thresh
  #################################################################
  easy_s_time <- Sys.time()
  easy_test_new <- pod_fda(
    test_data = test_data, 
    cutoff = 1 - out_rate 
  )
  
  # Total run time
  easy_time <- as.numeric(
    difftime(Sys.time(), easy_s_time, units = "secs")
  )
  
  # Grab the results
  easy_outliers_new_tuk <- easy_test_new$outliers_found
  #################################################################
  
  ## Run easy outlier method (new rule) with user specified thresh
  #################################################################
  easy_s_time_c1.5 <- Sys.time()
  easy_test_c1.5 <- pod_fda(
    test_data = test_data, 
    cutoff = "classical1.5" 
  )
  
  # Total run time
  easy_time_c1.5 <- as.numeric(
    difftime(Sys.time(), easy_s_time_c1.5, units = "secs")
  )
  
  # Grab the results
  easy_outliers_c1.5_tuk <- easy_test_c1.5$outliers_found
  #################################################################
  
  ## Outlier detection with MS-Plot
  #################################################################
  ms_s_time <- Sys.time()
  ms_test <- msplot(dts = test_data, 
                    return_mvdir = F,
                    plot = F)
  
  # Total run time
  ms_time <- difftime(Sys.time(), ms_s_time, units = "secs")
  
  # Grab the results
  ms_outliers <- paste0("Obs", 
                        str_pad(ms_test$outliers, 3, "left", "0"))
  #################################################################
  
  
  ## Run TVD on the data
  #################################################################
  tvd_s_time <- Sys.time()
  tvd_test <- tvdmss_mc(dts = test_data)
  
  # Total run time
  tvd_time <- difftime(Sys.time(), tvd_s_time, units = "secs")
  
  # Grab the results
  tvd_outliers <- paste0("Obs", 
                         str_pad(tvd_test$outliers, 3, "left", "0"))
  #################################################################
  
  
  ## Run the MUOD method
  #################################################################
  muod_s_time <- Sys.time()
  muod_test_tan <- muod(dts = test_data, cut_method = "tangent")
  muod_test_box <- muod(dts = test_data, cut_method = "boxplot")
  
  # Total run time
  muod_time <- difftime(Sys.time(), muod_s_time, units = "secs")
  
  # Grab the results
  muod_outliers_tan <- paste0(
    "Obs", str_pad(unique(unlist(muod_test_tan$outliers)), 
                   3, "left", "0")
  )
  
  muod_outliers_box <- paste0(
    "Obs", str_pad(unique(unlist(muod_test_box$outliers)), 
                   3, "left", "0")
  )
  #################################################################
  
  
  
  ### Make some calculations
  #################################################################
  
  ## First put all the outliers into a list
  test_outlier_list <- list(
    easy_outliers_new_tuk = easy_outliers_new_tuk, 
    easy_outliers_c1.5_tuk = easy_outliers_c1.5_tuk, 
    ms_outliers = ms_outliers,
    tvd_outliers = tvd_outliers,
    muod_outliers_tan = muod_outliers_tan,
    muod_outliers_box = muod_outliers_box
  )
  
  
  ## Next, calculate some things
  # Run through lapply over the list
  outlier_results <- lapply(test_outlier_list, FUN = function(x) {
    n_found <- length(x)
    n_tp <- sum(x %in% true_outs)
    n_fp <- sum(!(x %in% true_outs))
    n_fn <- sum(!(true_outs %in% x))
    
    return(c(n_found, n_tp, n_fp, n_fn))
  })
  
  # Combine into a matrix
  outlier_results_mat <- cbind(names(outlier_results),
                               rep(n_outliers, n_meths),
                               do.call(rbind, outlier_results),
                               round(c(easy_time,
                                       easy_time_c1.5,
                                       ms_time, tvd_time, 
                                       muod_time, muod_time), 3))
  ###################################################################
  
  # Print the index number to know how far the simulation made it
  cat(iter, ", ", sep = "")
  
  # Index numbers
  s_idx <- 1 + (iter - 1)*n_meths
  e_idx <- iter*n_meths
  
  # Update the data frame with the results
  sim_results[s_idx:e_idx, 7:13] <- outlier_results_mat
  
  p_bar$tick()
}
  
#   # Specify the time limit and what to do
# }, timeout = 85800, onTimeout = "warning")
#####################################################################

# cat("Run time: ", Sys.time() - sss_time)

# Post simulation computations
#####################################################################
# Calculate AUC
sim_results_sum <- sim_results %>%
  # Remove rows that didn't run before the timeout:
  dplyr::filter(iter_n <= iter) %>%
  dplyr::mutate(n_outliers = as.numeric(n_outliers), 
                n_found = as.numeric(n_found), 
                n_tp = as.numeric(n_tp), 
                n_fp = as.numeric(n_fp), 
                n_fn = as.numeric(n_fn), 
                run_times = as.numeric(run_times)) %>%
  dplyr::mutate(n_tn = n_obs - n_fn - n_tp - n_fp) %>%
  dplyr::mutate(sensitivity = n_tp/(n_tp + n_fn), 
                specificity = n_tn/(n_fp + n_tn), 
                accuracy = (n_tp + n_tn)/n_obs, 
                precision = n_tp/(n_tp + n_fp))

#####################################################################



## Save the results to a data frame
#####################################################################
# Create a new folder in the data directory for the new results
curr_date <- gsub("-", "_", Sys.Date())
dir.create(file.path(data_path, paste0("sim_", curr_date)), 
           showWarnings = FALSE)

# Save the new results
save_path <- file.path(data_path, 
                       paste0("sim_", curr_date), 
                       file_name)
base::save(sim_results_sum, file = save_path)
#####################################################################

