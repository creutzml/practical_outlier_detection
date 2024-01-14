#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### June 20th, 2023                                               ###
###                                                               ###
###   Simulation script to investigate and test the ability       ###
### of each method for classifying the type of functional         ###
### outlier present. Adding a simulation model that creates       ###
### different combinations of outliers that are both magnitude    ###
### and shape outliers.                                           ###
#####################################################################

## Libraries
library(tidyverse)
library(MASS)
library(fda)
library(fdaoutlier)
library(pracma)
library(R.utils)
options(dplyr.summarise.inform = FALSE)
library(here)
library(progress)

## Backend setup
#####################################################################
## Block below is commented out, since that was for use with bash 
## scripting to an HPC
# # Read in command arguments from bash files
# args <- commandArgs(trailingOnly=TRUE)
# if (length(args) != 6) {
#   stop("Looping parameters failed to initialize in 'commandArgs'.")
# }
# 
# # Make new objects from input
# sim_mod <- as.character(args[1])
# n_obs <- as.numeric(args[2])
# n_smpl_pts <- as.numeric(args[3])
# out_rate <- as.numeric(args[4])
# beta_cov <- as.numeric(args[5])

# Set the simulation parameters
{sim_mod = "1"; n_obs = 30; n_smpl_pts = 30; out_rate = 0.05;
beta_cov = .1}

# File name for saved results later
file_name <- paste0(
  "Simulation_Model", sim_mod, "_n_obs=", n_obs, "_n_smpl_pts=", 
  n_smpl_pts, "_out_rate=", out_rate, "_cov_beta=", beta_cov,
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
if (sim_mod == "B") {
  source(file.path(dir_path, "simulation_model_b.R"))
}
if (sim_mod == "B2") {
  source(file.path(dir_path, "simulation_model_b_only.R"))
}

# Specify the simulation function for getting new data
sim_mod_fct <- match.fun(paste0("simulation_model", sim_mod))
#####################################################################



## Start the simulation
#####################################################################
# Number of simulation iterations
n_iters <- 1

# Number of methods being used:
n_meths <- 5

# Number of total observations being recorded
n_length <- n_iters*n_meths

# Empty data frame to store results
sim_results <- data.frame(
  iter_n = rep(1:n_iters, each = n_meths),
  sim_model = rep(sim_mod, n_length),
  n_obs = rep(n_obs, n_length),
  n_smpl_pts = rep(n_smpl_pts, n_length),
  outlier_rate = rep(out_rate, n_length),
  cov_beta = rep(beta_cov, n_length),
  out_meth = rep("", n_length),
  n_outliers = rep(0, n_length),
  n_found = rep(0, n_length),
  n_tp = rep(0, n_length),
  n_fp = rep(0, n_length),
  n_fn = rep(0, n_length),
  n_shift = rep(0, n_length),
  n_found_shift = rep(0, n_length),
  n_tp_shift = rep(0, n_length),
  n_fp_shift = rep(0, n_length),
  n_fn_shift = rep(0, n_length),
  n_shape = rep(0, n_length),
  n_found_shape = rep(0, n_length),
  n_tp_shape = rep(0, n_length),
  n_fp_shape = rep(0, n_length),
  n_fn_shape = rep(0, n_length),
  run_times = rep(0, n_length),
  stringsAsFactors = F
)


# sss_time <- Sys.time()
# Progress bar for tracking the simulation progress
p_bar <- progress_bar$new(total = n_iters)

# Code structure from bash scripting
# # Add a stopper for the max wall time that slurm allows, in order to
# # save whatever iterations completed by that time
# withTimeout({
#   
# Loop over iterations
for (iter in 1:n_iters) {
  
  
  ## Simulate the data used for testing
  ###################################################################
  # Use the function to simulate data
  test_data_obj <- sim_mod_fct(
    n = n_obs,
    p = n_smpl_pts,
    outlier_rate = out_rate,
    cov_beta = beta_cov,
    plot = F
  )
  
  # Obtain just the matrix of data
  test_data <- test_data_obj[[1]]
  
  # Obtain the list of true outliers and a count
  if (sim_mod %in% c("B", "B2")) {
    true_outs <- test_data_obj[[2]][[2]]
    # test_data_obj$true_outliers$true_outliers
    n_outliers <- test_data_obj[[2]][[1]]
    # test_data_obj$true_outliers$n_outliers
  } else {
    true_outs <- test_data_obj[[2]]
    n_outliers <- length(true_outs)
  }
  
  # Empty vectors
  true_shift <- true_shape <- c()
  n_shift <- n_shape <- 0
  
  # Format the outlier names
  true_outs <- paste0("Obs", str_pad(true_outs, 3, "left", "0"))
  
  # Specify shape vs shift
  if (sim_mod == "B") {
    true_shift <- test_data_obj[[2]][[3]]
    # test_data_obj$true_outliers$shift_outliers
    true_shift <- paste0("Obs", str_pad(true_shift, 3, "left", "0"))
    n_shift <- length(true_shift)
    true_shape <- test_data_obj[[2]][[4]]
    # test_data_obj$true_outliers$shape_outliers
    true_shape <- paste0("Obs", str_pad(true_shape, 3, "left", "0"))
    n_shape <- length(true_shape)
  } else if (sim_mod == "B2") {
    true_shift <- test_data_obj[[2]][[2]]
    true_shift <- paste0("Obs", str_pad(true_shift, 3, "left", "0"))
    n_shift <- length(true_shift)
    true_shape <- test_data_obj[[2]][[2]]
    true_shape <- paste0("Obs", str_pad(true_shape, 3, "left", "0"))
    n_shape <- length(true_shape)
  } else if (sim_mod == "7") {
    true_shape <- true_outs
    n_shape <- n_outliers
  } else {
    true_shift <- true_outs
    n_shift <- n_outliers
  }
  ###################################################################
  
  
  
  ## Test the method using the simulated data
  ###################################################################
  
  ### Need to run the five methods, which results in 13 unique lists
  ### of outliers
  
  ## Run easy outlier method (new rule) with user specified thresh
  #################################################################
  temp_pod_user_s <- Sys.time()
  temp_pod_user <- pod_fda(test_data = test_data,
                           cutoff = 1 - out_rate)
  
  # Total run time
  temp_pod_user_time <- as.numeric(
    difftime(Sys.time(), temp_pod_user_s, units = "secs")
  )
  
  # Grab the results
  temp_pod_user_outliers <- temp_pod_user$outliers_found
  
  # Empty vectors
  temp_pod_user_mag <- temp_pod_user_shape <- c()
  
  if (length(temp_pod_user_outliers) != 0) {
    # Shape vs Shift outliers
    temp_pod_user_mag <- unique(unlist(c(
      temp_pod_user$outliers_class %>%
        dplyr::filter(class %in% c("Magnitude", "Both")) %>%
        dplyr::select(obs_idx)
    )))
    
    temp_pod_user_shape <- unique(unlist(c(
      temp_pod_user$outliers_class %>%
        dplyr::filter(class %in% c("Shape", "Both")) %>%
        dplyr::select(obs_idx)
    )))
  }
  #################################################################
  
  
  ## Run easy outlier method (new rule) with Tukey threshold
  #################################################################
  temp_pod_c1.5_s <- Sys.time()
  temp_pod_c1.5 <- pod_fda(test_data = test_data,
                           cutoff = "classical1.5")
  
  # Total run time
  temp_pod_c1.5_time <- as.numeric(
    difftime(Sys.time(), temp_pod_c1.5_s, units = "secs")
  )
  
  # Grab the results
  temp_pod_c1.5_outliers <- temp_pod_c1.5$outliers_found
  
  # Empty vectors
  temp_pod_c1.5_mag <- temp_pod_c1.5_shape <- c()
  
  if (length(temp_pod_c1.5_outliers) != 0) {
    # Shape vs Shift outliers
    temp_pod_c1.5_mag <- unique(unlist(c(
      temp_pod_c1.5$outliers_class %>%
        dplyr::filter(class %in% c("Magnitude", "Both")) %>%
        dplyr::select(obs_idx)
    )))
    
    temp_pod_c1.5_shape <- unique(unlist(c(
      temp_pod_c1.5$outliers_class %>%
        dplyr::filter(class %in% c("Shape", "Both")) %>%
        dplyr::select(obs_idx)
    )))
  }
  #################################################################
  
  
  ## Run TVD on the data
  #################################################################
  tvd_s_time <- Sys.time()
  tvd_test <- tvdmss_mc(dts = test_data)
  
  # Total run time
  tvd_time <- difftime(Sys.time(), tvd_s_time, units = "secs")
  
  # Grab the results
  tvd_outliers <- paste0("Obs",
                         str_pad(
                           tvd_test$outliers, 3, "left", "0"
                         ))
  
  # Identify shape and shift outliers
  tvd_shift <- paste0("Obs",
                      str_pad(
                        tvd_test$magnitude_outliers, 3, "left", "0"
                      ))
  tvd_shape <- paste0("Obs",
                      str_pad(
                        tvd_test$shape_outliers, 3, "left", "0"
                      ))
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
  
  ## something weird happening on Alpine, where there's times that
  ## muod returns a list for outliers, then returns a data frame
  ## other times
  
  if (is.list(muod_test_tan$outliers)) {
    # Identify shape vs shift outliers for both methods
    muod_outliers_tan_shift <- paste0(
      "Obs", str_pad(muod_test_tan$outliers$magnitude,
                     3, "left", "0")
    )
    muod_outliers_tan_shape <- paste0(
      "Obs", str_pad(unique(c(muod_test_tan$outliers$shape,
                              muod_test_tan$outliers$amplitude)),
                     3, "left", "0")
    )
    
  } else {
    # Identify shape vs shift outliers for both methods
    muod_outliers_tan_shift <- tryCatch({
      paste0("Obs", str_pad(muod_test_tan$outliers[,3], #$magnitude,
                            3, "left", "0")
      )
    }, error = function(cond) {
      return("")
    }, finally = "")
    
    muod_outliers_tan_shape <- tryCatch({
      paste0("Obs", str_pad(unique(c(muod_test_tan$outliers[,1], #$shape,
                                     muod_test_tan$outliers[,2])), #$amplitude)),
                            3, "left", "0")
      )
    }, error = function(cond) {
      return("")
    }, finally = "")
  }
  
  if (is.list(muod_test_box$outliers)) {
    muod_outliers_box_shift <- paste0(
      "Obs", str_pad(muod_test_box$outliers$magnitude,
                     3, "left", "0")
    )
    muod_outliers_box_shape <- paste0(
      "Obs", str_pad(unique(c(muod_test_box$outliers$shape,
                              muod_test_box$outliers$amplitude)),
                     3, "left", "0")
    )
    
  } else {
    muod_outliers_box_shift <- tryCatch({
      paste0("Obs", str_pad(muod_test_box$outliers[,3], #$magnitude,
                            3, "left", "0")
      )
    }, error = function(cond) {
      return("")
    }, finally = "")
    
    muod_outliers_box_shape <- tryCatch({
      paste0("Obs", str_pad(unique(c(muod_test_box$outliers[,1], #$shape,
                                     muod_test_box$outliers[,2])), #$amplitude)),
                            3, "left", "0")
      )
    }, error = function(cond) {
      return("")
    }, finally = "")
  }
  #################################################################
  
  
  
  ### Make some calculations
  #################################################################
  
  ## First put all the outliers into a list
  test_outlier_list <- list(
    pod_outliers_user = list(temp_pod_user_outliers,
                             temp_pod_user_mag,
                             temp_pod_user_shape),
    pod_outliers_c1.5 = list(temp_pod_c1.5_outliers,
                             temp_pod_c1.5_mag,
                             temp_pod_c1.5_shape),
    tvd_outliers = list(tvd_outliers,
                        tvd_shift,
                        tvd_shape),
    muod_outliers_tan = list(muod_outliers_tan,
                             muod_outliers_tan_shift,
                             muod_outliers_tan_shape),
    muod_outliers_box = list(muod_outliers_box,
                             muod_outliers_box_shift,
                             muod_outliers_box_shape)
  )
  
  
  ## Next, calculate some things
  # Run through lapply over the list
  outlier_results <- lapply(test_outlier_list, FUN = function(x) {
    # Outliers overall
    n_found <- length(x[[1]])
    n_tp <- sum(x[[1]] %in% true_outs)
    n_fp <- sum(!(x[[1]] %in% true_outs))
    n_fn <- sum(!(true_outs %in% x[[1]]))
    
    # Shift outliers
    n_found_shift <- length(x[[2]])
    n_tp_shift <- sum(x[[2]] %in% true_shift)
    n_fp_shift <- sum(!(x[[2]] %in% true_shift))
    n_fn_shift <- sum(!(true_shift %in% x[[2]]))
    
    # Shape outliers
    n_found_shape <- length(x[[3]])
    n_tp_shape <- sum(x[[3]] %in% true_shape)
    n_fp_shape <- sum(!(x[[3]] %in% true_shape))
    n_fn_shape <- sum(!(true_shape %in% x[[3]]))
    
    return(c(n_found, n_tp, n_fp, n_fn,
             n_shift, n_found_shift, n_tp_shift, n_fp_shift,
             n_fn_shift,
             n_shape, n_found_shape, n_tp_shape, n_fp_shape,
             n_fn_shape))
  })
  
  # Combine into a matrix
  outlier_results_mat <- cbind(names(outlier_results),
                               rep(n_outliers, n_meths),
                               do.call(rbind, outlier_results),
                               round(c(temp_pod_user_time,
                                       temp_pod_c1.5_time,
                                       tvd_time,
                                       muod_time,
                                       muod_time), 3))
  ###################################################################
  
  # Print the index number to know how far the simulation made it
  base::cat(iter, ", ", sep = "")
  
  # Index numbers
  s_idx <- 1 + (iter - 1)*n_meths
  e_idx <- iter*n_meths
  
  # Update the data frame with the results
  sim_results[s_idx:e_idx, 7:23] <- outlier_results_mat
  
  p_bar$tick()
}
  
# Bash script code
#   # Specify the time limit and what to do
# }, timeout = 80000, onTimeout = "warning")
#####################################################################

# cat("Run time: ", Sys.time() - sss_time)

# Post simulation computations
#####################################################################
# Calculate AUC
sim_results_sum <- sim_results %>%
  # Remove any iterations that did not complete before timeout
  dplyr::filter(iter_n <= iter) %>%
  dplyr::mutate(n_outliers = as.numeric(n_outliers), 
                n_found = as.numeric(n_found), 
                n_tp = as.numeric(n_tp), 
                n_fp = as.numeric(n_fp), 
                n_fn = as.numeric(n_fn),
                n_shift = as.numeric(n_shift), 
                n_found_shift = as.numeric(n_found_shift), 
                n_tp_shift = as.numeric(n_tp_shift), 
                n_fp_shift = as.numeric(n_fp_shift), 
                n_fn_shift = as.numeric(n_fn_shift),
                n_shape = as.numeric(n_shape), 
                n_found_shape = as.numeric(n_found_shape), 
                n_tp_shape = as.numeric(n_tp_shape), 
                n_fp_shape = as.numeric(n_fp_shape), 
                n_fn_shape = as.numeric(n_fn_shape),
                run_times = as.numeric(run_times)) %>%
  dplyr::mutate(n_tn = n_obs - n_fn - n_tp - n_fp, 
                n_tn_shift = n_obs - n_fn_shift - 
                  n_tp_shift - n_fp_shift, 
                n_tn_shape = n_obs - n_fn_shape - 
                  n_tp_shape - n_fp_shape) %>%
  dplyr::mutate(sensitivity = n_tp/(n_tp + n_fn), 
                specificity = n_tn/(n_fp + n_tn), 
                accuracy = (n_tp + n_tn)/n_obs, 
                precision = n_tp/(n_tp + n_fp), 
                sensitivity_shift = n_tp_shift/
                  (n_tp_shift + n_fn_shift), 
                specificity_shift = n_tn_shift/
                  (n_fp_shift + n_tn_shift), 
                accuracy_shift = (n_tp_shift + n_tn_shift)/n_obs, 
                precision_shift = n_tp_shift/
                  (n_tp_shift + n_fp_shift), 
                sensitivity_shape = n_tp_shape/
                  (n_tp_shape + n_fn_shape), 
                specificity_shape = n_tn_shape/
                  (n_fp_shape + n_tn_shape), 
                accuracy_shape = (n_tp_shape + n_tn_shape)/n_obs, 
                precision_shape = n_tp_shape/
                  (n_tp_shape + n_fp_shape))

#####################################################################



## Save the results to a data frame
#####################################################################
# Create a new folder in the data directory for the new results
curr_date <- gsub("-", "_", Sys.Date())
dir.create(file.path(data_path, paste0("sim_type_", curr_date)), 
           showWarnings = FALSE)

# Save the new results
save_path <- file.path(data_path, 
                       paste0("sim_type_", curr_date), 
                       file_name)

base::save(sim_results_sum, file = save_path)
#####################################################################

