#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### July 6th, 2022                                                ###
###                                                               ###
###   Using fast and fair prediction bands to develop a           ###
### functional outlier detection tool, through the use of         ###
### resampling techniques. This function is explained in detail   ###
### in Chapter 2, Creutzinger (2024+).                            ###
#####################################################################


# Packages
require(devtools)
install_github("creutzml/ffscbExtra")
library(tidyverse)
library(fdaoutlier)
library(fda)
library(ffscb)
library(ffscbExtra)



## Develop the function for ff_outlier
#####################################################################
pbod <- function(
    test_data, ho_pct = 0.5, n_steps1 = NULL, n_steps2 = NULL, 
    band_alpha = 0.05, band_n_int = 3, cutoff = 0.90
) {
  # test_data (matrix/df): matrix or data frame of data from which 
  #   outliers should be identified, where nrows = number of 
  #   observations and ncols = number of sampling points
  # ho_pct (num): proportion between 0 and 1 for how much data should
  #   be held out when creating the prediction band in each step
  # n_steps1 (num): how many iterations of the first resampling
  #   period should be done?
  # n_steps2 (num): how many iterations of the second resampling
  #   period should be done?
  # band_alpha (num): number between 0 and 1 denoting the level of 
  #   significance used for creating the simultaneous band
  # band_n_int (num): number of intervals used to create FF band. A smaller
  #   number of intervals results in a narrower band overall
  
  
  # Dimensions of data
  n_obs <- nrow(test_data)
  n_smpl_pts <- ncol(test_data)
  
  # Set the default values of n_obs and n_smpl_pts
  if (is.null(n_steps1)) n_steps1 <- n_obs
  if (is.null(n_steps2)) n_steps2 <- 2*n_obs
  
  # Choose t or Gaussian for user based on sample size
  if (n_obs*(1-ho_pct) >= 100) {
    band_type <- "FFSCB.z"
  } else {
    band_type <- "FFSCB.t"
  }
  
  # Vector of weights inversely related to the probability of an 
  # observation being an outlier
  obs_out_weights <- rep(1, n_obs)
  
  # Vector recording the average number of points along the curve 
  # that are marked as outlying
  obs_pts_out <- rep(NA, n_obs)
  
  # Start the first resampling procedure, randomly sampling a set of
  # observations without consideration of weights. This step will 
  # continuously update weights as it runs
  for (i in 1:n_steps1) {
    
    # Specify training and test sets
    test_set_idx <- sample(x = 1:n_obs, size = ceiling(ho_pct*n_obs))
    train_set <- test_data[-test_set_idx,]
    test_set <- test_data[test_set_idx,]
    train_df <- nrow(train_set) - 1
      
    # Calculate mean, covariance, and tau of training set
    train_mu <- colMeans(train_set)
    train_cov <- crossprod(train_set - train_mu)/train_df
    train_cov_mu <- train_cov/(train_df + 1)
    train_tau <- ffscb::tau_fun(t(train_set))
    
    # Create the prediction band
    train_band <- confidence_band(
      x = train_mu,
      cov.x = train_cov_mu, 
      tau = train_tau, 
      df = train_df,
      type = band_type,
      conf.level = 1 - band_alpha,
      n_int = band_n_int,
      int.type = "prediction",
      n.curves = n_obs
    )
    
    # Compare the test set to the band created, and make the weight
    # equal to `1 - (# of points outside band)/(# smpl pts)`
    new_weights <- apply(test_set, MARGIN = 1, FUN = function(x) {
      
      # compare vector to band
      n_pts_out <- sum(x < train_band[,3] | x > train_band[,2])
      
      # Return the new weight
      return(1 - n_pts_out/n_smpl_pts)
      
    })
    
    # Update the weight matrix
    obs_out_weights[test_set_idx] <- rowMeans(
      cbind(obs_out_weights[test_set_idx], new_weights)
    )
  }
  
  
  # Resampling procedure step 2: resample observaitons using the 
  # weights created in the first step
  for (i in 1:n_steps2) {
    
    # Specify training and test sets
    test_set_idx <- sample(
      x = 1:n_obs, 
      size = ceiling(ho_pct*n_obs), 
      prob = obs_out_weights
    )
    train_set <- test_data[-test_set_idx,]
    test_set <- test_data[test_set_idx,]
    train_df <- nrow(train_set) - 1
    
    # Calculate mean, covariance, and tau of training set
    train_mu <- colMeans(train_set)
    train_cov <- crossprod(train_set - train_mu)/train_df
    train_cov_mu <- train_cov/(train_df + 1)
    train_tau <- ffscb::tau_fun(t(train_set))
    
    # Create the prediction band
    train_band <- confidence_band(
      x = train_mu,
      cov.x = train_cov_mu, 
      tau = train_tau, 
      df = train_df,
      type = band_type,
      conf.level = 1 - band_alpha,
      n_int = band_n_int,
      int.type = "prediction",
      n.curves = n_obs
    )
    
    # Compare test set to train band and record number of points 
    # along the curve that are outlying
    pts_out <- apply(test_set, MARGIN = 1, FUN = function(x) {
      
      # compare vector to band
      n_pts_out <- sum(x < train_band[,3] | x > train_band[,2])
      
      # Return the new weight
      return(n_pts_out)
      
    })
    
    # Update vector recording number of points outlying
    obs_pts_out[test_set_idx] <- rowMeans(
      cbind((i - 1)*obs_pts_out[test_set_idx], pts_out),
      na.rm = T
    )
  }
  
  # Create a data frame of results
  obs_pts_out_df <- data.frame(
    obs_idx = paste0("Obs", str_pad(1:n_obs, 3, "left", "0")), 
    obs_pts_out = obs_pts_out
  )
  
  # Identify the outliers:
  out_names <- unique(unlist(c(
    obs_pts_out_df %>%
      dplyr::filter(obs_pts_out >= 
                      quantile(obs_pts_out, cutoff, na.rm = T)) %>%
      dplyr::select(obs_idx)
  )))
  
  # Return the vector of obs_pts_out 
  return(list(obs_pts_out_df = obs_pts_out_df, 
              outliers = out_names))
}




# ## Short test example:
# # Sample the data
# test_data_obj <- fdaoutlier::simulation_model9(n = 300, p = 50)
# test_data <- test_data_obj$data
# true_outs <- test_data_obj$true_outliers
# 
# # Implement PBOD
# test_run <- pbod(test_data = test_data, ho_pct = 0.05,
#                  n_steps1 = 10, n_steps2 = 10,
#                  band_alpha = 0.10)

