#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### July 19th, 2023                                               ###
###                                                               ###
###   Latest rendition of POD.                                    ###
#####################################################################

## Libraries
library(tidyverse)
library(MASS)
library(fda)
library(fdaoutlier)
library(pracma)
options(dplyr.summarise.inform = FALSE)





## Helper function and main method
#####################################################################
# Roughness calculation function
rough <- function(x) {
  l_x <- length(x)
  sum((x[3:l_x] - 2*x[2:(l_x - 1)] + x[1:(l_x - 2)])^2/4)
}

# Area under the curve using spline interpolation
auc_mc <- function(x, y) {
  
  # Create function to integrate over using spline interpolation
  int_fct <- stats::splinefun(x = x, y = y, method = "natural")
  
  # Use the integrate function
  int_val <- stats::integrate(
    int_fct, lower = min(x), upper = max(x), 
    subdivisions = 500, stop.on.error = F
  )
  
  # Return the value
  return(int_val$value)
}
#####################################################################



# Main function
#####################################################################
pod_fda <- function(test_data, cutoff = .9) {
  # test_data (mat): matrix of discrete functional observations, with
  #  the number of rows equivalent to the number of observations, and
  #  the number of columns equivalent to the number of sampling 
  #  points along the domain
  # plot_it (bool): should a plot be given of the curves, with the 
  #  outliers identified in a different color?
  # cutoff (num): value between 0 and 1 specifying the presumed 
  #  percentage of observations that are outliers

  
  ## Step 1: Split the sampling domain into intervals and calculate
  ## summary statistics over each interval for all observations
  ###################################################################
  # Acquire number of curves and sampling points
  n_sp <- ncol(test_data)
  n_curves <- nrow(test_data)
  
  # Need to determine interval sizes based on number of sampling pts
  if (n_sp >= 60) {
    int_size1 <- 15
    int_size2 <- 20
  } else if (n_sp >= 45) {
    int_size1 <- 3
    int_size2 <- 15
  } else {
    int_size1 <- 3
    int_size2 <- 8
  }
  
  # Set the number of intervals n_int_eigth <- n_sp*0.125
  int_cat1 <- cut(x = 1:n_sp, breaks = int_size1)
  int_cat2 <- cut(x = 1:n_sp, breaks = int_size2)
  
  # Now compute the summary statistics for each interval, for each 
  # functional observation
  test_profile_int <- test_data %>%
    as.data.frame() %>%
    dplyr::mutate(
      obs_idx = paste0("Obs", str_pad(1:n_curves, 3, "left", "0"))
    ) %>%
    tidyr::pivot_longer(cols = c(dplyr::everything(),-obs_idx)) %>%
    dplyr::mutate(`int1` = rep(int_cat1, n_curves), 
                  `int2` = rep(int_cat2, n_curves)) %>%
    tidyr::pivot_longer(cols = c(`int1`, `int2`), 
                        names_to = "n_intervals", 
                        values_to = "interval") %>%
    dplyr::group_by(obs_idx, n_intervals, interval) %>%
    dplyr::summarise(min_value = min(value), 
                     max_value = max(value), 
                     mean_value = mean(value),
                     med_value = median(value),
                     obs_range = max_value - min_value,
                     obs_rough = rough(value),
                     obs_auc = auc_mc(x = 1:length(value), 
                                      y = value),
                     obs_var = var(value),
                     obs_coef_var = sqrt(obs_var)/mean_value*100)
  ###################################################################
  
  
  ## Step 2: Use Tukey's classical boxplot to identify extreme 
  ## summary statistics relative to each interval and type
  ###################################################################
  # Find the upper and lower fence of each summary stat calculated, 
  # for each interval
  test_summary_grp <- test_profile_int %>%
    tidyr::pivot_longer(c(-obs_idx, -n_intervals, -interval), 
                        names_to = "summary_stat") %>%
    dplyr::group_by(summary_stat, n_intervals, interval)

  # Make the calculations for fences
  test_summary <- test_summary_grp %>%
    dplyr::summarize(
      q1_stat = quantile(value, probs = 0.25),
      q3_stat = quantile(value, probs = 0.75), 
      iqr_stat = q3_stat - q1_stat,
      u_fence_tuk = q3_stat + 1.5*iqr_stat, 
      l_fence_tuk = q1_stat - 1.5*iqr_stat
    )
  
  # Compare summary stat values to the fences created
  test_outliers <- test_summary_grp %>%
    dplyr::arrange(obs_idx, n_intervals, interval) %>%
    dplyr::left_join(
      test_summary, by = c("n_intervals", "interval", "summary_stat")
    ) %>%
    dplyr::mutate(out_tuk = (value < l_fence_tuk | 
                               value > u_fence_tuk))
  
  # Total number of outlying values per curve (adding up the whole 
  # domain and intervals)
  test_outliers_tot <- test_outliers %>%
    dplyr::group_by(obs_idx) %>%
    dplyr::summarize(n_tot_out_tuk = sum(out_tuk)) 
  ###################################################################
  
  
  ## Step 3: Identify outlying observations based on cutoff
  ###################################################################
  ## Determine the list of outliers based on the decided cutoff
  if (is.numeric(cutoff)) {
    # Acquire the names of observations identified as outliers, by 
    # specifying an outlier as those with a total number of outlying
    # values greater than or equal to the percentile
    outliers_tot_tuk <- unique(unlist(c(
      test_outliers_tot %>%
        dplyr::filter(
          n_tot_out_tuk >= quantile(n_tot_out_tuk, cutoff),
        ) %>%
        dplyr::select(obs_idx)
    )))
    
  } else {
    # Acquire the names of observations identified as outliers, by 
    # specifying an outlier as those with a total number of outlying
    # values outside the fences of classical boxplot (1.5)
    
    # Make the fences
    count_q1_tuk <- quantile(test_outliers_tot$n_tot_out_tuk, 
                             probs = .25)
    count_q3_tuk <- quantile(test_outliers_tot$n_tot_out_tuk, 
                             probs = .75)
    count_iqr_tuk <- count_q3_tuk - count_q1_tuk
    count_u_tuk <- count_q3_tuk + count_iqr_tuk*1.5
    
    # Filter on the fences
    outliers_tot_tuk <- unique(unlist(c(
      test_outliers_tot %>%
        dplyr::filter(n_tot_out_tuk >= count_u_tuk) %>%
        dplyr::select(obs_idx)
    )))
  }
  ###################################################################
  
  
  ## Step 4: Classify the type of functional outlier:
  ###################################################################
  # Calculate the number of intervals with five Location extreme stats
  # for the magnitude indicator, and the number of intervals with at
  # least one extreme Variation statistic. Then take the average of the
  # indicator to find the proportion of intervals over which each
  # criteria was met
  location_variation_counts <- test_outliers %>%
    dplyr::filter(obs_idx %in% outliers_tot_tuk) %>%
    group_by(obs_idx, n_intervals, interval) %>%
    dplyr::mutate(
      stat_type = dplyr::case_when(
        summary_stat %in% c("max_value", "mean_value",
                            "med_value", "min_value",
                            "obs_auc") ~ "Location_Stats",
        TRUE ~ "Variation_Stats"
      )) %>%
    group_by(obs_idx, n_intervals, interval, stat_type) %>%
    dplyr::summarize(n_out_tuk = sum(out_tuk)) %>%
    pivot_wider(names_from = stat_type, 
                values_from = n_out_tuk) %>%
    mutate(magnitude_ind = (Location_Stats > 2), 
           shape_ind = (Variation_Stats > 1)) 
  
  test_outliers_class <- location_variation_counts %>%
    dplyr::ungroup() %>%
    dplyr::group_by(obs_idx) %>%
    summarize(prop_magnitude_ind = mean(magnitude_ind), 
              prop_shape_ind = mean(shape_ind)) %>%
    # pivot_longer(c(prop_magnitude_ind, prop_shape_ind), 
    #              names_to = "stat_group", 
    #              values_to = "proportion_intervals") %>%
    dplyr::mutate(
      out_type_pred = dplyr::case_when(
        prop_magnitude_ind >= 0.33 & prop_shape_ind >= 0.20 ~ "Both",
        prop_magnitude_ind >= 0.33 ~ "Magnitude",
        prop_shape_ind >= 0.20 ~ "Shape", 
        TRUE ~ "Shape"
      ))
  ###################################################################
  
  
  ## Step 5: Return back the results
  ###################################################################
  outliers_class <- test_outliers_class %>%
    dplyr::select(obs_idx, out_type_pred) %>%
    dplyr::rename("class" = "out_type_pred")
  
  return(list(outliers_found = outliers_tot_tuk,
              outliers_class = outliers_class, 
              outliers_stats = location_variation_counts))
  ###################################################################
}
#####################################################################
# 
# test_data_obj <- simulation_model7()
# test_data <- test_data_obj$data
# 
# pod_test <- pod_fda(test_data, cutoff = "classical1.5")
