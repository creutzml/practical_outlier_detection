#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### November 27th, 2022                                           ###
###                                                               ###
###   Correcting function `tvdmss`, because current function      ###
### trips an error whenever there is a "positive" MSS outlier.    ###
#####################################################################


tvdmss_mc <- function(dts, 
                      emp_factor_mss = 1.5, 
                      emp_factor_tvd = 1.5, 
                      central_region_tvd = 0.5) {
  
  # Calculate the modified shape similarity index
  depths_mss <- total_variation_depth(dts = dts) 
  
  # Save the results
  tvd <- tvd_old <- depths_mss$tvd
  mss <- depths_mss$mss
  dta_dim <- dim(dts)
  n_curves <- dta_dim[1]
  n_points <- dta_dim[2]
  index <- (1:n_curves)
  n_central_obs <- ceiling(n_curves/2)
  
  # Find the (small) outlying MSS values
  shape_boxstats <- boxplot(mss, range = emp_factor_mss, plot = F)
  shape_outliers <- NULL
  
  # Check if any exist, and remove the shape outliers if they do
  if (length(shape_boxstats$out) != 0) { 
    if (any(shape_boxstats$out < mean(mss))) {
      shape_outliers <- which(
        mss %in% shape_boxstats$out[shape_boxstats$out < mean(mss)]
      )
      dts <- dts[-shape_outliers, ]
      tvd <- tvd[-shape_outliers]
      index <- index[-shape_outliers]
    }}
  
  # Next find the magnitude outliers
  magnitude_outliers <- NULL
  outliers <- functional_boxplot(
    dts, 
    depth_values = tvd, 
    emp_factor = emp_factor_tvd, 
    central_region = central_region_tvd * n_curves/nrow(dts)
  )$outliers
  
  # Save the list of outliers if any are found
  if (length(outliers) != 0) 
    magnitude_outliers = index[outliers]
  
  # Return back a list of all objects
  return(list(outliers = sort(c(magnitude_outliers, shape_outliers)), 
              shape_outliers = shape_outliers, 
              magnitude_outliers = magnitude_outliers, 
              tvd = tvd_old, mss = mss))
}