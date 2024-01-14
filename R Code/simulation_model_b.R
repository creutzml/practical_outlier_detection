#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### June 20th, 2023                                               ###
###                                                               ###
###   Adapts the sample generation function from Ojo et al.       ###
### (2021) to produce a random sample of data with multiple types ###
### of functional outliers present, including outliers that are a ###
### half and half combination of magnitude and shape outlier.     ###
#####################################################################


simulation_modelB <- function (
  n = 100, p = 50, outlier_rate = 0.05, mu = 4, q = 8, 
  kprob = 0.5, a = 0.1, b = 0.9, l = 0.05, cov_alpha = 1, 
  cov_beta = 1, cov_nu = 1, deterministic = TRUE, seed = NULL, 
  plot = F, plot_title = "Simulation Model Mix", title_cex = 1.5, 
  show_legend = T, ylabel = "", xlabel = "gridpoints"
) {
  
  tt <- seq(0, 1, length.out = p)
  covfun <- fdaoutlier:::covfunexp(tt, cov_alpha, cov_beta, cov_nu)
  muu <- mu * tt
  L <- chol(covfun)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  e <- matrix(rnorm(n * p), nrow = p, ncol = n)
  y <- muu + t(L) %*% e
  dtt <- fdaoutlier:::determine(deterministic, n, outlier_rate)
  n_outliers <- dtt$n_outliers
  true_outliers <- dtt$true_outliers
  
  # Split into shape and shift
  outlier_type <- cut(x = 1:length(true_outliers), 
                      breaks = 4, 
                      labels = c("Shift", 
                                 "Shape", 
                                 "Half and Half",
                                 "Combined"))
  
  true_shift <- true_outliers[outlier_type == "Shift"]
  n_shift <- length(true_shift)
  true_shape <- true_outliers[outlier_type == "Shape"]
  n_shape <- length(true_shape)
  true_hh <- true_outliers[outlier_type == "Half and Half"]
  n_hh <- length(true_hh)
  true_combined <- true_outliers[outlier_type == "Combined"]
  n_combined <- length(true_combined)
  
  # Save the indices to dtt
  dtt[["shift_outliers"]] <- true_shift
  dtt[["shape_outliers"]] <- true_shape
  dtt[["half_and_half_outliers"]] <- true_hh
  dtt[["combined_outliers"]] <- true_combined
  
  if (n_outliers > 0) {
    # Create matrix of random errors for the outliers
    e <- matrix(rnorm(p * n_outliers), nrow = p)
    
    # Split for shift versus shape outliers
    e_shift <- e[,1:n_shift]
    e_shape <- e[,(n_shift + 1):(n_shift + n_shape)]
    e_hh <- e[,(n_shift + n_shape + 1):(n_shift + n_shape + n_hh)]
    e_combined <- e[,(n_shift + n_shape + n_hh + 1):
                      (n_shift + n_shape + n_hh + n_combined)]
    
    # Make the shift outliers
    qcoeffk <- rbinom(n_shift, 1, kprob)
    qcoeffk[qcoeffk == 0] <- -1
    qcoeffk <- qcoeffk * q
    y[, true_shift] <- (muu + t(L) %*% e_shift) + 
      rep(qcoeffk, rep(p, n_shift))
    
    # Make the shape outliers
    theta <- rep(runif(n_shape, a, b), each = p)
    indicator <- 2 * sin(4 * pi * (tt + theta))
    y[, true_shape] <- (muu + t(L) %*% e_shape) + indicator

    # Combined outliers: magnitude and shift across the
    # sampling domain
    qcoeffk <- rbinom(n_combined, 1, kprob)
    qcoeffk[qcoeffk == 0] <- -1
    qcoeffk <- qcoeffk * q
    theta <- rep(runif(n_combined, a, b), each = p)
    indicator <- 2 * sin(4 * pi * (tt + theta))
    y[, true_combined] <- 
      (muu + t(L) %*% e_combined) + 
      rep(qcoeffk, rep(p, n_combined)) + 
      indicator[1:(p*n_combined)]
    
    # Half and half outliers: magnitude for part, then shape
    # for second half
    qcoeffk <- rbinom(n_hh, 1, kprob)
    qcoeffk[qcoeffk == 0] <- -1
    qcoeffk <- qcoeffk * q
    theta <- rep(runif(n_hh, a, b), each = p)
    indicator <- 2 * sin(4 * pi * (tt + theta))
    y[, true_hh] <- 
      (muu + t(L) %*% e_hh) + 
      c(rep(qcoeffk, rep(p, n_hh))[1:(p%/%2)], 
        rep(0, length((p %/% 2 + 1):p))) + 
      c(rep(0, length(1:(p%/%2))), 
        indicator[(p%/%2 + 1):p])
  }
  if (plot) {
    fdaoutlier:::plot_dtt(y, tt, p, true_outliers, show_legend, 
                          plot_title, title_cex, ylabel, xlabel)
  }
  return(list(data = t(y), true_outliers = dtt))
}

