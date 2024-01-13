simulation_modelB2 <- function (
  n = 100, p = 50, outlier_rate = 0.05, mu = 4, q = 8, 
  kprob = 0.5, a = 0.1, b = 0.9, l = 0.05, cov_alpha = 1, 
  cov_beta = 1, cov_nu = 1, deterministic = TRUE, seed = NULL, 
  plot = F, plot_title = "Simulation Model Combined", title_cex = 1.5, 
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
                      breaks = 2, 
                      labels = c("Half and Half",
                                 "Combined"))

  true_hh <- true_outliers[outlier_type == "Half and Half"]
  n_hh <- length(true_hh)
  true_combined <- true_outliers[outlier_type == "Combined"]
  n_combined <- length(true_combined)
  
  # Save the indices to dtt
  dtt[["half_and_half_outliers"]] <- true_hh
  dtt[["combined_outliers"]] <- true_combined
  
  if (n_outliers > 0) {
    # Create matrix of random errors for the outliers
    e <- matrix(rnorm(p * n_outliers), nrow = p)
    
    # Split for shift versus shape outliers
    e_hh <- e[,1:n_hh]
    e_combined <- e[,(n_hh + 1):(n_hh + n_combined)]

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

# simulation_modelB2(cov_beta = .9, plot = T)

