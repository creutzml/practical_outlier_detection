#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### November 15th, 2022                                           ###
### Project: Plot Creation                                        ###
###                                                               ###
###   Generating example figures for presentation/dissertation.   ###
#####################################################################

# Libraries
library(tidyverse)
library(ggplot2)
library(fda)
library(fdaoutlier)
library(ggthemes)
library(gridExtra)


## Ojo data models
model_list <- vector("list", 9)
plot_colors <- c("lightgray", "red")
for (i in 1:9) {
  sim_mod_fct <- match.fun(paste0("simulation_model", i))
  model_temp <- sim_mod_fct()
  model_data_temp <- t(model_temp$data) %>%
    as.data.frame() %>%
    dplyr::mutate(t = 1:50) %>%
    tidyr::pivot_longer(-t, names_to = "Obs", values_to = "Value") %>%
    dplyr::mutate(Out = case_when(
      Obs %in% paste0("V", model_temp$true_outliers) ~ "Outlier",
      TRUE ~ "None"))
  
  model_list[[i]] <- ggplot() + 
    geom_line(aes(x = t, y = Value, group = Obs), 
              color = "lightgrey",
              data = model_data_temp %>% filter(Out == "None")) +
    geom_line(aes(x = t, y = Value, group = Obs), 
              color = "darkorange2",
              data = model_data_temp %>% filter(Out == "Outlier")) +
    guides(color = "none") +
    labs(title = paste0("Simulation Model ", i)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
}  

do.call("grid.arrange", c(model_list, ncol = 3))


  