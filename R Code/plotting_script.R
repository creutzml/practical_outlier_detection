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

## Berkeley growth study
# data
data("growth")
growth_data <- growth$hgtf[,1:10] %>%
  data.frame() %>%
  mutate(age = growth$age) %>%
  pivot_longer(cols = -age, names_to = "girl", values_to = "height")

ggplot() + 
  geom_line(data = growth_data, 
            aes(x = age, y = height, color = girl)) +
  geom_point(data = growth_data, 
             aes(x = age, y = height, color = girl)) +
  geom_vline(xintercept = seq(3.12, 18, length.out = 8), linetype = "dashed") +
  guides(color = "none") +
  labs(x = "Age (years)", 
       y = "Height (cm)", 
       caption = paste0("Heights collected on 10 females at 31 ages",
                        ", from the Berkeley Growth Study")) +
  theme_bw() +
  theme(text=element_text(size=24))

## Gait cycle data
data("gait")
hip_angle <- gait[,,1] %>%
  data.frame() %>%
  mutate(time = as.numeric(rownames(gait[,,1]))) %>%
  pivot_longer(-time, names_to = "boy", values_to = "angle") %>%
  mutate(location = "Hip")
knee_angle <- gait[,,2] %>%
  data.frame() %>%
  mutate(time = as.numeric(rownames(gait[,,2]))) %>%
  pivot_longer(-time, names_to = "boy", values_to = "angle") %>%
  mutate(location = "Knee")

gait_data <- rbind.data.frame(hip_angle, knee_angle) %>%
  as.data.frame()

ggplot() +
  geom_line(data = gait_data, 
            aes(x = time, y = angle, color = boy)) +
  geom_point(data = gait_data, 
             aes(x = time, y = angle, color = boy)) +
  facet_grid(rows = vars(location)) +
  guides(color = "none") +
  labs(x = "Time (proportion of gait cycle)", 
       y = "Angle (degrees)",
       title = "Multivariate Random Sample of Gait Cycle Data",
       subtitle = "n = 39, p = 2, T = 20") +
  theme_bw() +
  theme(text=element_text(size=24))
  