#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### March 30th, 2023                                              ###
###                                                               ###
###   Functional Outlier Detection case study data example. The   ###
### data set is provided by `fdaoutlier` and is already pre-      ###
### processed. World Population Growth is given for 105 countries ###
### for the years 1950-2010 (T = 61). The countries were chosen   ###
### by including those with more than 1 million and less than     ###
### 15 million for the year of 1980.                              ###
#####################################################################

## Libraries
library(fda)
library(fdaoutlier)
library(tidyverse)
library(ggplot2)
library(ggrepel)


## Helper function and main method
#####################################################################
# Directories
dir_path <- file.path(here::here(), "R Code")
data_path <- file.path(here::here(), "Data")

# Source the functions
source(file.path(dir_path, "practical_outlier.R"))
source(file.path(dir_path, "tvdmss_mc.R"))
#####################################################################



## Load the data in
#####################################################################
# View(world_population)
# matplot(t(world_population), type = "l")

# Grab common values
n_obs <- nrow(world_population)
n_sp <- ncol(world_population)
#####################################################################

### Note: since my methodologies require a user specified cutoff, I 
### will use the results of other methods to determine the cutoff. 
### With that in mind, most choose anywhere between 20% and 25% of
### the countries as "outliers." So, I will set the cutoffs to be
### 20%
# thresh <- 0.9
# 
# 
# ## Run easy outlier method (user specified)
# #################################################################
# easy_test_user <- easy_out_fda(
#   test_data = world_population, 
#   cutoff = thresh
# )
# 
# # Grab the results
# easy_outliers_user <- easy_test_user$outliers_tot
# easy_outliers_user_summary <- easy_test_user$outliers_summary
# 
# 
# # What countries?
# easy_outliers_user_idx <- unlist(lapply(
#   1:length(easy_outliers_user),
#   FUN = function(x) {
#     as.numeric(substr(
#       easy_outliers_user[x], 
#       nchar(easy_outliers_user[x]) - 2,
#       nchar(easy_outliers_user[x])
#     ))
# }))
# 
# # Names of outlying countries
# easy_outliers_user_countries <- 
#   rownames(world_population)[easy_outliers_user_idx]
# #################################################################

## Run easy outlier method with Tukey 1.5*IQR rule for final
#################################################################
# Run with these interval sizes instead
easy_test_c1.5 <- pod_fda(
  test_data = world_population, 
  cutoff = "classical1.5"
)

# Grab the results
easy_outliers_c1.5 <- easy_test_c1.5$outliers_found
easy_outliers_c1.5_summary <- easy_test_c1.5$outliers_stats
easy_outliers_c1.5_shape <- easy_test_c1.5$outliers_class %>%
  dplyr::filter(class %in% c("Shape", "Both")) %>%
  dplyr::pull(obs_idx)
easy_outliers_c1.5_shift <- easy_test_c1.5$outliers_class %>%
  dplyr::filter(class %in% c("Magnitude", "Both")) %>%
  dplyr::pull(obs_idx)


# What countries?
easy_outliers_c1.5_idx <- unlist(lapply(
  1:length(easy_outliers_c1.5),
  FUN = function(x) {
    as.numeric(substr(
      easy_outliers_c1.5[x], 
      nchar(easy_outliers_c1.5[x]) - 2,
      nchar(easy_outliers_c1.5[x])
    ))
  }))

# Names of outlying countries
easy_outliers_c1.5_countries <- 
  rownames(world_population)[easy_outliers_c1.5_idx]
#################################################################



## Outlier detection with MS-Plot
#################################################################
ms_test <- msplot(dts = world_population, 
                  return_mvdir = F,
                  plot = F)

# Grab the results
ms_outliers <- paste0("Obs", 
                      str_pad(ms_test$outliers, 3, "left", "0"))

# Country names
ms_outliers_countries <- rownames(world_population)[ms_test$outliers]
#################################################################


## Run TVD on the data
#################################################################
tvd_test <- tvdmss_mc(dts = world_population)

# Grab the results
tvd_outliers <- paste0("Obs", 
                       str_pad(tvd_test$outliers, 3, "left", "0"))
tvd_outliers_shape <- paste0(
  "Obs", str_pad(tvd_test$shape_outliers, 3, "left", "0")
)
tvd_outliers_shift <- c() # there are none

# Country names
tvd_outliers_countries <- rownames(world_population)[tvd_test$outliers]
#################################################################


## Run the MUOD method
#################################################################
muod_test_tan <- muod(dts = world_population, cut_method = "tangent")
muod_test_box <- muod(dts = world_population, cut_method = "boxplot")

# Grab the results
muod_outliers_tan <- paste0(
  "Obs", str_pad(unique(unlist(muod_test_tan$outliers)), 
                 3, "left", "0")
)

muod_outliers_box <- paste0(
  "Obs", str_pad(unique(unlist(muod_test_box$outliers)), 
                 3, "left", "0")
)

muod_outliers_tan_shift <- paste0(
  "Obs", str_pad(unique(unlist(muod_test_tan$outliers$magnitude)), 
                 3, "left", "0")
)

muod_outliers_box_shift <- paste0(
  "Obs", str_pad(unique(unlist(muod_test_box$outliers$magnitude)), 
                 3, "left", "0")
)

muod_outliers_tan_shape <- paste0(
  "Obs", str_pad(unique(c(muod_test_tan$outliers$shape, 
                          muod_test_tan$outliers$amplitude)),
                 3, "left", "0")
)
  
muod_outliers_box_shape <- paste0(
  "Obs", str_pad(unique(c(muod_test_box$outliers$shape, 
                          muod_test_box$outliers$amplitude)),
                 3, "left", "0")
)

muod_test_box$outliers$magnitude

# Country names
muod_outliers_tan_countries <- rownames(world_population)[
  unique(unlist(muod_test_tan$outliers))
]

muod_outliers_box_countries <- rownames(world_population)[
  unique(unlist(muod_test_box$outliers))
]
#################################################################


### Plot of results for POD
#################################################################
world_population_long <- world_population %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Country") %>%
  dplyr::mutate(Outlier = Country %in% easy_outliers_c1.5_countries) %>%
  pivot_longer(-c(Country, Outlier), 
               names_to = "Year", 
               values_to = "Population") %>%
  dplyr::mutate(Year = as.numeric(substr(Year, 2, 5))) %>%
  dplyr::mutate(`Outlier Type` = case_when(
    Outlier & Country == "Netherlands" ~ "Magnitude", 
    Outlier & Country == "Sudan" ~ "Both",
    Outlier & !(Country %in% c("Sudan", "Netherlands")) ~ "Shape",
    TRUE ~ "None"
  )) %>%
  dplyr::mutate(`Outlier Type` = factor(`Outlier Type`, 
                                        levels = c("Magnitude", 
                                                   "Shape",
                                                   "Both",
                                                   "None")))

## Plot of world population outliers by POD
# Make colorblind palette
# The palette with grey:
cbPalette <- c("#E69F00", "#56B4E9", "#CC79A7", "#999999")

## Recreate Figure 5 in Section 4
# Make the plot
ggplot() + 
  geom_line(aes(x = Year, 
                y = Population, 
                color = `Outlier Type`, 
                group = Country,
                linetype = `Outlier Type`),
            size = 0.5,
            alpha = 0.5,
            data = world_population_long %>% 
              dplyr::filter(`Outlier Type` == "None")) +
  geom_line(aes(x = Year, 
                y = Population, 
                color = `Outlier Type`,
                linetype = `Outlier Type`,
                group = Country), 
            size = 2,
            data = world_population_long %>% 
              dplyr::filter(`Outlier Type` != "None")) +
  geom_label_repel(aes(x = Year, 
                       y = Population, 
                       label = Country), 
                   data = world_population_long %>% 
                     dplyr::filter(Outlier, Year == 2010), 
                   size = 6) +
  scale_color_manual(values = cbPalette, 
                     breaks = c("Magnitude", "Shape", "Both", "None")) +
  scale_linetype_manual(values = c("F1", "twodash", "dashed", "solid"),
                        breaks = c("Magnitude", "Shape", "Both", "None")) +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.15, 0.75), 
        legend.background = element_rect(fill = "white", 
                                         color = "black"), 
        legend.key.size = unit(3, "cm"),
        legend.key.height = unit(1.5, "line"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  labs(y = "Population ('000)", 
       color = "Outlier Type", 
       linetype = "Outlier Type")
#################################################################



### Plot of results for TVD
#################################################################
# world_population_long <- world_population %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "Country") %>%
#   dplyr::mutate(Outlier = Country %in% tvd_outliers_countries) %>%
#   pivot_longer(-c(Country, Outlier), 
#                names_to = "Year", 
#                values_to = "Population") %>%
#   dplyr::mutate(Year = as.numeric(substr(Year, 2, 5))) %>%
#   dplyr::mutate(`Outlier Type` = case_when(
#     Outlier ~ "Shape",
#     TRUE ~ "None"
#   )) %>%
#   dplyr::mutate(`Outlier Type` = factor(`Outlier Type`, 
#                                         levels = c("Magnitude", 
#                                                    "Shape", 
#                                                    "None")))
# 
# ## Plot of world population outliers by POD
# # Make colorblind palette
# # The palette with grey:
# cbPalette <- c("#56B4E9", "#999999")
# # Make the plot
# ggplot() + 
#   geom_line(aes(x = Year, 
#                 y = Population, 
#                 color = `Outlier Type`, 
#                 group = Country), 
#             data = world_population_long) +
#   geom_label_repel(aes(x = Year, 
#                        y = Population, 
#                        label = Country), 
#                    data = world_population_long %>% 
#                      dplyr::filter(Outlier, Year == 2010), 
#                    size = 6) +
#   scale_color_manual(values = cbPalette) +
#   theme_bw(base_size = 20) +
#   labs(y = "Population ('000)")
#################################################################


## Comparison of methods:
#################################################################
world_comp <- data.frame(
  Method = c("POD", "MS-Plot", "TVD", "MUOD (box)", "MUOD (tan)"), 
  `n (per)` = c(paste0(length(easy_outliers_c1.5), " (", 
                     format(round(length(easy_outliers_c1.5)/n_obs, 3),
                            nsmall = 3), 
                     ")"), 
              paste0(length(ms_outliers), " (", 
                     format(round(length(ms_outliers)/n_obs, 3),
                            nsmall = 3),
                     ")"),
              paste0(length(tvd_outliers), " (", 
                     format(round(length(tvd_outliers)/n_obs, 3), 
                            nsmall = 3),
                     ")"), 
              paste0(length(muod_outliers_box), " (", 
                     format(round(length(muod_outliers_box)/n_obs, 3),
                            nsmall = 3), 
                     ")"), 
              paste0(length(muod_outliers_tan), " (", 
                     format(round(length(muod_outliers_tan)/n_obs, 3),
                            nsmall = 3), 
                     ")")),
  `n Magnitude (per)` = c(
    paste0(length(easy_outliers_c1.5_shift), " (", 
           format(round(length(easy_outliers_c1.5_shift)/n_obs, 3), 
                  nsmall = 3), 
           ")"),
    NA,
    paste0(length(tvd_outliers_shift), " (", 
           format(round(length(tvd_outliers_shift)/n_obs, 3), 
                  nsmall = 3),
           ")"), 
    paste0(length(muod_outliers_box_shift), " (", 
           format(round(length(muod_outliers_box_shift)/n_obs, 3),
                  nsmall = 3), 
           ")"), 
    paste0(length(muod_outliers_tan_shift), " (", 
           format(round(length(muod_outliers_tan_shift)/n_obs, 3),
                  nsmall = 3), 
           ")")),
  `n Shape (per)` = c(
    paste0(length(easy_outliers_c1.5_shape), " (", 
           format(round(length(easy_outliers_c1.5_shape)/n_obs, 3), 
                  nsmall = 3), 
           ")"),
    NA,
    paste0(length(tvd_outliers_shape), " (", 
           format(round(length(tvd_outliers_shape)/n_obs, 3), 
                  nsmall = 3),
           ")"), 
    paste0(length(muod_outliers_box_shape), " (", 
           format(round(length(muod_outliers_box_shape)/n_obs, 3),
                  nsmall = 3), 
           ")"), 
    paste0(length(muod_outliers_tan_shape), " (", 
           format(round(length(muod_outliers_tan_shape)/n_obs, 3),
                  nsmall = 3), 
           ")")), 
  POD = c(NA, 
          sum(ms_outliers %in% easy_outliers_c1.5), 
          sum(tvd_outliers %in% easy_outliers_c1.5),
          sum(muod_outliers_box %in% easy_outliers_c1.5),
          sum(muod_outliers_tan %in% easy_outliers_c1.5)),
  `MS-Plot` = c(sum(ms_outliers %in% easy_outliers_c1.5), 
                NA,
                sum(tvd_outliers %in% ms_outliers),
                sum(muod_outliers_box %in% ms_outliers),
                sum(muod_outliers_tan %in% ms_outliers)),
  TVD = c(sum(tvd_outliers %in% easy_outliers_c1.5), 
          sum(tvd_outliers %in% ms_outliers),
          NA,
          sum(muod_outliers_box %in% tvd_outliers),
          sum(muod_outliers_tan %in% tvd_outliers)), 
  `MUOD (box)` = c(sum(muod_outliers_box %in% easy_outliers_c1.5), 
                   sum(muod_outliers_box %in% ms_outliers),
                   sum(muod_outliers_box %in% tvd_outliers),
                   NA,
                   sum(muod_outliers_tan %in% muod_outliers_box)), 
  `MUOD (tan)` = c(sum(muod_outliers_tan %in% easy_outliers_c1.5), 
                   sum(muod_outliers_tan %in% ms_outliers),
                   sum(muod_outliers_tan %in% tvd_outliers),
                   sum(muod_outliers_tan %in% muod_outliers_box),
                   NA)
)

## Reproduces Table 3 in Section 4
# Produce latex code
kableExtra::kbl(world_comp, 
                booktabs = TRUE, 
                format = "latex") %>%
  kableExtra::kable_styling(latex_options = c("scale_down", 
                                              "hold_position"))

## Describing the similar results of all methods applied to the
## World Population Growth case study
# Calculate how many were agreed on by every method:
all_agreed <- intersect(
  intersect(
    intersect(
      intersect(easy_outliers_c1.5, ms_outliers), tvd_outliers), 
    muod_outliers_box),
  muod_outliers_tan
)
length(all_agreed)
all_agreed_countries <- lapply(
  1:length(all_agreed), 
  FUN = function(x) {
    temp_idx <- as.numeric(substr(all_agreed[x], 
                                  nchar(all_agreed[x]) - 2, 
                                  nchar(all_agreed[x])))
    temp_idx <- unlist(temp_idx)                       
    rownames(world_population)[temp_idx]
  }
)

do.call(cat, all_agreed_countries)

#################################################################
