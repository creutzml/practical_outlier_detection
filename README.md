# practical_outlier_detection
Code developed to implement Practical Outlier Detection (POD), a method for identifying and classifying the type of functional outlier in a univariate, random sample of functional data.  

## General Structure of Repository
The contents of this repository are organized into three main folders: R Code, Data, and Figures. A short description of each folder is given below.

- **R Code**: folder that contains all R scripts for implementing POD, creating a random sample of functional data with outliers present, reproducing simulations, and reproducing the case study.
- **Data**: the saved simulation .RData files are provided here for easy reproducibility of the simulation figures and tables. Additionally, a copy of the World Population Growth study can also be found in this folder. 
- **Figures**: copies of the plots and figures presented in the manuscript.
\end{itemize}

## R Code: Details
The primary R script of interest is `practical_outlier.R`, which contains the function to implement my functional outlier detection method, Practical Outlier Detection (POD). 

**An Example**:
```
# Load the Ojo et al. (2021) library for easy data generation
library(fdaoutlier)

# Create an example data set
test_data_obj <- simulation_model7()
test_data <- test_data_obj$data

# Implement POD
pod_test <- pod_fda(test_data, cutoff = "classical1.5")
```

The remaining scripts in the folder allow one to recreate the simulation results and case study results found in the manuscript. Specifically, a description of each script is given below:
- `simulaton_model_b.R`: function to simulate a random sample of functional data with magnitude, shape, and two types of both magnitude and shape functional outliers. Adapted from Ojo et al.'s (2021) function `simulation_model1()` in the package `fdaoutlier`.
- `simulation_model_b_only.R`: function to simulate a random sample of functional data with two types of both magnitude and shape functional outliers. Adapted from Ojo et al.'s (2021) function `simulation_model1()` in the package `fdaoutlier`.
- `tvdmss_mc.R`: 
- `plotting_script.R`: allows one to recreate Figures 1, 6, and 7 in the manuscript
- `tvdmss_mc.R`: adapts the code originally created by Huang and Sun (2019) to correct for an error that occurred when there is a "positive" MSS outlier.
- `run_outlier_sim_for_type.R`: allows one to run a single parameter combination (at a time) of the simulation described in Section 3.3
- `results_outlier_sim_for_type.R`: summarizes and reproduces the simulation results of Section 3.3
- `world_population.R`: reproduces the results of the case study presented in Section 4


