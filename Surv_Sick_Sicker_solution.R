
############  Microsimulation Sick-Sicker model ##########
# Includes: 
# individual characteristicss: age
# age dependent mortality probabilities 
# individual treatment effect modifyer 
################################################################################

# Developed by the Decision Analysis in R for Technologies in Health (DARTH) group
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (1)	
# M.G. Myriam Hunink, MD, PhD (2,3)
# Hawre J. Jalal, MD, PhD (4) 
# Eline M. Krijkamp, MSc (2)
# Petros Pechlivanoglou, PhD (5) 

# In collaboration of: 		
# 1 University of Minnesota School of Public Health, Minneapolis, MN, USA
# 2 Erasmus MC, Rotterdam, The Netherlands
# 3 Harvard T.H. Chan School of Public Health, Boston, USA
# 4 University of Pittsburgh Graduate School of Public Health, Pittsburgh, PA, USA
# 5 The Hospital for Sick Children, Toronto and University of Toronto, Toronto ON, Canada

################################################################################
# Please cite our publications when using this code:
# darthworkgroup.com
# Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400-422.
# See GitHub for more information or code updates
# https://github.com/DARTH-git/Microsimulation-tutorial
################################################################################
# Copyright 2017, 
# THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide.  
# Copyright, trademarks, trade names and any and all associated intellectual property
# are exclusively owned by THE HOSPITAL FOR SICK CHILDREN and the collaborating 
# institutions and may not be used, reproduced, modified, distributed or adapted 
# in any way without appropriate citation.
################################################################################

rm(list =ls()) # clear memory (removes all the variables from the workspace)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the folder where the file is saved 

#### 01 Load packages ####
library(dplyr) # load dplyr including the useful inner_join function
library(dampack) # for CEA, PSA, and visualization
library(flexsurv)

#### 02 Load Functions ####
source("functions.R")


#### 03 Input Model Parameters ####
set.seed(1)                    # set the seed  
v_n   <- c("H","S1","S2","D")  # the model states names
# Model structure 
n_t     <- 30                    # time horizon, 30 cycles
c_l     <- 1
n_i   <- 100000                # number of simulated individuals
n_s   <- length(v_n)           # the number of health states
d_r   <- 0.03                  # discount rate of 3% per cycle
v_dwe <- v_dwc <- 1 / ((1 + d_r) ^ (0:n_t))   # discount weight 
v_names_str <- c("no treatment", "treatment") # strategy names

# Event probabilities (per cycle)
## Annual transition probabilities
p_HS1   <- 0.15                # probability of becoming sick when healthy

## Annual probabilities of death
p_mort   <- read.csv("mortProb_age.csv")        # load age dependent probability
dist_Age <- read.csv("MyPopulation-AgeDistribution.csv") # load age distribution

# Cost inputs
c_H     <- 2000             # cost of one cycle in the healthy state
c_S1    <- 4000             # cost of one cycle in the sick state
c_S2    <- 15000            # cost of one cycle in the sicker state
c_D     <- 0                # cost of one cycle in the dead state
c_Trt   <- 12000            # cost of treatment (per cycle)

# Utility inputs
u_H     <- 1                # utility when healthy 
u_S1    <- 0.75             # utility when sick 
u_S2    <- 0.5              # utility when sicker
u_D     <- 0                # utility when dead
u_Trt   <- 0.95             # utility when sick(er) and being treated

times     <- seq(0, n_t, c_l)  # the cycles in years
set.seed(2019)

#### 04 Sample individual level characteristics ####
#### 04.1 Static characteristics ####
v_x      <- runif(n_i, min = 0.95, max = 1.05) # treatment effect modifier at baseline  
v_age0   <- sample(x = dist_Age$age, prob = dist_Age$prop, size = n_i, replace = TRUE) 
# sample from age distribution an initial age for every individual
df_X     <- data.frame(ID = 1:n_i, x = v_x, Age = v_age0)

########################### Survival analysis component  ##################



# load the Sicker data 
data_long <- read.csv("data_long_Sicker.csv",row.names = 1)
head(data_long)


# Multistate models can be fitted independently for each transition.  This is more flexible!

# Create subsets for each transition
data_S1H  <- subset(data_long, trans == 1)
data_S1S2 <- subset(data_long, trans == 2)
data_S1D  <- subset(data_long, trans == 3)
data_S2D  <- subset(data_long, trans == 4)


# fit independent models for each transition and pick the one with the lowest AIC

fit_S1H <- fit.fun(time ="time", status = "status", data = data_S1H, times = times, extrapolate = F)
fit_S1S2 <- fit.fun(time ="time", status = "status", data = data_S1S2, times = times, extrapolate = F)
fit_S1D <- fit.fun(time ="time", status = "status", data = data_S1D, times = times, extrapolate = F)
fit_S2D <- fit.fun(time ="time", status = "status", data = data_S2D, times = times, extrapolate = F)


best.fit_S1H <-  fit_S1H [[which.min(fit_S1H$AIC)]]
best.fit_S1S2 <- fit_S1S2[[which.min(fit_S1S2$AIC)]]
best.fit_S1D <-  fit_S1D [[which.min(fit_S1D$AIC)]]
best.fit_S2D <-  fit_S2D [[which.min(fit_S2D$AIC)]]


# Extract transition probabilities from the best fitting models
p_S1H  <- flexsurvreg_prob(object = best.fit_S1H,t = times)
p_S1S2 <- flexsurvreg_prob(object = best.fit_S1S2,t = times)
p_S1D  <- flexsurvreg_prob(object = best.fit_S1D,t = times)
p_S2D  <- flexsurvreg_prob(object = best.fit_S2D,t = times)






#### 04.2 Dynamic characteristics 
# Specify the initial health state of the individuals 
# everyone begins in the healthy state (in this example)
v_M_init  <- rep("H", n_i)       # a vector with the initial health state for all individuals 
v_Ts1_init <- rep(0, n_i)         # a vector with the time of being sick at the start of the model 
v_Ts2_init <- rep(0, n_i)         # a vector with the time of being sick at the start of the model 

#### 05 Define Simulation Functions ####
#### 05.1 Probability function ####
# The Probs function that updates the transition probabilities of every cycle is shown below.

# By specifying all the arguments in the `MicroSim()` the simulation can be started
# In this example the outcomes are of the simulation are stored in the variables `outcomes_no_tr` and `outcomes_trt`.


### Run the simulation for both no treatment and treatment options
outcomes_no_trt  <- MicroSim(n_i, df_X, Trt = FALSE, seed = 1)
outcomes_trt     <- MicroSim(n_i, df_X, Trt = TRUE, seed = 1)


#### 07 Visualize results ####
options(scipen = 999)

### No treatment
plot(density(outcomes_no_trt$tc), main = paste("Total cost per person"), xlab = "Cost ($)")
plot(density(outcomes_no_trt$te), main = paste("Total QALYs per person"), xlab = "QALYs")

plot_m_TR(outcomes_no_trt$m_M)    # health state trace

### Treatment
plot(density(outcomes_trt$tc), main = paste("Total cost per person"), xlab = "Cost ($)")
plot(density(outcomes_trt$te), main = paste("Total QALYs per person"), xlab = "QALYs")

plot_m_TR(outcomes_trt$m_M)    # health state trace



#### 08 Cost Effectiveness Analysis ####
# store the mean costs of each strategy in a new variable C (vector of costs)
v_C <- c(outcomes_no_trt$tc_hat, outcomes_trt$tc_hat)
# store the mean QALYs of each strategy in a new variable E (vector of effects)
v_E <- c(outcomes_no_trt$te_hat, outcomes_trt$te_hat)

# use dampack to calculate the ICER
calculate_icers(cost       = v_C,
                effect     = v_E,
                strategies = v_names_str)






