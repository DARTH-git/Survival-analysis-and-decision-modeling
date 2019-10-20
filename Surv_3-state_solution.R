#####################################################################################
##########            Simple 3-state Partitioned Survival model in R ################
#####################################################################################

# Developed by the Decision Analysis in R for Technologies in Health (DARTH) workgroup
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (2)	
# M.G. Myriam Hunink, MD, PhD (3,4)
# Hawre J. Jalal, MD, PhD (5) 
# Eline M. Krijkamp, MSc (3)	
# Petros Pechlivanoglou, PhD (6) 

# In collaboration of: 		
# 1 Drug Policy Program, Center for Research and Teaching in Economics (CIDE) - CONACyT, 
#   Aguascalientes, Mexico
# 2 University of Minnesota School of Public Health, Minneapolis, MN, USA
# 3 Erasmus MC, Rotterdam, The Netherlands
# 4 Harvard T.H. Chan School of Public Health, Boston, USA
# 5 University of Pittsburgh Graduate School of Public Health, Pittsburgh, PA, USA
# 6 The Hospital for sick Children, Toronto and University of Toronto, Toronto ON, Canada

#####################################################################################
# Please cite our publications when using this code
# - Jalal H, Pechlivanoglou P, Krijkamp E, Alarid-Escudero F, Enns E, Hunink MG. 
# An Overview of R in Health Decision Sciences. Med Decis Making. 2017; 37(3): 735-746. 
# https://journals.sagepub.com/doi/abs/10.1177/0272989X16686559
# - Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400–22. 
# https://journals.sagepub.com/doi/abs/10.1177/0272989X18754513
# - Krijkamp EM, Alarid-Escudero F, Enns E, Pechlivanoglou P, Hunink MM, Jalal H. 
# A Multidimensional Array Representation of State-Transition Model Dynamics. 
# BioRxiv 670612 2019.https://www.biorxiv.org/content/10.1101/670612v1

#####################################################################################
# Copyright 2017, THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide. Copyright, 
# trademarks, trade names and any and all associated intellectual property are 
# exclusively owned by THE HOSPITAL FOR SICK CHILDREN and the collaborating 
# institutions. These materials may be used, reproduced, modified, distributed 
# and adapted with proper attribution.
#####################################################################################


#### 01 Load packages ####
if (!require('gems')) install.packages('gems'); library(gems)
if (!require('msm')) install.packages('msm'); library(msm)
if (!require('flexsurv')) install.packages('flexsurv'); library(flexsurv)
if (!require('dplyr')) install.packages('dplyr'); library(flexsurv)


#### 02 Load Functions ####
source("functions.R")

#### 03 Input Model Parameters ####

# number of states in the model
v_n       <- c("healthy", "sick", "dead") # state names
n_s       <- length(v_n)           # No of states 
n_i       <- 5000                  # number of simulations 
c_l       <- 1 / 12                # cycle length (a month)
n_t       <- 40                    # number of years (20 years)
times     <- seq(0, n_t, c_l)  # the cycles in years
set.seed(2009)                 # set the seed


# Create a transition probability matrix  with all transitions indicated and numbered

tmat <- matrix(NA, n_s, n_s, dimnames = list(v_n,v_n))
tmat["healthy", "sick"]  <- 1
tmat["healthy", "dead"]  <- 2
tmat["sick"   , "dead"]  <- 3

####################### Generate Data ################################
n_pat     <- 550                   # cohort size
n_years   <- 60                    # number of years 



generate <- gen_data(n_pat,n_years) # generates true, censored and OS/PFS data 
true_data <- generate$true_data
sim_data <- generate$sim_data
status   <- generate$status
OS_PFS_data <- generate$OS_PFS_data

head(true_data)
head(sim_data)
head(status)
head(OS_PFS_data)

##########################    Analysis #########################



## Showcasing the use of packages survival , flexsurv 

fit.KM <- survfit (Surv(time = OS_time, event= OS_status)~ 1 , data = OS_PFS_data, type ="fleming-harrington")
fit.weib <- flexsurvreg(Surv(time = OS_time, event= OS_status)~ 1 , data = OS_PFS_data, dist = "weibull")
plot(fit.KM,mark.time = T)
plot(fit.weib)

########################Partitioned Survival model ################

# fit all parametric models to the data and extract the AIC/BIC. 
# Select the one with the most appropriate fit
# Repeat for PFS and OS

fit.pfs  <- fit.fun(time ="PFS_time", status = "PFS_status", data = OS_PFS_data, times = times, extrapolate = T) 
fit.os   <- fit.fun(time = "OS_time", status = "OS_status" , data = OS_PFS_data, times = times,
                    extrapolate = T) 

best.pfs <- fit.pfs[["Weibull"]]  
best.os  <- fit.os [["Weibull"]]

# construct a partitioned survival model out of the fitted models
m_M_PSM <- partsurv(best.pfs,best.os, time = times)$trace


####################### MultiState modeling method 1 ##############

# The existing functions in R require the data in a long rather than a wide format
# convert the data in a way that flexsurv understands using the mstate package

data_long <- mstate::msprep(time = sim_data, status = status, trans = tmat )
data_long$trans <- as.factor(data_long$trans) # convert trans to a factor

data_long$from  <- case_when(data_long$from == 1 ~ "healthy", data_long$from ==2~"sick", data_long$from ==3~"dead")

data_long$to  <- case_when(data_long$to == 1 ~ "healthy", data_long$to ==2 ~ "sick", data_long$to ==3~"dead")


# fit all parametric multistate models simultaneously to the  data and extract the AIC/BIC. 
# Select the one with the lowest AIC

fits <- fit.mstate(time ="time", status = "status", trans, data = data_long, times = times, extrapolate = T )
best.fit <- fits[["Loglogistic"]]

# Construct a DES model out of the simultaneously fitted multistate model
DES_data <- sim.fmsm(best.fit, start = 1, t = n_years, trans = tmat, M = n_i)


m_M_DES  <- trace.DES(DES_data, n_i  = n_i , times = times, tmat = tmat)


####################### MultiState modeling method 2 ##############

# Multistate models can be fitted independently for each transition.  This is more flexible!
# Create subsets for each transition
data_HS <- subset(data_long, trans == 1)
data_HD <- subset(data_long, trans == 2)
data_SD <- subset(data_long, trans == 3)

# fit independent models for each transition and pick the one with the lowest AIC

fit_HS <- fit.fun(time = "time", status = "status", data = data_HS, times = times, extrapolate = T)
fit_HD <- fit.fun(time = "time", status = "status", data = data_HD, times = times, extrapolate = T)
fit_SD <- fit.fun(time = "time", status = "status", data = data_SD, times = times, extrapolate = T)

best.fit_HS <- fit_HS[["Gamma"]]
best.fit_HD <- fit_HD[["Weibull"]]
best.fit_SD <- fit_SD[["Lognormal"]]

# A microsimulation can be fitted instead of a DES. more computationally expensive
# but it provides more freeedom to the modeller. For the Microsimulation to be run
# we need transition probabilities per unit of time

# Extract transition probabilities from the best fitting models
p_HS <- flexsurvreg_prob(object = best.fit_HS, times = times)
p_HD <- flexsurvreg_prob(object = best.fit_HD, times = times)
p_SD <- flexsurvreg_prob(object = best.fit_SD, times = times)


# everyone starts in the "healthy" state and therefore has not spent time in "sick"
v_M_init  <- rep("healthy", times = n_i)   
v_Ts_init <- rep(0, n_i)  # a vector with the time of being sick at the start of the model 



Probs <- function(M_t, v_Ts, t) { 
  # Arguments:
  # M_t: health state occupied by at cycle t (character variable)
  # v_Ts: vector with the duration of being sick
  # t:     current cycle 
  # Returns: 
  #   transition probabilities for that cycle
  
  m_p_t           <- matrix(0, nrow = n_s, ncol = n_i)  # create matrix of state transition probabilities
  rownames(m_p_t) <-  v_n                               # give the state names to the rows
  
  # update m_p_t with the appropriate probabilities   
  m_p_t[, M_t == "healthy"] <- rbind(1 - p_HD[t] - p_HS[t], p_HS[t], p_HD[t])    # transition probabilities when healthy 
  m_p_t[, M_t == "sick"]    <- rbind(0, 1 - p_SD[v_Ts], p_SD[v_Ts])  # transition probabilities when sick 
  m_p_t[, M_t == "dead"]    <- c(0, 0, 1)                            # transition probabilities when dead     
  return(t(m_p_t))
}       


#### 06 Run Microsimulation ####
MicroSim <- function(n_i,  seed = 1) {
  # Arguments:  
  # n_i:     number of individuals
  # seed:    default is 1
  
  set.seed(seed) # set the seed
  
  # m_M is used to store the health state information over time for every individual
  times     <- seq(0, n_t, c_l)  # the cycles in years
  
  m_M  <-  matrix(nrow = n_i, ncol = length(times) + 1, 
                  dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                  paste("year", c(0,times), sep = " ")))  

  m_M[, 1] <- v_M_init          # initial health state for individual i
  v_Ts     <- v_Ts_init        # initialize time since illness onset for individual i
  
  # open a loop for time running cycles 1 to n_t 
  for (t in 1:length(times)) {
    v_p <- Probs(m_M[, t], v_Ts, t)  # calculate the transition probabilities for the cycle based on health state t
    m_M[, t + 1]  <- samplev(v_p, 1)       # sample the current health state and store that state in matrix m_M 
    
    v_Ts <- ifelse(m_M[, t + 1] == "sick", v_Ts + 1, 0) # update time since illness onset for t + 1 
    
    # Display simulation progress
    if(t/(n_t/10) == round(t/(n_t/10), 0)) { # display progress every 10%
      cat('\r', paste(t/n_t * 100, "% done", sep = " "))
    }
    
  } # close the loop for the time points 
  
  # store the results from the simulation in a list
  results <- list(m_M = m_M)   
  return(results)  # return the results
} # end of the MicroSim function  




# Run the simulation model
Micro_data <- MicroSim(n_i, seed = 1)

# create the microsimulation trace
m_M_Micro <- t(apply(Micro_data$m_M, 2, function(x) table(factor(x, levels = v_n, ordered = TRUE)))) 
m_M_Micro <- m_M_Micro / n_i                                  # calculate the proportion of individuals 
colnames(m_M_Micro) <- v_n                                    # name the rows of the matrix
rownames(m_M_Micro) <- paste("Cycle", c(0,times), sep = " ")       # name the columns of the matrix




# Calculate trace for the real data
m_M_data <- gems::transitionProbabilities(generate$cohort, times = times)@probabilities


# visually compare all methods
matplot(times,m_M_data, type='l', lty = 1, col = 1, ylab= "proportion of cohort", xlab = "Time",
        main = "Trace comparisons",xlim=c(0,25))
matlines(times, m_M_DES,  col = 2, lty = 1)
matlines(times, m_M_Micro[-1,], col = 3, lty = 1)
matlines(times, m_M_PSM,   col = 4, lty = 1)
legend("right", c("True Data", "DES","Microsim",  "PSM"),
        col = 1:4, lty = rep(1,4), bty= "n")

