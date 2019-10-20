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
# 1 Drug Policy Program, Center for Research and Teaching in Economics (CIDE) - CONACT, 
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
# Med Decis Making. 2018;38(3):400 
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
c_l       <- 1 / 12                # cycle length (a month)
n_t       <- 20                    # number of years (20 years)
d_r       <- 0.03                  # discount rate
times     <- seq(0, n_t, c_l)  # the cycles in years

c_H <- 200
c_S <- 500
c_D <- 0

v_c <- c(c_H, c_S, c_D  )

u_H <- 0.75
u_S <- 0.30
u_D <- 0
v_u <- c(u_H, u_S, u_D  )
v_dw <- 1 / (1 + d_r) ^ (times)   # discount weight 


####################### Digitized  Data ################################


# 1. use the function digitise() to translate the digitised OS and PFS data into patient level information

digitise("OS_Examp.txt","OS_Examp_AtRisk.txt", nevent_inp= 52, km_output="KMdata_OS.txt",ipd_output="IPDdata_OS.txt")


##2.Create IPD and KM data for the OS curves 
IPD.OS  <- make.ipd(ipd_files = c("IPDdata_OS.txt") ,ctr = 1,var.labs = c("time","event","arm"))

##3.Create IPD and KM data for the PFS curves


##########################    Analysis #########################

# 5.create KM plots for the OS and PFS using survfit() and plot the results 

# 6.fit a weibull  curve to the OS data and plot using the flexsurvreg() package



########################Partitioned Survival model ################

# 7.fit all parametric models to the PFS data and extract the AIC/BIC. 

# 8.Select the one with the most appropriate fit

# 9.Repeat for OS

# 10. construct a partitioned survival model out of the fitted models

###################################################################

#calculate total cost per cycle





