##################################### data generation ############################

v_n2   <- c("S1","H","S2","D")  # the model states names

tmat <- matrix(NA, n_s, n_s, dimnames = list(v_n2,v_n2))
tmat["S1", "H"]  <- 1
tmat["S1", "S2"] <- 2
tmat["S1", "D"]  <- 3
tmat["S2", "D"]  <- 4

n_pat   <- 400
n_years <- 30


# specification of hazard functions to generate data from
hazardf <- gems::generateHazardMatrix(n_s)
colnames(hazardf@list.matrix) <- 
  rownames(hazardf@list.matrix) <- v_n2

# specifying the transition hazard from healthy -> sick
hazardf[["S1","H"]] <- function (t, r1, r2){
  hweibull(t,r1, r2)
}

hazardf[["S1","D"]] <- function (t, r1, r2){
  hweibull(t,r1, r2)
}

hazardf[["S1","S2"]] <- function (t, r1, r2){
  hweibull(t,r1, r2)
}
hazardf[["S2","D"]] <- function (t, r1, r2){
  hweibull(t,r1, r2)
}




# list of parameters for the hazard functions defined above
mu        <- gems::generateParameterMatrix(hazardf) 
rownames(mu@list.matrix) <- 
  colnames(mu@list.matrix) <- v_n2

mu[["S1", "H"]]  <- list(0.8,  2) 
mu[["S1", "D"]]  <- list(1.2,  40)  
mu[["S1", "S2"]] <- list(1,   10) 
mu[["S2", "D"]]  <- list(1, 20.8) 




# simulate the cohort
cohort <- gems::simulateCohort(
  transitionFunctions = hazardf,
  parameters = mu,
  cohortSize = n_pat,
  to = n_years)

# extract the simulated true data 
true_data <- cohort@time.to.state
colnames(true_data) <- v_n2

true_status = (is.na(true_data)!=1)*1


true_data$D[(is.na(true_data$D) &is.na(true_data$S2))] <- apply(cbind(true_data[(is.na(true_data$D)&is.na(true_data$S2)),-1],n_years),1,min,na.rm=T)

true_data$D[is.na(true_data$D)] <- n_years

true_data$S2[is.na(true_data$S2)] <- apply(cbind(true_data[is.na(true_data$S2),-1],n_years),1,min,na.rm=T)


true_data$H[is.na(true_data$H)] <- apply(cbind(true_data[is.na(true_data$H),-1],n_years),1,min,na.rm=T)


data_long <- mstate::msprep(true_data,true_status,trans = tmat )
data_long$trans <- as.factor(data_long$trans) # convert trans to a factor

data_long$from  <-case_when(data_long$from == 1~"S1", data_long$from ==2~"H", data_long$from ==3~"S2", data_long $from ==4 ~"D")

data_long$to  <-case_when(data_long$to == 1~"S1", data_long$to ==2~"H", data_long$to ==3~"S2", data_long$to ==4 ~"D")



write.csv(data_long, "data_long_Sicker.csv")