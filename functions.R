## Function to fit multiple functional forms to survival data

fit.fun <- function(time, status, data = data , add = FALSE, extrapolate = FALSE, times)  
{
  #Extact the right data columns 
  data$time   <-   data[,   time]  
  data$status <-   data[, status]  

    if (extrapolate == TRUE)  {
    plot.times <- max(times)
  } else if  (extrapolate == FALSE) {
    plot.times <- max(data$time)
  }
  
  # Progression free survival  
  KM.fit     <-     survfit(Surv(time, status) ~ 1, data = data)                   # fit Kaplan-Meier curve 
  fit.llogis <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "llogis" ) # fit model with loglogistic distribution
  fit.weib   <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "weibull") # fit model with Weibull distribution
  fit.lnorm  <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "lnorm"  ) # fit model with lognormal distribution
  fit.gamma  <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "gamma"  ) # fit model with gamma distribution 
  fit.exp    <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "exp"    ) # fit model with exponential distribution
  fit.gengamma  <- flexsurvreg(Surv(time, status) ~ 1, data = data, dist = "gengamma"  ) # fit model with gamma distribution  
  
  
  # extarapolate all models beyond the KM curve
  if(add){ lines(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0, plot.times), conf.int= F)}
  if(!add){ plot(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0, plot.times), conf.int= F, mark.time= T)}
  lines(fit.llogis,   t = times, col = 2, ci = F)
  lines(fit.weib,     t = times, col = 3, ci = F)
  lines(fit.lnorm,    t = times, col = 4, ci = F)
  lines(fit.gamma,    t = times, col = 5, ci = F)
  lines(fit.gengamma,    t = times, col = 6, ci = F)
  lines(fit.exp,      t = times, col = 7, ci = F)
  legend("topright", cex = 0.7, c("Kaplan-Meier", "Loglogistic", "Weibull", "Lognormal", "Gamma","GenGamma", "Exponential"), col = 1:7, lty = rep(1, 7), bty="n")
  
  # compare AIC values
  AIC <- c(    Loglogistic = AIC(fit.llogis),                                         
               Weibull     = AIC(fit.weib), 
               Lognormal   = AIC(fit.lnorm), 
               Gamma       = AIC(fit.gamma),
               GenGamma       = AIC(fit.gengamma),
               Exponentail = AIC(fit.exp))
  AIC= round(AIC,3)
  
  # compare BIC values
  BIC <- c(    Loglogistic = BIC(fit.llogis),                                         
               Weibull     = BIC(fit.weib), 
               Lognormal   = BIC(fit.lnorm), 
               Gamma       = BIC(fit.gamma),
               GenGamma    = BIC(fit.gengamma),
               Exponential = BIC(fit.exp))
  
  BIC <- round(BIC,3)
  
  res <- list(Loglogistic = fit.llogis,
              Weibull     = fit.weib,
              Lognormal   = fit.lnorm, 
              Gamma       = fit.gamma,
              GenGamma       = fit.gengamma,
              Exponential = fit.exp, 
              AIC         = AIC,
              BIC         = BIC)
  res
}


fit.mstate <- function(time, status, trans,  data = data , add = FALSE, extrapolate = FALSE, times)  
{
  data$time  <- data[, time  ]
  data$tatus <- data[, status]

    if (extrapolate == TRUE)  {
    plot.times <- max(times)
  } else if  (extrapolate == FALSE) {
    plot.times <- max(data$time)
  }
  
  # Progression free survival  
  KM.fit     <-     survfit(Surv(time, status) ~ trans , data = data)                   # fit Kaplan-Meier curve 
  fit.llogis <- flexsurvreg(Surv(time, status) ~ trans + shape(trans), data = data, dist = "llogis" ) # fit model with loglogistic distribution
  fit.weib   <- flexsurvreg(Surv(time, status) ~ trans + shape(trans), data = data, dist = "weibull") # fit model with Weibull distribution
  fit.lnorm  <- flexsurvreg(Surv(time, status) ~ trans + sdlog(trans), data = data, dist = "lnorm"  ) # fit model with lognormal distribution
  fit.gamma  <- flexsurvreg(Surv(time, status) ~ trans + shape(trans), data = data, dist = "gamma"  ) # fit model with gamma distribution 
  fit.gengamma  <- flexsurvreg(Surv(time, status) ~ trans + Q(trans) + sigma(trans), data = data, dist = "gengamma"  ) # fit model with gamma distribution 
  
  
  # extarapolate all models beyond the KM curve
  if(add){ lines(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0, plot.times), conf.int= F)}
  if(!add){ plot(KM.fit, ylab = "Survival Probability", xlab = "Time", ylim = c(0,1), xlim = c(0, plot.times), conf.int= F)}
  lines(fit.llogis,   t = times, col = 2, ci = F)
  lines(fit.weib,     t = times, col = 3, ci = F)
  lines(fit.lnorm,    t = times, col = 4, ci = F)
  lines(fit.gamma,    t = times, col = 5, ci = F)
  lines(fit.gengamma,    t = times, col = 6, ci = F)
  legend("topright", cex = 0.7, c("Kaplan-Meier", "Loglogistic", "Weibull", "Lognormal", "Gen.Gamma"), col = 1:5, lty = rep(1, 5), bty="n")
  
  # compare AIC values
  AIC <- c(    Loglogistic = AIC(fit.llogis),                                         
               Weibull     = AIC(fit.weib), 
               Lognormal   = AIC(fit.lnorm), 
               Gamma       = AIC(fit.gamma),
               GenGamma    = AIC(fit.gengamma))
  AIC= round(AIC,3)
  
  # compare BIC values
  BIC <- c(    Loglogistic = BIC(fit.llogis),                                         
               Weibull     = BIC(fit.weib), 
               Lognormal   = BIC(fit.lnorm), 
               Gamma       = BIC(fit.gamma),
               GenGamma    = BIC(fit.gengamma))
  
  BIC <- round(BIC,3)
  
  res <- list(Loglogistic = fit.llogis,
              Weibull     = fit.weib,
              Lognormal   = fit.lnorm, 
              Gamma       = fit.gamma,
              GenGamma    = fit.gengamma,
              AIC         = AIC,
              BIC         = BIC)
  res
}



trace.DES = function(msm_sim = des_sim, tmat, n_i, times )
{
  # Restructure the data to extract markov trace
  data.mstate.sim <- data.frame(cbind(matrix(t(msm_sim$st), ncol=1),
                                      matrix(t(msm_sim$t) , ncol=1)))
  colnames(data.mstate.sim) <- c("state","time")
  data.mstate.sim$subject <- rep(1:n_i, each = ncol(msm_sim$st))
  
  data.mstate.sim = na.omit(data.mstate.sim)
  data.mstate.sim = data.mstate.sim[!duplicated(data.mstate.sim), ] # remove duplicate entries in the dataset
  
  # create transition intensitiy matrix with initial values based on the structure of tmat
  twoway7.q               <- tmat
  twoway7.q[!is.na(tmat)] <- 0.5
  twoway7.q[is.na(tmat)]  <- 0
  # fit msm model only so that we can extract the prevalence (i.e. trace) thrrough the prevalence.msm function
  
  fit.msm.sim <- msm(state ~ time,subject = subject, data = data.mstate.sim, qmatrix = twoway7.q, 
                     exacttimes = T, use.deriv = TRUE, analyticp = FALSE, fixedpars = TRUE, hessian = F)
 
  M.tr.des <- prevalence.msm(fit.msm.sim, times = times) # Markov trace when DES model is used
  
  
  return(M.tr.des[[3]]/100)
}



partsurv <- function(fit.pfs, fit.os, time = times){
  # Input
  # fit.pfs: flexsurv obj fitting pfs
  # fit.os: flexsurv obj fitting os
  # title:
  # time = numeric vector of time to estimate probabilities
  # output:
  #  res a list w/ one entry of a data frame w/ probabilities associated w/ stable ,prog and dead.
  
  pfs.surv <- summary(fit.pfs, t = time, ci = F)[[1]]$est
  os.surv  <- summary(fit.os,  t = time, ci = F)[[1]]$est
  sick                 <- os.surv - pfs.surv      # estimate the probability of remaining in the progressed state
  sick[sick < 0]       <- 0                       # in cases where the probability is negative replace with zero
  healthy               <- pfs.surv                # probability of remaining stable
  dead                 <- 1 - os.surv             # probability of being dead
  trace <- cbind(healthy, sick, dead)
  res   <- list(trace = trace)

  return(res)
}

flexsurvreg_prob <- function(object, newparams = NULL, times){

  if(is.null(newparams) == T ){
  params <- object$res[,1]  
  params <- as.matrix(t(params))
  }else {
    params <- newparams 
    params <- as.matrix(params)
    }

  if (ncol(params)== 1){
  surv <- object$dfns$p(times, params[,1], lower.tail = F)
  }else if (ncol(params)== 2){
    surv <- object$dfns$p(times,params[,1],params[,2], lower.tail = F)
   }else if (ncol(params)== 2){
     surv <- object$dfns$p(times,params[,1],params[,2],params[,3], lower.tail = F)
   } else{
     surv <- object$dfns$p(times,params[,1],params[,2],params[,3], lower.tail = F)
   }
   
  t.p <- 1- surv[-1]/(surv[-length(surv)])
return(t.p = t.p)
  }
# 
# flexsurvreg_prob <- function(object, newdata,t, cycle ){
#   # input:
#   # object: flexsurv object will be using
#   # newdata: a dataframe that includes the coefficients used in obejct or F for no covariates
#   # t: time at each transition
#   # cycle: cycle length
#   # This uses the main part of the flexsurv summary function
#   
#   # output:
#   # vector of length = length(t) probability for each t and newdata 
#   if(! is.data.frame(newdata)){
#     start  = t - cycle 
#     x <- object
#     dat <- x$data
#     fn <-  function(t,start,...) {
#       ret <- 1-(1 - x$dfns$p(t,...))/(1 - x$dfns$p(start,...))
#       ret[t<start] <- 1 # prob[t<start] was previously 0
#       ret
#     }
#     fncall <- list(t, start)
#     beta <- if (x$ncovs == 0) {0}
#     X <- as.matrix(0, nrow = 1, ncol = max(x$ncoveffs, 1))
#     
#     dlist <- x$dlist
#     ret <- vector(nrow(X), mode = "list")
#     
#     basepars.mat <- flexsurv:::add.covs(x, x$res.t[dlist$pars, "est"], 
#                                         beta, X[1, , drop = FALSE], transform = FALSE)
#     basepars <- as.list(as.data.frame(basepars.mat))
#     fncall[dlist$pars] <- basepars
#     y <- do.call(fn, fncall)   
#     
#     
#   } else {
#     
#     
#     # Save flexsurv object as x and saving coefficients as Xraw
#     x <- object
#     Xraw <-  model.frame(x)[,unique(attr(model.frame(x),"covnames.orig")),drop=FALSE]
#     
#     # Creating model matrix with the coefficients from newdata and extracting covariate names
#     X <- flexsurv:::form.model.matrix(object, as.data.frame(newdata))
#     
#     # Creating start which is equal to t-sycle
#     start <- t - cycle
#     #fn is the function that will be using to generate estimates
#     # Pr = 1 - (1-S(t+1))/(1-S(t))
#     
#     fn <- function(t,start,...) {
#       ret <- 1-(1 - x$dfns$p(t,...))/(1 - x$dfns$p(start,...))
#       ret[t<start] <- 1 # prob[t<start] was previously 0
#       ret
#     }
#     
#     fn <- flexsurv:::expand.summfn.args(fn)
#     
#     # fncall is a list for each argument that we will be calling
#     fncall <- list(t,start)
#     
#     # name of beta
#     beta <- if (x$ncovs==0) 0 else x$res[x$covpars,"est"]
#     dlist <- x$dlist
#     # returning a entry for each one
#     ret <- vector(nrow(X), mode="list")
#     covnames <- rownames(X)
#     names(ret) <- covnames
#     
#     # what this is doing
#     
#     # x, pars, beta, X, transform=FALSE)
#     transform= F
#     nres <- nrow(X)
#     pars <- x$res.t[dlist$pars,"est"]
#     pars <- matrix(pars, nrow=nres, ncol=length(pars), byrow=TRUE)
#     beta <- matrix(beta, nrow=1)
#     
#     
#     for (j in seq(along=x$dlist$pars)){
#       covinds <- x$mx[[x$dlist$pars[j]]]
#       if (length(covinds) > 0){
#         pars[,j] <- pars[,j] + beta[,covinds] %*% t(X[,covinds,drop=FALSE])
#       }
#       
#       if (!transform)
#         pars[,j] <- x$dlist$inv.transforms[[j]](pars[,j])
#     }
#     
#     basepars.mat <- flexsurv:::add.covs(x, x$res.t[dlist$pars,"est"], beta, X[,,drop=FALSE], transform=FALSE)
#     basepars <- as.list(as.data.frame(basepars.mat))
#     fncall[dlist$pars] <- basepars
#     
#     y <- do.call(fn, fncall)
#   }
#   return(y)
# }
# 

gen_data <- function(n_pat, n_years)
{
  # specification of hazard functions to generate data from
  hazardf <- gems::generateHazardMatrix(n_s)
  colnames(hazardf@list.matrix) <- 
    rownames(hazardf@list.matrix) <- v_n
  
  # specifying the transition hazard from healthy -> sick
  hazardf[["healthy","sick"]] <- function (t, r1, r2){
    hweibull(t,r1, r2)
  }
  
  # specifying the transition hazard from healthy -> dead 
  hazardf[["healthy","dead"]] <- function (t, r1, r2){
    flexsurv::hgompertz(t,r1, r2)
  }
  
  # specifying the transition hazard from sick -> dead 
  hazardf[["sick","dead"]] <- function (t, r1, r2){
    hweibull(t,r1, r2)
  }
  
  
  # list of parameters for the hazard functions defined above
  mu        <- gems::generateParameterMatrix(hazardf) 
  rownames(mu@list.matrix) <- 
    colnames(mu@list.matrix) <- v_n
  
  mu[["healthy", "sick"]] <- list(1.5, 6)   #  the Weibull parameters for H -> S
  mu[["healthy", "dead"]] <- list(0.25, 0.08)  # the Gompertz params for H -> D
  mu[["sick",    "dead"]] <- list(0.5,4)     #  the Weibull parameters for S -> D
  
  
  
  # simulate the cohort
  cohort <- gems::simulateCohort(
    transitionFunctions = hazardf,
    parameters = mu,
    cohortSize = n_pat,
    to = n_years)
  
  # extract the simulated true data 
  true_data <- cohort@time.to.state
  colnames(true_data) <- v_n
  
  true_data$dead[is.na(true_data$dead)] <- n_years
  true_data$sick[is.na(true_data$sick)] <- true_data$dead[is.na(true_data$sick)]
  
  
  # create a status variable that will capture the transition events
  true_status         <- matrix(NA, nrow = n_pat, ncol = n_s, dimnames = list(1:n_pat,v_n))
  true_status         <- as.data.frame(true_status)
  true_status$healthy <- ifelse(is.na(true_data$healthy),0,1)
  true_status$dead    <- ifelse(true_data$dead == n_years, 0, 1)
  true_status$sick    <- ifelse(true_data$dead == true_data$sick, 0, 1)
  
  
  censtime <- runif(n = n_pat, 0, n_years)
  
  censored_sick <- ifelse(censtime      <= true_data$sick |
                          true_data$sick >  5, 1, 0)
  censored_dead <- ifelse(censtime <= true_data$dead|
                    true_data$dead >5, 1, 0)

  sim_data <- true_data
  
  sim_data$sick[censored_sick == 1] <-  censtime[censored_sick == 1]
  sim_data$sick[sim_data$sick >5 ]  <-  5
  
  sim_data$dead[censored_dead == 1] <-  censtime[censored_dead == 1]
  sim_data$dead[sim_data$dead >5] <-  5
  
  status <- true_status
  
  status$sick[censored_sick == 1] = 0
  status$dead[censored_dead == 1] = 0
  
  # Usually trials report OS/PFS outcomes so we will recreate those
  
  OS_PFS_data <- data.frame(row.names = 1:n_pat)
  
  OS_PFS_data$PFS_time        <- apply(sim_data[, c("sick","dead")], 1, min) 
  OS_PFS_data$PFS_status      <- ifelse(status$dead == 1 | status$sick == 1, 1, 0 )
  
  OS_PFS_data$OS_time         <- sim_data$dead
  OS_PFS_data$OS_status       <- status$dead 
  list(cohort = cohort, true_data = true_data, true_status = true_status, 
        sim_data =  sim_data,      status = status, OS_PFS_data = OS_PFS_data)
}

samplev <- function(m.Probs, m) {
  # Arguments
  # m.Probs: matrix with probabilities (n.i * n.s)
  # m:       number of states than need to be sampled per individual  
  # Return
  # ran:    n.i x m matrix filled with sampled health state(s) per individual
  
  d <- dim(m.Probs)  # dimensions of the matrix filled with the multinomical probabilities for the health states 
  n <- d[1]          # first dimension - number of rows (number of individuals to sample for)
  k <- d[2]          # second dimension - number of columns (number of health states considered)
  lev <- dimnames(m.Probs)[[2]]  # extract the names of the health states considered for sampling
  if (!length(lev))  # in case names for the health states are missing, use numbers to specify the health states
    lev <- 1:k       # create a sequence from 1:k (number of health states considered)
  # create a matrix 
  ran <- matrix(lev[1], ncol = m, nrow = n) # create the matrix ran, filled with the first health state of the levels 
  U <- t(m.Probs)    # transposed m.Probs matrix n.i x n.s --> n.s x n.i 
  
  for(i in 2:k) {    # start loop, from the 2nd health states
    U[i, ] <- U[i, ] + U[i - 1, ] # start summing the probabilities of the different health states per individual 
  }
  if (any((U[k, ] - 1) > 1e-05))  # sum of all probs per individual - 1 should be 0 (use 1e-05 for rounding issues), else print the error statement
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) {   # start loop of the state that needs to be sampled (m)
    un <- rep(runif(n), rep(k, n))       # sample from a uniform distribution of length n*k
    ran[, j] <- lev[1 + colSums(un > U)] # store the health state at the jth column of the U matrix
  }
  ran # return the new health state per individual n.i x m
} # close the function 






#plot health state trace
plot_m_TR <- function(m_M) {
  # plot the distribution of the population across health states over time (trace)
  # count the number of individuals in each health state at each cycle
  m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_n, ordered = TRUE)))) 
  m_TR <- m_TR / n_i                                       # calculate the proportion of individuals 
  colnames(m_TR) <- v_n                                    # name the rows of the matrix
  rownames(m_TR) <- paste("Cycle", 0:n_t, sep = " ")       # name the columns of the matrix
  # Plot trace of first health state
  matplot(m_TR, type = "l", main = "Health state trace", col= 1:n_s,
          ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
  legend("topright", v_n, col = 1:n_s,    # add a legend to current plot
         lty = rep(1, 3), bty = "n", cex = 0.65)
  
}



digitise<-function (surv_inp, nrisk_inp, nevent_inp=NA,km_output = "KMdata.txt", ipd_output = "IPDdata.txt") 
{
  #################################################
  ###From survHE package
  ###Author: Patricia Guyot and Gianluca Baio
  #################################################
  working.dir <- dirname(surv_inp)
  if(is.na(nevent_inp)){tot.events <- "NA"}
  if(!is.na(nevent_inp)){tot.events<-nevent_inp}
  arm.id <- 1
  digizeit <- read.table(surv_inp, header = TRUE, row.names = NULL)
  t.S <- digizeit[, 2]
  S <- digizeit[, 3]
  pub.risk <- read.table(nrisk_inp, header = TRUE, row.names = NULL)
  pub.risk <- pub.risk[pub.risk[, 4] > 0, ]
  if (!(pub.risk[1, 3] == 1)) {
    pub.risk[1, 3] <- 1
  }
  t.risk <- pub.risk[, 2]
  lower <- pub.risk[, 3]
  upper <- pub.risk[, 4]
  n.risk <- pub.risk[, 5]
  n.int <- length(n.risk)
  n.t <- upper[n.int]
  arm <- rep(arm.id, n.risk[1])
  n.censor <- rep(0, (n.int - 1))
  n.hat <- rep(n.risk[1] + 1, n.t)
  cen <- d <- rep(0, n.t)
  KM.hat <- rep(1, n.t)
  last.i <- rep(1, n.int)
  sumdL <- 0
  if (n.int > 1) {
    for (i in 1:(n.int - 1)) {
      n.censor[i] <- round(n.risk[i] * S[lower[i + 1]]/S[lower[i]] - 
                             n.risk[i + 1])
      while ((n.hat[lower[i + 1]] > n.risk[i + 1]) || ((n.hat[lower[i + 
                                                                    1]] < n.risk[i + 1]) && (n.censor[i] > 0))) {
        if (n.censor[i] <= 0) {
          cen[lower[i]:upper[i]] <- 0
          n.censor[i] <- 0
        }
        if (n.censor[i] > 0) {
          cen.t <- rep(0, n.censor[i])
          for (j in 1:n.censor[i]) {
            cen.t[j] <- t.S[lower[i]] + j * (t.S[lower[(i + 
                                                          1)]] - t.S[lower[i]])/(n.censor[i] + 1)
          }
          cen[lower[i]:upper[i]] <- hist(cen.t, breaks = t.S[lower[i]:lower[(i + 
                                                                               1)]], plot = F)$counts
        }
        n.hat[lower[i]] <- n.risk[i]
        last <- last.i[i]
        for (k in lower[i]:upper[i]) {
          if (i == 1 & k == lower[i]) {
            d[k] <- 0
            KM.hat[k] <- 1
          }
          else {
            d[k] <- round(n.hat[k] * (1 - (S[k]/KM.hat[last])))
            KM.hat[k] <- KM.hat[last] * (1 - (d[k]/n.hat[k]))
          }
          n.hat[k + 1] <- n.hat[k] - d[k] - cen[k]
          if (d[k] != 0) 
            last <- k
        }
        n.censor[i] <- n.censor[i] + (n.hat[lower[i + 
                                                    1]] - n.risk[i + 1])
      }
      if (n.hat[lower[i + 1]] < n.risk[i + 1]) 
        n.risk[i + 1] <- n.hat[lower[i + 1]]
      last.i[(i + 1)] <- last
    }
  }
  if (n.int > 1) {
    n.censor[n.int] <- min(round(sum(n.censor[1:(n.int -1)]) * 
                                   (t.S[upper[n.int]] - t.S[lower[n.int]])/
                                   (t.S[upper[(n.int - 1)]] - t.S[lower[1]])), n.risk[n.int])
  }
  if (n.int == 1) {
    n.censor[n.int] <- 0
  }
  if (n.censor[n.int] <= 0) {
    cen[lower[n.int]:(upper[n.int] - 1)] <- 0
    n.censor[n.int] <- 0
  }
  if (n.censor[n.int] > 0) {
    cen.t <- rep(0, n.censor[n.int])
    for (j in 1:n.censor[n.int]) {
      cen.t[j] <- t.S[lower[n.int]] + j * (t.S[upper[n.int]] - 
                                             t.S[lower[n.int]])/(n.censor[n.int] + 1)
    }
    cen[lower[n.int]:(upper[n.int] - 1)] <- hist(cen.t, breaks = t.S[lower[n.int]:upper[n.int]], 
                                                 plot = F)$counts
  }
  n.hat[lower[n.int]] <- n.risk[n.int]
  last <- last.i[n.int]
  for (k in lower[n.int]:upper[n.int]) {
    if (KM.hat[last] != 0) {
      d[k] <- round(n.hat[k] * (1 - (S[k]/KM.hat[last])))
    }
    else {
      d[k] <- 0
    }
    KM.hat[k] <- KM.hat[last] * (1 - (d[k]/n.hat[k]))
    n.hat[k + 1] <- n.hat[k] - d[k] - cen[k]
    if (n.hat[k + 1] < 0) {
      n.hat[k + 1] <- 0
      cen[k] <- n.hat[k] - d[k]
    }
    if (d[k] != 0) 
      last <- k
  }
  if (tot.events != "NA") {
    if (n.int > 1) {
      sumdL <- sum(d[1:upper[(n.int - 1)]])
      if (sumdL >= tot.events) {
        d[lower[n.int]:upper[n.int]] <- rep(0, (upper[n.int] - 
                                                  lower[n.int] + 1))
        cen[lower[n.int]:(upper[n.int] - 1)] <- rep(0, 
                                                    (upper[n.int] - lower[n.int]))
        n.hat[(lower[n.int] + 1):(upper[n.int] + 1)] <- rep(n.risk[n.int], 
                                                            (upper[n.int] + 1 - lower[n.int]))
      }
    }
    if ((sumdL < tot.events) || (n.int == 1)) {
      sumd <- sum(d[1:upper[n.int]])
      while ((sumd > tot.events) || ((sumd < tot.events) && 
                                     (n.censor[n.int] > 0))) {
        n.censor[n.int] <- n.censor[n.int] + (sumd - 
                                                tot.events)
        if (n.censor[n.int] <= 0) {
          cen[lower[n.int]:(upper[n.int] - 1)] <- 0
          n.censor[n.int] <- 0
        }
        if (n.censor[n.int] > 0) {
          cen.t <- rep(0, n.censor[n.int])
          for (j in 1:n.censor[n.int]) {
            cen.t[j] <- t.S[lower[n.int]] + j * (t.S[upper[n.int]] - 
                                                   t.S[lower[n.int]])/(n.censor[n.int] + 1)
          }
          cen[lower[n.int]:(upper[n.int] - 1)] <- hist(cen.t, 
                                                       breaks = t.S[lower[n.int]:upper[n.int]], 
                                                       plot = F)$counts
        }
        n.hat[lower[n.int]] <- n.risk[n.int]
        last <- last.i[n.int]
        for (k in lower[n.int]:upper[n.int]) {
          d[k] <- round(n.hat[k] * (1 - (S[k]/KM.hat[last])))
          KM.hat[k] <- KM.hat[last] * (1 - (d[k]/n.hat[k]))
          if (k != upper[n.int]) {
            n.hat[k + 1] <- n.hat[k] - d[k] - cen[k]
            if (n.hat[k + 1] < 0) {
              n.hat[k + 1] <- 0
              cen[k] <- n.hat[k] - d[k]
            }
          }
          if (d[k] != 0) 
            last <- k
        }
        sumd <- sum(d[1:upper[n.int]])
      }
    }
  }
  KMdata <- data.frame(time = t.S, n.risk = n.hat[1:n.t], n.event = d, 
                       n.censored = cen)
  write.table(KMdata, km_output, sep = "\t", row.names = FALSE, 
              col.names = TRUE)
  t.IPD <- rep(t.S[n.t], n.risk[1])
  event.IPD <- rep(0, n.risk[1])
  k <- 1
  for (j in 1:n.t) {
    if (d[j] != 0) {
      t.IPD[k:(k + d[j] - 1)] <- rep(t.S[j], d[j])
      event.IPD[k:(k + d[j] - 1)] <- rep(1, d[j])
      k <- k + d[j]
    }
  }
  for (j in 1:(n.t - 1)) {
    if (cen[j] != 0) {
      t.IPD[k:(k + cen[j] - 1)] <- rep(((t.S[j] + t.S[j + 
                                                        1])/2), cen[j])
      event.IPD[k:(k + cen[j] - 1)] <- rep(0, cen[j])
      k <- k + cen[j]
    }
  }
  IPD <- data.frame(time = t.IPD, event = event.IPD, arm)
  write.table(IPD, ipd_output, sep = "\t", row.names = FALSE, 
              col.names = TRUE)
  if (dirname(km_output) == ".") {
    cat("\n")
    cat(paste0("Kaplan Meier data written to file: ", working.dir, 
               km_output))
  }
  else {
    cat("\n")
    cat(paste0("Kaplan Meier data written to file: ", km_output))
  }
  if (dirname(ipd_output) == ".") {
    cat("\n")
    cat(paste0("IPD data written to file: ", working.dir, 
               ipd_output))
    cat("\n")
  }
  else {
    cat("\n")
    cat(paste0("IPD data written to file: ", ipd_output))
    cat("\n")
  }
}
# 
# trace.DES = function(msm.sim = des.sim)
# {
#   # Restructure the data to extract markov trace
#   data.mstate.sim <- data.frame(cbind(matrix(t(msm.sim$st), ncol=1),
#                                       matrix(t(msm.sim$t) , ncol=1)))
#   colnames(data.mstate.sim) <- c("state","time")
#   data.mstate.sim$subject <- rep(1:n.i, each = ncol(msm.sim$st))
#   
#   data.mstate.sim = na.omit(data.mstate.sim)
#   data.mstate.sim = data.mstate.sim[!duplicated(data.mstate.sim), ] # remove duplicate entries in the dataset
#   
#   # create transition intensitiy matrix with initial values based on the structure of tmat
#   twoway7.q               <- tmat
#   twoway7.q[!is.na(tmat)] <- 0.5
#   twoway7.q[is.na(tmat)]  <- 0
#   # fit msm model only so that we can extract the prevalence (i.e. trace) thrrough the prevalence.msm function
#   
#   fit.msm.sim <- msm(state ~ time,subject = subject, data = data.mstate.sim, qmatrix = twoway7.q, 
#                      exacttimes = T, use.deriv = TRUE, analyticp = FALSE, fixedpars = TRUE, hessian = F)
#   
#   M.tr.des <- prevalence.msm(fit.msm.sim, times = times) # Markov trace when DES model is used
#   
#   
#   return(M.tr.des[[3]]/100)
# }
# 
# 

make.ipd <- function(ipd_files,ctr=1,var.labs=c("time","event","arm")) {
  ## Piles in the simulated IPD resulting from running digitise for more than one treatment arm  
  ## ipd_files = a list including the names of the IPD files created as output of digitise
  ## ctr = the index of the file associated with the control arm (default, the first file).
  ##       This will be coded as 0
  ## var.labs = a vector of labels for the column of the resulting data matrix. NB these
  ##            should match the arguments to the formula specified for fit.models. The
  ##            user can specify values. These should be 3 elements (TIME, EVENT, ARM)
  
  #################################################
  ###From survHE package
  ###Author: Gianluca Baio
  #################################################
  
  # Identifies the number of arms (= number of IPD files)
  n_arms <- length(ipd_files)
  index <- 1:n_arms
  active <- index[-ctr]
  data <- read.table(ipd_files[[ctr]],header=TRUE,row.names=NULL)
  data[,"arm"] <- 0 # sets the value of "arm" to 0, for the control group
  arm.ind <- 1
  for (i in active) {
    tmp <- read.table(ipd_files[[index[i]]],header=TRUE,row.names=NULL)
    tmp[,"arm"] <- arm.ind
    data <- rbind(data,tmp)
    arm.ind <- arm.ind+1
  }
  colnames(data) <- var.labs
  return(data)
}

plot_compare<-function(M.tr,M.tr.dig){
  par(mfrow=c(2,2))
  plot(M.tr[order(M.tr[,1]),1],M.tr.dig[order(M.tr[,1]),1],type="l",lwd=2,
       xlab="Partitioned Survival Using IPD",ylab="Estimated Partitioned Survival",
       main="Sick Free")
  abline(a=0,b=1,col="red")
  plot(M.tr[order(M.tr[,2]),2],M.tr.dig[order(M.tr[,2]),2],type="l",lwd=2,
       xlab="Partitioned Survival Using IPD",ylab="Estimated Partitioned Survival",
       main="Sick")
  abline(a=0,b=1,col="red")
  plot(M.tr[order(M.tr[,3]),3],M.tr.dig[order(M.tr[,3]),3],type="l",lwd=2,
       xlab="Partitioned Survival Using IPD",ylab="Estimated Partitioned Survival",
       main="Dead")
  abline(a=0,b=1,col="red")
  
  matplot(M.tr.dig,type='l', lwd=2,main="Estimated Partitioned Survival")
  par(mfrow=c(1,1))
}


plot_PSA<-function(psa.pfs_scale.dig,psa.pfs_scale,psa.pfs_shape.dig,psa.pfs_shape){
  par(mfrow=c(2,1))
  hist(psa.pfs_scale.dig,main="Density of PSA distribution for Scale",xlab="Scale",freq=FALSE,col=rgb(0,0,1))
  hist(psa.pfs_scale,add=TRUE,freq=FALSE,col=rgb(1,0,0,0.5))
  
  hist(psa.pfs_shape.dig,main="Density of PSA distribution for Shape",xlab="Shape",freq=FALSE,col=rgb(0,0,1))
  hist(psa.pfs_shape,add=TRUE,freq=FALSE,col=rgb(1,0,0,0.5))
  par(mfrow=c(1,1))
}




#################### Functions used in the microsimulation Sick Sicker example #######################


Probs <- function(M_it, df_X, v_Ts1, v_Ts2, t) { 
  # Arguments:
  # M_it: health state occupied by individual i at cycle t (character variable)
  # v_Ts: time an individual is sick
  # t:     current cycle 
  # Returns: 
  #   transition probabilities for that cycle
  
  m_p_it           <- matrix(0, nrow = n_s, ncol = n_i)  # create matrix of state transition probabilities
  rownames(m_p_it) <-  v_n                               # give the state names to the rows
  
  # lookup baseline probability and rate of dying based on individual characteristics
  p_HD_all <- inner_join(df_X, p_mort, by = c("Age"))
  p_HD     <- p_HD_all[M_it == "H", "p_HD"]
  
  # update the v_p with the appropriate probabilities   
  m_p_it[, M_it == "H"]  <- rbind(1 - p_HS1 - p_HD, p_HS1, 0, p_HD)                              # transition probabilities when healthy
  m_p_it[, M_it == "S1"] <- rbind(p_S1H[v_Ts1], 1 - p_S1H[v_Ts1] - p_S1S2[v_Ts1] - p_S1D[v_Ts1], p_S1S2[v_Ts1], p_S1D[v_Ts1])  # transition probabilities when sick
  m_p_it[, M_it == "S2"] <- rbind(0, 0, 1 - p_S2D[v_Ts2], p_S2D[v_Ts2])                                            # transition probabilities when sicker
  m_p_it[, M_it == "D"]  <- c(0, 0, 0, 1)                                                        # transition probabilities when dead   
  return(t(m_p_it))
}       

#### 05.2 Cost function ####
# The Costs function estimates the costs at every cycle.

Costs <- function (M_it, Trt = FALSE) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Trt:  is the individual being treated? (default is FALSE) 
  
  c_it <- 0                                  # by default the cost for everyone is zero 
  c_it[M_it == "H"]  <- c_H                  # update the cost if healthy
  c_it[M_it == "S1"] <- c_S1 + c_Trt * Trt   # update the cost if sick conditional on treatment
  c_it[M_it == "S2"] <- c_S2 + c_Trt * Trt   # update the cost if sicker conditional on treatment
  c_it[M_it == "D"]  <- c_D                  # update the cost if dead
  
  return(c_it)        		                   # return the costs
}

#### 05.3 Health outcome function ####
# The Effs function to update the utilities at every cycle.

Effs <- function (M_it, df_X, Trt = FALSE, cl = 1) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # df_X: individual characteristics including Age, Sex and the effect modifier of the treatment effect
  # Trt:  is the individual treated? (default is FALSE) 
  # cl:   cycle length (default is 1)
  
  u_it <- 0                                        # by default the utility for everyone is zero
  u_it[M_it == "H"]  <- u_H                        # update the utility if healthy
  u_it[M_it == "S1" & Trt == FALSE] <- u_S1        # update the utility if sick
  u_it[M_it == "S1" & Trt == TRUE]  <- u_Trt * df_X$x[M_it == "S1"]  # update the utility if sick but on treatment (adjust for individual effect modifier) 
  u_it[M_it == "S2"] <- u_S2                       # update the utility if sicker
  u_it[M_it == "D"]  <- u_D                        # update the utility if dead
  
  QALYs <-  u_it * cl            # calculate the QALYs during cycle t
  return(QALYs)                  # return the QALYs
}


#### 06 Run Microsimulation ####
MicroSim <- function(n_i, df_X , Trt = FALSE, seed = 1) {
  # Arguments:  
  # n_i:     number of individuals
  # df_X     data frame with individual data 
  ## Age      age of the individuals
  ## Sex      sex of the individuals 
  ## x        effect modifier  
  # Trt:     is this the individual receiving treatment? (default is FALSE)
  # seed:    default is 1
  
  set.seed(seed) # set the seed
  
  n_s <- length(v_n) # the number of health states
  
  # create three matrices called m_M, m_C and m_E
  # number of rows is equal to the n_i, the number of columns is equal to n_t  (the initial state and all the n_t cycles)
  # m_M is used to store the health state information over time for every individual
  # m_C is used to store the costs information over time for every individual
  # m_E is used to store the effects information over time for every individual
  
  m_M <- m_C <- m_E <- m_Ts <-  matrix(nrow = n_i, ncol = n_t + 1, 
                                       dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                                       paste("cycle", 0:n_t, sep = " ")))  
  
  m_M [, 1] <- v_M_init    # initial health state at cycle 0 for individual i
  v_Ts1      <- v_Ts1_init   # initialize time since illnes onset for individual i
  v_Ts2      <- v_Ts2_init   # initialize time since illnes onset for individual i
  
  m_C[, 1]  <- Costs(m_M[, 1], Trt)         # calculate costs per individual during cycle 0
  m_E[, 1]  <- Effs (m_M[, 1], df_X, Trt)   # calculate QALYs per individual during cycle 0
  
  # open a loop for time running cycles 1 to n_t 
  for (t in 1:n_t) {
    v_p <- Probs(m_M[, t], df_X, v_Ts1,v_Ts2, t)             # calculate the transition probabilities for the cycle based on  health state t
    m_M[, t + 1]  <- samplev(v_p, 1)                  # sample the current health state and store that state in matrix m_M 
    m_C[, t + 1]  <- Costs(m_M[, t + 1], Trt)         # calculate costs per individual during cycle t + 1
    m_E[, t + 1]  <- Effs(m_M[, t + 1], df_X, Trt)    # calculate QALYs per individual during cycle t + 1
    
    v_Ts1 <- if_else(m_M[, t + 1] == "S1", v_Ts1 + 1, 0)
    v_Ts2 <- if_else(m_M[, t + 1] == "S2", v_Ts2 + 1, 0) 
    df_X$Age[m_M[, t + 1] != "D"]  <- df_X$Age[m_M[, t + 1] != "D"] + 1
    
    # Display simulation progress
    if(t/(n_t/10) == round(t/(n_t/10), 0)) { # display progress every 10%
      cat('\r', paste(t/n_t * 100, "% done", sep = " "))
    }
    
  } # close the loop for the time points 
  
  # calculate  
  tc <- m_C %*% v_dwc    # total (discounted) cost per individual
  te <- m_E %*% v_dwe    # total (discounted) QALYs per individual 
  tc_hat <- mean(tc)     # average (discounted) cost 
  te_hat <- mean(te)     # average (discounted) QALYs
  
  # store the results from the simulation in a list
  results <- list(m_M = m_M, m_C = m_C, m_E = m_E, tc = tc , te = te, tc_hat = tc_hat, te_hat = te_hat)   
  return(results)  # return the results
} # end of the MicroSim function  


