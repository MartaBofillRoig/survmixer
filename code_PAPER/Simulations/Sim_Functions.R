
################################################################
# SIMULATIONS MDA Research
# Estimating the achieved power
# Last update 2019-01-10 (Eivissa); 
# Marta Bofill and Guadalupe Gómez
################################################################

# install.packages("survRM2")
require(survRM2)
require(survival)
# install.packages("FHtest")
# library(FHtest)
# The FHtest package requires package Icens from the Bioconductor:
# source("http://bioconductor.org/biocLite.R")
# biocLite("Icens")

# functions for the simulation study
fun_sim <- function(n0,n1,p0,p1,bet0,bet1,lambda_r0,lambda_r1,lambda_nr0,lambda_nr1,lambda_cens,truncated=F,tau=NULL){
  n0= round(n0)
  n1= round(n1)
  
  # group1
  BE1 = ifelse(runif(n1)<p1,1,0)
  n_r1 = sum(BE1)

  TE_r1 = rweibull(n_r1,scale=lambda_r1,shape=bet1)
  TE_nr1 = rweibull(n1-n_r1,scale=lambda_nr1,shape=bet1) 
  CE1 = rexp(n1,1/lambda_cens) 
  
  db1 = data.frame(resp=c(rep(1,n_r1), rep(0,n1-n_r1)), time_obs=c(TE_r1,TE_nr1), cens_obs=c(CE1), treat=c(rep(1,n1)))
  
  if(truncated==F){
    db1$time = ifelse(db1$time_obs<=db1$cens_obs, db1$time_obs,db1$cens_obs)
    db1$status = ifelse(db1$time_obs<=db1$cens_obs, 1,0)
    db1$time_obs <- NULL
    db1$cens_obs <- NULL
  }else{
    db1$time = ifelse(db1$time_obs<=pmin(db1$cens_obs,tau), db1$time_obs, pmin(db1$cens_obs,tau))
    db1$status = ifelse(db1$time_obs<=pmin(db1$cens_obs,tau), 1,0)
    db1$time_obs <- NULL
    db1$cens_obs <- NULL
  }
  
  
  # group 0
  BE0 = ifelse(runif(n0)<p0,1,0)
  n_r0 = sum(BE0)
  TE_r0 = rweibull(n_r0,scale=lambda_r0,shape=bet0)
  TE_nr0 = rweibull(n0-n_r0,scale=lambda_nr0,shape=bet0) 
  CE0 = rexp(n0,1/lambda_cens) 
  
  db0 = data.frame(resp=c(rep(1,n_r0), rep(0,n0-n_r0)), time_obs=c(TE_r0,TE_nr0), cens_obs=c(CE0), treat=c(rep(0,n0)))
  
  if(truncated==F){
    db0$time = ifelse(db0$time_obs<=db0$cens_obs, db0$time_obs,db0$cens_obs)
    db0$status = ifelse(db0$time_obs<=db0$cens_obs, 1,0)
    db0$time_obs <- NULL
    db0$cens_obs <- NULL
  }else{
    db0$time = ifelse(db0$time_obs<=pmin(db0$cens_obs,tau), db0$time_obs, pmin(db0$cens_obs,tau))
    db0$status = ifelse(db0$time_obs<=pmin(db0$cens_obs,tau), 1,0)
    db0$time_obs <- NULL
    db0$cens_obs <- NULL
  }
  
  
  data = rbind(db0,db1)
  
  return(data)
}

fun_test <- function(data,tau){
  # data = (resp,time,status,treat)
  
  n1 = dim(data[data$treat==1,])[1]
  n0 = dim(data[data$treat==0,])[1]
  
  tau1max= max(subset(data,data$treat==1)$time)
  tau0max= max(subset(data,data$treat==0)$time)
  
  if(tau> min(tau1max,tau0max)){
    tau = min(tau1max,tau0max)
  }
  
  # estimating the diff of restricted mean survival times directly from the whole population 
  k = rmst2(time=data$time,status=data$status,arm=data$treat,tau=tau) 
  num = k$RMST.arm1$rmst[1] - k$RMST.arm0$rmst[1]
  # den = sqrt((k$RMST.arm0$rmst.var)/n0 + (k$RMST.arm1$rmst.var)/n1)
  den = sqrt((k$RMST.arm0$rmst.var) + (k$RMST.arm1$rmst.var))
  test = num/den
  
  return(test)
  
}

fun_simtest <- function(n0,n1,p0,p1,bet0,bet1,lambda_r0,lambda_r1,lambda_nr0,lambda_nr1,lambda_cens,truncated,tau){
  db <- fun_sim(n0,n1,p0,p1,bet0,bet1,lambda_r0,lambda_r1,lambda_nr0,lambda_nr1,lambda_cens,truncated,tau)
  t <- fun_test(db,tau)
  return(t)
}

fun_simtest_LR <- function(n0,n1,p0,p1,bet0,bet1,lambda_r0,lambda_r1,lambda_nr0,lambda_nr1,lambda_cens,truncated,tau){
  db <- fun_sim(n0,n1,p0,p1,bet0,bet1,lambda_r0,lambda_r1,lambda_nr0,lambda_nr1,lambda_cens,truncated,tau)
  t <- survdiff(Surv(time,status)~treat,data=db)
  # fit <- with(db, Surv(time,status)~treat)
  # t <- FHtestrcc(fit)
  # return(t$statistic)
  return(t$chisq)
}

