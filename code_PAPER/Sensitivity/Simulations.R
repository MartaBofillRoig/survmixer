
################################################################
# SIMULATIONS MDA Research
# Marta Bofill and Guadalupe G�mez
################################################################

rm(list = ls())
setwd("C:/Users/marta.bofill/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Sensitivity")
source('C:/Users/marta.bofill/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Sensitivity/Sim_Functions.R')
path.results <- 'C:/Users/marta.bofill/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Sensitivity/results_sim/'

#####################################################################################
# PREAMBLE
#####################################################################################
library(survRM2)
library(survival)
library(xlsx)

data <- read.xlsx(file="scenarios/complete_scenarios.xls", sheetIndex=2, header=TRUE, colClasses=NA)

#####################################################################################
# Parameters
alpha=0.05
beta=0.2
z_alpha <- qnorm(1-alpha,0,1)
z_beta <-  qnorm(1-beta,0,1)
q_chi=qchisq(1-alpha, df=1)

# nsim: number of simulations
nsim=1000


#####################################################################################
# simulation seed
set.seed(9876)

t0=Sys.time()
data$Test_Reject=0
data$Test_Reject_LR=0

for(i in 1:dim(data)[1]){
# for(i in 1:1){
  data$Test_Reject[i] <- sum(replicate(nsim,
                                       fun_simtest(n0=data$os_samplesize[i]/2,n1=data$os_samplesize[i]/2,
                                                   p0=data$p0[i],p1=data$p1[i],
                                                   # bshape0=data$bshape0[i],bshape1=data$bshape1[i],
                                                   bshape0=1.5,bshape1=1.5,
                                                   ascale0_r=data$ascale0_r[i],ascale1_r=data$ascale1_r[i],
                                                   ascale0_nr=data$ascale0_nr[i],ascale1_nr=data$ascale1_nr[i],
                                                   ascale_cens=data$ascale_cens[i],
                                                   truncated=T,
                                                   tau=data$tau[i])) > z_alpha)/nsim

  data$Test_Reject_LR[i] <- sum(replicate(nsim,
                                          fun_simtest_LR(n0=data$os_samplesize[i]/2,n1=data$os_samplesize[i]/2,
                                                   p0=data$p0[i],p1=data$p1[i],
                                                   # bshape0=data$bshape0[i],bshape1=data$bshape1[i],
                                                   bshape0=1.5,bshape1=1.5,
                                                   ascale0_r=data$ascale0_r[i],ascale1_r=data$ascale1_r[i],
                                                   ascale0_nr=data$ascale0_nr[i],ascale1_nr=data$ascale1_nr[i],
                                                   ascale_cens=data$ascale_cens[i],
                                                   truncated=T,
                                                   tau=data$tau[i])) > q_chi)/nsim

  t1=Sys.time()-t0
  cat(i, "\t", data$Test_Reject[i], "\t", t1, "\n", file="results_sim/LOG_power.txt", append=TRUE)
}

t1=Sys.time()-t0
cat(t1, "\n", file="results_sim/LOG_power.txt", append=TRUE)
(t1)

rm(i)
save.image("C:/Users/marta.bofill/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Sensitivity/results_sim/RESULTS_sim.RData")


#####################################################################################
# simulation seed
set.seed(543767)

t0=Sys.time()
data$Test_Reject_size=0
data$Test_Reject_LR_size=0

for(i in 1:dim(data)[1]){
  data$Test_Reject_size[i] <- sum(replicate(nsim,
                                            fun_simtest(n0=data$os_samplesize[i]/2,n1=data$os_samplesize[i]/2,
                                                        p0=data$p0[i],p1=data$p0[i],
                                                        # bshape0=data$bshape0[i],bshape1=data$bshape0[i],
                                                        bshape0=2,bshape1=2,
                                                        ascale0_r=data$ascale0_r[i],ascale1_r=data$ascale0_r[i],
                                                        ascale0_nr=data$ascale0_nr[i],ascale1_nr=data$ascale0_nr[i],
                                                        ascale_cens=data$ascale_cens[i],
                                                        truncated=T,
                                                        tau=data$tau[i])) > z_alpha,na.rm = T)/nsim

  data$Test_Reject_LR_size[i] <- sum(replicate(nsim,
                                          fun_simtest_LR(n0=data$os_samplesize[i]/2,n1=data$os_samplesize[i]/2,
                                                         p0=data$p0[i],p1=data$p0[i],
                                                         # bshape0=data$bshape0[i],bshape1=data$bshape0[i],
                                                         bshape0=2,bshape1=2,
                                                         ascale0_r=data$ascale0_r[i],ascale1_r=data$ascale0_r[i],
                                                         ascale0_nr=data$ascale0_nr[i],ascale1_nr=data$ascale0_nr[i],
                                                         ascale_cens=data$ascale_cens[i],
                                                         truncated=T,
                                                         tau=data$tau[i])) > q_chi)/nsim

  t1=Sys.time()-t0
  cat(i, "\t", data$Test_Reject_size[i], "\t", t1, "\n", file="results_sim/LOG_size.txt", append=TRUE)

}

t1=Sys.time()-t0
cat(t1, "\n", file="results_sim/LOG_size.txt", append=TRUE)
(t1)

rm(i)
save.image("C:/Users/marta.bofill/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Sensitivity/results_sim/RESULTS_sim_inc.RData")

# write.xlsx(): append=FALSE when we are overwriting the sheet. Otherwise use append=TRUE
# write.xlsx(data, file="scenarios/complete_scenarios_results.xls", sheetName="complete_results_sim", append=FALSE, col.names=TRUE)
# write.xlsx(data, file="scenarios/complete_scenarios.xls", sheetName="complete_results_sim", append=TRUE, col.names=TRUE)

#####################################################################################
# simulation seed
set.seed(2231)

t0=Sys.time()
data$Test_Reject_size=0
data$Test_Reject_LR_size=0

for(i in 1:dim(data)[1]){
  data$Test_Reject_size[i] <- sum(replicate(nsim,
                                            fun_simtest(n0=data$os_samplesize[i]/2,n1=data$os_samplesize[i]/2,
                                                        p0=data$p0[i],p1=data$p0[i],
                                                        # bshape0=data$bshape0[i],bshape1=data$bshape0[i],
                                                        bshape0=2,bshape1=2,
                                                        ascale0_r=data$ascale0_r[i],ascale1_r=data$ascale0_r[i],
                                                        ascale0_nr=data$ascale0_nr[i],ascale1_nr=data$ascale0_nr[i],
                                                        ascale_cens=data$ascale_cens[i],
                                                        truncated=T,
                                                        tau=data$tau[i])) > z_alpha,na.rm = T)/nsim

  data$Test_Reject_LR_size[i] <- sum(replicate(nsim,
                                               fun_simtest_LR(n0=data$os_samplesize[i]/2,n1=data$os_samplesize[i]/2,
                                                              p0=data$p0[i],p1=data$p0[i],
                                                              # bshape0=data$bshape0[i],bshape1=data$bshape0[i],
                                                              bshape0=2,bshape1=2,
                                                              ascale0_r=data$ascale0_r[i],ascale1_r=data$ascale0_r[i],
                                                              ascale0_nr=data$ascale0_nr[i],ascale1_nr=data$ascale0_nr[i],
                                                              ascale_cens=data$ascale_cens[i],
                                                              truncated=T,
                                                              tau=data$tau[i])) > q_chi)/nsim

  t1=Sys.time()-t0
  cat(i, "\t", data$Test_Reject_size[i], "\t", t1, "\n", file="results_sim/LOG_size.txt", append=TRUE)

}

t1=Sys.time()-t0
cat(t1, "\n", file="results_sim/LOG_size.txt", append=TRUE)
(t1)

rm(i)
save.image("C:/Users/marta.bofill/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Sensitivity/results_sim/RESULTS_sim_inc.RData")

# write.xlsx(): append=FALSE when we are overwriting the sheet. Otherwise use append=TRUE
# write.xlsx(data, file="scenarios/complete_scenarios_results.xls", sheetName="complete_results_sim", append=FALSE, col.names=TRUE)
# write.xlsx(data, file="scenarios/complete_scenarios.xls", sheetName="complete_results_sim", append=TRUE, col.names=TRUE)
#####################################################################################
#####################################################################################

# DECREASING HAZARDS
#####################################################################################
# simulation seed
set.seed(9876)

t0=Sys.time()
data$Test_Reject=0
data$Test_Reject_LR=0

for(i in 1:dim(data)[1]){
  # for(i in 1:1){
  data$Test_Reject[i] <- sum(replicate(nsim,
                                       fun_simtest(n0=data$os_samplesize[i]/2,n1=data$os_samplesize[i]/2,
                                                   p0=data$p0[i],p1=data$p1[i],
                                                   # bshape0=data$bshape0[i],bshape1=data$bshape1[i],
                                                   bshape0=0.5,bshape1=0.5,
                                                   ascale0_r=data$ascale0_r[i],ascale1_r=data$ascale1_r[i],
                                                   ascale0_nr=data$ascale0_nr[i],ascale1_nr=data$ascale1_nr[i],
                                                   ascale_cens=data$ascale_cens[i],
                                                   truncated=T,
                                                   tau=data$tau[i])) > z_alpha)/nsim

  data$Test_Reject_LR[i] <- sum(replicate(nsim,
                                          fun_simtest_LR(n0=data$os_samplesize[i]/2,n1=data$os_samplesize[i]/2,
                                                         p0=data$p0[i],p1=data$p1[i],
                                                         # bshape0=data$bshape0[i],bshape1=data$bshape1[i],
                                                         bshape0=0.5,bshape1=0.5,
                                                         ascale0_r=data$ascale0_r[i],ascale1_r=data$ascale1_r[i],
                                                         ascale0_nr=data$ascale0_nr[i],ascale1_nr=data$ascale1_nr[i],
                                                         ascale_cens=data$ascale_cens[i],
                                                         truncated=T,
                                                         tau=data$tau[i])) > q_chi)/nsim

  t1=Sys.time()-t0
  cat(i, "\t", data$Test_Reject[i], "\t", t1, "\n", file="results_sim/LOG_power.txt", append=TRUE)
}

t1=Sys.time()-t0
cat(t1, "\n", file="results_sim/LOG_power.txt", append=TRUE)
(t1)

rm(i)
save.image("C:/Users/marta.bofill/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Sensitivity/results_sim/RESULTS_sim_dec.RData")


#####################################################################################
# simulation seed
set.seed(1903)

t0=Sys.time()
data$Test_Reject_size=0
data$Test_Reject_LR_size=0

for(i in 1:dim(data)[1]){
  data$Test_Reject_size[i] <- sum(replicate(nsim,
                                            fun_simtest(n0=data$os_samplesize[i]/2,n1=data$os_samplesize[i]/2,
                                                        p0=data$p0[i],p1=data$p0[i],
                                                        # bshape0=data$bshape0[i],bshape1=data$bshape0[i],
                                                        bshape0=0.5,bshape1=0.5,
                                                        ascale0_r=data$ascale0_r[i],ascale1_r=data$ascale0_r[i],
                                                        ascale0_nr=data$ascale0_nr[i],ascale1_nr=data$ascale0_nr[i],
                                                        ascale_cens=data$ascale_cens[i],
                                                        truncated=T,
                                                        tau=data$tau[i])) > z_alpha,na.rm = T)/nsim

  data$Test_Reject_LR_size[i] <- sum(replicate(nsim,
                                               fun_simtest_LR(n0=data$os_samplesize[i]/2,n1=data$os_samplesize[i]/2,
                                                              p0=data$p0[i],p1=data$p0[i],
                                                              # bshape0=data$bshape0[i],bshape1=data$bshape0[i],
                                                              bshape0=0.5,bshape1=0.5,
                                                              ascale0_r=data$ascale0_r[i],ascale1_r=data$ascale0_r[i],
                                                              ascale0_nr=data$ascale0_nr[i],ascale1_nr=data$ascale0_nr[i],
                                                              ascale_cens=data$ascale_cens[i],
                                                              truncated=T,
                                                              tau=data$tau[i])) > q_chi)/nsim

  t1=Sys.time()-t0
  cat(i, "\t", data$Test_Reject_size[i], "\t", t1, "\n", file="results_sim/LOG_size.txt", append=TRUE)

}

t1=Sys.time()-t0
cat(t1, "\n", file="results_sim/LOG_size.txt", append=TRUE)
(t1)

rm(i)
save.image("C:/Users/marta.bofill/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Sensitivity/results_sim/RESULTS_sim_dec.RData")

# write.xlsx(): append=FALSE when we are overwriting the sheet. Otherwise use append=TRUE
# write.xlsx(data, file="scenarios/complete_scenarios_results.xls", sheetName="complete_results_sim", append=FALSE, col.names=TRUE)
# write.xlsx(data, file="scenarios/complete_scenarios.xls", sheetName="complete_results_sim", append=TRUE, col.names=TRUE)

#####################################################################################
#####################################################################################
