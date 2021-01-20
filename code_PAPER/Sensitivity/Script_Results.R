
################################################################
# SIMULATIONS RESULTS - MDA Research
# Marta Bofill and Guadalupe Gómez
################################################################

rm(list = ls())
load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Sensitivity/results_sim/RESULTS_sim.RData")

library(ggplot2)
library(gridExtra)
library(ggpubr)


################################################################
################################################################

data$cases = 4
for(i in 1:dim(data)[1]){
  if(data$Delta_r[i]==0 & data$Delta_nr[i]==0) data$cases[i]= 1
  if(data$Delta_r[i]==0 & data$Delta_nr[i]!=0) data$cases[i]= 2
  if(data$Delta_r[i]!=0 & data$Delta_nr[i]==0) data$cases[i]=3
}
data$cases = as.factor(data$cases)
summary(data$cases)

data = subset(data, data$os_samplesize>100 & data$os_samplesize<5000)

summary(data)
summary(data$Test_Reject)
summary(data$Test_Reject_LR)
summary(data$Test_Reject_size)
summary(data$Test_Reject_LR_size)

################################################################
################################################################


load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Sensitivity/results_sim/RESULTS_sim_dec.RData")
data = subset(data, data$os_samplesize>100 & data$os_samplesize<5000)

summary(data)
summary(data$Test_Reject)
summary(data$Test_Reject_LR)
summary(data$Test_Reject_size)
summary(data$Test_Reject_LR_size)

################################################################
################################################################

load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Sensitivity/results_sim/RESULTS_sim_inc.RData")
data = subset(data, data$os_samplesize>100 & data$os_samplesize<5000)

summary(data)
summary(data$Test_Reject)
summary(data$Test_Reject_LR)
summary(data$Test_Reject_size)
summary(data$Test_Reject_LR_size)
