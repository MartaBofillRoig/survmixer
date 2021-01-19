
##########################################################################################
# Scenarios simulation MDA Research
# Authors: Marta Bofill Roig, Guadalupe Gómez Melis
##########################################################################################

# Preparing the scenarios
rm(list = ls())

# setwd("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Sensitivity")
setwd("C:/Users/marta.bofill/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Sensitivity")

# install.packages(c("xlsx", "readxl", "knitr", "dplyr", "tidyverse", "tidyr"))
library(xlsx)
library(readxl)
library(knitr)
library(plyr)
library(tidyverse)
library(tidyr)

# getwd()

#'
#' Reading the scenarios
## ------------------------------------------------------------------------
inputs_scenarios <- read_excel("scenarios/inputs_scenarios_exp.xls")
# sum(duplicated(inputs_scenarios))

inputs_scenarios <- inputs_scenarios[!duplicated(inputs_scenarios), ]
p0 = c(0.1,0.3)
# delta_p =c(0.1,0.3)
delta_p =c(0.1)
inputs_scenarios = expand_grid(inputs_scenarios,p0,delta_p)


rm(p0,delta_p)
# write.xlsx(inputs_scenarios, file="scenarios/inputs_scenarios2_test.xls", sheetName="inputs_scenarios", append=TRUE, col.names=TRUE)


##########################################################################################

purl("Functions_Vignette.Rmd", output = "Functions.R", documentation = 2)
source('Functions.R')
# getwd()

#'
#' Significance level and power
## ------------------------------------------------------------------------
alpha_error=0.05
beta_error=0.2
# z_alpha <- qnorm(1-alpha_error,0,1)
# z_beta <-  qnorm(1-beta_error,0,1)


#'
## Creating results table
## ------------------------------------------------------------------------

data <- data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
data <- plyr::rename(data,c('c(0)'='Scenario',
                      X0 = 'PH',
                      X0.1 = 'tau',
                      X0.2 = 'bshape0',
                      X0.3 = 'ascale0_r',
                      X0.4 = 'ascale0_nr',
                      X0.5 = 'bshape1',
                      X0.6 = 'ascale1_r',
                      X0.7 = 'ascale1_nr',
                      X0.8 = 'mean0_r',
                      X0.9 = 'mean1_r',
                      X0.10 = 'mean0_nr',
                      X0.11 = 'mean1_nr',
                      X0.12 = 'median0_r',
                      X0.13 = 'median1_r',
                      X0.14 = 'median0_nr',
                      X0.15 = 'median1_nr',
                      X0.16 = 'diffmean_r',
                      X0.17 = 'diffmean_nr',
                      X0.18 = 'diffmean_0',
                      X0.19 = 'diffmedian_r',
                      X0.20 = 'diffmedian_nr',
                      X0.21 = 'diffmedian_0',
                      X0.22 = 'ascale_cens',
                      X0.23 = 'p1',
                      X0.24 = 'p0',
                      X0.25 = 'Delta_r',
                      X0.26 = 'Delta_nr',
                      X0.27 = 'delta_p',
                      X0.28 = 'Delta_0',
                      X0.29 = 'k_1',
                      X0.30 = 'k_0',
                      X0.31 = 'surv0_r_tau',
                      X0.32 = 'surv0_nr_tau',
                      X0.33 = 'diffsurv0_tau',
                      X0.34 = 'surv1_r_tau',
                      X0.35 = 'surv1_nr_tau',
                      X0.36 = 'diffsurv1_tau',
                      X0.37 = 'surv0_tau',
                      X0.38 = 'surv1_tau',
                      X0.39 = 'diffsurv_tau',
                      X0.40 = 'var0',
                      X0.41 = 'var1',
                      X0.42 = 'os_effect',
                      X0.43 = 'os_samplesize'
                      ))


summary_data <- data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
summary_data <- plyr::rename(summary_data, c('c(0)'='Scenario',
                                       X0 =  'PH',
                                       X0.1 =  'cases',
                                       X0.2 = 'tau',
                                       X0.3 = 'bshape0',
                                       X0.4 = 'bshape1',
                                       X0.5 = 'mean0_r',
                                       X0.6 = 'mean0_nr',
                                       X0.7 = 'median0_r',
                                       X0.8 = 'median0_nr',
                                       X0.9 = 'diffmean_r',
                                       X0.10 = 'diffmean_nr',
                                       X0.11 = 'diffmean_0',
                                       X0.12 = 'diffmedian_r',
                                       X0.13 = 'diffmedian_nr',
                                       X0.14 = 'diffmedian_0',
                                       X0.15 = 'p0',
                                       X0.16 = 'delta_p',
                                       X0.17 = 'Delta_r',
                                       X0.18 = 'Delta_nr',
                                       X0.19 = 'Delta_0',
                                       X0.20 = 'k_0',
                                       X0.21 = 'surv0_tau',
                                       X0.22 = 'diffsurv_tau',
                                       X0.23 = 'os_effect',
                                       X0.24 = 'os_samplesize'))



it=1

for(i in 1:dim(inputs_scenarios)[1]){

  if(inputs_scenarios$bshape1[i] == inputs_scenarios$bshape0[i]){
    PH = TRUE
  }else{
    PH = FALSE
  }

  cases = inputs_scenarios$cases[i]

  p1 = inputs_scenarios$delta_p[i] +  inputs_scenarios$p0[i]
  Delta_0 = rmstw_f(ascale=inputs_scenarios$ascale0_r[i],bshape=inputs_scenarios$bshape0[i],tau=inputs_scenarios$tau[i]) - rmstw_f(ascale=inputs_scenarios$ascale0_nr[i],bshape=inputs_scenarios$bshape0[i],tau=inputs_scenarios$tau[i])

  Delta_r = rmstw_f(ascale=inputs_scenarios$ascale1_r[i],bshape=inputs_scenarios$bshape1[i],tau=inputs_scenarios$tau[i]) - rmstw_f(ascale=inputs_scenarios$ascale0_r[i],bshape=inputs_scenarios$bshape0[i],tau=inputs_scenarios$tau[i])
  Delta_nr = rmstw_f(ascale=inputs_scenarios$ascale1_nr[i],bshape=inputs_scenarios$bshape1[i],tau=inputs_scenarios$tau[i]) - rmstw_f(ascale=inputs_scenarios$ascale0_nr[i],bshape=inputs_scenarios$bshape0[i],tau=inputs_scenarios$tau[i])

  k_1 = p1*rmstw_f(ascale=inputs_scenarios$ascale1_r[i],bshape=inputs_scenarios$bshape1[i],tau=inputs_scenarios$tau[i]) + (1-p1)*rmstw_f(ascale=inputs_scenarios$ascale1_nr[i],bshape=inputs_scenarios$bshape1[i],tau=inputs_scenarios$tau[i])
  k_0 = inputs_scenarios$p0[i]*rmstw_f(ascale=inputs_scenarios$ascale0_r[i],bshape=inputs_scenarios$bshape0[i],tau=inputs_scenarios$tau[i]) + (1-inputs_scenarios$p0[i])*rmstw_f(ascale=inputs_scenarios$ascale0_nr[i],bshape=inputs_scenarios$bshape0[i],tau=inputs_scenarios$tau[i])

  os_effect = survw_effectsize(ascale0_r=inputs_scenarios$ascale0_r[i],
                                ascale0_nr=inputs_scenarios$ascale0_nr[i],
                                delta_p=inputs_scenarios$delta_p[i],
                                p0=inputs_scenarios$p0[i],
                                bshape0=inputs_scenarios$bshape0[i],
                                bshape1=inputs_scenarios$bshape1[i],
                                ascale1_r=inputs_scenarios$ascale1_r[i],
                                ascale1_nr=inputs_scenarios$ascale1_nr[i],
                                tau=inputs_scenarios$tau[i],
                                anticipated_effects=FALSE)

  # os_effect2 = survw_effectsize(delta_p=inputs_scenarios$delta_p[i],
  #                               p0=inputs_scenarios$p0[i],
  #                               Delta_r=Delta_r, Delta_0=Delta_0, Delta_nr=Delta_nr, anticipated_effects=TRUE)

  mean0_r = meanw_f(ascale=inputs_scenarios$ascale0_r[i],bshape=inputs_scenarios$bshape0[i])
  mean1_r = meanw_f(ascale=inputs_scenarios$ascale1_r[i],bshape=inputs_scenarios$bshape1[i])
  mean0_nr = meanw_f(ascale=inputs_scenarios$ascale0_nr[i],bshape=inputs_scenarios$bshape0[i])
  mean1_nr = meanw_f(ascale=inputs_scenarios$ascale1_nr[i],bshape=inputs_scenarios$bshape1[i])

  ascale_cens_value = 2*mean0_nr

  diffmean_r = mean1_r - mean0_r
  diffmean_nr = mean1_nr - mean0_nr
  diffmean_0 = mean0_r - mean0_nr

  median0_r = medianw_f(ascale=inputs_scenarios$ascale0_r[i],bshape=inputs_scenarios$bshape0[i])
  median1_r = medianw_f(ascale=inputs_scenarios$ascale1_r[i],bshape=inputs_scenarios$bshape1[i])
  median0_nr = medianw_f(ascale=inputs_scenarios$ascale0_nr[i],bshape=inputs_scenarios$bshape0[i])
  median1_nr = medianw_f(ascale=inputs_scenarios$ascale1_nr[i],bshape=inputs_scenarios$bshape1[i])

  diffmedian_r = median1_r - median0_r
  diffmedian_nr = median1_nr - median0_nr
  diffmedian_0 = median0_r - median0_nr

  # var0 = var_f(ascale_r=inputs_scenarios$ascale0_r[i],ascale_nr=inputs_scenarios$ascale0_nr[i],tau=inputs_scenarios$tau[i],bshape=inputs_scenarios$bshape0[i],ascale_cens=inputs_scenarios$ascale_cens[i],p=inputs_scenarios$p0[i])
  # var1 = var_f(ascale_r=inputs_scenarios$ascale1_r[i],ascale_nr=inputs_scenarios$ascale1_nr[i],tau=inputs_scenarios$tau[i],bshape=inputs_scenarios$bshape1[i],ascale_cens=inputs_scenarios$ascale_cens[i],p=p1)

  var0 = var_f(ascale_r=inputs_scenarios$ascale0_r[i],ascale_nr=inputs_scenarios$ascale0_nr[i],tau=inputs_scenarios$tau[i],bshape=inputs_scenarios$bshape0[i],ascale_cens=ascale_cens_value,p=inputs_scenarios$p0[i])
  var1 = var_f(ascale_r=inputs_scenarios$ascale1_r[i],ascale_nr=inputs_scenarios$ascale1_nr[i],tau=inputs_scenarios$tau[i],bshape=inputs_scenarios$bshape1[i],ascale_cens=ascale_cens_value,p=p1)

  # os_samplesize = ((z_alpha+z_beta)/(os_effect))^2*(var0 + var1)/0.5
  os_samplesize = survw_samplesize(ascale0_r=inputs_scenarios$ascale0_r[i],
                                    ascale0_nr=inputs_scenarios$ascale0_nr[i],
                                    delta_p=inputs_scenarios$delta_p[i],
                                    p0=inputs_scenarios$p0[i],
                                    bshape0=inputs_scenarios$bshape0[i],
                                    bshape1=inputs_scenarios$bshape1[i],
                                    ascale1_r=inputs_scenarios$ascale1_r[i],
                                    ascale1_nr=inputs_scenarios$ascale1_nr[i],
                                    ascale_cens=ascale_cens_value,
                                    tau=inputs_scenarios$tau[i],alpha=alpha_error,beta=beta_error)

  surv0_r_tau = survw_f(t=inputs_scenarios$tau[i],ascale=inputs_scenarios$ascale0_r[i],bshape=inputs_scenarios$bshape0[i])
  surv0_nr_tau = survw_f(t=inputs_scenarios$tau[i],ascale=inputs_scenarios$ascale0_nr[i],bshape=inputs_scenarios$bshape0[i])

  diffsurv0_tau =  surv0_r_tau - surv0_nr_tau

  surv1_r_tau = survw_f(t=inputs_scenarios$tau[i],ascale=inputs_scenarios$ascale1_r[i],bshape=inputs_scenarios$bshape1[i])
  surv1_nr_tau = survw_f(t=inputs_scenarios$tau[i],ascale=inputs_scenarios$ascale1_nr[i],bshape=inputs_scenarios$bshape1[i])

  diffsurv1_tau =  surv1_r_tau - surv1_nr_tau

  surv0_tau = survmixture_f(t=inputs_scenarios$tau[i], ascale_r=inputs_scenarios$ascale0_r[i],ascale_nr=inputs_scenarios$ascale0_nr[i], bshape=inputs_scenarios$bshape0[i], p=inputs_scenarios$p0[i])
  surv1_tau = survmixture_f(t=inputs_scenarios$tau[i], ascale_r=inputs_scenarios$ascale1_r[i],ascale_nr=inputs_scenarios$ascale1_nr[i], bshape=inputs_scenarios$bshape1[i], p=p1)

  diffsurv_tau =  surv1_tau - surv0_tau

  data[it,]<- c(PH,
                inputs_scenarios$tau[i],
                inputs_scenarios$bshape0[i],
                inputs_scenarios$ascale0_r[i],
                inputs_scenarios$ascale0_nr[i],
                inputs_scenarios$bshape1[i],
                inputs_scenarios$ascale1_r[i],
                inputs_scenarios$ascale1_nr[i],
                mean0_r,
                mean1_r,
                mean0_nr,
                mean1_nr,
                median0_r,
                median1_r,
                median0_nr,
                median1_nr,
                diffmean_r,
                diffmean_nr,
                diffmean_0,
                diffmedian_r,
                diffmedian_nr,
                diffmedian_0,
                ascale_cens_value,
                # inputs_scenarios$ascale_cens[i],
                p1,
                inputs_scenarios$p0[i],
                Delta_r,
                Delta_nr,
                inputs_scenarios$delta_p[i],
                Delta_0,
                k_1,
                k_0,
                surv0_r_tau,
                surv0_nr_tau,
                diffsurv0_tau,
                surv1_r_tau,
                surv1_nr_tau,
                diffsurv1_tau,
                surv0_tau,
                surv1_tau,
                diffsurv_tau,
                var0,
                var1,
                os_effect,
                os_samplesize)
                # os_effect2,os_effect3)


  summary_data[it,]<- c(PH, cases,
                inputs_scenarios$tau[i],
                inputs_scenarios$bshape0[i],
                inputs_scenarios$bshape1[i],
                mean0_r,
                mean0_nr,
                median0_r,
                median0_nr,
                diffmean_r,
                diffmean_nr,
                diffmean_0,
                diffmedian_r,
                diffmedian_nr,
                diffmedian_0,
                inputs_scenarios$p0[i],
                inputs_scenarios$delta_p[i],
                Delta_r,
                Delta_nr,
                Delta_0,
                k_0,
                surv0_tau,
                diffsurv_tau,
                os_effect,
                os_samplesize)
  it=it+1
}



#'
## ------------------------------------------------------------------------

# head(data,20)

rm(i,it,
   PH,
   cases,
   p1,
   Delta_r,
   Delta_nr,
   Delta_0,
   k_1,
   k_0,
   mean0_r,
   mean1_r,
   mean0_nr,
   mean1_nr,
   median0_r,
   median1_r,
   median0_nr,
   median1_nr,
   diffmean_r,
   diffmean_nr,
   diffmean_0,
   diffmedian_r,
   diffmedian_nr,
   diffmedian_0,
   surv0_tau,
   surv1_tau,
   diffsurv_tau,
   var0,
   var1,
   os_effect,
   os_samplesize,
   beta_error,alpha_error)


#' # save scenarios
write.xlsx(inputs_scenarios, file="scenarios/complete_scenarios.xls", sheetName="inputs_scenarios", col.names=TRUE)
write.xlsx(data, file="scenarios/complete_scenarios.xls", sheetName="complete_results", append=TRUE, col.names=TRUE)
write.xlsx(summary_data, file="scenarios/complete_scenarios.xls", sheetName="summary_results", append=TRUE, col.names=TRUE)
