################################################################
# EXAMPLE Simulation - MDA Research
# Marta Bofill and Guadalupe Gómez
################################################################

rm(list = ls())
setwd("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Example")
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Simulations/Functions.R')

library(ggplot2)
################################################################
# Functions
param_exp <- function(s,t){
  p = (-log(s))/t
  return(p)
}

################################################################
# From Figure 3: Subgroup analysis for event-free survival

s1_r = 0.87
s0_r = 0.55

s1_nr = 0.38
s0_nr = 0.41

lambda1_r = param_exp(s=s1_r, t=5)
lambda0_r = param_exp(s=s0_r, t=5)

lambda1_nr = param_exp(s=s1_nr, t=5)
lambda0_nr = param_exp(s=s0_nr, t=5)

# Censoring distribution
mean0_nr = meanw_f(lambda=lambda0_nr,bet=1)
lambda_c = 2*mean0_nr


# Means
wscale0_r = 1/lambda0_r
wscale1_r = 1/lambda1_r
wscale0_nr = 1/lambda0_nr
wscale1_nr = 1/lambda1_nr


# Response pCR

p1=45/117
p0=23/118
p= 68/235


################################################################
# Simulation study
library(survRM2)
library(survival)
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Simulations/Sim_Functions.R')
# library(xlsx)

alpha=0.05
beta=0.2
z_alpha <- qnorm(1-alpha,0,1)
z_beta <-  qnorm(1-beta,0,1)
q_chi=qchisq(1-alpha, df=1)

# nsim: number of simulations
nsim=10000

# Using our design
(n= survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, lambda1_r=wscale1_r,
                     lambda1_nr=wscale1_nr, lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2))


power_rmst <- sum(replicate(nsim, fun_simtest(n0=n/2,n1=n/2,
                                              p0=p0,p1=p1,
                                              bet0=1,bet1=1,
                                              lambda_r0=wscale0_r,lambda_r1=wscale1_r,
                                              lambda_nr0=wscale0_nr,lambda_nr1=wscale1_nr,
                                              lambda_cens=1/lambda_c,
                                              # lambda_cens=3,
                                              tau=5,truncated=T)) > z_alpha)/nsim

power_lr <- sum(replicate(nsim,  fun_simtest_LR(n0=n/2,n1=n/2,
                                                p0=p0,p1=p1,
                                                bet0=1,bet1=1,
                                                lambda_r0=wscale0_r,lambda_r1=wscale1_r,
                                                lambda_nr0=wscale0_nr,lambda_nr1=wscale1_nr,
                                                lambda_cens=1/lambda_c,
                                                # lambda_cens=3,
                                                tau=5,truncated=T)) > q_chi)/nsim


(c(power_rmst,power_lr))

###
# Using NOAH's design
powerNT_rmst <- sum(replicate(nsim, fun_simtest(n0=118,n1=117,
                                                p0=p0,p1=p1,
                                                bet0=1,bet1=1,
                                                lambda_r0=wscale0_r,lambda_r1=wscale1_r,
                                                lambda_nr0=wscale0_nr,lambda_nr1=wscale1_nr,
                                                lambda_cens=1/lambda_c,
                                                # lambda_cens=3,
                                                tau=5,truncated=T)) > z_alpha)/nsim

powerNT_lr <- sum(replicate(nsim,  fun_simtest_LR(n0=118,n1=117,
                                                  p0=p0,p1=p1,
                                                  bet0=1,bet1=1,
                                                  lambda_r0=wscale0_r,lambda_r1=wscale1_r,
                                                  lambda_nr0=wscale0_nr,lambda_nr1=wscale1_nr,
                                                  lambda_cens=1/lambda_c,
                                                  # lambda_cens=3,
                                                  tau=5,truncated=T)) > q_chi)/nsim


(c(powerNT_rmst,powerNT_lr))





################################################################
# Simulation EXAMPLE study

# install.packages("survminer")
library(survminer)
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Simulations/Sim_Functions.R')

set.seed(1425)

data_full <- fun_sim(n0=10000,n1=10000,
                     p0=p0,p1=p1,
                     bet0=1,bet1=1,
                     lambda_r0=wscale0_r,lambda_r1=wscale1_r,
                     lambda_nr0=wscale1_nr,lambda_nr1=wscale1_nr,
                     lambda_cens=1/lambda_c,
                     # lambda_cens=3,
                     tau=5,truncated=T)

data_r = subset(data_full,data_full$resp==1)
data_nr = subset(data_full,data_full$resp==0)

################################################################

fit <- survfit(Surv(time, status) ~ treat, data = data_full)
plot_full <- ggsurvplot(
  title    = "Survival curves (full population)",
  fit,                     # survfit object with calculated statistics.
  data = data_full,         # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for
  # point estimates of survival curves.
  xlim = c(0,5),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  # break.time.by = 100,     # break X axis in time intervals by 500.
  ggtheme = theme_grey(), # customize plot and risk table with a theme. theme_grey() theme_light()
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

# Responders
fit_r <- survfit(Surv(time, status) ~ treat, data = data_r)
plot_r <- ggsurvplot(
  title    = "Survival curves for responders",
  fit_r,                     # survfit object with calculated statistics.
  data = data_r,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for
  # point estimates of survival curves.
  xlim = c(0,5),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  # break.time.by = 100,     # break X axis in time intervals by 500.
  ggtheme = theme_grey(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

# nonResponders
fit_nr <- survfit(Surv(time, status) ~ treat, data = data_nr)
plot_nr <- ggsurvplot(
  title    = "Survival curves for non-responders",
  fit_nr,                     # survfit object with calculated statistics.
  data = data_nr,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for
  # point estimates of survival curves.
  xlim = c(0,5),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  # break.time.by = 100,     # break X axis in time intervals by 500.
  ggtheme = theme_grey(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

plots = list(plot_full, plot_r, plot_nr)
arrange_ggsurvplots(plots, ncol=3, nrow=1)

# arrange_ggsurvplots(list(plot_full, plot_r, plot_nr, gg_full), ncol=2, nrow=1)
# cowplot::plot_grid(plot_full, plot_r, plot_nr, gg_full, ncol=2)
################################################################
# HAZARD RATIO PLOTS
# install.packages("rms")
library(rms)

hzh <- hazard.ratio.plot(data_full$treat,Surv(time=data_full$time,event=data_full$status),  e=20,  xlim=c(0,5), legendloc='ll',antilog=T, xlab = "Time")
se <- hzh$se[1,]
lhr <- hzh$log.hazard.ratio[1,]
t <- hzh$time
hz <- exp(lhr)
LI <- exp(lhr - 2*se)
LS <- exp(lhr + 2*se)
data_hr <- data.frame(t,hz,LI,LS)
#

hzh <- hazard.ratio.plot(data_nr$treat,Surv(time=data_nr$time,event=data_nr$status),  e=20,  xlim=c(0,5), legendloc='ll',antilog=T, xlab = "Time")
se <- hzh$se[1,]
lhr <- hzh$log.hazard.ratio[1,]
t <- hzh$time
hz <- exp(lhr)
LI <- exp(lhr - 2*se)
LS <- exp(lhr + 2*se)
datanr_hr <- data.frame(t,hz,LI,LS)
#

hzh <- hazard.ratio.plot(data_r$treat,Surv(time=data_r$time,event=data_r$status),  e=20,  xlim=c(0,5), legendloc='ll',antilog=T, xlab = "Time")
se <- hzh$se[1,]
lhr <- hzh$log.hazard.ratio[1,]
t <- hzh$time
hz <- exp(lhr)
LI <- exp(lhr - 2*se)
LS <- exp(lhr + 2*se)
datar_hr <- data.frame(t,hz,LI,LS)

# PLOTS
# windows()
gg_full <- ggplot(data_hr,aes(x=t,y=hz)) + geom_smooth(colour="#000099") + xlab('Time') + ylab("Hazard Ratio")+ theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))

gg_full


################################################################
# ESTIMATED RMSTs
library(survRM2)
library(survival)

rmst2(time=data_full$time,status=data_full$status,arm=data_full$treat,tau=5)

rmst2(time=data_r$time,status=data_r$status,arm=data_r$treat,tau=5)

rmst2(time=data_nr$time,status=data_nr$status,arm=data_nr$treat,tau=5)








