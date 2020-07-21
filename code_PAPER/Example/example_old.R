################################################################
# EXAMPLE - MDA Research
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


survw_samplesize <- function(lambda0_r,lambda0_nr,delta_p,p0,beta0,beta1,lambda1_r,lambda1_nr,lambda_cens,tau,alpha=0.025,beta=0.2){

  z_alpha <- qnorm(1-alpha,0,1)
  z_beta <-  qnorm(1-beta,0,1)
  p1 = delta_p +  p0

  os_effect = survw_effectsize(lambda0_r,lambda0_nr,delta_p,p0,beta0,beta1,lambda1_r,lambda1_nr,tau)

  var0 <- var_f(lambda_r=lambda0_r,lambda_nr=lambda0_nr,tau=tau,bet=beta0,lambda_cens=lambda_cens,p=p0)
  var1 <- var_f(lambda_r=lambda1_r,lambda_nr=lambda1_nr,tau=tau,bet=beta1,lambda_cens=lambda_cens,p=p1)
  ss = ((z_alpha+z_beta)/(os_effect))^2*(var0 + var1)/0.5

  return(list(samplesize=ss,effectsize=os_effect,var0=var0,var1=var1))
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

HR_r =  lambda1_r/lambda0_r
HR_nr =  lambda1_nr/lambda0_nr
HR0 = lambda0_r/lambda0_nr
HR1 = lambda1_r/lambda1_nr

list(HR_r=HR_r,HR_nr=HR_nr,HR1=HR1,HR0=HR0)

# Censoring distribution
mean0_nr = meanw_f(lambda=lambda0_nr,bet=1)
lambda_c = 2*mean0_nr


################################################################
# Means

wscale0_r = 1/lambda0_r
wscale1_r = 1/lambda1_r
wscale0_nr = 1/lambda0_nr
wscale1_nr = 1/lambda1_nr

list(mean0_r=wscale0_r,mean0_nr=wscale0_nr,mean1_r=wscale1_r,mean1_nr=wscale1_nr)

################################################################
# Response pCR

p1=45/117
p0=23/118
p= 68/235

################################################################
# Effects

# rmstw_f(lambda=lambda1_nr_vector[5],bet=1,tau=5)-rmstw_f(lambda=lambda1_nr_vector[5],bet=1,tau=5)

Delta_0 = rmstw_f(lambda=wscale0_r,bet=1,tau=5) - rmstw_f(lambda=wscale0_nr,bet=1,tau=5)
Delta_r = rmstw_f(lambda=wscale1_r,bet=1,tau=5) - rmstw_f(lambda=wscale0_r,bet=1,tau=5)
Delta_nr = rmstw_f(lambda=wscale1_nr,bet=1,tau=5) - rmstw_f(lambda=wscale0_nr,bet=1,tau=5)
delta_p = p1-p0


list(Delta_0=round(Delta_0,2),Delta_r=round(Delta_r,2),Delta_nr=round(Delta_nr,2),delta_p=round(delta_p,2))

################################################################
# sample size NOAH TRIAL
z_alpha <- qnorm(1-0.05,0,1)
z_beta <-  qnorm(1-0.2,0,1)

(n_sch = 4*(z_alpha+z_alpha)^2/(0.5*log(0.545)^2))
(n_sch = 4*(z_alpha+z_alpha)^2/(0.5*log(0.77)^2))

################################################################
# Overall 5-year event-free survival

survmixture_f(t=5,lambda_r=wscale0_r, lambda_nr=wscale0_nr, bet=1, p=p0)
survmixture_f(t=5,lambda_r=wscale1_r, lambda_nr=wscale1_nr, bet=1, p=p1)

################################################################
# Overall 5-year event-free survival

survw_effectsize(lambda0_r=wscale0_r,lambda0_nr=wscale0_nr,delta_p=p1-p0,p0=p0,beta0=1,beta1=1,lambda1_r=wscale1_r,lambda1_nr=wscale1_nr,tau=1)
survw_effectsize(lambda0_r=wscale0_r,lambda0_nr=wscale0_nr,delta_p=p1-p0,p0=p0,beta0=1,beta1=1,lambda1_r=wscale1_r,lambda1_nr=wscale1_nr,tau=3)
survw_effectsize(lambda0_r=wscale0_r,lambda0_nr=wscale0_nr,delta_p=p1-p0,p0=p0,beta0=1,beta1=1,lambda1_r=wscale1_r,lambda1_nr=wscale1_nr,tau=5)

survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, lambda1_r=wscale1_r,
                 lambda1_nr=wscale1_nr, lambda_cens = 1/lambda_c, tau=1, alpha=0.05, beta=0.2)
survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, lambda1_r=wscale1_r,
                 lambda1_nr=wscale1_nr, lambda_cens = 1/lambda_c, tau=3, alpha=0.05, beta=0.2)
survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, lambda1_r=wscale1_r,
                 lambda1_nr=wscale1_nr, lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2)

n= survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, lambda1_r=wscale1_r,
                    lambda1_nr=wscale1_nr, lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2)


################################################################
# Plot sample size

lambda1_nr_vector = c(wscale1_nr,6,7,8,wscale1_nr*2)

fun.1 <- function(x) as.numeric(survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, lambda1_r=x,
                                      lambda1_nr=lambda1_nr_vector[1], lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2)[1])
fun.2 <- function(x) as.numeric(survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, lambda1_r=x,
                                                 lambda1_nr=lambda1_nr_vector[2], lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2)[1])
fun.3 <- function(x) as.numeric(survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, lambda1_r=x,
                                                 lambda1_nr=lambda1_nr_vector[3], lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2)[1])
fun.4 <- function(x) as.numeric(survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, lambda1_r=x,
                                                 lambda1_nr=lambda1_nr_vector[4], lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2)[1])
fun.5 <- function(x) as.numeric(survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, lambda1_r=x,
                                                 lambda1_nr=lambda1_nr_vector[5], lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2)[1])

lambda1_r_vector = seq(15,40,0.01)

samplesize_vector1 <- sapply(lambda1_r_vector,fun.1)
samplesize_vector2 <- sapply(lambda1_r_vector,fun.2)
samplesize_vector3 <- sapply(lambda1_r_vector,fun.3)
samplesize_vector4 <- sapply(lambda1_r_vector,fun.4)
samplesize_vector5 <- sapply(lambda1_r_vector,fun.5)


dataset = data.frame(lambda1_r_vector=c(lambda1_r_vector,lambda1_r_vector,lambda1_r_vector,lambda1_r_vector,lambda1_r_vector),
                     samplesize_vector=c(samplesize_vector1,samplesize_vector2,samplesize_vector3,samplesize_vector4,samplesize_vector5),
                     lambda1_nr_vector=c(rep(lambda1_nr_vector[1],length(lambda1_r_vector)),
                                         rep(lambda1_nr_vector[2],length(lambda1_r_vector)),
                                         rep(lambda1_nr_vector[3],length(lambda1_r_vector)),
                                         rep(lambda1_nr_vector[4],length(lambda1_r_vector)),
                                         rep(lambda1_nr_vector[5],length(lambda1_r_vector))
                                         )
)

windows(width = 8, height = 8)
ggplot(dataset,aes(x=lambda1_r_vector,y=lambda1_r_vector,col=as.factor(round(lambda1_nr_vector,2)))) + geom_point(aes(lambda1_r_vector,samplesize_vector))  +   geom_vline(xintercept = wscale1_r,linetype="dashed") +guides(col=guide_legend(title='Rate parameter non-responders in trastuzumab group', title.position = "top", title.theme = element_text(size=10))) + xlab('Rate parameter responders in trastuzumab group') + ylab('Sample size')  + theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12))
# axis.title=element_text(size=10,face="bold"))

################################################################
################################################################
# Plot survival distributions
survmixture_f <- function(t,lambda_r, lambda_nr, bet=1, p){
  s <- survw_f(t,lambda_r,bet)*p + survw_f(t,lambda_nr,bet)*(1-p)
  return(s)
}


funs1 <- function(x) survmixture_f(x,lambda_r=wscale1_r, lambda_nr=wscale1_nr, bet=1, p=p1)
funs0 <- function(x) survmixture_f(x,lambda_r=wscale0_r, lambda_nr=wscale0_nr, bet=1, p=p0)


time_vector = seq(0,5,0.01)

survival_vectors0 <- sapply(time_vector,funs0)
survival_vectors1 <- sapply(time_vector,funs1)


dataset = data.frame(times_vector=c(time_vector,time_vector),
                     survival_vector=c(survival_vectors0,survival_vectors1),
                     resp_vector=c(rep("Chemotherapy group",length(time_vector)),
                                   rep("Trastuzumab group",length(time_vector))
                     )
)

# windows()
ggf <- ggplot(dataset,aes(x=time_vector,y=survival_vector,col=as.factor(resp_vector))) + geom_point(aes(times_vector,survival_vector)) +scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,5)) +guides(col=guide_legend(title='Treatment group', title.position = "top", title.theme = element_text(size=10))) + xlab('Years') + ylab('Survival function')  + theme(legend.position="bottom", axis.text=element_text(size=10), axis.title=element_text(size=10))


################################################################
# In responders
funs1_r <- function(x) survw_f(x,lambda=wscale1_r,bet=1)
funs0_r <- function(x) survw_f(x,lambda=wscale0_r,bet=1)

time_vector = seq(0,5,0.01)

survival_vectors0r <- sapply(time_vector,funs0_r)
survival_vectors1r <- sapply(time_vector,funs1_r)

dataset = data.frame(times_vector=c(time_vector,time_vector),
                     survival_vector=c(survival_vectors0r,survival_vectors1r),
                     resp_vector=c(rep("Chemotherapy group",length(time_vector)),
                                   rep("Trastuzumab group",length(time_vector))
                     )
)

# windows()
ggr <- ggplot(dataset,aes(x=time_vector,y=survival_vector,col=as.factor(resp_vector))) + geom_point(aes(times_vector,survival_vector)) +scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,5)) +guides(col=guide_legend(title='Treatment group', title.position = "top", title.theme = element_text(size=10))) + xlab('Years') + ylab('Survival function for responders')  + theme(legend.position="bottom", axis.text=element_text(size=10), axis.title=element_text(size=10))

################################################################
# In non-responders
funs1_nr <- function(x) survw_f(x,lambda=wscale1_nr,bet=1)
funs0_nr <- function(x) survw_f(x,lambda=wscale0_nr,bet=1)

time_vector = seq(0,5,0.01)

survival_vectors0nr <- sapply(time_vector,funs0_nr)
survival_vectors1nr <- sapply(time_vector,funs1_nr)

dataset = data.frame(times_vector=c(time_vector,time_vector),
                     survival_vector=c(survival_vectors0nr,survival_vectors1nr),
                     resp_vector=c(rep("Chemotherapy group",length(time_vector)),
                                   rep("Trastuzumab group",length(time_vector))
                     )
)

# windows()
ggnr <- ggplot(dataset,aes(x=time_vector,y=survival_vector,col=as.factor(resp_vector))) + geom_point(aes(times_vector,survival_vector)) +scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,5)) +guides(col=guide_legend(title='Treatment group', title.position = "top", title.theme = element_text(size=10))) + xlab('Years') + ylab('Survival function for non-responders')  + theme(legend.position="bottom", axis.text=element_text(size=10), axis.title=element_text(size=10))

################################################################
# In group 0
funs0_r <- function(x) survw_f(x,lambda=wscale0_r,bet=1)
funs0_nr <- function(x) survw_f(x,lambda=wscale0_nr,bet=1)

time_vector = seq(0,5,0.01)

survival_vectors0r <- sapply(time_vector,funs0_r)
survival_vectors0nr <- sapply(time_vector,funs0_nr)
# samplesize_vector3 <- sapply(lambda1_r_vector,fun.3)
# samplesize_vector4 <- sapply(lambda1_r_vector,fun.4)
# samplesize_vector5 <- sapply(lambda1_r_vector,fun.5)


dataset = data.frame(times_vector=c(time_vector,time_vector),
                     survival_vector=c(survival_vectors0r,survival_vectors0nr),
                     resp_vector=c(rep("Responders",length(time_vector)),
                                  rep("Non-responders",length(time_vector))
                     )
)

# windows()
gg0 <- ggplot(dataset,aes(x=time_vector,y=survival_vector,col=as.factor(resp_vector))) + geom_point(aes(times_vector,survival_vector)) +scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,5)) +guides(col=guide_legend(title='Response pCR', title.position = "top", title.theme = element_text(size=10))) + xlab('Years') + ylab('Survival function')  + theme(legend.position="bottom", axis.text=element_text(size=10), axis.title=element_text(size=10))


################################################################
################################################################
# Plot hazard ratio

# denstity_f <- function(x){exp(-x/lambda)/lambda}
# funs1 <- function(x) survmixture_f(x,lambda_r=wscale1_r, lambda_nr=wscale1_nr, bet=1, p=p1)
# funs0 <- function(x) survmixture_f(x,lambda_r=wscale0_r, lambda_nr=wscale0_nr, bet=1, p=p0)

hazardratio_f <- function(x){
  return(
    survmixture_f(x,lambda_r=wscale0_r, lambda_nr=wscale0_nr, bet=1, p=p0)*(p1*(exp(-x/wscale1_r)/wscale1_r)+p0*(exp(-x/wscale1_nr)/wscale1_nr))/  (survmixture_f(x,lambda_r=wscale1_r, lambda_nr=wscale1_nr, bet=1, p=p1)*(p1*(exp(-x/wscale0_r)/wscale0_r)+p0*(exp(-x/wscale0_nr)/wscale0_nr)))
  )
}


time_vector = seq(0,5,0.01)

hazardratio_vector <- sapply(time_vector,hazardratio_f)


dataset = data.frame(times_vector=time_vector,
                     hazardratio_vector=hazardratio_vector
)

# windows()
gghr <- ggplot(dataset,aes(x=time_vector,y=hazardratio_vector)) + geom_point(aes(times_vector,hazardratio_vector)) +scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,5)) +guides(col=guide_legend(title='Treatment group', title.position = "top", title.theme = element_text(size=10))) + xlab('Years') + ylab('Hazard Ratio')  + theme(legend.position="bottom", axis.text=element_text(size=10), axis.title=element_text(size=10))
gghr



################################################################
# windows()
# cowplot::plot_grid(ggf, ggr, ggnr, gg0, nrow=1)

windows(width = 12, height = 12)
cowplot::plot_grid(ggf, gghr, ggr, ggnr, nrow=2)

