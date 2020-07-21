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
################################################################

# Weibull parametrization
# S(x) =  exp(- (x/a)^b)

# param_scale: returns the value of the scale parameter a given the survival (s) at time t
param_scale <- function(s,t,shape=1){
  scale = -t/((log(s))^(1/shape))
  return(scale)
}

# The function `survw_samplesize` calculates the sample size according to the distributional parameters of the responders and non-responders, also returns the effect size and the variances in each group.
#
survw_samplesize <- function(ascale0_r,ascale0_nr,delta_p,p0,bshape0,bshape1,ascale1_r,ascale1_nr,ascale_cens,tau,alpha=0.025,beta=0.2){

  z_alpha <- qnorm(1-alpha,0,1)
  z_beta <-  qnorm(1-beta,0,1)
  p1 = delta_p +  p0

  os_effect = survw_effectsize(ascale0_r,ascale0_nr,delta_p,p0,bshape0,bshape1,ascale1_r,ascale1_nr,tau)

  var0 <- var_f(ascale_r=ascale0_r,ascale_nr=ascale0_nr,tau=tau,bshape=bshape0,ascale_cens=ascale_cens,p=p0)
  var1 <- var_f(ascale_r=ascale1_r,ascale_nr=ascale1_nr,tau=tau,bshape=bshape1,ascale_cens=ascale_cens,p=p1)
  ss = ((z_alpha+z_beta)/(os_effect))^2*(var0 + var1)/0.5

  return(list(samplesize=ss,effectsize=os_effect,var0=var0,var1=var1))
}

################################################################
# From Figure 3: Subgroup analysis for event-free survival

s1_r = 0.87
s0_r = 0.55

s1_nr = 0.38
s0_nr = 0.41


ascale1_r = param_scale(s=s1_r, t=5)
ascale0_r = param_scale(s=s0_r, t=5)

ascale1_nr = param_scale(s=s1_nr, t=5)
ascale0_nr = param_scale(s=s0_nr, t=5)

HR_r =  (1/ascale1_r)/(1/ascale0_r)
HR_nr =  (1/ascale1_nr)/(1/ascale0_nr)
HR0 = (1/ascale0_r)/(1/ascale0_nr)
HR1 = (1/ascale1_r)/(1/ascale1_nr)

list(HR_r=HR_r,HR_nr=HR_nr,HR1=HR1,HR0=HR0)


################################################################
# Means

mean0_nr = meanw_f(ascale=ascale0_nr,bshape=1)
mean0_r = meanw_f(ascale=ascale0_r,bshape=1)
mean1_nr = meanw_f(ascale=ascale1_nr,bshape=1)
mean1_r = meanw_f(ascale=ascale1_r,bshape=1)


list(mean0_r=ascale0_r,mean0_nr=ascale0_nr,mean1_r=ascale1_r,mean1_nr=ascale1_nr)

################################################################
# Censoring distribution
# ascale_c = 2*mean0_nr
# mean_cens = meanw_f(ascale=ascale_c,bshape=1)
# (mean_cens)

ascale_c = 7

################################################################
# Response pCR

p1=45/117
p0=23/118
p= 68/235

################################################################
# Effects

Delta_0 = rmstw_f(ascale=ascale0_r,bshape=1,tau=5) - rmstw_f(ascale=ascale0_nr,bshape=1,tau=5)
Delta_r = rmstw_f(ascale=ascale1_r,bshape=1,tau=5) - rmstw_f(ascale=ascale0_r,bshape=1,tau=5)
Delta_nr = rmstw_f(ascale=ascale1_nr,bshape=1,tau=5) - rmstw_f(ascale=ascale0_nr,bshape=1,tau=5)
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

survmixture_f(t=5,ascale_r=ascale0_r, ascale_nr=ascale0_nr, bshape=1, p=p0)
survmixture_f(t=5,ascale_r=ascale1_r, ascale_nr=ascale1_nr, bshape=1, p=p1)

################################################################
# Overall 5-year event-free survival

# original design
survw_effectsize(ascale0_r=ascale0_r,ascale0_nr=ascale0_nr,delta_p=p1-p0,p0=p0,bshape0=1,bshape1=1,ascale1_r=ascale1_r,ascale1_nr=ascale1_nr,tau=5)

n= survw_samplesize(ascale0_r=ascale0_r, ascale0_nr=ascale0_nr, delta_p=p1-p0, p0=p0, bshape0=1, bshape1=1, ascale1_r=ascale1_r,
                    ascale1_nr=ascale1_nr, ascale_cens = ascale_c, tau=5, alpha=0.05, beta=0.2)
(n)

# case III
rmstw_f(ascale=ascale1_nr,bshape=1,tau=5) - rmstw_f(ascale=ascale1_nr,bshape=1,tau=5)
survw_samplesize(ascale0_r=ascale0_r, ascale0_nr=ascale1_nr, delta_p=p1-p0, p0=p0, bshape0=1, bshape1=1, ascale1_r=ascale1_r,
                 ascale1_nr=ascale1_nr, ascale_cens = ascale_c, tau=5, alpha=0.05, beta=0.2)

# case IV
rmstw_f(ascale=ascale1_nr,bshape=1,tau=5) - rmstw_f(ascale=ascale1_nr-0.5,bshape=1,tau=5)
survw_samplesize(ascale0_r=ascale0_r, ascale0_nr=ascale1_nr-0.5, delta_p=p1-p0, p0=p0, bshape0=1, bshape1=1, ascale1_r=ascale1_r,
                 ascale1_nr=ascale1_nr, ascale_cens = ascale_c, tau=5, alpha=0.05, beta=0.2)

################################################################
# Plot sample size

ascale1_nr_vector = c(ascale1_nr,6,7,8,ascale1_nr*2)

fun.1 <- function(x) as.numeric(survw_samplesize(ascale0_r=ascale0_r, ascale0_nr=ascale0_nr, delta_p=p1-p0, p0=p0, bshape0=1, bshape1=1, ascale1_r=x, ascale1_nr=ascale1_nr_vector[1], ascale_cens = ascale_c, tau=5, alpha=0.05, beta=0.2)[1])
fun.2 <- function(x) as.numeric(survw_samplesize(ascale0_r=ascale0_r, ascale0_nr=ascale0_nr, delta_p=p1-p0, p0=p0, bshape0=1, bshape1=1, ascale1_r=x, ascale1_nr=ascale1_nr_vector[2], ascale_cens = ascale_c, tau=5, alpha=0.05, beta=0.2)[1])
fun.3 <- function(x) as.numeric(survw_samplesize(ascale0_r=ascale0_r, ascale0_nr=ascale0_nr, delta_p=p1-p0, p0=p0, bshape0=1, bshape1=1, ascale1_r=x, ascale1_nr=ascale1_nr_vector[3], ascale_cens = ascale_c, tau=5, alpha=0.05, beta=0.2)[1])
fun.4 <- function(x) as.numeric(survw_samplesize(ascale0_r=ascale0_r, ascale0_nr=ascale0_nr, delta_p=p1-p0, p0=p0, bshape0=1, bshape1=1, ascale1_r=x, ascale1_nr=ascale1_nr_vector[4], ascale_cens = ascale_c, tau=5, alpha=0.05, beta=0.2)[1])
fun.5 <- function(x) as.numeric(survw_samplesize(ascale0_r=ascale0_r, ascale0_nr=ascale0_nr, delta_p=p1-p0, p0=p0, bshape0=1, bshape1=1, ascale1_r=x, ascale1_nr=ascale1_nr_vector[5], ascale_cens = ascale_c, tau=5, alpha=0.05, beta=0.2)[1])

ascale1_r_vector = seq(15,40,0.01)

samplesize_vector1 <- sapply(ascale1_r_vector,fun.1)
samplesize_vector2 <- sapply(ascale1_r_vector,fun.2)
samplesize_vector3 <- sapply(ascale1_r_vector,fun.3)
samplesize_vector4 <- sapply(ascale1_r_vector,fun.4)
samplesize_vector5 <- sapply(ascale1_r_vector,fun.5)


dataset = data.frame(ascale1_r_vector=c(ascale1_r_vector,ascale1_r_vector,ascale1_r_vector,ascale1_r_vector,ascale1_r_vector),
                     samplesize_vector=c(samplesize_vector1,samplesize_vector2,samplesize_vector3,samplesize_vector4,samplesize_vector5),
                     ascale1_nr_vector=c(rep(ascale1_nr_vector[1],length(ascale1_r_vector)),
                                         rep(ascale1_nr_vector[2],length(ascale1_r_vector)),
                                         rep(ascale1_nr_vector[3],length(ascale1_r_vector)),
                                         rep(ascale1_nr_vector[4],length(ascale1_r_vector)),
                                         rep(ascale1_nr_vector[5],length(ascale1_r_vector))
                                         )
)

windows(width = 8, height = 8)
ggplot(dataset,aes(x=ascale1_r_vector,y=ascale1_r_vector,col=as.factor(round(ascale1_nr_vector,2)))) + geom_point(aes(ascale1_r_vector,samplesize_vector))  +   geom_vline(xintercept = ascale1_r,linetype="dashed") +guides(col=guide_legend(title='Survival mean for non-responders in trastuzumab group', title.position = "top", title.theme = element_text(size=10))) + xlab('Survival mean for responders in trastuzumab group') + ylab('Sample size')  + theme(legend.position="bottom", axis.text=element_text(size=12), axis.title=element_text(size=12))
# axis.title=element_text(size=10,face="bold"))

################################################################
################################################################
# Plot survival distributions
survmixture_f <- function(t,ascale_r, ascale_nr, bshape=1, p){
  s <- survw_f(t,ascale_r,bshape)*p + survw_f(t,ascale_nr,bshape)*(1-p)
  return(s)
}


funs1 <- function(x) survmixture_f(x,ascale_r=ascale1_r, ascale_nr=ascale1_nr, bshape=1, p=p1)
funs0 <- function(x) survmixture_f(x,ascale_r=ascale0_r, ascale_nr=ascale0_nr, bshape=1, p=p0)


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
funs1_r <- function(x) survw_f(x,ascale=ascale1_r,bshape=1)
funs0_r <- function(x) survw_f(x,ascale=ascale0_r,bshape=1)

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
funs1_nr <- function(x) survw_f(x,ascale=ascale1_nr,bshape=1)
funs0_nr <- function(x) survw_f(x,ascale=ascale0_nr,bshape=1)

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
funs0_r <- function(x) survw_f(x,ascale=ascale0_r,bshape=1)
funs0_nr <- function(x) survw_f(x,ascale=ascale0_nr,bshape=1)

time_vector = seq(0,5,0.01)

survival_vectors0r <- sapply(time_vector,funs0_r)
survival_vectors0nr <- sapply(time_vector,funs0_nr)
# samplesize_vector3 <- sapply(ascale1_r_vector,fun.3)
# samplesize_vector4 <- sapply(ascale1_r_vector,fun.4)
# samplesize_vector5 <- sapply(ascale1_r_vector,fun.5)


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

# denstity_f <- function(x){exp(-x/ascale)/ascale}
# funs1 <- function(x) survmixture_f(x,ascale_r=ascale1_r, ascale_nr=ascale1_nr, bshape=1, p=p1)
# funs0 <- function(x) survmixture_f(x,ascale_r=ascale0_r, ascale_nr=ascale0_nr, bshape=1, p=p0)

hazardratio_f <- function(x){
  return(
    survmixture_f(x,ascale_r=ascale0_r, ascale_nr=ascale0_nr, bshape=1, p=p0)*(p1*(exp(-x/ascale1_r)/ascale1_r)+p0*(exp(-x/ascale1_nr)/ascale1_nr))/  (survmixture_f(x,ascale_r=ascale1_r, ascale_nr=ascale1_nr, bshape=1, p=p1)*(p1*(exp(-x/ascale0_r)/ascale0_r)+p0*(exp(-x/ascale0_nr)/ascale0_nr)))
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

