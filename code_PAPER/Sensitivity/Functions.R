#' ---
#' title: "survmixer: Sample size and effect size calculations for mixture survival  distributions"
#' author: "M. Bofill Roig"
#' date: "`r Sys.Date()`" 
#' output: 
#'   html_document:
#'     code_folding: hide
#' vignette: >
#'   %\VignetteIndexEntry{Computations}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' # html_document
#' ---
#' 
## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
# setwd("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Simulations")

#' 
#' 
#' ## Preamble 
#' 
#' Required packages:
## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------------------------------------------------------------------------
library(knitr) # purl function  
library(pracma)

#' 
#' ## Introduction
#' 
#' Survival mixture distributions are a type of survival distribution in which it is assumed that there are  a proportion of subjects  who had experienced an event that will behave differently from those subjects wo not experienced the event.
#' <!-- a proportion of subjects who will not experience the event.  -->
#' In our model, these ‘responders’ and ‘non-responders’ subjects are modeled separately  using a parametric survival distribution. 
#' 
#' ## Assumptions
#' 
#'   - Let $T_r^{(i)}$ and $T_{nr}^{(i)}$ the time to event for responders and non-responders in the group $i$, respectively.
#'   - Suppose that $T_r^{(i)}$ and $T_{nr}^{(i)}$ follow a weibull distribution with scale parameters $a_{r}^{(i)}$ and $a_{nr}^{(i)}$, respectively, and equal shape parameter  $b^{(i)}$, that is,  $T_r^{(i)} \sim Weibull(a_{r}^{(i)}, b^{(i)})$ and $T_{nr}^{(i)}\sim Weibull(a_{nr}^{(i)}, b^{(i)})$. 
#' 
#' ## Preliminary functions
#' 
#' 
#' The function `survw_f` computes the Weibull survival function. 
## ----survexp-----------------------------------------------------------------------------------------------------------------------------------------------------------
survw_f <- function(t,ascale,bshape){
  return(exp(-(t/ascale)^bshape))
} 

#' 
#' The function `rmstw_f` computes the restricted mean survival times (RMST) according to the Weibull survival function. 
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
rmstw_f <- function(ascale,bshape,tau,low=0){ 
  r <- integrate(survw_f, lower = low, upper = tau, ascale, bshape)$value
  return(r)
} 

#' 
#' 
#' 
#' The function `survw_derivf` computes the derivative of the survival distribution  `survw_f`.
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
survw_derivf <- function(t,ascale,bshape=1){
  return(-bshape*(t/ascale)^bshape*exp(-(t/ascale)^bshape)/t)
} 

#' 
#' The functions `meanw_f` and `medianw_f` calculate the mean and median for Weibull distributions, respectively.
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
medianw_f <- function(ascale,bshape){
  median = ascale*(log(2)^(1/bshape))
  return(median)
}

#' 
#' 
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
meanw_f <- function(ascale,bshape){
  mean = ascale*gamma(1+1/bshape)
  return(mean)
}

#' 
#' 
#' 
#' ## Mixture distribution 
#' 
#' 
#' The function `survmixture_f` computes the survival distribution as a mixture of  of responders and non-responders. The responders and non-responders distributions are assumed to be Weibull distributions.
## ----survmixt----------------------------------------------------------------------------------------------------------------------------------------------------------
survmixture_f <- function(t,ascale_r, ascale_nr, bshape=1, p){
  s <- survw_f(t,ascale_r,bshape)*p + survw_f(t,ascale_nr,bshape)*(1-p)
  return(s)
}

#' 
#' ##  Variance computation
#' 
#' The following three functions are used to calculate the variance of th difference of two RMSTs. `survw_integratef` is used for the integrations; `inside_var` calculates the expression inside the integral; finally, `var_f` computes the variance.
#' 
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
survw_integratef <- function(t,tau, ascale,bshape){
 int <- integrate(survw_f,lower=t, upper=tau,ascale,bshape)$value
 return(int) 
}

#' 
#' 
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
inside_var <- function(t,ascale_r,ascale_nr,tau,bshape,ascale_cens,p){
  
  num = p*sapply(t,survw_integratef,ascale=ascale_r,bshape=bshape,tau=tau)+(1-p)*sapply(t,survw_integratef,ascale=ascale_nr,bshape=bshape,tau=tau)
  den = survmixture_f(t,ascale_r, ascale_nr, bshape, p)
  dervS = p*survw_derivf(t,ascale_r,bshape) + (1-p)*survw_derivf(t,ascale_nr,bshape)

  inside_integral <- (num/den)^2*(1/survw_f(t,ascale_cens,bshape=1))*dervS

  return(-inside_integral)
}

#' 
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
var_f <- function(ascale_r,ascale_nr,tau,bshape,ascale_cens,p){
  integrate(inside_var,lower=0,upper=tau,ascale_r=ascale_r, ascale_nr= ascale_nr,tau=tau,bshape=bshape,ascale_cens=ascale_cens,p=p)$value
}

#' 
#' 
#' ## Effect size and Sample size calculation
#' 
#' The function `survw_effectsize` calculates the effect size according to the distributional parameters of the responders and non-responders.
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

survw_effectsize <- function(ascale0_r,ascale0_nr,delta_p,p0,bshape0,bshape1,ascale1_r,ascale1_nr,tau, 
                             Delta_r=NULL, Delta_0=NULL, Delta_nr=NULL, anticipated_effects=FALSE){ 

  if(anticipated_effects == TRUE){
    os_effect = (p0 + delta_p)*Delta_r + (1-p0-delta_p)*Delta_nr + delta_p*Delta_0  
  }
  if(anticipated_effects == FALSE){
    p1 = delta_p +  p0
    Delta_0 = rmstw_f(ascale=ascale0_r,bshape=bshape0,tau=tau) - rmstw_f(ascale=ascale0_nr,bshape=bshape0,tau=tau)
    
    Delta_r = rmstw_f(ascale=ascale1_r,bshape=bshape1,tau=tau) - rmstw_f(ascale=ascale0_r,bshape=bshape0,tau=tau) 
    Delta_nr = rmstw_f(ascale=ascale1_nr,bshape=bshape1,tau=tau) - rmstw_f(ascale=ascale0_nr,bshape=bshape0,tau=tau)  
    
    k_1 = p1*rmstw_f(ascale=ascale1_r,bshape=bshape1,tau=tau) + (1-p1)*rmstw_f(ascale=ascale1_nr,bshape=bshape1,tau=tau)
    k_0 = p0*rmstw_f(ascale=ascale0_r,bshape=bshape0,tau=tau) + (1-p0)*rmstw_f(ascale=ascale0_nr,bshape=bshape0,tau=tau)
    
    os_effect = (p0 + delta_p)*Delta_r + (1-p0-delta_p)*Delta_nr + delta_p*Delta_0  
  }  
  
  
 return(os_effect) 
}

#' 
#' 
#' 
#' The function `survw_samplesize` calculates the sample size according to the distributional parameters of the responders and non-responders.
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
survw_samplesize <- function(ascale0_r,ascale0_nr,delta_p,p0,bshape0,bshape1,ascale1_r,ascale1_nr,ascale_cens,tau,alpha=0.025,beta=0.2){ 

  z_alpha <- qnorm(1-alpha,0,1)  
  z_beta <-  qnorm(1-beta,0,1)
  p1 = delta_p +  p0
  
  os_effect = survw_effectsize(ascale0_r,ascale0_nr,delta_p,p0,bshape0,bshape1,ascale1_r,ascale1_nr,tau)
  
  var0 <- var_f(ascale_r=ascale0_r,ascale_nr=ascale0_nr,tau=tau,bshape=bshape0,ascale_cens=ascale_cens,p=p0)
  var1 <- var_f(ascale_r=ascale1_r,ascale_nr=ascale1_nr,tau=tau,bshape=bshape1,ascale_cens=ascale_cens,p=p1)
  ss = ((z_alpha+z_beta)/(os_effect))^2*(var0 + var1)/0.5
  
 return(ss) 
}


#' 
#' 
#' 
#' <!-- # Save R code -->
#' <!-- # ```{r} -->
#' <!-- # # require(knitr) # purl function -->
#' <!-- # # purl("Functions_Vignette.Rmd", output = "Functions.R", documentation = 2) -->
#' <!-- # ``` -->
#' 
