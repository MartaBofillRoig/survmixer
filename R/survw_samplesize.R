#' Sample size calculation for mixture survival distributions
#'
#' @description The function `survw_samplesize` calculates the sample size according to the distributional parameters of the responders and non-responders.
#'
#' @param ascale0_r scale parameter for the Weibull distribution in the control group for responders
#' @param ascale0_nr scale parameter for the Weibull distribution in the control group for non-responders
#' @param ascale1_r scale parameter for the Weibull distribution in the intervention group for responders
#' @param ascale1_nr scale parameter for the Weibull distribution in the intervention group for non-responders
#' @param delta_p effect size for the response rate
#' @param p0 event rate for the response
#' @param m0_r survival mean for responders in the control group
#' @param m0_nr survival mean for non-responders in the control group
#' @param diffm_r difference in survival means between groups for responders
#' @param diffm_nr difference in survival means between groups for responders
#' @param S0_r tau-year survival rates for responders in the control group
#' @param S0_nr tau-year survival rates for non-responders in the control group
#' @param diffS_r difference in tau-year survival rates for responders
#' @param diffS_nr difference in tau-year survival rates for non-responders
#' @param Delta_r restricted mean survival times (RMST) difference between intervention and control groups for responders
#' @param Delta_nr RMST difference between intervention and control groups for non-responders
#' @param ascale_cens distributional parameter for the exponential distribution for the censoring
#' @param tau follow-up
#' @param bshape0 shape parameter for the Weibull distribution in the control group
#' @param bshape1 shape parameter for the Weibull distribution in the intervention group
#' @param alpha type I error
#' @param beta type II error
#'
#' @export
#'
#' @return Sample size for overall survival
#' @author Marta Bofill Roig

survw_samplesize <- function(ascale0_r,ascale0_nr,ascale1_r,ascale1_nr,delta_p,p0,
                             m0_r, m0_nr, diffm_r, diffm_nr,
                             S0_r, S0_nr, diffS_r, diffS_nr,
                             Delta_r, Delta_nr,
                             ascale_cens,tau,
                             bshape0=1,bshape1=1,alpha=0.025,beta=0.2,
                             ss_strategy=0){

  z_alpha <- qnorm(1-alpha,0,1)
  z_beta <-  qnorm(1-beta,0,1)
  p1 = delta_p +  p0

  if(ss_strategy==1){

    m1_r = diffm_r+m0_r
    m1_nr = diffm_nr+m0_nr

    # note: mean = ascale*gamma(1+1/bshape)
    ascale0_r = m0_r/gamma(1+1/bshape0)
    ascale0_nr = m0_nr/gamma(1+1/bshape0)
    ascale1_r = m1_r/gamma(1+1/bshape1)
    ascale1_nr = m1_nr/gamma(1+1/bshape1)

  }

  if(ss_strategy==2){

    S1_r = diffS_r+ S0_r
    S1_nr = diffS_nr + S0_nr

    ascale0_r = param_scale(s=S0_r,t=tau,shape=bshape0)
    ascale0_nr = param_scale(s=S0_nr,t=tau,shape=bshape0)
    ascale1_r = param_scale(s=S1_r,t=tau,shape=bshape1)
    ascale1_nr = param_scale(s=S1_nr,t=tau,shape=bshape1)

  }

  if(ss_strategy==3){

    ascale0_r = param_scale(s=S0_r,t=tau,shape=bshape0)
    ascale0_nr = param_scale(s=S0_nr,t=tau,shape=bshape0)

    if(beta0==1 && beta1==1){
      ascale1_r = scale1_taylorf(ascale0=ascale0_r,Delta=Delta_r,tau=tau)
      ascale1_nr = scale1_taylorf(ascale0=ascale0_nr,Delta=Delta_nr,tau=tau)
    }
  }

  os_effect = survw_effectsize(ascale0_r,ascale0_nr,delta_p,p0,bshape0,bshape1,ascale1_r,ascale1_nr,tau)

  var0 <- var_f(ascale_r=ascale0_r,ascale_nr=ascale0_nr,tau=tau,bshape=bshape0,ascale_cens=ascale_cens,p=p0)
  var1 <- var_f(ascale_r=ascale1_r,ascale_nr=ascale1_nr,tau=tau,bshape=bshape1,ascale_cens=ascale_cens,p=p1)
  ss = ((z_alpha+z_beta)/(os_effect))^2*(var0 + var1)/0.5

  return(ss)
}

##################################################################################


