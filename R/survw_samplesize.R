#' Sample size calculation for mixture survival distributions
#'
#' @description The function `survw_samplesize` calculates the sample size according to the distributional parameters of the responders and non-responders.
#'
#' @param ascale0_r scale parameter for the Weibull distribution in the control group for responders
#' @param ascale0_nr scale parameter for the Weibull distribution in the control group for non-responders
#' @param delta_p effect size for the response rate
#' @param p0 event rate for the response
#' @param bshape0 shape parameter for the Weibull distribution in the control group
#' @param bshape1 shape parameter for the Weibull distribution in the intervention group
#' @param ascale1_r scale parameter for the Weibull distribution in the intervention group for responders
#' @param ascale1_nr scale parameter for the Weibull distribution in the intervention group for non-responders
#' @param ascale_cens distributional parameter for the exponential distribution for the censoring
#' @param tau follow-up
#' @param alpha type I error
#' @param beta type II error
#'
#' @export
#'
#' @return Sample size for overall survival
#' @author Marta Bofill Roig

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

##################################################################################


