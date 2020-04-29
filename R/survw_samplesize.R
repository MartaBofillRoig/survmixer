#' Sample size calculation for mixture survival distributions
#'
#' @description The function `survw_samplesize` calculates the sample size according to the distributional parameters of the responders and non-responders.
#'
#' @param lambda0_r scale parameter for the Weibull distribution in the control group for responders
#' @param lambda0_nr scale parameter for the Weibull distribution in the control group for non-responders
#' @param delta_p effect size for the response rate
#' @param p0 event rate for the response
#' @param beta0 shape parameter for the Weibull distribution in the control group
#' @param beta1 shape parameter for the Weibull distribution in the intervention group
#' @param lambda1_r scale parameter for the Weibull distribution in the intervention group for responders
#' @param lambda1_nr scale parameter for the Weibull distribution in the intervention group for non-responders
#' @param lambda_cens distributional parameter for the exponential distribution for the censoring
#' @param tau follow-up
#' @param alpha type I error
#' @param beta type II error
#'
#' @export
#'
#' @return Sample size for overall survival
#' @author Marta Bofill Roig

survw_samplesize <- function(lambda0_r,lambda0_nr,delta_p,p0,beta0,beta1,lambda1_r,lambda1_nr,lambda_cens,tau,alpha=0.025,beta=0.2){

  z_alpha <- qnorm(1-alpha,0,1)
  z_beta <-  qnorm(1-beta,0,1)
  p1 = delta_p +  p0

  os_effect = survw_effectsize(lambda0_r,lambda0_nr,delta_p,p0,beta0,beta1,lambda1_r,lambda1_nr,tau)

  var0 <- var_f(lambda_r=lambda0_r,lambda_nr=lambda0_nr,tau=tau,bet=beta0,lambda_cens=lambda_cens,p=p0)
  var1 <- var_f(lambda_r=lambda1_r,lambda_nr=lambda1_nr,tau=tau,bet=beta1,lambda_cens=lambda_cens,p=p1)
  ss = ((z_alpha+z_beta)/(os_effect))^2*(var0 + var1)/0.5

  return(ss)
}

##################################################################################


