#' Effect size calculation for mixture survival distributions
#'
#' @description The function `survw_effectsize` calculates the effect size  according to the information on responders and non-responders.
#'
#' @param lambda0_r scale parameter for the Weibull distribution in the control group for responders
#' @param lambda0_nr scale parameter for the Weibull distribution in the control group for non-responders
#' @param delta_p effect size for the response rate
#' @param p0 event rate for the response
#' @param beta0 shape parameter for the Weibull distribution in the control group
#' @param beta1 shape parameter for the Weibull distribution in the intervention group
#' @param lambda1_r scale parameter for the Weibull distribution in the intervention group for responders
#' @param lambda1_nr scale parameter for the Weibull distribution in the intervention group for non-responders
#' @param tau follow-up
#' @param Delta_r survival effect size between intervention and control groups for responders
#' @param Delta_0  survival effect size between responders and non-responders in the control group
#' @param Delta_nr survival effect size between intervention and control groups for non-responders
#' @param anticipated_effects Logical parameter. If it is TRUE then the effect size is computed based on previous information on the effect sizes on response rate and survival-by-responses (that is, based on Delta_r, Delta_0, Delta_nr); otherwise is based on the distributional parameters (lambda0_r, lambda0_nr, lambda1_r, lambda1_nr, beta0, beta1).
#'
#'
#' @export
#'
#' @return Effect size for overall survival
#' @author Marta Bofill Roig

survw_effectsize <- function(lambda0_r,lambda0_nr,delta_p,p0,beta0,beta1,lambda1_r,lambda1_nr,tau,
                             Delta_r=NULL, Delta_0=NULL, Delta_nr=NULL, anticipated_effects=FALSE){

  if(anticipated_effects == TRUE){
    os_effect = (p0 + delta_p)*Delta_r + (1-p0-delta_p)*Delta_nr + delta_p*Delta_0
  }
  if(anticipated_effects == FALSE){
    p1 = delta_p +  p0
    Delta_0 = rmstw_f(lambda=lambda0_r,bet=beta0,tau=tau) - rmstw_f(lambda=lambda0_nr,bet=beta0,tau=tau)

    Delta_r = rmstw_f(lambda=lambda1_r,bet=beta1,tau=tau) - rmstw_f(lambda=lambda0_r,bet=beta0,tau=tau)
    Delta_nr = rmstw_f(lambda=lambda1_nr,bet=beta1,tau=tau) - rmstw_f(lambda=lambda0_nr,bet=beta0,tau=tau)

    k_1 = p1*rmstw_f(lambda=lambda1_r,bet=beta1,tau=tau) + (1-p1)*rmstw_f(lambda=lambda1_nr,bet=beta1,tau=tau)
    k_0 = p0*rmstw_f(lambda=lambda0_r,bet=beta0,tau=tau) + (1-p0)*rmstw_f(lambda=lambda0_nr,bet=beta0,tau=tau)

    os_effect = (p0 + delta_p)*Delta_r + (1-p0-delta_p)*Delta_nr + delta_p*Delta_0
  }


  return(os_effect)
}

##################################################################################


