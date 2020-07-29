#' Effect size calculation for mixture survival distributions
#'
#' @description The function `survm_effectsize` calculates the effect size in terms of the difference of restricted mean survival times (RMST) according to the information on responders and non-responders.
#'
#' @param ascale0_r scale parameter for the Weibull distribution in the control group for responders
#' @param ascale0_nr scale parameter for the Weibull distribution in the control group for non-responders
#' @param delta_p effect size for the response rate
#' @param p0 event rate for the response
#' @param bshape0 shape parameter for the Weibull distribution in the control group
#' @param bshape1 shape parameter for the Weibull distribution in the intervention group
#' @param ascale1_r scale parameter for the Weibull distribution in the intervention group for responders
#' @param ascale1_nr scale parameter for the Weibull distribution in the intervention group for non-responders
#' @param tau follow-up
#' @param Delta_r RMST difference between intervention and control groups for responders
#' @param Delta_nr RMST difference between intervention and control groups for non-responders
#' @param Delta_0  RMST difference between responders and non-responders in the control group
#' @param anticipated_effects Logical parameter. If it is TRUE then the effect size is computed based on previous information on the effect sizes on response rate and survival-by-responses (that is, based on Delta_r, Delta_0, Delta_nr); otherwise is based on the distributional parameters (ascale0_r, ascale0_nr, ascale1_r, ascale1_nr, bshape0, bshape1).
#'
#'
#' @export
#'
#' @return Effect size for overall survival
#' @author Marta Bofill Roig

survm_effectsize <- function(ascale0_r,ascale0_nr,delta_p,p0,bshape0=1,bshape1=1,ascale1_r,ascale1_nr,tau,
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

  # output <- data.frame(Parameter=c("RMST difference","RMST difference responders","RMST difference non-responders","Response difference"),
  #                      Value=accounting(c(os_effect,Delta_r,Delta_nr,delta_p)))

  output <- data.frame(Parameter=c("RMST difference","RMST difference responders","RMST difference non-responders","Response difference"),
                       Value=c(os_effect,Delta_r,Delta_nr,delta_p))

  return(output)
}

##################################################################################


