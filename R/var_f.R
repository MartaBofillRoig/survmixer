#' Inside variance computation
#'
#' @description The following three functions are used to calculate the variance of th difference of two RMSTs. `survw_integratef` is used for the integrations; `inside_var` calculates the expression inside the integral; finally, `var_f` computes the variance.
#'
#'
#' @param t time at which the survival distribution is evaluated
#' @param ascale_r scale parameter for the Weibull distribution   for responders
#' @param ascale_nr scale parameter for the Weibull distribution  for non-responders
#' @param ascale_cens distributional parameter for the exponential distribution for the censoring
#' @param bshape shape parameter for the Weibull distribution
#' @param p event rate for the response
#' @param tau follow-up
#'
#'
#' @export
#' @keywords internal
#' @return Variance computation
#' @author Marta Bofill Roig
#'
#'
#'
inside_var <- function(t,ascale_r,ascale_nr,tau,bshape,ascale_cens,p){

  num = p*sapply(t,survw_integratef,ascale=ascale_r,bshape=bshape,tau=tau)+(1-p)*sapply(t,survw_integratef,ascale=ascale_nr,bshape=bshape,tau=tau)
  den = survmixture_f(t,ascale_r, ascale_nr, bshape, p)
  dervS = p*survw_derivf(t,ascale_r,bshape) + (1-p)*survw_derivf(t,ascale_nr,bshape)

  inside_integral <- (num/den)^2*(1/survw_f(t,ascale_cens,bshape=1))*dervS

  return(-inside_integral)
}


#' Variance computation
#' @description The   function  `var_f` computes the variance.
#'
#'
#' @param ascale_r scale parameter for the Weibull distribution   for responders
#' @param ascale_nr scale parameter for the Weibull distribution  for non-responders
#' @param bshape shape parameter for the Weibull distribution
#' @param p event rate for the response
#' @param tau follow-up
#' @param ascale_cens distributional parameter for the exponential distribution for the censoring
#'
#'
#' @export
#' @keywords internal
#' @return Variance computation
#' @author Marta Bofill Roig
#'
#'
#'
var_f <- function(ascale_r,ascale_nr,tau,bshape,ascale_cens,p){
  integrate(inside_var,lower=0,upper=tau,ascale_r=ascale_r, ascale_nr= ascale_nr,tau=tau,bshape=bshape,ascale_cens=ascale_cens,p=p)$value
}

#' Integrate function
#' @description the function `survw_integratef` is used for the integrations
#'
#'
#' @param ascale scale parameter for the Weibull distribution
#' @param bshape shape parameter for the Weibull distribution
#' @param tau follow-up
#' @param t time
#'
#'
#' @export
#' @keywords internal
#' @return Variance computation
#' @author Marta Bofill Roig
#'
#'
#'
survw_integratef <- function(t,tau, ascale,bshape){
  int <- integrate(survw_f,lower=t, upper=tau,ascale,bshape)$value
  return(int)
}
