#' Inside variance computation
#'
#' @description The following three functions are used to calculate the variance of th difference of two RMSTs. `survw_integratef` is used for the integrations; `inside_var` calculates the expression inside the integral; finally, `var_f` computes the variance.
#'
#'
#' @param t time at which the survival distribution is evaluated
#' @param lambda_r scale parameter for the Weibull distribution   for responders
#' @param lambda_nr scale parameter for the Weibull distribution  for non-responders
#' @param lambda_cens distributional parameter for the exponential distribution for the censoring
#' @param bet shape parameter for the Weibull distribution
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
inside_var <- function(t,lambda_r,lambda_nr,tau,bet,lambda_cens,p){

  num = p*sapply(t,survw_integratef,lambda=lambda_r,bet=bet,tau=tau)+(1-p)*sapply(t,survw_integratef,lambda=lambda_nr,bet=bet,tau=tau)
  den = survmixture_f(t,lambda_r, lambda_nr, bet, p)
  dervS = p*survw_derivf(t,lambda_r,bet) + (1-p)*survw_derivf(t,lambda_nr,bet)

  inside_integral <- (num/den)^2*(1/survw_f(t,lambda_cens,bet=1))*dervS

  return(-inside_integral)
}


#' Variance computation
#' @description The   function  `var_f` computes the variance.
#'
#'
#' @param lambda_r scale parameter for the Weibull distribution   for responders
#' @param lambda_nr scale parameter for the Weibull distribution  for non-responders
#' @param bet shape parameter for the Weibull distribution
#' @param p event rate for the response
#' @param tau follow-up
#' @param lambda_cens distributional parameter for the exponential distribution for the censoring
#'
#'
#' @export
#' @keywords internal
#' @return Variance computation
#' @author Marta Bofill Roig
#'
#'
#'
var_f <- function(lambda_r,lambda_nr,tau,bet,lambda_cens,p){
  integrate(inside_var,lower=0,upper=tau,lambda_r=lambda_r, lambda_nr= lambda_nr,tau=tau,bet=bet,lambda_cens=lambda_cens,p=p)$value
}

#' Integrate function
#' @description the function `survw_integratef` is used for the integrations
#'
#'
#' @param lambda scale parameter for the Weibull distribution
#' @param bet shape parameter for the Weibull distribution
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
survw_integratef <- function(t,tau, lambda,bet){
  int <- integrate(survw_f,lower=t, upper=tau,lambda,bet)$value
  return(int)
}
