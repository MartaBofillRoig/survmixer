#' Mixture survival function
#'
#' @description The function `survmixture_f` computes the survival distribution as a mixture of  of responders and non-responders. The responders and non-responders distributions are assumed to be Weibull distributions.
#'
#'
#' @param t time at which the survival distribution is evaluated
#' @param lambda_r scale parameter for the Weibull distribution   for responders
#' @param lambda_nr scale parameter for the Weibull distribution  for non-responders
#' @param bet shape parameter for the Weibull distribution
#' @param p event rate for the response
#'
#' @export
#'
#' @return Mixture survival function evaluated at t
#' @author Marta Bofill Roig
#'
#'
survmixture_f <- function(t,lambda_r, lambda_nr, bet=1, p){
  s <- survw_f(t,lambda_r,bet)*p + survw_f(t,lambda_nr,bet)*(1-p)
  return(s)
}
