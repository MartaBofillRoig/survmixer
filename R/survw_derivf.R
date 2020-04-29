#' Derivative Weibull survival function
#'
#' @description The function `survw_derivf` computes the derivative of the survival distribution  `survw_f`.
#'
#'
#' @param t time
#' @param lambda scale parameter for the Weibull distribution
#' @param bet shape parameter for the Weibull distribution
#'
#'
#' @export
#' @keywords internal
#' @return derivative
#' @author Marta Bofill Roig
#'
survw_derivf <- function(t,lambda,bet=1){
  return(-bet*(t/lambda)^bet*exp(-(t/lambda)^bet)/t)
}

