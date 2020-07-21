#' Derivative Weibull survival function
#'
#' @description The function `survw_derivf` computes the derivative of the survival distribution  `survw_f`.
#'
#'
#' @param t time
#' @param ascale scale parameter for the Weibull distribution
#' @param bshape shape parameter for the Weibull distribution
#'
#'
#' @export
#' @keywords internal
#' @return derivative
#' @author Marta Bofill Roig
#'
survw_derivf <- function(t,ascale,bshape=1){
  return(-bshape*(t/ascale)^bshape*exp(-(t/ascale)^bshape)/t)
}

