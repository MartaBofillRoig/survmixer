#' Restricted mean survival times Weibull distribution
#'
#' @description The function `rmstw_f` computes the restricted mean survival times (RMST) according to the Weibull survival function.
#'
#'
#' @param low RMST evaluated from low to tau
#' @param tau RMST evaluated from low to tau
#' @param ascale scale parameter for the Weibull distribution
#' @param bshape shape parameter for the Weibull distribution
#'
#'
#' @export
#' @keywords internal
#' @return rmst
#'
#' @author Marta Bofill Roig
#'

#'
rmstw_f <- function(ascale,bshape,tau,low=0){
  requireNamespace("stats")
  r <- integrate(survw_f, lower = low, upper = tau, ascale, bshape)$value
  return(r)
}

