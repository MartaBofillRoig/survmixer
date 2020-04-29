#' Restricted mean survival times Weibull distribution
#'
#' @description The function `rmstw_f` computes the restricted mean survival times (RMST) according to the Weibull survival function.
#'
#'
#' @param low rmst evaluated from low to tau
#' @param tau rmst evaluated from low to tau
#' @param lambda scale parameter for the Weibull distribution
#' @param bet shape parameter for the Weibull distribution
#'
#'
#' @export
#' @keywords internal
#' @return rmst
#' @author Marta Bofill Roig
#'

#'
rmstw_f <- function(lambda,bet,tau,low=0){
  r <- integrate(survw_f, lower = low, upper = tau, lambda, bet)$value
  return(r)
}

