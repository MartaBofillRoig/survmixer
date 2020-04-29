#' Weibull survival function
#'
#' @description The function `survw_f` computes the Weibull survival function.
#'
#'
#' @param t time
#' @param lambda scale parameter for the Weibull distribution
#' @param bet shape parameter for the Weibull distribution
#'
#'
#' @export
#' @keywords internal
#' @return survival function
#' @author Marta Bofill Roig
#'

#'
survw_f <- function(t,lambda,bet){
  return(exp(-(t/lambda)^bet))
}


