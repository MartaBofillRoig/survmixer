#' Weibull survival function
#'
#' @description The function `survw_f` computes the Weibull survival function.
#'
#'
#' @param t time
#' @param ascale scale parameter for the Weibull distribution
#' @param bshape shape parameter for the Weibull distribution
#'
#'
#' @export
#' @keywords internal
#' @return survival function
#' @author Marta Bofill Roig
#'

#'
survw_f <- function(t,ascale,bshape){
  return(exp(-(t/ascale)^bshape))
}


