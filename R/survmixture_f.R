#' Mixture survival function
#'
#' @description The function `survmixture_f` computes the survival distribution as a mixture of  responders and non-responders. The responders and non-responders distributions are assumed to be Weibull distributions.
#'
#'
#' @param t time at which the survival distribution is evaluated
#' @param ascale_r scale parameter for the Weibull distribution   for responders
#' @param ascale_nr scale parameter for the Weibull distribution  for non-responders
#' @param bshape shape parameter for the Weibull distribution
#' @param p event rate for the response
#'
#' @export
#'
#' @return This function returns the survival function evaluated at t based on a  mixture model of  responders and non-responders.
#' @author Marta Bofill Roig.
#' @references Design of phase III trials with long-term survival outcomes based on short-term binary results. Marta Bofill Roig, Yu Shen, Guadalupe Gomez Melis. 	arXiv:2008.12887
#'
#'
survmixture_f <- function(t,ascale_r, ascale_nr, bshape=1, p){

  if( 0 > p || p > 1){
    stop("Response probability must be a number between 0 and 1")
  }
  # if(t<0){
  #   stop("The time must be a positive number")
  # }
  if(ascale_r<0 || ascale_nr<0){
    stop("Scale parameter must be a positive number")
  }
  if(bshape<0 ){
    stop("Shape parameter must be a positive number")
  }

  s <- survw_f(t,ascale_r,bshape)*p + survw_f(t,ascale_nr,bshape)*(1-p)
  return(s)
}
