#' Scale parameter computation
#'
#' @description returns the value of the scale parameter in the intervention group using Taylor series
#'
#'
#' @param ascale0 scale parameter for the weibull distribution in the control group
#' @param Delta RMST difference between groups
#' @param tau end of follow-up
#'
#'
#' @export
#' @keywords internal
#' @note Weibull parametrization: S(x) =  exp(- (x/a)^b)
#' @return Variance computation
#' @author Marta Bofill Roig
#'
#'
#' 
 
scale1_taylorf <- function(ascale0,Delta,tau){
  c = Delta + ascale0*(1-exp(-tau/ascale0))
  
  sol1_ascale1 =  (-3*tau^2 + sqrt(3)*sqrt(8*c*tau^3-5*tau^4))/(12*(c-tau))
  sol2_ascale1 =  (-3*tau^2 - sqrt(3)*sqrt(8*c*tau^3-5*tau^4))/(12*(c-tau)) 
  
  return(list=c(sol1_ascale1,sol2_ascale1))
}