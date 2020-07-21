#' Scale parameter computation
#'
#' @description returns the value of the scale parameter a given the survival (s) at time t
#'
#'
#' @param t time at which the survival distribution is evaluated
#' @param s survival rate at time t 
#' @param shape shape parameter for the Weibull distribution 
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

# param_scale: returns the value of the scale parameter a given the survival (s) at time t
param_scale <- function(s,t,shape=1){
  scale = -t/((log(s))^(1/shape))
  return(scale)
}

