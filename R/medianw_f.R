#' Median Weibull survival function
#'
#' @description The functions `meanw_f` and `medianw_f` calculate the mean and median for Weibull distributions, respectively.
#'
#'
#' @param ascale scale parameter for the Weibull distribution
#' @param bshape shape parameter for the Weibull distribution
#'
#'
#' @export
#' @keywords internal
#' @return median
#' @author Marta Bofill Roig
#'
#' The functions `meanw_f` and `medianw_f` calculate the mean and median for Weibull distributions, respectively.
#'
#'
medianw_f <- function(ascale,bshape){
  median = ascale*(log(2)^(1/bshape))
  return(median)
}


#' Mean Weibull survival function
#'
#' @description The functions `meanw_f` and `medianw_f` calculate the mean and median for Weibull distributions, respectively.
#'
#'
#' @param ascale scale parameter for the Weibull distribution
#' @param bshape shape parameter for the Weibull distribution
#'
#'
#' @export
#' @keywords internal
#' @return mean
#' @author Marta Bofill Roig
#'
meanw_f <- function(ascale,bshape){
  mean = ascale*gamma(1+1/bshape)
  return(mean)
}

