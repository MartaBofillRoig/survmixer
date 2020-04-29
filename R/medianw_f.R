#' Median Weibull survival function
#'
#' @description The functions `meanw_f` and `medianw_f` calculate the mean and median for Weibull distributions, respectively.
#'
#'
#' @param lambda scale parameter for the Weibull distribution
#' @param bet shape parameter for the Weibull distribution
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
medianw_f <- function(lambda,bet){
  median = lambda*(log(2)^(1/bet))
  return(median)
}


#' Mean Weibull survival function
#'
#' @description The functions `meanw_f` and `medianw_f` calculate the mean and median for Weibull distributions, respectively.
#'
#'
#' @param lambda scale parameter for the Weibull distribution
#' @param bet shape parameter for the Weibull distribution
#'
#'
#' @export
#' @keywords internal
#' @return mean
#' @author Marta Bofill Roig
#'
meanw_f <- function(lambda,bet){
  mean = lambda*gamma(1+1/bet)
  return(mean)
}

