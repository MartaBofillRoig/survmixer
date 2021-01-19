#' Sample size calculation for mixture survival distributions
#'
#' @description The function `survm_samplesize` calculates the sample size according to the distributional parameters of the responders and non-responders.
#'
#' @param ascale0_r scale parameter for the Weibull distribution in the control group for responders
#' @param ascale0_nr scale parameter for the Weibull distribution in the control group for non-responders
#' @param ascale1_r scale parameter for the Weibull distribution in the intervention group for responders
#' @param ascale1_nr scale parameter for the Weibull distribution in the intervention group for non-responders
#' @param delta_p effect size for the response rate
#' @param p0 event rate for the response
#' @param m0_r survival mean for responders in the control group
#' @param m0_nr survival mean for non-responders in the control group
#' @param diffm_r difference in survival means between groups for responders
#' @param diffm_nr difference in survival means between groups for responders
#' @param S0_r tau-year survival rates for responders in the control group
#' @param S0_nr tau-year survival rates for non-responders in the control group
#' @param diffS_r difference in tau-year survival rates for responders
#' @param diffS_nr difference in tau-year survival rates for non-responders
#' @param Delta_r restricted mean survival times (RMST) difference between intervention and control groups for responders
#' @param Delta_nr RMST difference between intervention and control groups for non-responders
#' @param ascale_cens distributional parameter for the exponential distribution for the censoring
#' @param tau follow-up
#' @param bshape0 shape parameter for the Weibull distribution in the control group
#' @param bshape1 shape parameter for the Weibull distribution in the intervention group
#' @param all_ratio allocation ratio. The ratio of numbers of participants allocated in the control group. By default is assumed 1:1 (i.e., all_ratio=0.5)
#' @param alpha type I error
#' @param beta type II error
#' @param set_param Set of parameters to be used for the responders/non-responders survival functions If the set of parameters is =1, then the sample size is computed using the survival means (m0_r,m0_nr,diffm _r,diffm_nr); if set_param=2, it is computed using the tau-year survival rates (S0_r,S0_nr,diffS_r,diffS_nr); if set_param=2, it is computed using the RMSTs and survival rates (Delta_r,Delta_nr,S0_r,S0_nr). If set_param=0, the computation is based on the distributional parameters (ascale0_r, ascale0_nr, ascale1_r, ascale1_nr).
#'
#' @export
#' @import stats
#' @example survm_samplesize(S0_r=0.55,S0_nr=0.41,diffS_r=0.32,diffS_nr=0,p0=0.19,delta_p=0.19,ascale_cens=7,tau=5,alpha=0.05,beta=0.2,set_param=2)
#' @return This function returns the total sample size needed and the expected effect size for overall   survival  (RMST difference between groups).
#' @author Marta Bofill Roig.
#' @references Design of phase III trials with long-term survival outcomes based on short-term binary results. Marta Bofill Roig, Yu Shen, Guadalupe Gomez Melis. 	arXiv:2008.12887

survm_samplesize <- function(ascale0_r,ascale0_nr,ascale1_r,ascale1_nr,delta_p,p0,
                             m0_r, m0_nr, diffm_r, diffm_nr,
                             S0_r, S0_nr, diffS_r, diffS_nr,
                             Delta_r, Delta_nr,
                             ascale_cens,tau,
                             bshape0=1,bshape1=1,
                             all_ratio=0.5,alpha=0.025,beta=0.2,
                             set_param=0){

  requireNamespace("stats")
  if( 0 > p0 || p0 > 1){
    stop("Response probability must be a number between 0 and 1")
  }
  if( 0 > delta_p){
    stop("Effect size for the response rate must be a positive number")
  }
  if(ascale_cens<0){
    stop("Scale parameter must be a positive number")
  }
  # set of parameters
  if(set_param==0){
    if(bshape0<0 || bshape0<0){
      stop("Shape parameter must be a positive number")
    }else if(ascale0_r<0 || ascale0_nr<0 || ascale1_r<0 || ascale1_nr<0){
      stop("Scale parameter must be a positive number")
    }
  }
  if(set_param==1){
    if(bshape0<0 || bshape0<0){
      stop("Shape parameter must be a positive number")
    }else if(m0_r<0 || m0_nr<0){
      stop("Survival means must be a positive number")
    }else if(diffm_r<0 || diffm_nr<0){
      stop("Difference in survival means must be a positive number")
    }
  }
  if(set_param==2){
    if(bshape0<0 || bshape0<0){
      stop("Shape parameter must be a positive number")
    }else if(S0_r<0 || S0_nr<0){
      stop("tau-year survival rates must be a positive number")
    }else if(diffS_r<0 || diffS_nr<0){
      stop("Difference in tau-year survival rates must be a positive number")
    }
  }
  if(set_param==3){
    if(bshape0<0 || bshape0<0){
      stop("Shape parameter must be a positive number")
    }else if(S0_r<0 || S0_nr<0){
      stop("tau-year survival rates must be a positive number")
    }else if(Delta_r<0 || Delta_nr<0){
      stop("RMST difference between intervention and control groups must be a positive number")
    }
  }
  if(set_param>3||set_param<0){
    stop("The selected set of parameters is not valid")
  }
  if(tau<0){
    stop("The follow-up must be a positive number")
  }else if(0 > all_ratio || all_ratio > 1){
    stop("Allocation ratio must be number between 0 and 1")
  }else if( 0 > alpha || alpha > 1){
    stop("Alpha must be number between 0 and 1")
  }else if( 0 > beta || beta > 1){
    stop("Beta must be number between 0 and 1")
  }


  z_alpha <- qnorm(1-alpha,0,1)
  z_beta <-  qnorm(1-beta,0,1)
  p1 = delta_p +  p0

  if(set_param==1){

    m1_r = diffm_r+m0_r
    m1_nr = diffm_nr+m0_nr

    # note: mean = ascale*gamma(1+1/bshape)
    ascale0_r = m0_r/gamma(1+1/bshape0)
    ascale0_nr = m0_nr/gamma(1+1/bshape0)
    ascale1_r = m1_r/gamma(1+1/bshape1)
    ascale1_nr = m1_nr/gamma(1+1/bshape1)

  }

  if(set_param==2){

    S1_r = diffS_r+ S0_r
    S1_nr = diffS_nr + S0_nr

    ascale0_r = param_scale(s=S0_r,t=tau,shape=bshape0)
    ascale0_nr = param_scale(s=S0_nr,t=tau,shape=bshape0)
    ascale1_r = param_scale(s=S1_r,t=tau,shape=bshape1)
    ascale1_nr = param_scale(s=S1_nr,t=tau,shape=bshape1)

  }

  if(set_param==3){

    ascale0_r = param_scale(s=S0_r,t=tau,shape=bshape0)
    ascale0_nr = param_scale(s=S0_nr,t=tau,shape=bshape0)

    if(bshape0==1 && bshape1==1){
      ascale1_r = scale1_taylorf(ascale0=ascale0_r,Delta=Delta_r,tau=tau)[2]
      ascale1_nr = scale1_taylorf(ascale0=ascale0_nr,Delta=Delta_nr,tau=tau)[2]
    }
  }

  os_effect = survm_effectsize(ascale0_r,ascale0_nr,delta_p,p0,bshape0,bshape1,ascale1_r,ascale1_nr,tau)$Value[1]

  var0 <- var_f(ascale_r=ascale0_r,ascale_nr=ascale0_nr,tau=tau,bshape=bshape0,ascale_cens=ascale_cens,p=p0)
  var1 <- var_f(ascale_r=ascale1_r,ascale_nr=ascale1_nr,tau=tau,bshape=bshape1,ascale_cens=ascale_cens,p=p1)
  ss = ((z_alpha+z_beta)/(os_effect))^2*(var0/all_ratio + var1/(1-all_ratio))

  # output <- list(samplesize=ss,effectsize=os_effect)

  output <- data.frame(Parameter=c("Sample size","RMST difference"),
                   Value=c(ss, os_effect))

  return(output)
}

##################################################################################


