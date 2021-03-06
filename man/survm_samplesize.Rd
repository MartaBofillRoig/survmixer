% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survm_samplesize.R
\name{survm_samplesize}
\alias{survm_samplesize}
\title{Sample size calculation for mixture survival distributions}
\usage{
survm_samplesize(
  ascale0_r,
  ascale0_nr,
  ascale1_r,
  ascale1_nr,
  delta_p,
  p0,
  m0_r,
  m0_nr,
  diffm_r,
  diffm_nr,
  S0_r,
  S0_nr,
  diffS_r,
  diffS_nr,
  Delta_r,
  Delta_nr,
  ascale_cens,
  tau,
  bshape0 = 1,
  bshape1 = 1,
  all_ratio = 0.5,
  alpha = 0.025,
  beta = 0.2,
  set_param = 0
)
}
\arguments{
\item{ascale0_r}{scale parameter for the Weibull distribution in the control group for responders}

\item{ascale0_nr}{scale parameter for the Weibull distribution in the control group for non-responders}

\item{ascale1_r}{scale parameter for the Weibull distribution in the intervention group for responders}

\item{ascale1_nr}{scale parameter for the Weibull distribution in the intervention group for non-responders}

\item{delta_p}{effect size for the response rate}

\item{p0}{event rate for the response}

\item{m0_r}{survival mean for responders in the control group}

\item{m0_nr}{survival mean for non-responders in the control group}

\item{diffm_r}{difference in survival means between groups for responders}

\item{diffm_nr}{difference in survival means between groups for responders}

\item{S0_r}{tau-year survival rates for responders in the control group}

\item{S0_nr}{tau-year survival rates for non-responders in the control group}

\item{diffS_r}{difference in tau-year survival rates for responders}

\item{diffS_nr}{difference in tau-year survival rates for non-responders}

\item{Delta_r}{restricted mean survival times (RMST) difference between intervention and control groups for responders}

\item{Delta_nr}{RMST difference between intervention and control groups for non-responders}

\item{ascale_cens}{distributional parameter for the exponential distribution for the censoring}

\item{tau}{follow-up}

\item{bshape0}{shape parameter for the Weibull distribution in the control group}

\item{bshape1}{shape parameter for the Weibull distribution in the intervention group}

\item{all_ratio}{allocation ratio. The ratio of numbers of participants allocated in the control group. By default is assumed 1:1 (i.e., all_ratio=0.5)}

\item{alpha}{type I error}

\item{beta}{type II error}

\item{set_param}{Set of parameters to be used for the responders/non-responders survival functions If the set of parameters is =1, then the sample size is computed using the survival means (m0_r,m0_nr,diffm _r,diffm_nr); if set_param=2, it is computed using the tau-year survival rates (S0_r,S0_nr,diffS_r,diffS_nr); if set_param=2, it is computed using the RMSTs and survival rates (Delta_r,Delta_nr,S0_r,S0_nr). If set_param=0, the computation is based on the distributional parameters (ascale0_r, ascale0_nr, ascale1_r, ascale1_nr).}
}
\value{
This function returns the total sample size needed and the expected effect size for overall   survival  (RMST difference between groups).
}
\description{
The function `survm_samplesize` calculates the sample size according to the distributional parameters of the responders and non-responders.
}
\references{
Design of phase III trials with long-term survival outcomes based on short-term binary results. Marta Bofill Roig, Yu Shen, Guadalupe Gomez Melis. 	arXiv:2008.12887
}
\author{
Marta Bofill Roig.
}
