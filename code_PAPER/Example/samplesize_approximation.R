


lambda1_f <- function(lambda0,delta,tau){
  c = delta + lambda0*(1-exp(-tau/lambda0))
  
  lambda1 = tau^2/(2*(tau-c))
  return(lambda1)
}

lambda1_f2 <- function(lambda0,delta,tau){
  c = delta + lambda0*(1-exp(-tau/lambda0))
  
  sol1_lambda1 =  (-3*tau^2 + sqrt(3)*sqrt(8*c*tau^3-5*tau^4))/(12*(c-tau))
  sol2_lambda1 =  (-3*tau^2 - sqrt(3)*sqrt(8*c*tau^3-5*tau^4))/(12*(c-tau))
  
  # sol1_lambda1= 2*tau^3/(3*tau^2 + sqrt(3)*sqrt(tau^3*(11*tau-8*c)))
  # sol2_lambda1= 2*tau^3/(3*tau^2 - sqrt(3)*sqrt(tau^3*(11*tau-8*c)))
  
  # sol1_lambda1 = (tau^2+sqrt(tau^3*(5*tau-4*c)))/(2*(c-tau))
  # sol2_lambda1 = (tau^2-sqrt(tau^3*(5*tau-4*c)))/(2*(c-tau))  
    
  return(list=c(sol1_lambda1,sol2_lambda1))
}


# taylor(f = exp, x0 = 0, n = 4)


lambda1_f(lambda0=wscale0_r, delta=Delta_r, tau=5)
lambda1_f2(lambda0=wscale0_r, delta=Delta_r, tau=5)
wscale1_r


lambda1_f(lambda0=wscale0_nr, delta=Delta_nr, tau=5)
lambda1_f2(lambda0=wscale0_nr, delta=Delta_nr, tau=5)
wscale1_nr


###########

survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, lambda1_r=wscale1_r,
                 lambda1_nr=wscale1_nr, lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2)


survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, 
                 lambda1_r=lambda1_f(lambda0=wscale0_r, delta=Delta_r, tau=5), 
                 lambda1_nr=lambda1_f(lambda0=wscale0_nr, delta=Delta_nr, tau=5), lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2)


survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, 
                 lambda1_r=lambda1_f2(lambda0=wscale0_r, delta=Delta_r, tau=5)[2], 
                 lambda1_nr=lambda1_f2(lambda0=wscale0_nr, delta=Delta_nr, tau=5)[2], lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2)

 

###########



survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, lambda1_r=wscale1_r,
                 lambda1_nr=wscale1_nr, lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2)


survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, 
                 lambda1_r=lambda1_f(lambda0=wscale0_r, delta=Delta_r, tau=5), 
                 lambda1_nr=wscale1_nr, lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2)


survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, 
                 lambda1_r=lambda1_f2(lambda0=wscale0_r, delta=Delta_r, tau=5)[2], 
                 lambda1_nr=wscale1_nr, lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2)


###########





survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, lambda1_r=wscale1_r,
                 lambda1_nr=wscale1_nr, lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2)


survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, 
                 lambda1_r=wscale0_r, 
                 lambda1_nr=lambda1_f(lambda0=wscale0_nr, delta=Delta_nr, tau=5), lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2)


survw_samplesize(lambda0_r=wscale0_r, lambda0_nr=wscale0_nr, delta_p=p1-p0, p0=p0, beta0=1, beta1=1, 
                 lambda1_r=wscale0_r, 
                 lambda1_nr=lambda1_f2(lambda0=wscale0_nr, delta=Delta_nr, tau=5)[2], lambda_cens = 1/lambda_c, tau=5, alpha=0.05, beta=0.2)


