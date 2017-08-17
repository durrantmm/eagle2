
# Sampler should be a stanfit object. This creates a skeleton list of
# parameters with the correct dimensionality. 
get_skeleton=function(sampler) {
  m_pars <- sampler@model_pars
  idx_wo_lp <- which(m_pars != "lp__")
  m_pars <- m_pars[idx_wo_lp]
  p_dims <- sampler@par_dims[idx_wo_lp]
  rstan:::create_skeleton(m_pars, p_dims)
}

# Adagrad stochastic optimizer
adagrad=function(grad, x, master_stepsize=0.1, eps=1e-6, iterations=300, verbose=F) {
  historical_grad=x * 0.
  #remember_g2=x * 0.
  progress=list()
  for (i in 1:iterations) {
    g=grad(x)
    historical_grad=historical_grad+g^2
    #remember_g2=0.99 * remember_g2 + 0.01 * g
    x=x-master_stepsize*g/(eps+sqrt(historical_grad))
    progress[[i]]=attr(g,"log_prob")
    if (verbose) cat(i, attr(g,"log_prob"), "\n")
  }
  list(x=x,log_prob=unlist(progress))
}

robust_adagrad=function(...) {
  tryCatch( adagrad(...), error=function(g) NULL )
}

# pretty unstable
adamax=function(grad, x, master_stepsize=0.1, b1=0.1, b2=0.001, iterations=300, verbose=F) {
  momentum=x * 0.
  inf_norm=x * 0.
  progress=list()
  for (i in 1:iterations) {
    g=grad(x)
    momentum=(1. - b1) * momentum + b1 * g
    inf_norm=max( b2 * inf_norm, abs(g) + 1e-8 )
    x=x - master_stepsize * momentum / inf_norm
    progress[[i]]=attr(g,"log_prob")
    if (verbose) cat(i, attr(g,"log_prob"), "\n")
  }
  list(x=x,log_prob=unlist(progress))
}

# adadelta: works but goes very slowly
adadelta=function(grad, x, master_stepsize=1.0, oneMrho=0.05, eps=1e-8, iterations=300, verbose=F) {
  Eg2=x * 0.
  Edx2=x * 0.
  progress=list()
  for (i in 1:iterations) {
    g=grad(x)
    Eg2=(1. - oneMrho) * Eg2 + oneMrho * g^2
    deltax=sqrt(Edx2)/sqrt(Eg2 + eps) * g
    Edx2=(1. - oneMrho) * Edx2 + oneMrho * deltax^2
    x=x - master_stepsize * deltax
    progress[[i]]=attr(g,"log_prob")
    if (verbose) cat(i, attr(g,"log_prob"), "\n")
  }
  list(x=x,log_prob=unlist(progress))
}

#' Stochastic variational EM. 
#' 
#' Uses Stochastic Variational Inference (SVI) over parameters where to_optim==F and regular Stochastic Gradient (well Adagrad) over parameters for which to_optim==T. 
#' 
#' @param grad_log_prob Function which takes unconstrained parameters and returns the gradient of the log joint distribution. 
#' @param to_optim Logical vector denoting which parameters should be optimized over rather than integrated over using SVI. 
#' @param init Initialization. 
#' @param plot.elbo Whether to plot the progress in term of the ELBO (Evidence Lower Bound)
#' @param log_prob (optional) Function returning log joint distribution. 
#' 
#' @return List with
#' \item{m}{Mean of approximate posterior for each parameter (or point estimate)}
#' \item{s}{Standard deviation of approximate posterior for each param, or 0 or optimized params}
#' \item{elbo_func}{Function to calculate (noisy) estimate of ELBO using a single sample from rnorm(sum(!to_optim)).}
#' \item{elbo_progress}{Vector of estimated ELBO with algorithm iteration.}
svem=function(grad_log_prob, to_optim, init=NULL, plot.elbo=F, log_prob=NULL, master_stepsize=1., ...) {

  #to_optim=c(F,T,T)
  noptim=sum(to_optim)
  nintegrate=sum(!to_optim)
  P=length(to_optim)
  
  elbo_grad_mixed=function(temp0, eta=rnorm(nintegrate)) {
    x=numeric(P)
    x[to_optim]=temp0[1:noptim]
    temp=temp0[(noptim+1):length(temp0)]
    m=temp[1:nintegrate]
    logs=temp[(nintegrate+1):(2*nintegrate)]
    s=exp(logs)
    x[!to_optim]=m+s*eta
    g=grad_log_prob(x)
    stopifnot(!any(is.na(g)))
    grad_m=g[!to_optim]
    grad_logs=grad_m*eta*s+1 # plus 1 from entropy
    res=-c(g[to_optim],grad_m,grad_logs)
    attr(res,"log_prob")=attr(g,"log_prob") + sum(logs)
    res
  }
  
  if (is.null(init))
    init=numeric(noptim+2*nintegrate) else
      if (is.list(init)) 
        init=c( init$m[to_optim], init$m[!to_optim], log(init$s[!to_optim]))
  
  stopifnot(length(init)==noptim+2*nintegrate)
  
  adagrad_fit=robust_adagrad(elbo_grad_mixed, init, master_stepsize, ...) 
  while(is.null(adagrad_fit)) {
    master_stepsize=master_stepsize/2
    cat("Reducing learning rate to",master_stepsize,"\n")
    adagrad_fit=robust_adagrad(elbo_grad_mixed, init, master_stepsize, ...)
    if (master_stepsize < 1e-12) return(NULL)
  }
  
  elbo_func=function(eta=rnorm(nintegrate)) attr( elbo_grad_mixed(adagrad_fit$x, eta), "log_prob" )
  
  if (plot.elbo) {
    require(stats)
    ma <- function(x,n=5){stats::filter(x,rep(1/n,n), sides=2)}
    elbo_progression=adagrad_fit$log_prob
    elbo_ma=ma(elbo_progression,100)
    plot(elbo_progression,pch=16,col=rgb(.5,.5,.5,.5), ylim=c(min(elbo_progression[100:length(elbo_progression)]),max(elbo_progression)))
    lines(elbo_ma, col="red", lwd=2)
  }
  
  temp0=adagrad_fit$x
  
  if (!is.null(log_prob)) {
    x=numeric(P)
    x[to_optim]=temp0[1:noptim]
    temp=temp0[(noptim+1):length(temp0)]
    mh=temp[1:nintegrate]
    logs=temp[(nintegrate+1):(2*nintegrate)]
    sh=exp(logs)
    sumlogs=sum(logs)
    elbo_func=function(eta=rnorm(nintegrate)) {
      x[!to_optim]=mh+sh*eta
      log_prob(x) + sumlogs
    }
  }
  
  m=numeric(P)
  s=numeric(P)
  m[to_optim]=temp0[1:noptim]
  s[to_optim]=0
  temp=temp0[(noptim+1):length(temp0)]
  m[!to_optim]=temp[1:nintegrate]
  s[!to_optim]=exp(temp[(nintegrate+1):(2*nintegrate)])
  
  list(m=m, s=s, elbo_func= elbo_func, elbo_progress=adagrad_fit$log_prob)
}


svem_lbfgs=function(log_joint, grad_log_joint, to_optim, init=NULL, samples=100, trace=0, iterations=1000,factr=1e7) {
  
  #to_optim=c(F,T,T)
  noptim=sum(to_optim)
  nintegrate=sum(!to_optim)
  P=length(to_optim)
  
  eta_fixed=foreach(i=1:samples) %do% { rnorm(nintegrate) }
  
  elbo_grad_mixed=function(temp0) {
    x=numeric(P)
    x[to_optim]=temp0[1:noptim]
    temp=temp0[(noptim+1):length(temp0)]
    m=temp[1:nintegrate]
    logs=temp[(nintegrate+1):(2*nintegrate)]
    s=exp(logs)
    foreach(eta=eta_fixed, .combine = rbind) %do% { # could parallelize
      x[!to_optim]=m+s*eta
      g=grad_log_joint(x)
      stopifnot(!any(is.na(g)))
      grad_m=g[!to_optim]
      grad_logs=grad_m*eta*s+1 # plus 1 from entropy
      c(g[to_optim],grad_m,grad_logs)
    } %>% colMeans()
  }
  
  elbo_mixed=function(temp0) {
    x=numeric(P)
    x[to_optim]=temp0[1:noptim]
    temp=temp0[(noptim+1):length(temp0)]
    m=temp[1:nintegrate]
    logs=temp[(nintegrate+1):(2*nintegrate)]
    s=exp(logs)
    mean( foreach(eta=eta_fixed, .combine=c) %do% { # could parallelize
      x[!to_optim]=m+s*eta
      log_joint(x)
    } )  + sum(logs)
  }
  
  if (is.null(init))
    init=numeric(noptim+2*nintegrate) else
      if (is.list(init)) 
        init=c( init$m[to_optim], init$m[!to_optim], log(init$s[!to_optim]))
  
  stopifnot(length(init)==noptim+2*nintegrate)
  a=elbo_grad_mixed(init) 
 # stop(1)
  fit=optim(init, elbo_mixed, elbo_grad_mixed, method="L-BFGS-B", control=list(fnscale=-1, trace=trace,maxit=iterations,factr=factr) )
  
  temp0=fit$par
  
  m=numeric(P)
  s=numeric(P)
  m[to_optim]=temp0[1:noptim]
  s[to_optim]=0
  temp=temp0[(noptim+1):length(temp0)]
  m[!to_optim]=temp[1:nintegrate]
  s[!to_optim]=exp(temp[(nintegrate+1):(2*nintegrate)])
  
  list(m=m, s=s, elbo_opt=fit$value, elbo_mixed=function(x) { 
    elbo_mixed( c( x$m[to_optim], x$m[!to_optim], log(x$s[!to_optim])) )
  })
}
