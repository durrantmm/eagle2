
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
# Would rather use AdaMax instead? 
adagrad=function(grad, x, master_stepsize=0.1, eps=1e-6, iterations=300, verbose=F) {
  historical_grad=0
  progress=list()
  for (i in 1:iterations) {
    g=grad(x)
    historical_grad=historical_grad+g^2
    x=x-master_stepsize*g/(eps+sqrt(historical_grad))
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
#' @param samples_for_elbo How many samples to use to calculate the final ELBO. 
#' @param log_prob (optional) Function returning log joint distribution. 
#' 
#' @return List with
#' \item{m}{Mean of approximate posterior for each parameter (or point estimate)}
#' \item{s}{Standard deviation of approximate posterior for each param, or 0 or optimized params}
#' \item{elbo_func}{Function to calculate (noisy) estimate of ELBO using a single sample from rnorm(sum(!to_optim)).}
#' \item{elbo_progress}{Vector of estimated ELBO with algorithm iteration.}
svem=function(grad_log_prob, to_optim, init=NULL, plot.elbo=F, samples_for_elbo=10000, log_prob=NULL, ...) {

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
  
  adagrad_fit=adagrad(elbo_grad_mixed, init, ...) 
  
  elbo_func=function(eta=rnorm(nintegrate)) attr( elbo_grad_mixed(adagrad_fit$x, eta), "log_prob" )
  
  elbo_estimate=if (samples_for_elbo>=1) mean(unlist(foreach(i=1:samples_for_elbo) %do% { elbo_func() }), na.rm=T) else NA
  
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
  
  list(m=m, s=s, elbo=elbo_estimate, elbo_func= elbo_func, elbo_progress=adagrad_fit$log_prob)
}
