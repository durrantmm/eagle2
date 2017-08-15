

logit=function(g) log(g/(1.0-g))
inv_logit=function(g) { 1.0/(1.0+exp(-g)) }

#' Beta binomial GLMM with flips. Prior on concentration parameter is Gamma(concShape,concRate)
#'
#' Includes a per individual, per SNP random effect (shared across conditions) and uses stochastic 
#' variational inference to integrate over these.
#' 
#' @param ys numerator counts [n x T x K] where n are individuals, T are timepoints, K are SNPs
#' @param ns denominator counts [n x T x K]
#' @param concShape Shape of prior on concentration
#' @param concRate Rate of prior on concentration
#' @return List with likelihood ratio, p-value and fits
#' @importFrom rstan optimizing
#' @importFrom foreach foreach %do%
#' @importFrom abind abind
#' @export
eagle2_re=function(ys,ns,concShape=1.0001,concRate=1e-4,elbo_samples=3000,...) {
  
  N=dim(ys)[1]
  Ti=dim(ys)[2]
  K=dim(ys)[3]

  mode(ns)="integer"
  mode(ys)="integer"
  
  xNull=aperm( abind( foreach(i=seq_len(Ti)) %do% matrix(1,N,1), along=3 ), c(3,1,2) )
  
  dat=list(N=N,P=1L,T=Ti,K=K,ys=ys,ns=ns,x=xNull,concShape=concShape,concRate=concRate)
  
  # get stan_fit object so we can call sampler$grad_log_prob for the gradient
  # no sampling is done, just a hack to get a stanfit object
  sampler=sampling(stanmodels$bb_with_re, dat, iter=1, chains=0)
  
  # specific which parameters to optimize (rather than integrate over)
  sk=get_skeleton(sampler)
  sk$conc[]=T
  sk$p[]=T
  sk$log_rev=T
  sk$beta[]=T
  sk$re[]=F
  to_optim=as.logical(unlist(sk))
  
  # initialize using the model fit without the random effects
  o=optimizing(stanmodels$bb, dat, as_vector=F)
  init=list(m=get_skeleton(sampler), s=get_skeleton(sampler))
  init$m$beta=o$par$beta
  #init$m$p=logit(o$par$p)
  #init$m$p=pmin( pmax(init$m$p,-4), 4 )
  init$m$p=logit(o$par$p) 
  init$m$conc=log(o$par$conc)
  #init$m$conc=o$par$conc
  for (n in names(init$s)) init$s[[n]][]=1 # set all standard deviations to 1 initially
  for (n in names(init)) init[[n]]=unlist(init[[n]])
  
  # set the seed so we can use the same seed for the alternative model
  set.seed(1)
  
  null_gradient_function=function(g) grad_log_prob(sampler,g,T)
  #null_gradient_function(init$m)
  null_likelihood=function(g) log_prob(sampler,g,T,F)
  v_null=svem(null_gradient_function, to_optim, init, plot.elbo = F, samples_for_elbo = 0, log_prob = null_likelihood, ...)
  
  # Fit alternative model
  dat_full=dat
  temp=matrix(0,N,Ti)
  temp[,1]=1
  dat_full$x=aperm( abind( foreach(i=seq_len(Ti)) %do% { temp[,i]=1; temp }, along=3 ), c(3,1,2) )
  
  dat_full$P=Ti
  
  # specify which parameters to optimize
  sampler_full=sampling(stanmodels$bb_with_re, dat_full, iter=1, chains=0)
  sk_full=get_skeleton(sampler_full)
  sk_full$conc[]=T
  sk_full$p[]=T
  sk_full$log_rev=T
  sk_full$beta[]=T
  sk_full$re[]=F
  to_optim=as.logical(unlist(sk_full))
  
  # initialization choices: initializing from null fit seems best
  init=list( m=rstan:::rstan_relist(v_null$m, sk), 
    s=rstan:::rstan_relist(v_null$s, sk) )
  init$m$beta=c(init$m$beta, numeric(Ti-1) )
  init$s$beta=c(init$s$beta, numeric(Ti-1) ) # optimized anyway so this value will be ignored
  
  init=list(m=unlist(init$m),s=unlist(init$s))
  # fit alternative model
  
  full_gradient_function=function(g) grad_log_prob(sampler_full,g,T)
  full_likelihood=function(g) log_prob(sampler_full,g,T,F)
  
  set.seed(1)
  v_full=svem(full_gradient_function, to_optim, init=init, plot.elbo = F, samples_for_elbo = 0, log_prob = full_likelihood, ... )

  # using the same random draws to estimate the deviance is more statistically efficient
  nintegrate=sum(!to_optim)
  loglr=mean( foreach(g=1:elbo_samples, .combine=c) %do% {
    x=rnorm(nintegrate)
    v_full$elbo_func( x ) - v_null$elbo_func( x ) 
  }, na.rm=T )
  #loglr=mean(v_full$elbo_progress[(maxit/2):maxit] - v_null$elbo_progress[(maxit/2):maxit])
    
  list(loglr=loglr, df=Ti-1, lrtp=pchisq( 2.0*loglr, lower.tail = F , df=Ti-1 ), fit_full=rstan:::rstan_relist(v_full$m, sk_full) )
}
