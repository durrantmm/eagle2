

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
#' @importFrom rstan optimizing sampling grad_log_prob log_prob
#' @importFrom foreach foreach %do%
#' @importFrom abind abind
#' @export
eagle2_re=function(ys,ns,concShape=1.0001,concRate=1e-4,USE_LBFGS=T,burnin=3000,iterations=1000,elbo_samples=1000,learning_rate=1.,seed=1,...) {
  
  N=dim(ys)[1] # individuals
  Ti=dim(ys)[2] # conditions
  K=dim(ys)[3] # SNPs
  
  # gradient_scale_factor=2000 / (N*K) # sum(ns) / Ti # average number of reads per condition
  gradient_scale_factor=1.

  mode(ns)="integer"
  mode(ys)="integer"
  
  xNull=aperm( abind( foreach(i=seq_len(Ti)) %do% matrix(1,N,1), along=3 ), c(3,1,2) )
  
  dat=list(N=N,P=1L,T=Ti,K=K,ys=ys,ns=ns,x=xNull,concShape=concShape,concRate=concRate)
  
  # get stan_fit object so we can call sampler$grad_log_prob for the gradient
  # no sampling is done, just a hack to get a stanfit object
  sampler=sampling(stanmodels$bb_with_re, dat, iter=1, chains=0)
  
  # specific which parameters to optimize (rather than integrate over)
  make_skeleton=function(samp) {
    sk=get_skeleton(samp)
    sk$conc[]=T
    sk$p[]=T
    sk$log_rev=T
    sk$beta[]=T
    sk$re[]=F
    sk
  }
  
  sk=make_skeleton(sampler)
  
  to_optim=as.logical(unlist(sk))
  
  # initialize using the model fit without the random effects
  o=optimizing(stanmodels$bb, dat, as_vector=F, seed=seed)
  init=list(m=get_skeleton(sampler), s=get_skeleton(sampler))
  init$m$beta=o$par$beta
  #init$m$p=logit(o$par$p)
  #init$m$p=pmin( pmax(init$m$p,-4), 4 )
  init$m$p=logit(o$par$p) 
  init$m$conc=log(o$par$conc)
  #init$m$conc=o$par$conc
  for (n in names(init$s)) init$s[[n]][]=0.5 # set all standard deviations
  for (n in names(init)) init[[n]]=unlist(init[[n]])
  
  # set the seed so we can use the same seed for the alternative model
  set.seed(seed)
  
  null_gradient_function=function(g) { grad_log_prob(sampler,g,T) / gradient_scale_factor }
  #null_gradient_function(init$m)
  null_likelihood=function(g) { log_prob(sampler,g,T,F) }
  
  cat("Fitting initial null model\n")
  v_init=if (USE_LBFGS) {
    v_init=svem_lbfgs(null_likelihood, null_gradient_function, to_optim, init, samples=10, iterations=iterations, ...)
    set.seed(seed)
    svem_lbfgs(null_likelihood, null_gradient_function, to_optim, v_init, iterations=iterations, ...)
  } else svem(null_gradient_function, to_optim, init, plot.elbo = F, log_prob = null_likelihood, iterations=burnin, master_stepsize=learning_rate, ...)
  if (is.null(init)) return(NULL)
  cat("Re-fitting null model\n")
  set.seed(seed)
  v_null=if (USE_LBFGS) 
    svem_lbfgs(null_likelihood, null_gradient_function, to_optim, v_init, iterations=iterations, ...) else svem(null_gradient_function, to_optim, v_init, plot.elbo = F, log_prob = null_likelihood, iterations=iterations, master_stepsize=learning_rate, ...)
  if (is.null(v_null)) return(NULL)
  # Fit alternative model
  dat_full=dat
  temp=matrix(0,N,Ti)
  temp[,1]=1
  dat_full$x=aperm( abind( foreach(i=seq_len(Ti)) %do% { temp[,i]=1; temp }, along=3 ), c(3,1,2) )
  
  dat_full$P=Ti
  
  # specify which parameters to optimize
  sampler_full=sampling(stanmodels$bb_with_re, dat_full, iter=1, chains=0)
  
  sk_full = make_skeleton(sampler_full)
  to_optim_full=as.logical(unlist(sk_full))
  
  # initialization choices: initializing from null fit seems best
  init=list( m=rstan:::rstan_relist(v_init$m, sk), s=rstan:::rstan_relist(v_init$s, sk) )
  
  init$m$beta=c(init$m$beta, numeric(Ti-1) )
  init$s$beta=c(init$s$beta, numeric(Ti-1) ) # optimized anyway so this value will be ignored
  
  init=list(m=unlist(init$m),s=unlist(init$s))
  # fit alternative model
  
  full_gradient_function=function(g) { grad_log_prob(sampler_full,g,T)/gradient_scale_factor }
  full_likelihood=function(g) log_prob(sampler_full,g,T,F)
  
  set.seed(seed)
  cat("Fitting full model\n")
  v_full=if (USE_LBFGS) 
     svem_lbfgs(full_likelihood, full_gradient_function, to_optim_full, init=init, iterations=iterations, ...) else svem(full_gradient_function, to_optim_full, init=init, plot.elbo = F, log_prob = full_likelihood, iterations=iterations, master_stepsize=learning_rate, ... )
  if (is.null(v_full)) return(NULL)
  # using the same random draws to estimate the deviance is more statistically efficient
  nintegrate=sum(!to_optim)
  loglr=if (USE_LBFGS) (v_full$elbo_opt - v_null$elbo_opt) else { 
    mean( foreach(g=1:elbo_samples, .combine=c) %do% {
    x=rnorm(nintegrate)
    v_full$elbo_func( x ) - v_null$elbo_func( x ) 
    }, na.rm=T ) }
  
  list(loglr=loglr, df=Ti-1, lrtp=pchisq( 2.0*loglr, lower.tail = F , df=Ti-1 ), fit_full=rstan:::rstan_relist(v_full$m, sk_full) )
}
