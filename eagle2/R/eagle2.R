#' Robust inverse of PSD matrix using SVD
#' 
#' @param si Matrix to invert
#' @param eigenThreshold Threshold eigenvalues smaller than this
#' @return Inverted matrix
robustSolve=function(si, eigenThreshold=0.01){
  svdSigma=eigen(si)
  svdSigma$values=pmax(svdSigma$values,eigenThreshold)
  svdSigma$vectors %*% diag(1/svdSigma$values) %*% t(svdSigma$vectors)
}

#' Function to extract coefficients, standard errors, and Wald p-values. 
#' 
#' @param fit Fitted object returned from rstan::optimizing. Must have $par$beta and $hessian available
#' @return Data.frame with coefs, standard errors, Wald p-values
#' @export
get_coefs=function(fit) {
  variance=robustSolve(-fit$hessian)
  dimnames(variance)=dimnames(fit$hessian)
  P=length(fit$par$beta)
  betase=sqrt(diag(variance))[paste0("beta.",seq_len(P))]
  beta=fit$par$beta
  zscore=fit$par$beta / betase
  data.frame(coef=fit$par$beta,se=betase,wald_p=2.0*pnorm(-abs(zscore)))
}

#' Detect potential homozygote (or imprinted) SNP-individual pairs and zero them out. 
#' 
#' @param a Numerator counts [n x T x K] where n are individuals, T are timepoints, K are SNPs
#' @param nh denominator counts [n x T x K]
#' @param concShape Shape of prior on concentration
#' @param concRate Rate of prior on concentration
#' @param errorRate Assumed probability of erroneously observing a read with the alternative allele, due to sequencing or alignment error
#' @param posterior_threshold If posterior probability of being a het is below this the individual/SNP pair will be ignored. 
#' @return Filtered a and nh matrices, and which SNPs were kept. 
#' @import foreach
#' @importFrom rstan optimizing
#' @export
detect_homozygotes=function(a,nh,concShape=1.0001,concRate=1e-4,errorRate=0.01,posterior_theshold=0.95,verbose=T)
{
  
  homo=foreach(snp_index=seq_len(dim(nh)[3]), .combine=cbind) %do% {
    as=a[,,snp_index]
    nhh=nh[,,snp_index]
    ind_to_keep=rowSums(nhh)>0
    treat_to_keep=colSums(nhh)>0
    
    as=as[ind_to_keep,treat_to_keep,drop=F]
    nhh=nhh[ind_to_keep,treat_to_keep,drop=F]
    
    o=optimizing(stanmodels$is_het, dat=list(N=nrow(nhh), T=ncol(nhh), errorRate=errorRate, concShape=concShape, concRate=concRate, ys=as, ns=nhh), as_vector=F)
    eo=exp(o$par$probs)
    pr=sweep(eo, 1, rowSums(eo), "/")
    
    homo=pr[,1]<posterior_theshold
    
    if (verbose) cat("Removing",sum(homo),"individual(s) from SNP",dimnames(a)[[3]][snp_index],"\n")
    
    nh[which(ind_to_keep)[homo],treat_to_keep,snp_index]=0
    a[which(ind_to_keep)[homo],treat_to_keep,snp_index]=0
  }
  
  snp_to_keep=apply(nh>0, 3, any)
  ind_to_keep=apply(nh,1,sum) > 0
  if (sum(snp_to_keep)==0 | sum(ind_to_keep)==0) return(NULL)
  a=a[ind_to_keep,,snp_to_keep,drop=F]
  nh=nh[ind_to_keep,,snp_to_keep,drop=F]
  
  list(a=a,nh=nh,snp_to_keep=snp_to_keep)
}

#' Beta binomial GLM with flips. Prior on concentration parameter is Gamma(concShape,concRate)
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
eagle2=function(ys,ns,concShape=1.0001,concRate=1e-4,hessian=T,...) {
  
  N=dim(ys)[1]
  Ti=dim(ys)[2]
  K=dim(ys)[3]
  
  preprocess_x=function(g) aperm( abind(g,along=3), c(3,1,2) )
  
  xNull=preprocess_x( foreach(i=seq_len(Ti)) %do% matrix(1,N,1) )
  dat=list(N=N,P=1,T=Ti,K=K,ys=ys,ns=ns,x=xNull,concShape=concShape,concRate=concRate)

  # Fit null model
  fit_null <- optimizing(stanmodels$bb, data=dat, as_vector=F, hessian=T, ...)
  
  # Initialize alternative model using null model
  initFull=fit_null$par
  initFull$beta=c(fit_null$par$beta,numeric(Ti-1))
  
  # Fit alternative model
  datFull=dat
  temp=matrix(0,N,Ti)
  temp[,1]=1
  datFull$x=preprocess_x(foreach(i=seq_len(Ti)) %do% { temp[,i]=1; temp })
  datFull$P=Ti
  fit_full <- optimizing(stanmodels$bb, data=datFull, init=initFull, as_vector=F, hessian=T, ...)
  
  loglr=fit_full$value - fit_null$value
  
  list( loglr=loglr, lrtp=pchisq( 2.0*loglr, lower.tail = F , df=Ti-1 ), fit_full=fit_full, fit_null=fit_null )
}


