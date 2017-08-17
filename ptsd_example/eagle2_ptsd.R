# Data adapted from 
# Gene networks specific for innate immunity define post-traumatic stress disorder. M S Breen, A X Maihofer, S J Glatt, D S Tylee, S D Chandler, M T Tsuang, V B Risbrough, D G Baker, D T O'Connor, C M Nievergelt and C H Woelk. Molecular psychiatry 2015. 
# Accession number GSE64814. 

# Load required packages
require(data.table)
require(dplyr)
require(doMC)
registerDoMC(detectCores()-1)

USE_RANDOM_EFFECT=T

require(reshape2)

require(eagle2)
require(abind)

# Load metadata
meta=read.table("SraRunTable.txt", sep="\t", header=T, stringsAsFactors = F)
rownames(meta)=meta$Run_s
meta$id=sort(rep(1:94,2))

# Load counts 
counts_data=fread("zcat < ptsd_chr22.txt.gz", sep="\t", data.table = F)

# Load gene - SNP mapping
gene_snp=fread("zcat < gene_snps_chr22.txt.gz", sep="\t", data.table = F)

test_genes=unique(gene_snp$gene)

counts_data$cond=factor(meta[counts_data$srr,"time_point_s"])
counts_data$id=factor(meta[counts_data$srr,"id"])

gene_snp$chr_pos=paste( gene_snp$chr, gene_snp$pos, sep=":")
# chr22:-:300282:RP3-394A18.1
# Helper function to make [individuals x conditions x SNPs] count tensor
cast_me=function(snp_names, stat) { do.call(abind, c(foreach(snp=snp_names) %do% {
  res=dcast(counts_data[ counts_data$chr_pos==snp, ,drop=F], id ~ cond, value.var = stat, drop=F )
  res$id=NULL
  res[is.na(res)]=0
  res
}, along=3 )) }

# Run over genes (may take an hour or so)
results = foreach(gene=test_genes, .combine = bind_rows) %dopar% {
  
  # Get SNPs in this gene
  snp_names=gene_snp$chr_pos[ gene_snp$gene %in% gene ]
  #snp_indices = which(counts_data$chr_pos %in% snp_names)
  snp_names = intersect(counts_data$chr_pos, snp_names)

  # Get alternative allele counts
  a=cast_me( snp_names, "y" )
  # Get total counts
  nh=a+cast_me( snp_names, "r" )
  
  # Remove SNP-gene pairs that are probably homozygous
  filtered_data=detect_homozygotes(a,nh)
  if (is.null(filtered_data)) return(NULL)
  
  if (F)  {
    require(rstan)
    ys=filtered_data$a
    ns=filtered_data$nh
    concShape=1.0001
    concRate=1e-4
    elbo_samples=0
    stanmodels=eagle2:::stanmodels
    svem=eagle2:::svem
    get_skeleton=eagle2:::get_skeleton
  }
  
  # Run EAGLE2
  eagle_results = tryCatch( {
  	if (USE_RANDOM_EFFECT) eagle2_re( filtered_data$a, filtered_data$nh, iterations=1000, elbo_samples=3000, trace=2 ) else eagle2( filtered_data$a, filtered_data$nh ) 
	}, error=function(e) NULL )
  if (is.null(eagle_results)) return(NULL)
  # Store LRT p-value, coefficients, standard errors and Wald p-values [assumes just two conditions]
  if (USE_RANDOM_EFFECT) data.frame(gene=gene, loglr=eagle_results$loglr, lrtp=eagle_results$lrtp, coef=eagle_results$fit_full$beta[2]) else data.frame(gene=gene, loglr=eagle_results$loglr, lrtp=eagle_results$lrtp, get_coefs(eagle_results$fit_full)[2,])
} 

write.table(results, "ptsd_results.txt", quote=F, sep="\t", row.names=F)

# require(ggplot2)

# Do the LRT and Wald p-values agree
#theme_set(theme_bw(base_size = 14))
#ggplot(results, aes(lrtp, wald_p)) + geom_point()

# P-values are somewhat conservative. 
#hist(results$lrtp)

# Plot the top three associated genes. 
#pdf("plots.pdf", height=12, width=15)
for (i in 1:3) {
  print(i)
  gene=results$gene[ order(results$lrtp)[i] ]
  gs=gene_snp[ gene_snp$gene %in% gene, ]
  gs=paste( gs$chr, gs$pos, sep=":")
  
  snp_data=counts_data[counts_data$chr_pos %in% gs, ]
  colnames(snp_data)[colnames(snp_data)=="id"]="Individual"
  colnames(snp_data)[colnames(snp_data)=="ar"]="AllelicRatio"
  colnames(snp_data)[colnames(snp_data)=="cond"]="TimePoint"
  colnames(snp_data)[colnames(snp_data)=="chr_pos"]="SNP"
  snp_data$Coverage=snp_data$r+snp_data$y
  snp_data$TimePoint=as.numeric(as.factor(snp_data$TimePoint))
  snp_data$Individual=as.factor(snp_data$Individual)
  #print( ggplot(snp_data, aes( TimePoint, AllelicRatio, col=SNP, label=Individual, shape=Individual, linetype=Individual)) + geom_text( aes( size=Coverage)) + ylim(0,1) + theme_bw(base_size = 16) + geom_line() + scale_shape_manual(values=seq_along(levels(snp_data$Individual))) + ggtitle(gene) ) 
  snp_data$inter=interaction( snp_data$TimePoint, snp_data$SNP )
  snp_data$ii=as.numeric(snp_data$inter)
  print( ggplot(snp_data, aes( inter, AllelicRatio, col=SNP, label=Individual, shape=Individual)) + geom_text( aes( size=Coverage, label=Individual)) + ylim(0,1) + theme_bw(base_size = 16)  + scale_shape_manual(values=seq_along(levels(snp_data$Individual))) + ggtitle(gene) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_line(aes(ii,AllelicRatio),alpha=.5) + xlab("Condition.SNP") )
}
#dev.off()