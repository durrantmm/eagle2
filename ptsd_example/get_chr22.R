counts_data=fread("zcat < ptsd_counts.txt.gz", sep="\t")
setDF(counts_data)

# Load gene - SNP mapping
gene_snp=fread("zcat < gene_snps.txt.gz", sep="\t")
setDF(gene_snp)

counts_data=counts_data[counts_data$chr=="chr22",]
gene_snp=gene_snp[gene_snp$chr=="chr22",]

write.table(counts_data, file="ptsd_chr22.txt", quote = F, row.names = F, sep="\t")
write.table(gene_snp, file="gene_snps_chr22.txt", quote = F, row.names = F, sep="\t")
