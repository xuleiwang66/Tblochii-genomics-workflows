#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) stop('Usage: binary_candidate_effects.R model_data.tsv kinship.txt bed_prefix candidates.tsv output.tsv')
suppressPackageStartupMessages({library(data.table); library(GMMAT)})
d <- fread(args[1]); K <- as.matrix(fread(args[2],header=FALSE)); rownames(K)<-d$sample_id; colnames(K)<-d$sample_id; cand <- fread(args[4])
snps <- unique(cand$snp_id)
res <- GMMAT::glmm.wald(
  fixed=binary_status ~ body_weight_g + source_L + source_Q + PC1 + PC2 + PC3,
  data=as.data.frame(d), kins=K, id='sample_id', family=binomial(link='logit'),
  infile=args[3], snps=snps, method='REML', method.optim='AI', maxiter=500L,
  tol=1e-5, center=TRUE, select=NULL, missing.method='impute2mean', verbose=TRUE)
x <- as.data.table(res)
# Preserve package output and add portable effect columns when detectable.
beta_col <- intersect(c('BETA','beta','Estimate'),names(x))[1]
se_col <- intersect(c('SE','se','Std. Error'),names(x))[1]
p_col <- intersect(c('PVAL','P','p','pvalue'),names(x))[1]
snp_col <- intersect(c('SNP','snp','ID'),names(x))[1]
if (any(is.na(c(beta_col,se_col,p_col,snp_col)))) stop('Unexpected glmm.wald output columns: ',paste(names(x),collapse=','))
out <- data.table(snp_id=as.character(x[[snp_col]]), beta_log_or=as.numeric(x[[beta_col]]), se=as.numeric(x[[se_col]]), p_wald=as.numeric(x[[p_col]]))
out[, `:=`(OR=exp(beta_log_or), OR95_CI_lower=exp(beta_log_or-1.96*se), OR95_CI_upper=exp(beta_log_or+1.96*se), tolerance_effect=beta_log_or)]
fwrite(out,args[5],sep='\t',quote=FALSE,na='NA')
