#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop('Usage: fit_binary_null.R model_data.tsv kinship.txt output.rds')
suppressPackageStartupMessages({library(data.table); library(GMMAT)})
d <- fread(args[1]); K <- as.matrix(fread(args[2], header=FALSE)); storage.mode(K) <- 'double'
stopifnot(nrow(d)==nrow(K), ncol(K)==nrow(K), identical(as.character(packageVersion('GMMAT')), '1.4.1'))
rownames(K) <- d$sample_id; colnames(K) <- d$sample_id
fit <- GMMAT::glmmkin(
  fixed=binary_status ~ body_weight_g + source_L + source_Q + PC1 + PC2 + PC3,
  data=as.data.frame(d), kins=K, id='sample_id', family=binomial(link='logit'),
  method='REML', method.optim='AI', maxiter=500L, tol=1e-5,
  taumin=1e-5, taumax=1e5, tauregion=10L, verbose=TRUE)
saveRDS(fit,args[3],compress=FALSE)
