#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) stop('Usage: fit_cox_null.R model_data.tsv kinship.txt model.rds tau.tsv')
suppressPackageStartupMessages({library(data.table); library(coxmeg)})
d <- fread(args[1]); K <- as.matrix(fread(args[2],header=FALSE)); storage.mode(K)<-'double'
X <- as.matrix(d[,.(body_weight_g,source_L,source_Q,PC1,PC2,PC3)]); storage.mode(X)<-'double'
outcome <- cbind(as.numeric(d$time_to_event_h),as.integer(d$event))
fit <- coxmeg::coxmeg(outcome=outcome,corr=K,type='dense',X=X,eps=1e-6,min_tau=1e-4,max_tau=5,order=1,detap='exact',opt='bobyqa',solver=1,spd=FALSE,verbose=TRUE)
saveRDS(fit,args[3],compress=FALSE); fwrite(data.table(tau=fit$tau),args[4],sep='\t')
