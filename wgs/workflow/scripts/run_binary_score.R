#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop('Usage: run_binary_score.R null.rds bed_prefix output.tsv')
suppressPackageStartupMessages(library(GMMAT))
obj <- readRDS(args[1])
GMMAT::glmm.score(obj=obj, infile=args[2], outfile=args[3], center=TRUE,
  select=NULL, MAF.range=c(1e-7,0.5), miss.cutoff=1,
  missing.method='impute2mean', nperbatch=100L, tol=1e-5,
  ncores=1L, verbose=TRUE)
