#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: merge_cox_results.R output.tsv input1.tsv ...")
suppressPackageStartupMessages(library(data.table))

x <- rbindlist(lapply(args[-1], fread), use.names = TRUE, fill = TRUE)
chr_levels <- sprintf("Chr%02d", 1:24)
x[, chromosome_order := match(original_CHR, chr_levels)]
if (anyNA(x$chromosome_order)) {
  stop("Unexpected chromosome labels in Cox results: ", paste(unique(x[is.na(chromosome_order), original_CHR]), collapse = ","))
}
setorder(x, chromosome_order, original_POS, mapped_SNP)
x[, chromosome_order := NULL]
fwrite(x, args[1], sep = "\t", quote = FALSE, na = "NA")
