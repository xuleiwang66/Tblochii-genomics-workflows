#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("Usage: cox_candidate_effects.R candidates.tsv model_data.tsv kinship.txt null.rds bfile_prefix workdir output.tsv")
}
suppressPackageStartupMessages({library(data.table); library(coxmeg)})

candidate <- fread(args[1])
d <- fread(args[2])
K <- as.matrix(fread(args[3], header = FALSE))
null <- readRDS(args[4])
bfile <- args[5]
work <- args[6]
out <- args[7]
dir.create(work, recursive = TRUE, showWarnings = FALSE)

if (nrow(K) != nrow(d) || ncol(K) != nrow(d)) stop("Kinship dimension does not match Cox model data")
rownames(K) <- d$sample_id
colnames(K) <- d$sample_id

if (nrow(candidate) == 0L) {
  empty <- data.table(
    snp_id = character(), beta_log_hr = numeric(), se = numeric(), p_wald = numeric(),
    tau = numeric(), HR = numeric(), HR95_CI_lower = numeric(), HR95_CI_upper = numeric(),
    tolerance_effect = numeric()
  )
  fwrite(empty, out, sep = "\t", quote = FALSE, na = "NA")
  quit(save = "no", status = 0)
}

extract <- file.path(work, "snps.txt")
allele_file <- file.path(work, "alleles.txt")
fwrite(candidate[, .(snp_id)], extract, col.names = FALSE)
fwrite(candidate[, .(snp_id, effect_allele)], allele_file, col.names = FALSE, sep = "\t")

prefix <- file.path(work, "candidate_genotypes")
status <- system2(
  "plink2",
  c(
    "--bfile", bfile, "--chr-set", "24", "no-xy", "no-mt", "--extract", extract,
    "--export", "Av", "--export-allele", allele_file, "--out", prefix, "--threads", "1"
  )
)
if (status != 0) stop("PLINK dosage export failed")

traw <- fread(paste0(prefix, ".traw"))
fam <- fread(paste0(bfile, ".fam"), header = FALSE)
metadata_columns <- c("CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT")
sample_columns <- setdiff(names(traw), metadata_columns)
expected_fid_iid <- paste(fam$V1, fam$V2, sep = "_")
if (!identical(sample_columns, expected_fid_iid) && !identical(sample_columns, as.character(fam$V2))) {
  stop("PLINK Av sample columns do not match FAM order")
}

variant_order <- match(candidate$snp_id, traw$SNP)
if (anyNA(variant_order)) stop("One or more candidate SNPs are absent from PLINK Av output")
traw <- traw[variant_order]

dosage_variant_by_sample <- as.matrix(traw[, ..sample_columns])
storage.mode(dosage_variant_by_sample) <- "double"
genotype <- t(dosage_variant_by_sample)
rownames(genotype) <- fam$V2
colnames(genotype) <- candidate$snp_id
if (!identical(rownames(genotype), d$sample_id)) stop("Candidate genotype sample order does not match Cox model data")
if (any(!is.finite(genotype))) stop("Candidate genotype matrix contains missing or non-finite dosages")

outcome <- cbind(d$time_to_event_h, d$event)
covariates <- as.matrix(d[, .(body_weight_g, source_L, source_Q, PC1, PC2, PC3)])
storage.mode(covariates) <- "double"

fit <- coxmeg::coxmeg_m(
  X = genotype, outcome = outcome, corr = K, type = "dense", cov = covariates,
  tau = null$tau, min_tau = 1e-4, max_tau = 5, eps = 1e-6, order = 1,
  detap = "exact", opt = "bobyqa", score = FALSE, threshold = 1,
  solver = 1, spd = FALSE, verbose = TRUE
)
summary <- as.data.table(fit$summary)
required <- c("tau", "beta_exact", "sd_beta_exact", "p_exact")
missing <- setdiff(required, names(summary))
if (length(missing)) stop("Unexpected coxmeg_m schema; missing: ", paste(missing, collapse = ","))
if (nrow(summary) != nrow(candidate)) stop("coxmeg_m result count does not match candidate count")

result <- data.table(
  snp_id = candidate$snp_id,
  beta_log_hr = as.numeric(summary$beta_exact),
  se = as.numeric(summary$sd_beta_exact),
  p_wald = as.numeric(summary$p_exact),
  tau = as.numeric(summary$tau)
)
result[, `:=`(
  HR = exp(beta_log_hr),
  HR95_CI_lower = exp(beta_log_hr - 1.96 * se),
  HR95_CI_upper = exp(beta_log_hr + 1.96 * se),
  tolerance_effect = -beta_log_hr
)]
fwrite(result, out, sep = "\t", quote = FALSE, na = "NA")
