#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop("Usage: run_cox_score_by_chr.R chromosome bed_prefix model_data.tsv kinship.txt null.rds workdir out.tsv plink2")
}
suppressPackageStartupMessages({library(data.table); library(coxmeg)})

chr <- args[1]
source <- args[2]
d <- fread(args[3])
K <- as.matrix(fread(args[4], header = FALSE))
null <- readRDS(args[5])
work <- args[6]
out <- args[7]
plink2 <- args[8]
dir.create(work, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

if (nrow(K) != nrow(d) || ncol(K) != nrow(d)) stop("Kinship dimension does not match Cox model data")
rownames(K) <- d$sample_id
colnames(K) <- d$sample_id

canonical_chromosome <- function(value) {
  text <- sub("^chr", "", sub("^Chr", "", as.character(value)))
  number <- suppressWarnings(as.integer(text))
  ifelse(is.na(number), as.character(value), sprintf("Chr%02d", number))
}

prefix <- file.path(work, "genotypes")
full_bim <- fread(paste0(source, ".bim"), header = FALSE)
keep_ids <- full_bim[canonical_chromosome(V1) == chr, V2]
if (length(keep_ids) == 0L) stop("No variants found for ", chr)
extract <- file.path(work, "extract.snps.txt")
fwrite(data.table(keep_ids), extract, col.names = FALSE)

status <- system2(
  plink2,
  c(
    "--bfile", source, "--chr-set", "24", "no-xy", "no-mt", "--allow-extra-chr",
    "--extract", extract, "--make-bed", "--out", prefix, "--threads", "1"
  )
)
if (status != 0) stop("PLINK chromosome extraction failed")

original_bim <- fread(paste0(prefix, ".bim"), header = FALSE)
recode_bim <- copy(original_bim)
recode_bim[[1]] <- 1L
fwrite(recode_bim, paste0(prefix, ".bim"), sep = "\t", col.names = FALSE, quote = FALSE)

pheno <- data.table(FID = d$sample_id, IID = d$sample_id, time = d$time_to_event_h, status = d$event)
covariates <- data.table(
  FID = d$sample_id, IID = d$sample_id, body_weight_g = d$body_weight_g,
  source_L = d$source_L, source_Q = d$source_Q, PC1 = d$PC1, PC2 = d$PC2, PC3 = d$PC3
)
pheno_file <- file.path(work, "pheno.tsv")
covariate_file <- file.path(work, "covariates.tsv")
fwrite(pheno, pheno_file, sep = "\t")
fwrite(covariates, covariate_file, sep = "\t")

result <- coxmeg::coxmeg_plink(
  pheno = pheno_file, corr = K, type = "dense", bed = prefix, tmp_dir = work,
  cov_file = covariate_file, tau = null$tau, eps = 1e-6, min_tau = 1e-4,
  max_tau = 5, order = 1, detap = "exact", solver = 1, spd = FALSE,
  maf = 1e-7, score = TRUE, verbose = TRUE
)
x <- as.data.table(result$summary)
required <- c("snp.id", "chromosome", "position", "allele", "afreq", "afreq_inc", "index", "score", "score_test", "p")
missing <- setdiff(required, names(x))
if (length(missing)) stop("Unexpected coxmeg_plink schema; missing: ", paste(missing, collapse = ","))

map <- match(as.character(x$snp.id), as.character(original_bim[[2]]))
if (anyNA(map)) stop("Cannot map coxmeg rows to source BIM by SNP ID")
alleles <- tstrsplit(as.character(x$allele), "[/|]", perl = TRUE)
if (length(alleles) < 2) stop("Cannot parse coxmeg allele field")

x[, `:=`(
  original_CHR = chr,
  original_POS = original_bim[[4]][map],
  mapped_SNP = original_bim[[2]][map],
  counted_allele = alleles[[1]],
  other_allele = alleles[[2]],
  counted_allele_frequency = afreq_inc,
  cox_sensitivity_score = score,
  cox_tolerance_score = -score,
  tested_N = result$nsam,
  tau_reused = null$tau
)]
fwrite(x, out, sep = "\t", quote = FALSE, na = "NA")
