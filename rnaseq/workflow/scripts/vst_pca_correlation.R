suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
})

robust_z <- function(x) {
  med <- median(x, na.rm=TRUE)
  sc <- mad(x, constant=1.4826, na.rm=TRUE)
  if (!is.finite(sc) || sc == 0) sc <- sd(x, na.rm=TRUE)
  if (!is.finite(sc) || sc == 0) return(rep(0, length(x)))
  (x-med)/sc
}

counts_path <- snakemake@input[["counts"]]
samples_path <- snakemake@input[["samples"]]
tissue_name <- snakemake@wildcards[["tissue"]]
min_count <- as.integer(snakemake@params[["min_count"]])
min_samples <- as.integer(snakemake@params[["min_samples"]])
source_levels <- unlist(snakemake@params[["source_levels"]])
treatment_levels <- unlist(snakemake@params[["treatment_levels"]])
th <- snakemake@params[["thresholds"]]

meta <- fread(samples_path, sep="\t")
meta <- meta[tissue == tissue_name]
if (nrow(meta) == 0) stop("No samples found for tissue: ", tissue_name)
meta[, source_cohort := factor(source_cohort, levels=source_levels)]
meta[, treatment := factor(treatment, levels=treatment_levels)]
setorder(meta, source_cohort, treatment, biological_replicate, sample_id)

count_dt <- fread(counts_path, sep="\t")
sample_ids <- meta$sample_id
count_mat <- as.matrix(count_dt[, ..sample_ids])
rownames(count_mat) <- count_dt$gene_id
storage.mode(count_mat) <- "numeric"
if (anyNA(count_mat) || any(count_mat < 0) || any(abs(count_mat-round(count_mat)) > 1e-8)) stop("Invalid raw counts")
count_mat <- round(count_mat); storage.mode(count_mat) <- "integer"
keep <- rowSums(count_mat >= min_count) >= min_samples
filtered <- count_mat[keep, , drop=FALSE]
if (nrow(filtered) == 0) stop("No genes retained for VST")

coldata <- as.data.frame(meta); rownames(coldata) <- coldata$sample_id
dds <- DESeqDataSetFromMatrix(filtered, coldata, design=~ source_cohort + treatment)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind=FALSE)
vst_mat <- assay(vsd)
if (anyNA(vst_mat) || any(!is.finite(vst_mat))) stop("Invalid VST matrix")

vst_dt <- data.table(gene_id=rownames(vst_mat))
vst_dt <- cbind(vst_dt, as.data.table(vst_mat))
fwrite(vst_dt, snakemake@output[["vst"]], sep="\t", quote=FALSE)

pca <- prcomp(t(vst_mat), center=TRUE, scale.=FALSE)
percent <- pca$sdev^2/sum(pca$sdev^2)*100
coord <- data.table(
  sample_id=rownames(pca$x), tissue=tissue_name,
  source_cohort=as.character(meta[match(rownames(pca$x), sample_id), source_cohort]),
  treatment=as.character(meta[match(rownames(pca$x), sample_id), treatment]),
  biological_replicate=as.character(meta[match(rownames(pca$x), sample_id), biological_replicate]),
  replacement_pituitary=as.character(meta[match(rownames(pca$x), sample_id), replacement_pituitary]),
  PC1=pca$x[,1], PC2=pca$x[,2], PC3=if(ncol(pca$x)>=3) pca$x[,3] else NA_real_,
  PC1_percent_variance=percent[1], PC2_percent_variance=percent[2], PC3_percent_variance=if(length(percent)>=3) percent[3] else NA_real_
)
fwrite(coord, snakemake@output[["pca"]], sep="\t", quote=FALSE, na="NA")

cor_mat <- cor(vst_mat, method="pearson")
cor_long <- as.data.table(as.table(cor_mat))
setnames(cor_long, c("sample_id_1", "sample_id_2", "pearson_correlation"))
cor_long[, tissue := tissue_name]
setcolorder(cor_long, c("tissue", "sample_id_1", "sample_id_2", "pearson_correlation"))
fwrite(cor_long, snakemake@output[["correlations"]], sep="\t", quote=FALSE)

cor_no_self <- cor_mat; diag(cor_no_self) <- NA_real_
z1 <- robust_z(coord$PC1); z2 <- robust_z(coord$PC2)
robust_distance <- sqrt(z1^2+z2^2)
warning_rows <- list()
for (i in seq_len(nrow(meta))) {
  sid <- meta$sample_id[i]
  peers <- meta[source_cohort == meta$source_cohort[i] & treatment == meta$treatment[i] & sample_id != sid, sample_id]
  within <- cor_mat[sid, peers]
  median_within <- median(within, na.rm=TRUE)
  level <- "OK"
  if (median_within < as.numeric(th[["review_within_group_correlation"]]) || robust_distance[i] > as.numeric(th[["review_robust_pc_distance"]])) level <- "REVIEW"
  if (median_within < as.numeric(th[["flag_within_group_correlation"]]) || robust_distance[i] > as.numeric(th[["flag_robust_pc_distance"]])) level <- "FLAG"
  if (level != "OK") {
    warning_rows[[length(warning_rows)+1]] <- data.table(
      sample_id=sid, tissue=tissue_name,
      source_cohort=as.character(meta$source_cohort[i]), treatment=as.character(meta$treatment[i]),
      biological_replicate=as.character(meta$biological_replicate[i]), replacement_pituitary=as.character(meta$replacement_pituitary[i]),
      median_correlation_within_group=median_within,
      robust_PC1_PC2_distance=robust_distance[i], warning_level=level,
      action="review_only_no_automatic_exclusion"
    )
  }
}
warning_dt <- if(length(warning_rows)) rbindlist(warning_rows) else data.table(
  sample_id=character(), tissue=character(), source_cohort=character(), treatment=character(), biological_replicate=character(), replacement_pituitary=character(),
  median_correlation_within_group=numeric(), robust_PC1_PC2_distance=numeric(), warning_level=character(), action=character()
)
fwrite(warning_dt, snakemake@output[["warnings"]], sep="\t", quote=FALSE)
