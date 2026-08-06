suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
})

counts_path <- snakemake@input[["counts"]]
samples_path <- snakemake@input[["samples"]]
tissue_name <- snakemake@wildcards[["tissue"]]
min_count <- as.integer(snakemake@params[["min_count"]])
min_samples <- as.integer(snakemake@params[["min_samples"]])
padj_cutoff <- as.numeric(snakemake@params[["padj"]])
lfc_cutoff <- as.numeric(snakemake@params[["lfc"]])
independent_filtering <- as.logical(snakemake@params[["independent"]])
source_levels <- unlist(snakemake@params[["source_levels"]])
treatment_levels <- unlist(snakemake@params[["treatment_levels"]])
contrast_names <- unlist(snakemake@params[["contrast_names"]])
contrast_numerators <- unlist(snakemake@params[["contrast_numerators"]])
contrast_denominators <- unlist(snakemake@params[["contrast_denominators"]])

meta <- fread(samples_path, sep="\t")
meta <- meta[tissue == tissue_name]
if (nrow(meta) == 0) stop("No samples found for tissue: ", tissue_name)
meta[, source_cohort := factor(source_cohort, levels=source_levels)]
meta[, treatment := factor(treatment, levels=treatment_levels)]
if (anyNA(meta$source_cohort) || anyNA(meta$treatment)) stop("Invalid design level in metadata")

count_dt <- fread(counts_path, sep="\t")
if (anyDuplicated(count_dt$gene_id)) stop("Duplicated gene_id in raw count matrix")
sample_ids <- meta$sample_id
if (!all(sample_ids %in% names(count_dt))) stop("Count matrix missing tissue samples")
count_mat <- as.matrix(count_dt[, ..sample_ids])
rownames(count_mat) <- count_dt$gene_id
storage.mode(count_mat) <- "numeric"
if (anyNA(count_mat) || any(count_mat < 0) || any(abs(count_mat-round(count_mat)) > 1e-8)) stop("Raw counts must be non-negative integers")
count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"

keep <- rowSums(count_mat >= min_count) >= min_samples
filtered <- count_mat[keep, , drop=FALSE]
if (nrow(filtered) == 0) stop("No genes retained after low-expression filtering")

coldata <- as.data.frame(meta)
rownames(coldata) <- coldata$sample_id
design_matrix <- model.matrix(~ source_cohort + treatment, data=coldata)
if (qr(design_matrix)$rank != ncol(design_matrix)) stop("DESeq2 design matrix is not full rank")

dds <- DESeqDataSetFromMatrix(countData=filtered, colData=coldata, design=~ source_cohort + treatment)
dds <- DESeq(dds)

size_factors <- data.table(
  sample_id=colnames(dds), tissue=tissue_name,
  source_cohort=as.character(colData(dds)$source_cohort),
  treatment=as.character(colData(dds)$treatment),
  biological_replicate=as.character(colData(dds)$biological_replicate),
  replacement_pituitary=as.character(colData(dds)$replacement_pituitary),
  size_factor=as.numeric(sizeFactors(dds))
)
fwrite(size_factors, snakemake@output[["size_factors"]], sep="\t", quote=FALSE)

norm <- counts(dds, normalized=TRUE)
norm_dt <- data.table(gene_id=rownames(norm))
norm_dt <- cbind(norm_dt, as.data.table(norm))
fwrite(norm_dt, snakemake@output[["normalized"]], sep="\t", quote=FALSE)

filter_summary <- data.table(
  tissue=tissue_name,
  genes_before_filter=nrow(count_mat),
  genes_after_filter=nrow(filtered),
  filter_rule=sprintf("rowSums(counts >= %d) >= %d", min_count, min_samples)
)
fwrite(filter_summary, snakemake@output[["filter_summary"]], sep="\t", quote=FALSE)

complete_outputs <- unlist(snakemake@output[["complete"]])
significant_outputs <- unlist(snakemake@output[["significant"]])

for (i in seq_along(contrast_names)) {
  name <- as.character(contrast_names[i])
  numerator <- as.character(contrast_numerators[i])
  denominator <- as.character(contrast_denominators[i])
  result <- results(
    dds,
    contrast=c("treatment", numerator, denominator),
    alpha=padj_cutoff,
    independentFiltering=independent_filtering
  )
  dt <- as.data.table(as.data.frame(result), keep.rownames="gene_id")
  dt[, tissue := tissue_name]
  dt[, contrast := name]
  dt[, contrast_numerator := numerator]
  dt[, contrast_denominator := denominator]
  dt[, significant := fifelse(!is.na(padj) & padj < padj_cutoff & abs(log2FoldChange) >= lfc_cutoff, "yes", "no")]
  dt[, direction := fifelse(significant == "yes" & log2FoldChange >= lfc_cutoff, "up", fifelse(significant == "yes" & log2FoldChange <= -lfc_cutoff, "down", fifelse(is.na(padj), "not_tested_padj_NA", "not_significant")))]
  setcolorder(dt, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "tissue", "contrast", "contrast_numerator", "contrast_denominator", "significant", "direction"))
  setorder(dt, is.na(padj), padj, is.na(pvalue), pvalue)
  complete_path <- complete_outputs[grepl(paste0("/", name, "\\.complete\\.tsv$"), complete_outputs)]
  significant_path <- significant_outputs[grepl(paste0("/", name, "\\.significant\\.tsv$"), significant_outputs)]
  if (length(complete_path) != 1 || length(significant_path) != 1) stop("Cannot resolve Snakemake outputs for contrast: ", name)
  fwrite(dt, complete_path, sep="\t", quote=FALSE, na="NA")
  fwrite(dt[significant == "yes"], significant_path, sep="\t", quote=FALSE, na="NA")
}
