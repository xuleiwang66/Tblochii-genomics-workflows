suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

pca <- fread(snakemake@input[["pca"]], sep="\t")
samples <- fread(snakemake@input[["samples"]], sep="\t")
tissue_name <- snakemake@wildcards[["tissue"]]
p <- ggplot(pca, aes(PC1, PC2, color=treatment, shape=source_cohort)) +
  geom_point(size=3, alpha=0.9) + theme_bw() +
  labs(x=sprintf("PC1 (%.2f%%)", pca$PC1_percent_variance[1]), y=sprintf("PC2 (%.2f%%)", pca$PC2_percent_variance[1]), title=tissue_name)
ggsave(snakemake@output[["pca"]], p, width=7.5, height=6.2)

vst <- fread(snakemake@input[["vst"]], sep="\t")
genes <- vst$gene_id; vst[, gene_id := NULL]
mat <- as.matrix(vst)
cor_mat <- cor(mat, method="pearson")
ann <- as.data.frame(samples[tissue == tissue_name, .(source_cohort, treatment, biological_replicate, replacement_pituitary)])
rownames(ann) <- samples[tissue == tissue_name, sample_id]
ann <- ann[colnames(cor_mat), , drop=FALSE]
pheatmap(cor_mat, annotation_col=ann, annotation_row=ann, filename=snakemake@output[["correlation"]], border_color=NA, width=9, height=8)
