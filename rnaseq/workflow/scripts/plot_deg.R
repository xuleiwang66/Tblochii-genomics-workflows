suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
dt <- fread(snakemake@input[[1]], sep="\t")
padj_cutoff <- as.numeric(snakemake@params[["padj"]])
lfc_cutoff <- as.numeric(snakemake@params[["lfc"]])
dt[, neg_log10_padj := ifelse(is.na(padj), NA_real_, -log10(pmax(padj, 1e-300)))]
dt[, category := factor(ifelse(direction %in% c("up", "down"), direction, "not_significant"), levels=c("up", "down", "not_significant"))]
p <- ggplot(dt, aes(log2FoldChange, neg_log10_padj, color=category)) + geom_point(alpha=0.65, size=0.75, na.rm=TRUE) +
  geom_vline(xintercept=c(-lfc_cutoff, lfc_cutoff), linetype="dashed") + geom_hline(yintercept=-log10(padj_cutoff), linetype="dashed") +
  theme_bw() + labs(x="log2 fold change", y="-log10 adjusted P-value")
ggsave(snakemake@output[[1]], p, width=7.5, height=6.2)
