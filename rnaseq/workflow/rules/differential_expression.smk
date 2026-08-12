rule deseq2_by_tissue:
    input:
        counts="results/counts/raw_gene_counts.tsv",
        samples="results/reference/samples.validated.tsv"
    output:
        complete=expand("results/deseq2/{{tissue}}/{contrast}.complete.tsv", contrast=CONTRAST_NAMES),
        significant=expand("results/deseq2/{{tissue}}/{contrast}.significant.tsv", contrast=CONTRAST_NAMES),
        size_factors="results/deseq2/{tissue}/size_factors.tsv",
        normalized="results/deseq2/{tissue}/normalized_counts.tsv.gz",
        filter_summary="results/deseq2/{tissue}/gene_filter_summary.tsv"
    params:
        source_levels=config["experimental_design"]["source_levels"],
        treatment_levels=config["experimental_design"]["treatment_levels"],
        contrast_names=CONTRAST_NAMES,
        contrast_numerators=[x["numerator"] for x in CONTRASTS],
        contrast_denominators=[x["denominator"] for x in CONTRASTS],
        min_count=int(config["deseq2"]["min_count"]),
        min_samples=int(config["deseq2"]["min_samples"]),
        padj=float(config["deseq2"]["padj_cutoff"]),
        lfc=float(config["deseq2"]["absolute_log2fc_cutoff"]),
        independent=bool(config["deseq2"].get("independent_filtering", True))
    conda:
        "../envs/r_deseq2.yaml"
    script:
        "../scripts/deseq2_by_tissue.R"


rule summarize_deseq2:
    input:
        complete=[
            f"results/deseq2/{tissue}/{contrast}.complete.tsv"
            for tissue in TISSUES for contrast in CONTRAST_NAMES
        ]
    output:
        "results/qc/DESeq2_DEG_summary.tsv"
    params:
        tissues=TISSUES,
        contrasts=CONTRAST_NAMES
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize_deseq2.py"
