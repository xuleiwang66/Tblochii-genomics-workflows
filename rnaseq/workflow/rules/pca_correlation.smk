rule vst_pca_correlation_by_tissue:
    input:
        counts="results/counts/raw_gene_counts.tsv",
        samples="results/reference/samples.validated.tsv"
    output:
        vst="results/expression/vst/{tissue}.tsv.gz",
        pca="results/qc/pca/{tissue}.tsv",
        correlations="results/qc/correlations/{tissue}.tsv.gz",
        warnings="results/qc/sample_warnings/{tissue}.tsv"
    params:
        source_levels=config["experimental_design"]["source_levels"],
        treatment_levels=config["experimental_design"]["treatment_levels"],
        min_count=int(config["deseq2"]["min_count"]),
        min_samples=int(config["deseq2"]["min_samples"]),
        thresholds=config["qc_warnings"]
    conda:
        "../envs/r_deseq2.yaml"
    script:
        "../scripts/vst_pca_correlation.R"


rule combine_pca_coordinates:
    input:
        expand("results/qc/pca/{tissue}.tsv", tissue=TISSUES)
    output:
        "results/qc/RNAseq_PCA_coordinates.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/combine_tables.py"


rule combine_sample_correlations:
    input:
        expand("results/qc/correlations/{tissue}.tsv.gz", tissue=TISSUES)
    output:
        "results/qc/RNAseq_sample_correlations.tsv.gz"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/combine_tables.py"


rule combine_sample_warnings:
    input:
        expand("results/qc/sample_warnings/{tissue}.tsv", tissue=TISSUES)
    output:
        "results/qc/RNAseq_sample_QC_warnings.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/combine_tables.py"


rule plot_pca_correlation:
    input:
        vst="results/expression/vst/{tissue}.tsv.gz",
        pca="results/qc/pca/{tissue}.tsv",
        samples="results/reference/samples.validated.tsv"
    output:
        pca="results/plots/pca/{tissue}.pdf",
        correlation="results/plots/correlation/{tissue}.pdf"
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_pca_correlation.R"
