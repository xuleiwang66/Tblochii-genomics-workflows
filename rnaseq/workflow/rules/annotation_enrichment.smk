rule annotate_deseq2_result:
    input:
        complete="results/deseq2/{tissue}/{contrast}.complete.tsv",
        annotation=lambda wc: FUNCTIONAL_ANNOTATION
    output:
        complete="results/deseq2_annotated/{tissue}/{contrast}.complete.tsv",
        significant="results/deseq2_annotated/{tissue}/{contrast}.significant.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/annotate_deseq2_results.py"


rule run_ora:
    input:
        annotation=lambda wc: FUNCTIONAL_ANNOTATION,
        complete=[
            f"results/deseq2/{tissue}/{contrast}.complete.tsv"
            for tissue in TISSUES for contrast in CONTRAST_NAMES
        ],
        significant=[
            f"results/deseq2/{tissue}/{contrast}.significant.tsv"
            for tissue in TISSUES for contrast in CONTRAST_NAMES
        ]
    output:
        results="results/enrichment/ORA_complete_results.tsv.gz",
        summary="results/enrichment/ORA_summary.tsv",
        go="results/enrichment/TERM2GENE_GO.tsv.gz",
        ko="results/enrichment/TERM2GENE_KEGG_ko.tsv.gz",
        pathway="results/enrichment/TERM2GENE_KEGG_Pathway.tsv.gz"
    params:
        tissues=TISSUES,
        contrasts=CONTRAST_NAMES,
        padj=float(config["enrichment"]["padj_cutoff"])
    conda:
        "../envs/scipy.yaml"
    script:
        "../scripts/run_ora.py"


rule plot_deg_result:
    input:
        "results/deseq2/{tissue}/{contrast}.complete.tsv"
    output:
        "results/plots/deg/{tissue}.{contrast}.volcano.pdf"
    params:
        padj=float(config["deseq2"]["padj_cutoff"]),
        lfc=float(config["deseq2"]["absolute_log2fc_cutoff"])
    conda:
        "../envs/r_plots.yaml"
    script:
        "../scripts/plot_deg.R"
