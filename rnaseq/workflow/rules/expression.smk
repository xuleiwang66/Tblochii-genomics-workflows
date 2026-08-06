rule calculate_tpm:
    input:
        counts="results/counts/raw_gene_counts.tsv",
        lengths="results/reference/gene_length.tsv"
    output:
        matrix="results/expression/tpm.tsv",
        sums="results/qc/tpm_sample_sums.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/calculate_tpm.py"
