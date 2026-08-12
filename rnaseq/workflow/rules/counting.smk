rule featurecounts:
    input:
        saf="results/reference/gene_exon_union.saf",
        bams=expand("results/alignment/{sample}.sorted.bam", sample=SAMPLES),
        strandedness="results/qc/strandedness_summary.tsv" if config["run"].get("strandedness_qc", True) else []
    output:
        counts="results/counts/featureCounts.txt",
        summary="results/counts/featureCounts.txt.summary"
    threads:
        int(config["featurecounts"]["threads"])
    params:
        s=int(config["featurecounts"]["strandedness"]),
        both="-B" if config["featurecounts"].get("require_both_ends_mapped", True) else "",
        chimeric="" if config["featurecounts"].get("count_chimeric_fragments", True) else "-C",
        multimapping="-M" if config["featurecounts"].get("count_multimapping", False) else "",
        overlap="-O" if config["featurecounts"].get("allow_multi_overlap", False) else ""
    conda:
        "../envs/counting.yaml"
    shell:
        r'''
        featureCounts -F SAF -a {input.saf} -o {output.counts} \
          -s {params.s} -p --countReadPairs {params.both} {params.chimeric} \
          {params.multimapping} {params.overlap} -T {threads} {input.bams}
        test -s {output.summary}
        '''


rule collect_featurecounts:
    input:
        samples="results/reference/samples.validated.tsv",
        counts="results/counts/featureCounts.txt",
        summary="results/counts/featureCounts.txt.summary"
    output:
        matrix="results/counts/raw_gene_counts.tsv",
        by_sample="results/qc/featurecounts_by_sample.tsv",
        overall="results/qc/featurecounts_summary.tsv",
        warnings="results/qc/featurecounts_warnings.tsv"
    params:
        thresholds=config["qc_warnings"]
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/collect_featurecounts.py"
