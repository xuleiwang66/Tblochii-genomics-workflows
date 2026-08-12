rule star_align:
    input:
        index="results/reference/star_index",
        r1="results/clean/{sample}_R1.fastq.gz",
        r2="results/clean/{sample}_R2.fastq.gz"
    output:
        bam="results/alignment/{sample}.sorted.bam",
        bai="results/alignment/{sample}.sorted.bam.bai",
        log="results/alignment/logs/{sample}.Log.final.out",
        sj="results/alignment/logs/{sample}.SJ.out.tab"
    threads:
        int(config["star_alignment"]["threads"])
    params:
        prefix=lambda wc: f"results/alignment/star_tmp/{wc.sample}.",
        multimap=int(config["star_alignment"]["out_filter_multimap_nmax"]),
        sj=int(config["star_alignment"]["align_sj_overhang_min"]),
        sjdb=int(config["star_alignment"]["align_sjdb_overhang_min"]),
        mismatch=float(config["star_alignment"]["out_filter_mismatch_nover_lmax"]),
        sort_ram=int(config["star_alignment"]["limit_bam_sort_ram"])
    conda:
        "../envs/mapping.yaml"
    shell:
        r'''
        mkdir -p results/alignment/logs results/alignment/star_tmp
        STAR --runThreadN {threads} \
          --genomeDir {input.index} \
          --readFilesIn {input.r1} {input.r2} \
          --readFilesCommand zcat \
          --outFileNamePrefix {params.prefix} \
          --outSAMtype BAM SortedByCoordinate \
          --outSAMattributes NH HI AS nM MD \
          --outFilterMultimapNmax {params.multimap} \
          --alignSJoverhangMin {params.sj} \
          --alignSJDBoverhangMin {params.sjdb} \
          --outFilterMismatchNoverLmax {params.mismatch} \
          --limitBAMsortRAM {params.sort_ram}
        mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
        mv {params.prefix}Log.final.out {output.log}
        mv {params.prefix}SJ.out.tab {output.sj}
        samtools index -@ {threads} {output.bam} {output.bai}
        samtools quickcheck -v {output.bam}
        '''


rule summarize_star:
    input:
        samples="results/reference/samples.validated.tsv",
        logs=expand("results/alignment/logs/{sample}.Log.final.out", sample=SAMPLES)
    output:
        by_sample="results/qc/star_mapping_by_sample.tsv",
        overall="results/qc/star_mapping_summary.tsv",
        warnings="results/qc/star_mapping_warnings.tsv"
    params:
        thresholds=config["qc_warnings"]
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize_star.py"


rule infer_strandedness:
    input:
        bam="results/alignment/{sample}.sorted.bam",
        bed="results/reference/annotation.rseqc.bed12"
    output:
        "results/qc/strandedness/{sample}.txt"
    params:
        sample_size=int(config["strandedness"]["sample_size"]),
        mapq=int(config["strandedness"]["min_mapq"])
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} -s {params.sample_size} -q {params.mapq} > {output}"


rule summarize_strandedness:
    input:
        samples="results/reference/samples.validated.tsv",
        reports=expand("results/qc/strandedness/{sample}.txt", sample=SAMPLES)
    output:
        by_sample="results/qc/strandedness_by_sample.tsv",
        summary="results/qc/strandedness_summary.tsv"
    params:
        expected=int(config["strandedness"]["configured_featurecounts_s"]),
        minimum_consensus=float(config["strandedness"]["minimum_consensus"])
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize_strandedness.py"
