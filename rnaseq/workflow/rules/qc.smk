rule fastp:
    input:
        r1=lambda wc: fastq(wc, "r1"),
        r2=lambda wc: fastq(wc, "r2")
    output:
        r1="results/clean/{sample}_R1.fastq.gz",
        r2="results/clean/{sample}_R2.fastq.gz",
        json="results/qc/fastp/{sample}.json",
        html="results/qc/fastp/{sample}.html"
    threads:
        int(config["fastp"]["threads"])
    params:
        adapter_r1=config["fastp"]["adapter_sequence_r1"],
        adapter_r2=config["fastp"]["adapter_sequence_r2"],
        q=int(config["fastp"]["qualified_quality_phred"]),
        unqualified=int(config["fastp"]["unqualified_percent_limit"]),
        n_limit=int(config["fastp"]["n_base_limit"]),
        length=int(config["fastp"]["length_required"]),
        compression=int(config["fastp"]["compression"]),
        detect="--detect_adapter_for_pe" if config["fastp"].get("detect_adapter_for_pe", True) else ""
    conda:
        "../envs/qc.yaml"
    shell:
        r'''
        fastp --in1 {input.r1} --in2 {input.r2} \
          --out1 {output.r1} --out2 {output.r2} \
          --json {output.json} --html {output.html} \
          {params.detect} \
          --adapter_sequence {params.adapter_r1} \
          --adapter_sequence_r2 {params.adapter_r2} \
          --qualified_quality_phred {params.q} \
          --unqualified_percent_limit {params.unqualified} \
          --n_base_limit {params.n_limit} \
          --length_required {params.length} \
          --compression {params.compression} \
          --thread {threads}
        '''


rule fastqc_clean:
    input:
        r1="results/clean/{sample}_R1.fastq.gz",
        r2="results/clean/{sample}_R2.fastq.gz"
    output:
        html1="results/qc/fastqc/{sample}/{sample}_R1_fastqc.html",
        zip1="results/qc/fastqc/{sample}/{sample}_R1_fastqc.zip",
        html2="results/qc/fastqc/{sample}/{sample}_R2_fastqc.html",
        zip2="results/qc/fastqc/{sample}/{sample}_R2_fastqc.zip"
    threads:
        int(config["fastqc"]["threads"])
    conda:
        "../envs/qc.yaml"
    params:
        outdir=lambda wc: f"results/qc/fastqc/{wc.sample}"
    shell:
        "mkdir -p {params.outdir} && fastqc --threads {threads} --outdir {params.outdir} {input.r1} {input.r2}"


rule multiqc:
    input:
        fastp=expand("results/qc/fastp/{sample}.json", sample=SAMPLES),
        fastqc=expand("results/qc/fastqc/{sample}/{sample}_R1_fastqc.zip", sample=SAMPLES)
    output:
        "results/qc/multiqc/rnaseq_qc.html"
    conda:
        "../envs/qc.yaml"
    shell:
        "multiqc results/qc/fastp results/qc/fastqc --filename rnaseq_qc.html --outdir results/qc/multiqc --force"


rule summarize_fastp:
    input:
        samples="results/reference/samples.validated.tsv",
        jsons=expand("results/qc/fastp/{sample}.json", sample=SAMPLES)
    output:
        by_sample="results/qc/fastp_qc_by_sample.tsv",
        overall="results/qc/fastp_qc_summary.tsv",
        warnings="results/qc/fastp_qc_warnings.tsv"
    params:
        thresholds=config["qc_warnings"]
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize_fastp.py"
