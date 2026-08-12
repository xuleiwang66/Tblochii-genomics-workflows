rule fastqc:
    input:
        lambda wc: fastq_for(wc)
    output:
        html="results/fastqc/{sample}_{read}_fastqc.html",
        zip="results/fastqc/{sample}_{read}_fastqc.zip"
    wildcard_constraints:
        read="r1|r2"
    threads: config["resources"]["fastqc_threads"]
    conda: "../envs/qc.yaml"
    params:
        outdir="results/fastqc",
        link=lambda wc: f"results/fastqc/{wc.sample}_{wc.read}.fastq.gz"
    shell:
        r"""
        ln -sfn $(realpath {input}) {params.link}
        fastqc --threads {threads} --outdir {params.outdir} {params.link}
        rm -f {params.link}
        """


rule bwa_mem_sort:
    input:
        r1=lambda wc: SAMPLE_INFO[wc.sample]["r1"],
        r2=lambda wc: SAMPLE_INFO[wc.sample]["r2"],
        ref=REF,
        idx=rules.bwa_index.output
    output:
        bam=temp("results/bam/sorted/{sample}.bam"),
        bai=temp("results/bam/sorted/{sample}.bam.bai")
    threads: config["resources"]["alignment_threads"]
    resources:
        mem_mb=16000
    conda: "../envs/alignment.yaml"
    params:
        rg=lambda wc: f"@RG\\tID:{wc.sample}\\tSM:{wc.sample}\\tPL:DNBSEQ\\tLB:{wc.sample}\\tPU:{wc.sample}"
    shell:
        r"""
        bwa mem -M -t {threads} -R '{params.rg}' {input.ref} {input.r1} {input.r2} |
          samtools sort -@ 1 -m 3G -o {output.bam} -
        samtools index {output.bam}
        samtools quickcheck {output.bam}
        """

rule mark_duplicates:
    input:
        bam="results/bam/sorted/{sample}.bam",
        bai="results/bam/sorted/{sample}.bam.bai"
    output:
        bam=temp("results/bam/markdup/{sample}.markdup.bam"),
        bai=temp("results/bam/markdup/{sample}.markdup.bam.bai"),
        metrics="results/bam/markdup/{sample}.dup_metrics.txt"
    resources:
        mem_mb=24000
    conda: "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options '-Xmx20g' MarkDuplicates \
          -I {input.bam} -O {output.bam} -M {output.metrics} \
          --READ_NAME_REGEX null --REMOVE_DUPLICATES false \
          --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX false
        samtools index {output.bam}
        """

rule set_nm_md_uq_tags:
    input:
        bam="results/bam/markdup/{sample}.markdup.bam",
        bai="results/bam/markdup/{sample}.markdup.bam.bai",
        ref=REF,
        fai=REF+".fai",
        dictionary=rules.gatk_dict.output
    output:
        bam="results/bam/final/{sample}.bam",
        bai="results/bam/final/{sample}.bam.bai"
    resources:
        mem_mb=24000
    conda: "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options '-Xmx20g' SetNmMdAndUqTags \
          -R {input.ref} -I {input.bam} -O {output.bam} \
          --CREATE_INDEX false --VALIDATION_STRINGENCY LENIENT
        samtools index {output.bam}
        samtools quickcheck {output.bam}
        """

rule bam_qc:
    input:
        bam="results/bam/final/{sample}.bam",
        bai="results/bam/final/{sample}.bam.bai",
        bed="results/reference/chromosomes.bed",
        dup="results/bam/markdup/{sample}.dup_metrics.txt"
    output:
        tsv="results/qc/bam/{sample}.tsv",
        flagstat=temp("results/qc/bam/{sample}.flagstat.txt"),
        depth=temp("results/qc/bam/{sample}.depth_summary.tsv")
    threads: 2
    conda: "../envs/alignment.yaml"
    shell:
        r"""
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat}
        samtools depth -aa -b {input.bed} {input.bam} |
          awk 'BEGIN{{OFS="\t"; n=0; s=0; x1=0; x3=0; x5=0; x10=0}}
               {{n++; d=$3; s+=d; if(d>=1)x1++; if(d>=3)x3++; if(d>=5)x5++; if(d>=10)x10++}}
               END{{print "mean_depth_chr01_to_chr24",s/n;
                    print "breadth_1x_chr01_to_chr24_pct",100*x1/n;
                    print "breadth_3x_chr01_to_chr24_pct",100*x3/n;
                    print "breadth_5x_chr01_to_chr24_pct",100*x5/n;
                    print "breadth_10x_chr01_to_chr24_pct",100*x10/n}}' > {output.depth}
        python workflow/scripts/bam_qc.py --sample {wildcards.sample} --flagstat {output.flagstat} --depth {output.depth} --dup-metrics {input.dup} --out {output.tsv}
        """

rule summarize_bam_qc:
    input:
        expand("results/qc/bam/{sample}.tsv", sample=SAMPLES)
    output:
        "results/qc/wgs_sample_qc_summary.tsv"
    conda: "../envs/python.yaml"
    shell:
        "python workflow/scripts/concat_tables.py --inputs {input} --out {output}"
