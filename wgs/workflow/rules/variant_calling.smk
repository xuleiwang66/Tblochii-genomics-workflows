rule haplotypecaller_gvcf:
    input:
        bam="results/bam/final/{sample}.bam",
        bai="results/bam/final/{sample}.bam.bai",
        ref=REF,
        fai=REF+".fai",
        dictionary=rules.gatk_dict.output,
        chromosomes=config["chromosomes"]
    output:
        gvcf="results/gvcf/{sample}.g.vcf.gz",
        csi="results/gvcf/{sample}.g.vcf.gz.csi"
    threads: config["resources"]["haplotypecaller_threads"]
    resources:
        mem_mb=16000
    conda: "../envs/gatk.yaml"
    params:
        intervals=" ".join(f"-L {c}" for c in CHROMOSOMES),
        ploidy=PLOIDY
    shell:
        r"""
        gatk --java-options '-Xmx12g' HaplotypeCaller \
          -R {input.ref} -I {input.bam} -O {output.gvcf} \
          -ERC GVCF --sample-ploidy {params.ploidy} \
          --native-pair-hmm-threads {threads} \
          {params.intervals}
        bcftools index --force --csi {output.gvcf}
        """

rule make_gvcf_sample_map:
    input:
        samples=config["samples"],
        gvcfs=expand("results/gvcf/{sample}.g.vcf.gz", sample=SAMPLES)
    output:
        "results/joint/gvcf_sample_map.tsv"
    conda: "../envs/python.yaml"
    shell:
        "python workflow/scripts/make_gvcf_map.py --samples {input.samples} --gvcf-dir results/gvcf --out {output}"

rule genomicsdb_import:
    input:
        sample_map="results/joint/gvcf_sample_map.tsv",
        gvcfs=expand("results/gvcf/{sample}.g.vcf.gz", sample=SAMPLES)
    output:
        workspace=directory("results/joint/genomicsdb/{chrom}")
    threads: config["resources"]["genomicsdb_threads"]
    resources:
        mem_mb=64000
    conda: "../envs/gatk.yaml"
    shell:
        r"""
        rm -rf {output.workspace}
        gatk --java-options '-Xmx32g' GenomicsDBImport \
          --sample-name-map {input.sample_map} \
          --genomicsdb-workspace-path {output.workspace} \
          -L {wildcards.chrom} --reader-threads {threads} --batch-size 50
        """

rule genotype_gvcfs:
    input:
        workspace="results/joint/genomicsdb/{chrom}",
        ref=REF,
        fai=REF+".fai",
        dictionary=rules.gatk_dict.output
    output:
        vcf="results/joint/raw/{chrom}.raw.vcf.gz",
        csi="results/joint/raw/{chrom}.raw.vcf.gz.csi"
    threads: 4
    resources:
        mem_mb=160000
    conda: "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options '-Xmx12g' GenotypeGVCFs \
          -R {input.ref} -V gendb://{input.workspace} -L {wildcards.chrom} \
          -O {output.vcf}
        bcftools index --force --csi --threads {threads} {output.vcf}
        """
