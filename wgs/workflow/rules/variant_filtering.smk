rule select_biallelic_snps:
    input:
        vcf="results/joint/raw/{chrom}.raw.vcf.gz",
        index="results/joint/raw/{chrom}.raw.vcf.gz.csi",
        ref=REF
    output:
        vcf=temp("results/variants/biallelic/{chrom}.vcf.gz"),
        csi=temp("results/variants/biallelic/{chrom}.vcf.gz.csi")
    resources: mem_mb=64000
    conda: "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options '-Xmx48g' SelectVariants -R {input.ref} -V {input.vcf} -L {wildcards.chrom} \
          --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O {output.vcf}
        bcftools index --force --csi {output.vcf}
        """

rule hard_filter_snps:
    input:
        vcf="results/variants/biallelic/{chrom}.vcf.gz",
        csi="results/variants/biallelic/{chrom}.vcf.gz.csi",
        ref=REF
    output:
        vcf=temp("results/variants/hard_filtered/{chrom}.vcf.gz"),
        csi=temp("results/variants/hard_filtered/{chrom}.vcf.gz.csi")
    resources: mem_mb=48000
    conda: "../envs/gatk.yaml"
    params:
        qd=config["hard_filter"]["QD"], fs=config["hard_filter"]["FS"], sor=config["hard_filter"]["SOR"],
        mq=config["hard_filter"]["MQ"], mqrs=config["hard_filter"]["MQRankSum"], rprs=config["hard_filter"]["ReadPosRankSum"]
    shell:
        r"""
        gatk --java-options '-Xmx32g' VariantFiltration -R {input.ref} -V {input.vcf} -L {wildcards.chrom} -O {output.vcf} \
          --missing-values-evaluate-as-failing false \
          --filter-name QD_lt_2 --filter-expression 'QD < {params.qd}' \
          --filter-name FS_gt_60 --filter-expression 'FS > {params.fs}' \
          --filter-name SOR_gt_3 --filter-expression 'SOR > {params.sor}' \
          --filter-name MQ_lt_40 --filter-expression 'MQ < {params.mq}' \
          --filter-name MQRankSum_lt_neg12_5 --filter-expression 'MQRankSum < {params.mqrs}' \
          --filter-name ReadPosRankSum_lt_neg8 --filter-expression 'ReadPosRankSum < {params.rprs}'
        bcftools index --force --csi {output.vcf}
        """

rule select_pass_snps:
    input:
        vcf="results/variants/hard_filtered/{chrom}.vcf.gz",
        csi="results/variants/hard_filtered/{chrom}.vcf.gz.csi",
        ref=REF
    output:
        vcf=temp("results/variants/pass/{chrom}.vcf.gz"),
        csi=temp("results/variants/pass/{chrom}.vcf.gz.csi")
    resources: mem_mb=48000
    conda: "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options '-Xmx16g' SelectVariants -R {input.ref} -V {input.vcf} -L {wildcards.chrom} --exclude-filtered -O {output.vcf}
        bcftools index --force --csi {output.vcf}
        """

rule filter_site_call_rate:
    input:
        vcf="results/variants/pass/{chrom}.vcf.gz",
        csi="results/variants/pass/{chrom}.vcf.gz.csi"
    output:
        vcf="results/variants/filtered/{chrom}.vcf.gz",
        csi="results/variants/filtered/{chrom}.vcf.gz.csi"
    threads: 2
    resources: mem_mb=16000
    conda: "../envs/vcf.yaml"
    params:
        min_an=MIN_AN
    shell:
        r"""
        test {params.min_an} -gt 0 || (echo 'samples.tsv contains no samples' >&2; exit 1)
        bcftools view -i 'INFO/AN>={params.min_an}' -Oz -o {output.vcf} {input.vcf}
        bcftools index --force --csi --threads {threads} {output.vcf}
        """

rule summarize_variant_filtering:
    input:
        raw=expand("results/joint/raw/{chrom}.raw.vcf.gz",chrom=CHROMOSOMES),
        biallelic=expand("results/variants/biallelic/{chrom}.vcf.gz",chrom=CHROMOSOMES),
        passed=expand("results/variants/pass/{chrom}.vcf.gz",chrom=CHROMOSOMES),
        final=expand("results/variants/filtered/{chrom}.vcf.gz",chrom=CHROMOSOMES),
        chromosomes=config["chromosomes"]
    output:
        by_chr="results/qc/variant_filtering_by_chromosome.tsv",
        summary="results/qc/variant_filtering_summary.tsv"
    conda: "../envs/vcf.yaml"
    shell:
        r"""
        python workflow/scripts/summarize_variant_filtering.py --chromosomes {input.chromosomes} \
          --raw-pattern 'results/joint/raw/{{chrom}}.raw.vcf.gz' \
          --biallelic 'results/variants/biallelic/{{chrom}}.vcf.gz' \
          --pass-pattern 'results/variants/pass/{{chrom}}.vcf.gz' \
          --final 'results/variants/filtered/{{chrom}}.vcf.gz' \
          --out-by-chr {output.by_chr} --out-summary {output.summary}
        """
