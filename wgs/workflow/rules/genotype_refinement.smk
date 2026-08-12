rule subset_analysis_samples:
    input:
        vcf="results/variants/filtered/{chrom}.vcf.gz",
        csi="results/variants/filtered/{chrom}.vcf.gz.csi",
        samples="results/metadata/analysis_samples.txt"
    output:
        vcf=temp("results/genotype_refinement/input/{chrom}.vcf.gz"),
        csi=temp("results/genotype_refinement/input/{chrom}.vcf.gz.csi")
    threads: 4
    conda: "../envs/vcf.yaml"
    shell:
        r"""
        bcftools view --threads {threads} -S {input.samples} -Oz -o {output.vcf} {input.vcf}
        bcftools index --force --csi --threads {threads} {output.vcf}
        """

rule beagle_phase_and_impute:
    input:
        vcf="results/genotype_refinement/input/{chrom}.vcf.gz",
        csi="results/genotype_refinement/input/{chrom}.vcf.gz.csi",
        jar=config["beagle"]["jar"]
    output:
        vcf="results/genotype_refinement/beagle/{chrom}.vcf.gz",
        csi="results/genotype_refinement/beagle/{chrom}.vcf.gz.csi",
        log="results/genotype_refinement/beagle/{chrom}.log"
    threads: config["beagle"]["threads"]
    resources: mem_mb=32000
    conda: "../envs/beagle.yaml"
    params:
        prefix=lambda wc: f"results/genotype_refinement/beagle/{wc.chrom}",
        mem=config["beagle"]["java_mem_gb"],
        window=config["beagle"]["window"],
        overlap=config["beagle"]["overlap"]
    shell:
        r"""
        java -Xmx{params.mem}g -jar {input.jar} gt={input.vcf} out={params.prefix} nthreads={threads} window={params.window} overlap={params.overlap}
        bcftools index --force --csi --threads {threads} {output.vcf}
        """
