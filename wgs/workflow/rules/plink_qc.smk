rule vcf_to_plink2:
    input:
        vcf="results/genotype_refinement/beagle/{chrom}.vcf.gz",
        csi="results/genotype_refinement/beagle/{chrom}.vcf.gz.csi"
    output:
        pgen="results/plink/by_chromosome/{chrom}.pgen",
        pvar="results/plink/by_chromosome/{chrom}.pvar",
        psam="results/plink/by_chromosome/{chrom}.psam"
    threads: 2
    resources: mem_mb=16000
    conda: "../envs/plink.yaml"
    params:
        prefix=lambda wc:f"results/plink/by_chromosome/{wc.chrom}"
    shell:
        r"""
        plink2 --vcf {input.vcf} --double-id --allow-extra-chr \
          --chr-set 24 no-x no-y no-mt \
          --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 missing \
          --make-pgen --threads {threads} --out {params.prefix}
        """

rule make_pmerge_list:
    input:
        pgen=expand("results/plink/by_chromosome/{chrom}.pgen",chrom=CHROMOSOMES),
        chromosomes=config["chromosomes"]
    output:
        "results/plink/pmerge_list.txt"
    conda: "../envs/python.yaml"
    shell:
        "python workflow/scripts/make_pmerge_list.py --chromosomes {input.chromosomes} --prefix-pattern 'results/plink/by_chromosome/{{chrom}}' --out {output}"

rule merge_plink2:
    input:
        merge_list="results/plink/pmerge_list.txt",
        pgen=expand("results/plink/by_chromosome/{chrom}.pgen",chrom=CHROMOSOMES),
        pvar=expand("results/plink/by_chromosome/{chrom}.pvar",chrom=CHROMOSOMES),
        psam=expand("results/plink/by_chromosome/{chrom}.psam",chrom=CHROMOSOMES)
    output:
        pgen="results/plink/raw.pgen",
        pvar="results/plink/raw.pvar",
        psam="results/plink/raw.psam"
    threads: config["resources"]["plink_threads"]
    resources: mem_mb=128000
    conda: "../envs/plink.yaml"
    params: prefix="results/plink/raw"
    shell:
        "plink2 --pmerge-list {input.merge_list} pfile --sample-inner-join --chr-set 24 no-x no-y no-mt --sort-vars --make-pgen --threads {threads} --out {params.prefix}"

rule plink_common_qc:
    input:
        pgen="results/plink/raw.pgen", pvar="results/plink/raw.pvar", psam="results/plink/raw.psam"
    output:
        missing_pgen=temp("results/plink/missingness_qc.pgen"),
        missing_pvar=temp("results/plink/missingness_qc.pvar"),
        missing_psam=temp("results/plink/missingness_qc.psam"),
        pgen="results/plink/common_qc.pgen",
        pvar="results/plink/common_qc.pvar",
        psam="results/plink/common_qc.psam"
    threads: config["resources"]["plink_threads"]
    resources: mem_mb=128000
    conda: "../envs/plink.yaml"
    params:
        mind=config["plink_qc"]["mind"], geno=config["plink_qc"]["geno"], maf=config["plink_qc"]["maf"]
    shell:
        r"""
        plink2 --pfile results/plink/raw --chr-set 24 no-x no-y no-mt \
          --mind {params.mind} --geno {params.geno} --make-pgen --threads {threads} --out results/plink/missingness_qc
        plink2 --pfile results/plink/missingness_qc --chr-set 24 no-x no-y no-mt \
          --maf {params.maf} --make-pgen --threads {threads} --out results/plink/common_qc
        """

rule ld_prune:
    input:
        pgen="results/plink/common_qc.pgen", pvar="results/plink/common_qc.pvar", psam="results/plink/common_qc.psam"
    output:
        prune_in="results/plink/ld_pruning.prune.in",
        prune_out="results/plink/ld_pruning.prune.out",
        pgen="results/plink/ld_pruned.pgen",
        pvar="results/plink/ld_pruned.pvar",
        psam="results/plink/ld_pruned.psam"
    threads: config["resources"]["plink_threads"]
    resources: mem_mb=128000
    conda: "../envs/plink.yaml"
    params:
        window=config["plink_qc"]["ld_window"], step=config["plink_qc"]["ld_step"], r2=config["plink_qc"]["ld_r2"]
    shell:
        r"""
        plink2 --pfile results/plink/common_qc --chr-set 24 no-x no-y no-mt \
          --indep-pairwise {params.window} {params.step} {params.r2} --threads {threads} --out results/plink/ld_pruning
        plink2 --pfile results/plink/common_qc --chr-set 24 no-x no-y no-mt \
          --extract {output.prune_in} --make-pgen --threads {threads} --out results/plink/ld_pruned
        """

rule export_ld_pruned_bed:
    input:
        pgen="results/plink/ld_pruned.pgen",pvar="results/plink/ld_pruned.pvar",psam="results/plink/ld_pruned.psam"
    output:
        bed="results/plink/ld_pruned.bed",bim="results/plink/ld_pruned.bim",fam="results/plink/ld_pruned.fam"
    threads: 4
    conda: "../envs/plink.yaml"
    shell:
        "plink2 --pfile results/plink/ld_pruned --chr-set 24 no-x no-y no-mt --make-bed --threads {threads} --out results/plink/ld_pruned"

rule export_common_bed:
    input:
        pgen="results/plink/common_qc.pgen",pvar="results/plink/common_qc.pvar",psam="results/plink/common_qc.psam"
    output:
        bed="results/plink/common_qc.bed",bim="results/plink/common_qc.bim",fam="results/plink/common_qc.fam"
    threads: 4
    conda: "../envs/plink.yaml"
    shell:
        r"""
        plink2 --pfile results/plink/common_qc --chr-set 24 no-x no-y no-mt --make-bed --threads {threads} --out results/plink/common_qc
        """

rule summarize_plink_qc:
    input:
        raw="results/plink/raw.pgen", missing="results/plink/missingness_qc.pgen", common="results/plink/common_qc.pgen", pruned="results/plink/ld_pruned.pgen"
    output:
        "results/qc/plink_qc_summary.tsv"
    conda: "../envs/python.yaml"
    shell:
        "python workflow/scripts/summarize_plink_qc.py --raw results/plink/raw --missingness results/plink/missingness_qc --common results/plink/common_qc --pruned results/plink/ld_pruned --out {output}"
