rule prepare_continuous_gwas_inputs:
    input:
        psam="results/plink/common_qc.psam",
        metadata="results/metadata/phenotypes.tsv",
        pca="results/population_structure/pca/pca10.scores.tsv",
        kinship_order="results/population_structure/kinship/sample_order.tsv"
    output:
        phenotype="results/gwas/continuous/input/phenotype.txt",
        covariates="results/gwas/continuous/input/covariates.txt",
        annotated="results/gwas/continuous/input/annotated.tsv",
        order="results/gwas/continuous/input/sample_order.tsv"
    conda: "../envs/scipy.yaml"
    shell:
        r"""
        python workflow/scripts/prepare_continuous_gwas_inputs.py --psam {input.psam} --metadata {input.metadata} --pca {input.pca} --kinship-order {input.kinship_order} \
          --phenotype {output.phenotype} --covariates {output.covariates} --annotated {output.annotated} --sample-order {output.order}
        """

rule make_continuous_gwas_bed:
    input:
        bed="results/plink/common_qc.bed",bim="results/plink/common_qc.bim",fam="results/plink/common_qc.fam"
    output:
        bed="results/gwas/continuous/input/genotypes.bed",
        bim="results/gwas/continuous/input/genotypes.bim",
        fam="results/gwas/continuous/input/genotypes.fam"
    shell:
        r"""
        cp {input.bed} {output.bed}; cp {input.bim} {output.bim}
        awk 'BEGIN{{OFS="\t"}} {{$6=1; print}}' {input.fam} > {output.fam}
        """

rule continuous_gwas_lmm:
    input:
        bed="results/gwas/continuous/input/genotypes.bed",
        bim="results/gwas/continuous/input/genotypes.bim",
        fam="results/gwas/continuous/input/genotypes.fam",
        phenotype="results/gwas/continuous/input/phenotype.txt",
        covariates="results/gwas/continuous/input/covariates.txt",
        kinship="results/population_structure/kinship/centered_kinship.cXX.txt"
    output:
        assoc="results/gwas/continuous/main_lmm.assoc.txt",
        log="results/gwas/continuous/main_lmm.log.txt"
    resources: mem_mb=128000
    conda: "../envs/gwas.yaml"
    shell:
        r"""
        gemma -bfile results/gwas/continuous/input/genotypes \
          -p {input.phenotype} -n 1 -c {input.covariates} -k {input.kinship} -km 1 -lmm 4 \
          -outdir results/gwas/continuous -o main_lmm
        """

rule continuous_gwas_qc:
    input:
        assoc="results/gwas/continuous/main_lmm.assoc.txt"
    output:
        standardized="results/gwas/continuous/main_lmm.standardized.tsv.gz",
        candidates="results/gwas/continuous/candidate_variants.tsv",
        summary="results/gwas/continuous/qc_summary.tsv"
    conda: "../envs/scipy.yaml"
    params: threshold=config["analysis"]["candidate_p"]
    shell:
        "python workflow/scripts/gwas_qc.py --model continuous --input {input.assoc} --out {output.standardized} --candidates {output.candidates} --summary {output.summary} --threshold {params.threshold}"

rule continuous_clump_input:
    input: "results/gwas/continuous/candidate_variants.tsv"
    output: "results/gwas/continuous/clump_input.tsv"
    conda: "../envs/python.yaml"
    shell: "python workflow/scripts/make_clump_input.py --candidates {input} --out {output}"

rule continuous_clump:
    input:
        bed="results/gwas/continuous/input/genotypes.bed",bim="results/gwas/continuous/input/genotypes.bim",fam="results/gwas/continuous/input/genotypes.fam",assoc="results/gwas/continuous/clump_input.tsv"
    output:
        clumped="results/gwas/continuous/clumps.clumped"
    conda: "../envs/plink.yaml"
    params:
        p1=config["analysis"]["clump_p1"],p2=config["analysis"]["clump_p2"],r2=config["analysis"]["clump_r2"],kb=config["analysis"]["clump_kb"]
    shell:
        "plink --bfile results/gwas/continuous/input/genotypes --allow-extra-chr --clump {input.assoc} --clump-snp-field SNP --clump-field P --clump-p1 {params.p1} --clump-p2 {params.p2} --clump-r2 {params.r2} --clump-kb {params.kb} --out results/gwas/continuous/clumps"

rule parse_continuous_clumps:
    input: "results/gwas/continuous/clumps.clumped"
    output:
        leads="results/gwas/continuous/lead_snps.tsv",
        members="results/gwas/continuous/locus_members.tsv"
    conda: "../envs/python.yaml"
    params: kb=config["analysis"]["clump_kb"]
    shell: "python workflow/scripts/parse_clumps.py --clumped {input} --model continuous --window-kb {params.kb} --leads {output.leads} --members {output.members}"
