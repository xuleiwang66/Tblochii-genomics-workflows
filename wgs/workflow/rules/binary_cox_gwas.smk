rule prepare_binary_cox_inputs:
    input:
        metadata="results/metadata/phenotypes.tsv",
        pca="results/population_structure/pca/pca10.scores.tsv",
        psam="results/plink/common_qc.psam",
        kinship="results/population_structure/kinship/centered_kinship.cXX.txt"
    output:
        assignment="results/gwas/binary/input/assignment.tsv",
        counts="results/gwas/binary/input/group_counts.tsv",
        keep="results/gwas/binary/input/keep.tsv",
        binary_data="results/gwas/binary/input/model_data.tsv",
        binary_kinship="results/gwas/binary/input/kinship.txt",
        cox_data="results/gwas/cox/input/model_data.tsv",
        cox_kinship="results/gwas/cox/input/kinship.txt"
    conda: "../envs/scipy.yaml"
    params: fraction=config["analysis"]["binary_tail_fraction"]
    shell:
        r"""
        python workflow/scripts/prepare_binary_cox_inputs.py --metadata {input.metadata} --pca {input.pca} --psam {input.psam} --kinship {input.kinship} \
          --tail-fraction {params.fraction} --assignment {output.assignment} --counts {output.counts} --keep {output.keep} \
          --binary-model-data {output.binary_data} --binary-kinship {output.binary_kinship} --cox-model-data {output.cox_data} --cox-kinship {output.cox_kinship}
        """

rule subset_binary_genotypes:
    input:
        pgen="results/plink/common_qc.pgen",pvar="results/plink/common_qc.pvar",psam="results/plink/common_qc.psam",keep="results/gwas/binary/input/keep.tsv"
    output:
        pgen="results/gwas/binary/input/genotypes.pgen",pvar="results/gwas/binary/input/genotypes.pvar",psam="results/gwas/binary/input/genotypes.psam",
        bed="results/gwas/binary/input/genotypes.bed",bim="results/gwas/binary/input/genotypes.bim",fam="results/gwas/binary/input/genotypes.fam"
    threads: 8
    resources: mem_mb=64000
    conda: "../envs/plink.yaml"
    shell:
        r"""
        plink2 --pfile results/plink/common_qc --keep {input.keep} --make-pgen --threads {threads} --out results/gwas/binary/input/genotypes
        plink2 --pfile results/gwas/binary/input/genotypes --chr-set 24 no-x no-y no-mt --make-bed --threads {threads} --out results/gwas/binary/input/genotypes
        """

rule fit_binary_null:
    input:
        model="results/gwas/binary/input/model_data.tsv",kinship="results/gwas/binary/input/kinship.txt"
    output: "results/gwas/binary/null_model.rds"
    resources: mem_mb=64000
    conda: "../envs/gwas.yaml"
    shell: "Rscript workflow/scripts/fit_binary_null.R {input.model} {input.kinship} {output}"

rule binary_score_gwas:
    input:
        null="results/gwas/binary/null_model.rds",bed="results/gwas/binary/input/genotypes.bed",bim="results/gwas/binary/input/genotypes.bim",fam="results/gwas/binary/input/genotypes.fam"
    output: "results/gwas/binary/score.raw.tsv"
    resources: mem_mb=96000
    conda: "../envs/gwas.yaml"
    shell: "Rscript workflow/scripts/run_binary_score.R {input.null} results/gwas/binary/input/genotypes {output}"

rule binary_gwas_qc:
    input: "results/gwas/binary/score.raw.tsv"
    output:
        standardized="results/gwas/binary/score.standardized.tsv.gz",candidates="results/gwas/binary/candidate_variants.tsv",summary="results/gwas/binary/qc_summary.tsv"
    conda: "../envs/scipy.yaml"
    params: threshold=config["analysis"]["candidate_p"]
    shell: "python workflow/scripts/gwas_qc.py --model binary --input {input} --out {output.standardized} --candidates {output.candidates} --summary {output.summary} --threshold {params.threshold}"

rule binary_candidate_effects:
    input:
        candidates="results/gwas/binary/candidate_variants.tsv",model="results/gwas/binary/input/model_data.tsv",kinship="results/gwas/binary/input/kinship.txt",bed="results/gwas/binary/input/genotypes.bed"
    output: "results/gwas/binary/candidate_effects.tsv"
    conda: "../envs/gwas.yaml"
    shell: "Rscript workflow/scripts/binary_candidate_effects.R {input.model} {input.kinship} results/gwas/binary/input/genotypes {input.candidates} {output}"

rule binary_clump_input:
    input: "results/gwas/binary/candidate_variants.tsv"
    output: "results/gwas/binary/clump_input.tsv"
    conda: "../envs/python.yaml"
    shell: "python workflow/scripts/make_clump_input.py --candidates {input} --out {output}"

rule binary_clump:
    input:
        bed="results/gwas/binary/input/genotypes.bed",bim="results/gwas/binary/input/genotypes.bim",fam="results/gwas/binary/input/genotypes.fam",assoc="results/gwas/binary/clump_input.tsv"
    output: "results/gwas/binary/clumps.clumped"
    conda: "../envs/plink.yaml"
    params: p1=config["analysis"]["clump_p1"],p2=config["analysis"]["clump_p2"],r2=config["analysis"]["clump_r2"],kb=config["analysis"]["clump_kb"]
    shell: "plink --bfile results/gwas/binary/input/genotypes --allow-extra-chr --clump {input.assoc} --clump-snp-field SNP --clump-field P --clump-p1 {params.p1} --clump-p2 {params.p2} --clump-r2 {params.r2} --clump-kb {params.kb} --out results/gwas/binary/clumps"

rule parse_binary_clumps:
    input: "results/gwas/binary/clumps.clumped"
    output: leads="results/gwas/binary/lead_snps.tsv",members="results/gwas/binary/locus_members.tsv"
    conda: "../envs/python.yaml"
    params: kb=config["analysis"]["clump_kb"]
    shell: "python workflow/scripts/parse_clumps.py --clumped {input} --model binary --window-kb {params.kb} --leads {output.leads} --members {output.members}"

rule fit_cox_null:
    input: model="results/gwas/cox/input/model_data.tsv",kinship="results/gwas/cox/input/kinship.txt"
    output: model="results/gwas/cox/null_model.rds",tau="results/gwas/cox/null_tau.tsv"
    resources: mem_mb=96000
    conda: "../envs/gwas.yaml"
    shell: "Rscript workflow/scripts/fit_cox_null.R {input.model} {input.kinship} {output.model} {output.tau}"

rule cox_score_by_chromosome:
    input:
        bed="results/plink/common_qc.bed",bim="results/plink/common_qc.bim",fam="results/plink/common_qc.fam",model="results/gwas/cox/input/model_data.tsv",kinship="results/gwas/cox/input/kinship.txt",null="results/gwas/cox/null_model.rds"
    output: "results/gwas/cox/by_chromosome/{chrom}.tsv"
    resources: mem_mb=64000
    conda: "../envs/gwas.yaml"
    params: work=lambda wc:f"results/gwas/cox/work/{wc.chrom}"
    shell: "Rscript workflow/scripts/run_cox_score_by_chr.R {wildcards.chrom} results/plink/common_qc {input.model} {input.kinship} {input.null} {params.work} {output} $(command -v plink2)"

rule merge_cox_score:
    input: expand("results/gwas/cox/by_chromosome/{chrom}.tsv",chrom=CHROMOSOMES)
    output: "results/gwas/cox/score.raw.tsv"
    conda: "../envs/gwas.yaml"
    shell: "Rscript workflow/scripts/merge_cox_results.R {output} {input}"

rule cox_gwas_qc:
    input: "results/gwas/cox/score.raw.tsv"
    output: standardized="results/gwas/cox/score.standardized.tsv.gz",candidates="results/gwas/cox/candidate_variants.tsv",summary="results/gwas/cox/qc_summary.tsv"
    conda: "../envs/scipy.yaml"
    params: threshold=config["analysis"]["candidate_p"]
    shell: "python workflow/scripts/gwas_qc.py --model cox --input {input} --out {output.standardized} --candidates {output.candidates} --summary {output.summary} --threshold {params.threshold}"

rule cox_candidate_effects:
    input: candidates="results/gwas/cox/candidate_variants.tsv",model="results/gwas/cox/input/model_data.tsv",kinship="results/gwas/cox/input/kinship.txt",null="results/gwas/cox/null_model.rds",bed="results/plink/common_qc.bed"
    output: "results/gwas/cox/candidate_effects.tsv"
    resources: mem_mb=96000
    conda: "../envs/gwas.yaml"
    params: work="results/gwas/cox/candidate_effect_work"
    shell: "Rscript workflow/scripts/cox_candidate_effects.R {input.candidates} {input.model} {input.kinship} {input.null} results/plink/common_qc {params.work} {output}"

rule cox_clump_input:
    input: "results/gwas/cox/candidate_variants.tsv"
    output: "results/gwas/cox/clump_input.tsv"
    conda: "../envs/python.yaml"
    shell: "python workflow/scripts/make_clump_input.py --candidates {input} --out {output}"

rule cox_clump:
    input: bed="results/plink/common_qc.bed",bim="results/plink/common_qc.bim",fam="results/plink/common_qc.fam",assoc="results/gwas/cox/clump_input.tsv"
    output: "results/gwas/cox/clumps.clumped"
    conda: "../envs/plink.yaml"
    params: p1=config["analysis"]["clump_p1"],p2=config["analysis"]["clump_p2"],r2=config["analysis"]["clump_r2"],kb=config["analysis"]["clump_kb"]
    shell: "plink --bfile results/plink/common_qc --allow-extra-chr --clump {input.assoc} --clump-snp-field SNP --clump-field P --clump-p1 {params.p1} --clump-p2 {params.p2} --clump-r2 {params.r2} --clump-kb {params.kb} --out results/gwas/cox/clumps"

rule parse_cox_clumps:
    input: "results/gwas/cox/clumps.clumped"
    output: leads="results/gwas/cox/lead_snps.tsv",members="results/gwas/cox/locus_members.tsv"
    conda: "../envs/python.yaml"
    params: kb=config["analysis"]["clump_kb"]
    shell: "python workflow/scripts/parse_clumps.py --clumped {input} --model cox --window-kb {params.kb} --leads {output.leads} --members {output.members}"
