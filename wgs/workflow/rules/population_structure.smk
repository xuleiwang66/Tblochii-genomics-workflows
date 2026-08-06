rule pca:
    input:
        pgen="results/plink/ld_pruned.pgen",pvar="results/plink/ld_pruned.pvar",psam="results/plink/ld_pruned.psam"
    output:
        eigenvec="results/population_structure/pca/pca10.eigenvec",
        eigenval="results/population_structure/pca/pca10.eigenval"
    threads: config["resources"]["pca_threads"]
    resources: mem_mb=64000
    conda: "../envs/plink.yaml"
    params:
        n=config["population_structure"]["pca_components"]
    shell:
        "plink2 --pfile results/plink/ld_pruned --chr-set 24 no-x no-y no-mt --pca {params.n} --threads {threads} --out results/population_structure/pca/pca10"

rule attach_pca_metadata:
    input:
        eigenvec="results/population_structure/pca/pca10.eigenvec",
        metadata="results/metadata/phenotypes.tsv"
    output:
        "results/population_structure/pca/pca10.scores.tsv"
    conda: "../envs/python.yaml"
    shell:
        "python workflow/scripts/attach_pca_metadata.py --eigenvec {input.eigenvec} --metadata {input.metadata} --out {output}"

rule make_kinship_bed:
    input:
        bed="results/plink/common_qc.bed",bim="results/plink/common_qc.bim",fam="results/plink/common_qc.fam"
    output:
        bed="results/population_structure/kinship/genotypes.bed",
        bim="results/population_structure/kinship/genotypes.bim",
        fam="results/population_structure/kinship/genotypes.fam",
        order="results/population_structure/kinship/sample_order.tsv"
    conda: "../envs/plink.yaml"
    shell:
        r"""
        cp {input.bed} {output.bed}; cp {input.bim} {output.bim}
        awk 'BEGIN{{OFS="\t"}} {{$6=1; print}}' {input.fam} > {output.fam}
        awk 'BEGIN{{OFS="\t"; print "order_index","FID","IID"}} {{print NR,$1,$2}}' {output.fam} > {output.order}
        """

rule gemma_kinship:
    input:
        bed="results/population_structure/kinship/genotypes.bed",
        bim="results/population_structure/kinship/genotypes.bim",
        fam="results/population_structure/kinship/genotypes.fam"
    output:
        matrix="results/population_structure/kinship/centered_kinship.cXX.txt",
        log="results/population_structure/kinship/centered_kinship.log.txt"
    resources: mem_mb=96000
    conda: "../envs/population.yaml"
    shell:
        r"""
        gemma -bfile results/population_structure/kinship/genotypes -gk 1 \
          -outdir results/population_structure/kinship -o centered_kinship
        """

rule prepare_admixture_input:
    input:
        bed="results/plink/ld_pruned.bed",bim="results/plink/ld_pruned.bim",fam="results/plink/ld_pruned.fam"
    output:
        bed="results/population_structure/admixture/input/genotypes.bed",
        bim="results/population_structure/admixture/input/genotypes.bim",
        fam="results/population_structure/admixture/input/genotypes.fam",
        order="results/population_structure/admixture/sample_order.tsv"
    shell:
        r"""
        mkdir -p results/population_structure/admixture/input
        ln -sfn $(realpath {input.bed}) {output.bed}; ln -sfn $(realpath {input.bim}) {output.bim}; ln -sfn $(realpath {input.fam}) {output.fam}
        awk 'BEGIN{{OFS="\t"; print "order_index","FID","IID"}} {{print NR,$1,$2}}' {input.fam} > {output.order}
        """

rule admixture:
    input:
        bed="results/population_structure/admixture/input/genotypes.bed",
        bim="results/population_structure/admixture/input/genotypes.bim",
        fam="results/population_structure/admixture/input/genotypes.fam"
    output:
        q="results/population_structure/admixture/input/genotypes.{K}.Q",
        p=temp("results/population_structure/admixture/input/genotypes.{K}.P"),
        log="results/population_structure/admixture/logs/K{K}.log"
    threads: config["resources"]["admixture_threads"]
    resources: mem_mb=32000
    conda: "../envs/population.yaml"
    params: cv=config["population_structure"]["admixture_cv_folds"]
    shell:
        r"""
        mkdir -p results/population_structure/admixture/logs
        cd results/population_structure/admixture/input
        admixture --cv={params.cv} -j{threads} genotypes.bed {wildcards.K} > ../logs/K{wildcards.K}.log 2>&1
        """

rule summarize_admixture:
    input:
        logs=expand("results/population_structure/admixture/logs/K{K}.log",K=K_VALUES)
    output:
        "results/population_structure/admixture/cv_error_by_k.tsv"
    conda: "../envs/python.yaml"
    shell:
        "python workflow/scripts/summarize_admixture.py --logs {input.logs} --out {output}"
