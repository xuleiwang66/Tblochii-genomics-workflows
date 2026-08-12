rule harmonize_gwas_candidates:
    input:
        bim="results/plink/common_qc.bim",
        continuous="results/gwas/continuous/candidate_variants.tsv",
        binary="results/gwas/binary/candidate_variants.tsv",
        cox="results/gwas/cox/candidate_variants.tsv",
        binary_effects="results/gwas/binary/candidate_effects.tsv",
        cox_effects="results/gwas/cox/candidate_effects.tsv"
    output:
        "results/gwas/integrated/harmonized_candidate_variants.tsv"
    conda: "../envs/python.yaml"
    shell:
        "python workflow/scripts/harmonize_candidates.py --bim {input.bim} --candidate-files {input.continuous} {input.binary} {input.cox} --effect-files {input.continuous} {input.binary_effects} {input.cox_effects} --models continuous binary cox --out {output}"

rule integrate_cross_model_loci:
    input:
        continuous_leads="results/gwas/continuous/lead_snps.tsv",
        binary_leads="results/gwas/binary/lead_snps.tsv",
        cox_leads="results/gwas/cox/lead_snps.tsv",
        continuous_members="results/gwas/continuous/locus_members.tsv",
        binary_members="results/gwas/binary/locus_members.tsv",
        cox_members="results/gwas/cox/locus_members.tsv"
    output:
        support="results/gwas/integrated/union_lead_support_matrix.tsv",
        cluster_members="results/gwas/integrated/cross_model_locus_cluster_members.tsv",
        clusters="results/gwas/integrated/cross_model_locus_clusters.tsv"
    conda: "../envs/python.yaml"
    params:
        kb=config["analysis"]["clump_kb"]
    shell:
        r"""
        python workflow/scripts/integrate_gwas.py \
          --lead-files {input.continuous_leads} {input.binary_leads} {input.cox_leads} \
          --member-files {input.continuous_members} {input.binary_members} {input.cox_members} \
          --models continuous binary cox --window-kb {params.kb} \
          --support {output.support} --cluster-members {output.cluster_members} --clusters {output.clusters}
        """

rule annotate_cross_model_candidates:
    input:
        clusters="results/gwas/integrated/cross_model_locus_clusters.tsv",
        cluster_members="results/gwas/integrated/cross_model_locus_cluster_members.tsv",
        gff=GFF
    output:
        "results/gwas/integrated/cross_model_candidate_genes.tsv"
    conda: "../envs/python.yaml"
    params:
        functional=FUNCTIONAL_ANNOTATION
    shell:
        r"""
        python workflow/scripts/annotate_clusters.py --clusters {input.clusters} --cluster-members {input.cluster_members} \
          --gff3 {input.gff} --functional '{params.functional}' --out {output}
        """
