rule prepare_metadata:
    input:
        phenotypes=config["phenotypes"],
        samples=config["samples"]
    output:
        clean="results/metadata/phenotypes.tsv",
        analysis="results/metadata/analysis_samples.txt",
        excluded="results/metadata/pre_target_death_samples.txt"
    conda:
        "../envs/python.yaml"
    params:
        censor=config["analysis"]["censor_time_h"],
        exclude=" ".join(str(x) for x in config["analysis"]["downstream_status_exclude"])
    shell:
        """
        python workflow/scripts/prepare_metadata.py \
          --phenotypes {input.phenotypes} \
          --samples {input.samples} \
          --clean {output.clean} \
          --analysis-samples {output.analysis} \
          --pre-target-deaths {output.excluded} \
          --censor-time {params.censor} \
          --exclude-status {params.exclude}
        """

rule samtools_faidx:
    input: REF
    output: REF + ".fai"
    conda: "../envs/alignment.yaml"
    shell: "samtools faidx {input}"

rule gatk_dict:
    input: REF
    output: str(Path(REF).with_suffix(".dict"))
    conda: "../envs/gatk.yaml"
    shell: "gatk CreateSequenceDictionary -R {input} -O {output}"

rule bwa_index:
    input: REF
    output:
        amb=REF+".amb", ann=REF+".ann", bwt=REF+".bwt", pac=REF+".pac", sa=REF+".sa"
    conda: "../envs/alignment.yaml"
    shell: "bwa index {input}"

rule reference_intervals:
    input:
        fai=REF+".fai",
        chromosomes=config["chromosomes"]
    output:
        bed="results/reference/chromosomes.bed",
        validated="results/reference/reference_contigs.validated.tsv"
    conda: "../envs/python.yaml"
    shell:
        "python workflow/scripts/make_reference_intervals.py --fai {input.fai} --chromosomes {input.chromosomes} --bed {output.bed} --validated {output.validated}"
