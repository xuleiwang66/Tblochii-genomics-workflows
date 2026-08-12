rule validate_samples:
    input:
        config["samples"]
    output:
        "results/reference/samples.validated.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/validate_samples.py"


rule fasta_index:
    input:
        REF
    output:
        REF + ".fai"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools faidx {input}"


rule gff3_to_gtf:
    input:
        gff3=GFF3,
        chromosomes=config["chromosomes"]
    output:
        "results/reference/annotation.for_star.gtf"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/gff3_to_gtf.py"


rule star_index:
    input:
        fasta=REF,
        gtf="results/reference/annotation.for_star.gtf"
    output:
        directory("results/reference/star_index")
    threads:
        int(config["star_index"]["threads"])
    params:
        sjdb_overhang=int(config["star_index"]["sjdb_overhang"]),
        sa_index=int(config["star_index"]["genome_sa_index_nbases"]),
        chr_bin=int(config["star_index"]["genome_chr_bin_nbits"]),
        ram=int(config["star_index"]["limit_genome_generate_ram"])
    conda:
        "../envs/mapping.yaml"
    shell:
        r'''
        mkdir -p {output}
        STAR --runThreadN {threads} \
          --runMode genomeGenerate \
          --genomeDir {output} \
          --genomeFastaFiles {input.fasta} \
          --sjdbGTFfile {input.gtf} \
          --sjdbOverhang {params.sjdb_overhang} \
          --genomeSAindexNbases {params.sa_index} \
          --genomeChrBinNbits {params.chr_bin} \
          --limitGenomeGenerateRAM {params.ram}
        '''


rule prepare_gene_exon_saf:
    input:
        gff3=GFF3,
        fai=REF + ".fai",
        chromosomes=config["chromosomes"]
    output:
        saf="results/reference/gene_exon_union.saf",
        lengths="results/reference/gene_length.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/prepare_gene_exon_saf.py"


rule annotation_to_bed12:
    input:
        gtf="results/reference/annotation.for_star.gtf",
        fai=REF + ".fai"
    output:
        "results/reference/annotation.rseqc.bed12"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/annotation_to_bed12.py"
