# RNA-seq Snakemake workflow

Minimal workflow for the *Trachinotus blochii* RNA-seq analysis: fastp/FastQC/MultiQC, STAR alignment, RSeQC strandedness QC, exon-union featureCounts, TPM, tissue-specific DESeq2, VST, PCA, sample correlations and optional GO/KEGG ORA.

## Inputs

- `config/samples.tsv`: `sample_id`, `r1`, `r2`, `source_cohort`, `treatment`, `tissue`, `biological_replicate`, `replacement_pituitary`
- `resources/Tov.fa`
- `resources/Tov_rename.gff3`
- optional eggNOG-mapper annotation table

Labels: cohorts `C/L/Q`; treatments `control/early/late`; tissues `gill/hypothalamus/liver/pituitary`.

## Run

```bash
snakemake --use-conda --cores 1 --dry-run
snakemake --use-conda --cores 32
```

## Main parameters

- fastp: quality 6; <=50% low-quality bases; <=15 N bases; minimum length 15
- STAR: maximum 20 loci/read; mismatch proportion <=0.04
- RSeQC: `infer_experiment.py -s 200000 -q 30`
- featureCounts: paired fragments, `-s 0`, `-B`, no multimapping or multi-overlap counting
- gene length: gene-level exon union on `Chr01–Chr24`
- DESeq2: per tissue, `~ source_cohort + treatment`
- expression filter: count >=10 in >=3 samples
- contrasts: early/control, late/control and early/late
- DEG: `padj<0.05` and `|log2FC|>=1`
- VST: `blind=FALSE`; PCA centered/unscaled; Pearson sample correlations
- optional ORA: one-sided hypergeometric test with BH correction by tissue × contrast × term type

Confirmed original-analysis versions are listed in `software_versions.tsv`.
