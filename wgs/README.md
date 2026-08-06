# WGS Snakemake workflow

Minimal workflow for the *Trachinotus blochii* WGS analysis: FASTQ QC, BWA-MEM alignment, GATK variant calling/filtering, Beagle refinement, PLINK QC, PCA, ADMIXTURE, kinship, continuous/binary/Cox GWAS and cross-model candidate integration.

## Inputs

- `config/samples.tsv`: `sample_id`, `r1`, `r2`
- `config/phenotypes.tsv`: `sample_id`, `fish_id`, `source_cohort`, `LOE_time_h`, `LOE_status`, `body_weight_g`
- `resources/Tov.fa`
- `resources/Tov_rename.gff3`
- `resources/beagle.5.5.jar`
- optional functional annotation table

`LOE_status`: `1` event, `0` right-censored at 24 h, `2` pre-target death excluded from downstream genotype/GWAS analyses.

## Run

```bash
snakemake --use-conda --cores 1 --dry-run
snakemake --use-conda --cores 32
```

## Main parameters

- GATK filters: `QD<2`, `FS>60`, `SOR>3`, `MQ<40`, `MQRankSum<-12.5`, `ReadPosRankSum<-8`
- site call rate: `>=0.80`
- Beagle: `window=2`, `overlap=0.2`, no external panel, no DR2 filter
- PLINK: `--mind 0.1 --geno 0.1 --maf 0.05`
- LD pruning: `--indep-pairwise 50 10 0.2`
- ADMIXTURE: `K=2–60`, 10-fold CV
- GWAS covariates: body weight, source cohort, PC1–PC3 and kinship
- candidate threshold: `P<=1e-5`; clumping: `p1=p2=1e-5`, `r2=0.2`, `500 kb`

Confirmed original-analysis versions are listed in `software_versions.tsv`.
