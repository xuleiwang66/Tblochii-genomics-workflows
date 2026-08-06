import csv
from pathlib import Path


def summarize(paths, output):
    rows = []
    for path in paths:
        p = Path(path)
        tissue = p.parent.name
        contrast = p.name.replace(".complete.tsv", "")
        tested = significant = up = down = padj_na = 0
        with open(path, encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                tested += 1
                if row.get("padj", "NA") in {"NA", "", "nan"}:
                    padj_na += 1
                if row.get("significant") == "yes":
                    significant += 1
                    if row.get("direction") == "up": up += 1
                    if row.get("direction") == "down": down += 1
        rows.append([tissue, contrast, tested, significant, up, down, padj_na])
    rows.sort()
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t"); writer.writerow(["tissue", "contrast", "tested_genes", "significant_DEGs", "up_DEGs", "down_DEGs", "padj_NA"]); writer.writerows(rows)


def main_snakemake(sm):
    summarize(list(map(str, sm.input.complete)), str(sm.output[0]))


if "snakemake" in globals():
    main_snakemake(snakemake)
