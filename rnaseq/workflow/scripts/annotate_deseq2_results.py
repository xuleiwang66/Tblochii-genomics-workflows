import csv
import math
from pathlib import Path
from annotation_utils import ANNOTATION_COLUMNS, read_eggnog


def annotate(input_path, annotation_path, complete_out, significant_out):
    annotations = read_eggnog(annotation_path)
    Path(complete_out).parent.mkdir(parents=True, exist_ok=True)
    with open(input_path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        base_fields = list(reader.fieldnames or [])
        fields = base_fields + ["plot_category", "neg_log10_padj", "log10_baseMean_plus1"] + ANNOTATION_COLUMNS + ["annotation_status"]
        with open(complete_out, "w", encoding="utf-8", newline="") as complete, open(significant_out, "w", encoding="utf-8", newline="") as sig:
            wc = csv.DictWriter(complete, fieldnames=fields, delimiter="\t"); ws = csv.DictWriter(sig, fieldnames=fields, delimiter="\t")
            wc.writeheader(); ws.writeheader()
            seen = set()
            for row in reader:
                gene = row["gene_id"]
                if gene in seen:
                    raise ValueError(f"Duplicated gene_id: {gene}")
                seen.add(gene)
                padj = None if row.get("padj") in {None, "", "NA", "nan"} else float(row["padj"])
                base_mean = float(row.get("baseMean") or 0)
                row["plot_category"] = row.get("direction", "not_significant")
                row["neg_log10_padj"] = "NA" if padj is None else f"{-math.log10(max(padj, 1e-300)):.8g}"
                row["log10_baseMean_plus1"] = f"{math.log10(base_mean + 1.0):.8g}"
                ann = annotations.get(gene, {})
                for col in ANNOTATION_COLUMNS:
                    row[col] = ann.get(col, "")
                row["annotation_status"] = "annotated" if any(row[col] for col in ANNOTATION_COLUMNS) else "unannotated"
                wc.writerow(row)
                if row.get("significant") == "yes":
                    ws.writerow(row)


def main_snakemake(sm):
    annotate(str(sm.input.complete), str(sm.input.annotation), str(sm.output.complete), str(sm.output.significant))


if "snakemake" in globals():
    main_snakemake(snakemake)
