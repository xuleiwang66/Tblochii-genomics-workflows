import csv
import statistics
from pathlib import Path


def read_samples(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_log(path):
    values = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if "|" in line:
                key, value = line.split("|", 1)
                values[key.strip()] = value.strip()
    def number(key, percent=False):
        raw = values.get(key, "0").replace("%", "")
        val = float(raw)
        return val / 100.0 if percent else val
    return {
        "input_reads": int(number("Number of input reads")),
        "uniquely_mapped_reads": int(number("Uniquely mapped reads number")),
        "uniquely_mapped_rate": number("Uniquely mapped reads %", True),
        "multimapped_reads": int(number("Number of reads mapped to multiple loci")),
        "multimapped_rate": number("% of reads mapped to multiple loci", True),
        "too_many_loci_reads": int(number("Number of reads mapped to too many loci")),
        "too_many_loci_rate": number("% of reads mapped to too many loci", True),
        "unmapped_mismatch_rate": number("% of reads unmapped: too many mismatches", True),
        "unmapped_too_short_rate": number("% of reads unmapped: too short", True),
        "unmapped_other_rate": number("% of reads unmapped: other", True),
        "mismatch_rate_per_base": number("Mismatch rate per base, %", True),
    }


def summarize(samples_path, log_paths, outputs, thresholds):
    samples = read_samples(samples_path)
    log_by_id = {Path(p).name.replace(".Log.final.out", ""): p for p in log_paths}
    rows = []
    warnings = []
    for meta in samples:
        sid = meta["sample_id"]
        if sid not in log_by_id:
            raise ValueError(f"Missing STAR Log.final.out for {sid}")
        row = {"sample_id": sid, "source_cohort": meta["source_cohort"], "treatment": meta["treatment"], "tissue": meta["tissue"]}
        row.update(parse_log(log_by_id[sid]))
        rows.append(row)
        u = row["uniquely_mapped_rate"]
        m = row["multimapped_rate"]
        short = row["unmapped_too_short_rate"]
        if u < float(thresholds["min_unique_mapping_rate"]):
            warnings.append((sid, "low_unique_mapping_rate", u, thresholds["min_unique_mapping_rate"]))
        elif u < float(thresholds["review_unique_mapping_rate"]):
            warnings.append((sid, "review_unique_mapping_rate", u, thresholds["review_unique_mapping_rate"]))
        if m > float(thresholds["max_multimapping_rate"]):
            warnings.append((sid, "high_multimapping_rate", m, thresholds["max_multimapping_rate"]))
        if short > float(thresholds["max_unmapped_too_short_rate"]):
            warnings.append((sid, "high_unmapped_too_short_rate", short, thresholds["max_unmapped_too_short_rate"]))
    fields = list(rows[0])
    Path(outputs[0]).parent.mkdir(parents=True, exist_ok=True)
    with open(outputs[0], "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows(rows)
    metrics = []
    for field in ["input_reads", "uniquely_mapped_rate", "multimapped_rate", "unmapped_too_short_rate", "mismatch_rate_per_base"]:
        vals = [float(r[field]) for r in rows]
        metrics.extend([(field + "_median", statistics.median(vals)), (field + "_min", min(vals)), (field + "_max", max(vals))])
    with open(outputs[1], "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t"); writer.writerow(["metric", "value"]); writer.writerows(metrics)
    with open(outputs[2], "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t"); writer.writerow(["sample_id", "warning", "value", "threshold"]); writer.writerows(warnings)


def main_snakemake(sm):
    summarize(str(sm.input.samples), list(map(str, sm.input.logs)), [str(sm.output.by_sample), str(sm.output.overall), str(sm.output.warnings)], dict(sm.params.thresholds))


if "snakemake" in globals():
    main_snakemake(snakemake)
