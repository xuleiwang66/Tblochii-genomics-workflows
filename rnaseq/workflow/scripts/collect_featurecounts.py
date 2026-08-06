import csv
import os
import statistics
from collections import defaultdict
from pathlib import Path


def read_samples(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sample_from_column(value, valid_ids):
    base = os.path.basename(value)
    for suffix in [".sorted.bam", ".bam"]:
        if base.endswith(suffix):
            base = base[:-len(suffix)]
    if base in valid_ids:
        return base
    raise ValueError(f"Cannot map featureCounts column to sample_id: {value}")


def collect(samples_path, count_path, summary_path, outputs, thresholds):
    samples = read_samples(samples_path)
    order = [x["sample_id"] for x in samples]
    valid = set(order)
    header = None
    rows = []
    with open(count_path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0] in {"Geneid", "GeneID"}:
                header = parts
            elif header and parts:
                rows.append(parts)
    if header is None:
        raise ValueError("featureCounts header not found")
    fc_ids = [sample_from_column(x, valid) for x in header[6:]]
    if set(fc_ids) != valid or len(fc_ids) != len(set(fc_ids)):
        raise ValueError("featureCounts sample columns do not match samples.tsv")
    index = {sid: i for i, sid in enumerate(fc_ids)}
    gene_ids = []
    counts = {sid: [] for sid in order}
    seen = set()
    for row in rows:
        gene = row[0]
        if gene in seen:
            raise ValueError(f"Duplicated gene_id in featureCounts output: {gene}")
        seen.add(gene); gene_ids.append(gene)
        values = row[6:]
        for sid in order:
            value = float(values[index[sid]])
            if value < 0 or not value.is_integer():
                raise ValueError(f"Invalid raw count for {gene}/{sid}: {value}")
            counts[sid].append(int(value))
    Path(outputs[0]).parent.mkdir(parents=True, exist_ok=True)
    with open(outputs[0], "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t"); writer.writerow(["gene_id"] + order)
        for i, gene in enumerate(gene_ids):
            writer.writerow([gene] + [counts[sid][i] for sid in order])
    with open(summary_path, encoding="utf-8", newline="") as handle:
        summary_rows = list(csv.reader(handle, delimiter="\t"))
    summary_ids = [sample_from_column(x, valid) for x in summary_rows[0][1:]]
    sidx = {sid: i for i, sid in enumerate(summary_ids)}
    categories = [row[0] for row in summary_rows[1:] if row]
    per_sample = {sid: {} for sid in order}
    for row in summary_rows[1:]:
        if not row:
            continue
        for sid in order:
            per_sample[sid][row[0]] = int(float(row[1 + sidx[sid]]))
    group_totals = defaultdict(list)
    group_assigned = defaultdict(list)
    output_rows = []
    meta_by_id = {x["sample_id"]: x for x in samples}
    for sid in order:
        meta = meta_by_id[sid]
        total = sum(per_sample[sid].values())
        assigned = per_sample[sid].get("Assigned", 0)
        key = (meta["source_cohort"], meta["treatment"], meta["tissue"])
        group_totals[key].append(total); group_assigned[key].append(assigned)
        output_rows.append({"sample_id": sid, "source_cohort": key[0], "treatment": key[1], "tissue": key[2], **per_sample[sid], "total_fragments": total, "assigned_rate": assigned / total if total else 0.0})
    med_total = {k: statistics.median(v) for k, v in group_totals.items()}
    med_assigned = {k: statistics.median(v) for k, v in group_assigned.items()}
    warnings = []
    for row in output_rows:
        key = (row["source_cohort"], row["treatment"], row["tissue"])
        row["total_fragments_ratio_to_group_median"] = row["total_fragments"] / med_total[key] if med_total[key] else 0.0
        row["assigned_ratio_to_group_median"] = row.get("Assigned", 0) / med_assigned[key] if med_assigned[key] else 0.0
        rate = row["assigned_rate"]
        ratio = row["total_fragments_ratio_to_group_median"]
        if rate < float(thresholds["min_count_assigned_rate"]):
            warnings.append((row["sample_id"], "low_assigned_rate", rate, thresholds["min_count_assigned_rate"]))
        elif rate < float(thresholds["review_count_assigned_rate"]):
            warnings.append((row["sample_id"], "review_assigned_rate", rate, thresholds["review_count_assigned_rate"]))
        if ratio < float(thresholds["min_group_library_ratio"]) or ratio > float(thresholds["max_group_library_ratio"]):
            warnings.append((row["sample_id"], "library_size_ratio_outside_group_range", ratio, f"{thresholds['min_group_library_ratio']}-{thresholds['max_group_library_ratio']}"))
    fields = ["sample_id", "source_cohort", "treatment", "tissue"] + categories + ["total_fragments", "assigned_rate", "total_fragments_ratio_to_group_median", "assigned_ratio_to_group_median"]
    with open(outputs[1], "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", extrasaction="ignore"); writer.writeheader(); writer.writerows(output_rows)
    rates = [row["assigned_rate"] for row in output_rows]
    with open(outputs[2], "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t"); writer.writerow(["metric", "value"]); writer.writerows([("genes", len(gene_ids)), ("samples", len(order)), ("assigned_rate_median", statistics.median(rates)), ("assigned_rate_min", min(rates)), ("assigned_rate_max", max(rates))])
    with open(outputs[3], "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t"); writer.writerow(["sample_id", "warning", "value", "threshold"]); writer.writerows(warnings)


def main_snakemake(sm):
    collect(str(sm.input.samples), str(sm.input.counts), str(sm.input.summary), [str(sm.output.matrix), str(sm.output.by_sample), str(sm.output.overall), str(sm.output.warnings)], dict(sm.params.thresholds))


if "snakemake" in globals():
    main_snakemake(snakemake)
