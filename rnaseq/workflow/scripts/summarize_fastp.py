import csv
import json
import statistics
from pathlib import Path


def read_samples(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def nested(data, *keys, default=0):
    cur = data
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def summarize(samples_path, json_paths, outputs, thresholds):
    samples = read_samples(samples_path)
    json_by_id = {Path(p).stem: p for p in json_paths}
    rows = []
    warnings = []
    for meta in samples:
        sid = meta["sample_id"]
        path = json_by_id.get(sid)
        if not path:
            raise ValueError(f"Missing fastp JSON for {sid}")
        with open(path, encoding="utf-8") as handle:
            data = json.load(handle)
        before_reads = int(nested(data, "summary", "before_filtering", "total_reads"))
        after_reads = int(nested(data, "summary", "after_filtering", "total_reads"))
        pass_rate = after_reads / before_reads if before_reads else 0.0
        q20 = float(nested(data, "summary", "after_filtering", "q20_rate"))
        q30 = float(nested(data, "summary", "after_filtering", "q30_rate"))
        gc = float(nested(data, "summary", "after_filtering", "gc_content"))
        row = {
            "sample_id": sid,
            "source_cohort": meta["source_cohort"],
            "treatment": meta["treatment"],
            "tissue": meta["tissue"],
            "before_reads": before_reads,
            "after_reads": after_reads,
            "pass_rate": pass_rate,
            "q20_rate": q20,
            "q30_rate": q30,
            "gc_content": gc,
            "duplication_rate": float(nested(data, "duplication", "rate", default=0.0)),
            "low_quality_reads": int(nested(data, "filtering_result", "low_quality_reads")),
            "too_many_N_reads": int(nested(data, "filtering_result", "too_many_N_reads")),
            "too_short_reads": int(nested(data, "filtering_result", "too_short_reads")),
            "adapter_trimmed_reads": int(nested(data, "adapter_cutting", "adapter_trimmed_reads")),
            "insert_size_peak": int(nested(data, "insert_size", "peak", default=0)),
        }
        rows.append(row)
        if pass_rate < float(thresholds["min_fastp_pass_rate"]):
            warnings.append((sid, "low_fastp_pass_rate", pass_rate, thresholds["min_fastp_pass_rate"]))
        if q30 < float(thresholds["min_fastp_q30_rate"]):
            warnings.append((sid, "low_fastp_q30_rate", q30, thresholds["min_fastp_q30_rate"]))
        if gc < float(thresholds["min_fastp_gc_content"]) or gc > float(thresholds["max_fastp_gc_content"]):
            warnings.append((sid, "fastp_gc_outside_review_range", gc, f"{thresholds['min_fastp_gc_content']}-{thresholds['max_fastp_gc_content']}"))
    Path(outputs[0]).parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with open(outputs[0], "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    metrics = []
    for field in ["before_reads", "after_reads", "pass_rate", "q20_rate", "q30_rate", "gc_content", "duplication_rate"]:
        values = [float(r[field]) for r in rows]
        metrics.extend([(field + "_median", statistics.median(values)), (field + "_min", min(values)), (field + "_max", max(values))])
    with open(outputs[1], "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t"); writer.writerow(["metric", "value"]); writer.writerows(metrics)
    with open(outputs[2], "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t"); writer.writerow(["sample_id", "warning", "value", "threshold"]); writer.writerows(warnings)


def main_snakemake(sm):
    summarize(str(sm.input.samples), list(map(str, sm.input.jsons)), [str(sm.output.by_sample), str(sm.output.overall), str(sm.output.warnings)], dict(sm.params.thresholds))


if "snakemake" in globals():
    main_snakemake(snakemake)
