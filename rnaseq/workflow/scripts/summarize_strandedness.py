import csv
import re
from pathlib import Path

PATTERN = re.compile(r'Fraction of reads explained by "([^"]+)":\s*([0-9.eE+-]+)')
FAILED = re.compile(r'Fraction of reads failed to determine:\s*([0-9.eE+-]+)')


def classify(path):
    text = Path(path).read_text(encoding="utf-8")
    failed_match = FAILED.search(text)
    groups = PATTERN.findall(text)
    failed = float(failed_match.group(1)) if failed_match else 1.0
    if len(groups) < 2:
        return failed, 0.0, 0.0, "failed", None
    first = float(groups[0][1]); second = float(groups[1][1])
    dominant = max(first, second)
    difference = abs(first - second)
    if failed > 0.50:
        return failed, first, second, "high_failed_to_determine", None
    if dominant < 0.60 or difference < 0.20:
        return failed, first, second, "unstranded_or_ambiguous", 0
    if first > second:
        return failed, first, second, "stranded_group_1", 1
    return failed, first, second, "stranded_group_2", 2


def summarize(samples_path, reports, by_sample, summary, expected, minimum_consensus):
    with open(samples_path, encoding="utf-8-sig", newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    report_by_id = {Path(p).stem: p for p in reports}
    rows = []
    recommendations = []
    for meta in samples:
        sid = meta["sample_id"]
        if sid not in report_by_id:
            raise ValueError(f"Missing strandedness report for {sid}")
        failed, group1, group2, category, recommendation = classify(report_by_id[sid])
        rows.append([sid, meta["tissue"], failed, group1, group2, category, "NA" if recommendation is None else recommendation])
        if recommendation is not None:
            recommendations.append(recommendation)
    counts = {x: recommendations.count(x) for x in [0, 1, 2]}
    consensus_value = max(counts, key=counts.get) if recommendations else 0
    consensus_fraction = counts[consensus_value] / len(recommendations) if recommendations else 0.0
    global_recommendation = consensus_value if consensus_fraction >= minimum_consensus else 0
    if global_recommendation != expected:
        raise ValueError(f"Configured featureCounts -s {expected} conflicts with RSeQC recommendation {global_recommendation}")
    Path(by_sample).parent.mkdir(parents=True, exist_ok=True)
    with open(by_sample, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t"); writer.writerow(["sample_id", "tissue", "failed_fraction", "group1_fraction", "group2_fraction", "classification", "recommended_featureCounts_s"]); writer.writerows(rows)
    with open(summary, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t"); writer.writerow(["metric", "value"])
        writer.writerows([
            ("configured_featureCounts_s", expected),
            ("recommended_featureCounts_s", global_recommendation),
            ("consensus_fraction", consensus_fraction),
            ("samples_recommending_0", counts[0]),
            ("samples_recommending_1", counts[1]),
            ("samples_recommending_2", counts[2]),
        ])


def main_snakemake(sm):
    summarize(str(sm.input.samples), list(map(str, sm.input.reports)), str(sm.output.by_sample), str(sm.output.summary), int(sm.params.expected), float(sm.params.minimum_consensus))


if "snakemake" in globals():
    main_snakemake(snakemake)
