#!/usr/bin/env python3
"""Create final within-cohort binary tails and full-cohort Cox inputs."""
import argparse
import csv
import math
from pathlib import Path

import numpy as np


TOLERANCE = 1e-10


def read_rows(path):
    with open(path, encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise SystemExit(f"Cannot write empty table: {path}")
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


parser = argparse.ArgumentParser()
parser.add_argument("--metadata", required=True)
parser.add_argument("--pca", required=True)
parser.add_argument("--psam", required=True)
parser.add_argument("--kinship", required=True)
parser.add_argument("--tail-fraction", type=float, default=0.30)
parser.add_argument("--assignment", required=True)
parser.add_argument("--counts", required=True)
parser.add_argument("--keep", required=True)
parser.add_argument("--binary-model-data", required=True)
parser.add_argument("--binary-kinship", required=True)
parser.add_argument("--cox-model-data", required=True)
parser.add_argument("--cox-kinship", required=True)
args = parser.parse_args()

psam = read_rows(args.psam)
ids = [row.get("IID") or row.get("#IID") for row in psam]
if any(sample_id is None for sample_id in ids) or len(ids) != len(set(ids)):
    raise SystemExit("Invalid or duplicate IID values in PSAM")
metadata = {row["sample_id"]: row for row in read_rows(args.metadata)}
pca = {row["sample_id"]: row for row in read_rows(args.pca)}
missing = [sample_id for sample_id in ids if sample_id not in metadata or sample_id not in pca]
if missing:
    raise SystemExit(f"Missing metadata/PCA samples: {missing[:10]}")

kinship = np.loadtxt(args.kinship)
if kinship.shape != (len(ids), len(ids)):
    raise SystemExit("Kinship dimension mismatch")
if not np.all(np.isfinite(kinship)):
    raise SystemExit("Kinship contains non-finite values")
if not np.allclose(kinship, kinship.T, atol=1e-8, rtol=0):
    raise SystemExit("Kinship is not symmetric")

assignment_by_id = {}
count_rows = []
selected_set = set()

cohort_order = [cohort for cohort in ["C", "L", "Q"] if cohort in {metadata[s]["source_cohort"] for s in ids}]
cohort_order += sorted({metadata[s]["source_cohort"] for s in ids} - set(cohort_order))

for cohort in cohort_order:
    current = [sample_id for sample_id in ids if metadata[sample_id]["source_cohort"] == cohort]
    current.sort(key=lambda sample_id: (float(metadata[sample_id]["LOE_time_h"]), sample_id))
    cohort_n = len(current)
    tail_size = math.floor(args.tail_fraction * cohort_n)
    if tail_size < 1 or 2 * tail_size > cohort_n:
        raise SystemExit(f"Invalid binary tail size for cohort {cohort}: n={cohort_n}, tail={tail_size}")

    status = {sample_id: None for sample_id in current}
    reason = {sample_id: "middle_excluded" for sample_id in current}
    for sample_id in current[:tail_size]:
        status[sample_id] = 0
        reason[sample_id] = "within_cohort_early_tail"
    for sample_id in current[cohort_n - tail_size:]:
        status[sample_id] = 1
        reason[sample_id] = "within_cohort_late_tail"

    early_cutoff = float(metadata[current[tail_size - 1]]["LOE_time_h"])
    late_cutoff = float(metadata[current[cohort_n - tail_size]]["LOE_time_h"])

    early_at_boundary = [sample_id for sample_id in current if abs(float(metadata[sample_id]["LOE_time_h"]) - early_cutoff) <= TOLERANCE]
    late_at_boundary = [sample_id for sample_id in current if abs(float(metadata[sample_id]["LOE_time_h"]) - late_cutoff) <= TOLERANCE]
    early_tie_crosses = any(status[sample_id] == 0 for sample_id in early_at_boundary) and any(status[sample_id] is None for sample_id in early_at_boundary)
    late_tie_crosses = any(status[sample_id] == 1 for sample_id in late_at_boundary) and any(status[sample_id] is None for sample_id in late_at_boundary)

    if early_tie_crosses:
        for sample_id in early_at_boundary:
            status[sample_id] = None
            reason[sample_id] = "early_boundary_tie_excluded"
    if late_tie_crosses:
        for sample_id in late_at_boundary:
            status[sample_id] = None
            reason[sample_id] = "late_boundary_tie_excluded"

    early_n = sum(value == 0 for value in status.values())
    late_n = sum(value == 1 for value in status.values())
    excluded_n = cohort_n - early_n - late_n
    count_rows.append({
        "source_cohort": cohort,
        "downstream_n": cohort_n,
        "target_tail_size": tail_size,
        "early_n": early_n,
        "late_n": late_n,
        "excluded_n": excluded_n,
        "early_cutoff_h": f"{early_cutoff:.12g}",
        "late_cutoff_h": f"{late_cutoff:.12g}",
        "early_boundary_tie_crosses": "YES" if early_tie_crosses else "NO",
        "late_boundary_tie_crosses": "YES" if late_tie_crosses else "NO",
    })

    for rank, sample_id in enumerate(current, 1):
        binary_status = status[sample_id]
        label = "early_sensitive" if binary_status == 0 else ("late_tolerant" if binary_status == 1 else "middle_excluded")
        assignment_by_id[sample_id] = {
            "sample_id": sample_id,
            "source_cohort": cohort,
            "LOE_time_h": metadata[sample_id]["LOE_time_h"],
            "binary_status": "NA" if binary_status is None else binary_status,
            "binary_label": label,
            "within_cohort_rank": rank,
            "tail_size": tail_size,
            "selection_reason": reason[sample_id],
            "full_genotype_order": ids.index(sample_id) + 1,
        }
        if binary_status is not None:
            selected_set.add(sample_id)

assignment_rows = [assignment_by_id[sample_id] for sample_id in ids]
selected_ids = [sample_id for sample_id in ids if sample_id in selected_set]
count_rows.append({
    "source_cohort": "ALL",
    "downstream_n": len(ids),
    "target_tail_size": sum(int(row["target_tail_size"]) for row in count_rows),
    "early_n": sum(row["binary_status"] == 0 for row in assignment_rows),
    "late_n": sum(row["binary_status"] == 1 for row in assignment_rows),
    "excluded_n": sum(row["binary_status"] == "NA" for row in assignment_rows),
    "early_cutoff_h": "NA",
    "late_cutoff_h": "NA",
    "early_boundary_tie_crosses": "YES" if any(row["early_boundary_tie_crosses"] == "YES" for row in count_rows) else "NO",
    "late_boundary_tie_crosses": "YES" if any(row["late_boundary_tie_crosses"] == "YES" for row in count_rows) else "NO",
})

selected_indices = [ids.index(sample_id) for sample_id in selected_ids]
binary_kinship = kinship[np.ix_(selected_indices, selected_indices)]


def model_rows(use_ids, binary=False):
    output = []
    for sample_id in use_ids:
        m = metadata[sample_id]
        q = pca[sample_id]
        covariates = {
            "sample_id": sample_id,
            "body_weight_g": m["body_weight_g"],
            "source_L": m["source_L"],
            "source_Q": m["source_Q"],
            "PC1": q["PC1"],
            "PC2": q["PC2"],
            "PC3": q["PC3"],
        }
        if binary:
            output.append({"sample_id": sample_id, "binary_status": assignment_by_id[sample_id]["binary_status"], **{k: v for k, v in covariates.items() if k != "sample_id"}})
        else:
            output.append({"sample_id": sample_id, "time_to_event_h": m["time_to_event_h"], "event": m["event"], **{k: v for k, v in covariates.items() if k != "sample_id"}})
    return output

write_rows(args.assignment, assignment_rows)
write_rows(args.counts, count_rows)
write_rows(args.binary_model_data, model_rows(selected_ids, binary=True))
write_rows(args.cox_model_data, model_rows(ids, binary=False))
Path(args.keep).parent.mkdir(parents=True, exist_ok=True)
with open(args.keep, "w", encoding="utf-8") as handle:
    for sample_id in selected_ids:
        handle.write(f"{sample_id}\t{sample_id}\n")
np.savetxt(args.binary_kinship, binary_kinship, delimiter="\t", fmt="%.12g")
np.savetxt(args.cox_kinship, kinship, delimiter="\t", fmt="%.12g")
