#!/usr/bin/env python3
"""Integrate model-specific PLINK clumps using fixed lead-centered windows."""
import argparse
import csv
from pathlib import Path


def read_rows(path):
    with open(path, encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def chromosome_key(value):
    text = str(value)
    stripped = text.lower().replace("chr", "")
    try:
        return (0, int(stripped))
    except ValueError:
        return (1, text)


parser = argparse.ArgumentParser()
parser.add_argument("--lead-files", nargs="+", required=True)
parser.add_argument("--member-files", nargs="+", required=True)
parser.add_argument("--models", nargs="+", required=True)
parser.add_argument("--window-kb", type=int, default=500)
parser.add_argument("--support", required=True)
parser.add_argument("--cluster-members", required=True)
parser.add_argument("--clusters", required=True)
args = parser.parse_args()

if not (len(args.lead_files) == len(args.member_files) == len(args.models)):
    raise SystemExit("lead-files, member-files and models must have equal lengths")

window_bp = args.window_kb * 1000
loci = []
all_members = []
for model, lead_file, member_file in zip(args.models, args.lead_files, args.member_files):
    for row in read_rows(lead_file):
        position = int(float(row["position"]))
        loci.append({
            "model": model,
            "source_locus_id": row.get("source_locus_id") or f"{model}:{row['lead_snp']}",
            "chromosome": row["chromosome"],
            "position": position,
            "lead_snp": row["lead_snp"],
            "lead_p": row["lead_p"],
            "window_start": int(float(row.get("window_start") or max(1, position - window_bp))),
            "window_end": int(float(row.get("window_end") or position + window_bp)),
            "clump_start": int(float(row.get("clump_start") or position)),
            "clump_end": int(float(row.get("clump_end") or position)),
            "member_count": int(float(row.get("member_count") or 1)),
        })
    for row in read_rows(member_file):
        row["model"] = model
        all_members.append(row)

# Exact lead-SNP support.
by_lead = {}
for locus in loci:
    by_lead.setdefault(locus["lead_snp"], []).append(locus)
support_rows = []
for snp, records in sorted(by_lead.items(), key=lambda item: (chromosome_key(item[1][0]["chromosome"]), item[1][0]["position"], item[0])):
    models = sorted({record["model"] for record in records})
    support_rows.append({
        "lead_snp": snp,
        "chromosome": records[0]["chromosome"],
        "position": records[0]["position"],
        "models": ",".join(models),
        "model_count": len(models),
    })

# Connected components of overlapping fixed lead-centered windows.
cluster_id = 0
cluster_members = []
cluster_rows = []
for chromosome in sorted({locus["chromosome"] for locus in loci}, key=chromosome_key):
    chromosome_loci = sorted(
        [locus for locus in loci if locus["chromosome"] == chromosome],
        key=lambda row: (row["window_start"], row["window_end"], row["model"], row["lead_snp"]),
    )
    current = []
    current_end = None
    components = []
    for locus in chromosome_loci:
        if not current or locus["window_start"] <= current_end:
            current.append(locus)
            current_end = locus["window_end"] if current_end is None else max(current_end, locus["window_end"])
        else:
            components.append(current)
            current = [locus]
            current_end = locus["window_end"]
    if current:
        components.append(current)

    for component in components:
        cluster_id += 1
        cid = f"L{cluster_id:04d}"
        models = sorted({locus["model"] for locus in component})
        leads = sorted({locus["lead_snp"] for locus in component})
        cluster_rows.append({
            "cluster_id": cid,
            "chromosome": chromosome,
            "start": min(locus["window_start"] for locus in component),
            "end": max(locus["window_end"] for locus in component),
            "models": ",".join(models),
            "model_count": len(models),
            "locus_count": len(component),
            "lead_snp_count": len(leads),
            "lead_snps": ",".join(leads),
            "minimum_lead_p": min(float(locus["lead_p"]) for locus in component),
        })
        for locus in component:
            cluster_members.append({"cluster_id": cid, **locus})

outputs = [
    (args.support, support_rows, ["lead_snp", "chromosome", "position", "models", "model_count"]),
    (
        args.cluster_members,
        cluster_members,
        [
            "cluster_id", "model", "source_locus_id", "chromosome", "position", "lead_snp",
            "lead_p", "window_start", "window_end", "clump_start", "clump_end", "member_count",
        ],
    ),
    (
        args.clusters,
        cluster_rows,
        [
            "cluster_id", "chromosome", "start", "end", "models", "model_count",
            "locus_count", "lead_snp_count", "lead_snps", "minimum_lead_p",
        ],
    ),
]

for output_path, rows, fields in outputs:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
