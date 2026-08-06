#!/usr/bin/env python3
"""Annotate each model-specific locus, then aggregate gene support by cluster."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_attributes(text):
    attributes = {}
    for token in text.split(";"):
        if "=" in token:
            key, value = token.split("=", 1)
            attributes[key] = value
    return attributes


def distance_to_gene(position, gene):
    if gene["start"] <= position <= gene["end"]:
        return 0
    return min(abs(position - gene["start"]), abs(position - gene["end"]))


parser = argparse.ArgumentParser()
parser.add_argument("--clusters", required=True)
parser.add_argument("--cluster-members", required=True)
parser.add_argument("--gff3", required=True)
parser.add_argument("--functional", default="")
parser.add_argument("--out", required=True)
args = parser.parse_args()

genes = []
with open(args.gff3, encoding="utf-8", errors="replace") as handle:
    for line in handle:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] not in {"gene", "pseudogene"}:
            continue
        attributes = parse_attributes(fields[8])
        gene_id = attributes.get("ID") or attributes.get("gene_id") or attributes.get("Name")
        if gene_id:
            genes.append({
                "chromosome": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "gene_id": gene_id,
                "gene_name": attributes.get("Name", gene_id),
            })

genes_by_chromosome = defaultdict(list)
for gene in genes:
    genes_by_chromosome[gene["chromosome"]].append(gene)

functional = {}
functional_fields = []
if args.functional and Path(args.functional).exists():
    with open(args.functional, encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames:
            key = next((name for name in ["query", "gene_id", "ID", "#query"] if name in reader.fieldnames), reader.fieldnames[0])
            functional_fields = [name for name in reader.fieldnames if name != key]
            for row in reader:
                functional[row.get(key, "")] = row

with open(args.clusters, encoding="utf-8-sig") as handle:
    clusters = {row["cluster_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
with open(args.cluster_members, encoding="utf-8-sig") as handle:
    loci = list(csv.DictReader(handle, delimiter="\t"))

# Support is first evaluated per model-specific locus, matching the analysis definition.
records = {}
for locus in loci:
    chromosome = locus["chromosome"]
    lead_position = int(float(locus["position"]))
    window_start = int(float(locus["window_start"]))
    window_end = int(float(locus["window_end"]))
    clump_start = int(float(locus["clump_start"]))
    clump_end = int(float(locus["clump_end"]))
    chromosome_genes = genes_by_chromosome.get(chromosome, [])

    window_genes = [gene for gene in chromosome_genes if gene["start"] <= window_end and gene["end"] >= window_start]
    if chromosome_genes:
        distances = [distance_to_gene(lead_position, gene) for gene in chromosome_genes]
        minimum_distance = min(distances)
        nearest_genes = [gene for gene, distance in zip(chromosome_genes, distances) if distance == minimum_distance]
    else:
        nearest_genes = []

    selected = {gene["gene_id"]: gene for gene in window_genes + nearest_genes}
    for gene_id, gene in selected.items():
        key = (locus["cluster_id"], gene_id)
        record = records.setdefault(key, {
            "cluster_id": locus["cluster_id"],
            "chromosome": chromosome,
            "cluster_start": clusters[locus["cluster_id"]]["start"],
            "cluster_end": clusters[locus["cluster_id"]]["end"],
            "locus_models": clusters[locus["cluster_id"]]["models"],
            "gene_id": gene_id,
            "gene_name": gene["gene_name"],
            "gene_start": gene["start"],
            "gene_end": gene["end"],
            "window_models": set(),
            "clump_span_models": set(),
            "nearest_models": set(),
            "nearest_distance_bp": None,
            "supporting_loci": set(),
        })
        record["supporting_loci"].add(locus["source_locus_id"])
        if gene in window_genes:
            record["window_models"].add(locus["model"])
        if gene["start"] <= clump_end and gene["end"] >= clump_start:
            record["clump_span_models"].add(locus["model"])
        if gene in nearest_genes:
            record["nearest_models"].add(locus["model"])
            distance = distance_to_gene(lead_position, gene)
            if record["nearest_distance_bp"] is None or distance < record["nearest_distance_bp"]:
                record["nearest_distance_bp"] = distance

rows = []
for (_, gene_id), record in sorted(records.items(), key=lambda item: (item[0][0], item[1]["gene_start"], item[0][1])):
    row = {
        "cluster_id": record["cluster_id"],
        "chromosome": record["chromosome"],
        "cluster_start": record["cluster_start"],
        "cluster_end": record["cluster_end"],
        "locus_models": record["locus_models"],
        "gene_id": gene_id,
        "gene_name": record["gene_name"],
        "gene_start": record["gene_start"],
        "gene_end": record["gene_end"],
        "window_support_models": ",".join(sorted(record["window_models"])) or "NA",
        "clump_span_support_models": ",".join(sorted(record["clump_span_models"])) or "NA",
        "nearest_gene_support_models": ",".join(sorted(record["nearest_models"])) or "NA",
        "minimum_nearest_distance_bp": "NA" if record["nearest_distance_bp"] is None else record["nearest_distance_bp"],
        "supporting_model_loci": ",".join(sorted(record["supporting_loci"])),
    }
    annotation = functional.get(gene_id, {})
    for field in functional_fields:
        row[field] = annotation.get(field, "NA") or "NA"
    rows.append(row)

base_fields = [
    "cluster_id", "chromosome", "cluster_start", "cluster_end", "locus_models",
    "gene_id", "gene_name", "gene_start", "gene_end", "window_support_models",
    "clump_span_support_models", "nearest_gene_support_models", "minimum_nearest_distance_bp",
    "supporting_model_loci",
]
fields = base_fields + functional_fields
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
with open(args.out, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
