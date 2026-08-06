#!/usr/bin/env python3
"""Parse PLINK 1.9 clumping output into lead- and member-level tables."""
import argparse
import csv
import re
from pathlib import Path


def parse_variant_id(variant_id: str):
    """Parse IDs produced by PLINK2 --set-all-var-ids '@:#:$r:$a'."""
    match = re.match(r"^([^:]+):(\d+):", variant_id or "")
    if not match:
        return None, None
    return match.group(1), int(match.group(2))


parser = argparse.ArgumentParser()
parser.add_argument("--clumped", required=True)
parser.add_argument("--model", required=True)
parser.add_argument("--window-kb", type=int, default=500)
parser.add_argument("--leads", required=True)
parser.add_argument("--members", required=True)
args = parser.parse_args()

lead_rows = []
member_rows = []
window_bp = args.window_kb * 1000
path = Path(args.clumped)

if path.exists() and path.stat().st_size:
    lines = [line.split() for line in path.open(encoding="utf-8", errors="replace") if line.strip()]
    if lines:
        header = lines[0]
        for values in lines[1:]:
            row = dict(zip(header, values))
            lead_snp = row["SNP"]
            chromosome = row["CHR"]
            lead_position = int(float(row["BP"]))
            lead_p = row["P"]
            locus_id = f"{args.model}:{lead_snp}"

            members = [(lead_snp, "INDEX", chromosome, lead_position)]
            sp2 = row.get("SP2", "NONE")
            if sp2 not in {"NONE", "NA", "."}:
                for token in sp2.split(","):
                    member_snp = token.split("(", 1)[0]
                    member_chromosome, member_position = parse_variant_id(member_snp)
                    members.append((member_snp, "MEMBER", member_chromosome, member_position))

            valid_positions = [position for _, _, _, position in members if position is not None]
            clump_start = min(valid_positions) if valid_positions else lead_position
            clump_end = max(valid_positions) if valid_positions else lead_position

            lead_rows.append({
                "model": args.model,
                "source_locus_id": locus_id,
                "chromosome": chromosome,
                "position": lead_position,
                "lead_snp": lead_snp,
                "lead_p": lead_p,
                "window_start": max(1, lead_position - window_bp),
                "window_end": lead_position + window_bp,
                "clump_start": clump_start,
                "clump_end": clump_end,
                "member_count": len(members),
            })

            for member_snp, member_type, member_chromosome, member_position in members:
                member_rows.append({
                    "model": args.model,
                    "source_locus_id": locus_id,
                    "lead_snp": lead_snp,
                    "member_snp": member_snp,
                    "member_type": member_type,
                    "chromosome": member_chromosome or chromosome,
                    "position": "NA" if member_position is None else member_position,
                })

outputs = [
    (
        args.leads,
        lead_rows,
        [
            "model", "source_locus_id", "chromosome", "position", "lead_snp", "lead_p",
            "window_start", "window_end", "clump_start", "clump_end", "member_count",
        ],
    ),
    (
        args.members,
        member_rows,
        [
            "model", "source_locus_id", "lead_snp", "member_snp", "member_type",
            "chromosome", "position",
        ],
    ),
]

for output_path, rows, fields in outputs:
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
