import re
from collections import defaultdict
from pathlib import Path
from gff_utils import read_fai, merge_intervals

ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')


def convert(gtf, fai, output):
    chrom_lengths = read_fai(fai)
    tx = defaultdict(list)
    with open(gtf, encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2].lower() != "exon":
                continue
            chrom, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
            attrs = dict(ATTR_RE.findall(parts[8]))
            transcript = attrs.get("transcript_id")
            if not transcript:
                raise ValueError(f"Missing transcript_id in GTF line {line_no}")
            if chrom not in chrom_lengths or end > chrom_lengths[chrom]:
                raise ValueError(f"GTF exon outside FASTA index at line {line_no}")
            tx[(chrom, strand, transcript)].append((start, end))
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    with open(output, "w", encoding="utf-8") as out:
        for (chrom, strand, transcript), intervals in sorted(tx.items()):
            merged = merge_intervals(intervals)
            chrom_start = merged[0][0] - 1
            chrom_end = merged[-1][1]
            block_sizes = [end - start + 1 for start, end in merged]
            block_starts = [start - 1 - chrom_start for start, _ in merged]
            fields = [
                chrom, str(chrom_start), str(chrom_end), transcript, "0", strand,
                str(chrom_start), str(chrom_end), "0", str(len(merged)),
                ",".join(map(str, block_sizes)) + ",",
                ",".join(map(str, block_starts)) + ",",
            ]
            out.write("\t".join(fields) + "\n")
            rows += 1
    if rows == 0:
        raise ValueError("No BED12 records generated")


def main_snakemake(sm):
    convert(str(sm.input.gtf), str(sm.input.fai), str(sm.output[0]))


if "snakemake" in globals():
    main_snakemake(snakemake)
