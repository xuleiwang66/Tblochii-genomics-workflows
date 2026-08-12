import csv
from collections import defaultdict
from pathlib import Path
from gff_utils import parse_attributes, split_parent, read_chromosomes, read_fai, merge_intervals, scan_gene_transcript_maps, resolve_exon


def build(gff3, fai, chromosomes_file, saf_out, lengths_out):
    chromosomes = read_chromosomes(chromosomes_file)
    allowed = set(chromosomes)
    rank = {c: i for i, c in enumerate(chromosomes)}
    lengths = read_fai(fai)
    missing = [c for c in chromosomes if c not in lengths]
    if missing:
        raise ValueError(f"Configured chromosomes absent from FASTA index: {missing}")
    genes, tx_to_genes = scan_gene_transcript_maps(gff3)
    grouped = defaultdict(lambda: defaultdict(list))
    with open(gff3, encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                raise ValueError(f"Malformed GFF3 line {line_no}")
            seqid, _, feature, start, end, _, strand, _, attrs_text = parts
            if feature.lower() != "exon" or seqid not in allowed:
                continue
            s, e = int(start), int(end)
            if s < 1 or e < s or e > lengths[seqid]:
                raise ValueError(f"Invalid exon coordinates at GFF3 line {line_no}")
            parents = split_parent(parse_attributes(attrs_text).get("Parent"))
            if not parents:
                raise ValueError(f"Exon without Parent at GFF3 line {line_no}")
            resolved = resolve_exon(parents, genes, tx_to_genes)
            if not resolved:
                raise ValueError(f"Exon Parent cannot resolve to gene at GFF3 line {line_no}")
            for _, gene_id in resolved:
                grouped[gene_id][(seqid, strand)].append((s, e))
    saf_rows = []
    gene_rows = []
    for gene_id, by_location in grouped.items():
        total = 0
        intervals_n = 0
        chroms = []
        strands = []
        for (chrom, strand), intervals in by_location.items():
            chroms.append(chrom)
            strands.append(strand)
            for start, end in merge_intervals(intervals):
                saf_rows.append((gene_id, chrom, start, end, strand))
                total += end - start + 1
                intervals_n += 1
        gene_rows.append((gene_id, total, intervals_n, ",".join(sorted(set(chroms), key=lambda x: rank[x])), ",".join(sorted(set(strands)))))
    if not saf_rows:
        raise ValueError("No SAF rows generated")
    saf_rows.sort(key=lambda x: (rank[x[1]], x[2], x[3], x[0]))
    gene_rows.sort(key=lambda x: x[0])
    Path(saf_out).parent.mkdir(parents=True, exist_ok=True)
    with open(saf_out, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["GeneID", "Chr", "Start", "End", "Strand"])
        writer.writerows(saf_rows)
    with open(lengths_out, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["gene_id", "gene_length_bp", "n_union_intervals", "chromosomes", "strands"])
        writer.writerows(gene_rows)


def main_snakemake(sm):
    build(str(sm.input.gff3), str(sm.input.fai), str(sm.input.chromosomes), str(sm.output.saf), str(sm.output.lengths))


if "snakemake" in globals():
    main_snakemake(snakemake)
