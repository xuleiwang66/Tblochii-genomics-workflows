import urllib.parse
from collections import defaultdict

TRANSCRIPT_FEATURES = {
    "mrna", "transcript", "lnc_rna", "ncrna", "rrna", "trna",
    "snrna", "snorna", "mirna", "primary_transcript"
}


def parse_attributes(text):
    out = {}
    for item in text.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
        else:
            continue
        out[key.strip()] = urllib.parse.unquote(value.strip().strip('"'))
    return out


def split_parent(value):
    if not value:
        return []
    return [x.strip() for x in str(value).split(",") if x.strip()]


def read_chromosomes(path):
    with open(path, encoding="utf-8") as handle:
        return [x.strip() for x in handle if x.strip() and not x.startswith("#")]


def read_fai(path):
    lengths = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                chrom, length, *_ = line.rstrip("\n").split("\t")
                lengths[chrom] = int(length)
    return lengths


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(x) for x in merged]


def scan_gene_transcript_maps(gff3):
    genes = set()
    transcript_to_genes = defaultdict(set)
    with open(gff3, encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                raise ValueError(f"Malformed GFF3 line {line_no}: expected 9 columns")
            feature = parts[2].lower()
            attrs = parse_attributes(parts[8])
            if feature == "gene":
                gene_id = attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name")
                if gene_id:
                    genes.add(gene_id)
            elif feature in TRANSCRIPT_FEATURES:
                tx_id = attrs.get("ID") or attrs.get("transcript_id")
                parents = split_parent(attrs.get("Parent") or attrs.get("gene_id"))
                if tx_id:
                    transcript_to_genes[tx_id].update(parents)
    return genes, transcript_to_genes


def resolve_exon(parents, genes, transcript_to_genes):
    resolved = []
    for parent in parents:
        if parent in transcript_to_genes:
            for gene in sorted(transcript_to_genes[parent]):
                resolved.append((parent, gene))
        elif parent in genes:
            resolved.append((parent, parent))
    return sorted(set(resolved))
