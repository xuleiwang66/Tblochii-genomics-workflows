from pathlib import Path
from gff_utils import parse_attributes, split_parent, read_chromosomes, scan_gene_transcript_maps, resolve_exon


def convert(gff3, chromosomes_file, output):
    allowed = set(read_chromosomes(chromosomes_file))
    genes, tx_to_genes = scan_gene_transcript_maps(gff3)
    emitted = 0
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(gff3, encoding="utf-8") as handle, open(output, "w", encoding="utf-8") as out:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                raise ValueError(f"Malformed GFF3 line {line_no}")
            seqid, source, feature, start, end, score, strand, phase, attrs_text = parts
            if feature.lower() != "exon" or seqid not in allowed:
                continue
            attrs = parse_attributes(attrs_text)
            parents = split_parent(attrs.get("Parent"))
            if not parents:
                raise ValueError(f"Exon without Parent at GFF3 line {line_no}")
            resolved = resolve_exon(parents, genes, tx_to_genes)
            if not resolved:
                raise ValueError(f"Exon Parent cannot resolve to gene at GFF3 line {line_no}")
            s, e = int(start), int(end)
            if s < 1 or e < s:
                raise ValueError(f"Invalid exon coordinates at GFF3 line {line_no}")
            for transcript_id, gene_id in resolved:
                attributes = f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
                out.write("\t".join([seqid, source, "exon", str(s), str(e), score, strand, phase, attributes]) + "\n")
                emitted += 1
    if emitted == 0:
        raise ValueError("No exon records were written to STAR GTF")


def main_snakemake(sm):
    convert(str(sm.input.gff3), str(sm.input.chromosomes), str(sm.output[0]))


if "snakemake" in globals():
    main_snakemake(snakemake)
