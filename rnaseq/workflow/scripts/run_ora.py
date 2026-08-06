import csv
import gzip
import math
import re
from collections import defaultdict
from pathlib import Path
from scipy.stats import hypergeom
from annotation_utils import read_eggnog, split_terms


def bh_adjust(pvalues):
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    adjusted = [1.0] * n
    running = 1.0
    for rank_from_end, i in enumerate(reversed(order), 1):
        rank = n - rank_from_end + 1
        running = min(running, pvalues[i] * n / rank)
        adjusted[i] = min(1.0, running)
    return adjusted


def normalize_terms(column, value):
    terms = split_terms(value)
    out = []
    for raw in terms:
        term = raw
        if column == "GOs":
            if not re.fullmatch(r"GO:\d+", term):
                continue
        elif column == "KEGG_ko":
            term = re.sub(r"^ko:", "", term)
            if not re.fullmatch(r"K\d{5}", term):
                continue
        else:
            term = re.sub(r"^path:", "", term)
            if re.fullmatch(r"(?:map|ko)\d{5}", term):
                term = "map" + term[-5:]
        if term:
            out.append((term, raw))
    return sorted(set(out))


def build_term_maps(annotations):
    definitions = [("GO", "GOs"), ("KEGG_ko", "KEGG_ko"), ("KEGG_Pathway", "KEGG_Pathway")]
    maps = {}
    raw_sources = {}
    for term_type, column in definitions:
        gene_to_terms = defaultdict(set)
        source_ids = defaultdict(set)
        for gene, ann in annotations.items():
            for term, raw in normalize_terms(column, ann.get(column, "")):
                gene_to_terms[gene].add(term)
                source_ids[term].add(raw)
        maps[term_type] = gene_to_terms
        raw_sources[term_type] = source_ids
    return maps, raw_sources


def read_gene_ids(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return [row["gene_id"] for row in csv.DictReader(handle, delimiter="\t")]


def write_gzip_tsv(path, fields, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows(rows)


def run(annotation_path, complete_paths, significant_paths, outputs, padj_cutoff):
    annotations = read_eggnog(annotation_path)
    maps, raw_sources = build_term_maps(annotations)
    complete_index = {(Path(p).parent.name, Path(p).name.replace(".complete.tsv", "")): p for p in complete_paths}
    significant_index = {(Path(p).parent.name, Path(p).name.replace(".significant.tsv", "")): p for p in significant_paths}
    all_results = []
    summary = []
    term2gene_outputs = {"GO": outputs[2], "KEGG_ko": outputs[3], "KEGG_Pathway": outputs[4]}
    for term_type, outpath in term2gene_outputs.items():
        rows = []
        for gene, terms in maps[term_type].items():
            for term in sorted(terms):
                rows.append({"term_type": term_type, "term_id": term, "raw_source_term_ids": ";".join(sorted(raw_sources[term_type][term])), "gene_id": gene})
        write_gzip_tsv(outpath, ["term_type", "term_id", "raw_source_term_ids", "gene_id"], rows)
    for key in sorted(complete_index):
        tissue, contrast = key
        all_genes = set(read_gene_ids(complete_index[key]))
        deg_genes = set(read_gene_ids(significant_index[key]))
        if not deg_genes.issubset(all_genes):
            raise ValueError(f"DEG set is not a subset of tested genes: {tissue}/{contrast}")
        for term_type in ["GO", "KEGG_ko", "KEGG_Pathway"]:
            gene_to_terms = maps[term_type]
            background = sorted(all_genes & set(gene_to_terms))
            deg_annotated = sorted(deg_genes & set(background))
            N = len(background); n = len(deg_annotated)
            term_bg = defaultdict(set); term_deg = defaultdict(set)
            for gene in background:
                for term in gene_to_terms[gene]: term_bg[term].add(gene)
            for gene in deg_annotated:
                for term in gene_to_terms[gene]: term_deg[term].add(gene)
            tested_terms = sorted(term_deg)  # preserves the original Stage 9 behavior: k > 0 terms only
            pvalues = []
            temp = []
            for term in tested_terms:
                M = len(term_bg[term]); k = len(term_deg[term])
                pvalue = float(hypergeom.sf(k - 1, N, M, n)) if N and n else 1.0
                pvalues.append(pvalue)
                temp.append((term, M, k, pvalue))
            adjusted = bh_adjust(pvalues) if pvalues else []
            significant_terms = 0
            for (term, M, k, pvalue), padj in zip(temp, adjusted):
                significant = padj < padj_cutoff
                significant_terms += int(significant)
                gene_ratio = k / n if n else 0.0
                background_ratio = M / N if N else 0.0
                all_results.append({
                    "tissue": tissue, "contrast": contrast, "term_type": term_type,
                    "term_id": term, "raw_source_term_ids": ";".join(sorted(raw_sources[term_type][term])),
                    "DEG_n_total": len(deg_genes), "DEG_n_annotated": n,
                    "background_n_total": len(all_genes), "background_n_annotated": N,
                    "term_gene_n_in_background": M, "term_DEG_n": k,
                    "gene_ratio": gene_ratio, "background_ratio": background_ratio,
                    "fold_enrichment": gene_ratio / background_ratio if background_ratio else 0.0,
                    "pvalue": pvalue, "padj": padj, "significant": "yes" if significant else "no",
                    "DEG_gene_ids": ";".join(sorted(term_deg[term])),
                    "background_term_gene_ids": ";".join(sorted(term_bg[term])),
                })
            summary.append({
                "tissue": tissue, "contrast": contrast, "term_type": term_type,
                "background_n_total": len(all_genes), "background_n_annotated": N,
                "background_annotation_coverage": N / len(all_genes) if all_genes else 0.0,
                "DEG_n_total": len(deg_genes), "DEG_n_annotated": n,
                "DEG_annotation_coverage": n / len(deg_genes) if deg_genes else 0.0,
                "terms_tested_with_DEG_overlap": len(tested_terms),
                "significant_terms_padj_lt_cutoff": significant_terms,
            })
    result_fields = ["tissue", "contrast", "term_type", "term_id", "raw_source_term_ids", "DEG_n_total", "DEG_n_annotated", "background_n_total", "background_n_annotated", "term_gene_n_in_background", "term_DEG_n", "gene_ratio", "background_ratio", "fold_enrichment", "pvalue", "padj", "significant", "DEG_gene_ids", "background_term_gene_ids"]
    write_gzip_tsv(outputs[0], result_fields, all_results)
    Path(outputs[1]).parent.mkdir(parents=True, exist_ok=True)
    with open(outputs[1], "w", encoding="utf-8", newline="") as out:
        fields = list(summary[0]) if summary else ["tissue", "contrast", "term_type"]
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows(summary)


def main_snakemake(sm):
    run(str(sm.input.annotation), list(map(str, sm.input.complete)), list(map(str, sm.input.significant)), [str(sm.output.results), str(sm.output.summary), str(sm.output.go), str(sm.output.ko), str(sm.output.pathway)], float(sm.params.padj))


if "snakemake" in globals():
    main_snakemake(snakemake)
