import csv
import re

ANNOTATION_COLUMNS = ["Preferred_name", "Description", "GOs", "KEGG_ko", "KEGG_Pathway", "PFAMs"]


def read_eggnog(path):
    with open(path, encoding="utf-8") as handle:
        lines = [line for line in handle if not line.startswith("##")]
    header_index = None
    for i, line in enumerate(lines):
        if line.startswith("#query\t") or line.startswith("query\t"):
            header_index = i
            break
    if header_index is None:
        raise ValueError("eggNOG annotation header containing query was not found")
    lines[header_index] = lines[header_index].lstrip("#")
    reader = csv.DictReader(lines[header_index:], delimiter="\t")
    if "query" not in (reader.fieldnames or []):
        raise ValueError("eggNOG annotation lacks query column")
    annotations = {}
    for row in reader:
        query = (row.get("query") or "").strip()
        if not query or query in annotations:
            continue
        annotations[query] = {col: (row.get(col) or "").strip() for col in ANNOTATION_COLUMNS}
    return annotations


def split_terms(value):
    if not value or value in {"-", "NA", "."}:
        return []
    return sorted({x.strip() for x in re.split(r"[,;|]+", value) if x.strip() and x.strip() not in {"-", "NA", "."}})
