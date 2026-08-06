import csv
from pathlib import Path

REQUIRED = [
    "sample_id", "r1", "r2", "source_cohort", "treatment", "tissue",
    "biological_replicate", "replacement_pituitary"
]


def validate(input_path, output_path, source_levels, treatment_levels, tissue_levels):
    with open(input_path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [x for x in REQUIRED if x not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Missing samples.tsv columns: {missing}")
        rows = []
        seen = set()
        for line_no, row in enumerate(reader, 2):
            clean = {k: (row.get(k) or "").strip() for k in REQUIRED}
            sid = clean["sample_id"]
            if not sid:
                raise ValueError(f"Empty sample_id at line {line_no}")
            if sid in seen:
                raise ValueError(f"Duplicated sample_id: {sid}")
            seen.add(sid)
            if not clean["r1"] or not clean["r2"]:
                raise ValueError(f"Missing R1/R2 path for {sid}")
            if clean["source_cohort"] not in source_levels:
                raise ValueError(f"Invalid source_cohort for {sid}: {clean['source_cohort']}")
            if clean["treatment"] not in treatment_levels:
                raise ValueError(f"Invalid treatment for {sid}: {clean['treatment']}")
            if clean["tissue"] not in tissue_levels:
                raise ValueError(f"Invalid tissue for {sid}: {clean['tissue']}")
            if clean["replacement_pituitary"].lower() not in {"yes", "no"}:
                raise ValueError(f"replacement_pituitary must be yes/no for {sid}")
            clean["replacement_pituitary"] = clean["replacement_pituitary"].lower()
            rows.append(clean)
    if not rows:
        raise ValueError("samples.tsv contains no samples")
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=REQUIRED, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main_snakemake(sm):
    design = sm.config["experimental_design"]
    validate(
        str(sm.input[0]), str(sm.output[0]),
        set(design["source_levels"]), set(design["treatment_levels"]),
        set(design["tissue_levels"]),
    )


if "snakemake" in globals():
    main_snakemake(snakemake)
