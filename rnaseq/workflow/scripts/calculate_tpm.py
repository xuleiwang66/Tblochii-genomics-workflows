import csv
from pathlib import Path


def calculate(counts_path, lengths_path, matrix_out, sums_out):
    lengths = {}
    with open(lengths_path, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = row["gene_id"]
            length = float(row["gene_length_bp"])
            if length <= 0 or gene in lengths:
                raise ValueError(f"Invalid or duplicated gene length: {gene}")
            lengths[gene] = length
    with open(counts_path, encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        samples = header[1:]
        rows = []
        denominators = [0.0] * len(samples)
        seen = set()
        for row in reader:
            gene = row[0]
            if gene in seen or gene not in lengths:
                raise ValueError(f"Count/gene-length mismatch or duplicate: {gene}")
            seen.add(gene)
            values = [float(x) for x in row[1:]]
            rpk = [value / (lengths[gene] / 1000.0) for value in values]
            denominators = [a + b for a, b in zip(denominators, rpk)]
            rows.append((gene, rpk))
    if seen != set(lengths):
        raise ValueError("Count matrix and gene_length.tsv gene sets differ")
    if any(x <= 0 for x in denominators):
        raise ValueError("TPM denominator is zero for at least one sample")
    Path(matrix_out).parent.mkdir(parents=True, exist_ok=True)
    sums = [0.0] * len(samples)
    with open(matrix_out, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t"); writer.writerow(["gene_id"] + samples)
        for gene, rpk in rows:
            tpm = [value / denom * 1_000_000.0 for value, denom in zip(rpk, denominators)]
            sums = [a + b for a, b in zip(sums, tpm)]
            writer.writerow([gene] + [f"{x:.6f}" for x in tpm])
    with open(sums_out, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t"); writer.writerow(["sample_id", "tpm_sum", "absolute_error_from_1e6"])
        for sid, total in zip(samples, sums):
            writer.writerow([sid, f"{total:.6f}", f"{abs(total - 1_000_000.0):.6f}"])


def main_snakemake(sm):
    calculate(str(sm.input.counts), str(sm.input.lengths), str(sm.output.matrix), str(sm.output.sums))


if "snakemake" in globals():
    main_snakemake(snakemake)
