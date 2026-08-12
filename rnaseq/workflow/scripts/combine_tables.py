import gzip
from pathlib import Path


def opener(path, mode):
    return gzip.open(path, mode + "t", encoding="utf-8") if str(path).endswith(".gz") else open(path, mode, encoding="utf-8")


def combine(inputs, output):
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    expected_header = None
    with opener(output, "w") as out:
        for path in inputs:
            with opener(path, "r") as handle:
                header = handle.readline()
                if not header:
                    raise ValueError(f"Empty table: {path}")
                if expected_header is None:
                    expected_header = header
                    out.write(header)
                elif header != expected_header:
                    raise ValueError(f"Header mismatch while combining {path}")
                for line in handle:
                    out.write(line)


def main_snakemake(sm):
    combine(list(map(str, sm.input)), str(sm.output[0]))


if "snakemake" in globals():
    main_snakemake(snakemake)
