"""merge StringTie gene abundance tables into Coverage, FPKM and TPM matrices.

usage: python3 merge_gene_expression.py sample_gtf.tsv stringtie_dir output_dir
the sample/GTF list has two tab-separated columns without a header.
columns follow sample-list order; genes absent from a sample table receive 0.
use complete StringTie abundance tables from the same reference annotation.
"""

import csv
import math
from pathlib import Path
import sys


def merge_expression(sample_list, quant_dir, output_dir):
    samples = []
    with open(sample_list) as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 2 or not fields[0] or not fields[1]:
                raise ValueError("Expected sample ID<TAB>GTF path, without a header")
            samples.append(fields[0])
    if not samples or len(samples) != len(set(samples)):
        raise ValueError("Sample list must be nonempty and contain unique IDs")
    if any(len(sample.split()) != 1 or "/" in sample or sample in {".", ".."}
           for sample in samples):
        raise ValueError("Sample IDs must not contain whitespace or path components")

    metrics = ("Coverage", "FPKM", "TPM")
    data = {}
    for sample in samples:
        path = Path(quant_dir) / sample / f"{sample}.gene_abund.tsv"
        seen = set()
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"Gene ID", *metrics}.issubset(reader.fieldnames or []):
                raise ValueError(f"Missing Gene ID/Coverage/FPKM/TPM columns: {path}")
            for row in reader:
                gene = (row["Gene ID"] or "").strip()
                if not gene or gene in seen:
                    raise ValueError(f"Empty or duplicate Gene ID {gene!r}: {path}")
                seen.add(gene)
                values = tuple(row[metric] for metric in metrics)
                if any(value is None or not value.strip() for value in values):
                    raise ValueError(f"Missing abundance value for {gene}: {path}")
                try:
                    valid = all(math.isfinite(float(value)) and float(value) >= 0
                                for value in values)
                except ValueError:
                    valid = False
                if not valid:
                    raise ValueError(f"Expected finite, nonnegative abundances for {gene}: {path}")
                data.setdefault(gene, {})[sample] = values

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for index, metric in enumerate(metrics):
        path = output_dir / f"merged_gene_{metric}_matrix.csv"
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["gene_id", *samples])
            for gene in sorted(data):
                writer.writerow([gene] + [data[gene].get(sample, ("0",) * 3)[index]
                                          for sample in samples])


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    merge_expression(*sys.argv[1:])
