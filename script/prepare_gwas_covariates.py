"""write GEMMA PC covariates in the order of a PLINK FAM file.

usage: python3 prepare_gwas_covariates.py genotype.fam pca.eigenvec trait.txt output.txt pc_count
the PLINK 1.9 eigenvec file has no header: FID IID PC1 PC2 ...
the phenotype file has one value per FAM row, no header; missing values are NA.
phenotype row counts and numeric values are checked; its sample order must
already match FAM because the phenotype file contains no sample identifiers.
output has PC columns only, with no header or sample identifiers.
"""

import math
from pathlib import Path
import sys


def prepare_covariates(fam_path, eigenvec_path, phenotype_path, output_path, pc_count):
    if pc_count < 1:
        raise ValueError("PC count must be positive")
    with open(fam_path) as handle:
        fam = [line.split() for line in handle]
    if not fam or any(len(row) != 6 for row in fam):
        raise ValueError("Expected a nonempty six-column FAM file")
    samples = [tuple(row[:2]) for row in fam]
    if len(samples) != len(set(samples)):
        raise ValueError("Duplicate FID/IID pairs in FAM")

    pcs = {}
    with open(eigenvec_path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) != pc_count + 2:
                raise ValueError(f"Expected FID IID followed by {pc_count} PC values")
            sample = tuple(fields[:2])
            if sample in pcs:
                raise ValueError(f"Duplicate FID/IID pair in eigenvec: {sample}")
            if not all(math.isfinite(float(value)) for value in fields[2:]):
                raise ValueError(f"Non-finite PC values for {sample}")
            pcs[sample] = fields[2:]
    if set(pcs) != set(samples):
        raise ValueError("FAM and eigenvec sample sets differ")

    with open(phenotype_path) as handle:
        phenotype = [line.split() for line in handle]
    if len(phenotype) != len(samples) or any(len(row) != 1 for row in phenotype):
        raise ValueError("Expected one phenotype value per FAM row, without a header")
    observed = []
    for row in phenotype:
        if row[0] != "NA":
            value = float(row[0])
            if not math.isfinite(value):
                raise ValueError("Phenotype values must be finite or NA")
            observed.append(value)
    if len(set(observed)) < 2:
        raise ValueError("Phenotype must have at least two distinct observed values")

    output = Path(output_path)
    for source in (fam_path, eigenvec_path, phenotype_path):
        if output.exists() and output.samefile(source):
            raise ValueError("Covariate output must not overwrite an input file")
    with output.open("w") as handle:
        for sample in samples:
            handle.write("\t".join(pcs[sample]) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 6:
        raise SystemExit(__doc__)
    prepare_covariates(*sys.argv[1:5], pc_count=int(sys.argv[5]))
