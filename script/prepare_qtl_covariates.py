"""merge PEER factors with genotype PCs for molecular QTL mapping.

requires pandas and numpy. Reuses the phenotype workflow's sample-list reader.
reads a PEER matrix and a headerless PLINK 1.9 eigenvec file, retaining
sample ID strings and writing PEER rows followed by PC rows in sample-list order.
both inputs must contain exactly the requested samples; missing, ambiguous or
non-finite values are rejected. The merged covariates must leave positive
residual degrees of freedom and have full rank together with an intercept.
the intercept is checked but is not written to the output matrix.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from prepare_qtl_phenotypes import MISSING, read_samples, read_table


def read_eigenvec(path):
    table = pd.read_csv(path, sep=r"\s+", header=None, dtype=str, keep_default_na=False)
    if table.empty or table.shape[1] < 3:
        raise ValueError("Expected FID IID followed by PC columns in eigenvec")
    if (table.iloc[:, :2].isna().any().any() or table.iloc[:, :2].eq("").any().any()
            or table.iloc[:, 1].duplicated().any()):
        raise ValueError(f"FID/IID fields must be nonempty and each IID must be unique: {path}")
    return table.set_index(1, drop=False)


def numeric_matrix(table):
    values = table.apply(pd.to_numeric, errors="raise").astype(float)
    if not np.isfinite(values.to_numpy()).all():
        raise ValueError("Covariates must contain finite numbers, without missing values")
    return values


def check_samples(actual, expected, source):
    if set(actual) != set(expected):
        missing = sorted(set(expected) - set(actual))
        extra = sorted(set(actual) - set(expected))
        raise ValueError(f"Sample IDs differ in {source}; missing: {missing}; extra: {extra}")


def merge_covariates(peer_path, eigenvec_path, samples, pc_count):
    if pc_count < 1:
        raise ValueError("PC count must be positive")
    peer = read_table(peer_path)
    if peer.columns[0] != "id" or peer.empty:
        raise ValueError("PEER matrix must start with id and contain factor rows")
    ids = peer["id"]
    if ids.isin(MISSING).any() or ids.duplicated().any() or ids.str.contains(r"\s", regex=True).any():
        raise ValueError("PEER factor IDs must be nonempty, distinct and contain no whitespace")
    check_samples(peer.columns[1:], samples, "PEER matrix")
    peer = numeric_matrix(peer.set_index("id")[samples])
    eigenvec = read_eigenvec(eigenvec_path)
    check_samples(eigenvec.index, samples, "PCA eigenvec")
    if eigenvec.shape[1] < pc_count + 2:
        raise ValueError(f"PCA eigenvec contains fewer than {pc_count} PCs")
    pcs = numeric_matrix(eigenvec.loc[samples].iloc[:, 2:pc_count + 2]).T
    pcs.index = [f"PC{i}" for i in range(1, pc_count + 1)]
    merged = pd.concat([peer, pcs])
    if merged.index.duplicated().any():
        raise ValueError("PEER and PC row IDs overlap")
    if len(samples) <= len(merged) + 2:
        raise ValueError("Insufficient samples for covariates, an intercept and a tested variant")
    values = merged.to_numpy().T
    magnitudes = np.abs(values).max(axis=0)
    centered = values / np.where(magnitudes == 0, 1, magnitudes)
    centered -= centered.mean(axis=0)
    scales = np.linalg.norm(centered, axis=0)
    if np.any(scales == 0):
        raise ValueError("Constant covariates cannot be fitted with an intercept")
    design = np.column_stack([np.ones(len(samples)), centered / scales])
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("PEER factors and genotype PCs are collinear with an intercept")
    merged.index.name = "id"
    return merged


def prepare_covariates(args):
    inputs = [args.samples, args.peer, args.eigenvec]
    output = Path(args.output)
    for source in map(Path, inputs):
        if (output.resolve() == source.resolve()
                or (output.exists() and source.exists() and output.samefile(source))):
            raise ValueError("Output must not overwrite an input file")
    samples = read_samples(args.samples)
    table = merge_covariates(args.peer, args.eigenvec, samples, args.pc_count)
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, sep="\t", index_label="id")
    print(f"Saved {output}: {len(samples)} samples, {len(table)} covariates")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True,
                        help="Same sample list used for phenotype preparation")
    parser.add_argument("--output", required=True)
    parser.add_argument("--peer", required=True, help="PEER factor TSV: id followed by sample columns")
    parser.add_argument("--eigenvec", required=True,
                        help="Headerless PLINK 1.9 eigenvec: FID IID PC1 PC2 ...")
    parser.add_argument("--pc-count", type=int, default=5)
    try:
        prepare_covariates(parser.parse_args())
    except (OSError, ValueError, pd.errors.ParserError) as error:
        parser.exit(1, f"Error: {error}\n")
