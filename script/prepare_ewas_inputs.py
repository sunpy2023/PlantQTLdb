"""prepare trait-specific OSCA inputs from filtered FPKM, PSI or PDUI matrices.

requires pandas. The molecular CSV/TSV has feature IDs in its first column and
sample IDs in the remaining columns. The one-column agronomic phenotype file
must follow FAM row order, without a header; NA and -9 denote missing traits.
samples are matched by IID and written in FAM order. Use NA for missing molecular
values; features with fewer than two distinct observed values are removed.
outputs: <prefix>.omics.tsv and <prefix>.phenotype.tsv (FID, IID, trait; no header).
"""

import argparse
from pathlib import Path

import pandas as pd

from prepare_qtl_phenotypes import MISSING, read_table


def validate_ids(ids, label):
    if (ids.isin(MISSING | {"None"}).any() or ids.duplicated().any()
            or ids.str.contains(r"\s", regex=True).any()):
        raise ValueError(f"{label} must be nonempty, distinct and contain no whitespace")


def numeric_values(values, missing):
    text = values.str.strip()
    numeric = pd.to_numeric(text.mask(text.isin(missing)), errors="raise")
    if numeric.isin([float("inf"), float("-inf")]).any():
        raise ValueError("Numeric values must be finite or missing")
    return numeric


def prepare_inputs(args):
    inputs = [Path(value) for value in [args.matrix, args.fam, args.phenotype]]
    outputs = [Path(f"{args.out_prefix}.{suffix}.tsv") for suffix in ["omics", "phenotype"]]
    for index, output in enumerate(outputs):
        for other in inputs + outputs[index + 1:]:
            if (output.resolve() == other.resolve()
                    or (output.exists() and other.exists() and output.samefile(other))):
                raise ValueError("Outputs must be distinct and must not overwrite input files")
    if not args.trait.strip() or any(char in args.trait for char in "\t\r\n"):
        raise ValueError("Trait name must be nonempty and contain no tabs or line breaks")

    rows = [line.split() for line in inputs[1].read_text().splitlines()]
    if not rows or any(len(row) != 6 for row in rows):
        raise ValueError("FAM requires six fields per row, without a header")
    samples = pd.DataFrame([row[:2] for row in rows], columns=["FID", "IID"])
    validate_ids(samples["IID"], "FAM individual IDs")
    if samples["FID"].isin(MISSING | {"None"}).any():
        raise ValueError("FAM family IDs must be nonempty")

    rows = [line.split() for line in inputs[2].read_text().splitlines()]
    if len(rows) != len(samples) or any(len(row) != 1 for row in rows):
        raise ValueError("Phenotype requires one value per FAM row, without a header; use NA for missing")
    phenotype = numeric_values(pd.Series([row[0] for row in rows]), MISSING | {"None"})
    samples["trait"] = phenotype.mask(phenotype.eq(-9))

    matrix = read_table(inputs[0])
    if matrix.eq("").any().any():
        raise ValueError("Empty molecular fields or incomplete rows; use NA for missing values")
    validate_ids(matrix.iloc[:, 0], "Molecular feature IDs")
    if matrix.iloc[:, 0].isin(["FID", "IID"]).any():
        raise ValueError("Molecular feature IDs must not use the reserved names FID or IID")
    validate_ids(pd.Series(matrix.columns[1:], dtype=str), "Molecular sample IDs")
    matrix = matrix.set_index(matrix.columns[0])
    samples = samples.loc[samples["IID"].isin(matrix.columns) & samples["trait"].notna()].copy()
    if len(samples) < 3 or samples["trait"].nunique() < 2:
        raise ValueError("At least three matched samples and a variable observed phenotype are required")
    values = matrix[samples["IID"]].apply(numeric_values, missing=MISSING | {"None"})
    n_features = len(values)
    values = values.loc[values.nunique(axis=1, dropna=True).ge(2)]
    if values.empty:
        raise ValueError("No variable molecular features remain in the matched samples")

    omics = values.T.reset_index(drop=True)
    omics.columns.name = None
    omics = pd.concat([samples[["FID", "IID"]].reset_index(drop=True), omics], axis=1)
    outputs[0].parent.mkdir(parents=True, exist_ok=True)
    omics.to_csv(outputs[0], sep="\t", index=False, na_rep="NA")
    samples.to_csv(outputs[1], sep="\t", index=False, header=False)
    print(f"{args.trait}: {len(phenotype)} FAM samples; {len(samples)} matched with observed trait; "
          f"{len(values)}/{n_features} molecular features retained")
    for output in outputs:
        print(f"Saved {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True, help="Filtered molecular CSV/TSV, before normalization")
    parser.add_argument("--fam", required=True, help="FAM defining sample IDs and phenotype row order")
    parser.add_argument("--phenotype", required=True, help="One trait column in FAM order, without a header")
    parser.add_argument("--trait", required=True, help="Agronomic trait name")
    parser.add_argument("--out-prefix", required=True, help="Output path prefix")
    try:
        prepare_inputs(parser.parse_args())
    except (OSError, ValueError, pd.errors.ParserError) as error:
        parser.exit(1, f"Error: {error}\n")
