"""prepare FUSION inputs from normalized molecular phenotypes and full GEMMA results.

requires pandas. Matrix columns are feature ID followed by sample IDs; locations
use geneid, chr, s1 and s2 from prepare_qtl_phenotypes.py. Samples are matched by
IID in FAM order. Missing molecular values are omitted separately for each feature;
at least five observed samples and two distinct values are required for fivefold CV.
stable file keys preserve original feature IDs without putting them in shell paths.
outputs: features.tsv (no header), pheno/<key>.pheno (FID, IID, value; no header),
and <trait>.sumstats.tsv with allele1 as A1 and Z = beta / se.
"""

import argparse
import hashlib
from pathlib import Path
import re

import pandas as pd

from prepare_ewas_inputs import numeric_values, validate_ids
from prepare_qtl_phenotypes import MISSING, align_locations, read_table


def check_outputs(outputs, inputs):
    def identity(path):
        path = Path(path)
        if path.exists():
            info = path.stat()
            return info.st_dev, info.st_ino
        return "path", str(path.resolve())

    keys = [identity(path) for path in outputs]
    if len(keys) != len(set(keys)) or set(keys) & {identity(path) for path in inputs}:
        raise ValueError("Outputs must be distinct and must not overwrite input files")


def gwas_statistics(path):
    columns = ["rs", "allele1", "allele0", "beta", "se", "p_wald", "chr", "ps"]
    table = read_table(path, columns).apply(lambda column: column.str.strip())
    validate_ids(table["rs"], "GWAS variant IDs")
    for column in ["beta", "se", "p_wald"]:
        table[column] = numeric_values(table[column], MISSING)
    p_values = table["p_wald"]
    if (~(p_values.isna() | p_values.between(0, 1))).any():
        raise ValueError("GWAS P-values must be between 0 and 1, or missing")
    valid = table[["beta", "se", "p_wald"]].notna().all(axis=1) & table["se"].gt(0)
    valid &= ~table["allele1"].isin(MISSING) & ~table["allele0"].isin(MISSING)
    table = table.loc[valid].copy()
    if table.empty:
        raise ValueError("No GWAS records with observed beta, positive SE and valid P-values")
    for column in ["chr", "ps"]:
        values = pd.to_numeric(table[column], errors="raise")
        if not (values.gt(0) & values.lt(float("inf")) & values.mod(1).eq(0)).all():
            raise ValueError(f"GWAS {column} must contain positive integers")
        table[column] = values.astype("int64")
    z_scores = table["beta"] / table["se"]
    if z_scores.abs().eq(float("inf")).any():
        raise ValueError("GWAS beta / se produced a nonfinite Z-score")
    a1, a2 = table["allele1"].str.upper(), table["allele0"].str.upper()
    if (a1.eq(a2) | a1.str.contains(r"\s") | a2.str.contains(r"\s")).any():
        raise ValueError("GWAS alleles must be distinct and contain no whitespace")
    return pd.DataFrame({"SNP": table["rs"], "A1": a1, "A2": a2, "Z": z_scores,
                         "BETA": table["beta"], "SE": table["se"], "P": table["p_wald"],
                         "CHR": table["chr"], "BP": table["ps"]})


def prepare_inputs(args):
    out_dir = Path(args.out_dir).resolve()
    if not re.fullmatch(r"[A-Za-z0-9_/.-]+", str(out_dir)):
        raise ValueError("Use a shell-safe output path without spaces for FUSION's internal commands")
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", args.trait):
        raise ValueError("Trait name must contain only letters, numbers, underscores, dots or hyphens")

    rows = [line.split() for line in Path(args.fam).read_text().splitlines()]
    if not rows or any(len(row) != 6 for row in rows):
        raise ValueError("FAM requires six fields per row, without a header")
    samples = pd.DataFrame([row[:2] for row in rows], columns=["FID", "IID"])
    validate_ids(samples["IID"], "FAM individual IDs")
    if samples["FID"].isin(MISSING | {"None"}).any():
        raise ValueError("FAM family IDs must be nonempty")

    matrix = read_table(args.matrix)
    if matrix.eq("").any().any():
        raise ValueError("Empty molecular fields or incomplete rows; use NA for missing values")
    validate_ids(matrix.iloc[:, 0], "Molecular feature IDs")
    validate_ids(pd.Series(matrix.columns[1:], dtype=str), "Molecular sample IDs")
    matrix = matrix.set_index(matrix.columns[0])
    samples = samples.loc[samples["IID"].isin(matrix.columns)].set_index("IID", drop=False)
    if len(samples) < 5:
        raise ValueError("At least five samples must match between the matrix and FAM")
    values = matrix[samples["IID"]].apply(numeric_values, missing=MISSING | {"None"})
    keep = values.count(axis=1).ge(5) & values.nunique(axis=1, dropna=True).ge(2)
    values = values.loc[keep]
    if values.empty:
        raise ValueError("No variable molecular features have at least five observed samples")

    locations = read_table(args.locations, ["geneid", "chr", "s1", "s2"])
    locations = align_locations(locations, values.index, None, True)
    if len(locations) != len(values):
        raise ValueError("Every retained feature requires valid, matching coordinates")
    chromosomes = locations["chr"]
    if not chromosomes.str.fullmatch(r"[1-9][0-9]*").all():
        raise ValueError("Location chromosomes must use the same positive integer codes as genotype QC")
    keys = ["f_" + hashlib.sha256(name.encode()).hexdigest()[:24] for name in values.index]
    if len(keys) != len(set(keys)):
        raise ValueError("Feature file keys are not unique")
    features = pd.DataFrame({"key": keys, "id": values.index, "chr": chromosomes.to_numpy(),
                             "start": locations["s1"].to_numpy(), "end": locations["s2"].to_numpy()})
    statistics = gwas_statistics(args.gwas)
    outputs = [out_dir / "features.tsv", out_dir / f"{args.trait}.sumstats.tsv"]
    outputs += [out_dir / "pheno" / f"{key}.pheno" for key in keys]
    check_outputs(outputs, [args.matrix, args.locations, args.fam, args.gwas])

    (out_dir / "pheno").mkdir(parents=True, exist_ok=True)
    for key, (_, row) in zip(keys, values.iterrows()):
        observed = row.dropna()
        phenotype = samples.loc[observed.index, ["FID", "IID"]].copy()
        phenotype["value"] = observed.to_numpy()
        phenotype.to_csv(out_dir / "pheno" / f"{key}.pheno", sep="\t", index=False, header=False)
    features.to_csv(outputs[0], sep="\t", index=False, header=False)
    statistics.to_csv(outputs[1], sep="\t", index=False)
    print(f"{args.trait}: {len(samples)} matched samples; {len(values)}/{len(matrix)} features retained; "
          f"{len(statistics)} valid GWAS records")
    print(f"Saved training inputs and GWAS summary statistics in {out_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True, help="Normalized molecular CSV/TSV")
    parser.add_argument("--locations", required=True, help="TSV with geneid, chr, s1 and s2")
    parser.add_argument("--fam", required=True, help="FAM defining sample IDs and output order")
    parser.add_argument("--gwas", required=True, help="Full GEMMA association table")
    parser.add_argument("--trait", required=True, help="Agronomic trait name")
    parser.add_argument("--out-dir", required=True, help="FUSION training input/output directory")
    try:
        prepare_inputs(parser.parse_args())
    except (OSError, ValueError, pd.errors.ParserError) as error:
        parser.exit(1, f"Error: {error}\n")
