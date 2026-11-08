"""merge FUSION chromosome results and apply feature-level Bonferroni correction.

requires pandas. Weights must have WGT, ID, CHR, P0 and P1 columns. Each indexed
chromosome requires <trait>.chr<chr>.dat; an accompanying .dat.MHC is merged too.
row counts must match indexed weights; duplicate or unexpected feature IDs fail.
unidentified placeholder rows from failed prediction models are omitted. Features
with missing TWAS.P remain in the full table but do not enter the tested count.
the threshold is alpha / n_tested for this dataset, trait and molecular type.
outputs: <prefix>.all.tsv and <prefix>.significant.tsv.
"""

import argparse
import math
from pathlib import Path

import pandas as pd

from prepare_ewas_inputs import numeric_values, validate_ids
from prepare_qtl_phenotypes import MISSING, read_table
from prepare_twas_inputs import check_outputs


def summarize_results(args):
    if not math.isfinite(args.alpha) or not 0 < args.alpha <= 1:
        raise ValueError("Alpha must be greater than 0 and at most 1")
    weights = read_table(args.weights)
    required = {"WGT", "ID", "CHR", "P0", "P1"}
    if not required.issubset(weights.columns) or weights.empty:
        raise ValueError("Weights require a nonempty table with WGT, ID, CHR, P0 and P1 columns")
    validate_ids(weights["ID"], "Indexed feature IDs")
    if not weights["CHR"].str.fullmatch(r"[1-9][0-9]*").all():
        raise ValueError("Weights must use positive integer chromosome codes")
    tables, inputs, columns = [], [args.weights], None
    failed = 0
    for chromosome, group in weights.groupby("CHR", sort=False):
        source = Path(args.results_dir) / f"{args.trait}.chr{chromosome}.dat"
        paths = [source]
        if Path(f"{source}.MHC").exists():
            paths.append(Path(f"{source}.MHC"))
        pieces = []
        for path in paths:
            table = read_table(path).apply(lambda column: column.str.strip())
            if not {"ID", "CHR", "TWAS.P"}.issubset(table.columns):
                raise ValueError(f"Missing ID, CHR or TWAS.P columns: {path}")
            if columns is None:
                columns = table.columns.tolist()
            if table.columns.tolist() != columns:
                raise ValueError(f"Inconsistent FUSION result columns: {path}")
            pieces.append(table)
            inputs.append(path)
        table = pd.concat(pieces, ignore_index=True)
        if len(table) != len(group):
            raise ValueError(f"Result row count does not match indexed weights on chromosome {chromosome}")
        identified = ~table["ID"].isin(MISSING)
        failed += int((~identified).sum())
        table = table.loc[identified].copy()
        if (not table["ID"].isin(group["ID"]).all() or not table["CHR"].eq(chromosome).all()
                or table[["CHR", "TWAS.P"]].eq("").any().any()):
            raise ValueError(f"Invalid feature IDs, chromosomes or incomplete rows in {source}")
        tables.append(table)

    results = pd.concat(tables, ignore_index=True)
    validate_ids(results["ID"], "Result feature IDs")
    results["TWAS.P"] = numeric_values(results["TWAS.P"], MISSING)
    p_values = results["TWAS.P"]
    if (~(p_values.isna() | p_values.between(0, 1))).any():
        raise ValueError("TWAS P-values must be between 0 and 1, or missing")
    n_tested = int(p_values.notna().sum())
    if not n_tested:
        raise ValueError("No valid TWAS P-values; no summary outputs were written")
    threshold = args.alpha / n_tested
    results["TWAS_Bonferroni_threshold"] = threshold
    results["TWAS_significant"] = p_values.le(threshold)
    results = results.sort_values("TWAS.P", kind="stable", na_position="last")
    significant = results.loc[results["TWAS_significant"]]
    outputs = [Path(f"{args.out_prefix}.{suffix}.tsv") for suffix in ["all", "significant"]]
    check_outputs(outputs, inputs)
    outputs[0].parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(outputs[0], sep="\t", index=False, na_rep="NA")
    significant.to_csv(outputs[1], sep="\t", index=False, na_rep="NA")
    print(f"{args.trait}: {n_tested} tested features; {failed} unidentified failed-model rows omitted; "
          f"threshold={threshold:.17g}; {len(significant)} significant associations")
    for output in outputs:
        print(f"Saved {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, help="Directory containing chromosome .dat files")
    parser.add_argument("--weights", required=True, help="FUSION weight index for this molecular type")
    parser.add_argument("--trait", required=True, help="Agronomic trait name matching chromosome filenames")
    parser.add_argument("--alpha", type=float, default=0.05, help="Family-wise significance level")
    parser.add_argument("--out-prefix", required=True, help="Output path prefix")
    try:
        summarize_results(parser.parse_args())
    except (OSError, ValueError, pd.errors.ParserError) as error:
        parser.exit(1, f"Error: {error}\n")
