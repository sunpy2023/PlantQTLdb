"""apply feature-level Bonferroni correction and annotate SMR HEIDI results.

requires pandas. Each input contains one dataset, trait and molecular type.
features with observed p_SMR define n_tested; threshold = alpha / n_tested.
HEIDI_pass means p_HEIDI >= the specified threshold; missing HEIDI remains NA.
all Bonferroni-significant SMR associations are retained with their HEIDI results.
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
    for name, value in [("Alpha", args.alpha), ("HEIDI threshold", args.heidi_threshold)]:
        if not math.isfinite(value) or not 0 < value <= 1:
            raise ValueError(f"{name} must be greater than 0 and at most 1")
    results = read_table(args.input).apply(lambda column: column.str.strip())
    if not {"probeID", "p_SMR", "p_HEIDI"}.issubset(results.columns):
        raise ValueError("SMR results require probeID, p_SMR and p_HEIDI columns")
    if results[["probeID", "p_SMR", "p_HEIDI"]].eq("").any().any():
        raise ValueError("Incomplete SMR rows; use NA for missing results")
    validate_ids(results["probeID"], "SMR feature IDs")
    for column in ["p_SMR", "p_HEIDI"]:
        results[column] = numeric_values(results[column], MISSING)
        values = results[column]
        if (~(values.isna() | values.between(0, 1))).any():
            raise ValueError(f"{column} must be between 0 and 1, or missing")
    n_tested = int(results["p_SMR"].notna().sum())
    if not n_tested:
        raise ValueError("No valid SMR P-values; no summary outputs were written")
    threshold = args.alpha / n_tested
    results["SMR_Bonferroni_threshold"] = threshold
    results["SMR_significant"] = results["p_SMR"].le(threshold)
    results["HEIDI_pass"] = results["p_HEIDI"].ge(args.heidi_threshold).astype("boolean")
    results.loc[results["p_HEIDI"].isna(), "HEIDI_pass"] = pd.NA
    results = results.sort_values("p_SMR", kind="stable", na_position="last")
    outputs = [Path(f"{args.out_prefix}.{suffix}.tsv") for suffix in ["all", "significant"]]
    check_outputs(outputs, [args.input])
    outputs[0].parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(outputs[0], sep="\t", index=False, na_rep="NA")
    significant = results.loc[results["SMR_significant"]]
    significant.to_csv(outputs[1], sep="\t", index=False, na_rep="NA")
    print(f"{n_tested} tested features; threshold={threshold:.17g}; "
          f"{len(significant)} significant SMR associations; "
          f"{int(significant['p_HEIDI'].isna().sum())} with missing HEIDI")
    for output in outputs:
        print(f"Saved {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="SMR .smr result file")
    parser.add_argument("--alpha", type=float, default=0.05, help="Family-wise significance level")
    parser.add_argument("--heidi-threshold", type=float, default=0.05,
                        help="HEIDI heterogeneity P-value cutoff")
    parser.add_argument("--out-prefix", required=True, help="Output path prefix")
    try:
        summarize_results(parser.parse_args())
    except (OSError, ValueError, pd.errors.ParserError) as error:
        parser.exit(1, f"Error: {error}\n")
