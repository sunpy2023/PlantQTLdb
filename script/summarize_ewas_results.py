"""label OSCA MOMENT or converged MOA fallback results and filter unadjusted P values.

requires pandas. The model is inferred from the matching OSCA log, since both
models can produce a .moment file. The final model must converge and finish
without errors; its reported probe count must match the result table.
outputs: <prefix>.all.tsv and <prefix>.significant.tsv, including MODEL.
missing P values remain in the full table and never enter the significant table.
"""

import argparse
import math
from pathlib import Path
import re

import pandas as pd

from prepare_ewas_inputs import numeric_values, validate_ids
from prepare_qtl_phenotypes import MISSING, read_table


def summarize_results(args):
    if not math.isfinite(args.p_threshold) or not 0 < args.p_threshold <= 1:
        raise ValueError("P-value threshold must be greater than 0 and at most 1")
    inputs = [Path(args.input), Path(args.log)]
    outputs = [Path(f"{args.out_prefix}.{suffix}.tsv") for suffix in ["all", "significant"]]
    for index, output in enumerate(outputs):
        for other in inputs + outputs[index + 1:]:
            if (output.resolve() == other.resolve()
                    or (output.exists() and other.exists() and output.samefile(other))):
                raise ValueError("Outputs must be distinct and must not overwrite input files")

    log = inputs[1].read_text()
    model = "MOA" if "Degrade to MOA" in log else "MOMENT"
    final_log = log.rsplit("Degrade to MOA", 1)[-1]
    if ("Error:" in final_log or "Log-likelihood ratio converged." not in final_log
            or not re.search(r"^Analysis finished:", final_log, re.MULTILINE)):
        raise ValueError(f"{model} did not converge and finish successfully; inspect {inputs[1]}")
    counts = re.findall(r"Results of (\d+) probes have been saved in file", final_log)
    if not counts:
        raise ValueError("Log does not confirm that association results were saved")

    results = read_table(inputs[0])
    columns = ["Chr", "Probe", "bp", "Gene", "Orientation", "b", "se", "p"]
    if results.columns.tolist() != columns or results.eq("").any().any():
        raise ValueError("Expected complete OSCA rows with Chr, Probe, bp, Gene, Orientation, b, se and p")
    if len(results) != int(counts[-1]):
        raise ValueError("Result row count does not match the final model's reported probe count")
    validate_ids(results["Probe"], "OSCA probe IDs")
    p_values = numeric_values(results["p"], MISSING)
    p_values = p_values.mask(p_values.eq(-9))
    if (~(p_values.isna() | p_values.between(0, 1))).any():
        raise ValueError("OSCA P-values must be between 0 and 1, missing or -9")
    n_tested = int(p_values.notna().sum())
    if not n_tested:
        raise ValueError("No valid EWAS P-values; no summary outputs were written")

    results["MODEL"] = model
    significant = results.loc[p_values.le(args.p_threshold)]
    outputs[0].parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(outputs[0], sep="\t", index=False)
    significant.to_csv(outputs[1], sep="\t", index=False)
    print(f"{model}: {n_tested} tested features; P <= {args.p_threshold:g}; "
          f"{len(significant)} exploratory associations")
    for output in outputs:
        print(f"Saved {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="OSCA .moment result table")
    parser.add_argument("--log", required=True, help="Matching OSCA MOMENT run log")
    parser.add_argument("--p-threshold", type=float, default=0.001, help="Unadjusted P-value cutoff")
    parser.add_argument("--out-prefix", required=True, help="Output path prefix")
    try:
        summarize_results(parser.parse_args())
    except (OSError, ValueError, pd.errors.ParserError) as error:
        parser.exit(1, f"Error: {error}\n")
