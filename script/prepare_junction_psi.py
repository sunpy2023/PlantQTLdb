"""extract a junction PSI matrix from a VOILA modulize TSV/CSV report.

usage: python3 prepare_junction_psi.py -i junctions.tsv -o junction_psi.tsv
requires pandas. Leading '#' metadata lines are skipped.
the id joins lsv_id, event_id, junction_name and junction_coord with underscores.
sample order follows the report; .sorted_median_psi/_median_psi suffixes are removed.
plain sample columns following the standard 15 annotation columns are also supported.
rows without lsv_id and exact duplicate output rows are removed; conflicting
duplicate ids are rejected. Missing PSI values are written as NA, never zero.
semicolon-delimited PSI values are checked and retained without averaging;
these are not scalar phenotypes and require handling during QTL preprocessing.
"""

import argparse
import csv
from pathlib import Path
import re

import pandas as pd


ANNOTATION_COLUMNS = [
    "module_id", "gene_id", "gene_name", "seqid", "strand", "lsv_id", "event_id",
    "complex", "module_event_combination", "denovo", "junction_name", "junction_coord",
    "event_non_changing", "event_changing", "junction_changing",
]
ID_COLUMNS = ["lsv_id", "event_id", "junction_name", "junction_coord"]
MISSING = {"", "NA", "NaN", "nan", "N/A", "."}


def prepare_matrix(input_path, output_path):
    input_path, output_path = Path(input_path), Path(output_path)
    if (input_path.resolve() == output_path.resolve()
            or (output_path.exists() and output_path.samefile(input_path))):
        raise ValueError("Input and output must be different files")
    delimiter = "," if input_path.suffix.lower() == ".csv" else "\t"
    with input_path.open(newline="") as handle:
        first = next((line for line in handle if line.strip() and not line.startswith("#")), None)
    if first is None:
        raise ValueError("Empty junction report")
    header = next(csv.reader([first], delimiter=delimiter))
    if any(not name.strip() for name in header) or len(header) != len(set(header)):
        raise ValueError("Empty or duplicate input column names")
    if not set(ID_COLUMNS).issubset(header):
        raise ValueError("Missing lsv_id/event_id/junction_name/junction_coord columns")
    psi_columns = [name for name in header if name.endswith("_median_psi")]
    if psi_columns:
        samples = [re.sub(r"(?:\.sorted)?_median_psi$", "", name) for name in psi_columns]
    else:
        if header[:len(ANNOTATION_COLUMNS)] != ANNOTATION_COLUMNS:
            raise ValueError("No median_psi columns or standard annotation prefix found")
        psi_columns = header[len(ANNOTATION_COLUMNS):]
        if any(name.endswith(("_median_reads", "_var_psi")) for name in psi_columns):
            raise ValueError("Expected PSI sample columns after annotations")
        samples = psi_columns
    if (not samples or len(samples) != len(set(samples))
            or any(not name or name == "id" or any(char.isspace() for char in name)
                   for name in samples)):
        raise ValueError("Missing, invalid or duplicate output sample names")

    data = pd.read_csv(input_path, sep=delimiter, comment="#", dtype=str,
                       keep_default_na=False, usecols=ID_COLUMNS + psi_columns)
    for column in ID_COLUMNS:
        data[column] = data[column].str.strip()
    n_missing_lsv = int(data["lsv_id"].isin(MISSING).sum())
    data = data.loc[~data["lsv_id"].isin(MISSING)].copy()
    if data.empty:
        raise ValueError("No junction rows with lsv_id to export")
    ids = data[ID_COLUMNS].copy()
    for column in ID_COLUMNS:
        ids[column] = ids[column].mask(ids[column].isin(MISSING), "NA")
    multivalue_count = 0
    for column in psi_columns:
        values = data[column].str.strip()
        multivalue_count += int(values.str.contains(";", regex=False).sum())
        tokens = values.str.split(";").explode().str.strip()
        numeric = pd.to_numeric(tokens.mask(tokens.isin(MISSING)), errors="raise")
        if not (numeric.isna() | numeric.between(0, 1)).all():
            raise ValueError(f"PSI must be within [0, 1] or missing: {column}")
        data[column] = values.mask(values.isin(MISSING), "NA")
    output = data[psi_columns].copy()
    output.columns = samples
    output.insert(0, "id", ids.agg("_".join, axis=1))
    n_duplicates = len(output) - len(output.drop_duplicates())
    output = output.drop_duplicates()
    if output["id"].duplicated().any():
        example = output.loc[output["id"].duplicated(), "id"].iloc[0]
        raise ValueError(f"Conflicting PSI records share the same id: {example}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, sep="\t", index=False)
    print(f"Exported {len(output)} junction traits and {len(samples)} samples to {output_path}")
    print(f"Removed {n_missing_lsv} rows without lsv_id and {n_duplicates} exact duplicate rows")
    if multivalue_count:
        print(f"Found {multivalue_count} semicolon-delimited PSI fields before deduplication; "
              "handle these during QTL preprocessing")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, help="VOILA junction report")
    parser.add_argument("-o", "--output", required=True, help="Output PSI matrix TSV")
    args = parser.parse_args()
    prepare_matrix(args.input, args.output)
