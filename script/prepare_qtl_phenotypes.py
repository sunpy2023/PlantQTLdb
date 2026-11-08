"""select QTL samples, filter raw molecular phenotypes and prepare matching locations.

requires pandas. Input matrices may be CSV/TSV, with trait IDs in the first
column and sample IDs in the remaining columns. Sample IDs must already match
the genotype IDs; the sample list determines the subset and output order.
expression uses GFF3 gene coordinates, PSI uses the original VOILA report,
and PDUI uses the locations produced by prepare_apa_matrix.py.
chromosome mappings are applied exactly; unmapped names are retained.
coordinate values are retained from their sources; rows without valid bounds
are removed from both the matrix and locations, without shifting endpoints.
missing values stay NA. Valid semicolon-delimited PSI fields become NA,
matching the original filter_psi.py behavior, without averaging their values.
outputs: <prefix>.filtered.tsv and <prefix>.locations.tsv, in matching row order.
"""

import argparse
import csv
import gzip
import math
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

from prepare_junction_psi import ID_COLUMNS, MISSING


def read_table(path, columns=None):
    with Path(path).open() as handle:
        first = next((line for line in handle if line.strip() and not line.startswith("#")), "")
    delimiter = "\t" if "\t" in first else ","
    header = next(csv.reader([first], delimiter=delimiter), [])
    if (not header or len(header) != len(set(header))
            or any(not name.strip() for name in header)):
        raise ValueError(f"Empty or duplicate table column names: {path}")
    if columns is not None and not set(columns).issubset(header):
        raise ValueError(f"Missing required columns {columns}: {path}")
    table = pd.read_csv(path, sep=delimiter, dtype=str, keep_default_na=False,
                        comment="#", usecols=columns)
    if not isinstance(table.index, pd.RangeIndex):
        raise ValueError(f"Unexpected number of fields in table: {path}")
    if table.isna().any().any():
        raise ValueError(f"Incomplete table rows: {path}")
    return table if columns is None else table[columns]


def read_samples(path):
    samples = [line.strip() for line in Path(path).read_text().splitlines() if line.strip()]
    if (len(samples) < 2 or len(samples) != len(set(samples))
            or any(name == "id" or any(char.isspace() for char in name) for name in samples)):
        raise ValueError("Sample list requires at least two distinct IDs, one per line, without a header")
    return samples


def read_values(path, samples, phenotype_type):
    table = read_table(path)
    ids = table.iloc[:, 0]
    if (ids.isin(MISSING).any() or ids.duplicated().any()
            or ids.str.contains(r"\s", regex=True).any()):
        raise ValueError("Phenotype IDs must be nonempty, distinct and contain no whitespace")
    table = table.set_index(table.columns[0])
    if not set(samples).issubset(table.columns):
        absent = sorted(set(samples) - set(table.columns))
        raise ValueError(f"Requested samples absent from phenotype matrix: {absent}")
    values = table[samples].copy()
    multivalue_count = 0
    for column in samples:
        text = values[column].str.strip()
        multi = text.str.contains(";", regex=False)
        tokens = text.str.split(";").explode().str.strip() if phenotype_type == "psi" else text
        numeric = pd.to_numeric(tokens.mask(tokens.isin(MISSING)), errors="raise")
        valid = numeric.isna() | (numeric.ge(0) & numeric.lt(float("inf")))
        if phenotype_type != "expression":
            valid &= numeric.isna() | numeric.le(1)
        if not valid.all():
            raise ValueError(f"Invalid {phenotype_type} value in sample {column}")
        multivalue_count += int(multi.sum())
        values[column] = pd.to_numeric(text.mask(text.isin(MISSING) | multi), errors="raise")
    values.index.name = "id"
    return values, multivalue_count


def gene_locations(path, ids):
    wanted, rows = set(ids), []
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rt") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Expected nine GFF3 columns: {path}:{number}")
            if fields[2] != "gene":
                continue
            attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            gene = unquote(attrs.get("ID", ""))
            aliases = [gene]
            if gene.startswith("gene-ZEAMMB73_"):
                aliases.append(gene[len("gene-ZEAMMB73_"):])
            rows.extend([name, fields[0], fields[3], fields[4]] for name in aliases if name in wanted)
    return pd.DataFrame(rows, columns=["geneid", "chr", "s1", "s2"])


def junction_locations(path, ids):
    meta = read_table(path, [*ID_COLUMNS, "seqid"])
    for column in ID_COLUMNS:
        text = meta[column].str.strip()
        meta[column] = text.mask(text.isin(MISSING), "NA")
    meta["geneid"] = meta[ID_COLUMNS].agg("_".join, axis=1)
    meta = meta.loc[meta["geneid"].isin(ids)].copy()
    coordinates = meta["junction_coord"].str.extract(r"^(?P<s1>\d+)-(?P<s2>\d+)$")
    return pd.DataFrame({"geneid": meta["geneid"], "chr": meta["seqid"],
                         "s1": coordinates["s1"], "s2": coordinates["s2"]})


def align_locations(locations, ids, chrom_map, require_all):
    locations = locations.loc[locations["geneid"].isin(ids)].drop_duplicates()
    if locations["geneid"].duplicated().any():
        raise ValueError("A phenotype ID maps to conflicting locations")
    locations = locations.set_index("geneid")
    if require_all and not set(ids).issubset(locations.index):
        raise ValueError("Phenotype IDs are absent from the supplied location table or junction report")
    locations = locations.reindex(ids)
    if chrom_map:
        with Path(chrom_map).open(newline="") as handle:
            rows = [row for row in csv.reader(handle, delimiter="\t") if row]
        if (not rows or any(len(row) != 2 or any(not value or any(c.isspace() for c in value)
                                                for value in row) for row in rows)
                or len({row[0] for row in rows}) != len(rows)):
            raise ValueError("Chromosome map requires distinct old IDs and two columns, without a header")
        mapping = dict(rows)
        locations["chr"] = locations["chr"].map(lambda name: mapping.get(name, name))
    for column in ["s1", "s2"]:
        locations[column] = pd.to_numeric(locations[column], errors="coerce")
    valid = (locations["chr"].notna() & ~locations["chr"].isin(MISSING)
             & ~locations["chr"].str.contains(r"\s", regex=True, na=True)
             & locations["s1"].gt(0) & locations["s1"].le(locations["s2"])
             & locations["s1"].mod(1).eq(0) & locations["s2"].mod(1).eq(0))
    locations = locations.loc[valid].copy()
    locations[["s1", "s2"]] = locations[["s1", "s2"]].astype("int64")
    locations.index.name = "geneid"
    return locations


def prepare_phenotypes(args):
    thresholds = [args.max_missing_frac, args.min_mean_expression, args.expression_cutoff,
                  args.min_expressed_frac, args.min_intermediate_frac, args.intermediate_low,
                  args.intermediate_high, args.min_range, args.min_variance]
    bounded = [args.max_missing_frac, args.min_expressed_frac, args.min_intermediate_frac,
               args.intermediate_low, args.intermediate_high, args.min_range, args.min_variance]
    if (not all(math.isfinite(value) and value >= 0 for value in thresholds)
            or any(value > 1 for value in bounded)
            or args.intermediate_low > args.intermediate_high):
        raise ValueError("Invalid phenotype filtering thresholds")
    source = {"expression": args.annotation, "psi": args.junction_meta, "pdui": args.locations}[args.type]
    if not source:
        raise ValueError("Provide --annotation for expression, --junction-meta for PSI "
                         "or --locations for PDUI")
    inputs = [Path(value) for value in [args.input, args.samples, source, args.chrom_map] if value]
    outputs = [Path(f"{args.out_prefix}.{suffix}.tsv") for suffix in ["filtered", "locations"]]
    for index, output in enumerate(outputs):
        for other in inputs + outputs[index + 1:]:
            if (output.resolve() == other.resolve()
                    or (output.exists() and other.exists() and output.samefile(other))):
                raise ValueError("Outputs must be distinct and must not overwrite input files")
    samples = read_samples(args.samples)
    values, multivalue_count = read_values(args.input, samples, args.type)
    n_input = len(values)
    observed = values.notna().sum(axis=1)
    keep = values.isna().mean(axis=1).le(args.max_missing_frac)
    if args.type == "expression":
        expressed = values.gt(args.expression_cutoff).sum(axis=1).div(observed)
        keep &= values.mean(axis=1).gt(args.min_mean_expression) & expressed.ge(args.min_expressed_frac)
    else:
        intermediate = (values.ge(args.intermediate_low) & values.le(args.intermediate_high)).sum(axis=1)
        keep &= (intermediate.div(observed).ge(args.min_intermediate_frac)
                 & (values.max(axis=1) - values.min(axis=1)).ge(args.min_range)
                 & values.var(axis=1, ddof=1).ge(args.min_variance))
    values = values.loc[keep]
    n_filtered = len(values)
    if args.type == "expression":
        locations = gene_locations(source, values.index)
    elif args.type == "psi":
        locations = junction_locations(source, values.index)
    else:
        locations = read_table(source, ["geneid", "chr", "s1", "s2"])
    locations = align_locations(locations, values.index, args.chrom_map, args.type != "expression")
    values = values.loc[locations.index]
    if values.empty:
        raise ValueError("No phenotypes remain after filtering and location matching; "
                         "no outputs were written")
    outputs[0].parent.mkdir(parents=True, exist_ok=True)
    values.to_csv(outputs[0], sep="\t", na_rep="NA", index_label="id")
    locations.to_csv(outputs[1], sep="\t", index_label="geneid")
    print(f"{args.type}: {len(samples)} samples; {n_input} input traits; {n_filtered} pass phenotype QC; "
          f"{len(values)} have valid locations")
    if multivalue_count:
        print(f"Treated {multivalue_count} semicolon-delimited PSI fields as NA, without averaging")
    for output in outputs:
        print(f"Saved {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--type", choices=["expression", "psi", "pdui"], required=True)
    parser.add_argument("--input", required=True, help="Raw phenotype CSV/TSV")
    parser.add_argument("--samples", required=True,
                        help="Matched sample IDs in output order, without a header")
    parser.add_argument("--annotation", help="GFF3 gene annotation, optionally gzipped, for expression")
    parser.add_argument("--junction-meta", help="Original VOILA junctions.tsv/CSV for PSI")
    parser.add_argument("--locations", help="APA location TSV with geneid, chr, s1 and s2 columns")
    parser.add_argument("--chrom-map", help="Old chromosome ID<TAB>standard ID, without a header")
    parser.add_argument("--max-missing-frac", type=float, default=0.2)
    parser.add_argument("--min-mean-expression", type=float, default=0.1)
    parser.add_argument("--expression-cutoff", type=float, default=0.1)
    parser.add_argument("--min-expressed-frac", type=float, default=0.1)
    parser.add_argument("--min-intermediate-frac", type=float, default=0.1)
    parser.add_argument("--intermediate-low", type=float, default=0.1)
    parser.add_argument("--intermediate-high", type=float, default=0.9)
    parser.add_argument("--min-range", type=float, default=0.1)
    parser.add_argument("--min-variance", type=float, default=0.001)
    parser.add_argument("--out-prefix", required=True, help="Output path prefix")
    try:
        prepare_phenotypes(parser.parse_args())
    except (OSError, ValueError, pd.errors.ParserError) as error:
        parser.exit(1, f"Error: {error}\n")
