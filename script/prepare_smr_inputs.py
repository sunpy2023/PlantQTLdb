"""prepare SMR query and GWAS files from the public molQTL and GEMMA outputs.

requires pandas. QTL IDs must be CHROM:POS:REF:ALT, as written by the mapping
scripts; BIM chromosomes must use the same positive integer codes. Variants are
matched by chromosome, position and alleles, then renamed to BIM IDs. QTL beta
is aligned to BIM A1; SE = abs(beta / t). No strand complementation is inferred.
all available cis pairs are retained for HEIDI, including nonsignificant pairs.
feature IDs and starts are used as probe labels and anchors; strand is unknown.
gene labels follow the gene prefixes in splicing/APA feature IDs.
per-feature N/df, when present, are saved separately, without assigning a global N.
GWAS A1 and its frequency follow GEMMA; N = observed trait samples - n_miss.
outputs: qtl.query.tsv, qtl.features.tsv and <trait>.gwas.ma.
"""

import argparse
import csv
import math
from pathlib import Path
import re
from tempfile import TemporaryDirectory

import pandas as pd

from prepare_ewas_inputs import numeric_values, validate_ids
from prepare_qtl_phenotypes import MISSING, read_table
from prepare_twas_inputs import check_outputs


QTL_FILES = {"expression": "eqtl/cis_egenes.all_cis.tsv",
             "splicing": "sqtl/cis_splicing_traits.all_cis.tsv",
             "apa": "apaqtl/cis_apa_traits.all_cis.tsv"}
QUERY_COLUMNS = ["SNP", "Chr", "BP", "A1", "A2", "Freq", "Probe", "Probe_Chr",
                 "Probe_bp", "Gene", "Orientation", "b", "se", "p"]


def table_rows(path, required, whitespace=False):
    with Path(path).open(newline="") as handle:
        rows = (line.split() for line in handle) if whitespace else csv.reader(handle, delimiter="\t")
        header = next(rows, [])
        if len(header) != len(set(header)) or not set(required).issubset(header):
            raise ValueError(f"Missing required or duplicate columns in {path}: {required}")
        for number, fields in enumerate(rows, 2):
            if not fields:
                continue
            if len(fields) != len(header) or any(not value.strip() for value in fields):
                raise ValueError(f"Incomplete or malformed row: {path}:{number}")
            yield dict(zip(header, (value.strip() for value in fields)))


def number(value):
    if value in MISSING:
        return None
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("Numeric values must be finite or missing")
    return value


def positive_integer(value):
    if not re.fullmatch(r"[1-9][0-9]*", value):
        raise ValueError(f"Expected a positive integer chromosome or coordinate: {value}")
    return int(value)


def genotype_metadata(bim, frequency):
    variants, positions = {}, {}
    with Path(bim).open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) != 6:
                raise ValueError("BIM requires six fields per row")
            chrom, snp, _, bp, a1, a2 = fields
            positive_integer(chrom)
            positive_integer(bp)
            a1, a2 = a1.upper(), a2.upper()
            if snp in variants or snp in MISSING or a1 == a2 or {a1, a2} & (MISSING | {"0"}):
                raise ValueError(f"Duplicate BIM ID or invalid alleles: {snp}")
            key = (chrom, bp, tuple(sorted((a1, a2))))
            if key in positions:
                raise ValueError(f"Ambiguous BIM chromosome, position and alleles: {snp}")
            record = [snp, chrom, bp, a1, a2, None]
            variants[snp] = positions[key] = record
    if not variants:
        raise ValueError("BIM is empty")
    for row in table_rows(frequency, ["CHR", "SNP", "A1", "A2", "MAF"], whitespace=True):
        snp = row["SNP"]
        if snp not in variants:
            raise ValueError(f"Frequency variant absent from BIM: {snp}")
        record = variants[snp]
        freq = number(row["MAF"])
        if (record[5] is not None or row["CHR"] != record[1] or freq is None or not 0 <= freq <= 1
                or (row["A1"].upper(), row["A2"].upper()) != tuple(record[3:5])):
            raise ValueError(f"Inconsistent frequency metadata; use --keep-allele-order: {snp}")
        record[5] = freq
    if any(record[5] is None for record in variants.values()):
        raise ValueError("Every BIM variant requires a matching allele frequency")
    return variants, positions


def observed_samples(fam, phenotype):
    rows = [line.split() for line in Path(fam).read_text().splitlines()]
    if not rows or any(len(row) != 6 for row in rows):
        raise ValueError("FAM requires six fields per row, without a header")
    validate_ids(pd.Series([row[1] for row in rows]), "FAM individual IDs")
    values = [line.split() for line in Path(phenotype).read_text().splitlines()]
    if len(values) != len(rows) or any(len(row) != 1 for row in values):
        raise ValueError("Phenotype requires one value per FAM row, without a header; use NA for missing")
    values = numeric_values(pd.Series([row[0] for row in values]), MISSING | {"None"})
    values = values.mask(values.eq(-9)).dropna()
    if len(values) < 3 or values.nunique() < 2:
        raise ValueError("At least three observed samples and a variable trait are required")
    return len(values)


def write_qtl(args, source, positions, query, feature_file):
    id_column = {"expression": "gene", "splicing": "junction_id", "apa": "apa_id"}[args.molecular_type]
    locations = read_table(args.locations, ["geneid", "chr", "s1", "s2"])
    validate_ids(locations["geneid"], "Location feature IDs")
    locations = locations.set_index("geneid").to_dict("index")
    features, skipped = {}, 0
    with query.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(QUERY_COLUMNS)
        for row in table_rows(source, ["SNP", id_column, "beta", "t-stat", "p-value"]):
            feature = row[id_column]
            if feature in MISSING or any(char.isspace() for char in feature) or feature not in locations:
                raise ValueError(f"QTL feature requires an unambiguous matching location: {feature}")
            location = locations[feature]
            chrom, anchor = location["chr"], location["s1"]
            positive_integer(chrom)
            if positive_integer(anchor) > positive_integer(location["s2"]):
                raise ValueError(f"Feature end precedes start: {feature}")
            n, df = "NA", "NA"
            if "n_samples" in row or "df" in row:
                n, df = number(row.get("n_samples", "NA")), number(row.get("df", "NA"))
                if n is None or df is None or n != int(n) or df != int(df) or not 0 < df <= n - 2:
                    raise ValueError(f"Invalid QTL n_samples/df: {feature}")
                n, df = int(n), int(df)
            counts = features.setdefault(feature, [n, df, 0, 0])
            if counts[:2] != [n, df]:
                raise ValueError(f"Inconsistent per-feature n_samples/df: {feature}")
            beta, t, p = (number(row[column]) for column in ["beta", "t-stat", "p-value"])
            if p is not None and not 0 <= p <= 1:
                raise ValueError(f"QTL P-value outside [0, 1]: {feature}")
            if beta is None or t is None or p is None or beta == 0 or t == 0:
                skipped += 1
                continue
            se = abs(beta / t)
            if not math.isfinite(se) or se <= 0:
                skipped += 1
                continue
            parts = row["SNP"].rsplit(":", 3)
            if len(parts) != 4:
                raise ValueError(f"QTL SNP ID must be CHROM:POS:REF:ALT: {row['SNP']}")
            snp_chrom, bp, ref, alt = parts
            positive_integer(snp_chrom)
            positive_integer(bp)
            ref, alt = ref.upper(), alt.upper()
            key = (snp_chrom, bp, tuple(sorted((ref, alt))))
            if ref == alt or key not in positions or snp_chrom != chrom:
                raise ValueError(f"QTL variant mismatches BIM alleles or feature chromosome: {row['SNP']}")
            snp, _, _, a1, a2, freq = positions[key]
            counted = ref if args.counted_allele == "REF" else alt
            flipped = counted != a1
            gene = feature
            if args.molecular_type == "splicing":
                gene = feature.split(":", 1)[0]
            elif args.molecular_type == "apa":
                gene = feature.split("|", 1)[0]
            writer.writerow([snp, snp_chrom, bp, a1, a2, freq, feature, chrom, anchor,
                             gene, "NA", -beta if flipped else beta, se, p])
            counts[2] += 1
            counts[3] += int(flipped)
    retained = {feature: counts for feature, counts in features.items() if counts[2]}
    if not retained:
        raise ValueError("No QTL pairs have recoverable SE and matching genotype alleles")
    with feature_file.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Probe", "n_samples", "df", "n_variants", "n_flipped"])
        writer.writerows([feature, *counts] for feature, counts in retained.items())
    print(f"QTL: {len(retained)} features; {sum(row[2] for row in retained.values())} cis pairs; "
          f"{sum(row[3] for row in retained.values())} effects flipped; {skipped} unrecoverable rows omitted")


def write_gwas(source, variants, observed, output):
    required = ["rs", "chr", "ps", "allele1", "allele0", "af", "beta", "se", "p_wald", "n_miss"]
    seen, kept, skipped = set(), 0, 0
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "A1", "A2", "freq", "b", "se", "p", "N"])
        for row in table_rows(source, required):
            snp = row["rs"]
            if snp in seen or snp not in variants:
                raise ValueError(f"Duplicate GWAS ID or variant absent from BIM: {snp}")
            seen.add(snp)
            record = variants[snp]
            a1, a2 = row["allele1"].upper(), row["allele0"].upper()
            if (row["chr"] != record[1] or row["ps"] != record[2]
                    or a1 == a2 or {a1, a2} != set(record[3:5])):
                raise ValueError(f"GWAS and BIM positions or alleles disagree: {snp}")
            beta, se, p, freq, missing = (number(row[column])
                                         for column in ["beta", "se", "p_wald", "af", "n_miss"])
            if ((p is not None and not 0 <= p <= 1)
                    or (freq is not None and not 0 <= freq <= 1)):
                raise ValueError(f"GWAS P-value or allele frequency outside [0, 1]: {snp}")
            if any(value is None for value in [beta, se, p, freq, missing]) or se <= 0:
                skipped += 1
                continue
            if missing != int(missing) or not 0 <= missing < observed:
                raise ValueError(f"Invalid GWAS n_miss for {observed} observed trait samples: {snp}")
            writer.writerow([snp, a1, a2, freq, beta, se, p, observed - int(missing)])
            kept += 1
    if not kept:
        raise ValueError("No valid GWAS associations remain")
    print(f"GWAS: {observed} observed trait samples; {kept} records retained; "
          f"{skipped} invalid tests omitted")


def prepare_inputs(args):
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", args.trait):
        raise ValueError("Trait name must contain only letters, numbers, underscores, dots or hyphens")
    source = Path(args.qtl_root) / QTL_FILES[args.molecular_type]
    out_dir = Path(args.out_dir)
    names = ["qtl.query.tsv", "qtl.features.tsv", f"{args.trait}.gwas.ma"]
    outputs = [out_dir / name for name in names]
    check_outputs(outputs, [source, args.locations, args.bim, args.frequency,
                            args.fam, args.phenotype, args.gwas])
    variants, positions = genotype_metadata(args.bim, args.frequency)
    observed = observed_samples(args.fam, args.phenotype)
    out_dir.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".prepare-smr-", dir=out_dir) as temporary:
        paths = [Path(temporary) / name for name in names]
        write_qtl(args, source, positions, paths[0], paths[1])
        write_gwas(args.gwas, variants, observed, paths[2])
        for path, output in zip(paths, outputs):
            path.replace(output)
    for output in outputs:
        print(f"Saved {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecular-type", choices=QTL_FILES, required=True)
    parser.add_argument("--qtl-root", required=True, help="molQTL directory containing eqtl/sqtl/apaqtl")
    parser.add_argument("--locations", required=True, help="TSV with geneid, chr, s1 and s2")
    parser.add_argument("--bim", required=True, help="QC genotype BIM defining variant IDs and alleles")
    parser.add_argument("--fam", required=True, help="FAM defining phenotype row order")
    parser.add_argument("--frequency", required=True, help="PLINK --keep-allele-order --freq output")
    parser.add_argument("--counted-allele", choices=["REF", "ALT"], required=True,
                        help="Allele counted in the corresponding QTL mapping run")
    parser.add_argument("--gwas", required=True, help="Full GEMMA association table, including af and n_miss")
    parser.add_argument("--phenotype", required=True, help="One trait column in FAM order; NA/-9 for missing")
    parser.add_argument("--trait", required=True, help="Agronomic trait name")
    parser.add_argument("--out-dir", required=True, help="SMR working directory")
    try:
        prepare_inputs(parser.parse_args())
    except (OSError, ValueError, pd.errors.ParserError) as error:
        parser.exit(1, f"Error: {error}\n")
