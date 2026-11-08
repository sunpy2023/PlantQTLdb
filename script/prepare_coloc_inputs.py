"""prepare coloc inputs from public cis-QTL tables, a QC VCF and GEMMA results.

requires pandas and pysam. Matrix columns identify the QTL samples; the trait
file has one value per FAM row, with NA/-9 for missing values. VCF and FAM sample
sets must match. Chromosomes use the same positive integer codes in all inputs.
QTL IDs are CHROM:POS:REF:ALT; variants are checked against VCF/BIM alleles and
renamed to BIM IDs. Different biallelic variants may share the same position.
one VCF scan calculates MAF in the QTL cohort and in observed GWAS trait samples.
QTL varbeta = (beta / t)^2; GWAS varbeta = se^2. Monomorphic or unobserved variants
and tests without recoverable positive variance are omitted, without a P filter.
eQTL N comes from matrix columns; sQTL/APA N and df must match the observed matrix
counts. GWAS N is the observed trait count, matching the original coloc workflow.
outputs: qtl.coloc.tsv.gz, gwas.coloc.tsv.gz and gwas.clump.tsv (SNP, P).
"""

import argparse
import csv
import gzip
import math
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
import pysam

from prepare_ewas_inputs import numeric_values, validate_ids
from prepare_qtl_phenotypes import MISSING, read_table
from prepare_smr_inputs import QTL_FILES, number, positive_integer, table_rows
from prepare_twas_inputs import check_outputs


def sample_metadata(args):
    rows = [line.split() for line in Path(args.fam).read_text().splitlines()]
    if not rows or any(len(row) != 6 for row in rows):
        raise ValueError("FAM requires six fields per row, without a header")
    samples = pd.Series([row[1] for row in rows])
    validate_ids(samples, "FAM individual IDs")
    rows = [line.split() for line in Path(args.phenotype).read_text().splitlines()]
    if len(rows) != len(samples) or any(len(row) != 1 for row in rows):
        raise ValueError("Phenotype requires one value per FAM row; no header, NA for missing")
    trait = numeric_values(pd.Series([row[0] for row in rows]), MISSING | {"None"})
    trait = trait.mask(trait.eq(-9))
    if trait.count() < 3 or trait.nunique() < 2:
        raise ValueError("At least three observed samples and a variable trait are required")

    matrix = read_table(args.matrix)
    validate_ids(matrix.iloc[:, 0], "Molecular feature IDs")
    qtl_samples = pd.Series(matrix.columns[1:], dtype=str)
    validate_ids(qtl_samples, "Molecular sample IDs")
    if len(qtl_samples) < 3 or not qtl_samples.isin(samples).all():
        raise ValueError("At least three QTL samples must match FAM IDs; unmatched columns are not dropped")
    if matrix.eq("").any().any():
        raise ValueError("Incomplete molecular matrix; use NA for missing values")
    values = matrix.iloc[:, 1:].apply(numeric_values, missing=MISSING | {"None"})
    if args.molecular_type == "expression" and values.isna().any().any():
        raise ValueError("Expression input must match the complete matrix used for eQTL mapping")
    counts = dict(zip(matrix.iloc[:, 0], values.count(axis=1)))
    return set(samples), set(qtl_samples), set(samples[trait.notna()]), counts


def variant_metadata(args, samples, qtl_samples, gwas_samples):
    variants, positions = {}, {}
    with Path(args.bim).open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) != 6:
                raise ValueError("BIM requires six fields per row")
            chrom, snp, _, bp, a1, a2 = fields
            positive_integer(chrom)
            positive_integer(bp)
            a1, a2 = a1.upper(), a2.upper()
            key = (chrom, bp, tuple(sorted((a1, a2))))
            if (snp in variants or key in positions or snp in MISSING
                    or a1 == a2 or {a1, a2} & (MISSING | {"0"})):
                raise ValueError(f"Duplicate BIM ID/variant or invalid alleles: {snp}")
            record = [snp, chrom, bp, a1, a2, None, None, None, None]
            variants[snp] = positions[key] = record
    if not variants:
        raise ValueError("BIM is empty")
    with pysam.VariantFile(args.vcf) as vcf:
        if set(vcf.header.samples) != samples:
            raise ValueError("VCF and FAM must contain the same sample IDs")
        selected = qtl_samples | gwas_samples
        vcf.subset_samples([sample for sample in vcf.header.samples if sample in selected])
        for rec in vcf:
            ref = rec.ref.upper()
            alt = rec.alts[0].upper() if len(rec.alts or []) == 1 else None
            if alt is None or ref == alt:
                raise ValueError(f"Expected a biallelic VCF variant: {rec.chrom}:{rec.pos}")
            record = positions.get((rec.chrom, str(rec.pos), tuple(sorted((ref, alt)))))
            if record is None:
                continue
            if record[5] is not None:
                raise ValueError(f"Duplicate VCF variant: {record[0]}")
            q_alt = q_total = g_alt = g_total = 0
            for sample, call in rec.samples.items():
                gt = call.get("GT") or ()
                if any(allele not in (None, 0, 1) for allele in gt):
                    raise ValueError(f"Non-biallelic genotype: {record[0]}")
                alt_count = gt.count(1)
                total = len(gt) - gt.count(None)
                if sample in qtl_samples:
                    q_alt += alt_count
                    q_total += total
                if sample in gwas_samples:
                    g_alt += alt_count
                    g_total += total
            q_maf = min(q_alt, q_total - q_alt) / q_total if q_total else None
            g_maf = min(g_alt, g_total - g_alt) / g_total if g_total else None
            record[5:] = [ref, alt, q_maf, g_maf]
    if any(record[5] is None for record in variants.values()):
        raise ValueError("Every BIM variant must have matching VCF coordinates and alleles")
    return variants, positions


def write_qtl(args, source, positions, counts, output):
    id_column = {"expression": "gene", "splicing": "junction_id", "apa": "apa_id"}[args.molecular_type]
    required = ["SNP", id_column, "beta", "t-stat", "p-value"]
    if args.molecular_type != "expression":
        required += ["n_samples", "df"]
    kept = skipped = 0
    degrees = {}
    with gzip.open(output, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["feature", "SNP", "chr", "position", "beta", "varbeta", "MAF", "p",
                         "n_samples", "df"])
        for row in table_rows(source, required):
            feature = row[id_column]
            if feature not in counts or counts[feature] <= 2:
                raise ValueError(f"QTL feature needs at least three observed matrix values: {feature}")
            n, df = counts[feature], "NA"
            if "n_samples" in row:
                n, df = number(row["n_samples"]), number(row.get("df", "NA"))
                if (n != counts[feature] or df is None or df != int(df) or not 0 < df <= n - 2):
                    raise ValueError(f"QTL n_samples/df disagree with the molecular matrix: {feature}")
                n, df = int(n), int(df)
                if degrees.setdefault(feature, df) != df:
                    raise ValueError(f"Inconsistent QTL df within feature: {feature}")
            parts = row["SNP"].rsplit(":", 3)
            if len(parts) != 4:
                raise ValueError(f"Expected CHROM:POS:REF:ALT QTL ID: {row['SNP']}")
            chrom, bp, ref, alt = parts
            ref, alt = ref.upper(), alt.upper()
            record = positions.get((chrom, bp, tuple(sorted((ref, alt)))))
            if record is None or [ref, alt] != record[5:7]:
                raise ValueError(f"QTL variant does not match the QC VCF: {row['SNP']}")
            beta, t, p = (number(row[column]) for column in ["beta", "t-stat", "p-value"])
            if p is not None and not 0 <= p <= 1:
                raise ValueError(f"Invalid QTL P-value: {row['SNP']}")
            if beta is None or t is None or p is None or t == 0:
                skipped += 1
                continue
            se = abs(beta / t)
            variance = se * se
            if not math.isfinite(variance) or variance <= 0 or not record[7]:
                skipped += 1
                continue
            writer.writerow([feature, record[0], chrom, bp, beta, variance, record[7], p, n, df])
            kept += 1
    print(f"QTL: {kept} cis pairs retained; {skipped} unrecoverable or monomorphic records omitted")


def write_gwas(source, variants, n_gwas, output, clump_file):
    required = ["rs", "chr", "ps", "allele1", "allele0", "beta", "se", "p_wald"]
    seen, kept, skipped = set(), 0, 0
    with gzip.open(output, "wt", newline="") as handle, clump_file.open("w", newline="") as clump:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        lead_writer = csv.writer(clump, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "chr", "position", "beta", "varbeta", "MAF", "p", "n_samples"])
        lead_writer.writerow(["SNP", "P"])
        for row in table_rows(source, required):
            snp = row["rs"]
            if snp in seen or snp not in variants:
                raise ValueError(f"Duplicate GWAS ID or variant absent from BIM: {snp}")
            seen.add(snp)
            record = variants[snp]
            if (row["chr"] != record[1] or row["ps"] != record[2]
                    or (row["allele1"].upper() == row["allele0"].upper())
                    or {row["allele1"].upper(), row["allele0"].upper()} != set(record[3:5])):
                raise ValueError(f"GWAS/BIM coordinates or alleles disagree: {snp}")
            beta, se, p = (number(row[column]) for column in ["beta", "se", "p_wald"])
            if p is not None and not 0 <= p <= 1:
                raise ValueError(f"Invalid GWAS P-value: {snp}")
            if beta is None or se is None or p is None or se <= 0 or not record[8]:
                skipped += 1
                continue
            variance = se * se
            if not math.isfinite(variance) or variance <= 0:
                skipped += 1
                continue
            writer.writerow([snp, record[1], record[2], beta, variance, record[8], p, n_gwas])
            lead_writer.writerow([snp, p])
            kept += 1
    if not kept:
        raise ValueError("No valid GWAS summary statistics remain")
    print(f"GWAS: {n_gwas} observed samples; {kept} records retained; {skipped} invalid records omitted")


def prepare_inputs(args):
    source = Path(args.qtl_root) / QTL_FILES[args.molecular_type]
    out_dir = Path(args.out_dir)
    names = ["qtl.coloc.tsv.gz", "gwas.coloc.tsv.gz", "gwas.clump.tsv"]
    outputs = [out_dir / name for name in names]
    check_outputs(outputs, [source, args.matrix, args.vcf, args.bim, args.fam, args.gwas, args.phenotype])
    samples, qtl_samples, gwas_samples, counts = sample_metadata(args)
    variants, positions = variant_metadata(args, samples, qtl_samples, gwas_samples)
    out_dir.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".prepare-coloc-", dir=out_dir) as temporary:
        paths = [Path(temporary) / name for name in names]
        write_qtl(args, source, positions, counts, paths[0])
        write_gwas(args.gwas, variants, len(gwas_samples), paths[1], paths[2])
        for path, output in zip(paths, outputs):
            path.replace(output)
    print(f"Prepared coloc inputs: {len(qtl_samples)} QTL samples; {len(variants)} QC variants; {out_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecular-type", choices=QTL_FILES, required=True)
    parser.add_argument("--qtl-root", required=True, help="molQTL directory containing eqtl/sqtl/apaqtl")
    parser.add_argument("--matrix", required=True, help="Normalized molecular matrix used for QTL mapping")
    parser.add_argument("--vcf", required=True, help="QC VCF matching the BIM/FAM data")
    parser.add_argument("--bim", required=True, help="QC BIM defining unique variant IDs, coordinates and alleles")
    parser.add_argument("--fam", required=True, help="FAM defining sample IDs and phenotype row order")
    parser.add_argument("--gwas", required=True, help="Full GEMMA association table")
    parser.add_argument("--phenotype", required=True, help="One trait column in FAM order; NA/-9 for missing")
    parser.add_argument("--out-dir", required=True, help="Working directory for coloc inputs")
    try:
        prepare_inputs(parser.parse_args())
    except (OSError, ValueError, pd.errors.ParserError) as error:
        parser.exit(1, f"Error: {error}\n")
