"""train FUSION models and run population-matched TWAS association tests.

inputs-dir is a fresh run directory from prepare_twas_inputs.py, containing
features.tsv, pheno/<key>.pheno and <trait>.sumstats.tsv. Failed or ineligible
features are omitted from weights.pos. A failed chromosome stops association;
summarize_twas_results.py also checks that every indexed chromosome is complete.
model settings and tool paths are supplied by 02_twas_analysis.sh.
"""

import argparse
import csv
from pathlib import Path
import re
import shlex
import shutil
import subprocess

from prepare_twas_inputs import check_outputs


def run_command(command, log, cwd=None):
    with log.open("w") as handle:
        subprocess.run([str(value) for value in command], cwd=cwd, stdout=handle,
                       stderr=subprocess.STDOUT, check=True)


def read_features(path, chromosome_count):
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    keys, ids = set(), set()
    for row in rows:
        if len(row) != 5:
            raise ValueError("Features require key, ID, chromosome, start and end, without a header")
        key, feature, chromosome, start, end = row
        if (not re.fullmatch(r"[A-Za-z0-9_-]+", key) or key in keys or feature in ids
                or not feature or any(char.isspace() for char in feature)):
            raise ValueError("Feature keys and IDs must be distinct; keys must be safe filenames")
        if not all(re.fullmatch(r"[1-9][0-9]*", value) for value in [chromosome, start, end]):
            raise ValueError("Feature coordinates must be positive integers")
        if int(chromosome) > chromosome_count or int(start) > int(end):
            raise ValueError("Feature coordinates must match the genotype chromosome set")
        keys.add(key)
        ids.add(feature)
    if not rows:
        raise ValueError("No molecular features available for training")
    return rows


def run_twas(args):
    if args.chr_count < 1 or args.cis_window < 0 or args.crossval < 2 or args.min_variants < 1:
        raise ValueError("Invalid chromosome count, cis window, CV folds or minimum variant count")
    models = args.models.split(",")
    if len(models) != len(set(models)) or not set(models).issubset({"top1", "lasso", "enet"}):
        raise ValueError("Models must be distinct choices from top1, lasso and enet")
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", args.trait):
        raise ValueError("Use a filename-safe agronomic trait name")
    inputs, output = Path(args.inputs_dir).resolve(), Path(args.out_dir).resolve()
    fusion, ldref = Path(args.fusion_dir).resolve(), Path(args.ldref_dir).resolve()
    genotype = Path(args.geno_prefix).resolve()
    plink, rscript = shutil.which("plink"), shutil.which("Rscript")
    if not plink or not rscript:
        raise ValueError("PLINK and Rscript must be available on PATH")
    gcta = shutil.which(args.gcta)
    if not gcta:
        raise ValueError("GCTA must be an executable file path or a command available on PATH")
    gcta = Path(gcta).resolve()
    features = read_features(inputs / "features.tsv", args.chr_count)
    statistics = inputs / f"{args.trait}.sumstats.tsv"
    required = [statistics, gcta, fusion / "FUSION.compute_weights.R", fusion / "FUSION.assoc_test.R"]
    required += [Path(f"{genotype}.{suffix}") for suffix in ["bed", "bim", "fam"]]
    required += [inputs / "pheno" / f"{row[0]}.pheno" for row in features]
    if any(not path.is_file() for path in required):
        raise ValueError("Missing genotype, phenotype, GWAS input, GCTA executable or FUSION script")
    weights, logs = output / "weights", output / "logs"
    results = output / "results" / args.trait
    index = output / "weights.pos"
    check_outputs([index], [inputs / "features.tsv", *required])
    for directory in [weights, logs, results, ldref]:
        directory.mkdir(parents=True, exist_ok=True)
    plink_options = ["--chr-set", args.chr_count, "--allow-extra-chr"]
    fusion_plink = shlex.join([plink, *map(str, plink_options)])
    trained = []

    with index.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["WGT", "ID", "CHR", "P0", "P1"])
        for key, feature, chromosome, start, end in features:
            cis = inputs / f"{key}.cis"
            weight = inputs / f"{key}.wgt.RDat"
            weight.unlink(missing_ok=True)
            try:
                run_command([plink, "--bfile", genotype, *plink_options, "--chr", chromosome,
                             "--from-bp", max(1, int(start) - args.cis_window),
                             "--to-bp", int(end) + args.cis_window, "--make-bed", "--out", cis],
                            logs / f"{key}.plink.log")
                with Path(f"{cis}.bim").open() as variants:
                    n_variants = sum(bool(line.strip()) for line in variants)
                if n_variants < args.min_variants:
                    continue
                run_command([rscript, "FUSION.compute_weights.R", "--bfile", cis,
                             "--pheno", inputs / "pheno" / f"{key}.pheno",
                             "--tmp", inputs / f"{key}.fusion", "--out", inputs / key,
                             "--models", args.models, "--crossval", args.crossval,
                             "--PATH_plink", fusion_plink, "--PATH_gcta", shlex.quote(str(gcta))],
                            logs / f"{key}.training.log", cwd=fusion)
            except subprocess.CalledProcessError:
                print(f"Skipped {feature}: external command failed; see {logs}")
                continue
            if weight.is_file() and weight.stat().st_size:
                weight.replace(weights / weight.name)
                writer.writerow([weight.name, feature, chromosome, start, end])
                trained.append(chromosome)
    if not trained:
        raise ValueError("No prediction weights generated; inspect the feature training logs")

    for chromosome in sorted(set(trained), key=int):
        run_command([plink, "--bfile", genotype, *plink_options, "--chr", chromosome,
                     "--make-bed", "--out", ldref / f"chr{chromosome}"],
                    logs / f"ldref.chr{chromosome}.log")
        run_command([rscript, "FUSION.assoc_test.R", "--sumstats", statistics,
                     "--weights", index, "--weights_dir", weights, "--ref_ld_chr", ldref / "chr",
                     "--chr", chromosome, "--out", inputs / f"{args.trait}.chr{chromosome}.dat"],
                    results / f"{args.trait}.chr{chromosome}.log", cwd=fusion)
    print(f"{args.trait}: {len(trained)}/{len(features)} features trained; "
          f"{len(set(trained))} chromosomes analyzed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs-dir", required=True, help="Fresh directory of prepared TWAS inputs")
    parser.add_argument("--geno-prefix", required=True, help="Population PLINK BED/BIM/FAM prefix")
    parser.add_argument("--fusion-dir", required=True, help="Directory containing the FUSION R scripts")
    parser.add_argument("--gcta", required=True, help="FUSION-compatible GCTA executable")
    parser.add_argument("--trait", required=True, help="Agronomic trait name")
    parser.add_argument("--out-dir", required=True, help="Weight and result output directory")
    parser.add_argument("--ldref-dir", required=True, help="Population LD reference directory")
    parser.add_argument("--chr-count", type=int, required=True)
    parser.add_argument("--cis-window", type=int, default=200000)
    parser.add_argument("--models", default="top1,lasso,enet")
    parser.add_argument("--crossval", type=int, default=5)
    parser.add_argument("--min-variants", type=int, default=5)
    try:
        run_twas(parser.parse_args())
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"Error: {error}\n")
