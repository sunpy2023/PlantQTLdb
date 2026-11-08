"""prepare APA event tables, a gene-level PDUI matrix and matching locations.

usage: python3 prepare_apa_matrix.py --result-prefix /path/to/apa/results \
    --chrom-list chromosomes.txt --tx2gene tx2gene.tsv --out-dir matrices
reads results_<chromosome>/dapars2_result_temp.<chromosome>.txt in chromosome-list order.
sample order is preserved; paths and terminal _PDUI suffixes are removed from names.
gene events are identified by gene ID, analysis interval and proximal poly(A) site.
duplicate gene events are merged only when their PDUI values and strands agree;
all transcript IDs are retained and the lowest fitting error is reported.
missing PDUI stays NA. No phenotype filtering, imputation or normalization is applied.
locations retain the 1-based inclusive analysis intervals from DaPars2 Loci.
"""

import argparse
import csv
import math
from pathlib import Path
import re

from extract_apa_regions import load_mapping


RESULT_FIELDS = ["Gene", "fit_value", "Predicted_Proximal_APA", "Loci"]
EVENT_FIELDS = ["event_id", "transcript_id", "gene_id", "chr", "strand",
                "s1", "s2", "proximal_site", "fit_value"]
GENE_FIELDS = ["event_id", "gene_id", "transcript_ids", "chr", "strand",
               "s1", "s2", "proximal_site", "fit_value"]
MISSING = {"", "NA", "nan", "NaN", "N/A", "."}


def sample_names(header):
    if len(header) < 5 or header[:4] != RESULT_FIELDS:
        raise ValueError("Expected DaPars2 metadata columns followed by sample_PDUI columns")
    columns = [Path(value).name for value in header[4:]]
    if not all(value.endswith("_PDUI") for value in columns):
        raise ValueError("Sample columns must end in _PDUI")
    samples = [value[:-5] for value in columns]
    if (len(samples) != len(set(samples))
            or any(not name or any(c.isspace() for c in name)
                   or name in EVENT_FIELDS + GENE_FIELDS for name in samples)):
        raise ValueError("Invalid, duplicate or reserved sample names")
    return samples


def parse_event(row, chromosome, mapping, context):
    parts = row[0].rsplit("|", 3)
    locus = re.fullmatch(r"([^:\s]+):(\d+)-(\d+)", row[3])
    if len(parts) != 4 or locus is None:
        raise ValueError(f"Invalid DaPars2 event ID or Loci: {context}")
    tx, raw_gene, raw_chrom, strand = parts
    chrom, start, end = locus.group(1), int(locus.group(2)), int(locus.group(3))
    proximal, fit = int(row[2]), float(row[1])
    if (chrom != chromosome or raw_chrom != chrom or strand not in {"+", "-"}
            or start < 1 or not start <= proximal <= end or not math.isfinite(fit) or fit < 0):
        raise ValueError(f"Invalid chromosome, coordinates, strand or fitting error: {context}")
    if tx not in mapping:
        raise ValueError(f"Unmapped transcript {tx}: {context}")
    gene = mapping[tx]
    if raw_gene not in MISSING and raw_gene != gene:
        raise ValueError(f"Conflicting gene mapping for {tx}: {context}")
    values, signature = [], []
    for value in row[4:]:
        value = value.strip()
        numeric = None if value in MISSING else float(value)
        if numeric is not None and (not math.isfinite(numeric) or not 0 <= numeric <= 1):
            raise ValueError(f"PDUI must be within [0, 1] or missing: {context}")
        values.append("NA" if numeric is None else value)
        signature.append(numeric)
    locus_id = f"{chrom}:{start}-{end}"
    return {"event_id": f"{tx}|{locus_id}|APA{proximal}", "transcript_id": tx,
            "gene_id": gene, "chr": chrom, "strand": strand, "s1": start, "s2": end,
            "proximal_site": proximal, "fit_value": row[1],
            "values": values, "signature": tuple(signature)}


def read_results(paths, chromosomes, mapping):
    samples, records, seen = None, [], set()
    for path, chromosome in zip(paths, chromosomes):
        with path.open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, [])
            current_samples = sample_names(header)
            if samples is None:
                samples = current_samples
            elif samples != current_samples:
                raise ValueError(f"Sample names or order differ across chromosome files: {path}")
            for number, row in enumerate(reader, 2):
                if not row:
                    continue
                if len(row) != len(header):
                    raise ValueError(f"Unexpected number of columns: {path}:{number}")
                record = parse_event(row, chromosome, mapping, f"{path}:{number}")
                if record["event_id"] in seen:
                    raise ValueError(f"Duplicate transcript event ID: {record['event_id']}")
                seen.add(record["event_id"])
                records.append(record)
    return samples, records


def collapse_events(records):
    groups = {}
    for record in records:
        key = tuple(record[name] for name in ("gene_id", "chr", "s1", "s2", "proximal_site"))
        if key not in groups:
            group = record.copy()
            gene, chrom, start, end, proximal = key
            group["event_id"] = f"{gene}|{chrom}:{start}-{end}|APA{proximal}"
            group["transcript_ids"] = set()
            groups[key] = group
        group = groups[key]
        if group["signature"] != record["signature"] or group["strand"] != record["strand"]:
            raise ValueError(f"Conflicting transcript records for gene event: {group['event_id']}")
        group["transcript_ids"].add(record["transcript_id"])
        if float(record["fit_value"]) < float(group["fit_value"]):
            group["fit_value"] = record["fit_value"]
    for group in groups.values():
        group["transcript_ids"] = ",".join(sorted(group["transcript_ids"]))
    return list(groups.values())


def write_table(path, header, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def prepare_matrix(result_prefix, chrom_list, tx2gene, out_dir):
    result_prefix, chrom_list, tx2gene, out_dir = map(Path, (result_prefix, chrom_list, tx2gene, out_dir))
    chromosomes = [line.strip() for line in chrom_list.read_text().splitlines() if line.strip()]
    if (not chromosomes or len(chromosomes) != len(set(chromosomes))
            or any("/" in chrom or any(c.isspace() for c in chrom) for chrom in chromosomes)):
        raise ValueError("Chromosome list must contain distinct chromosome names, one per line")
    paths = [Path(f"{result_prefix}_{chrom}") / f"dapars2_result_temp.{chrom}.txt" for chrom in chromosomes]
    inputs = [chrom_list, tx2gene, *paths]
    names = ["apa_transcript_events.tsv", "apa_gene_events.tsv", "apa_pdui.tsv", "apa_locations.tsv"]
    outputs = [out_dir / name for name in names]
    for index, output in enumerate(outputs):
        for other in inputs + outputs[index + 1:]:
            if (output.resolve() == other.resolve()
                    or (output.exists() and other.exists() and output.samefile(other))):
                raise ValueError("Output files must be distinct and must not overwrite inputs")
    mapping = load_mapping(tx2gene)
    samples, records = read_results(paths, chromosomes, mapping)
    genes = collapse_events(records)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_table(outputs[0], EVENT_FIELDS + samples,
                ([row[field] for field in EVENT_FIELDS] + row["values"] for row in records))
    write_table(outputs[1], GENE_FIELDS + samples,
                ([row[field] for field in GENE_FIELDS] + row["values"] for row in genes))
    write_table(outputs[2], ["event_id", *samples],
                ([row["event_id"], *row["values"]] for row in genes))
    write_table(outputs[3], ["geneid", "chr", "s1", "s2"],
                ([row["event_id"], row["chr"], row["s1"], row["s2"]] for row in genes))
    print(f"Prepared {len(records)} transcript events, {len(genes)} gene events and {len(samples)} samples")
    for output in outputs:
        print(f"Saved {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-prefix", required=True, help="Chromosome result directory prefix")
    parser.add_argument("--chrom-list", required=True, help="Chromosome list used by DaPars2")
    parser.add_argument("--tx2gene", required=True, help="Transcript-to-gene mapping TSV")
    parser.add_argument("--out-dir", required=True, help="Output matrix directory")
    args = parser.parse_args()
    prepare_matrix(args.result_prefix, args.chrom_list, args.tx2gene, args.out_dir)
