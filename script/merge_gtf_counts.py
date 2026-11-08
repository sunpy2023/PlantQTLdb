"""estimate gene and transcript counts from StringTie -e GTF coverage.

usage: python3 merge_gtf_counts.py sample_gtf.tsv gene_counts.csv transcript_counts.csv
the sample/GTF list has two tab-separated columns without a header.
READ_LEN is the average read length in bp (default: 75).
each transcript row must precede its exon rows and contain gene_id,
transcript_id and cov attributes. Use the same reference annotation for all samples.
counts are ceil(coverage * total exon length / read length), summed per gene.
columns follow sample-list order; absent genes/transcripts receive 0.
"""

import csv
import math
import os
from pathlib import Path
import re
import sys


def process_gtf(gtf_path, read_len):
    transcripts = {}
    current_tid = None
    with open(gtf_path) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("#") or not line.strip():
                continue
            location = f"{gtf_path}:{line_number}"
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Expected nine GTF columns: {location}")
            if fields[2] not in {"transcript", "exon"}:
                continue
            attributes = dict(re.findall(r'(\w+)\s+"([^"]*)"', fields[8]))
            tid = attributes.get("transcript_id")
            gene = attributes.get("gene_id")
            if not tid or not gene:
                raise ValueError(f"Missing transcript_id or gene_id: {location}")
            if fields[2] == "transcript":
                if tid in transcripts:
                    raise ValueError(f"Duplicate transcript_id {tid!r}: {location}")
                try:
                    coverage = float(attributes["cov"])
                except (KeyError, ValueError):
                    raise ValueError(f"Missing or invalid transcript coverage: {location}") from None
                if not math.isfinite(coverage) or coverage < 0:
                    raise ValueError(f"Expected finite, nonnegative coverage: {location}")
                transcripts[tid] = {"gene": gene, "coverage": coverage, "length": 0}
                current_tid = tid
            else:
                if tid != current_tid or transcripts[tid]["gene"] != gene:
                    raise ValueError(f"Exon does not match its preceding transcript: {location}")
                start, end = int(fields[3]), int(fields[4])
                if start < 1 or end < start:
                    raise ValueError(f"Invalid exon coordinates: {location}")
                transcripts[tid]["length"] += end - start + 1

    counts = {}
    for tid, transcript in transcripts.items():
        if transcript["length"] == 0:
            raise ValueError(f"Transcript {tid!r} has no exons: {gtf_path}")
        counts[tid] = math.ceil(transcript["coverage"] * transcript["length"] / read_len)
    return {tid: transcript["gene"] for tid, transcript in transcripts.items()}, counts


def write_matrix(path, data, samples, id_column):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([id_column, *samples])
        for feature in sorted(data):
            writer.writerow([feature] + [data[feature].get(sample, 0) for sample in samples])


def merge_counts(sample_list, gene_output, transcript_output, read_len):
    if not math.isfinite(read_len) or read_len <= 0:
        raise ValueError("READ_LEN must be finite and positive")
    samples = []
    with open(sample_list) as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 2 or not fields[0] or not fields[1]:
                raise ValueError("Expected sample ID<TAB>GTF path, without a header")
            samples.append(fields)
    sample_names = [sample for sample, _ in samples]
    if not samples or len(sample_names) != len(set(sample_names)):
        raise ValueError("Sample list must be nonempty and contain unique IDs")

    gene_counts, transcript_counts, transcript_genes = {}, {}, {}
    for sample, gtf_path in samples:
        genes, counts = process_gtf(gtf_path, read_len)
        for tid, count in counts.items():
            gene = genes[tid]
            if tid in transcript_genes and transcript_genes[tid] != gene:
                raise ValueError(f"Conflicting gene_id for transcript {tid!r}: {gtf_path}")
            transcript_genes[tid] = gene
            transcript_counts.setdefault(tid, {})[sample] = count
            values = gene_counts.setdefault(gene, {})
            values[sample] = values.get(sample, 0) + count

    if Path(gene_output).resolve() == Path(transcript_output).resolve():
        raise ValueError("Gene and transcript output paths must differ")
    write_matrix(gene_output, gene_counts, sample_names, "gene_id")
    write_matrix(transcript_output, transcript_counts, sample_names, "transcript_id")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    merge_counts(*sys.argv[1:], read_len=float(os.environ.get("READ_LEN", 75)))
