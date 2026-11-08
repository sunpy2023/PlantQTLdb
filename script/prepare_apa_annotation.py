"""prepare transcript BED12 and transcript-to-gene mapping for APA analysis.

usage: python3 prepare_apa_annotation.py --input annotation.gtf --format gtf \
    --out-bed12 transcripts.bed12 --out-tx2gene tx2gene.tsv
supports GTF/GFF3, optionally gzipped, without requiring transcript/exon row order.
BED12 blocks use sorted, distinct exons and zero-based, half-open coordinates.
coding boundaries are not inferred; thickStart/thickEnd span the transcript.
the mapping has transcript_id, gene_id and gene_name columns with a header.
when no gene identifier is provided, the transcript ID is used as the fallback.
"""

import argparse
import csv
import gzip
from pathlib import Path
import re
from urllib.parse import unquote


def read_models(path, annotation_format):
    models = {}
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Expected nine annotation columns at line {number}")
            chrom, _, feature, start, end, _, strand, _, text = fields
            kind = feature.lower()
            is_transcript = (kind == "transcript" or kind.endswith("rna")
                             or kind.endswith("_transcript"))
            if kind != "exon" and not is_transcript:
                continue
            start, end = int(start), int(end)
            if start < 1 or end < start or strand not in {"+", "-"}:
                raise ValueError(f"Invalid coordinates or strand at line {number}")
            if annotation_format == "gtf":
                attrs = dict(re.findall(r'([^\s;]+)\s+"([^"]*)"', text))
                transcripts = [attrs.get("transcript_id", "")]
                gene = attrs.get("gene_id", "")
                name = attrs.get("gene_name", "")
            else:
                attrs = dict(part.split("=", 1) for part in text.strip().split(";") if "=" in part)
                parents = attrs.get("Parent", "")
                if kind == "exon":
                    transcripts = [unquote(value) for value in parents.split(",") if value]
                    gene = unquote(attrs.get("gene_id") or attrs.get("gene", ""))
                    name = unquote(attrs.get("gene_name", ""))
                else:
                    transcripts = [unquote(attrs.get("ID", ""))]
                    if "," in parents:
                        raise ValueError(f"Transcript has multiple parent genes at line {number}")
                    gene = unquote(parents or attrs.get("gene_id") or attrs.get("gene", ""))
                    name = unquote(attrs.get("Name", ""))
            if not transcripts or any(not tx or any(c.isspace() for c in tx) for tx in transcripts):
                raise ValueError(f"Missing or invalid transcript ID at line {number}")
            if any(c.isspace() for c in gene):
                raise ValueError(f"Invalid gene ID at line {number}")
            for tx in transcripts:
                model = models.setdefault(
                    tx, {"chrom": chrom, "strand": strand, "gene": "", "name": "", "exons": set()}
                )
                if (model["chrom"], model["strand"]) != (chrom, strand):
                    raise ValueError(f"Conflicting chromosome or strand for {tx}")
                if gene and model["gene"] and gene != model["gene"]:
                    raise ValueError(f"Conflicting gene IDs for {tx}")
                if gene:
                    model["gene"] = gene
                if name:
                    model["name"] = name
                if kind == "exon":
                    model["exons"].add((start, end))
    return models


def prepare_annotation(input_path, bed_path, mapping_path, annotation_format=None):
    input_path, bed_path, mapping_path = map(Path, (input_path, bed_path, mapping_path))
    paths = [input_path, bed_path, mapping_path]
    for index, path in enumerate(paths):
        for other in paths[index + 1:]:
            if (path.resolve() == other.resolve()
                    or (path.exists() and other.exists() and path.samefile(other))):
                raise ValueError("Input, BED12 and mapping must be different files")
    if annotation_format is None:
        name = input_path.name.lower()
        annotation_format = "gtf" if name.endswith((".gtf", ".gtf.gz")) else "gff3"
    models = read_models(input_path, annotation_format)
    bed_rows, mapping_rows = [], []
    for tx, model in models.items():
        exons = sorted(model["exons"])
        if not exons:
            continue
        if any(left[1] >= right[0] for left, right in zip(exons, exons[1:])):
            raise ValueError(f"Overlapping exons within transcript {tx}")
        start, end = exons[0][0] - 1, exons[-1][1]
        sizes = ",".join(str(right - left + 1) for left, right in exons) + ","
        offsets = ",".join(str(left - 1 - start) for left, _ in exons) + ","
        bed_rows.append([model["chrom"], start, end, tx, 0, model["strand"],
                         start, end, 0, len(exons), sizes, offsets])
        gene = model["gene"] or tx
        mapping_rows.append([tx, gene, model["name"] or gene])
    if not bed_rows:
        raise ValueError("No transcript exon models found in the annotation")
    bed_rows.sort(key=lambda row: (row[0], row[1], row[2], row[3]))
    mapping_rows.sort(key=lambda row: row[0])
    for path, rows in ((bed_path, bed_rows), (mapping_path, mapping_rows)):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            if path == mapping_path:
                writer.writerow(["transcript_id", "gene_id", "gene_name"])
            writer.writerows(rows)
    print(f"Prepared {len(bed_rows)} transcripts: {bed_path}; {mapping_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Reference GTF/GFF3 annotation")
    parser.add_argument("--format", choices=["gtf", "gff3"], help="Annotation format (otherwise inferred)")
    parser.add_argument("--out-bed12", required=True, help="Output transcript BED12")
    parser.add_argument("--out-tx2gene", required=True, help="Output transcript-to-gene table")
    args = parser.parse_args()
    prepare_annotation(args.input, args.out_bed12, args.out_tx2gene, args.format)
