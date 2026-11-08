"""extract DaPars2 analysis regions from transcript BED12 and a gene mapping.

usage: python3 extract_apa_regions.py -b transcripts.bed12 -s tx2gene.tsv -o apa_regions.bed
requires bedtools. Terminal exons are used as candidate APA regions, not only
annotated untranslated bases. The shifted start coordinates, coordinate-based
deduplication and opposite-strand subtraction follow the existing DaPars2 workflow.
after subtraction, keep the leftmost fragment on + and the rightmost on -.
output names are transcript_id|gene_id|chromosome|strand; transcript IDs may
contain pipes. Contigs containing underscores are retained for plant genomes.
"""

import argparse
from pathlib import Path
import subprocess
import tempfile


def load_mapping(path):
    mapping = {}
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if fields[0] == "transcript_id":
                continue
            if len(fields) < 2 or not fields[0] or not fields[1]:
                raise ValueError(f"Invalid transcript-to-gene mapping at line {number}")
            tx, gene = fields[:2]
            if any(c.isspace() for c in tx + gene) or "|" in gene:
                raise ValueError(f"Unsupported transcript/gene ID at line {number}")
            if tx in mapping and mapping[tx] != gene:
                raise ValueError(f"Conflicting gene IDs for transcript {tx}")
            mapping[tx] = gene
    if not mapping:
        raise ValueError("Empty transcript-to-gene mapping")
    return mapping


def write_candidates(bed_path, mapping, output_path):
    seen_regions, seen_transcripts = set(), set()
    with bed_path.open() as source, output_path.open("w") as output:
        for number, line in enumerate(source, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Expected BED12 fields at line {number}")
            chrom, start, end, tx = fields[:4]
            start, end = int(start), int(end)
            strand, count = fields[5], int(fields[9])
            sizes = [int(value) for value in fields[10].rstrip(",").split(",")]
            offsets = [int(value) for value in fields[11].rstrip(",").split(",")]
            if (start < 0 or end <= start or strand not in {"+", "-"}
                    or count < 1 or count != len(sizes) or count != len(offsets)
                    or any(size <= 0 or offset < 0 or offset + size > end - start
                           for size, offset in zip(sizes, offsets))
                    or any(left + size > right for left, size, right
                           in zip(offsets, sizes, offsets[1:]))):
                raise ValueError(f"Invalid BED12 coordinates or blocks at line {number}")
            if tx in seen_transcripts or tx not in mapping:
                raise ValueError(f"Duplicate or unmapped transcript ID: {tx}")
            seen_transcripts.add(tx)
            if any(c.isspace() for c in chrom + tx) or "|" in chrom:
                raise ValueError(f"Unsupported chromosome/transcript ID at line {number}")
            # preserve the shifted start coordinates used by the DaPars2 helper
            left = start + offsets[-1] + 1 if strand == "+" else start + 1
            right = end if strand == "+" else start + sizes[0]
            key = (chrom, left, right, strand)
            if key in seen_regions:
                continue
            seen_regions.add(key)
            name = f"{tx}|{mapping[tx]}|{chrom}|{strand}"
            output.write(f"{chrom}\t{left}\t{right}\t{name}\t0\t{strand}\n")
    if not seen_regions:
        raise ValueError("No candidate APA regions found")
    return len(seen_regions)


def extract_regions(bed_path, mapping_path, output_path):
    bed_path, mapping_path, output_path = map(Path, (bed_path, mapping_path, output_path))
    for source in (bed_path, mapping_path):
        if (output_path.resolve() == source.resolve()
                or (output_path.exists() and output_path.samefile(source))):
            raise ValueError("Output must not overwrite an input file")
    mapping = load_mapping(mapping_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".apa_regions.", dir=output_path.parent) as folder:
        folder = Path(folder)
        candidates = folder / "candidates.bed"
        subtracted = folder / "subtracted.bed"
        final = folder / "final.bed"
        candidate_count = write_candidates(bed_path, mapping, candidates)
        with subtracted.open("w") as handle:
            subprocess.run(["bedtools", "subtract", "-a", str(candidates),
                            "-b", str(candidates), "-S"], stdout=handle, check=True)
        selected = {}
        with subtracted.open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                tx = fields[3].rsplit("|", 3)[0]
                previous = selected.get(tx)
                if (previous is None
                        or (fields[5] == "+" and int(fields[1]) < int(previous[1]))
                        or (fields[5] == "-" and int(fields[2]) > int(previous[2]))):
                    selected[tx] = fields
        if not selected:
            raise ValueError("No APA regions remain after opposite-strand subtraction")
        with final.open("w") as handle:
            for fields in selected.values():
                handle.write("\t".join(fields) + "\n")
        final.replace(output_path)
    print(f"Prepared {len(selected)} APA regions from {candidate_count} candidates: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-b", "--bed", required=True, help="Transcript BED12")
    parser.add_argument("-s", "--symbol", required=True, help="Transcript-to-gene mapping TSV")
    parser.add_argument("-o", "--output", required=True, help="Output DaPars2 annotation")
    args = parser.parse_args()
    extract_regions(args.bed, args.symbol, args.output)
