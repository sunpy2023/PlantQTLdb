"""collapse RNA and WGS/capture sample columns into one diploid GT per material.

read a coordinate-sorted, biallelic SNP VCF from stdin and write VCF to stdout.
sample names must end in __RNA and __WGS or __CAPTURE. Missing calls in one
source are filled from the other; phased calls are compared without phase.
the output contains unphased GT only; source INFO values are cleared.

for --conflict rna, RNA records must carry INFO/RNA_SOURCE=1. At positions
with different ALT alleles, retain the RNA record and discard the WGS record.
for --conflict missing, discard different-ALT positions upstream, as in the
S028 bcftools merge -m snps / view -m2 -M2 commands.

usage examples are in 01_genotype/03_genotype_merging.sh.
run the built-in check with: python3 vcf_assay_consensus.py --self-test
"""

import argparse
from itertools import groupby
import sys


def canonical_gt(value):
    alleles = value.replace("|", "/").split("/")
    if "." in alleles:
        return None
    if len(alleles) != 2 or any(allele not in {"0", "1"} for allele in alleles):
        raise ValueError(f"Expected a diploid biallelic GT, found {value!r}")
    return "/".join(sorted(alleles))


def consensus_gt(rna_gt, other_gt, conflict):
    if rna_gt is None:
        return other_gt or "./."
    if other_gt is None or rna_gt == other_gt or conflict == "rna":
        return rna_gt
    return "./."


def merge_records(source, output, materials, other_source, conflict):
    if not materials or len(materials) != len(set(materials)):
        raise ValueError("Material list must be nonempty and contain unique IDs")
    if any(len(material.split()) != 1 for material in materials):
        raise ValueError("Expected one material ID per line, without whitespace")

    rna_flag_header = False
    for line in source:
        if line.startswith("##INFO=<ID=RNA_SOURCE,"):
            rna_flag_header = True
        elif line.startswith("##"):
            output.write(line)
        elif line.startswith("#CHROM"):
            names = line.rstrip("\n").split("\t")[9:]
            break
        else:
            raise ValueError("VCF #CHROM header not found before records")
    else:
        raise ValueError("VCF #CHROM header not found")

    if len(names) != len(set(names)):
        raise ValueError("Duplicate sample columns in merged VCF")
    suffixes = ("__RNA", f"__{other_source}")
    if any(not name.endswith(suffixes) for name in names):
        raise ValueError(f"Sample names must end in {suffixes}; check rename maps")
    if conflict == "rna" and not rna_flag_header:
        raise ValueError("RNA priority requires the INFO/RNA_SOURCE header and site tags")
    lookup = {name: index for index, name in enumerate(names)}
    indices = [(lookup.get(material + suffixes[0]), lookup.get(material + suffixes[1]))
               for material in materials]
    if any(rna is None and other is None for rna, other in indices):
        raise ValueError("A requested material is absent from both sources")
    output.write(f'##assayConsensus="RNA/{other_source}; discordant calls: {conflict}"\n')
    output.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                            "FILTER", "INFO", "FORMAT"] + materials) + "\n")

    def records():
        for line in source:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 + len(names):
                raise ValueError("VCF record column count differs from header")
            if int(fields[1]) < 1 or len(fields[3]) != 1 or len(fields[4]) != 1:
                raise ValueError(f"Expected a biallelic SNP at {fields[0]}:{fields[1]}")
            yield fields

    n_records = n_discordant = n_dropped = 0
    for coordinate, group in groupby(records(), key=lambda row: (row[0], row[1])):
        rows = list(group)
        if conflict == "rna":
            rna_rows = [row for row in rows if "RNA_SOURCE=1" in row[7].split(";")]
            if rna_rows:
                n_dropped += len(rows) - len(rna_rows)
                rows = rna_rows
        if len(rows) != 1:
            raise ValueError(f"Ambiguous SNP records at {':'.join(coordinate)}")
        fields = rows[0]
        fmt = fields[8].split(":")
        if "GT" not in fmt:
            raise ValueError("Input record lacks FORMAT/GT")
        gt_index = fmt.index("GT")
        sample_values = fields[9:]

        def get_gt(index):
            if index is None:
                return None
            values = sample_values[index].split(":")
            if values == ["."]:
                return None
            if gt_index >= len(values):
                raise ValueError("Sample value lacks the declared FORMAT/GT field")
            return canonical_gt(values[gt_index])

        output_gt = []
        for rna_index, other_index in indices:
            rna_gt, other_gt = get_gt(rna_index), get_gt(other_index)
            if rna_gt is not None and other_gt is not None and rna_gt != other_gt:
                n_discordant += 1
            output_gt.append(consensus_gt(rna_gt, other_gt, conflict))
        fields[7:9] = [".", "GT"]
        output.write("\t".join(fields[:9] + output_gt) + "\n")
        n_records += 1
    return n_records, n_discordant, n_dropped


def self_test():
    from io import StringIO

    header = (
        '##fileformat=VCFv4.2\n'
        '##INFO=<ID=RNA_SOURCE,Number=1,Type=Integer,Description="RNA site">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
        'a__RNA\ta__WGS\tb__RNA\tb__WGS\tc__WGS\n'
    )
    data = (
        'chr1\t10\t.\tA\tC\t.\tPASS\tRNA_SOURCE=1\tDP:GT\t'
        '5:0|1\t5:1/0\t5:0/0\t5:1/1\t5:1/1\n'
        'chr1\t20\t.\tA\tG\t.\tPASS\tRNA_SOURCE=1\tGT\t'
        '0/1\t./.\t./.\t./.\t./.\n'
        'chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t'
        './.\t1/1\t./.\t0/1\t0/0\n'
    )
    output = StringIO()
    counts = merge_records(StringIO(header + data), output, ["a", "b", "c"], "WGS", "rna")
    rows = [line.split("\t") for line in output.getvalue().splitlines() if not line.startswith("#")]
    assert counts == (2, 1, 1)
    assert rows[0][9:] == ["0/1", "0/0", "1/1"]
    assert rows[1][4] == "G" and rows[1][9:] == ["0/1", "./.", "./."]
    output = StringIO()
    source = (header + data.splitlines(keepends=True)[0]).replace("__WGS", "__CAPTURE")
    merge_records(StringIO(source), output, ["c", "b", "a"], "CAPTURE", "missing")
    assert output.getvalue().splitlines()[-1].split("\t")[9:] == ["1/1", "./.", "0/1"]
    assert canonical_gt("1|0") == "0/1" and canonical_gt("0/.") is None
    assert consensus_gt(None, None, "missing") == "./."
    assert consensus_gt(None, "0/0", "rna") == "0/0"
    try:
        merge_records(StringIO(header + data), StringIO(), ["absent"], "WGS", "rna")
    except ValueError:
        pass
    else:
        raise AssertionError("Missing material was silently accepted")
    print("Self-test passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--materials", help="One material ID per line in output order, no header")
    parser.add_argument("--other-source", choices=("WGS", "CAPTURE"), default="CAPTURE")
    parser.add_argument("--conflict", choices=("rna", "missing"), default="missing")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if not args.materials:
        parser.error("--materials is required")
    with open(args.materials, encoding="utf-8") as handle:
        materials = [line.strip() for line in handle if line.strip()]
    counts = merge_records(sys.stdin, sys.stdout, materials, args.other_source, args.conflict)
    sys.stderr.write(f"records={counts[0]}\tdiscordant={counts[1]}\tdifferent_alt_dropped={counts[2]}\n")


if __name__ == "__main__":
    main()
