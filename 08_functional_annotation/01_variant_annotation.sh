# variant annotation workflow

INPUT_VCF="/path/to/genotype_qc/genotype.qc.vcf.gz"
CHR_MAP="/path/to/chromosome_map.tsv"  # reference ID<TAB>analysis ID; same mapping as genotype QC, no header
REF="/path/to/reference.fa"  # reference genome FASTA
GFF="/path/to/annotation.gff3"  # gene annotation matching the reference genome
CDS="/path/to/cds.fa"  # transcript IDs matching the gene annotation
PROTEIN="/path/to/protein.fa"  # transcript IDs matching the gene annotation
GENOME="reference_build"  # SnpEff database identifier
DB_DIR="/path/to/reference/snpeff"  # database directory for this reference and annotation
OUT_DIR="/path/to/functional_annotation"
SNPEFF_JAR="/path/to/snpEff/snpEff.jar"
BCFTOOLS_THREADS=12

CONFIG="$DB_DIR/snpEff.config"
GENOME_DIR="$DB_DIR/data/$GENOME"
PREFIX="$OUT_DIR/genotype.annotated"
REVERSE_MAP="$OUT_DIR/chromosome_map.reverse.tsv"
REFERENCE_VCF="$PREFIX.reference.vcf.gz"

mkdir -p "$GENOME_DIR" "$OUT_DIR"

# 1. build the annotation database (SnpEff v5.4c; once per reference and annotation)
# reuse a matching database by running steps 2-4
# use genes.gtf and -gtf22 for GTF input
cp "$REF" "$GENOME_DIR/sequences.fa"
cp "$GFF" "$GENOME_DIR/genes.gff"
cp "$CDS" "$GENOME_DIR/cds.fa"
cp "$PROTEIN" "$GENOME_DIR/protein.fa"

printf 'data.dir = ./data\n%s.genome : %s\n' "$GENOME" "$GENOME" > "$CONFIG"

java -Xmx24g -jar "$SNPEFF_JAR" build \
    -gff3 -v -noLog \
    -c "$CONFIG" \
    "$GENOME" \
    > "$DB_DIR/build.log" 2>&1

# 2. restore reference chromosome names and check alleles (BCFtools v1.21)
# the mapping must be one-to-one and use the same reference assembly
awk 'BEGIN { OFS="\t" }
    NF != 2 || reference[$1]++ || analysis[$2]++ {
        print "Invalid or non-unique chromosome mapping at line " NR > "/dev/stderr"
        exit 1
    }
    { print $2, $1 }
    END { if (NR == 0) exit 1 }
' "$CHR_MAP" > "$REVERSE_MAP" || exit 1

bcftools annotate \
    --rename-chrs "$REVERSE_MAP" \
    --threads "$BCFTOOLS_THREADS" \
    -Oz -o "$REFERENCE_VCF" \
    "$INPUT_VCF" || exit 1

# reference FASTA and SnpEff database must use the restored chromosome names
bcftools norm \
    --check-ref e \
    -f "$REF" \
    -Ou "$REFERENCE_VCF" \
    > /dev/null 2> "$PREFIX.reference_check.log" || exit 1

# 3. annotate SNPs and indels
java -Xmx24g -jar "$SNPEFF_JAR" ann \
    -v -noLog -nodownload \
    -c "$CONFIG" \
    -stats "$PREFIX.summary.html" \
    -csvStats "$PREFIX.summary.csv" \
    "$GENOME" "$REFERENCE_VCF" \
    > "$PREFIX.reference.annotated.vcf" 2> "$PREFIX.log" || exit 1

# 4. restore analysis chromosome names, compress and index
# retain the variant IDs assigned during genotype QC
bcftools annotate \
    --rename-chrs "$CHR_MAP" \
    --threads "$BCFTOOLS_THREADS" \
    -Oz -o "$PREFIX.vcf.gz" \
    "$PREFIX.reference.annotated.vcf" || exit 1

bcftools index --threads "$BCFTOOLS_THREADS" "$PREFIX.vcf.gz"
