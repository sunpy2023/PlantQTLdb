# alternative polyadenylation preprocessing workflow (BEDTools v2.30.0; SAMtools v1.21)

# 1. prepare transcript annotation
GTF="/path/to/annotation.gtf"  # gene annotation matching the BAM reference genome
OUT_DIR="/path/to/apa"
ANNOTATION_DIR="$OUT_DIR/annotation"
TRANSCRIPT_BED="$ANNOTATION_DIR/transcripts.bed12"
TX2GENE="$ANNOTATION_DIR/tx2gene.tsv"  # transcript_id<TAB>gene_id<TAB>gene_name, with header
PREPARE_ANNOTATION_SCRIPT="script/prepare_apa_annotation.py"

mkdir -p "$ANNOTATION_DIR"

python3 "$PREPARE_ANNOTATION_SCRIPT" \
    --input "$GTF" \
    --format gtf \
    --out-bed12 "$TRANSCRIPT_BED" \
    --out-tx2gene "$TX2GENE"

# 2. prepare APA analysis regions
APA_REGIONS="$ANNOTATION_DIR/apa_regions.bed"
CHROM_LIST="$ANNOTATION_DIR/chromosomes.txt"
EXTRACT_REGIONS_SCRIPT="script/extract_apa_regions.py"

python3 "$EXTRACT_REGIONS_SCRIPT" \
    -b "$TRANSCRIPT_BED" \
    -s "$TX2GENE" \
    -o "$APA_REGIONS"

cut -f1 "$APA_REGIONS" |
    LC_ALL=C sort -u > "$CHROM_LIST"

# 3. generate bedGraph coverage files
BAM_DIR="/path/to/bam"
WIG_DIR="$OUT_DIR/wig"

mkdir -p "$WIG_DIR"

for BAM in "$BAM_DIR"/*.sorted.bam; do
    SAMPLE=$(basename "$BAM" .sorted.bam)

    bedtools genomecov \
        -ibam "$BAM" \
        -bg \
        -split |
        LC_ALL=C sort -k1,1 -k2,2n > "$WIG_DIR/${SAMPLE}.wig"
done

# 4. calculate mapped-read counts
DEPTH_FILE="$OUT_DIR/wig_depth.tsv"  # WIG path<TAB>mapped-read count, no header

for BAM in "$BAM_DIR"/*.sorted.bam; do
    SAMPLE=$(basename "$BAM" .sorted.bam)

    MAPPED_READS=$(samtools idxstats "$BAM" |
        awk '{sum += $3} END {printf "%.0f\n", sum}')

    printf '%s\t%s\n' \
        "$WIG_DIR/${SAMPLE}.wig" \
        "$MAPPED_READS"
done > "$DEPTH_FILE"

LC_ALL=C sort -k1,1 -o "$DEPTH_FILE" "$DEPTH_FILE"
