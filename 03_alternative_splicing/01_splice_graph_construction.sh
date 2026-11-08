# splice graph construction workflow (MAJIQ v3.0.22)

# 1. prepare annotation splice graph
GFF3="/path/to/annotation.gff3"  # gene annotation matching the BAM reference genome
OUT_DIR="/path/to/majiq"
ANNOTATION_SG="$OUT_DIR/annotations/sg.zarr"

mkdir -p "$OUT_DIR/annotations"

majiq-v3 gff3 \
    "$GFF3" \
    "$ANNOTATION_SG"

# 2. extract junction and intron coverage
BAM_DIR="/path/to/bam"
SJ_DIR="$OUT_DIR/sj"
SJ_THREADS=4

mkdir -p "$SJ_DIR"

for BAM in "$BAM_DIR"/*.sorted.bam; do
    SAMPLE=$(basename "$BAM" .sorted.bam)

    majiq-v3 sj \
        --nthreads "$SJ_THREADS" \
        "$BAM" \
        "$ANNOTATION_SG" \
        "$SJ_DIR/${SAMPLE}.sj"
done

# 3. build cohort splice graph
BUILD_SG="$OUT_DIR/build/sg.zarr"

mkdir -p "$OUT_DIR/build"

majiq-v3 build \
    "$ANNOTATION_SG" \
    "$BUILD_SG" \
    --sjs "$SJ_DIR"/*.sj
