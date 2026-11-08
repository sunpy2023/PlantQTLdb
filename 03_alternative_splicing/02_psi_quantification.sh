# PSI quantification workflow (MAJIQ v3.0.22)

OUT_DIR="/path/to/majiq"
BUILD_SG="$OUT_DIR/build/sg.zarr"
SJ_DIR="$OUT_DIR/sj"

# 1. calculate splice graph coverage
SGC_DIR="$OUT_DIR/build/sgc_per_sample"

mkdir -p "$SGC_DIR"

for SJ in "$SJ_DIR"/*.sj; do
    SAMPLE=$(basename "$SJ" .sj)

    majiq-v3 sg-coverage \
        "$BUILD_SG" \
        "$SGC_DIR/${SAMPLE}.sgc" \
        "$SJ"
done

# 2. prepare coverage for PSI quantification
PSI_DIR="$OUT_DIR/psi"

mkdir -p "$PSI_DIR"

for SJ in "$SJ_DIR"/*.sj; do
    SAMPLE=$(basename "$SJ" .sj)

    majiq-v3 psi-coverage \
        "$BUILD_SG" \
        "$PSI_DIR/${SAMPLE}.psicov" \
        "$SJ"
done

# 3. quantify PSI and export per-sample tables
PSI_TSV_DIR="$OUT_DIR/psi_tsv"

mkdir -p "$PSI_TSV_DIR"

for PSICOV in "$PSI_DIR"/*.psicov; do
    SAMPLE=$(basename "$PSICOV" .psicov)

    majiq-v3 quantify \
        "$PSICOV" \
        --splicegraph "$BUILD_SG" \
        --output-tsv "$PSI_TSV_DIR/${SAMPLE}.tsv"
done
