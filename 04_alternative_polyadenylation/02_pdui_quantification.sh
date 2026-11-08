# PDUI quantification workflow (DaPars2 v2.1)

OUT_DIR="/path/to/apa"

# 1. prepare the DaPars2 configuration
APA_REGIONS="$OUT_DIR/annotation/apa_regions.bed"
CHROM_LIST="$OUT_DIR/annotation/chromosomes.txt"
DEPTH_FILE="$OUT_DIR/wig_depth.tsv"
CONFIG_FILE="$OUT_DIR/dapars2_config.txt"
RESULT_PREFIX="$OUT_DIR/results"  # output directories: results_<chromosome>
COVERAGE_THRESHOLD=20
DAPARS_THREADS=4

WIG_FILES=$(cut -f1 "$DEPTH_FILE" | paste -sd, -)

mkdir -p "$OUT_DIR"

cat > "$CONFIG_FILE" <<EOF
Annotated_3UTR=$APA_REGIONS
Aligned_Wig_files=$WIG_FILES
Output_directory=$RESULT_PREFIX
Output_result_file=dapars2
Coverage_threshold=$COVERAGE_THRESHOLD
Num_Threads=$DAPARS_THREADS
sequencing_depth_file=$DEPTH_FILE
EOF

# 2. infer proximal poly(A) sites and quantify PDUI
DAPARS_SCRIPT="/path/to/DaPars2/src/DaPars2_Multi_Sample_Multi_Chr.py"

python3 "$DAPARS_SCRIPT" \
    "$CONFIG_FILE" \
    "$CHROM_LIST"
