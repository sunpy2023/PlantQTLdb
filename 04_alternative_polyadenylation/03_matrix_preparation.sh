# APA matrix preparation workflow

OUT_DIR="/path/to/apa"

# 1. merge chromosome results and prepare APA matrices
RESULT_PREFIX="$OUT_DIR/results"
CHROM_LIST="$OUT_DIR/annotation/chromosomes.txt"
TX2GENE="$OUT_DIR/annotation/tx2gene.tsv"
MATRIX_DIR="$OUT_DIR/matrices"
PREPARE_MATRIX_SCRIPT="script/prepare_apa_matrix.py"

mkdir -p "$MATRIX_DIR"

python3 "$PREPARE_MATRIX_SCRIPT" \
    --result-prefix "$RESULT_PREFIX" \
    --chrom-list "$CHROM_LIST" \
    --tx2gene "$TX2GENE" \
    --out-dir "$MATRIX_DIR"
