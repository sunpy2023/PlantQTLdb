# splicing event classification workflow (VOILA v3.0.22)

OUT_DIR="/path/to/majiq"
BUILD_SG="$OUT_DIR/build/sg.zarr"
PSI_DIR="$OUT_DIR/psi"
SGC_DIR="$OUT_DIR/build/sgc_per_sample"

# 1. classify splicing events
EVENT_DIR="$OUT_DIR/modulizer/ALL"
VOILA_THREADS=12
DECOMPLEXIFY_READS=10  # minimum read support for decomplexification

mkdir -p "$EVENT_DIR"

voila modulize \
    -d "$EVENT_DIR" \
    "$BUILD_SG" \
    "$PSI_DIR"/*.psicov \
    "$SGC_DIR"/*.sgc \
    --show-read-counts \
    --decomplexify-reads-threshold "$DECOMPLEXIFY_READS" \
    -j "$VOILA_THREADS"

# 2. prepare the junction PSI matrix
MATRIX_DIR="$OUT_DIR/matrices"
PREPARE_PSI_SCRIPT="script/prepare_junction_psi.py"

mkdir -p "$MATRIX_DIR"

python3 "$PREPARE_PSI_SCRIPT" \
    -i "$EVENT_DIR/junctions.tsv" \
    -o "$MATRIX_DIR/junction_psi.tsv"
