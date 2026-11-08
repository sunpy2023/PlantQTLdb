# EWAS workflow

MOLECULAR_TYPE="expression"  # expression (FPKM), splicing (PSI) or apa (PDUI)
PHENO_DIR="/path/to/molqtl/phenotypes"
MOLECULAR_MATRIX="$PHENO_DIR/${MOLECULAR_TYPE}.filtered.tsv"  # before normalization
FAM="/path/to/genotype_qc/genotype.qc.fam"
TRAIT="trait_name"
PHENOTYPE="/path/to/trait.txt"  # one column in FAM order, no header; NA for missing values
OUT_DIR="/path/to/integrative_analysis/ewas/$MOLECULAR_TYPE"
OSCA="osca"
P_THRESHOLD=0.001  # unadjusted P value for exploratory associations
PREPARE_SCRIPT="script/prepare_ewas_inputs.py"
SUMMARIZE_SCRIPT="script/summarize_ewas_results.py"

mkdir -p "$OUT_DIR/runs" "$OUT_DIR/results"
# retain inputs, binary omics data and logs separately for each run
RUN_DIR=$(mktemp -d "$OUT_DIR/runs/${TRAIT}.XXXXXX")
PREFIX="$RUN_DIR/$TRAIT"
EFILE="$PREFIX.omics.tsv"
PHENO_FILE="$PREFIX.phenotype.tsv"

# 1. match samples and prepare trait-specific OSCA inputs
python3 "$PREPARE_SCRIPT" \
    --matrix "$MOLECULAR_MATRIX" \
    --fam "$FAM" \
    --phenotype "$PHENOTYPE" \
    --trait "$TRAIT" \
    --out-prefix "$PREFIX" \
    > "$PREFIX.prepare.log" 2>&1

# 2. prepare binary omics data (OSCA v1.20)
"$OSCA" \
    --efile "$EFILE" \
    --make-bod \
    --out "$PREFIX" \
    > "$PREFIX.bod.log" 2>&1

# 3. MOMENT association analysis
# OSCA switches to MOA when MOMENT does not converge
"$OSCA" \
    --moment \
    --befile "$PREFIX" \
    --pheno "$PHENO_FILE" \
    --out "$PREFIX" \
    > "$PREFIX.moment.log" 2>&1

# 4. label MOMENT or MOA results and filter associations
python3 "$SUMMARIZE_SCRIPT" \
    --input "$PREFIX.moment" \
    --log "$PREFIX.moment.log" \
    --p-threshold "$P_THRESHOLD" \
    --out-prefix "$OUT_DIR/results/$TRAIT"
