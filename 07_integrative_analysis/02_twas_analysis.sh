# TWAS workflow

MOLECULAR_TYPE="expression"  # expression (FPKM), splicing (PSI) or apa (PDUI)
PHENO_DIR="/path/to/molqtl/phenotypes"
MOLECULAR_MATRIX="$PHENO_DIR/${MOLECULAR_TYPE}.normalized.tsv"  # rank-normalized matrix
LOCATIONS="$PHENO_DIR/${MOLECULAR_TYPE}.locations.tsv"  # feature coordinates
GENO_PREFIX="/path/to/genotype.qc"  # PLINK BED/BIM/FAM prefix; chromosomes 1..CHR_COUNT
TRAIT="trait_name"
ASSOC_FILE="/path/to/gwas/${TRAIT}.assoc.txt"  # full GEMMA results
OUT_DIR="/path/to/integrative_analysis/twas/$MOLECULAR_TYPE"
LDREF_DIR="/path/to/integrative_analysis/twas/ldref"  # same population; shared across molecular types
FUSION_DIR="/path/to/fusion_twas"
GCTA="/path/to/gcta64"
CHR_COUNT=20  # nuclear chromosome count matching genotype QC
CIS_WINDOW=200000  # bp on each side of the molecular feature interval
MODELS="top1,lasso,enet"  # top1, LASSO and elastic net
CROSSVAL=5  # cross-validation folds
MIN_VARIANTS=5  # minimum number of cis variants for training
PREPARE_SCRIPT="script/prepare_twas_inputs.py"
TWAS_SCRIPT="script/run_twas.py"
SUMMARIZE_SCRIPT="script/summarize_twas_results.py"

mkdir -p "$OUT_DIR/runs" "$OUT_DIR/results"
# retain inputs, prediction weights and logs separately for each run
RUN_DIR=$(mktemp -d "$OUT_DIR/runs/${TRAIT}.XXXXXX")

# 1. prepare TWAS inputs
python3 "$PREPARE_SCRIPT" \
    --matrix "$MOLECULAR_MATRIX" \
    --locations "$LOCATIONS" \
    --fam "${GENO_PREFIX}.fam" \
    --gwas "$ASSOC_FILE" \
    --trait "$TRAIT" \
    --out-dir "$RUN_DIR" \
    > "$RUN_DIR/prepare.log" 2>&1

# 2. train prediction models and test TWAS associations (FUSION; PLINK 1.9; GCTA v1.95.1)
python3 "$TWAS_SCRIPT" \
    --inputs-dir "$RUN_DIR" \
    --geno-prefix "$GENO_PREFIX" \
    --fusion-dir "$FUSION_DIR" \
    --gcta "$GCTA" \
    --trait "$TRAIT" \
    --out-dir "$RUN_DIR" \
    --ldref-dir "$LDREF_DIR" \
    --chr-count "$CHR_COUNT" \
    --cis-window "$CIS_WINDOW" \
    --models "$MODELS" \
    --crossval "$CROSSVAL" \
    --min-variants "$MIN_VARIANTS" \
    > "$RUN_DIR/fusion.log" 2>&1

# 3. merge results and extract significant associations (Bonferroni correction)
python3 "$SUMMARIZE_SCRIPT" \
    --results-dir "$RUN_DIR" \
    --weights "$RUN_DIR/weights.pos" \
    --trait "$TRAIT" \
    --alpha 0.05 \
    --out-prefix "$OUT_DIR/results/$TRAIT"
