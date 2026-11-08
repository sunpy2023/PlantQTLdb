# SMR workflow

MOLECULAR_TYPE="expression"  # expression (eQTL), splicing (sQTL) or apa (apaQTL)
QTL_ROOT="/path/to/molqtl"
LOCATIONS="$QTL_ROOT/phenotypes/${MOLECULAR_TYPE}.locations.tsv"  # feature coordinates
GENO_PREFIX="/path/to/genotype.qc"  # matched PLINK BED/BIM/FAM prefix
COUNTED_ALLELE="ALT"  # match QTL mapping; current defaults: expression ALT, splicing/apa REF
TRAIT="trait_name"
ASSOC_FILE="/path/to/gwas/${TRAIT}.assoc.txt"  # full GEMMA results
PHENOTYPE="/path/to/trait.txt"  # one column in FAM order; no header, missing = NA
OUT_DIR="/path/to/integrative_analysis/smr/$MOLECULAR_TYPE"
SMR="smr"
CHR_COUNT=20  # nuclear chromosome count matching genotype QC
SMR_THREADS=10
PREPARE_SCRIPT="script/prepare_smr_inputs.py"
SUMMARIZE_SCRIPT="script/summarize_smr_results.py"

mkdir -p "$OUT_DIR/tmp" "$OUT_DIR/logs" "$OUT_DIR/results"
# create a separate working directory for each run; retain intermediate files
RUN_DIR=$(mktemp -d "$OUT_DIR/tmp/run.XXXXXX")
QTL_QUERY="$RUN_DIR/qtl.query.tsv"
GWAS_SUMMARY="$RUN_DIR/${TRAIT}.gwas.ma"
BESD_PREFIX="$RUN_DIR/qtl"
RESULT_PREFIX="$RUN_DIR/$TRAIT"

# 1. prepare allele frequencies and SMR inputs (PLINK 1.9)
plink \
    --bfile "$GENO_PREFIX" \
    --chr-set "$CHR_COUNT" \
    --allow-extra-chr \
    --keep-allele-order \
    --freq \
    --out "$RUN_DIR/genotype" \
    > "$OUT_DIR/logs/${TRAIT}.frequency.log" 2>&1

python3 "$PREPARE_SCRIPT" \
    --molecular-type "$MOLECULAR_TYPE" \
    --qtl-root "$QTL_ROOT" \
    --locations "$LOCATIONS" \
    --bim "${GENO_PREFIX}.bim" \
    --fam "${GENO_PREFIX}.fam" \
    --frequency "$RUN_DIR/genotype.frq" \
    --counted-allele "$COUNTED_ALLELE" \
    --gwas "$ASSOC_FILE" \
    --phenotype "$PHENOTYPE" \
    --trait "$TRAIT" \
    --out-dir "$RUN_DIR" \
    > "$OUT_DIR/logs/${TRAIT}.prepare.log" 2>&1

# 2. create BESD files (SMR v1.4.1)
"$SMR" \
    --qfile "$QTL_QUERY" \
    --make-besd \
    --out "$BESD_PREFIX" \
    > "$OUT_DIR/logs/${TRAIT}.besd.log" 2>&1

# 3. run SMR and HEIDI tests
"$SMR" \
    --bfile "$GENO_PREFIX" \
    --gwas-summary "$GWAS_SUMMARY" \
    --beqtl-summary "$BESD_PREFIX" \
    --thread-num "$SMR_THREADS" \
    --out "$RESULT_PREFIX" \
    > "$OUT_DIR/logs/${TRAIT}.smr.log" 2>&1

# 4. extract significant SMR associations and annotate HEIDI results
python3 "$SUMMARIZE_SCRIPT" \
    --input "${RESULT_PREFIX}.smr" \
    --alpha 0.05 \
    --heidi-threshold 0.05 \
    --out-prefix "$OUT_DIR/results/$TRAIT"
