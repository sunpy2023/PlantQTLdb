# colocalization workflow

MOLECULAR_TYPE="expression"  # expression (eQTL), splicing (sQTL) or apa (apaQTL)
QTL_ROOT="/path/to/molqtl"
MOLECULAR_MATRIX="$QTL_ROOT/phenotypes/${MOLECULAR_TYPE}.normalized.tsv"  # QTL samples and per-feature counts
GENO_PREFIX="/path/to/genotype.qc"  # matched PLINK BED/BIM/FAM and VCF prefix
TRAIT="trait_name"
ASSOC_FILE="/path/to/gwas/${TRAIT}.assoc.txt"  # full GEMMA results
PHENOTYPE="/path/to/trait.txt"  # one column in FAM order; no header, missing = NA
OUT_DIR="/path/to/integrative_analysis/coloc/$MOLECULAR_TYPE"
CHR_COUNT=20  # nuclear chromosome count matching genotype QC
PREPARE_SCRIPT="script/prepare_coloc_inputs.py"
COLOC_SCRIPT="script/run_colocalization.R"

mkdir -p "$OUT_DIR/tmp" "$OUT_DIR/logs" "$OUT_DIR/results"
# create a separate working directory for each run; retain intermediate files
RUN_DIR=$(mktemp -d "$OUT_DIR/tmp/run.XXXXXX")

# 1. prepare molecular QTL and GWAS summary statistics
python3 "$PREPARE_SCRIPT" \
    --molecular-type "$MOLECULAR_TYPE" \
    --qtl-root "$QTL_ROOT" \
    --matrix "$MOLECULAR_MATRIX" \
    --vcf "${GENO_PREFIX}.vcf.gz" \
    --bim "${GENO_PREFIX}.bim" \
    --fam "${GENO_PREFIX}.fam" \
    --gwas "$ASSOC_FILE" \
    --phenotype "$PHENOTYPE" \
    --out-dir "$RUN_DIR" \
    > "$OUT_DIR/logs/${TRAIT}.prepare.log" 2>&1

# 2. select GWAS lead variants (PLINK 1.9)
N_QC=$(awk 'END { print NR }' "${GENO_PREFIX}.bim")  # variants retained after genotype QC
CLUMP_P1=$(awk -v n="$N_QC" 'BEGIN { if (n > 0) printf "%.17g", 0.05 / n }')
CLUMP_P2=1e-4
CLUMP_R2=0.2
CLUMP_KB=200  # kb on each side of the index variant

plink \
    --bfile "$GENO_PREFIX" \
    --chr-set "$CHR_COUNT" \
    --allow-extra-chr \
    --clump "$RUN_DIR/gwas.clump.tsv" \
    --clump-p1 "$CLUMP_P1" \
    --clump-p2 "$CLUMP_P2" \
    --clump-r2 "$CLUMP_R2" \
    --clump-kb "$CLUMP_KB" \
    --out "$RUN_DIR/gwas" \
    > "$OUT_DIR/logs/${TRAIT}.clump.log" 2>&1

# 3. run Bayesian colocalization and retain supported associations (coloc v5.2.3)
WINDOW_BP=200000  # bp on each side of a GWAS lead variant for feature selection
MIN_SHARED=50  # minimum number of shared variants per molecular feature
P1=1e-4
P2=1e-4
P12=1e-5
PPH4_MIN=0.5
PPH4_STRONG=0.8

Rscript "$COLOC_SCRIPT" \
    --qtl "$RUN_DIR/qtl.coloc.tsv.gz" \
    --gwas "$RUN_DIR/gwas.coloc.tsv.gz" \
    --clump-prefix "$RUN_DIR/gwas" \
    --window-bp "$WINDOW_BP" \
    --min-shared "$MIN_SHARED" \
    --p1 "$P1" \
    --p2 "$P2" \
    --p12 "$P12" \
    --pph4-min "$PPH4_MIN" \
    --pph4-strong "$PPH4_STRONG" \
    --out-prefix "$OUT_DIR/results/$TRAIT" \
    > "$OUT_DIR/logs/${TRAIT}.coloc.log" 2>&1
