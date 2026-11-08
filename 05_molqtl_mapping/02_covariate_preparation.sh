# molecular QTL covariate preparation workflow

GENO_PREFIX="/path/to/genotype.qc"  # PLINK BED/BIM/FAM prefix; chromosomes 1..CHR_COUNT
PCA_DIR="/path/to/genotype_pca"  # calculate once; share PCs between GWAS and molQTL
PHENO_DIR="/path/to/molqtl/phenotypes"
SAMPLE_LIST="/path/to/qtl_samples.txt"  # same sample list used in phenotype preparation
OUT_DIR="/path/to/molqtl/covariates"
CHR_COUNT=20  # nuclear chromosome count matching genotype QC
PLINK_THREADS=12
PREPARE_COVARIATES_SCRIPT="script/prepare_qtl_covariates.py"
PEER_SCRIPT="script/estimate_peer_factors.R"

mkdir -p "$PCA_DIR" "$OUT_DIR"

# 1. LD pruning (PLINK 1.9)
LD_WINDOW=50  # window size in variants
LD_STEP=5  # window step in variants
LD_R2=0.2
PRUNE_PREFIX="$PCA_DIR/genotype"

plink \
    --bfile "$GENO_PREFIX" \
    --chr-set "$CHR_COUNT" \
    --allow-extra-chr \
    --indep-pairwise "$LD_WINDOW" "$LD_STEP" "$LD_R2" \
    --threads "$PLINK_THREADS" \
    --out "$PRUNE_PREFIX"

# 2. calculate genotype principal components (PLINK 1.9)
PC_COUNT=5
PCA_PREFIX="$PCA_DIR/genotype.pca"

plink \
    --bfile "$GENO_PREFIX" \
    --chr-set "$CHR_COUNT" \
    --allow-extra-chr \
    --extract "${PRUNE_PREFIX}.prune.in" \
    --pca "$PC_COUNT" \
    --threads "$PLINK_THREADS" \
    --out "$PCA_PREFIX"

# 3. estimate PEER factors (PEER R package v2.0)
EXPRESSION="$PHENO_DIR/expression.filtered.tsv"  # FPKM input for log2(FPKM + 1)
PEER_COUNT=20
PEER_FILE="$OUT_DIR/peer_factors.tsv"

Rscript "$PEER_SCRIPT" \
    --input "$EXPRESSION" \
    --factor-count "$PEER_COUNT" \
    --max-gene-missing 0.05 \
    --max-sample-missing 0.05 \
    --output "$PEER_FILE"

# 4. merge PEER factors and genotype principal components
COVARIATES="$OUT_DIR/covariates.tsv"

python3 "$PREPARE_COVARIATES_SCRIPT" \
    --peer "$PEER_FILE" \
    --eigenvec "${PCA_PREFIX}.eigenvec" \
    --samples "$SAMPLE_LIST" \
    --pc-count "$PC_COUNT" \
    --output "$COVARIATES"
