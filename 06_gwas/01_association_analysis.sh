# genome-wide association analysis workflow

GENO_PREFIX="/path/to/genotype.qc"  # PLINK BED/BIM/FAM prefix; chromosomes 1..CHR_COUNT
PCA_DIR="/path/to/genotype_pca"  # calculate once; share PCs between GWAS and molQTL
PHENO="/path/to/trait.txt"  # one column in FAM sample order; no header, missing = NA
TRAIT="trait_name"
OUT_DIR="/path/to/gwas"
CHR_COUNT=20  # nuclear chromosome count matching genotype QC
# chromosome counts: Arabidopsis 5, soybean 20, maize 10, rice 12, wheat 21, cotton 26,
# Brassica napus 19, tomato 12, sorghum 10
PLINK_THREADS=12

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

# prepare GEMMA covariates in FAM sample order
COVAR="$OUT_DIR/covariates.gemma.txt"
COVARIATE_SCRIPT="script/prepare_gwas_covariates.py"

python3 "$COVARIATE_SCRIPT" \
    "${GENO_PREFIX}.fam" \
    "${PCA_PREFIX}.eigenvec" \
    "$PHENO" \
    "$COVAR" \
    "$PC_COUNT"

# 3. calculate the relatedness matrix (GEMMA v0.98.5)
KINSHIP_PREFIX="${TRAIT}.kinship"

gemma \
    -bfile "$GENO_PREFIX" \
    -p "$PHENO" \
    -gk 2 \
    -outdir "$OUT_DIR" \
    -o "$KINSHIP_PREFIX"

# 4. genome-wide association analysis (GEMMA v0.98.5)
gemma \
    -bfile "$GENO_PREFIX" \
    -p "$PHENO" \
    -c "$COVAR" \
    -k "$OUT_DIR/${KINSHIP_PREFIX}.sXX.txt" \
    -lmm 1 \
    -outdir "$OUT_DIR" \
    -o "$TRAIT"

# 5. extract significant associations (Bonferroni correction)
ASSOC_FILE="$OUT_DIR/${TRAIT}.assoc.txt"

awk '
    BEGIN { FS = OFS = "\t" }
    NR == FNR {
        if (FNR > 1 && NF) n_tested++  # number of tested variants
        next
    }
    FNR == 1 {
        for (i = 1; i <= NF; i++) {
            if ($i == "p_wald") p_column = i
        }
        if (!p_column || !n_tested) {
            print "Empty GEMMA results or missing p_wald column" > "/dev/stderr"
            exit 1
        }
        threshold = 0.05 / n_tested
        print
        next
    }
    $p_column ~ /^[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?$/ {
        if ($p_column + 0 <= threshold) print
    }
    END { if (!n_tested) exit 1 }
' "$ASSOC_FILE" "$ASSOC_FILE" > "$OUT_DIR/${TRAIT}.significant.tsv"
