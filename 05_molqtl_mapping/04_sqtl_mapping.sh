# sQTL mapping workflow

GENOTYPE_VCF="/path/to/genotype_qc/genotype.qc.vcf.gz"
PHENO_DIR="/path/to/molqtl/phenotypes"
COVARIATES="/path/to/molqtl/covariates/covariates.tsv"
OUT_DIR="/path/to/molqtl/sqtl"
SQTL_SCRIPT="script/run_sqtl_mapping.R"

mkdir -p "$OUT_DIR"

# 1. cis/trans-sQTL mapping (MatrixEQTL v2.3; BCFtools v1.21)
PSI="$PHENO_DIR/splicing.normalized.tsv"  # rank-normalized PSI
JUNCTION_LOCATIONS="$PHENO_DIR/splicing.locations.tsv"
CIS_DISTANCE=200000  # bp flanking the junction interval
MIN_MAF=0.05
MAX_MISSING=0.10  # maximum missing genotype fraction per variant
COUNTED_ALLELE="REF"  # count REF or ALT copies; beta is per additional copy

# cis/trans significance threshold: 0.05 / retained variants, calculated in run_sqtl_mapping.R
Rscript "$SQTL_SCRIPT" \
    --vcf "$GENOTYPE_VCF" \
    --psi "$PSI" \
    --junction-locations "$JUNCTION_LOCATIONS" \
    --covariates "$COVARIATES" \
    --cis-distance "$CIS_DISTANCE" \
    --min-maf "$MIN_MAF" \
    --max-variant-missing "$MAX_MISSING" \
    --counted-allele "$COUNTED_ALLELE" \
    --out-dir "$OUT_DIR" \
    > "$OUT_DIR/sqtl.log" 2>&1
