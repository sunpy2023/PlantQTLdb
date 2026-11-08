# APAQTL mapping workflow

GENOTYPE_VCF="/path/to/genotype_qc/genotype.qc.vcf.gz"
PHENO_DIR="/path/to/molqtl/phenotypes"
COVARIATES="/path/to/molqtl/covariates/covariates.tsv"
OUT_DIR="/path/to/molqtl/apaqtl"
APAQTL_SCRIPT="script/run_apaqtl_mapping.R"

mkdir -p "$OUT_DIR"

# 1. cis/trans-APAQTL mapping (MatrixEQTL v2.3; BCFtools v1.21)
PDUI="$PHENO_DIR/apa.normalized.tsv"  # rank-normalized PDUI
APA_LOCATIONS="$PHENO_DIR/apa.locations.tsv"
CIS_DISTANCE=200000  # bp flanking the APA interval
MIN_MAF=0.05
MAX_MISSING=0.10  # maximum missing genotype fraction per variant
COUNTED_ALLELE="REF"  # count REF or ALT copies; beta is per additional copy

# cis/trans significance threshold: 0.05 / retained variants, calculated in run_apaqtl_mapping.R
Rscript "$APAQTL_SCRIPT" \
    --vcf "$GENOTYPE_VCF" \
    --pdui "$PDUI" \
    --apa-locations "$APA_LOCATIONS" \
    --covariates "$COVARIATES" \
    --cis-distance "$CIS_DISTANCE" \
    --min-maf "$MIN_MAF" \
    --max-variant-missing "$MAX_MISSING" \
    --counted-allele "$COUNTED_ALLELE" \
    --out-dir "$OUT_DIR" \
    > "$OUT_DIR/apaqtl.log" 2>&1
