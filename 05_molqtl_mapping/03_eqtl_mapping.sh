# eQTL mapping workflow

GENOTYPE_VCF="/path/to/genotype_qc/genotype.qc.vcf.gz"
PHENO_DIR="/path/to/molqtl/phenotypes"
COVARIATES="/path/to/molqtl/covariates/covariates.tsv"
OUT_DIR="/path/to/molqtl/eqtl"
EQTL_SCRIPT="script/run_eqtl_mapping.R"

mkdir -p "$OUT_DIR"

# 1. cis/trans-eQTL mapping (MatrixEQTL v2.3; BCFtools v1.21)
EXPRESSION="$PHENO_DIR/expression.normalized.tsv"
GENE_LOCATIONS="$PHENO_DIR/expression.locations.tsv"
CIS_DISTANCE=200000  # bp flanking the gene interval
MIN_MAF=0.05
MAX_MISSING=0.10  # maximum missing genotype fraction per variant
COUNTED_ALLELE="ALT"  # count REF or ALT copies; beta is per additional copy

# cis/trans significance threshold: 0.05 / retained variants, calculated in run_eqtl_mapping.R
Rscript "$EQTL_SCRIPT" \
    --vcf "$GENOTYPE_VCF" \
    --expression "$EXPRESSION" \
    --gene-locations "$GENE_LOCATIONS" \
    --covariates "$COVARIATES" \
    --cis-distance "$CIS_DISTANCE" \
    --min-maf "$MIN_MAF" \
    --max-variant-missing "$MAX_MISSING" \
    --counted-allele "$COUNTED_ALLELE" \
    --out-dir "$OUT_DIR" \
    > "$OUT_DIR/eqtl.log" 2>&1
