# molecular phenotype preparation workflow

OUT_DIR="/path/to/molqtl/phenotypes"
SAMPLE_LIST="/path/to/qtl_samples.txt"  # matched sample IDs, one per line in output order, no header
CHR_MAP="/path/to/chromosome_map.tsv"  # chromosome mapping used in genotype QC
PREPARE_PHENOTYPES_SCRIPT="script/prepare_qtl_phenotypes.py"
NORMALIZE_SCRIPT="script/normalize_phenotypes.R"

mkdir -p "$OUT_DIR"

# 1. prepare gene expression phenotypes
FPKM="/path/to/expression_matrices/merged_gene_FPKM_matrix.csv"
GFF="/path/to/annotation.gff3"  # reference annotation with gene IDs matching the expression matrix

python3 "$PREPARE_PHENOTYPES_SCRIPT" \
    --type expression \
    --input "$FPKM" \
    --samples "$SAMPLE_LIST" \
    --annotation "$GFF" \
    --chrom-map "$CHR_MAP" \
    --max-missing-frac 0.2 \
    --min-mean-expression 0.1 \
    --expression-cutoff 0.1 \
    --min-expressed-frac 0.1 \
    --out-prefix "$OUT_DIR/expression"

# 2. prepare splicing phenotypes
PSI="/path/to/majiq/matrices/junction_psi.tsv"  # junction IDs and sample PSI values
JUNCTION_META="/path/to/majiq/modulizer/ALL/junctions.tsv"  # VOILA report with junction coordinates

python3 "$PREPARE_PHENOTYPES_SCRIPT" \
    --type psi \
    --input "$PSI" \
    --samples "$SAMPLE_LIST" \
    --junction-meta "$JUNCTION_META" \
    --chrom-map "$CHR_MAP" \
    --max-missing-frac 0.2 \
    --min-intermediate-frac 0.1 \
    --intermediate-low 0.1 \
    --intermediate-high 0.9 \
    --min-range 0.1 \
    --min-variance 0.001 \
    --out-prefix "$OUT_DIR/splicing"

# 3. prepare APA phenotypes
PDUI="/path/to/apa/matrices/apa_pdui.tsv"
APA_LOCATIONS="/path/to/apa/matrices/apa_locations.tsv"

python3 "$PREPARE_PHENOTYPES_SCRIPT" \
    --type pdui \
    --input "$PDUI" \
    --samples "$SAMPLE_LIST" \
    --locations "$APA_LOCATIONS" \
    --chrom-map "$CHR_MAP" \
    --max-missing-frac 0.2 \
    --min-intermediate-frac 0.1 \
    --intermediate-low 0.1 \
    --intermediate-high 0.9 \
    --min-range 0.1 \
    --min-variance 0.001 \
    --out-prefix "$OUT_DIR/apa"

# 4. apply rank-based inverse normal transformation
for PHENOTYPE in expression splicing apa; do
    Rscript "$NORMALIZE_SCRIPT" \
        --input "$OUT_DIR/${PHENOTYPE}.filtered.tsv" \
        --output "$OUT_DIR/${PHENOTYPE}.normalized.tsv"
done
