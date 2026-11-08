# RNA-seq expression quantification workflow
# run step 1 for each sample, then steps 2-3 for all samples in the same cohort

# 1. expression quantification (StringTie v3.0.0)
SAMPLE="sample_id"
BAM_DIR="/path/to/bam"
GTF="/path/to/annotation.gtf"  # gene annotation matching the reference genome
QUANT_DIR="/path/to/stringtie"
THREADS=10
STRAND=""  # empty: unstranded; --rf: fr-firststrand; --fr: fr-secondstrand

mkdir -p "$QUANT_DIR/$SAMPLE"

stringtie \
    "$BAM_DIR/${SAMPLE}.sorted.bam" \
    -e \
    -B \
    -G "$GTF" \
    -p "$THREADS" \
    $STRAND \
    -A "$QUANT_DIR/$SAMPLE/${SAMPLE}.gene_abund.tsv" \
    -o "$QUANT_DIR/$SAMPLE/${SAMPLE}.gtf" \
    > "$QUANT_DIR/$SAMPLE/${SAMPLE}.stringtie.log" 2>&1

# 2. merge gene expression matrices (Python 3)
SAMPLE_LIST="/path/to/samples.txt"  # sample IDs, one per line in output order, no header
MATRIX_DIR="/path/to/expression_matrices"
MERGE_EXPRESSION_SCRIPT="script/merge_gene_expression.py"
GTF_LIST="$MATRIX_DIR/sample_gtf.tsv"  # sample ID<TAB>StringTie GTF path, no header

mkdir -p "$MATRIX_DIR"

while IFS= read -r SAMPLE_ID || [ -n "$SAMPLE_ID" ]; do
    printf '%s\t%s\n' "$SAMPLE_ID" "$QUANT_DIR/$SAMPLE_ID/${SAMPLE_ID}.gtf"
done < "$SAMPLE_LIST" > "$GTF_LIST"

python3 "$MERGE_EXPRESSION_SCRIPT" \
    "$GTF_LIST" \
    "$QUANT_DIR" \
    "$MATRIX_DIR"

# 3. generate estimated count matrices (Python 3)
MERGE_COUNTS_SCRIPT="script/merge_gtf_counts.py"
READ_LEN=75  # read length (bp) used to estimate counts

READ_LEN="$READ_LEN" python3 "$MERGE_COUNTS_SCRIPT" \
    "$GTF_LIST" \
    "$MATRIX_DIR/merged_gene_count_matrix.csv" \
    "$MATRIX_DIR/merged_mRNA_count_matrix.csv"
