# RNA-seq preprocessing workflow

# 1. clean data (fastp v1.0.1)
SAMPLE="sample_id"
LAYOUT="PE"  # PE: paired-end; SE: single-end
CLEAN_DIR="/path/to/clean_fastq"
FASTP_THREADS=4

mkdir -p "$CLEAN_DIR"

if [ "$LAYOUT" = "PE" ]; then
    R1="/path/to/${SAMPLE}_1.fastq.gz"
    R2="/path/to/${SAMPLE}_2.fastq.gz"

    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$CLEAN_DIR/${SAMPLE}_clean_1.fastq.gz" \
        -O "$CLEAN_DIR/${SAMPLE}_clean_2.fastq.gz" \
        --json "$CLEAN_DIR/${SAMPLE}.fastp.json" \
        --html "$CLEAN_DIR/${SAMPLE}.fastp.html" \
        -w "$FASTP_THREADS" \
        > "$CLEAN_DIR/${SAMPLE}.fastp.log" 2>&1
else
    READ="/path/to/${SAMPLE}.fastq.gz"

    fastp \
        -i "$READ" \
        -o "$CLEAN_DIR/${SAMPLE}_clean.fastq.gz" \
        --json "$CLEAN_DIR/${SAMPLE}.fastp.json" \
        --html "$CLEAN_DIR/${SAMPLE}.fastp.html" \
        -w "$FASTP_THREADS" \
        > "$CLEAN_DIR/${SAMPLE}.fastp.log" 2>&1
fi

# 2. prepare reference index (HISAT2 v2.2.1)
REF="/path/to/reference.fa"  # reference genome FASTA file
GTF="/path/to/annotation.gtf"  # gene annotation GTF file
HISAT2_INDEX="/path/to/hisat2_index/reference"  # HISAT2 index prefix
EXON_FILE="${HISAT2_INDEX}.exon.txt"
SS_FILE="${HISAT2_INDEX}.ss.txt"
INDEX_THREADS=10

mkdir -p "$(dirname "$HISAT2_INDEX")"

hisat2_extract_exons.py "$GTF" > "$EXON_FILE"
hisat2_extract_splice_sites.py "$GTF" > "$SS_FILE"

hisat2-build \
    -p "$INDEX_THREADS" \
    --exon "$EXON_FILE" \
    --ss "$SS_FILE" \
    "$REF" \
    "$HISAT2_INDEX"

# 3. alignment, sorting and indexing (HISAT2 v2.2.1; SAMtools v1.21)
BAM_DIR="/path/to/bam"
HISAT2_THREADS=12
SORT_THREADS=2
MIN_INTRON=20     # bp: Arabidopsis/soybean/maize/rice/wheat/Brassica napus/tomato/sorghum 20; cotton 40
MAX_INTRON=50000  # bp: Arabidopsis 12000; soybean 50000; maize/rice 20000; wheat 70000; cotton 100000; Brassica napus/tomato/sorghum 500000

mkdir -p "$BAM_DIR"

if [ "$LAYOUT" = "PE" ]; then
    hisat2 \
        -x "$HISAT2_INDEX" \
        -p "$HISAT2_THREADS" \
        --min-intronlen "$MIN_INTRON" \
        --max-intronlen "$MAX_INTRON" \
        --dta \
        -1 "$CLEAN_DIR/${SAMPLE}_clean_1.fastq.gz" \
        -2 "$CLEAN_DIR/${SAMPLE}_clean_2.fastq.gz" \
        2> "$BAM_DIR/${SAMPLE}.hisat2.log"
else
    hisat2 \
        -x "$HISAT2_INDEX" \
        -p "$HISAT2_THREADS" \
        --min-intronlen "$MIN_INTRON" \
        --max-intronlen "$MAX_INTRON" \
        --dta \
        -U "$CLEAN_DIR/${SAMPLE}_clean.fastq.gz" \
        2> "$BAM_DIR/${SAMPLE}.hisat2.log"
fi | \
    samtools sort \
        -@ "$SORT_THREADS" \
        -o "$BAM_DIR/${SAMPLE}.sorted.bam" \
        - \
        2> "$BAM_DIR/${SAMPLE}.sort.log"

samtools index -c -@ "$SORT_THREADS" "$BAM_DIR/${SAMPLE}.sorted.bam"
