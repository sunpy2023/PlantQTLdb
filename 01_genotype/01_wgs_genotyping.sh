# whole-genome resequencing genotyping workflow
# run steps 1-4 per sample, steps 5-6 per chromosome, then step 7

# prepare reference indexes
REF="/path/to/reference.fa"  # reference genome FASTA file

bwa-mem2 index "$REF"
samtools faidx "$REF"
gatk CreateSequenceDictionary -R "$REF" -O "${REF%.*}.dict"

# 1. clean data (fastp v1.0.1)
SAMPLE="sample_id"
R1="/path/to/${SAMPLE}_1.fastq.gz"
R2="/path/to/${SAMPLE}_2.fastq.gz"
CLEAN_DIR="/path/to/clean_fastq"
FASTP_THREADS=4

mkdir -p "$CLEAN_DIR"

fastp \
    -i "$R1" \
    -I "$R2" \
    -o "$CLEAN_DIR/${SAMPLE}_clean_1.fastq.gz" \
    -O "$CLEAN_DIR/${SAMPLE}_clean_2.fastq.gz" \
    --json "$CLEAN_DIR/${SAMPLE}.fastp.json" \
    --html "$CLEAN_DIR/${SAMPLE}.fastp.html" \
    -w "$FASTP_THREADS" \
    > "$CLEAN_DIR/${SAMPLE}.fastp.log" 2>&1

# 2. alignment, sorting and indexing (BWA-MEM2 v2.2.1; SAMtools v1.21)
BAM_DIR="/path/to/bam"
BWA_THREADS=12
SORT_THREADS=2
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}\tPU:${SAMPLE}"

mkdir -p "$BAM_DIR"

bwa-mem2 mem \
    -M \
    -R "$RG" \
    -t "$BWA_THREADS" \
    "$REF" \
    "$CLEAN_DIR/${SAMPLE}_clean_1.fastq.gz" \
    "$CLEAN_DIR/${SAMPLE}_clean_2.fastq.gz" \
    2> "$BAM_DIR/${SAMPLE}.bwa_mem2.log" | \
    samtools sort \
        -@ "$SORT_THREADS" \
        -o "$BAM_DIR/${SAMPLE}.sorted.bam" \
        - \
        2> "$BAM_DIR/${SAMPLE}.sort.log"

samtools index -c -@ "$SORT_THREADS" "$BAM_DIR/${SAMPLE}.sorted.bam"

# 3. mark duplicates (Picard v3.4.0)
PICARD="/path/to/picard.jar"
MARKDUP_DIR="/path/to/markdup"
JAVA_MEM="20g"
MARKDUP_BAM="$MARKDUP_DIR/${SAMPLE}.markdup.bam"
DUP_METRICS="$MARKDUP_DIR/${SAMPLE}.dup_metrics.txt"
TMP_DIR="$MARKDUP_DIR/tmp/${SAMPLE}"

mkdir -p "$MARKDUP_DIR" "$TMP_DIR"

java \
    -Xmx"$JAVA_MEM" \
    -Djava.io.tmpdir="$TMP_DIR" \
    -jar "$PICARD" \
    MarkDuplicates \
    I="$BAM_DIR/${SAMPLE}.sorted.bam" \
    O="$MARKDUP_BAM" \
    M="$DUP_METRICS" \
    TMP_DIR="$TMP_DIR" \
    CREATE_INDEX=false \
    VALIDATION_STRINGENCY=SILENT \
    > "$MARKDUP_DIR/${SAMPLE}.markdup.log" 2>&1

samtools index -c -@ "$SORT_THREADS" "$MARKDUP_BAM"

# 4. per-sample variant calling (GATK HaplotypeCaller v4.6.2.0)
GVCF_DIR="/path/to/gvcf"
HC_THREADS=4

mkdir -p "$GVCF_DIR"

gatk \
    --java-options "-Xmx${JAVA_MEM}" \
    HaplotypeCaller \
    -R "$REF" \
    -I "$MARKDUP_BAM" \
    -O "$GVCF_DIR/${SAMPLE}.g.vcf.gz" \
    -ERC GVCF \
    --native-pair-hmm-threads "$HC_THREADS" \
    --create-output-variant-index false \
    --tmp-dir "$TMP_DIR" \
    > "$GVCF_DIR/${SAMPLE}.HaplotypeCaller.log" 2>&1

# create a CSI index (BCFtools v1.21)
bcftools index "$GVCF_DIR/${SAMPLE}.g.vcf.gz"

# 5. import cohort gVCFs (GATK GenomicsDBImport v4.6.2.0)
SAMPLE_MAP="/path/to/sample_map.tsv"  # tab-delimited sample IDs and absolute gVCF paths
GENOMICSDB_DIR="/path/to/genomicsdb"
INTERVAL="chr1"  # chromosome name in the reference genome
BATCH_SIZE=20
DB_PATH="$GENOMICSDB_DIR/${INTERVAL}"
TMP_DIR="$GENOMICSDB_DIR/tmp/${INTERVAL}"
GVCF_HEADER="$GENOMICSDB_DIR/gvcf.header.vcf"

mkdir -p "$GENOMICSDB_DIR" "$TMP_DIR"

# use gVCFs with matching reference sequences and header definitions
FIRST_GVCF=$(awk -F '\t' 'NR == 1 {print $2; exit}' "$SAMPLE_MAP")
bcftools view -h "$FIRST_GVCF" | cut -f1-8 > "$GVCF_HEADER"

gatk \
    --java-options "-Xmx${JAVA_MEM}" \
    GenomicsDBImport \
    --sample-name-map "$SAMPLE_MAP" \
    --genomicsdb-workspace-path "$DB_PATH" \
    -L "$INTERVAL" \
    --batch-size "$BATCH_SIZE" \
    --bypass-feature-reader true \
    --avoid-nio true \
    --header "$GVCF_HEADER" \
    --tmp-dir "$TMP_DIR" \
    > "$GENOMICSDB_DIR/${INTERVAL}.GenomicsDBImport.log" 2>&1

# 6. joint genotyping (GATK GenotypeGVCFs v4.6.2.0)
VCF_DIR="/path/to/vcf"

mkdir -p "$VCF_DIR"

gatk \
    --java-options "-Xmx${JAVA_MEM}" \
    GenotypeGVCFs \
    -R "$REF" \
    -V "gendb://${DB_PATH}" \
    -L "$INTERVAL" \
    -O "$VCF_DIR/${INTERVAL}.raw.vcf.gz" \
    --create-output-variant-index false \
    --tmp-dir "$TMP_DIR" \
    > "$VCF_DIR/${INTERVAL}.GenotypeGVCFs.log" 2>&1

bcftools index "$VCF_DIR/${INTERVAL}.raw.vcf.gz"

# 7. concatenate and index VCF files
VCF_LIST="/path/to/vcf_list.txt"  # one VCF path per line, in reference chromosome order
FINAL_VCF="$VCF_DIR/cohort.raw.vcf.gz"
BCFTOOLS_THREADS=12

bcftools concat \
    -f "$VCF_LIST" \
    -Oz \
    -o "$FINAL_VCF" \
    --threads "$BCFTOOLS_THREADS" \
    > "$VCF_DIR/cohort.concat.log" 2>&1

bcftools index --threads "$BCFTOOLS_THREADS" "$FINAL_VCF"
