# RNA-seq genotyping workflow
# run steps 1-6 per sample, steps 7-8 per chromosome, then step 9

# prepare reference indexes
REF="/path/to/reference.fa"  # reference genome FASTA file
HISAT2_INDEX="/path/to/hisat2_index/reference"  # HISAT2 index prefix

mkdir -p "$(dirname "$HISAT2_INDEX")"

hisat2-build "$REF" "$HISAT2_INDEX"
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

# 2. alignment, sorting and indexing (HISAT2 v2.2.1; SAMtools v1.21)
BAM_DIR="/path/to/bam"
HISAT2_THREADS=12
SORT_THREADS=2
MIN_INTRON=20     # bp: Arabidopsis/soybean/maize/rice/wheat/Brassica napus/tomato/sorghum 20; cotton 40
MAX_INTRON=50000  # bp: Arabidopsis 12000; soybean 50000; maize/rice 20000; wheat 70000; cotton 100000; Brassica napus/tomato/sorghum 500000

mkdir -p "$BAM_DIR"

hisat2 \
    -x "$HISAT2_INDEX" \
    -p "$HISAT2_THREADS" \
    --min-intronlen "$MIN_INTRON" \
    --max-intronlen "$MAX_INTRON" \
    --dta \
    -1 "$CLEAN_DIR/${SAMPLE}_clean_1.fastq.gz" \
    -2 "$CLEAN_DIR/${SAMPLE}_clean_2.fastq.gz" \
    2> "$BAM_DIR/${SAMPLE}.hisat2.log" | \
    samtools sort \
        -@ "$SORT_THREADS" \
        -o "$BAM_DIR/${SAMPLE}.sorted.bam" \
        - \
        2> "$BAM_DIR/${SAMPLE}.sort.log"

samtools index -c -@ "$SORT_THREADS" "$BAM_DIR/${SAMPLE}.sorted.bam"

# 3. retain uniquely mapped reads (SAMtools v1.21)
UNIQUE_DIR="/path/to/unique_bam"
UNIQUE_BAM="$UNIQUE_DIR/${SAMPLE}.uniq.bam"

mkdir -p "$UNIQUE_DIR"

samtools view \
    -@ "$SORT_THREADS" \
    -b \
    -F 2308 \
    -d NH:1 \
    -o "$UNIQUE_BAM" \
    "$BAM_DIR/${SAMPLE}.sorted.bam" \
    > "$UNIQUE_DIR/${SAMPLE}.unique.log" 2>&1

samtools index -c -@ "$SORT_THREADS" "$UNIQUE_BAM"

# 4. add read groups and mark duplicates (GATK v4.6.2.0)
MARKDUP_DIR="/path/to/markdup"
JAVA_MEM="20g"
RG_BAM="$MARKDUP_DIR/${SAMPLE}.uniq.rg.bam"
MARKDUP_BAM="$MARKDUP_DIR/${SAMPLE}.uniq.markdup.bam"
DUP_METRICS="$MARKDUP_DIR/${SAMPLE}.dup_metrics.txt"
TMP_DIR="$MARKDUP_DIR/tmp/${SAMPLE}"

mkdir -p "$MARKDUP_DIR" "$TMP_DIR"

gatk \
    --java-options "-Xmx${JAVA_MEM}" \
    AddOrReplaceReadGroups \
    -I "$UNIQUE_BAM" \
    -O "$RG_BAM" \
    -RGID "$SAMPLE" \
    -RGLB "${SAMPLE}_lib" \
    -RGPL ILLUMINA \
    -RGPU "${SAMPLE}_unit" \
    -RGSM "$SAMPLE" \
    --SORT_ORDER coordinate \
    --TMP_DIR "$TMP_DIR" \
    --VALIDATION_STRINGENCY SILENT \
    > "$MARKDUP_DIR/${SAMPLE}.AddOrReplaceReadGroups.log" 2>&1

samtools index -c -@ "$SORT_THREADS" "$RG_BAM"

gatk \
    --java-options "-Xmx${JAVA_MEM}" \
    MarkDuplicates \
    -I "$RG_BAM" \
    -O "$MARKDUP_BAM" \
    -M "$DUP_METRICS" \
    --TMP_DIR "$TMP_DIR" \
    --VALIDATION_STRINGENCY SILENT \
    > "$MARKDUP_DIR/${SAMPLE}.markdup.log" 2>&1

samtools index -c -@ "$SORT_THREADS" "$MARKDUP_BAM"

# 5. split reads spanning introns (GATK SplitNCigarReads v4.6.2.0)
SPLIT_DIR="/path/to/split_bam"
SPLIT_BAM="$SPLIT_DIR/${SAMPLE}.uniq.markdup.split.bam"

mkdir -p "$SPLIT_DIR"

gatk \
    --java-options "-Xmx${JAVA_MEM}" \
    SplitNCigarReads \
    -R "$REF" \
    -I "$MARKDUP_BAM" \
    -O "$SPLIT_BAM" \
    --create-output-bam-index false \
    --tmp-dir "$TMP_DIR" \
    > "$SPLIT_DIR/${SAMPLE}.SplitNCigarReads.log" 2>&1

samtools index -c -@ "$SORT_THREADS" "$SPLIT_BAM"

# 6. per-sample variant calling (GATK HaplotypeCaller v4.6.2.0)
GVCF_DIR="/path/to/gvcf"
HC_THREADS=4

mkdir -p "$GVCF_DIR"

gatk \
    --java-options "-Xmx${JAVA_MEM}" \
    HaplotypeCaller \
    -R "$REF" \
    -I "$SPLIT_BAM" \
    -O "$GVCF_DIR/${SAMPLE}.g.vcf.gz" \
    -ERC GVCF \
    --dont-use-soft-clipped-bases true \
    --native-pair-hmm-threads "$HC_THREADS" \
    --create-output-variant-index false \
    --tmp-dir "$TMP_DIR" \
    > "$GVCF_DIR/${SAMPLE}.HaplotypeCaller.log" 2>&1

bcftools index "$GVCF_DIR/${SAMPLE}.g.vcf.gz"

# 7. import cohort gVCFs (GATK GenomicsDBImport v4.6.2.0)
SAMPLE_MAP="/path/to/sample_map.tsv"  # tab-delimited sample IDs and absolute gVCF paths
GENOMICSDB_DIR="/path/to/genomicsdb"
INTERVAL="chr1"  # chromosome name in the reference genome
BATCH_SIZE=20
DB_PATH="$GENOMICSDB_DIR/${INTERVAL}"
TMP_DIR="$GENOMICSDB_DIR/tmp/${INTERVAL}"
GVCF_HEADER="$GENOMICSDB_DIR/gvcf.header.vcf"

mkdir -p "$GENOMICSDB_DIR" "$TMP_DIR"

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

# 8. joint genotyping (GATK GenotypeGVCFs v4.6.2.0)
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

# 9. concatenate and index VCF files (BCFtools v1.21)
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
