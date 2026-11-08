# genotype merging workflow

CONSENSUS_SCRIPT="script/vcf_assay_consensus.py"

# S006: RNA-seq and whole-genome resequencing
# run step 1 once, repeat steps 2-4 for each chromosome, then run step 5

# 1. prepare sample identifiers (BCFtools v1.21)
REF="/path/to/soybean_reference.fa"  # reference genome FASTA file
RNA_VCF="/path/to/S006/rnaseq.raw.vcf.gz"
WGS_VCF="/path/to/S006/wgs.raw.vcf.gz"
RNA_MAP="/path/to/S006/rna.rename.tsv"  # old ID<TAB>material_id__RNA, no header
WGS_MAP="/path/to/S006/wgs.rename.tsv"  # old ID<TAB>material_id__WGS, no header
MATERIALS="/path/to/S006/materials.txt"  # material IDs to retain, one per line in output order
WORK_DIR="/path/to/S006/genotype_merging"
BCFTOOLS_THREADS=12

mkdir -p "$WORK_DIR"

bcftools reheader -f "$REF.fai" -s "$RNA_MAP" -o "$WORK_DIR/RNA.renamed.vcf.gz" "$RNA_VCF"
bcftools reheader -f "$REF.fai" -s "$WGS_MAP" -o "$WORK_DIR/WGS.renamed.vcf.gz" "$WGS_VCF"
bcftools index --threads "$BCFTOOLS_THREADS" "$WORK_DIR/RNA.renamed.vcf.gz"
bcftools index --threads "$BCFTOOLS_THREADS" "$WORK_DIR/WGS.renamed.vcf.gz"

# 2. filter RNA-seq and WGS SNPs independently (BCFtools v1.21)
CHROM="Chr01"  # reference chromosome name, excluding scaffolds
QC_TAGS='AC,AN,MAF,F_MISSING,N_HET:1=int(N_PASS(GT="het")),N_HOM_REF:1=int(N_PASS(GT="RR")),N_HOM_ALT:1=int(N_PASS(GT="AA"))'
S006_FILTER='INFO/F_MISSING<=0.5 && INFO/MAF>=0.05 && INFO/N_HET*5<N_SAMPLES &&
    ((INFO/AC[0]*2<=INFO/AN && INFO/N_HOM_ALT>=10) ||
    (INFO/AC[0]*2>INFO/AN && INFO/N_HOM_REF>=10))'

for SOURCE in RNA WGS; do
    bcftools view \
        -r "$CHROM" -v snps -m2 -M2 \
        "$WORK_DIR/${SOURCE}.renamed.vcf.gz" -Ou | \
        bcftools norm -f "$REF" -c e -d exact -Ou | \
        bcftools +setGT -Ou -- -t ./x -n . | \
        bcftools +fill-tags -Ou -- -t "$QC_TAGS" | \
        bcftools view -i "$S006_FILTER" -Ou | \
        bcftools annotate \
            -x INFO,^FORMAT/GT \
            --threads "$BCFTOOLS_THREADS" \
            -Oz -o "$WORK_DIR/${SOURCE}.${CHROM}.qc.vcf.gz"
done

# 3. impute RNA-seq and WGS genotypes independently (Beagle v5.0)
BEAGLE="/path/to/beagle.12Jul19.0df.jar"
JAVA_MEM="20g"
BEAGLE_THREADS=12
SEED=12345

for SOURCE in RNA WGS; do
    java \
        -Xmx"$JAVA_MEM" \
        -jar "$BEAGLE" \
        gt="$WORK_DIR/${SOURCE}.${CHROM}.qc.vcf.gz" \
        out="$WORK_DIR/${SOURCE}.${CHROM}.beagle" \
        nthreads="$BEAGLE_THREADS" \
        seed="$SEED" \
        > "$WORK_DIR/${SOURCE}.${CHROM}.beagle.stdout.log" 2>&1

    bcftools reheader \
        -f "$REF.fai" \
        -o "$WORK_DIR/${SOURCE}.${CHROM}.imputed.vcf.gz" \
        "$WORK_DIR/${SOURCE}.${CHROM}.beagle.vcf.gz"
done

# 4. merge genotypes with RNA-seq priority (BCFtools v1.21; Python 3)
# retain RNA alleles at positions with different ALT alleles between sources
bcftools +fill-tags \
    "$WORK_DIR/RNA.${CHROM}.imputed.vcf.gz" \
    -Oz -o "$WORK_DIR/RNA.${CHROM}.tagged.vcf.gz" \
    -- -t 'RNA_SOURCE:1=int(1)'

bcftools index --threads "$BCFTOOLS_THREADS" "$WORK_DIR/RNA.${CHROM}.tagged.vcf.gz"
bcftools index --threads "$BCFTOOLS_THREADS" "$WORK_DIR/WGS.${CHROM}.imputed.vcf.gz"

bcftools merge \
    -m none -i RNA_SOURCE:max \
    "$WORK_DIR/RNA.${CHROM}.tagged.vcf.gz" \
    "$WORK_DIR/WGS.${CHROM}.imputed.vcf.gz" -Ov | \
    python3 "$CONSENSUS_SCRIPT" \
        --materials "$MATERIALS" \
        --other-source WGS \
        --conflict rna | \
    bcftools +fill-tags \
        --threads "$BCFTOOLS_THREADS" \
        -Oz -o "$WORK_DIR/${CHROM}.merged.vcf.gz" \
        -- -t AC,AN,AF,MAF,NS,F_MISSING

bcftools index --threads "$BCFTOOLS_THREADS" "$WORK_DIR/${CHROM}.merged.vcf.gz"

# 5. concatenate and index merged VCF files (BCFtools v1.21)
VCF_LIST="/path/to/S006/merged_vcf.list"  # one chromosome VCF path per line in reference order
FINAL_VCF="$WORK_DIR/S006.merged.vcf.gz"

bcftools concat \
    --threads "$BCFTOOLS_THREADS" \
    -f "$VCF_LIST" \
    -Oz -o "$FINAL_VCF"

bcftools index --threads "$BCFTOOLS_THREADS" "$FINAL_VCF"

# S028: RNA-seq and regulatory-region capture
# run step 6 once, repeat steps 7-9 for each chromosome, then run step 10

# 6. prepare sample identifiers (BCFtools v1.21)
REF="/path/to/wheat_reference.fa"  # reference genome FASTA file
RNA_VCF="/path/to/S028/rnaseq.raw.vcf.gz"
CAPTURE_VCF="/path/to/S028/capture.raw.vcf.gz"
RNA_MAP="/path/to/S028/rna.rename.tsv"  # old ID<TAB>material_id__RNA, no header
CAPTURE_MAP="/path/to/S028/capture.rename.tsv"  # old ID<TAB>material_id__CAPTURE, no header
MATERIALS="/path/to/S028/materials.txt"  # material IDs for merging and imputation, one per line in output order
WORK_DIR="/path/to/S028/genotype_merging"
BCFTOOLS_THREADS=12

mkdir -p "$WORK_DIR"

bcftools reheader -f "$REF.fai" -s "$RNA_MAP" -o "$WORK_DIR/RNA.renamed.vcf.gz" "$RNA_VCF"
bcftools reheader -f "$REF.fai" -s "$CAPTURE_MAP" -o "$WORK_DIR/CAPTURE.renamed.vcf.gz" "$CAPTURE_VCF"
bcftools index --threads "$BCFTOOLS_THREADS" "$WORK_DIR/RNA.renamed.vcf.gz"
bcftools index --threads "$BCFTOOLS_THREADS" "$WORK_DIR/CAPTURE.renamed.vcf.gz"

# 7. filter RNA-seq and capture SNPs independently (BCFtools v1.21)
CHROM="chr1A"  # reference chromosome name, excluding unplaced scaffolds
S028_FILTER='INFO/NS>0 && INFO/F_MISSING<=0.5 && INFO/AC_Het<=0.03*INFO/NS'

for SOURCE in RNA CAPTURE; do
    SOURCE_FILTER="$S028_FILTER"
    if [ "$SOURCE" = "CAPTURE" ]; then
        SOURCE_FILTER="$SOURCE_FILTER && INFO/MAF>=0.05"
    fi

    bcftools view \
        -r "$CHROM" -v snps -m2 -M2 \
        "$WORK_DIR/${SOURCE}.renamed.vcf.gz" -Ou | \
        bcftools norm -f "$REF" -c e -d exact -Ou | \
        bcftools +setGT -Ou -- -t ./x -n . | \
        bcftools +setGT -Ou -- -t q -n . -i 'FMT/DP<3 | FMT/DP="."' | \
        bcftools +fill-tags -Ou -- -t AC,AN,AF,MAF,F_MISSING,AC_Het,NS | \
        bcftools view \
            -i "$SOURCE_FILTER" \
            --threads "$BCFTOOLS_THREADS" \
            -Oz -o "$WORK_DIR/${SOURCE}.${CHROM}.qc.vcf.gz"
done

bcftools index --threads "$BCFTOOLS_THREADS" "$WORK_DIR/RNA.${CHROM}.qc.vcf.gz"
bcftools index --threads "$BCFTOOLS_THREADS" "$WORK_DIR/CAPTURE.${CHROM}.qc.vcf.gz"

# 8. merge genotypes and set discordant calls to missing (BCFtools v1.21; Python 3)
# exclude positions with different ALT alleles between sources
bcftools merge \
    -m snps \
    "$WORK_DIR/RNA.${CHROM}.qc.vcf.gz" \
    "$WORK_DIR/CAPTURE.${CHROM}.qc.vcf.gz" -Ou | \
    bcftools view -m2 -M2 -Ov | \
    python3 "$CONSENSUS_SCRIPT" \
        --materials "$MATERIALS" \
        --other-source CAPTURE \
        --conflict missing | \
    bcftools +fill-tags -Ou -- -t AC,AN,AF,MAF,NS,F_MISSING | \
    bcftools view \
        -i 'INFO/MAF>0' \
        --threads "$BCFTOOLS_THREADS" \
        -Oz -o "$WORK_DIR/${CHROM}.consensus.vcf.gz"

# 9. impute merged genotypes (Beagle v5.0)
BEAGLE="/path/to/beagle.12Jul19.0df.jar"
JAVA_MEM="80g"
BEAGLE_THREADS=8
BEAGLE_NE=12000  # effective population size used for S028

java \
    -Xmx"$JAVA_MEM" \
    -jar "$BEAGLE" \
    gt="$WORK_DIR/${CHROM}.consensus.vcf.gz" \
    out="$WORK_DIR/${CHROM}.beagle" \
    ne="$BEAGLE_NE" \
    nthreads="$BEAGLE_THREADS" \
    > "$WORK_DIR/${CHROM}.beagle.stdout.log" 2>&1

bcftools reheader \
    -f "$REF.fai" \
    -o "$WORK_DIR/${CHROM}.merged.vcf.gz" \
    "$WORK_DIR/${CHROM}.beagle.vcf.gz"

bcftools index --threads "$BCFTOOLS_THREADS" "$WORK_DIR/${CHROM}.merged.vcf.gz"

# 10. concatenate and index merged VCF files (BCFtools v1.21)
VCF_LIST="/path/to/S028/merged_vcf.list"  # one chromosome VCF path per line in reference order
FINAL_VCF="$WORK_DIR/S028.merged.vcf.gz"

bcftools concat \
    --threads "$BCFTOOLS_THREADS" \
    -f "$VCF_LIST" \
    -Oz -o "$FINAL_VCF"

bcftools index --threads "$BCFTOOLS_THREADS" "$FINAL_VCF"
