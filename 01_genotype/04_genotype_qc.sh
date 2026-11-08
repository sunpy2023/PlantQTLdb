# genotype quality control workflow

# 1. standardize chromosome names (BCFtools v1.21)
INPUT_VCF="/path/to/genotype.vcf.gz"
CHR_MAP="/path/to/chromosome_map.tsv"  # reference ID<TAB>numeric chromosome ID (1..CHR_COUNT), no header
# use the same chromosome mapping for molecular feature coordinates
OUT_DIR="/path/to/genotype_qc"
BCFTOOLS_THREADS=12

mkdir -p "$OUT_DIR"

bcftools annotate \
    --rename-chrs "$CHR_MAP" \
    --threads "$BCFTOOLS_THREADS" \
    -Oz -o "$OUT_DIR/genotype.chrom.vcf.gz" \
    "$INPUT_VCF"

# 2. assign allele-specific variant identifiers (BCFtools v1.21)
bcftools annotate \
    --set-id '%CHROM:%POS:%REF:%ALT' \
    --threads "$BCFTOOLS_THREADS" \
    -Oz -o "$OUT_DIR/genotype.id.vcf.gz" \
    "$OUT_DIR/genotype.chrom.vcf.gz"

# 3. filter genotypes and export files (PLINK 2.0)
QC_PREFIX="$OUT_DIR/genotype.qc"
CHR_COUNT=20  # chromosome count: soybean 20; wheat 21; cotton 26
PLINK_THREADS=12

# remove matching duplicate records; stop if alleles, genotypes or annotations differ
plink2 \
    --vcf "$OUT_DIR/genotype.id.vcf.gz" \
    --chr-set "$CHR_COUNT" no-x no-y no-xy no-mt \
    --allow-extra-chr \
    --max-alleles 2 \
    --maf 0.05 \
    --geno 0.1 \
    --rm-dup error \
    --make-bed \
    --export vcf bgz id-paste=iid \
    --threads "$PLINK_THREADS" \
    --out "$QC_PREFIX" || exit 1

# 4. update variant statistics and index VCF files (BCFtools v1.21)
TAGGED_VCF="${QC_PREFIX}.tags.vcf.gz"

bcftools +fill-tags \
    "${QC_PREFIX}.vcf.gz" \
    --threads "$BCFTOOLS_THREADS" \
    -Oz -o "$TAGGED_VCF" \
    -- -t AC,AN,AF,MAF,NS,F_MISSING

mv "$TAGGED_VCF" "${QC_PREFIX}.vcf.gz"
bcftools index --threads "$BCFTOOLS_THREADS" "${QC_PREFIX}.vcf.gz"
