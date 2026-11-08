# map cis/trans-APAQTLs from a QC VCF, normalized PDUI and complete covariates
# run from the pipeline root; reuse the observed-sample regression and genotype QC helpers
# raw PDUI QC and rank normalization are performed in phenotype preparation
# PDUI NA are not imputed; each APA trait uses its observed samples only
# counted allele defaults to REF for the original TXT matrices; use ALT for original VCF GT dosage
# outputs: cis.significant.tsv, trans.significant.tsv, cis_apa_traits.normalized.tsv,
# cis_apa_traits.all_cis.tsv; association tables include n_samples and df

source("script/run_sqtl_mapping.R")

if (sys.nframe() == 0L) main("apa")
