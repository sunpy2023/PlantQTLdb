# GWAS result visualization

ASSOC_FILE="/path/to/gwas/trait_name.assoc.txt"
OUT_PREFIX="/path/to/gwas/plots/trait_name"
P_COLUMN="p_wald"  # GEMMA Wald-test P-value column
PLOT_SCRIPT="script/plot_gemma_result.py"

# 1. plot Manhattan and QQ plots
python3 "$PLOT_SCRIPT" \
    -i "$ASSOC_FILE" \
    -p "$P_COLUMN" \
    -o "$OUT_PREFIX"
