"""plot GEMMA association results as Manhattan and QQ plots.

usage: python3 plot_gemma_result.py -i trait.assoc.txt -p p_wald -o plots/trait
requires NumPy, pandas, Matplotlib and SciPy.
input columns: chr (positive integer), ps (1-based position), and the chosen P column.
the Bonferroni threshold is 0.05 divided by all input variant records, before
omitting missing P-values. P=0 is retained and capped only for display.
output: <prefix>.manhattan.png/.pdf and <prefix>.qq.png/.pdf, at 300 dpi.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2


MIN_P = np.nextafter(0.0, 1.0)


def load_results(path, p_column):
    data = pd.read_csv(path, sep="\t", usecols=["chr", "ps", p_column],
                       float_precision="round_trip")
    n_tested = len(data)  # number of tested variants before omitting missing P-values
    if not n_tested:
        raise ValueError("No tested variants in the association file")
    p_values = pd.to_numeric(data[p_column], errors="raise").to_numpy(dtype=float)
    missing = np.isnan(p_values)
    if np.any(~missing & (~np.isfinite(p_values) | (p_values < 0) | (p_values > 1))):
        raise ValueError("P-values must be between 0 and 1, or missing")
    data = data.loc[~missing].copy()
    if data.empty:
        raise ValueError("No nonmissing P-values to plot")
    data[p_column] = p_values[~missing]
    for column in ("chr", "ps"):
        values = pd.to_numeric(data[column], errors="raise").to_numpy(dtype=float)
        if np.any(~np.isfinite(values) | (values < 1) | (values != np.floor(values))):
            raise ValueError(f"{column} must contain positive integers for plotted variants")
        data[column] = values.astype(np.int64)
    data["logp"] = -np.log10(data[p_column].clip(lower=MIN_P))
    return data, n_tested


def plot_manhattan(data, n_tested):
    fig, ax = plt.subplots(figsize=(14, 5))
    colors = ("#4e79a7", "#f28e2b")
    offset = 0
    ticks, labels = [], []
    for index, (chromosome, group) in enumerate(data.groupby("chr", sort=True)):
        positions = group["ps"].to_numpy() + offset
        ax.scatter(positions, group["logp"], color=colors[index % 2],
                   s=5, linewidths=0, rasterized=True)
        ticks.append((positions.min() + positions.max()) / 2)
        labels.append(str(chromosome))
        offset += int(group["ps"].max())
    threshold = 0.05 / n_tested
    ax.axhline(-np.log10(threshold), color="#c44e52", linestyle="--", linewidth=1,
               label=f"P = 0.05/{n_tested:,} = {threshold:.2e}")
    ax.set_xticks(ticks, labels)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$-\log_{10}(P)$")
    ax.set_title("Manhattan plot")
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1, 1.01), fontsize=9)
    return fig


def plot_qq(data, p_column):
    p_values = np.sort(data[p_column].to_numpy())
    observed = -np.log10(np.maximum(p_values, MIN_P))
    expected = -np.log10(np.arange(1, len(p_values) + 1) / (len(p_values) + 1))
    lambda_gc = np.median(chi2.isf(p_values, df=1)) / chi2.isf(0.5, df=1)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(expected, observed, s=7, color="#4e79a7", linewidths=0, rasterized=True)
    x_limit = expected.max() * 1.05
    y_limit = max(expected.max(), observed.max()) * 1.05
    ax.plot([0, x_limit], [0, x_limit], color="#c44e52", linestyle="--", linewidth=1)
    ax.set_xlim(0, x_limit)
    ax.set_ylim(0, y_limit)
    ax.set_xlabel(r"Expected $-\log_{10}(P)$")
    ax.set_ylabel(r"Observed $-\log_{10}(P)$")
    ax.set_title(rf"QQ plot ($\lambda_{{GC}}$ = {lambda_gc:.3f})")
    return fig


def save_figure(fig, prefix, kind, zero_count):
    for ax in fig.axes:
        ax.spines[["top", "right"]].set_visible(False)
    if zero_count:
        fig.text(0.01, 0.01, f"P=0 shown at -log10(P)={-np.log10(MIN_P):.1f}", fontsize=8)
    fig.tight_layout(rect=(0, 0.035 if zero_count else 0, 1, 1))
    for extension in ("png", "pdf"):
        path = f"{prefix}.{kind}.{extension}"
        fig.savefig(path, dpi=300)
        print(f"Saved {path}", flush=True)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, help="Complete GEMMA association table")
    parser.add_argument("-p", "--p_col", default="p_wald", help="P-value column (default: p_wald)")
    parser.add_argument("-o", "--out", required=True, help="Output filename prefix")
    args = parser.parse_args()
    data, n_tested = load_results(args.input, args.p_col)
    zero_count = int((data[args.p_col] == 0).sum())
    print(f"Tested variants: {n_tested}; plotted: {len(data)}; P=0: {zero_count}; "
          f"Bonferroni threshold: {0.05 / n_tested:.17g}", flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "pdf.fonttype": 42})
    save_figure(plot_manhattan(data, n_tested), args.out, "manhattan", zero_count)
    save_figure(plot_qq(data, args.p_col), args.out, "qq", zero_count)


if __name__ == "__main__":
    main()
