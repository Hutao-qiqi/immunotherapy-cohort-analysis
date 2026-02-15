#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
from dataclasses import dataclass
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from lifelines import CoxPHFitter


def _default_base_dir() -> Path:
    return Path(__file__).resolve().parent


def _default_survival_file(base_dir: Path) -> Path:
    return (base_dir / ".." / "survival_curves" / "updated_survival_data.txt").resolve()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Forest-plot-like summary for selected genes")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=_default_base_dir(),
        help="Directory containing expression + annotation files",
    )
    parser.add_argument(
        "--expr-file",
        type=Path,
        default=None,
        help="Expression matrix file (default: <base-dir>/combined_expression_combat_corrected.txt)",
    )
    parser.add_argument(
        "--survival-file",
        type=Path,
        default=None,
        help="Survival file (default: ../survival_curves/updated_survival_data.txt)",
    )
    parser.add_argument(
        "--annotation-file",
        type=Path,
        default=None,
        help="Annotation file (default: <base-dir>/MYC_PVT1_annotation.txt)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=_default_base_dir() / "outputs",
        help="Output directory (default: ./outputs)",
    )
    return parser.parse_args()

# Interaction-significant (q<0.05) and immune / immunotherapy-adjacent candidates
GENES = ["CDHR2", "PROZ", "SERPINA1", "VTN", "CLEC1B", "CEACAM5", "HRG", "FCN3"]


@dataclass
class CoxInputs:
    base: pd.DataFrame
    expr: pd.DataFrame


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def load_inputs(expr_file: Path, survival_file: Path, annotation_file: Path) -> CoxInputs:
    expr = pd.read_csv(expr_file, sep="\t", index_col=0)

    surv = pd.read_csv(survival_file, sep=r"\s+", engine="python")
    if "Sample_ID" not in surv.columns:
        raise ValueError("Survival file must contain column: Sample_ID")
    surv = surv.set_index("Sample_ID")

    annot = pd.read_csv(annotation_file, sep="\t")
    if "Sample" not in annot.columns or "MYC_PVT1_Status" not in annot.columns:
        raise ValueError("MYC_PVT1_annotation.txt must contain columns: Sample, MYC_PVT1_Status")
    annot = annot.set_index("Sample")
    annot = annot[annot["MYC_PVT1_Status"].isin(["hi_hi", "lo_lo"])]

    # Harmonize sample IDs
    expr_cols = expr.columns.astype(str)
    expr.columns = expr_cols
    surv.index = surv.index.astype(str)
    annot.index = annot.index.astype(str)

    common = sorted(set(expr.columns) & set(surv.index) & set(annot.index))
    if len(common) < 50:
        raise ValueError(f"Too few common samples: {len(common)}")

    surv = surv.loc[common].copy()
    annot = annot.loc[common].copy()

    for col in ["OS_months", "OS_event"]:
        if col not in surv.columns:
            raise ValueError(f"Survival file must contain column: {col}")
        surv[col] = pd.to_numeric(surv[col], errors="coerce")

    base = pd.DataFrame(index=common)
    base["OS_months"] = surv["OS_months"]
    base["OS_event"] = surv["OS_event"]

    if "Dataset" in surv.columns:
        base["Dataset"] = surv["Dataset"].astype(str)

    base = base.dropna(subset=["OS_months", "OS_event"])
    base = base[base["OS_months"] > 0]

    base["status"] = (annot.loc[base.index, "MYC_PVT1_Status"] == "hi_hi").astype(int)

    return CoxInputs(base=base, expr=expr)


def split_high_within_status(df: pd.DataFrame, expr_col: str) -> pd.Series:
    gene_high = pd.Series(0, index=df.index, dtype=int)
    for status_value in [0, 1]:
        idx = df.index[df["status"] == status_value]
        if len(idx) < 10:
            continue
        ordered = df.loc[idx, expr_col].sort_values(kind="mergesort")
        half = len(ordered) // 2
        high_idx = ordered.index[half:]
        gene_high.loc[high_idx] = 1
    return gene_high


def fit_interaction_for_gene(inputs: CoxInputs, gene: str) -> dict:
    if gene not in inputs.expr.index:
        raise KeyError(f"Gene not found in expression matrix: {gene}")

    df = inputs.base.copy()
    df["expr"] = inputs.expr.loc[gene, df.index].astype(float)
    df = df.dropna(subset=["expr"]).copy()

    df["gene_high"] = split_high_within_status(df, "expr")
    df["gene_high_status"] = df["gene_high"].astype(float) * df["status"].astype(float)

    cols = ["OS_months", "OS_event", "gene_high", "status", "gene_high_status"]
    strata = None
    if "Dataset" in df.columns:
        cols.append("Dataset")
        strata = ["Dataset"]

    cph = CoxPHFitter()
    cph.fit(df[cols], duration_col="OS_months", event_col="OS_event", strata=strata)

    summ = cph.summary
    beta_g = float(summ.loc["gene_high", "coef"])
    beta_i = float(summ.loc["gene_high_status", "coef"])
    p_i = float(summ.loc["gene_high_status", "p"])

    v = cph.variance_matrix_
    var_g = float(v.loc["gene_high", "gene_high"])
    var_i = float(v.loc["gene_high_status", "gene_high_status"])
    cov_gi = float(v.loc["gene_high", "gene_high_status"])

    def ci_from_beta(beta: float, var: float) -> tuple[float, float]:
        se = float(np.sqrt(max(var, 0.0)))
        lo = beta - 1.96 * se
        hi = beta + 1.96 * se
        return float(np.exp(lo)), float(np.exp(hi))

    # lo_lo: beta_g
    hr_lo = float(np.exp(beta_g))
    ci_lo_l, ci_lo_u = ci_from_beta(beta_g, var_g)

    # hi_hi: beta_g + beta_i
    beta_hi = beta_g + beta_i
    var_hi = var_g + var_i + 2.0 * cov_gi
    hr_hi = float(np.exp(beta_hi))
    ci_hi_l, ci_hi_u = ci_from_beta(beta_hi, var_hi)

    # Vulnerability Index from interaction: exp(beta_interaction)
    vi = float(np.exp(beta_i))
    vi_l, vi_u = ci_from_beta(beta_i, var_i)

    return {
        "Gene": gene,
        "n": int(df.shape[0]),
        "events": int(df["OS_event"].sum()),
        "HR_lo_lo": hr_lo,
        "CI_lo_lo_low": ci_lo_l,
        "CI_lo_lo_high": ci_lo_u,
        "HR_hi_hi": hr_hi,
        "CI_hi_hi_low": ci_hi_l,
        "CI_hi_hi_high": ci_hi_u,
        "VI": vi,
        "VI_low": vi_l,
        "VI_high": vi_u,
        "beta_interaction": beta_i,
        "p_interaction": p_i,
    }


def plot_forest(results: pd.DataFrame, out_pdf: Path, out_png: Path) -> None:
    plt.rcParams["figure.dpi"] = 300
    plt.rcParams["savefig.dpi"] = 300
    plt.rcParams["axes.unicode_minus"] = False

    n_genes = results.shape[0]
    fig_h = max(3.5, 0.55 * n_genes)
    fig_w = 11

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.2, 2.6, 2.0], wspace=0.02)

    ax_gene = fig.add_subplot(gs[0, 0])
    ax_plot = fig.add_subplot(gs[0, 1])
    ax_txt = fig.add_subplot(gs[0, 2])

    y = np.arange(n_genes)[::-1]

    # left: gene names
    ax_gene.axis("off")
    for i, gene in enumerate(results["Gene"].tolist()):
        ax_gene.text(0.0, y[i], gene, ha="left", va="center", fontsize=12, fontweight="bold")

    # middle: forest plot (two groups)
    ax_plot.set_xscale("log")
    ax_plot.axvline(1.0, color="black", linestyle="--", linewidth=1.0)

    color_lo = "#377eb8"  # lo_lo
    color_hi = "#D95F02"  # hi_hi
    dy = 0.15

    for i, row in results.iterrows():
        yi = y[i]
        # lo_lo
        ax_plot.plot([row["CI_lo_lo_low"], row["CI_lo_lo_high"]], [yi - dy, yi - dy], color=color_lo, lw=2)
        ax_plot.scatter(row["HR_lo_lo"], yi - dy, color=color_lo, marker="s", s=40, edgecolor="white", lw=0.8, zorder=3)
        # hi_hi
        ax_plot.plot([row["CI_hi_hi_low"], row["CI_hi_hi_high"]], [yi + dy, yi + dy], color=color_hi, lw=2)
        ax_plot.scatter(row["HR_hi_hi"], yi + dy, color=color_hi, marker="s", s=40, edgecolor="white", lw=0.8, zorder=3)

    ax_plot.set_yticks(y)
    ax_plot.set_yticklabels([""] * n_genes)
    ax_plot.set_xlabel("Hazard Ratio (High vs Low within status)", fontsize=11, fontweight="bold")

    # dynamic xlim
    x_min = float(np.nanmin([results["CI_lo_lo_low"].min(), results["CI_hi_hi_low"].min(), 0.6]))
    x_max = float(np.nanmax([results["CI_lo_lo_high"].max(), results["CI_hi_hi_high"].max(), 3.0]))
    ax_plot.set_xlim(max(0.2, x_min * 0.9), x_max * 1.08)

    # right: text columns
    ax_txt.axis("off")
    ax_txt.text(0.00, 1.02, "HR (lo_lo)", transform=ax_txt.transAxes, fontsize=11, fontweight="bold", ha="left")
    ax_txt.text(0.40, 1.02, "HR (hi_hi)", transform=ax_txt.transAxes, fontsize=11, fontweight="bold", ha="left")
    ax_txt.text(0.78, 1.02, "VI / P(int)", transform=ax_txt.transAxes, fontsize=11, fontweight="bold", ha="left")

    def fmt_hr(hr: float, lo: float, hi: float) -> str:
        return f"{hr:.2f} ({lo:.2f}, {hi:.2f})"

    def fmt_p(p: float) -> str:
        if p < 1e-3:
            return "<0.001"
        return f"{p:.3f}"

    for i, row in results.iterrows():
        yi = y[i]
        ax_txt.text(0.00, yi, fmt_hr(row["HR_lo_lo"], row["CI_lo_lo_low"], row["CI_lo_lo_high"]), ha="left", va="center", fontsize=10)
        ax_txt.text(0.40, yi, fmt_hr(row["HR_hi_hi"], row["CI_hi_hi_low"], row["CI_hi_hi_high"]), ha="left", va="center", fontsize=10)
        ax_txt.text(0.78, yi, f"{row['VI']:.2f} / {fmt_p(row['p_interaction'])}", ha="left", va="center", fontsize=10, fontweight="bold" if row["q_interaction"] < 0.05 else "normal")

    # legend
    ax_plot.scatter([], [], color=color_lo, marker="s", s=40, label="lo_lo")
    ax_plot.scatter([], [], color=color_hi, marker="s", s=40, label="hi_hi")
    ax_plot.legend(loc="lower right", frameon=False, fontsize=10)

    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    fig.savefig(out_png, bbox_inches="tight", facecolor="white")


def main() -> None:
    args = parse_args()
    base_dir = args.base_dir.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    expr_file = (args.expr_file.resolve() if args.expr_file else (base_dir / "combined_expression_combat_corrected.txt"))
    survival_file = (args.survival_file.resolve() if args.survival_file else _default_survival_file(base_dir))
    annotation_file = (args.annotation_file.resolve() if args.annotation_file else (base_dir / "MYC_PVT1_annotation.txt"))

    out_pdf = out_dir / "forest_plot_like.pdf"
    out_png = out_dir / f"forest_plot_like_{int(time.time())}.png"
    out_tsv = out_dir / "forest_plot_like_interaction_stats.tsv"

    inputs = load_inputs(expr_file=expr_file, survival_file=survival_file, annotation_file=annotation_file)

    rows = []
    for g in GENES:
        rows.append(fit_interaction_for_gene(inputs, g))

    res = pd.DataFrame(rows)
    res["q_interaction"] = bh_fdr(res["p_interaction"].values)
    res = res.sort_values(["q_interaction", "p_interaction"], ascending=[True, True]).reset_index(drop=True)

    res.to_csv(out_tsv, sep="\t", index=False)
    plot_forest(res, out_pdf=out_pdf, out_png=out_png)

    print("✓ forest_plot_like regenerated")
    print("  - PDF:", out_pdf)
    print("  - PNG:", out_png)
    print("  - TSV:", out_tsv)


if __name__ == "__main__":
    main()
