#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Continuous Cox predicted survival curves (4-curve style).

Goal: match the look of `plot_survival.py` (KM + log-rank), but curves are
predicted from a continuous Cox model:

  Surv(OS_time, OS_event) ~ gene_z + status + gene_z:status

Where:
  - gene_z is z-scored expression across the full matched cohort (per gene)
  - status: lo_lo=0, hi_hi=1

To keep the same 4-curve layout and n= labels, we still define "high/low"
within each status by median-rank split (equal halves), but the curves are
Cox-predicted at the mean gene_z of each half.

P-values shown are Wald tests for the continuous gene effect within each
status group:
  - lo_lo: p(gene_z)
  - hi_hi: p(gene_z + interaction) via linear-combination SE
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter


BASE_DIR = r"E:/data/changyuan/免疫队列/单基因生存分析"
EXPR_FILE = os.path.join(BASE_DIR, "combined_expression_combat_corrected.txt")
ANNOT_FILE = os.path.join(BASE_DIR, "MYC_PVT1_annotation.txt")
SURV_FILE = os.path.join(os.path.dirname(BASE_DIR), "生存曲线", "updated_survival_data.txt")

OUT_DIR = os.path.join(BASE_DIR, "survival_plots_cox_pred")

TARGET_GENES = ["AGT", "APOA2", "APOC1", "FGFR4", "LGALS4", "TRIM6", "OASL", "ARG1"]

COLORS = {
    "hihi_high": "#7f0000",  # dark red
    "hihi_low": "#00008b",  # dark blue
    "lolo_high": "#ffa500",  # orange
    "lolo_low": "#87cefa",  # light blue
}

STYLE_SINGLE = {
    "title_fs": 44,
    "label_fs": 24,
    "tick_fs": 18,
    "spine_lw": 4.0,
    "tick_w": 3.0,
    "tick_l": 8,
    "curve_lw": 6.0,
    "legend_fs": 22,
    "p_fs": 22,
    "legend_anchor": (0.34, 0.90),
    "p_x": 0.79,
    "p1_y": (0.83, 0.75),
    "p2_y": (0.67, 0.59),
}

STYLE_GRID = {
    "title_fs": 24,
    "label_fs": 14,
    "tick_fs": 11,
    "spine_lw": 2.4,
    "tick_w": 1.8,
    "tick_l": 5,
    "curve_lw": 3.2,
    "legend_fs": 11,
    "p_fs": 11,
    "legend_anchor": (0.36, 0.90),
    "p_x": 0.79,
    "p1_y": (0.83, 0.75),
    "p2_y": (0.67, 0.59),
}


@dataclass(frozen=True)
class GroupDef:
    n_high: int
    n_low: int
    z_high: float
    z_low: float
    dataset_weights: Dict[str, float]


def _median_rank_split(values: pd.Series) -> pd.Series:
    s = values.sort_values(kind="mergesort")
    half = len(s) // 2
    out = pd.Series(0, index=s.index, dtype=int)
    out.loc[s.index[half:]] = 1
    return out.reindex(values.index)


def _p_from_z(z: float) -> float:
    if not np.isfinite(z):
        return float("nan")
    return float(math.erfc(abs(z) / math.sqrt(2.0)))


def _style_axes(ax: plt.Axes, title: str, style: Dict[str, float]) -> None:
    ax.set_title(title, fontsize=style["title_fs"], fontweight="bold", pad=12)
    ax.set_xlabel("Overall survival time (Months)", fontsize=style["label_fs"], fontweight="bold")
    ax.set_ylabel("Probability of Survival", fontsize=style["label_fs"], fontweight="bold")

    ax.set_ylim(0, 1.02)
    ax.tick_params(
        axis="both",
        which="both",
        labelsize=style["tick_fs"],
        width=style["tick_w"],
        length=style["tick_l"],
    )

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(style["spine_lw"])
        spine.set_color("black")

    ax.grid(False)


def _draw_p_bracket(ax: plt.Axes, x: float, y_top: float, y_bottom: float, p: float, fontsize: int) -> None:
    if not np.isfinite(p):
        p_txt = "NA"
    else:
        p_txt = f"{p:.4f}"

    ax.plot([x, x], [y_bottom, y_top], transform=ax.transAxes, color="black", lw=2.0, clip_on=False)
    ax.plot([x - 0.02, x], [y_top, y_top], transform=ax.transAxes, color="black", lw=2.0, clip_on=False)
    ax.plot([x - 0.02, x], [y_bottom, y_bottom], transform=ax.transAxes, color="black", lw=2.0, clip_on=False)
    ax.text(
        x + 0.01,
        (y_top + y_bottom) / 2,
        f"p = {p_txt}",
        transform=ax.transAxes,
        fontsize=fontsize,
        fontweight="bold",
        va="center",
        ha="left",
    )


def load_inputs() -> Tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    if not os.path.exists(EXPR_FILE):
        raise FileNotFoundError(EXPR_FILE)
    if not os.path.exists(ANNOT_FILE):
        raise FileNotFoundError(ANNOT_FILE)
    if not os.path.exists(SURV_FILE):
        raise FileNotFoundError(SURV_FILE)

    expr = pd.read_csv(EXPR_FILE, sep="\t", index_col=0)
    annot = pd.read_csv(ANNOT_FILE, sep="\t")
    status = pd.Series(annot["MYC_PVT1_Status"].values, index=annot["Sample"].astype(str).values)
    status = status[status.isin(["hi_hi", "lo_lo"])]

    surv = pd.read_csv(SURV_FILE, sep=r"\s+", engine="python").set_index("Sample_ID")
    return expr, status, surv


def build_base_df(expr: pd.DataFrame, status: pd.Series, surv: pd.DataFrame, gene: str) -> pd.DataFrame:
    if gene not in expr.index:
        raise KeyError(f"Gene not found in expression matrix: {gene}")

    common = sorted(set(expr.columns) & set(status.index) & set(surv.index))
    df = pd.DataFrame(index=common)
    df["status_str"] = status.loc[common].astype(str).values
    df["status"] = (df["status_str"] == "hi_hi").astype(int)
    df["OS_time"] = pd.to_numeric(surv.loc[common, "OS_months"], errors="coerce")
    df["OS_event"] = pd.to_numeric(surv.loc[common, "OS_event"], errors="coerce")
    if "Dataset" in surv.columns:
        df["Dataset"] = surv.loc[common, "Dataset"].astype(str).values
    else:
        df["Dataset"] = "all"
    df["expr"] = pd.to_numeric(expr.loc[gene, common], errors="coerce").values

    df = df.dropna(subset=["status", "OS_time", "OS_event", "Dataset", "expr"]).copy()
    df = df[(df["OS_time"] > 0) & (df["OS_event"].isin([0, 1]))].copy()

    x = df["expr"].to_numpy(dtype=float)
    df["gene_z"] = (x - np.nanmean(x)) / np.nanstd(x)
    df["interaction"] = df["gene_z"] * df["status"]
    return df


def define_groups(df: pd.DataFrame, status_value: int) -> GroupDef:
    sub = df[df["status"] == status_value].copy()
    if sub.empty:
        return GroupDef(0, 0, float("nan"), float("nan"), {})

    sub["Expression_Group"] = _median_rank_split(sub["expr"]).astype(int)
    high = sub[sub["Expression_Group"] == 1]
    low = sub[sub["Expression_Group"] == 0]

    # representative z values for predictions
    z_high = float(np.nanmean(high["gene_z"])) if len(high) else float("nan")
    z_low = float(np.nanmean(low["gene_z"])) if len(low) else float("nan")

    # dataset mixture weights within this status group
    counts = sub["Dataset"].value_counts().to_dict()
    total = float(sum(counts.values())) if counts else 0.0
    weights = {k: (v / total) for k, v in counts.items()} if total > 0 else {}

    return GroupDef(
        n_high=int(len(high)),
        n_low=int(len(low)),
        z_high=z_high,
        z_low=z_low,
        dataset_weights=weights,
    )


def _subgroup_df(df: pd.DataFrame, status_value: int, is_high: bool) -> pd.DataFrame:
    sub = df[df["status"] == status_value].copy()
    if sub.empty:
        return sub
    sub["Expression_Group"] = _median_rank_split(sub["expr"]).astype(int)
    return sub[sub["Expression_Group"] == (1 if is_high else 0)].copy()


def _plot_censor_marks(ax: plt.Axes, surv: pd.Series, censor_times: np.ndarray, color: str, style: Dict[str, float]) -> None:
    """
    Plot censor tick-marks at the predicted survival probability at censor time.
    """
    if censor_times.size == 0:
        return

    # Stepwise lookup: y(t) = S(t_last <= t)
    idx = surv.index.to_numpy(dtype=float)
    y = surv.to_numpy(dtype=float)
    order = np.argsort(idx)
    idx = idx[order]
    y = y[order]

    ct = np.asarray(censor_times, dtype=float)
    ct = ct[np.isfinite(ct)]
    if ct.size == 0:
        return

    # For each censor time, find rightmost index <= t.
    pos = np.searchsorted(idx, ct, side="right") - 1
    pos = np.clip(pos, 0, len(idx) - 1)
    cy = y[pos]

    ax.plot(
        ct,
        cy,
        linestyle="None",
        marker="|",
        markersize=max(8.0, style["curve_lw"] * 2.2),
        markeredgewidth=max(1.2, style["curve_lw"] * 0.35),
        color=color,
        alpha=0.9,
        zorder=5,
    )


def _baseline_survival_by_stratum(cph: CoxPHFitter) -> Dict[str, pd.Series]:
    """
    lifelines baseline_survival_ for stratified models is typically a DataFrame
    with one column per stratum. Return {stratum_value: Series(time->S0)}.
    """
    bs = cph.baseline_survival_
    if isinstance(bs, pd.Series):
        return {"all": bs}
    if not isinstance(bs, pd.DataFrame) or bs.shape[1] == 0:
        raise RuntimeError("Unexpected baseline_survival_ shape.")

    out: Dict[str, pd.Series] = {}
    for col in bs.columns:
        out[str(col)] = bs[col].copy()
    return out


def predict_mixture_survival(
    baseline_by_dataset: Dict[str, pd.Series],
    dataset_weights: Dict[str, float],
    linpred: float,
) -> pd.Series:
    """
    Compute weighted-average survival across dataset strata:
      S_mix(t) = sum_w S0_d(t) ^ exp(linpred)
    """
    if not dataset_weights:
        # fallback: equal weights across available baselines
        keys = list(baseline_by_dataset.keys())
        dataset_weights = {k: 1.0 / len(keys) for k in keys} if keys else {"all": 1.0}

    # common time grid
    all_times = sorted(set().union(*[s.index.tolist() for s in baseline_by_dataset.values()]))
    idx = pd.Index(all_times, name="timeline")

    out = pd.Series(0.0, index=idx)
    scale = float(np.exp(linpred))
    for ds, w in dataset_weights.items():
        if ds not in baseline_by_dataset:
            continue
        s0 = baseline_by_dataset[ds].reindex(idx).ffill().bfill()
        out += float(w) * (s0 ** scale)
    return out.clip(lower=0.0, upper=1.0)


def fit_interaction_cox(df: pd.DataFrame) -> CoxPHFitter:
    fit_df = df[["OS_time", "OS_event", "gene_z", "status", "interaction", "Dataset"]].copy()
    cph = CoxPHFitter()
    cph.fit(
        fit_df,
        duration_col="OS_time",
        event_col="OS_event",
        strata=["Dataset"],
        robust=True,
    )
    return cph


def gene_effect_pvalues(cph: CoxPHFitter) -> Tuple[float, float]:
    """
    Return (p_lo_lo, p_hi_hi) for gene_z effect within each status group.
    """
    beta_g = float(cph.params_.loc["gene_z"])
    beta_i = float(cph.params_.loc["interaction"])
    p_lo = float(cph.summary.loc["gene_z", "p"])

    v = (
        cph.variance_matrix_
        .loc[["gene_z", "interaction"], ["gene_z", "interaction"]]
        .to_numpy(dtype=float)
    )
    var_gene = float(v[0, 0])
    var_int = float(v[1, 1])
    cov = float(v[0, 1])
    var_sum = var_gene + var_int + 2.0 * cov
    se = float(np.sqrt(max(var_sum, 0.0)))
    p_hi = _p_from_z((beta_g + beta_i) / se) if se > 0 else float("nan")
    return p_lo, p_hi


def plot_one_gene(
    ax: plt.Axes,
    expr: pd.DataFrame,
    status: pd.Series,
    surv: pd.DataFrame,
    gene: str,
) -> None:
    style = STYLE_GRID if getattr(ax, "_compact", False) else STYLE_SINGLE

    df = build_base_df(expr, status, surv, gene)
    cph = fit_interaction_cox(df)

    p_lo, p_hi = gene_effect_pvalues(cph)
    baseline = _baseline_survival_by_stratum(cph)

    # define display groups (for labels + representative z values)
    hihi = define_groups(df, 1)
    lolo = define_groups(df, 0)
    hihi_high_df = _subgroup_df(df, 1, True)
    hihi_low_df = _subgroup_df(df, 1, False)
    lolo_high_df = _subgroup_df(df, 0, True)
    lolo_low_df = _subgroup_df(df, 0, False)

    beta_g = float(cph.params_.loc["gene_z"])
    beta_s = float(cph.params_.loc["status"])
    beta_i = float(cph.params_.loc["interaction"])

    def linpred(status_val: int, z_val: float) -> float:
        return beta_g * z_val + beta_s * status_val + beta_i * z_val * status_val

    # Predicted curves
    s_hihi_high = predict_mixture_survival(baseline, hihi.dataset_weights, linpred(1, hihi.z_high))
    s_hihi_low = predict_mixture_survival(baseline, hihi.dataset_weights, linpred(1, hihi.z_low))
    s_lolo_high = predict_mixture_survival(baseline, lolo.dataset_weights, linpred(0, lolo.z_high))
    s_lolo_low = predict_mixture_survival(baseline, lolo.dataset_weights, linpred(0, lolo.z_low))

    # Plot
    ax.plot(s_hihi_high.index, s_hihi_high.values, color=COLORS["hihi_high"], lw=style["curve_lw"], ls="-",
            label=f"hi_hi & {gene}_high (n={hihi.n_high})")
    ax.plot(s_hihi_low.index, s_hihi_low.values, color=COLORS["hihi_low"], lw=style["curve_lw"], ls="-",
            label=f"hi_hi & {gene}_low (n={hihi.n_low})")

    dash = (0, (10, 7))
    ax.plot(s_lolo_high.index, s_lolo_high.values, color=COLORS["lolo_high"], lw=style["curve_lw"], ls=dash,
            label=f"lo_lo & {gene}_high (n={lolo.n_high})")
    ax.plot(s_lolo_low.index, s_lolo_low.values, color=COLORS["lolo_low"], lw=style["curve_lw"], ls=dash,
            label=f"lo_lo & {gene}_low (n={lolo.n_low})")

    # Censor marks (OS_event == 0) for each subgroup
    _plot_censor_marks(
        ax,
        s_hihi_high,
        hihi_high_df.loc[hihi_high_df["OS_event"] == 0, "OS_time"].to_numpy(),
        COLORS["hihi_high"],
        style,
    )
    _plot_censor_marks(
        ax,
        s_hihi_low,
        hihi_low_df.loc[hihi_low_df["OS_event"] == 0, "OS_time"].to_numpy(),
        COLORS["hihi_low"],
        style,
    )
    _plot_censor_marks(
        ax,
        s_lolo_high,
        lolo_high_df.loc[lolo_high_df["OS_event"] == 0, "OS_time"].to_numpy(),
        COLORS["lolo_high"],
        style,
    )
    _plot_censor_marks(
        ax,
        s_lolo_low,
        lolo_low_df.loc[lolo_low_df["OS_event"] == 0, "OS_time"].to_numpy(),
        COLORS["lolo_low"],
        style,
    )

    _style_axes(ax, gene, style)

    legend = ax.legend(
        loc="upper left",
        bbox_to_anchor=style["legend_anchor"],
        frameon=False,
        fontsize=style["legend_fs"],
        handlelength=3.0,
    )
    for text in legend.get_texts():
        text.set_fontweight("bold")

    # p-values: hi_hi (continuous effect), lo_lo (continuous effect)
    _draw_p_bracket(ax, x=style["p_x"], y_top=style["p1_y"][0], y_bottom=style["p1_y"][1], p=p_hi, fontsize=style["p_fs"])
    _draw_p_bracket(ax, x=style["p_x"], y_top=style["p2_y"][0], y_bottom=style["p2_y"][1], p=p_lo, fontsize=style["p_fs"])


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    plt.rcParams["font.family"] = "Arial"
    # Keep text editable in Illustrator (embed TrueType)
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    expr, status, surv = load_inputs()

    # Individual plots
    for gene in TARGET_GENES:
        fig, ax = plt.subplots(figsize=(12.5, 8.5))
        plot_one_gene(ax, expr, status, surv, gene)
        fig.tight_layout()
        out_png = os.path.join(OUT_DIR, f"{gene}_cox_pred.png")
        out_pdf = os.path.join(OUT_DIR, f"{gene}_cox_pred.pdf")
        fig.savefig(out_png, dpi=300)
        fig.savefig(out_pdf)
        plt.close(fig)

    # Combined grid (2x4)
    fig, axes = plt.subplots(2, 4, figsize=(22, 14))
    for ax, gene in zip(axes.ravel(), TARGET_GENES):
        setattr(ax, "_compact", True)
        plot_one_gene(ax, expr, status, surv, gene)
    fig.tight_layout()
    grid_png = os.path.join(OUT_DIR, "cox_pred_8genes_grid.png")
    grid_pdf = os.path.join(OUT_DIR, "cox_pred_8genes_grid.pdf")
    fig.savefig(grid_png, dpi=300)
    fig.savefig(grid_pdf)
    plt.close(fig)

    print(f"Saved 8 individual Cox-predicted plots (PNG+PDF) + grid (PNG+PDF) to: {OUT_DIR}")


if __name__ == "__main__":
    main()
