#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot survival curves for selected genes in MYC/PVT1 hi_hi vs lo_lo groups.

Reproduces the 4-curve style shown in the example figure:
  - hi_hi & gene_high / gene_low (solid)
  - lo_lo & gene_high / gene_low (dashed)
and annotates per-group Cox p-values (high vs low within each status).

Split strategy: median-rank within each status group (equal halves).
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test


BASE_DIR = r"E:/data/changyuan/免疫队列/单基因生存分析"
EXPR_FILE = os.path.join(BASE_DIR, "combined_expression_combat_corrected.txt")
ANNOT_FILE = os.path.join(BASE_DIR, "MYC_PVT1_annotation.txt")
SURV_FILE = os.path.join(os.path.dirname(BASE_DIR), "生存曲线", "updated_survival_data.txt")

OUT_DIR = os.path.join(BASE_DIR, "survival_plots_median_rank")

TARGET_GENES = ["AGT", "APOA2", "APOC1", "FGFR4", "LGALS4", "TRIM6", "HAPLN1", "ARG1"]


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
class GroupStats:
    n_high: int
    n_low: int
    p_value: float


def _median_rank_split(values: pd.Series) -> pd.Series:
    s = values.sort_values(kind="mergesort")
    half = len(s) // 2
    out = pd.Series(0, index=s.index, dtype=int)
    out.loc[s.index[half:]] = 1
    return out.reindex(values.index)


def _prepare_gene_df(expr: pd.DataFrame, status: pd.Series, surv: pd.DataFrame, gene: str, group: str) -> pd.DataFrame:
    if gene not in expr.index:
        raise KeyError(f"Gene not found in expression matrix: {gene}")

    common = sorted(set(expr.columns) & set(status.index) & set(surv.index))
    df = pd.DataFrame(index=common)
    df["status"] = status.loc[common].astype(str).values
    df["OS_time"] = pd.to_numeric(surv.loc[common, "OS_months"], errors="coerce")
    df["OS_event"] = pd.to_numeric(surv.loc[common, "OS_event"], errors="coerce")
    df["expr"] = pd.to_numeric(expr.loc[gene, common], errors="coerce").values

    df = df.dropna(subset=["status", "OS_time", "OS_event", "expr"]).copy()
    df = df[(df["OS_time"] > 0) & (df["OS_event"].isin([0, 1]))].copy()
    df = df[df["status"] == group].copy()
    return df


def _fit_group_pvalue(df: pd.DataFrame) -> Tuple[pd.DataFrame, GroupStats]:
    if df.empty:
        return df, GroupStats(n_high=0, n_low=0, p_value=float("nan"))

    df = df.copy()
    df["Expression_Group"] = _median_rank_split(df["expr"]).astype(int)  # 1=high, 0=low
    n_high = int((df["Expression_Group"] == 1).sum())
    n_low = int((df["Expression_Group"] == 0).sum())

    # Log-rank p-value for high vs low within this group (matches example figure)
    try:
        high = df[df["Expression_Group"] == 1]
        low = df[df["Expression_Group"] == 0]
        lr = logrank_test(
            durations_A=high["OS_time"],
            durations_B=low["OS_time"],
            event_observed_A=high["OS_event"],
            event_observed_B=low["OS_event"],
        )
        p_value = float(lr.p_value)
    except Exception:
        p_value = float("nan")

    return df, GroupStats(n_high=n_high, n_low=n_low, p_value=p_value)


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

    # Full frame, thick spines
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

    # bracket in axes coordinates
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


def plot_one_gene(ax: plt.Axes, expr: pd.DataFrame, status: pd.Series, surv: pd.DataFrame, gene: str) -> Dict[str, GroupStats]:
    style = STYLE_GRID if getattr(ax, "_compact", False) else STYLE_SINGLE

    hihi_raw = _prepare_gene_df(expr, status, surv, gene, "hi_hi")
    lolo_raw = _prepare_gene_df(expr, status, surv, gene, "lo_lo")

    hihi_df, hihi_stats = _fit_group_pvalue(hihi_raw)
    lolo_df, lolo_stats = _fit_group_pvalue(lolo_raw)

    kmf = KaplanMeierFitter()

    # hi_hi curves (solid)
    for grp_val, key in [(1, "hihi_high"), (0, "hihi_low")]:
        sub = hihi_df[hihi_df["Expression_Group"] == grp_val]
        label = f"hi_hi & {gene}_{'high' if grp_val == 1 else 'low'} (n={len(sub)})"
        if len(sub) > 0:
            kmf.fit(sub["OS_time"], sub["OS_event"], label=label)
            kmf.plot_survival_function(ax=ax, ci_show=False, color=COLORS[key], lw=style["curve_lw"], ls="-")

    # lo_lo curves (dashed)
    for grp_val, key in [(1, "lolo_high"), (0, "lolo_low")]:
        sub = lolo_df[lolo_df["Expression_Group"] == grp_val]
        label = f"lo_lo & {gene}_{'high' if grp_val == 1 else 'low'} (n={len(sub)})"
        if len(sub) > 0:
            kmf.fit(sub["OS_time"], sub["OS_event"], label=label)
            kmf.plot_survival_function(
                ax=ax,
                ci_show=False,
                color=COLORS[key],
                lw=style["curve_lw"],
                ls=(0, (10, 7)),
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

    # Brackets + p-values aligned to legend entries (fixed positions)
    _draw_p_bracket(ax, x=style["p_x"], y_top=style["p1_y"][0], y_bottom=style["p1_y"][1], p=hihi_stats.p_value, fontsize=style["p_fs"])
    _draw_p_bracket(ax, x=style["p_x"], y_top=style["p2_y"][0], y_bottom=style["p2_y"][1], p=lolo_stats.p_value, fontsize=style["p_fs"])

    return {"hi_hi": hihi_stats, "lo_lo": lolo_stats}


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


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)

    plt.rcParams["font.family"] = "Arial"

    expr, status, surv = load_inputs()

    # Individual plots
    for gene in TARGET_GENES:
        fig, ax = plt.subplots(figsize=(12.5, 8.5))
        plot_one_gene(ax, expr, status, surv, gene)
        fig.tight_layout()
        out_path = os.path.join(OUT_DIR, f"{gene}_km_median_rank.png")
        fig.savefig(out_path, dpi=300)
        plt.close(fig)

    # Combined grid (2x4) for quick viewing
    fig, axes = plt.subplots(2, 4, figsize=(22, 14))
    for ax, gene in zip(axes.ravel(), TARGET_GENES):
        setattr(ax, "_compact", True)
        plot_one_gene(ax, expr, status, surv, gene)
    fig.tight_layout()
    grid_path = os.path.join(OUT_DIR, "km_8genes_median_rank_grid.png")
    fig.savefig(grid_path, dpi=300)
    plt.close(fig)

    print(f"Saved 8 individual plots + grid to: {OUT_DIR}")


if __name__ == "__main__":
    main()
