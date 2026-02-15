#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Interaction-term Cox test for a small gene list, using within-status high/low split.

This matches the HRR figure logic more closely:
- define gene_high within each status group (hi_hi / lo_lo) by median-rank (top half vs bottom half)
- fit a single interaction Cox model:
    Surv ~ gene_high + status + gene_high:status
- report interaction p-value.

Notes
- status is MYC/PVT1 proxy: lo_lo=0, hi_hi=1
- stratify by Dataset when available
"""

from __future__ import annotations

import os
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

BASE_DIR = r"E:\data\changyuan\免疫队列\单基因生存分析"
EXPR_FILE = os.path.join(BASE_DIR, "combined_expression_combat_corrected.txt")
ANNOT_FILE = os.path.join(BASE_DIR, "MYC_PVT1_annotation.txt")
SURV_FILE = os.path.join(os.path.dirname(BASE_DIR), "生存曲线", "updated_survival_data.txt")

GENES = [
    "FOLH1",
    "ZMYND8",
    "F2RL1",
    "ARG1",
    "ZHX1",
    "SOX9",
    "HAPLN1",
    "CDK6",
]

OUT_TSV = os.path.join(BASE_DIR, "cox_interaction_hrr_genes_hilo_within_status.tsv")

MIN_SAMPLES = 40
MIN_EVENTS = 8


def load_data() -> Tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    expr = pd.read_csv(EXPR_FILE, sep="\t", index_col=0)

    annot = pd.read_csv(ANNOT_FILE, sep="\t")
    status = pd.Series(annot["MYC_PVT1_Status"].values, index=annot["Sample"].astype(str).values)
    status = status[status.isin(["hi_hi", "lo_lo"])]

    surv = pd.read_csv(SURV_FILE, sep=r"\s+", engine="python").set_index("Sample_ID")

    common = sorted(set(expr.columns) & set(status.index) & set(surv.index))
    base = pd.DataFrame(index=common)
    base["OS_time"] = pd.to_numeric(surv.loc[common, "OS_months"], errors="coerce")
    base["OS_event"] = pd.to_numeric(surv.loc[common, "OS_event"], errors="coerce")
    base["status"] = (status.loc[common] == "hi_hi").astype(int).values

    if "Dataset" in surv.columns:
        base["Dataset"] = surv.loc[common, "Dataset"].astype(str).values

    base = base.dropna(subset=["OS_time", "OS_event", "status"]).copy()
    base = base[(base["OS_time"] > 0) & (base["OS_event"].isin([0, 1]))].copy()

    expr = expr.loc[:, base.index.tolist()]
    return expr, status.loc[base.index], base


def median_rank_split(values: pd.Series) -> pd.Series:
    # even split by rank: bottom half 0, top half 1
    s = values.sort_values(kind="mergesort")
    half = len(s) // 2
    low_idx = s.index[:half]
    high_idx = s.index[half:]
    out = pd.Series(0, index=s.index, dtype=int)
    out.loc[high_idx] = 1
    return out.reindex(values.index)


def fit_gene(expr: pd.DataFrame, base: pd.DataFrame, gene: str) -> Dict[str, object]:
    row: Dict[str, object] = {"Gene": gene}

    if gene not in expr.index:
        row.update({"ok": False, "reason": "gene_not_found"})
        return row

    x = expr.loc[gene]

    df = base.copy()
    df["expr"] = pd.to_numeric(x, errors="coerce")
    df = df.dropna(subset=["expr"]).copy()

    if df.shape[0] < MIN_SAMPLES or int(df["OS_event"].sum()) < MIN_EVENTS:
        row.update({"ok": False, "reason": "too_few_samples_or_events", "n": int(df.shape[0]), "events": int(df["OS_event"].sum())})
        return row

    # define gene_high within each status group
    gene_high = pd.Series(index=df.index, dtype=int)
    for s in [0, 1]:
        idx = df.index[df["status"] == s]
        if len(idx) < 10:
            row.update({"ok": False, "reason": "too_few_in_status_group"})
            return row
        gene_high.loc[idx] = median_rank_split(df.loc[idx, "expr"]).astype(int)

    df["gene_high"] = gene_high.astype(int)
    df["interaction"] = df["gene_high"] * df["status"]

    cols = ["OS_time", "OS_event", "gene_high", "status", "interaction"]
    fit_df = df[cols + (["Dataset"] if "Dataset" in df.columns else [])].copy()

    cph = CoxPHFitter()
    try:
        if "Dataset" in fit_df.columns:
            cph.fit(fit_df, duration_col="OS_time", event_col="OS_event", strata=["Dataset"], robust=True)
        else:
            cph.fit(fit_df, duration_col="OS_time", event_col="OS_event", robust=True)
    except Exception as e:
        row.update({"ok": False, "reason": f"fit_failed: {type(e).__name__}"})
        return row

    if "interaction" not in cph.summary.index or "gene_high" not in cph.summary.index:
        row.update({"ok": False, "reason": "missing_terms"})
        return row

    beta_g = float(cph.params_.loc["gene_high"])
    beta_i = float(cph.params_.loc["interaction"])

    row.update(
        {
            "ok": True,
            "n": int(fit_df.shape[0]),
            "events": int(fit_df["OS_event"].sum()),
            "p_interaction": float(cph.summary.loc["interaction", "p"]),
            "beta_gene": beta_g,
            "beta_interaction": beta_i,
            "HR_lo_lo_high_vs_low": float(np.exp(beta_g)),
            "HR_hi_hi_high_vs_low": float(np.exp(beta_g + beta_i)),
            "VI_from_interaction": float(np.exp(beta_i)),  # reviewer-acceptable replacement for HRR
        }
    )
    return row


def main() -> None:
    expr, _status_series, base = load_data()

    rows = [fit_gene(expr, base, g) for g in GENES]
    out = pd.DataFrame(rows)

    out.to_csv(OUT_TSV, sep="\t", index=False)
    print("Saved:", OUT_TSV)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
