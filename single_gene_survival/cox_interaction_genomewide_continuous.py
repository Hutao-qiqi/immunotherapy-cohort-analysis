#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Genome-wide continuous-expression Cox interaction scan
=====================================================

Model (reviewer-recommended):
    Surv(OS_time, OS_event) ~ gene_z + status + gene_z:status

Where:
  - gene_z: per-gene z-score across the full matched cohort
  - status: lo_lo=0, hi_hi=1 (MYC/PVT1 dual-expression groups)

This script extends `cox_interaction_term.py` by also reporting:
  - p_lo_lo: p-value for gene_z effect within lo_lo (beta_gene)
  - p_hi_hi: p-value for gene_z effect within hi_hi (beta_gene + beta_interaction)
  - VI: exp(beta_interaction)

It also outputs a hit list using the user's 3 criteria (raw p-values):
  (i)  p_interaction < 0.05
  (ii) within hi_hi: HR_hi_hi > 1 and p_hi_hi < 0.05
  (iii) within lo_lo: p_lo_lo >= 0.05
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from lifelines import CoxPHFitter
from tqdm import tqdm


BASE_DIR = r"E:\data\changyuan\免疫队列\单基因生存分析"
EXPR_FILE = os.path.join(BASE_DIR, "combined_expression_combat_corrected.txt")
ANNOT_FILE = os.path.join(BASE_DIR, "MYC_PVT1_annotation.txt")
SURV_FILE = os.path.join(os.path.dirname(BASE_DIR), "生存曲线", "updated_survival_data.txt")

OUT_ALL = os.path.join(BASE_DIR, "cox_interaction_continuous_genomewide.tsv")
OUT_HITS = os.path.join(BASE_DIR, "selected_continuous_3criteria_p0.05.tsv")

N_JOBS = min(8, os.cpu_count() or 1)
MIN_SAMPLES = 60
MIN_EVENTS = 10


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)

    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q

    p_ok = p[ok]
    order = np.argsort(p_ok)
    ranked = p_ok[order]

    m = float(len(ranked))
    q_ranked = ranked * m / (np.arange(1, len(ranked) + 1))
    q_ranked = np.minimum.accumulate(q_ranked[::-1])[::-1]
    q_ranked = np.clip(q_ranked, 0.0, 1.0)

    q_ok = np.empty_like(p_ok)
    q_ok[order] = q_ranked
    q[ok] = q_ok
    return q


@dataclass(frozen=True)
class InputData:
    samples: list[str]
    base_df: pd.DataFrame
    gene_names: list[str]
    expr_values: np.ndarray
    has_dataset: bool


def load_inputs() -> InputData:
    for p in (EXPR_FILE, ANNOT_FILE, SURV_FILE):
        if not os.path.exists(p):
            raise FileNotFoundError(p)

    expr = pd.read_csv(EXPR_FILE, sep="\t", index_col=0)

    annot = pd.read_csv(ANNOT_FILE, sep="\t")
    status_series = pd.Series(annot["MYC_PVT1_Status"].values, index=annot["Sample"].astype(str).values)
    status_series = status_series[status_series.isin(["hi_hi", "lo_lo"])]

    surv = pd.read_csv(SURV_FILE, sep=r"\s+", engine="python").set_index("Sample_ID")

    common = sorted(set(expr.columns) & set(status_series.index) & set(surv.index))
    if len(common) < MIN_SAMPLES:
        raise RuntimeError(f"Too few matched samples: {len(common)}")

    base_df = pd.DataFrame(index=common)
    base_df["OS_time"] = pd.to_numeric(surv.loc[common, "OS_months"], errors="coerce")
    base_df["OS_event"] = pd.to_numeric(surv.loc[common, "OS_event"], errors="coerce")
    base_df["status"] = (status_series.loc[common] == "hi_hi").astype(int).values

    has_dataset = "Dataset" in surv.columns
    if has_dataset:
        base_df["Dataset"] = surv.loc[common, "Dataset"].astype(str).values

    base_df = base_df.dropna(subset=["OS_time", "OS_event", "status"]).copy()
    base_df = base_df[(base_df["OS_time"] > 0) & (base_df["OS_event"].isin([0, 1]))].copy()
    if len(base_df) < MIN_SAMPLES:
        raise RuntimeError(f"Too few usable samples after QC: {len(base_df)}")
    if base_df["OS_event"].sum() < MIN_EVENTS:
        raise RuntimeError(f"Too few events after QC: {int(base_df['OS_event'].sum())}")

    samples = base_df.index.tolist()
    expr_sub = expr.loc[:, samples]
    expr_values = expr_sub.to_numpy(dtype=float, copy=True)
    gene_names = expr_sub.index.astype(str).tolist()

    return InputData(
        samples=samples,
        base_df=base_df,
        gene_names=gene_names,
        expr_values=expr_values,
        has_dataset=has_dataset,
    )


def _p_from_z(z: float) -> float:
    if not np.isfinite(z):
        return float("nan")
    return float(math.erfc(abs(z) / math.sqrt(2.0)))


def fit_one_gene(i: int, data: InputData) -> Optional[dict]:
    gene = data.gene_names[i]
    x = data.expr_values[i, :]

    x_mean = np.nanmean(x)
    x_std = np.nanstd(x)
    if not np.isfinite(x_std) or x_std == 0:
        return None

    x_z = (x - x_mean) / x_std

    df = data.base_df.copy()
    df["gene_z"] = x_z
    df = df.dropna(subset=["gene_z"]).copy()

    if len(df) < MIN_SAMPLES or df["OS_event"].sum() < MIN_EVENTS:
        return None

    df["interaction"] = df["gene_z"] * df["status"]

    cols = ["OS_time", "OS_event", "gene_z", "status", "interaction"]
    fit_df = df[cols + (["Dataset"] if data.has_dataset else [])].copy()

    cph = CoxPHFitter()
    try:
        if data.has_dataset:
            cph.fit(
                fit_df,
                duration_col="OS_time",
                event_col="OS_event",
                strata=["Dataset"],
                robust=True,
            )
        else:
            cph.fit(
                fit_df,
                duration_col="OS_time",
                event_col="OS_event",
                robust=True,
            )
    except Exception:
        return None

    if "interaction" not in cph.summary.index or "gene_z" not in cph.params_.index:
        return None

    beta_gene = float(cph.params_.loc["gene_z"])
    beta_int = float(cph.params_.loc["interaction"])
    p_int = float(cph.summary.loc["interaction", "p"])

    hr_lo = float(np.exp(beta_gene))
    hr_hi = float(np.exp(beta_gene + beta_int))
    vi = float(np.exp(beta_int))

    v = (
        cph.variance_matrix_
        .loc[["gene_z", "interaction"], ["gene_z", "interaction"]]
        .to_numpy(dtype=float)
    )
    var_gene = float(v[0, 0])
    var_int = float(v[1, 1])
    cov = float(v[0, 1])

    se_lo = float(np.sqrt(max(var_gene, 0.0)))
    se_hi = float(np.sqrt(max(var_gene + var_int + 2.0 * cov, 0.0)))

    p_lo = _p_from_z(beta_gene / se_lo) if se_lo > 0 else float("nan")
    p_hi = _p_from_z((beta_gene + beta_int) / se_hi) if se_hi > 0 else float("nan")

    # approximate 95% CI for HRs
    ci_lo_l = float(np.exp(beta_gene - 1.96 * se_lo)) if se_lo > 0 else float("nan")
    ci_lo_u = float(np.exp(beta_gene + 1.96 * se_lo)) if se_lo > 0 else float("nan")
    ci_hi_l = float(np.exp((beta_gene + beta_int) - 1.96 * se_hi)) if se_hi > 0 else float("nan")
    ci_hi_u = float(np.exp((beta_gene + beta_int) + 1.96 * se_hi)) if se_hi > 0 else float("nan")

    return {
        "Gene": gene,
        "n": int(len(fit_df)),
        "events": int(fit_df["OS_event"].sum()),
        "beta_gene": beta_gene,
        "beta_interaction": beta_int,
        "p_interaction": p_int,
        "q_interaction": np.nan,  # filled after aggregation
        "HR_gene_lo_lo": hr_lo,
        "p_lo_lo": p_lo,
        "CI95_lo_lo": f"{ci_lo_l:.3g}-{ci_lo_u:.3g}" if np.isfinite(ci_lo_l) and np.isfinite(ci_lo_u) else np.nan,
        "HR_gene_hi_hi": hr_hi,
        "p_hi_hi": p_hi,
        "CI95_hi_hi": f"{ci_hi_l:.3g}-{ci_hi_u:.3g}" if np.isfinite(ci_hi_l) and np.isfinite(ci_hi_u) else np.nan,
        "VI": vi,
    }


def main() -> None:
    print("Loading inputs...")
    data = load_inputs()
    print(f"Matched samples: {len(data.samples)}")
    print(f"Genes: {len(data.gene_names)}")
    print(f"Using dataset stratification: {data.has_dataset}")
    print(f"Parallel jobs: {N_JOBS}")

    results = Parallel(n_jobs=N_JOBS)(
        delayed(fit_one_gene)(i, data)
        for i in tqdm(range(len(data.gene_names)), desc="Cox interaction (continuous)")
    )

    rows = [r for r in results if r is not None]
    if not rows:
        raise RuntimeError("No valid model fits produced.")

    out = pd.DataFrame(rows)
    out["q_interaction"] = bh_fdr(out["p_interaction"].to_numpy())
    out = out.sort_values(["q_interaction", "p_interaction"], ascending=[True, True])
    out.to_csv(OUT_ALL, sep="\t", index=False)
    print(f"Saved: {OUT_ALL}")

    crit = (
        (out["p_interaction"] < 0.05)
        & (out["HR_gene_hi_hi"] > 1.0)
        & (out["p_hi_hi"] < 0.05)
        & (out["p_lo_lo"] >= 0.05)
    )
    hits = out.loc[crit, [
        "Gene",
        "n",
        "events",
        "p_interaction",
        "q_interaction",
        "VI",
        "HR_gene_hi_hi",
        "p_hi_hi",
        "HR_gene_lo_lo",
        "p_lo_lo",
        "beta_gene",
        "beta_interaction",
    ]].copy().sort_values(["p_interaction", "VI"], ascending=[True, False])

    hits.to_csv(OUT_HITS, sep="\t", index=False)
    print(f"Saved: {OUT_HITS}")
    print(f"Hits (3 criteria, raw p): {len(hits)}")
    if len(hits) > 0:
        print(hits.head(20).to_string(index=False))


if __name__ == "__main__":
    main()

