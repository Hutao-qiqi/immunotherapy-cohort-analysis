#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Cox interaction-term reanalysis

Implements the reviewer-recommended model to identify context-dependent
prognostic factors via an interaction term:

    Surv(OS, event) ~ Gene + 8q24_status + Gene:8q24_status

Here, 8q24_status is proxied by MYC/PVT1 dual-expression groups:
  - lo_lo -> 0
  - hi_hi -> 1

Outputs a TSV with interaction P-values and BH-FDR (q-values).
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from lifelines import CoxPHFitter
from tqdm import tqdm


def _default_base_dir() -> Path:
    return Path(__file__).resolve().parent


def _default_survival_file(base_dir: Path) -> Path:
    return (base_dir / ".." / "survival_curves" / "updated_survival_data.txt").resolve()


def parse_args() -> argparse.Namespace:
    base_dir = _default_base_dir()
    parser = argparse.ArgumentParser(description="Cox interaction-term reanalysis")
    parser.add_argument("--base-dir", type=Path, default=base_dir, help="Directory containing expr + annotation")
    parser.add_argument(
        "--expr-file",
        type=Path,
        default=None,
        help="Expression matrix file (default: <base-dir>/combined_expression_combat_corrected.txt)",
    )
    parser.add_argument(
        "--annotation-file",
        type=Path,
        default=None,
        help="Annotation file (default: <base-dir>/MYC_PVT1_annotation.txt)",
    )
    parser.add_argument(
        "--survival-file",
        type=Path,
        default=None,
        help="Survival file (default: ../survival_curves/updated_survival_data.txt)",
    )
    parser.add_argument(
        "--out-tsv",
        type=Path,
        default=base_dir / "outputs" / "cox_interaction_results.tsv",
        help="Output TSV (default: ./outputs/cox_interaction_results.tsv)",
    )
    parser.add_argument("--n-jobs", type=int, default=min(8, os.cpu_count() or 1), help="Parallel jobs")
    return parser.parse_args()


BASE_DIR = _default_base_dir()
EXPR_FILE = (BASE_DIR / "combined_expression_combat_corrected.txt")
ANNOT_FILE = (BASE_DIR / "MYC_PVT1_annotation.txt")
SURV_FILE = _default_survival_file(BASE_DIR)
OUT_TSV = (BASE_DIR / "outputs" / "cox_interaction_results.tsv")

N_JOBS = min(8, os.cpu_count() or 1)
MIN_SAMPLES = 60
MIN_EVENTS = 10


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR for a 1D array; preserves NaNs."""
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
    base_df: pd.DataFrame  # indexed by samples in the same order
    gene_names: list[str]
    expr_values: np.ndarray  # shape: (genes, samples)
    has_dataset: bool


def load_inputs() -> InputData:
    if not os.path.exists(EXPR_FILE):
        raise FileNotFoundError(str(EXPR_FILE))
    if not os.path.exists(ANNOT_FILE):
        raise FileNotFoundError(str(ANNOT_FILE))
    if not os.path.exists(SURV_FILE):
        raise FileNotFoundError(str(SURV_FILE))

    expr = pd.read_csv(EXPR_FILE, sep="\t", index_col=0)

    annot = pd.read_csv(ANNOT_FILE, sep="\t")
    if "Sample" not in annot.columns or "MYC_PVT1_Status" not in annot.columns:
        raise ValueError("MYC_PVT1_annotation.txt must contain columns: Sample, MYC_PVT1_Status")

    status_series = pd.Series(annot["MYC_PVT1_Status"].values, index=annot["Sample"].astype(str).values)
    status_series = status_series[status_series.isin(["hi_hi", "lo_lo"])]

    surv = pd.read_csv(SURV_FILE, sep=r"\s+", engine="python")
    if "Sample_ID" not in surv.columns:
        raise ValueError("updated_survival_data.txt must contain column: Sample_ID")
    surv = surv.set_index("Sample_ID")

    for col in ("OS_months", "OS_event"):
        if col not in surv.columns:
            raise ValueError(f"updated_survival_data.txt missing required column: {col}")

    common = sorted(set(expr.columns) & set(status_series.index) & set(surv.index))
    if len(common) < MIN_SAMPLES:
        raise RuntimeError(f"Too few matched samples: {len(common)}")

    # base survival + status
    base_df = pd.DataFrame(index=common)
    base_df["OS_time"] = pd.to_numeric(surv.loc[common, "OS_months"], errors="coerce")
    base_df["OS_event"] = pd.to_numeric(surv.loc[common, "OS_event"], errors="coerce")
    base_df["status"] = (status_series.loc[common] == "hi_hi").astype(int).values

    has_dataset = "Dataset" in surv.columns
    if has_dataset:
        base_df["Dataset"] = surv.loc[common, "Dataset"].astype(str).values

    # basic QC
    base_df = base_df.dropna(subset=["OS_time", "OS_event", "status"]).copy()
    base_df = base_df[(base_df["OS_time"] > 0) & (base_df["OS_event"].isin([0, 1]))].copy()

    if len(base_df) < MIN_SAMPLES:
        raise RuntimeError(f"Too few usable samples after QC: {len(base_df)}")
    if base_df["OS_event"].sum() < MIN_EVENTS:
        raise RuntimeError(f"Too few events after QC: {int(base_df['OS_event'].sum())}")

    # align expression to base_df sample order
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


def fit_one_gene(i: int, data: InputData) -> Optional[dict]:
    gene = data.gene_names[i]
    x = data.expr_values[i, :]

    if not np.all(np.isfinite(x)):
        # allow partial missingness by dropping rows
        pass

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

    # HR per 1-SD gene increase, within each status stratum
    hr_lo = float(np.exp(beta_gene))
    hr_hi = float(np.exp(beta_gene + beta_int))

    # approximate 95% CI using variance-covariance
    v = cph.variance_matrix_.loc[["gene_z", "interaction"], ["gene_z", "interaction"]].to_numpy(dtype=float)
    var_gene = float(v[0, 0])
    var_int = float(v[1, 1])
    cov = float(v[0, 1])

    se_lo = np.sqrt(max(var_gene, 0.0))
    se_hi = np.sqrt(max(var_gene + var_int + 2.0 * cov, 0.0))

    ci_lo_l = float(np.exp(beta_gene - 1.96 * se_lo))
    ci_lo_u = float(np.exp(beta_gene + 1.96 * se_lo))

    ci_hi_l = float(np.exp((beta_gene + beta_int) - 1.96 * se_hi))
    ci_hi_u = float(np.exp((beta_gene + beta_int) + 1.96 * se_hi))

    return {
        "Gene": gene,
        "n": int(len(fit_df)),
        "events": int(fit_df["OS_event"].sum()),
        "beta_gene": beta_gene,
        "beta_interaction": beta_int,
        "p_interaction": p_int,
        "HR_gene_lo_lo": hr_lo,
        "HR_gene_hi_hi": hr_hi,
        "CI95_lo_lo": f"{ci_lo_l:.3g}-{ci_lo_u:.3g}",
        "CI95_hi_hi": f"{ci_hi_l:.3g}-{ci_hi_u:.3g}",
    }


def main() -> None:
    args = parse_args()
    global BASE_DIR, EXPR_FILE, ANNOT_FILE, SURV_FILE, OUT_TSV, N_JOBS

    BASE_DIR = args.base_dir.resolve()
    EXPR_FILE = (args.expr_file.resolve() if args.expr_file else (BASE_DIR / "combined_expression_combat_corrected.txt"))
    ANNOT_FILE = (args.annotation_file.resolve() if args.annotation_file else (BASE_DIR / "MYC_PVT1_annotation.txt"))
    SURV_FILE = (args.survival_file.resolve() if args.survival_file else _default_survival_file(BASE_DIR))
    OUT_TSV = args.out_tsv.resolve()
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    N_JOBS = int(args.n_jobs)

    print("Loading inputs...")
    data = load_inputs()
    print(f"Matched samples: {len(data.samples)}")
    print(f"Genes: {len(data.gene_names)}")
    print(f"Using dataset stratification: {data.has_dataset}")
    print(f"Parallel jobs: {N_JOBS}")

    results = Parallel(n_jobs=N_JOBS)(
        delayed(fit_one_gene)(i, data)
        for i in tqdm(range(len(data.gene_names)), desc="Cox interaction")
    )

    rows = [r for r in results if r is not None]
    if not rows:
        raise RuntimeError("No valid model fits produced.")

    out = pd.DataFrame(rows)
    out["q_interaction"] = bh_fdr(out["p_interaction"].to_numpy())

    out = out.sort_values(["q_interaction", "p_interaction"], ascending=[True, True])

    out.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"Saved: {OUT_TSV}")

    sig = out[out["q_interaction"] < 0.05]
    print(f"Significant interactions (q<0.05): {len(sig)}")
    print(sig.head(10)[["Gene", "p_interaction", "q_interaction", "HR_gene_lo_lo", "HR_gene_hi_hi"]])


if __name__ == "__main__":
    main()
