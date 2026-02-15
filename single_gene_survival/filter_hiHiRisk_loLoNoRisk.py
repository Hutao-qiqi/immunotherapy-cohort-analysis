#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import argparse
from pathlib import Path

def _default_base_dir() -> Path:
    return Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    base_dir = _default_base_dir()
    parser = argparse.ArgumentParser(description="Filter interaction hits: hi_hi risk + lo_lo no-risk")
    parser.add_argument(
        "--in-path",
        type=Path,
        default=base_dir / "cox_interaction_results.tsv",
        help="Input TSV (default: ./cox_interaction_results.tsv)",
    )
    parser.add_argument(
        "--out-path",
        type=Path,
        default=(base_dir / "outputs" / "cox_interaction_hiHiRisk_loLoNoRisk.tsv"),
        help="Output TSV (default: ./outputs/cox_interaction_hiHiRisk_loLoNoRisk.tsv)",
    )
    return parser.parse_args()


def parse_ci(value):
    if pd.isna(value):
        return (np.nan, np.nan)
    s = str(value).strip()
    if "-" not in s:
        return (np.nan, np.nan)
    left, right = s.split("-", 1)
    try:
        return (float(left), float(right))
    except Exception:
        return (np.nan, np.nan)


def main():
    args = parse_args()
    in_path = args.in_path.resolve()
    out_path = args.out_path.resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(in_path, sep="\t")

    lo_ci = df["CI95_lo_lo"].apply(parse_ci)
    hi_ci = df["CI95_hi_hi"].apply(parse_ci)

    df["lo_ci_l"] = [x[0] for x in lo_ci]
    df["lo_ci_u"] = [x[1] for x in lo_ci]
    df["hi_ci_l"] = [x[0] for x in hi_ci]
    df["hi_ci_u"] = [x[1] for x in hi_ci]

    # Definition aligning with your goal + reviewer requirement:
    # - interaction significant (q<0.05)
    # - hi_hi has risk: HR_hi_hi>1 and 95%CI lower>1
    # - lo_lo no risk: 95%CI includes 1
    # - optional guard: HR_lo_lo <= 1.1 (avoid mild risk in lo_lo)
    crit = (
        (df["q_interaction"] < 0.05)
        & (df["HR_gene_hi_hi"] > 1.0)
        & (df["hi_ci_l"] > 1.0)
        & (df["lo_ci_l"] <= 1.0)
        & (df["lo_ci_u"] >= 1.0)
        & (df["HR_gene_lo_lo"] <= 1.1)
    )

    hits = (
        df.loc[
            crit,
            [
                "Gene",
                "n",
                "events",
                "p_interaction",
                "q_interaction",
                "HR_gene_lo_lo",
                "CI95_lo_lo",
                "HR_gene_hi_hi",
                "CI95_hi_hi",
                "beta_gene",
                "beta_interaction",
            ],
        ]
        .copy()
        .sort_values(["q_interaction", "p_interaction"], ascending=[True, True])
    )

    hits.to_csv(out_path, sep="\t", index=False)

    print(f"Input:  {in_path}")
    print(f"Output: {out_path}")
    print(f"Hits: {hits.shape[0]}")
    print("Top 10:")
    if hits.shape[0] == 0:
        print("(none)")
    else:
        print(hits.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
