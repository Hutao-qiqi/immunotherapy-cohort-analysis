#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

IN_PATH = r"E:\data\changyuan\免疫队列\单基因生存分析\cox_interaction_results.tsv"
OUT_PATH = r"E:\data\changyuan\免疫队列\单基因生存分析\cox_interaction_hiHiRisk_loLoNoRisk.tsv"


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
    df = pd.read_csv(IN_PATH, sep="\t")

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

    hits.to_csv(OUT_PATH, sep="\t", index=False)

    print(f"Input:  {IN_PATH}")
    print(f"Output: {OUT_PATH}")
    print(f"Hits: {hits.shape[0]}")
    print("Top 10:")
    if hits.shape[0] == 0:
        print("(none)")
    else:
        print(hits.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
