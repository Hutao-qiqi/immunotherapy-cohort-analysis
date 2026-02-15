#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd

IN_PATH = r"E:\data\changyuan\免疫队列\单基因生存分析\cox_interaction_results.tsv"

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


def main():
    df = pd.read_csv(IN_PATH, sep="\t")
    sub = df[df["Gene"].isin(GENES)].copy()

    # Keep order same as figure
    sub["_order"] = sub["Gene"].map({g: i for i, g in enumerate(GENES)})
    sub = sub.sort_values("_order").drop(columns=["_order"])

    cols = [
        "Gene",
        "p_interaction",
        "q_interaction",
        "HR_gene_lo_lo",
        "CI95_lo_lo",
        "HR_gene_hi_hi",
        "CI95_hi_hi",
        "beta_interaction",
    ]

    # Some genes might be missing if not fit
    missing = [g for g in GENES if g not in set(sub["Gene"])]

    print("Interaction-term Cox results for HRR figure genes")
    print("Input:", IN_PATH)
    if missing:
        print("Missing (no valid model fit):", ", ".join(missing))
    print(sub[cols].to_string(index=False))


if __name__ == "__main__":
    main()
