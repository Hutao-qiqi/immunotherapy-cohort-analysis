#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import argparse
from pathlib import Path

def _default_in_path() -> Path:
    return (Path(__file__).resolve().parent / "outputs" / "cox_interaction_results.tsv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Query interaction-term Cox results for a fixed gene list")
    parser.add_argument(
        "--in-path",
        type=Path,
        default=_default_in_path(),
        help="Input TSV (default: ./outputs/cox_interaction_results.tsv)",
    )
    return parser.parse_args()

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
    args = parse_args()
    in_path = args.in_path.resolve()
    df = pd.read_csv(in_path, sep="\t")
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
    print("Input:", in_path)
    if missing:
        print("Missing (no valid model fit):", ", ".join(missing))
    print(sub[cols].to_string(index=False))


if __name__ == "__main__":
    main()
