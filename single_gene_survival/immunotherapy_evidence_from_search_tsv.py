#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Build an "agent + experimental cues" immunotherapy evidence table from an
existing PubMed search TSV (immunotherapy_search_195_genes.tsv).

This avoids re-running per-gene ESearch; it only EFetch'es the PMIDs listed
in the search TSV and triages evidence based on titles/abstracts.
"""

from __future__ import annotations

import csv
import os
import re
import time
import urllib.parse
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional

import requests


IN_SEARCH = "immunotherapy_search_195_genes.tsv"
OUT_AGENT_EXP = "immunotherapy_evidence_195_agent_exp.tsv"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SLEEP_SEC = 0.34

AGENT_CUES = [
    "anti-pd-1",
    "anti-pd-l1",
    "anti-ctla-4",
    "pembrolizumab",
    "nivolumab",
    "atezolizumab",
    "durvalumab",
    "ipilimumab",
    "car-t",
    "car t",
    "bcg",
]

EXPERIMENT_CUES = [
    "in vivo",
    "mouse",
    "mice",
    "xenograft",
    "syngeneic",
    "knockdown",
    "silencing",
    "overexpression",
    "crispr",
    "shrna",
    "sirna",
    "knockout",
    "blocking antibody",
    "blockade",
]


@dataclass(frozen=True)
class PubMedRecord:
    pmid: str
    title: str
    abstract: str


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", (s or "")).strip().lower()


def _has_any(text: str, cues: List[str]) -> bool:
    t = _norm(text)
    return any(c.lower() in t for c in cues)


def _chunked(items: List[str], n: int) -> Iterable[List[str]]:
    for i in range(0, len(items), n):
        yield items[i : i + n]


def _efetch_pmids(pmids: List[str]) -> Dict[str, PubMedRecord]:
    if not pmids:
        return {}
    params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
        "tool": "codex-cli",
    }
    url = f"{EUTILS}/efetch.fcgi?{urllib.parse.urlencode(params)}"
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()

    root = ET.fromstring(resp.text)
    out: Dict[str, PubMedRecord] = {}
    for article in root.findall(".//PubmedArticle"):
        pmid_el = article.find(".//PMID")
        if pmid_el is None or not pmid_el.text:
            continue
        pmid = pmid_el.text.strip()

        title_el = article.find(".//ArticleTitle")
        title = "".join(title_el.itertext()).strip() if title_el is not None else ""

        abs_el = article.find(".//Abstract")
        if abs_el is None:
            abstract = ""
        else:
            parts = []
            for at in abs_el.findall(".//AbstractText"):
                txt = "".join(at.itertext()).strip()
                if txt:
                    parts.append(txt)
            abstract = "\n".join(parts).strip()

        out[pmid] = PubMedRecord(pmid=pmid, title=title, abstract=abstract)
    return out


def main() -> None:
    if not os.path.exists(IN_SEARCH):
        raise FileNotFoundError(IN_SEARCH)

    gene_to_pmids: Dict[str, List[str]] = {}
    all_pmids: List[str] = []

    with open(IN_SEARCH, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            gene = (row.get("Gene") or "").strip()
            pmids = [(x or "").strip() for x in (row.get("PMIDs") or "").split(",") if (x or "").strip()]
            if not gene:
                continue
            gene_to_pmids[gene] = pmids
            all_pmids.extend(pmids)

    uniq_pmids = list(dict.fromkeys(all_pmids))
    pmid_to_rec: Dict[str, PubMedRecord] = {}
    for batch in _chunked(uniq_pmids, 200):
        try:
            pmid_to_rec.update(_efetch_pmids(batch))
        except Exception:
            # retry smaller
            for sub in _chunked(batch, 50):
                try:
                    pmid_to_rec.update(_efetch_pmids(sub))
                except Exception:
                    continue
                time.sleep(SLEEP_SEC)
        time.sleep(SLEEP_SEC)

    out_rows: List[Dict[str, str]] = []
    for gene, pmids in gene_to_pmids.items():
        best: Optional[PubMedRecord] = None
        for p in pmids:
            rec = pmid_to_rec.get(p)
            if not rec:
                continue
            txt = rec.title + "\n" + rec.abstract
            if _has_any(txt, AGENT_CUES) and _has_any(rec.abstract, EXPERIMENT_CUES):
                best = rec
                break
        if best is None:
            continue
        out_rows.append(
            {
                "Gene": gene,
                "PMID": best.pmid,
                "Title": best.title,
            }
        )

    out_rows.sort(key=lambda d: (d["Gene"]))
    with open(OUT_AGENT_EXP, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["Gene", "PMID", "Title"], delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    print(f"Wrote: {OUT_AGENT_EXP}")
    print(f"Agent+EXP hits: {len(out_rows)}")


if __name__ == "__main__":
    main()

