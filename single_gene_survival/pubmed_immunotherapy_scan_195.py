#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PubMed immunotherapy literature scan for the 195 continuous-interaction hits.

Inputs:
  - selected_continuous_3criteria_p0.05.tsv  (must contain column: Gene)

Outputs:
  - immunotherapy_search_195_genes.tsv
  - immunotherapy_evidence_195_strict.tsv

Method (lightweight, automated):
  1) For each gene, PubMed ESearch for immunotherapy/ICI-related keywords.
  2) Fetch titles + abstracts for the top PMIDs (retmax=5).
  3) Heuristically flag "experimental evidence" if abstract mentions ICI/BCG etc
     AND contains experiment cues (mouse/in vivo/xenograft/CRISPR/knockdown/...).

Note: This is an automated triage. You should manually verify the strict hits
before using them in a manuscript.
"""

from __future__ import annotations

import csv
import itertools
import os
import re
import time
import urllib.parse
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import requests


IN_TSV = "selected_continuous_3criteria_p0.05.tsv"
OUT_SEARCH = "immunotherapy_search_195_genes.tsv"
OUT_STRICT = "immunotherapy_evidence_195_strict.tsv"
OUT_EXPERIMENTAL = "immunotherapy_evidence_195_experimental.tsv"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Keep this modest to avoid rate-limit issues.
SLEEP_SEC = 0.34  # ~3 requests/sec
RETMAX_PMIDS = 5


IMMUNO_TERMS = [
    "immunotherapy",
    "\"immune checkpoint\"",
    "\"immune checkpoint blockade\"",
    "\"PD-1\"",
    "\"PD-L1\"",
    "\"CTLA-4\"",
    "\"anti-PD-1\"",
    "\"anti-PD-L1\"",
    "\"anti-CTLA-4\"",
    "ICI",
    "ICB",
    "pembrolizumab",
    "nivolumab",
    "atezolizumab",
    "durvalumab",
    "ipilimumab",
    "BCG",
]

CANCER_TERMS = [
    "cancer",
    "tumor",
    "tumour",
    "carcinoma",
    "melanoma",
    "sarcoma",
    "leukemia",
    "lymphoma",
]

# Experimental evidence cues (very rough; used only for triage).
EXPERIMENT_CUES = [
    "in vivo",
    "mouse",
    "mice",
    "xenograft",
    "syngeneic",
    "knockdown",
    "silencing",
    "overexpression",
    "CRISPR",
    "sgRNA",
    "shRNA",
    "siRNA",
    "KO",
    "knockout",
    "blocking antibody",
    "blockade",
]

ICI_CUES = [
    "anti-pd-1",
    "anti-pd-l1",
    "pd-1 blockade",
    "pd-l1 blockade",
    "checkpoint blockade",
    "immune checkpoint",
    "pembrolizumab",
    "nivolumab",
    "atezolizumab",
    "durvalumab",
    "ipilimumab",
    "bcg",
]

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

# Common non-experimental / purely computational phrases in abstracts.
EXCLUDE_CUES = [
    "bioinformatics",
    "in silico",
    "signature",
    "risk model",
    "prognostic model",
    "nomogram",
    "machine learning",
    "tcga",
    "geo dataset",
    "data mining",
    "algorithm",
]

@dataclass(frozen=True)
class PubMedRecord:
    pmid: str
    title: str
    abstract: str


def _read_genes(path: str) -> List[str]:
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        genes = []
        for row in r:
            g = (row.get("Gene") or "").strip()
            if g:
                genes.append(g)
    # unique but preserve order
    seen = set()
    out = []
    for g in genes:
        if g not in seen:
            out.append(g)
            seen.add(g)
    return out


def _chunked(items: List[str], n: int) -> Iterable[List[str]]:
    for i in range(0, len(items), n):
        yield items[i : i + n]


def _esearch_query(gene: str) -> str:
    gene_term = f"\"{gene}\"[Title/Abstract]"
    immuno = "(" + " OR ".join(f"{t}[Title/Abstract]" if not t.startswith("\"") else f"{t}[Title/Abstract]" for t in IMMUNO_TERMS) + ")"
    cancer = "(" + " OR ".join(f"{t}[Title/Abstract]" for t in CANCER_TERMS) + ")"
    # Keep query simple; title/abstract constraints reduce noise.
    return f"({gene_term}) AND {immuno} AND {cancer}"


def _esearch(gene: str) -> Tuple[int, List[str]]:
    term = _esearch_query(gene)
    params = {
        "db": "pubmed",
        "term": term,
        "retmax": str(RETMAX_PMIDS),
        "retmode": "json",
        "sort": "relevance",
        "tool": "codex-cli",
    }
    url = f"{EUTILS}/esearch.fcgi?{urllib.parse.urlencode(params)}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()["esearchresult"]
    count = int(data.get("count", "0") or "0")
    pmids = [str(x) for x in data.get("idlist", [])]
    return count, pmids


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
                label = at.attrib.get("Label")
                txt = "".join(at.itertext()).strip()
                if not txt:
                    continue
                parts.append(f"{label}: {txt}" if label else txt)
            abstract = "\n".join(parts).strip()

        out[pmid] = PubMedRecord(pmid=pmid, title=title, abstract=abstract)

    return out


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", (s or "")).strip().lower()


def _has_any(text: str, cues: List[str]) -> bool:
    t = _norm(text)
    return any(c.lower() in t for c in cues)


def _evidence_type(text: str) -> str:
    t = _norm(text)
    if any(x in t for x in ["patient", "patients", "clinical trial", "treated with", "cohort"]):
        return "clinical"
    if any(x in t for x in ["mouse", "mice", "in vivo", "xenograft", "syngeneic"]):
        return "in vivo"
    if any(x in t for x in ["cell line", "in vitro", "knockdown", "overexpression", "crisper", "crispr", "shrna", "sirna"]):
        return "in vitro"
    return "unspecified"

def _efetch_with_retry(pmids: List[str], max_batch: int = 200) -> Dict[str, PubMedRecord]:
    """
    Robust efetch: if a large batch fails, retry with smaller batches.
    """
    out: Dict[str, PubMedRecord] = {}
    batches = list(_chunked(pmids, max_batch))
    for batch in batches:
        try:
            out.update(_efetch_pmids(batch))
        except Exception:
            # Split once and retry.
            if len(batch) <= 10:
                continue
            mid = len(batch) // 2
            for sub in (batch[:mid], batch[mid:]):
                try:
                    out.update(_efetch_pmids(sub))
                except Exception:
                    continue
        time.sleep(SLEEP_SEC)
    return out


def main() -> None:
    if not os.path.exists(IN_TSV):
        raise FileNotFoundError(IN_TSV)

    genes = _read_genes(IN_TSV)
    print(f"Genes: {len(genes)}")

    gene_hits: List[Dict[str, object]] = []
    gene_to_pmids: Dict[str, List[str]] = {}
    all_pmids: List[str] = []

    for i, gene in enumerate(genes, 1):
        try:
            count, pmids = _esearch(gene)
        except Exception as e:
            print(f"[{i}/{len(genes)}] {gene}: ESearch failed: {e}")
            count, pmids = 0, []
        time.sleep(SLEEP_SEC)

        gene_to_pmids[gene] = pmids
        all_pmids.extend(pmids)
        gene_hits.append(
            {
                "Gene": gene,
                "Count": count,
                "PMIDs": ",".join(pmids),
            }
        )

        if i % 25 == 0:
            print(f"  processed {i}/{len(genes)}")

    # Fetch details for unique PMIDs in batches.
    uniq_pmids = list(dict.fromkeys([p for p in all_pmids if p]))
    print(f"Unique PMIDs to fetch: {len(uniq_pmids)}")

    pmid_to_record = _efetch_with_retry(uniq_pmids, max_batch=200)

    # Build enriched search table.
    out_rows: List[Dict[str, object]] = []
    strict_rows: List[Dict[str, object]] = []
    experimental_rows: List[Dict[str, object]] = []

    for row in gene_hits:
        gene = str(row["Gene"])
        pmids = gene_to_pmids.get(gene, [])
        recs = [pmid_to_record.get(p) for p in pmids if p in pmid_to_record]
        titles = " || ".join([r.title for r in recs if r and r.title][:RETMAX_PMIDS])
        abstracts = " || ".join([_norm(r.abstract) for r in recs if r and r.abstract][:RETMAX_PMIDS])

        has_ici = _has_any(titles + "\n" + abstracts, ICI_CUES)
        has_agent = _has_any(titles + "\n" + abstracts, AGENT_CUES)
        has_exp = _has_any(abstracts, EXPERIMENT_CUES)
        excluded = _has_any(titles + "\n" + abstracts, EXCLUDE_CUES)
        strict = bool(has_ici and has_exp and pmids)
        experimental = bool(has_agent and has_exp and not excluded and pmids)

        out_rows.append(
            {
                "Gene": gene,
                "Count": row["Count"],
                "PMIDs": row["PMIDs"],
                "Top_titles": titles,
                "Flags": ";".join(
                    [
                        "ICI" if has_ici else "",
                        "AGENT" if has_agent else "",
                        "EXP" if has_exp else "",
                        "EXCL" if excluded else "",
                        "STRICT" if strict else "",
                        "EXP_ONLY" if experimental else "",
                    ]
                ).strip(";"),
            }
        )

        if strict:
            # Keep the most "obvious" record as the evidence line.
            best = None
            for r in recs:
                if r and _has_any(r.title + "\n" + r.abstract, ICI_CUES) and _has_any(r.abstract, EXPERIMENT_CUES):
                    best = r
                    break
            if best is None and recs:
                best = recs[0]

            if best is not None:
                strict_rows.append(
                    {
                        "Gene": gene,
                        "Evidence_type": _evidence_type(best.abstract),
                        "Immunotherapy_context": "ICI/BCG-related (auto triage)",
                        "PMID": best.pmid,
                        "Title": best.title,
                    }
                )

        if experimental:
            best = None
            for r in recs:
                if not r:
                    continue
                txt = r.title + "\n" + r.abstract
                if _has_any(txt, AGENT_CUES) and _has_any(r.abstract, EXPERIMENT_CUES) and not _has_any(txt, EXCLUDE_CUES):
                    best = r
                    break
            if best is None and recs:
                best = recs[0]
            if best is not None:
                experimental_rows.append(
                    {
                        "Gene": gene,
                        "Evidence_type": _evidence_type(best.abstract),
                        "Immunotherapy_context": "explicit agent + experimental cues (auto triage)",
                        "PMID": best.pmid,
                        "Title": best.title,
                    }
                )

    # Write outputs.
    with open(OUT_SEARCH, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["Gene", "Count", "PMIDs", "Top_titles", "Flags"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(out_rows)

    with open(OUT_STRICT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["Gene", "Evidence_type", "Immunotherapy_context", "PMID", "Title"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(strict_rows)

    with open(OUT_EXPERIMENTAL, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["Gene", "Evidence_type", "Immunotherapy_context", "PMID", "Title"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(experimental_rows)

    print(f"Wrote: {OUT_SEARCH}")
    print(f"Wrote: {OUT_STRICT}")
    print(f"Wrote: {OUT_EXPERIMENTAL}")
    print(f"Strict hits: {len(strict_rows)}")
    print(f"Experimental hits: {len(experimental_rows)}")


if __name__ == "__main__":
    main()
