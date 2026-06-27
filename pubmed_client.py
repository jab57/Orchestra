"""
PubMed novelty assessment via NCBI E-utilities.

Queries PubMed for a gene (or gene pair) in a cancer context and returns
a structured novelty verdict: established / emerging / novel.

Rate limits: 3 req/s without API key, 10 req/s with NCBI_API_KEY in .env.
All HTTP calls are synchronous (requests); wrapped in asyncio.to_thread for
use from async callers.
"""

import asyncio
import os
import xml.etree.ElementTree as ET

import requests
from dotenv import load_dotenv

load_dotenv()

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_REQUEST_TIMEOUT = 15  # seconds
_SSL_VERIFY = os.getenv("ORCHESTRA_SSL_NO_VERIFY") != "1"

# Synonym expansion for each cancer context string produced by _TCGA_TO_CANCER_CONTEXT.
# Each value is a self-contained OR clause ready to embed in a PubMed boolean query.
# Unknown contexts fall back to a plain [tiab] literal in _build_base_query.
_CANCER_SYNONYMS: dict[str, str] = {
    "cervical cancer": (
        '("cervical cancer"[tiab] OR "cervical carcinoma"[tiab] OR '
        '"cervical squamous cell carcinoma"[tiab] OR "CESC"[tiab])'
    ),
    "head and neck squamous": (
        '("head and neck squamous"[tiab] OR "HNSCC"[tiab] OR '
        '"head and neck cancer"[tiab] OR "head and neck squamous cell carcinoma"[tiab])'
    ),
    "colorectal": (
        '("colorectal"[tiab] OR "colorectal cancer"[tiab] OR '
        '"colon cancer"[tiab] OR "rectal cancer"[tiab])'
    ),
    "breast cancer": '("breast cancer"[tiab] OR "breast carcinoma"[tiab])',
    "lung adenocarcinoma": '("lung adenocarcinoma"[tiab] OR "LUAD"[tiab])',
    "lung squamous cell carcinoma": (
        '("lung squamous cell carcinoma"[tiab] OR "LUSC"[tiab] OR '
        '"lung squamous carcinoma"[tiab])'
    ),
    "hepatocellular carcinoma": (
        '("hepatocellular carcinoma"[tiab] OR "liver cancer"[tiab] OR "HCC"[tiab])'
    ),
    "ovarian cancer": '("ovarian cancer"[tiab] OR "ovarian carcinoma"[tiab])',
    "bladder cancer": '("bladder cancer"[tiab] OR "urothelial carcinoma"[tiab])',
    "kidney clear cell carcinoma": (
        '("kidney clear cell carcinoma"[tiab] OR "clear cell renal cell carcinoma"[tiab] OR '
        '"ccRCC"[tiab])'
    ),
    "gastric cancer": '("gastric cancer"[tiab] OR "stomach cancer"[tiab])',
    "pancreatic cancer": (
        '("pancreatic cancer"[tiab] OR "pancreatic ductal adenocarcinoma"[tiab] OR "PDAC"[tiab])'
    ),
    "prostate cancer": '("prostate cancer"[tiab] OR "prostate carcinoma"[tiab])',
    "endometrial cancer": (
        '("endometrial cancer"[tiab] OR "endometrial carcinoma"[tiab] OR "uterine cancer"[tiab])'
    ),
}

# Pair queries are far more specific than single-gene queries, so a lower hit count
# still indicates a well-characterised relationship.
_VERDICT_THRESHOLDS: dict[str, dict[str, int]] = {
    "single": {"established": 20, "emerging": 5},
    "pair":   {"established": 10, "emerging": 3},
}


def _ncbi_params(extra: dict) -> dict:
    params = {"db": "pubmed", "retmode": "xml", **extra}
    api_key = os.getenv("NCBI_API_KEY")
    if api_key:
        params["api_key"] = api_key
    return params


def _build_base_query(gene: str, cancer_context: str, gene2: str | None) -> str:
    ctx = _CANCER_SYNONYMS.get(cancer_context, f'"{cancer_context}"[tiab]')
    if gene2:
        return f'"{gene}"[tiab] AND "{gene2}"[tiab] AND {ctx}'
    return f'"{gene}"[tiab] AND {ctx}'


def _build_experimental_query(base: str) -> str:
    return (
        f'({base}) AND '
        f'(("in vitro"[tiab] OR "in vivo"[tiab] OR "experimental"[tiab]) '
        f'NOT ("computational"[tiab] OR "bioinformatic"[tiab] OR "bioinformatics"[tiab]))'
    )


def _esearch(query: str, retmax: int = 0, sort: str = "relevance") -> tuple[int, list[str]]:
    """Return (total_count, id_list) for a PubMed query."""
    params = _ncbi_params({"term": query, "retmax": str(retmax), "sort": sort})
    resp = requests.get(f"{NCBI_BASE}/esearch.fcgi", params=params, timeout=_REQUEST_TIMEOUT, verify=_SSL_VERIFY)
    resp.raise_for_status()
    root = ET.fromstring(resp.text)
    count_el = root.find("Count")
    count = int(count_el.text) if count_el is not None and count_el.text else 0
    ids = [el.text for el in root.findall(".//Id") if el.text]
    return count, ids


def _efetch_year(pmid: str) -> int | None:
    """Return the publication year for a PMID."""
    params = _ncbi_params({"id": pmid, "rettype": "abstract"})
    resp = requests.get(f"{NCBI_BASE}/efetch.fcgi", params=params, timeout=_REQUEST_TIMEOUT, verify=_SSL_VERIFY)
    resp.raise_for_status()
    root = ET.fromstring(resp.text)
    year_el = root.find(".//PubDate/Year")
    if year_el is not None and year_el.text:
        return int(year_el.text)
    medline_el = root.find(".//PubDate/MedlineDate")
    if medline_el is not None and medline_el.text:
        try:
            return int(medline_el.text[:4])
        except ValueError:
            return None
    return None


def _verdict(pubmed_hits: int, query_type: str = "single") -> str:
    t = _VERDICT_THRESHOLDS.get(query_type, _VERDICT_THRESHOLDS["single"])
    if pubmed_hits > t["established"]:
        return "established"
    if pubmed_hits >= t["emerging"]:
        return "emerging"
    return "novel"


def _rationale(pubmed_hits: int, experimental_hits: int, verdict: str) -> str:
    if pubmed_hits == 0:
        return "No prior papers found in this gene-cancer context"
    exp_desc = f"{experimental_hits} experimental" if experimental_hits > 0 else "no experimental"
    comp_hits = max(0, pubmed_hits - experimental_hits)
    comp_desc = f"{comp_hits} computational" if comp_hits > 0 else "no computational"
    base = f"{pubmed_hits} prior paper{'s' if pubmed_hits != 1 else ''}; {exp_desc}, {comp_desc}"
    if verdict == "novel":
        return f"{base} — limited prior characterization in this context"
    if verdict == "emerging":
        return f"{base} — active area with supporting evidence"
    return f"{base} — well-characterized in this context"


def _novelty_assessment_sync(
    gene: str,
    cancer_context: str,
    gene2: str | None,
    return_pmids: bool = False,
) -> dict:
    query_type = "pair" if gene2 else "single"
    base_query = _build_base_query(gene, cancer_context, gene2)
    exp_query = _build_experimental_query(base_query)

    # Call 1: total count + most recent PMID (single esearch with retmax=1)
    total, recent_ids = _esearch(base_query, retmax=1, sort="pub_date")

    # Call 2: experimental subset count
    experimental = 0
    if total > 0:
        experimental, _ = _esearch(exp_query, retmax=0)

    # Call 3: publication year of most recent paper
    most_recent_year = None
    if recent_ids:
        most_recent_year = _efetch_year(recent_ids[0])

    # Call 4 (opt-in): top-10 PMIDs by relevance for independent spot-checking
    top_pmids: list[str] = []
    if return_pmids and total > 0:
        _, top_pmids = _esearch(base_query, retmax=10, sort="relevance")

    verdict = _verdict(total, query_type)
    return {
        "gene": gene,
        "gene2": gene2,
        "cancer_context": cancer_context,
        "pubmed_hits": total,
        "experimental_hits": experimental,
        "computational_hits": max(0, total - experimental),
        "most_recent_year": most_recent_year,
        "novelty_verdict": verdict,
        "verdict_rationale": _rationale(total, experimental, verdict),
        "top_pmids": top_pmids,
    }


async def novelty_assessment(
    gene: str,
    cancer_context: str,
    gene2: str | None = None,
    return_pmids: bool = False,
) -> dict:
    """Query PubMed and return a structured novelty verdict for a gene (or gene pair) in a cancer context."""
    return await asyncio.to_thread(_novelty_assessment_sync, gene, cancer_context, gene2, return_pmids)


def format_novelty_report(result: dict) -> str:
    """Format a novelty_assessment result as a human-readable text report."""
    gene = result["gene"]
    gene2 = result.get("gene2")
    context = result["cancer_context"]
    subject = f"{gene}/{gene2}" if gene2 else gene

    lines = [
        f"Novelty Assessment: {subject} in {context}",
        "",
        f"  PubMed hits:       {result['pubmed_hits']}",
        f"    Experimental:    {result['experimental_hits']}",
        f"    Computational:   {result['computational_hits']}",
        f"  Most recent:       {result['most_recent_year'] or 'N/A'}",
        f"  Verdict:           {result['novelty_verdict'].upper()}",
        "",
        f"  {result['verdict_rationale']}",
    ]
    return "\n".join(lines)
