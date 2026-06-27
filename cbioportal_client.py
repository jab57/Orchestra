"""
cBioPortal REST API client for TCGA expression + methylation data.

Used by the fetch_tcga_methylation_correlation Orchestra tool to compute
Spearman correlations between a regulator's RNA-seq expression and target
gene promoter methylation beta values across TCGA tumour samples.

No authentication required. SSL bypass via ORCHESTRA_SSL_NO_VERIFY=1.
"""

import asyncio
import math
import os

import requests
from scipy.stats import spearmanr

_BASE = "https://www.cbioportal.org/api"
_TIMEOUT = 30
_SSL_VERIFY = os.getenv("ORCHESTRA_SSL_NO_VERIFY") != "1"

_TCGA_STUDY_IDS: dict[str, str] = {
    code: f"tcga_{code}" for code in [
        "blca", "brca", "cesc", "coad", "hnsc", "kirc",
        "lihc", "luad", "lusc", "ov", "paad", "prad", "stad", "ucec",
    ]
}


def _get(path: str, **params) -> list | dict:
    r = requests.get(f"{_BASE}{path}", params=params, timeout=_TIMEOUT, verify=_SSL_VERIFY)
    r.raise_for_status()
    return r.json()


def _post(path: str, body: dict | list, **params) -> list:
    r = requests.post(f"{_BASE}{path}", json=body, params=params, timeout=_TIMEOUT, verify=_SSL_VERIFY)
    r.raise_for_status()
    return r.json()


def _entrez_id(gene: str) -> int | None:
    """Resolve a Hugo gene symbol to its Entrez gene ID via cBioPortal."""
    try:
        data = _get(f"/genes/{gene}", geneIdType="HUGO_GENE_SYMBOL")
        if isinstance(data, dict):
            return data.get("entrezGeneId")
        return None
    except Exception:
        return None


def _entrez_ids(genes: list[str]) -> dict[str, int]:
    """
    Batch Hugo symbol → Entrez ID via single POST to /genes/fetch.
    Falls back to serial GET calls if the batch endpoint fails.
    """
    if not genes:
        return {}
    try:
        data = _post("/genes/fetch", body=genes, geneIdType="HUGO_GENE_SYMBOL")
        if isinstance(data, list):
            return {
                item["hugoGeneSymbol"]: item["entrezGeneId"]
                for item in data
                if isinstance(item, dict)
                and "hugoGeneSymbol" in item
                and "entrezGeneId" in item
            }
    except Exception:
        pass
    # Serial fallback (batch endpoint unavailable or returned unexpected shape)
    result = {}
    for gene in genes:
        eid = _entrez_id(gene)
        if eid is not None:
            result[gene] = eid
    return result


def _discover_profiles(study_id: str) -> tuple[str | None, str | None]:
    """
    Return (rna_profile_id, methylation_profile_id) for a TCGA study.
    Prefers rna_seq_v2_mrna and methylation_hm450; falls back to any
    MRNA_EXPRESSION / METHYLATION profile if the preferred ID is absent.
    """
    try:
        profiles = _get("/molecular-profiles", studyId=study_id)
    except Exception:
        return None, None

    by_id = {p["molecularProfileId"]: p for p in profiles}

    # RNA-seq: prefer rna_seq_v2_mrna, then any MRNA_EXPRESSION profile
    rna_id: str | None = None
    preferred_rna = f"{study_id}_rna_seq_v2_mrna"
    if preferred_rna in by_id:
        rna_id = preferred_rna
    else:
        for p in profiles:
            if p.get("molecularAlterationType") == "MRNA_EXPRESSION":
                rna_id = p["molecularProfileId"]
                break

    # Methylation: prefer HM450 over HM27, then any METHYLATION profile
    meth_id: str | None = None
    for suffix in ["methylation_hm450", "methylation_hm27"]:
        candidate = f"{study_id}_{suffix}"
        if candidate in by_id:
            meth_id = candidate
            break
    if meth_id is None:
        for p in profiles:
            if p.get("molecularAlterationType") == "METHYLATION":
                meth_id = p["molecularProfileId"]
                break

    return rna_id, meth_id


def _get_sample_ids(study_id: str) -> list[str]:
    """Return all sample IDs for a study, handling pagination."""
    ids: list[str] = []
    page = 0
    while True:
        batch = _get(f"/studies/{study_id}/samples", pageSize=500, pageNumber=page)
        if not batch:
            break
        ids.extend(s["sampleId"] for s in batch)
        if len(batch) < 500:
            break
        page += 1
    return ids


def _fetch_values(profile_id: str, entrez_id: int, sample_ids: list[str]) -> dict[str, float]:
    """
    Return {sampleId: float_value} for one gene in one molecular profile.
    Processes samples in batches of 500 (API limit). Skips failed batches silently.
    """
    values: dict[str, float] = {}
    for i in range(0, len(sample_ids), 500):
        batch = sample_ids[i:i + 500]
        try:
            rows = _post(
                "/molecular-data/fetch",
                body={"entrezGeneIds": [entrez_id], "sampleIds": batch},
                molecularProfileId=profile_id,
                projection="SUMMARY",
            )
            for row in rows:
                sid = row.get("sampleId")
                val = row.get("value")
                if sid and val is not None:
                    try:
                        values[sid] = float(val)
                    except (ValueError, TypeError):
                        pass
        except Exception:
            pass
    return values


def _correlation_sync(
    regulator: str,
    target_genes: list[str],
    tcga_network: str,
) -> dict:
    study_id = _TCGA_STUDY_IDS.get(tcga_network.lower())
    if not study_id:
        valid = ", ".join(sorted(_TCGA_STUDY_IDS))
        return {"error": f"Unknown TCGA code {tcga_network!r}. Valid: {valid}"}

    rna_profile, meth_profile = _discover_profiles(study_id)
    if not rna_profile:
        return {"error": f"No RNA-seq expression profile found for {study_id}"}
    if not meth_profile:
        return {"error": f"No methylation profile found for {study_id}"}

    all_genes = [regulator] + list(target_genes)
    entrez_map = _entrez_ids(all_genes)

    if regulator not in entrez_map:
        return {"error": f"Could not resolve Entrez ID for regulator {regulator!r}"}

    sample_ids = _get_sample_ids(study_id)
    if not sample_ids:
        return {"error": f"No samples found for {study_id}"}

    expr = _fetch_values(rna_profile, entrez_map[regulator], sample_ids)
    if not expr:
        return {"error": f"No expression data returned for {regulator} in {rna_profile}"}

    correlations = []
    for target in target_genes:
        if target not in entrez_map:
            correlations.append({
                "target_gene": target,
                "n_samples": 0,
                "rho": None,
                "p_value": None,
                "direction": "error",
                "note": f"Could not resolve Entrez ID for {target}",
            })
            continue

        meth = _fetch_values(meth_profile, entrez_map[target], sample_ids)
        shared = sorted(set(expr) & set(meth))

        if len(shared) < 5:
            correlations.append({
                "target_gene": target,
                "n_samples": len(shared),
                "rho": None,
                "p_value": None,
                "direction": "insufficient_data",
                "note": f"Only {len(shared)} samples with matched data (need ≥ 5)",
            })
            continue

        x = [expr[s] for s in shared]
        y = [meth[s] for s in shared]
        # spearmanr returns SpearmanrResult(statistic, pvalue); unpack as tuple
        # scipy stubs are incomplete in Pyright so we access by position.
        _rr = spearmanr(x, y)
        rho_raw: float = float(_rr[0])  # type: ignore[arg-type]
        p_raw: float = float(_rr[1])  # type: ignore[arg-type]

        rho_f = rho_raw if not math.isnan(rho_raw) else None
        p_f = p_raw if not math.isnan(p_raw) else None

        if rho_f is None:
            direction = "undefined"
        elif rho_f < -0.2:
            direction = "inverse (high expression → low methylation)"
        elif rho_f > 0.2:
            direction = "concordant (high expression → high methylation)"
        else:
            direction = "no clear trend (|ρ| < 0.2)"

        correlations.append({
            "target_gene": target,
            "n_samples": len(shared),
            "rho": round(rho_f, 3) if rho_f is not None else None,
            "p_value": p_f,
            "direction": direction,
        })

    return {
        "regulator": regulator,
        "tcga_network": tcga_network,
        "study_id": study_id,
        "rna_profile": rna_profile,
        "methylation_profile": meth_profile,
        "n_total_samples": len(sample_ids),
        "correlations": correlations,
    }


async def methylation_expression_correlation(
    regulator: str,
    target_genes: list[str],
    tcga_network: str,
) -> dict:
    """
    Spearman correlation between a regulator's RNA-seq expression and target gene
    methylation beta values across a TCGA cohort. All HTTP calls run in a thread
    so the asyncio event loop is not blocked.
    """
    return await asyncio.to_thread(_correlation_sync, regulator, target_genes, tcga_network)


def format_correlation_report(result: dict) -> str:
    """Format a methylation_expression_correlation result as a markdown report."""
    if "error" in result:
        return f"## TCGA Methylation-Expression Correlation\n\n**Error:** {result['error']}"

    regulator = result["regulator"]
    network = result["tcga_network"].upper()
    n_samples = result["n_total_samples"]
    rna_p = result["rna_profile"]
    meth_p = result["methylation_profile"]
    correlations = result["correlations"]

    lines = [
        f"## TCGA Methylation-Expression Correlation",
        f"**Regulator:** {regulator} (RNA-seq expression)  "
        f"**Cohort:** {network}  |  n = {n_samples} samples",
        f"**Profiles:** `{rna_p}` · `{meth_p}`",
        "",
        "| Target Gene | n | Spearman ρ | p-value | Direction |",
        "|-------------|---|------------|---------|-----------|",
    ]

    for c in correlations:
        target = c["target_gene"]
        n = c["n_samples"]
        rho = f"{c['rho']:+.3f}" if c["rho"] is not None else "—"
        p = f"{c['p_value']:.2e}" if c["p_value"] is not None else "—"
        direction = c.get("direction", "")
        note = c.get("note", "")
        direction_cell = note if note else direction
        lines.append(f"| {target} | {n} | {rho} | {p} | {direction_cell} |")

    lines += [
        "",
        "_Negative ρ (high regulator expression → low methylation beta) supports the hypothesis_",
        "_that regulator activity antagonises epigenetic silencing of these targets._",
        "_Thresholds: |ρ| > 0.2 = directional trend · p < 0.05 = nominally significant_",
        "",
        "> ⚠️ Computationally derived from bulk tumour RNA-seq + array methylation. "
        "Requires experimental validation.",
    ]
    return "\n".join(lines)
