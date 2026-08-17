"""
Independent replication check for the corroboration-threshold result, take 2 --
STAD instead of HNSC. Chosen specifically because CASCADE's own paper already
validated its downstream perturbation-prediction reliability in STAD (MYC: 85.7%
concordance; CCNE1: 98%; TOP2A: 100%; CCND3: 79.6%), so a non-replication here
cannot be explained away by "CASCADE just doesn't work well in this tissue" the
way the HNSC attempt could -- STAD isolates the corroboration-threshold question
more cleanly.

Focal panel is the 11 genes CONFIRMED to have GREmLN coverage in a pre-screen
(stad_coverage_check.py) -- CASCADE-validated-in-STAD genes (MYC, AURKA, TOP2A,
CCND3) plus TCGA STAD landscape drivers (TP53, PIK3CA, KRAS, SMAD4, ERBB2, APC,
CTNNB1). Four candidates (CCNE1, ARID1A, CDH1, RHOA) were excluded up front for
lacking GREmLN coverage -- not reported as in-experiment failures.
"""
import asyncio
import sys
import io
import json
from pathlib import Path

if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).parent.parent))
from dotenv import load_dotenv
load_dotenv(Path(__file__).parent.parent / ".env")

sys.path.insert(0, str(Path(__file__).parent))
from gene_match_guard import assert_gene_match, GeneMismatchError

from orchestra_langgraph_workflow import OrchestraWorkflow
from mcp_client import make_cascade_client, make_regnetagents_client

CELL_TYPE = "epithelial_cell"

# Confirmed GREmLN-covered subset only (see stad_coverage_check.py / stad_coverage_results.json)
STAD_FOCAL = ["MYC", "AURKA", "TOP2A", "CCND3", "TP53", "PIK3CA", "KRAS", "SMAD4", "ERBB2", "APC", "CTNNB1"]
HOUSEKEEPING = ["ACTB", "GAPDH", "HPRT1", "LDHA", "TUBB"]
NON_DRIVER = ["FASN", "PCNA", "PKM", "PABPC1", "VIM"]

QUERIES = (
    [{"gene": g, "cancer_type": "stad", "category": "focal"} for g in STAD_FOCAL]
    + [{"gene": g, "cancer_type": "stad", "category": "housekeeping"} for g in HOUSEKEEPING]
    + [{"gene": g, "cancer_type": "stad", "category": "non_driver"} for g in NON_DRIVER]
)

OUT_PATH = Path(__file__).parent.parent / "outputs" / "corroboration_replication_stad_raw.json"


async def run_one(workflow: OrchestraWorkflow, gene: str, cancer_type: str, category: str) -> dict:
    print(f"\n=== {gene} / {cancer_type} ({category}) ===", flush=True)
    try:
        result = await asyncio.wait_for(
            workflow.run_analysis(
                gene=gene, cell_type=CELL_TYPE, analysis_type="network_comparison",
                cancer_type=cancer_type, validate_tumor_acquired=True,
            ),
            timeout=180,
        )
    except Exception as e:
        print(f"  ERROR (workflow call): {e!r}")
        return {"gene": gene, "cancer_type": cancer_type, "category": category, "error": repr(e)}

    nc = result.get("network_comparison")
    if not nc:
        err = (result.get("errors") or {}).get("network_comparison", "no network_comparison in result")
        print(f"  ERROR (no network_comparison): {err}")
        return {"gene": gene, "cancer_type": cancer_type, "category": category, "error": err}

    try:
        assert_gene_match(nc, gene)
    except GeneMismatchError as e:
        print(f"  GENE MISMATCH -- discarding: {e}")
        return {"gene": gene, "cancer_type": cancer_type, "category": category, "error": f"gene_mismatch: {e}"}

    tumor_acquired = nc.get("regulators", {}).get("tumor_state_only", [])[:10]
    validation = nc.get("tumor_acquired_cascade_validation", {})

    scored = []
    for reg in tumor_acquired:
        cascade_res = validation.get(reg) or {}
        if "corroboration_count" not in cascade_res:
            scored.append({"gene": reg, "cascade_error": cascade_res.get("error", "no result")})
            continue
        scored.append({"gene": reg, **cascade_res})
        print(f"  {reg}: corroboration_count={cascade_res['corroboration_count']}/4")

    return {
        "gene": gene,
        "cancer_type": cancer_type,
        "category": category,
        "n_tumor_acquired_total": len(nc.get("regulators", {}).get("tumor_state_only", [])),
        "n_scored": len(tumor_acquired),
        "scored_candidates": scored,
        "rewiring": nc.get("interpretation", {}).get("regulatory_rewiring"),
    }


async def main():
    workflow = OrchestraWorkflow()
    OUT_PATH.parent.mkdir(exist_ok=True)

    async with make_cascade_client() as cascade, make_regnetagents_client() as regnetagents:
        workflow._persistent_cascade = cascade
        workflow._persistent_regnetagents = regnetagents

        all_results = []
        for q in QUERIES:
            r = await run_one(workflow, q["gene"], q["cancer_type"], q["category"])
            all_results.append(r)
            with open(OUT_PATH, "w") as f:
                json.dump(all_results, f, indent=2)

    print(f"\n\nDONE. {len(all_results)} queries. Results written to {OUT_PATH}")


if __name__ == "__main__":
    asyncio.run(main())
