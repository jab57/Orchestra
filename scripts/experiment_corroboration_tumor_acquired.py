"""
Orchestra corroboration-signal experiment: does CASCADE corroboration on the TCGA
tumor-acquired regulator tier predict OncoKB cancer-gene membership better than
uncorroborated candidates?

Calls OrchestraWorkflow.run_analysis() directly (in-process, real Orchestra code --
not bypassing it), with validate_tumor_acquired=True. Mirrors the pattern CASCADE's
own arXiv paper used (invoking CascadeWorkflow.run() directly as the real agentic
entry point, not round-tripping through the MCP stdio protocol).

Panel: reuses RegNetAgents' own published focal-gene panel exactly (BRCA + COAD),
plus its housekeeping/non-driver negative-control panels, for direct comparability
with RegNetAgents' own OncoKB validation methodology.

Every result is checked with a gene-match guard before being trusted -- see
gene_match_guard.py in the same directory.
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

# RegNetAgents' own published focal-gene panel (RegNetAgents_research_paper.tex,
# Section 2.3). Two overlapping genes tested in both cancer types intentionally.
BRCA_FOCAL = ["TP53", "MYC", "CTNNB1", "CCND1", "BRCA2", "PIK3CA", "PTEN", "RB1", "ERBB2", "ESR1", "GATA3"]
COAD_FOCAL = ["TP53", "MYC", "CTNNB1", "CCND1", "KRAS", "APC", "SMAD4", "BRAF", "PIK3CA", "PTEN", "FBXW7", "TCF7L2"]

# RegNetAgents' own negative-control panels (housekeeping genes, tumor-expressed
# non-driver genes) -- tested against both cancer types for full comparability.
HOUSEKEEPING = ["ACTB", "GAPDH", "HPRT1", "LDHA", "TUBB"]
NON_DRIVER = ["FASN", "PCNA", "PKM", "PABPC1", "VIM"]

QUERIES = (
    [{"gene": g, "cancer_type": "brca", "category": "focal"} for g in BRCA_FOCAL]
    + [{"gene": g, "cancer_type": "coad", "category": "focal"} for g in COAD_FOCAL]
    + [
        {"gene": g, "cancer_type": ct, "category": "housekeeping"}
        for g in HOUSEKEEPING
        for ct in ("brca", "coad")
    ]
    + [
        {"gene": g, "cancer_type": ct, "category": "non_driver"}
        for g in NON_DRIVER
        for ct in ("brca", "coad")
    ]
)

OUT_PATH = Path(__file__).parent.parent / "outputs" / "corroboration_tumor_acquired_raw.json"


async def run_one(workflow: OrchestraWorkflow, gene: str, cancer_type: str, category: str) -> dict:
    print(f"\n=== {gene} / {cancer_type} ({category}) ===", flush=True)
    try:
        result = await asyncio.wait_for(
            workflow.run_analysis(
                gene=gene,
                cell_type=CELL_TYPE,
                analysis_type="network_comparison",
                cancer_type=cancer_type,
                validate_tumor_acquired=True,
            ),
            timeout=180,  # lightweight CASCADE tools (LINCS, super-enhancer, DoRothEA, DepMap) — fast
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

    # corroboration_count is now precomputed inline by Orchestra's lightweight-tool
    # validation (LINCS, super-enhancer, DoRothEA, DepMap) -- no separate scoring step
    # needed here, just read it off each candidate's result.
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

    # Open both child-server connections ONCE and reuse them for all 43 queries via
    # run_analysis()'s existing persistent-connection code path (the same path the real
    # MCP server uses once warmed up), instead of the per-call path, which would otherwise
    # respawn + re-warm both subprocesses (model load, DoRothEA/DepMap/LINCS data) on every
    # single query -- ~20-30s of pure overhead x 43, on top of the actual analysis time.
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
