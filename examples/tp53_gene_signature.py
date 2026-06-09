"""
TP53 Gene Signature Analysis

Demonstrates Orchestra's gene signature path using a set of canonical TP53
transcriptional targets in epithelial_cell.

Given a differentially expressed gene list, Orchestra:
  1. RegNetAgents find_master_regulators — Fisher enrichment ranks TFs by
     overlap with the input gene set in the ARACNe regulon
  2. Parallel CASCADE comprehensive_perturbation_analysis on top 3 TFs —
     LINCS knockdown, DepMap essentiality, DoRothEA TF confidence, etc.
  3. Synthesis — regulators supported by both enrichment statistics and
     experimental evidence are reported as high-confidence drivers

Expected: TP53 ranks highly in Fisher enrichment (canonical TP53 targets) and
is corroborated by CASCADE (DoRothEA + LINCS + DepMap evidence).

Requirements: RegNetAgents and CASCADE must be installed and accessible.
See docs/installation.md for setup instructions.
"""

import asyncio
import sys
import io
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

from dotenv import load_dotenv
load_dotenv(Path(__file__).parent.parent / ".env")

from orchestra_langgraph_workflow import OrchestraWorkflow

# Canonical TP53 transcriptional target genes (epithelial context)
TP53_SIGNATURE = [
    "CDKN1A", "MDM2",  "BAX",       "BBC3",  "GADD45A",
    "FAS",    "DDB2",  "RRM2B",     "APAF1", "SESN2",
    "PLK3",   "GDF15", "TNFRSF10B", "PERP",  "PMAIP1",
    "ZMAT3",  "FDXR",  "TP53I3",    "CCNG1", "SFN",
]


async def main() -> None:
    print("Orchestra: TP53 Gene Signature Analysis")
    print(f"Input: {len(TP53_SIGNATURE)}-gene signature | Cell type: epithelial_cell | Path: signature")
    print(f"Genes: {', '.join(TP53_SIGNATURE)}")
    print("-" * 60)

    workflow = OrchestraWorkflow()
    result = await workflow.run_analysis(
        gene="",
        cell_type="epithelial_cell",
        analysis_type="gene_signature",
        gene_signature=TP53_SIGNATURE,
    )

    print(result.get("final_report", "(no report generated)"))

    synthesis = result.get("synthesis", {})
    drivers = synthesis.get("ranked_drivers", [])
    if drivers:
        print(f"\nTop driver(s) with cross-system support:")
        for d in drivers[:5]:
            fisher_p = d.get("fisher_pvalue", "n/a")
            overlap = d.get("overlap_count", "n/a")
            corroboration = d.get("corroboration_count", 0)
            cascade_sources = ", ".join(d.get("cascade_sources", [])) or "none"
            print(
                f"  {d['symbol']:12s}  Fisher p={fisher_p}  overlap={overlap}  "
                f"CASCADE {corroboration}/7 ({cascade_sources})"
            )

    mr_result = result.get("master_regulators") or {}
    query_summary = mr_result.get("query_summary", {})
    if query_summary:
        print(f"\nQuery: {query_summary.get('gene_set_size')} input genes, "
              f"{query_summary.get('genes_found_in_network')} found in network")

    errors = result.get("errors", {})
    if errors:
        print(f"\nPartial data — errors: {errors}")
    else:
        print("\nErrors: {} (both systems returned data)")


if __name__ == "__main__":
    asyncio.run(main())
