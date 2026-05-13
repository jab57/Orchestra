"""
TP53 Causal Chain Analysis

Demonstrates Orchestra's TF path with TP53 — a canonical master regulator.

TP53 is classified as a master regulator by CASCADE. Orchestra runs:
  1. CASCADE get_gene_metadata → gene_role: master_regulator
  2. Parallel:
     - RegNetAgents comprehensive_gene_analysis → network rank, pathway enrichment
     - CASCADE comprehensive_perturbation_analysis → LINCS, DepMap, DoRothEA, STRING
  3. Cross-system synthesis → genes corroborated by both network topology and
     experimental evidence

Expected: CDKN1A with 3 CASCADE sources (DoRothEA + LINCS + STRING) — a canonical
TP53 target validated in both systems, confirming the cross-system scoring works.

Requirements: RegNetAgents and CASCADE must be installed and accessible.
See docs/installation.md for setup instructions.
"""

import asyncio
import sys
import io
from pathlib import Path

if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

from dotenv import load_dotenv
load_dotenv(Path(__file__).parent.parent / ".env")

from orchestra_langgraph_workflow import OrchestraWorkflow


async def main() -> None:
    print("Orchestra: TP53 Causal Chain Analysis")
    print("Gene: TP53 | Cell type: epithelial_cell | Path: TF")
    print("-" * 60)

    workflow = OrchestraWorkflow()
    result = await workflow.run_analysis(
        gene="TP53",
        cell_type="epithelial_cell",
        analysis_type="causal_chain",
    )

    print(result.get("final_report", "(no report generated)"))

    synthesis = result.get("synthesis", {})
    cross_hits = synthesis.get("cross_system_hits", [])
    if cross_hits:
        print(f"\nCross-system hits ({len(cross_hits)} genes corroborated by both systems):")
        for g in cross_hits[:5]:
            sources = ", ".join(g.get("sources", []))
            print(f"  {g['symbol']}: CASCADE {g['source_count']} sources ({sources}) + RegNetAgents target")

    errors = result.get("errors", {})
    if errors:
        print(f"\nPartial data — errors: {errors}")
    else:
        print("\nErrors: {} (both systems returned data)")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).parent.parent))
    asyncio.run(main())
