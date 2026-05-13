"""
APC Effector Analysis

Demonstrates Orchestra's effector/scaffold path with the canonical APC→CTNNB1 use case.

APC is a tumor suppressor scaffold protein in the Wnt/β-catenin pathway. It has no
direct transcriptional targets, so a standard perturbation query returns empty results.
Orchestra detects this automatically and reroutes:

  1. CASCADE get_gene_metadata → APC classified as effector/scaffold
  2. CASCADE get_protein_interactions → STRING PPI; CTNNB1 identified (score 0.999)
  3. CASCADE comprehensive_perturbation_analysis(CTNNB1) → downstream genes
  4. RegNetAgents comprehensive_gene_analysis(CTNNB1) → network context, pathway enrichment

This is the strict-necessity case: neither child server alone completes the analysis.
CASCADE dead-ends on empty perturbation; RegNetAgents has no PPI to bridge APC → CTNNB1.

Requirements: RegNetAgents and CASCADE must be installed and accessible.
See docs/installation.md for setup instructions.
"""

import asyncio
import sys
import io
from pathlib import Path

# UTF-8 stdout for Windows terminals that default to cp1252
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

from dotenv import load_dotenv
load_dotenv(Path(__file__).parent.parent / ".env")

from orchestra_langgraph_workflow import OrchestraWorkflow


async def main() -> None:
    print("Orchestra: APC Effector Analysis")
    print("Gene: APC | Cell type: epithelial_cell | Path: effector")
    print("-" * 60)

    workflow = OrchestraWorkflow()
    result = await workflow.run_analysis(
        gene="APC",
        cell_type="epithelial_cell",
        analysis_type="causal_chain",
    )

    print(result.get("final_report", "(no report generated)"))

    synthesis = result.get("synthesis", {})
    tf_partner = synthesis.get("tf_partner")
    if tf_partner:
        print(f"\nTF partner identified: {tf_partner}")

    errors = result.get("errors", {})
    if errors:
        print(f"\nPartial data — errors: {errors}")
    else:
        print("\nErrors: {} (both systems returned data)")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).parent.parent))
    asyncio.run(main())
