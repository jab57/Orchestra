"""
MYC Cross-Cell-Type Context Comparison

Validates Orchestra's compare_cell_contexts tool using MYC across three
cell types: cd4_t_cells, cd8_t_cells, epithelial_cell.

Expected result: MYC is a hub_regulator in all three cell types (RegNetAgents)
and essential in all three (DepMap/LINCS via CASCADE) — demonstrating conserved
evidence for a canonical oncogene across immune and epithelial contexts.

This case is the biological validation for Issue #11 (compare_cell_contexts).

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

sys.path.insert(0, str(Path(__file__).parent.parent))
from orchestra_langgraph_workflow import OrchestraWorkflow


async def main() -> None:
    cell_types = ["cd4_t_cells", "cd8_t_cells", "epithelial_cell"]

    print("Orchestra: MYC Cross-Cell-Type Context Comparison")
    print(f"Gene: MYC | Cell types: {', '.join(cell_types)}")
    print("-" * 60)

    workflow = OrchestraWorkflow()
    result = await workflow.run_analysis(
        gene="MYC",
        cell_type="",
        analysis_type="cell_context_comparison",
        cell_types=cell_types,
    )

    print(result.get("final_report", "(no report generated)"))

    synthesis = result.get("synthesis", {})
    conservation_scores = synthesis.get("conservation_scores", {})
    if conservation_scores:
        print("\nConservation summary:")
        for src, cons in conservation_scores.items():
            print(f"  {src}: {cons['label']} ({cons['count']}/{cons['n']})")

    errors = result.get("errors", {})
    if errors:
        print(f"\nPartial data — errors: {errors}")
    else:
        print("\nErrors: {} (all systems returned data)")


if __name__ == "__main__":
    asyncio.run(main())
