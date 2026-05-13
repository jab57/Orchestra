"""
MYC Therapeutic Target Validation (BRD4 complementarity case)

Demonstrates Orchestra's validate_therapeutic_targets path using MYC in cd4_t_cells.

Orchestra queries three independent candidate sources:
  1. RegNetAgents: upstream TF regulators ranked by PageRank centrality
  2. CASCADE: therapeutic_target_discovery (super-enhancers, drug database)
  3. CASCADE: STRING PPI top interactors

Top candidates are validated via CASCADE comprehensive_perturbation_analysis and
scored against 7 independent evidence sources (corroboration table).

Key result: BRD4 scores 1/7 — found by CASCADE super-enhancer analysis, absent from
the ARACNe regulatory network. This is biologically expected: BRD4 is a BET bromodomain
co-activator that acts through chromatin-level super-enhancer binding, not direct
transcriptional regulation. ARACNe cannot detect chromatin-level mechanisms.

Orchestra presents both views together:
  - CASCADE: "MYC has super-enhancers in 32 cell types → BET inhibitor sensitivity"
  - RegNetAgents: "BRD4 absent from the ARACNe TF network"
The absence is as informative as the presence — it tells you why the mechanism is
epigenetic, not transcriptional.

BRD4→MYC via super-enhancers is published (Lovén et al. 2013, Cell).
BET inhibitors (JQ1, OTX015) are in clinical trials.

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
    print("Orchestra: MYC Therapeutic Target Validation")
    print("Gene: MYC | Cell type: cd4_t_cells | Path: therapeutic_validation")
    print("-" * 60)

    workflow = OrchestraWorkflow()
    result = await workflow.run_analysis(
        gene="MYC",
        cell_type="cd4_t_cells",
        analysis_type="therapeutic_validation",
    )

    print(result.get("final_report", "(no report generated)"))

    synthesis = result.get("synthesis", {})
    evidence_table = synthesis.get("evidence_table", [])
    if evidence_table:
        print("\nCorroboration summary (top candidates):")
        for row in evidence_table[:5]:
            score = f"{row['corroboration_count']}/{row['corroboration_denominator']}"
            se = "✓" if row.get("super_enhancer") else "-"
            pagerank = "✓" if row.get("pagerank_rank") else "-"
            print(f"  {row['gene']}: score={score}, super_enhancer={se}, pagerank={pagerank}")

    errors = result.get("errors", {})
    if errors:
        print(f"\nPartial data — errors: {errors}")
    else:
        print("\nErrors: {} (both systems returned data)")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).parent.parent))
    asyncio.run(main())
