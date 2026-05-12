"""
Integration tests for the TP53 TF causal chain path.

Requires CASCADE and RegNetAgents running as child servers.
Skip unless ORCHESTRA_INTEGRATION_TESTS=1.

Validates:
  - TP53 routes to the TF path (master_regulator classification)
  - Both systems return data without errors
  - CDKN1A appears in CASCADE corroborated targets (canonical TP53 target)
  - Cross-system corroboration hits are detected
  - The final report is well-formed and contains expected content
"""

import os
import pytest

from orchestra_langgraph_workflow import OrchestraWorkflow

pytestmark = pytest.mark.skipif(
    not os.getenv("ORCHESTRA_INTEGRATION_TESTS"),
    reason="requires running child servers; set ORCHESTRA_INTEGRATION_TESTS=1",
)


async def test_tp53_routes_to_tf_path():
    """TP53 is a canonical master regulator — must route to tf_path."""
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("TP53", "epithelial_cell")
    assert result["synthesis"]["routing"] == "tf"


async def test_tp53_gene_role_is_tf_type():
    """CASCADE must classify TP53 as a TF-type gene."""
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("TP53", "epithelial_cell")
    assert result["gene_role"] in ("master_regulator", "transcription_factor")


async def test_tp53_no_errors_from_either_server():
    """Both child servers must complete without error for TP53."""
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("TP53", "epithelial_cell")
    assert result["errors"] == {}


async def test_tp53_completed_steps_include_tf_path():
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("TP53", "epithelial_cell")
    assert "run_tf_path" in result["completed_steps"]


async def test_tp53_cdkn1a_in_corroborated_targets():
    """CDKN1A is the canonical TP53 transcriptional target; must appear in CASCADE corroborated targets."""
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("TP53", "epithelial_cell")
    corroborated = result["synthesis"].get("corroborated_targets", [])
    symbols = [g["symbol"] for g in corroborated]
    assert "CDKN1A" in symbols


async def test_tp53_cross_system_hits_detected():
    """At least one gene must appear in both RegNetAgents network targets and CASCADE
    experimental data — this is the core cross-system corroboration result."""
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("TP53", "epithelial_cell")
    assert len(result["synthesis"]["cross_system_hits"]) > 0


async def test_tp53_final_report_is_non_empty_string():
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("TP53", "epithelial_cell")
    assert isinstance(result["final_report"], str)
    assert len(result["final_report"]) > 0


async def test_tp53_report_mentions_tp53_and_cell_type():
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("TP53", "epithelial_cell")
    assert "TP53" in result["final_report"]
    assert "epithelial_cell" in result["final_report"]


async def test_tp53_synthesis_has_both_systems_available():
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("TP53", "epithelial_cell")
    synthesis = result["synthesis"]
    assert synthesis["cascade_available"] is True
    assert synthesis["regnetagents_available"] is True
