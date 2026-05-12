"""
Integration tests for the APC effector path.

Requires CASCADE and RegNetAgents running as child servers.
Skip unless ORCHESTRA_INTEGRATION_TESTS=1.

Validates:
  - APC routes to effector path (not a TF, no direct targets)
  - CTNNB1 is identified as the TF partner via PPI
  - Both systems return data without errors
  - The final report names CTNNB1 and references the effector routing
"""

import os
import pytest

from orchestra_langgraph_workflow import OrchestraWorkflow

pytestmark = pytest.mark.skipif(
    not os.getenv("ORCHESTRA_INTEGRATION_TESTS"),
    reason="requires running child servers; set ORCHESTRA_INTEGRATION_TESTS=1",
)


async def test_apc_routes_to_effector_path():
    """APC is a scaffold protein — it must route to the effector path, not tf_path."""
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("APC", "epithelial_cell")
    assert result["synthesis"]["routing"] == "effector"


async def test_apc_identifies_ctnnb1_as_tf_partner():
    """CTNNB1 is APC's primary PPI partner and should be selected as TF partner."""
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("APC", "epithelial_cell")
    assert result["tf_partner"] == "CTNNB1"


async def test_apc_no_errors_from_either_server():
    """Both child servers must complete without error for the APC use case."""
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("APC", "epithelial_cell")
    assert result["errors"] == {}


async def test_apc_completed_steps_include_effector_path():
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("APC", "epithelial_cell")
    assert "run_effector_path" in result["completed_steps"]


async def test_apc_final_report_is_non_empty_string():
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("APC", "epithelial_cell")
    assert isinstance(result["final_report"], str)
    assert len(result["final_report"]) > 0


async def test_apc_report_mentions_ctnnb1():
    """CTNNB1 must appear in the final report as the TF partner."""
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("APC", "epithelial_cell")
    assert "CTNNB1" in result["final_report"]


async def test_apc_report_mentions_effector_routing():
    """The report header must identify this as an effector/scaffold analysis."""
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("APC", "epithelial_cell")
    assert "effector" in result["final_report"].lower()


async def test_apc_synthesis_has_cascade_and_regnetagents_available():
    """Both systems must be flagged as available for the strict-necessity argument."""
    wf = OrchestraWorkflow()
    result = await wf.run_analysis("APC", "epithelial_cell")
    synthesis = result["synthesis"]
    assert synthesis["cascade_available"] is True
    assert synthesis["regnetagents_available"] is True
