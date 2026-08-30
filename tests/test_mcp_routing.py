"""
Regression tests for the MCP tool name -> analysis_type routing bug.

orchestra_mcp_server.call_tool() dispatches by MCP tool name. Most tools have an
explicit elif branch that sets the correct internal analysis_type string, but
validate_therapeutic_targets and effector_analysis previously fell through to the
generic handler, which passed the raw tool name straight through as analysis_type.
_routing_decision() checks for "therapeutic_validation" (not
"validate_therapeutic_targets"), so validate_therapeutic_targets silently routed
to gene-role-based default routing (tf_path/effector_path) instead of
validation_path -- since project inception, per `git log -S`. effector_analysis had
no matching check in _routing_decision() at all, so it also silently fell through to
gene-role-based default routing instead of always taking the effector path, contrary
to its documented "scaffold/effector routing" purpose.

These tests call orchestra_mcp_server.call_tool() directly (not just
_routing_decision() in isolation) to catch a regression at the actual dispatch
layer, not just in the routing helper.
"""
import sys
import os
from unittest.mock import AsyncMock, MagicMock

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import mcp.server.lowlevel.server as _mcp_server_module
import orchestra_mcp_server


def _install_request_context():
    """MCP tool handlers read app.request_context (a contextvar) for session/request_id."""
    fake_session = MagicMock()
    fake_session.send_log_message = AsyncMock()
    ctx = _mcp_server_module.RequestContext(
        request_id="test-request",
        meta=None,
        session=fake_session,
        lifespan_context=None,
    )
    token = _mcp_server_module.request_ctx.set(ctx)
    return token


@pytest.fixture
def mock_run_analysis(monkeypatch):
    """Replace the module-level `workflow` with a mock capturing run_analysis() kwargs."""
    captured = {}

    async def fake_run_analysis(**kwargs):
        captured.update(kwargs)
        return {"final_report": "mocked report"}

    mock_workflow = MagicMock()
    mock_workflow.run_analysis = AsyncMock(side_effect=fake_run_analysis)
    monkeypatch.setattr(orchestra_mcp_server, "workflow", mock_workflow)
    return captured


@pytest.mark.asyncio
async def test_validate_therapeutic_targets_routes_to_therapeutic_validation(mock_run_analysis):
    token = _install_request_context()
    try:
        await orchestra_mcp_server.call_tool(
            "validate_therapeutic_targets", {"gene": "MYC", "cell_type": "epithelial_cell"}
        )
    finally:
        _mcp_server_module.request_ctx.reset(token)

    assert mock_run_analysis["analysis_type"] == "therapeutic_validation"
    assert mock_run_analysis["gene"] == "MYC"


@pytest.mark.asyncio
async def test_effector_analysis_routes_to_effector_analysis_type(mock_run_analysis):
    token = _install_request_context()
    try:
        await orchestra_mcp_server.call_tool(
            "effector_analysis", {"gene": "APC", "cell_type": "epithelial_cell"}
        )
    finally:
        _mcp_server_module.request_ctx.reset(token)

    assert mock_run_analysis["analysis_type"] == "effector_analysis"
    assert mock_run_analysis["gene"] == "APC"


@pytest.mark.asyncio
async def test_causal_chain_analysis_unaffected(mock_run_analysis):
    """causal_chain_analysis intentionally relies on default gene-role routing --
    confirm this fix didn't change its (already-correct) behavior."""
    token = _install_request_context()
    try:
        await orchestra_mcp_server.call_tool(
            "causal_chain_analysis", {"gene": "TP53", "cell_type": "epithelial_cell"}
        )
    finally:
        _mcp_server_module.request_ctx.reset(token)

    assert mock_run_analysis["analysis_type"] == "causal_chain_analysis"
    assert mock_run_analysis["gene"] == "TP53"


class TestRoutingDecisionEffectorAnalysis:
    """Direct _routing_decision unit coverage for the effector_analysis fix."""

    def test_effector_analysis_forces_effector_path_regardless_of_role(self):
        from orchestra_langgraph_workflow import OrchestraWorkflow

        wf = object.__new__(OrchestraWorkflow)
        for role in ("master_regulator", "transcription_factor", "minor_regulator", "effector", "isolated", None):
            decision = wf._routing_decision({"analysis_type": "effector_analysis", "gene_role": role})
            assert decision == "effector_path", f"role={role} -> {decision}"
