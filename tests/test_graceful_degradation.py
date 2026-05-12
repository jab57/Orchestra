"""
Unit tests for graceful degradation — behavior when one child server fails.

Uses mock clients injected directly into the workflow; no running child servers
required. Validates that the workflow correctly:
  - continues and populates errors{} when one gather arm raises an exception
  - sets the appropriate *_available flag in synthesis
  - emits the corresponding unavailability warning in the formatted report
  - handles effector path PPI failure without crashing
  - handles validation path RegNetAgents failure while still returning CASCADE candidates
"""

import asyncio
import pytest
from unittest.mock import AsyncMock, MagicMock

from orchestra_langgraph_workflow import OrchestraWorkflow, OrchestraState


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_workflow() -> OrchestraWorkflow:
    wf = object.__new__(OrchestraWorkflow)
    wf.use_llm = False
    wf.llm_available = False
    wf.llm_client = None
    wf.ollama_model = "llama3.1:8b"
    wf.ollama_temperature = 0.3
    wf.ollama_max_tokens = 2000
    wf._cascade = None
    wf._regnetagents = None
    wf.graph = None
    return wf


def _base_state(**overrides) -> dict:
    base = {
        "gene": "TP53",
        "cell_type": "epithelial_cell",
        "analysis_type": "causal_chain",
        "analysis_depth": "comprehensive",
        "gene_role": "master_regulator",
        "tf_partner": None,
        "network_analysis": None,
        "perturbation_result": None,
        "ppi_interactions": None,
        "validated_targets": None,
        "completed_steps": [],
        "errors": {},
        "final_report": None,
        "synthesis": None,
    }
    base.update(overrides)
    return base


_MOCK_NETWORK = {
    "target_analysis": {
        "cascade_targets": [
            {"gene_symbol": "CDKN1A"},
            {"gene_symbol": "MDM2"},
        ]
    }
}

_MOCK_PERTURBATION = {
    "evidence_synthesis": {
        "key_findings": ["CDKN1A downregulated in 3 sources"],
        "multi_source_genes": [
            {"symbol": "CDKN1A", "source_count": 3, "sources": ["lincs", "depmap", "dorothea"]}
        ],
        "source_agreements": ["CDKN1A: network + LINCS"],
        "source_disagreements": [],
    }
}


# ---------------------------------------------------------------------------
# TF path — CASCADE fails
# ---------------------------------------------------------------------------

class TestTfPathCascadeDown:
    async def test_errors_dict_contains_perturbation_key(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(side_effect=TimeoutError("cascade timeout"))
        wf._regnetagents = MagicMock()
        wf._regnetagents.call_tool = AsyncMock(return_value=_MOCK_NETWORK)

        state = await wf._run_tf_path(_base_state())
        assert "perturbation" in state["errors"]

    async def test_network_analysis_populated_from_regnetagents(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(side_effect=TimeoutError("cascade timeout"))
        wf._regnetagents = MagicMock()
        wf._regnetagents.call_tool = AsyncMock(return_value=_MOCK_NETWORK)

        state = await wf._run_tf_path(_base_state())
        assert state["network_analysis"] == _MOCK_NETWORK

    async def test_perturbation_result_absent(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(side_effect=TimeoutError("cascade timeout"))
        wf._regnetagents = MagicMock()
        wf._regnetagents.call_tool = AsyncMock(return_value=_MOCK_NETWORK)

        state = await wf._run_tf_path(_base_state())
        assert not state["perturbation_result"]

    async def test_synthesis_cascade_available_false(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(side_effect=TimeoutError("cascade timeout"))
        wf._regnetagents = MagicMock()
        wf._regnetagents.call_tool = AsyncMock(return_value=_MOCK_NETWORK)

        state = await wf._run_tf_path(_base_state())
        state = wf._synthesize_tf_path(state)
        assert state["synthesis"]["cascade_available"] is False

    async def test_report_warns_cascade_unavailable(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(side_effect=TimeoutError("cascade timeout"))
        wf._regnetagents = MagicMock()
        wf._regnetagents.call_tool = AsyncMock(return_value=_MOCK_NETWORK)

        state = await wf._run_tf_path(_base_state())
        state = wf._synthesize_tf_path(state)
        report = "\n".join(wf._format_tf_report(state["synthesis"]))
        assert "CASCADE unavailable" in report


# ---------------------------------------------------------------------------
# TF path — RegNetAgents fails
# ---------------------------------------------------------------------------

class TestTfPathRegNetAgentsDown:
    async def test_errors_dict_contains_network_key(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(return_value=_MOCK_PERTURBATION)
        wf._regnetagents = MagicMock()
        wf._regnetagents.call_tool = AsyncMock(side_effect=ConnectionError("regnetagents down"))

        state = await wf._run_tf_path(_base_state())
        assert "network" in state["errors"]

    async def test_perturbation_result_populated_from_cascade(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(return_value=_MOCK_PERTURBATION)
        wf._regnetagents = MagicMock()
        wf._regnetagents.call_tool = AsyncMock(side_effect=ConnectionError("regnetagents down"))

        state = await wf._run_tf_path(_base_state())
        assert state["perturbation_result"] == _MOCK_PERTURBATION

    async def test_synthesis_regnetagents_available_false(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(return_value=_MOCK_PERTURBATION)
        wf._regnetagents = MagicMock()
        wf._regnetagents.call_tool = AsyncMock(side_effect=ConnectionError("regnetagents down"))

        state = await wf._run_tf_path(_base_state())
        state = wf._synthesize_tf_path(state)
        assert state["synthesis"]["regnetagents_available"] is False

    async def test_report_warns_regnetagents_unavailable(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(return_value=_MOCK_PERTURBATION)
        wf._regnetagents = MagicMock()
        wf._regnetagents.call_tool = AsyncMock(side_effect=ConnectionError("regnetagents down"))

        state = await wf._run_tf_path(_base_state())
        state = wf._synthesize_tf_path(state)
        report = "\n".join(wf._format_tf_report(state["synthesis"]))
        assert "RegNetAgents unavailable" in report

    async def test_cascade_key_findings_still_present_in_report(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(return_value=_MOCK_PERTURBATION)
        wf._regnetagents = MagicMock()
        wf._regnetagents.call_tool = AsyncMock(side_effect=ConnectionError("regnetagents down"))

        state = await wf._run_tf_path(_base_state())
        state = wf._synthesize_tf_path(state)
        report = "\n".join(wf._format_tf_report(state["synthesis"]))
        assert "CDKN1A" in report


# ---------------------------------------------------------------------------
# Effector path — CASCADE PPI fails
# ---------------------------------------------------------------------------

class TestEffectorPathPpiDown:
    async def test_errors_dict_contains_ppi_key(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(side_effect=TimeoutError("ppi timeout"))
        wf._regnetagents = MagicMock()

        state = _base_state(gene="APC", gene_role="effector")
        state = await wf._run_effector_path(state)
        assert "ppi" in state["errors"]

    async def test_completed_steps_still_updated(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(side_effect=TimeoutError("ppi timeout"))
        wf._regnetagents = MagicMock()

        state = _base_state(gene="APC", gene_role="effector")
        state = await wf._run_effector_path(state)
        assert "run_effector_path" in state["completed_steps"]

    async def test_no_tf_partner_set(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(side_effect=TimeoutError("ppi timeout"))
        wf._regnetagents = MagicMock()

        state = _base_state(gene="APC", gene_role="effector")
        state = await wf._run_effector_path(state)
        assert not state.get("tf_partner")


# ---------------------------------------------------------------------------
# Validation path — RegNetAgents fails, CASCADE still contributes candidates
# ---------------------------------------------------------------------------

class TestValidationPathRegNetAgentsDown:
    def _cascade_side_effects(self):
        """Return three sequential CASCADE responses for the validation path."""
        return [
            # therapeutic_target_discovery
            {
                "therapeutic_targets": [
                    {"target": "BRD4", "priority": "high", "reason": "super-enhancer"}
                ]
            },
            # get_protein_interactions
            {"interactions": []},
            # comprehensive_perturbation_analysis for BRD4 (top-1 validation)
            {
                "evidence_synthesis": {
                    "key_findings": ["super-enhancer present"],
                    "multi_source_genes": [],
                    "source_agreements": [],
                    "source_disagreements": [],
                }
            },
        ]

    async def test_errors_dict_contains_network_key(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(side_effect=self._cascade_side_effects())
        wf._regnetagents = MagicMock()
        wf._regnetagents.call_tool = AsyncMock(side_effect=ConnectionError("regnetagents down"))

        state = _base_state(analysis_type="therapeutic_validation", gene="MYC",
                            cell_type="cd4_t_cells")
        state = await wf._run_validation_path(state)
        assert "network" in state["errors"]

    async def test_cascade_candidates_still_returned(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(side_effect=self._cascade_side_effects())
        wf._regnetagents = MagicMock()
        wf._regnetagents.call_tool = AsyncMock(side_effect=ConnectionError("regnetagents down"))

        state = _base_state(analysis_type="therapeutic_validation", gene="MYC",
                            cell_type="cd4_t_cells")
        state = await wf._run_validation_path(state)
        targets = state.get("validated_targets") or []
        assert any(t["gene"] == "BRD4" for t in targets)

    async def test_synthesis_regnetagents_available_false(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(side_effect=self._cascade_side_effects())
        wf._regnetagents = MagicMock()
        wf._regnetagents.call_tool = AsyncMock(side_effect=ConnectionError("regnetagents down"))

        state = _base_state(analysis_type="therapeutic_validation", gene="MYC",
                            cell_type="cd4_t_cells")
        state = await wf._run_validation_path(state)
        state = wf._synthesize_validation_path(state)
        assert state["synthesis"]["regnetagents_available"] is False

    async def test_validation_report_warns_regnetagents_unavailable(self):
        wf = _make_workflow()
        wf._cascade = MagicMock()
        wf._cascade.call_tool = AsyncMock(side_effect=self._cascade_side_effects())
        wf._regnetagents = MagicMock()
        wf._regnetagents.call_tool = AsyncMock(side_effect=ConnectionError("regnetagents down"))

        state = _base_state(analysis_type="therapeutic_validation", gene="MYC",
                            cell_type="cd4_t_cells")
        state = await wf._run_validation_path(state)
        state = wf._synthesize_validation_path(state)
        report = "\n".join(wf._format_validation_report(state["synthesis"]))
        assert "RegNetAgents unavailable" in report
