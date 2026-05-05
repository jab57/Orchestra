"""
Orchestra test suite.

Unit tests cover routing logic, synthesis, and report formatting — all
runnable without live child servers. Integration tests (requiring running
RegNetAgents and CASCADE servers) are skipped unless
ORCHESTRA_INTEGRATION_TESTS=1 is set in the environment.
"""

import os
import pytest
from orchestra_langgraph_workflow import OrchestraWorkflow, OrchestraState


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_workflow() -> OrchestraWorkflow:
    """Instantiate OrchestraWorkflow bypassing LLM init side-effects."""
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


def _tf_state(**overrides) -> dict:
    base = {
        "gene": "TP53",
        "cell_type": "epithelial_cell",
        "analysis_type": "causal_chain",
        "analysis_depth": "comprehensive",
        "gene_role": "master_regulator",
        "tf_partner": None,
        "network_analysis": {},
        "perturbation_result": {},
        "ppi_interactions": None,
        "validated_targets": None,
        "completed_steps": [],
        "errors": {},
        "final_report": None,
        "synthesis": None,
    }
    base.update(overrides)
    return base


# ---------------------------------------------------------------------------
# State schema
# ---------------------------------------------------------------------------

class TestStateSchema:
    def test_required_fields_present(self):
        required = [
            "gene", "cell_type", "analysis_type", "analysis_depth",
            "gene_role", "network_analysis", "perturbation_result",
            "completed_steps", "errors", "final_report",
        ]
        annotations = OrchestraState.__annotations__
        for field in required:
            assert field in annotations, f"Missing state field: {field}"


# ---------------------------------------------------------------------------
# Workflow instantiation
# ---------------------------------------------------------------------------

class TestWorkflowInstantiation:
    def test_instantiation(self):
        wf = OrchestraWorkflow()
        assert wf is not None

    def test_llm_disabled_by_default(self):
        wf = OrchestraWorkflow()
        assert wf.use_llm is False


# ---------------------------------------------------------------------------
# Routing decision
# ---------------------------------------------------------------------------

class TestRoutingDecision:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def test_master_regulator_routes_tf(self, wf):
        state = _tf_state(gene_role="master_regulator")
        assert wf._routing_decision(state) == "tf_path"

    def test_transcription_factor_routes_tf(self, wf):
        state = _tf_state(gene_role="transcription_factor")
        assert wf._routing_decision(state) == "tf_path"

    def test_minor_regulator_routes_tf(self, wf):
        state = _tf_state(gene_role="minor_regulator")
        assert wf._routing_decision(state) == "tf_path"

    def test_effector_routes_effector(self, wf):
        state = _tf_state(gene_role="effector")
        assert wf._routing_decision(state) == "effector_path"

    def test_isolated_routes_effector(self, wf):
        state = _tf_state(gene_role="isolated")
        assert wf._routing_decision(state) == "effector_path"

    def test_none_role_defaults_effector(self, wf):
        state = _tf_state(gene_role=None)
        assert wf._routing_decision(state) == "effector_path"

    def test_therapeutic_validation_overrides_gene_role(self, wf):
        """therapeutic_validation analysis_type always routes to validation_path
        regardless of gene role."""
        state = _tf_state(analysis_type="therapeutic_validation",
                          gene_role="master_regulator")
        assert wf._routing_decision(state) == "validation_path"

    def test_therapeutic_validation_with_effector_role(self, wf):
        state = _tf_state(analysis_type="therapeutic_validation",
                          gene_role="effector")
        assert wf._routing_decision(state) == "validation_path"


# ---------------------------------------------------------------------------
# _extract_regnetagents_targets
# ---------------------------------------------------------------------------

class TestExtractRegNetagentsTargets:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def test_primary_path_cascade_targets(self, wf):
        network = {
            "target_analysis": {
                "cascade_targets": [
                    {"gene_symbol": "MYC"},
                    {"gene_symbol": "CDKN1A"},
                ]
            }
        }
        result = wf._extract_regnetagents_targets(network)
        assert "MYC" in result
        assert "CDKN1A" in result

    def test_fallback_top_level_targets(self, wf):
        network = {"targets": [{"gene_symbol": "TP53"}, {"symbol": "BRCA1"}]}
        result = wf._extract_regnetagents_targets(network)
        assert "TP53" in result
        assert "BRCA1" in result

    def test_fallback_string_list(self, wf):
        network = {"targets": ["MYC", "BCL2"]}
        result = wf._extract_regnetagents_targets(network)
        assert "MYC" in result
        assert "BCL2" in result

    def test_fallback_nested_block(self, wf):
        network = {
            "network_analysis": {
                "downstream_targets": [{"gene_symbol": "IRF4"}]
            }
        }
        result = wf._extract_regnetagents_targets(network)
        assert "IRF4" in result

    def test_empty_network(self, wf):
        assert wf._extract_regnetagents_targets({}) == set()

    def test_deduplicates_across_fields(self, wf):
        network = {
            "target_analysis": {"cascade_targets": [{"gene_symbol": "MYC"}]},
            "targets": [{"gene_symbol": "MYC"}],
        }
        result = wf._extract_regnetagents_targets(network)
        assert result == {"MYC"}

    def test_ignores_none_symbols(self, wf):
        network = {"targets": [{"gene_symbol": None}, {"gene_symbol": "BCL2"}]}
        result = wf._extract_regnetagents_targets(network)
        assert None not in result
        assert "BCL2" in result


# ---------------------------------------------------------------------------
# _synthesize_tf_path
# ---------------------------------------------------------------------------

class TestSynthesizeTfPath:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def _perturbation_with_findings(self, genes=None):
        genes = genes or ["MYC", "BCL2"]
        return {
            "evidence_synthesis": {
                "key_findings": ["MYC downregulated in 3 sources"],
                "multi_source_genes": [
                    {"symbol": g, "source_count": 2, "sources": ["lincs", "depmap"]}
                    for g in genes
                ],
                "source_agreements": ["MYC: network down + LINCS down"],
                "source_disagreements": [],
            }
        }

    def test_routing_key_is_tf(self, wf):
        state = _tf_state(
            perturbation_result=self._perturbation_with_findings(),
            network_analysis={},
        )
        result = wf._synthesize_tf_path(state)
        assert result["synthesis"]["routing"] == "tf"

    def test_cross_system_hits_when_overlap(self, wf):
        """Genes in both CASCADE multi_source_genes AND RegNetAgents targets
        should appear in cross_system_hits."""
        state = _tf_state(
            perturbation_result=self._perturbation_with_findings(["MYC", "BCL2"]),
            network_analysis={
                "target_analysis": {
                    "cascade_targets": [{"gene_symbol": "MYC"}]
                }
            },
        )
        result = wf._synthesize_tf_path(state)
        hits = result["synthesis"]["cross_system_hits"]
        assert any(h["symbol"] == "MYC" for h in hits)
        assert not any(h["symbol"] == "BCL2" for h in hits)

    def test_no_cross_system_hits_when_no_overlap(self, wf):
        state = _tf_state(
            perturbation_result=self._perturbation_with_findings(["CDKN1A"]),
            network_analysis={
                "target_analysis": {
                    "cascade_targets": [{"gene_symbol": "MYC"}]
                }
            },
        )
        result = wf._synthesize_tf_path(state)
        assert result["synthesis"]["cross_system_hits"] == []

    def test_errors_preserved(self, wf):
        state = _tf_state(
            errors={"network": "timeout"},
            perturbation_result={},
            network_analysis={},
        )
        result = wf._synthesize_tf_path(state)
        assert result["synthesis"]["errors"]["network"] == "timeout"

    def test_completed_steps_updated(self, wf):
        state = _tf_state(perturbation_result={}, network_analysis={})
        result = wf._synthesize_tf_path(state)
        assert "synthesize" in result["completed_steps"]


# ---------------------------------------------------------------------------
# _synthesize_effector_path
# ---------------------------------------------------------------------------

class TestSynthesizeEffectorPath:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def test_routing_key_is_effector(self, wf):
        state = _tf_state(
            gene_role="effector",
            tf_partner="CTNNB1",
            perturbation_result={},
            network_analysis={},
        )
        result = wf._synthesize_effector_path(state)
        assert result["synthesis"]["routing"] == "effector"

    def test_tf_partner_preserved(self, wf):
        state = _tf_state(
            gene_role="effector",
            tf_partner="CTNNB1",
            perturbation_result={},
            network_analysis={},
        )
        result = wf._synthesize_effector_path(state)
        assert result["synthesis"]["tf_partner"] == "CTNNB1"

    def test_gene_and_cell_type_preserved(self, wf):
        state = _tf_state(
            gene="APC",
            cell_type="epithelial_cell",
            gene_role="effector",
            tf_partner="CTNNB1",
            perturbation_result={},
            network_analysis={},
        )
        result = wf._synthesize_effector_path(state)
        assert result["synthesis"]["gene"] == "APC"
        assert result["synthesis"]["cell_type"] == "epithelial_cell"


# ---------------------------------------------------------------------------
# _synthesize_validation_path
# ---------------------------------------------------------------------------

class TestSynthesizeValidationPath:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def _state_with_candidates(self, candidates):
        return _tf_state(
            analysis_type="therapeutic_validation",
            validated_targets=candidates,
            network_analysis={},
        )

    def test_routing_key_is_validation(self, wf):
        state = self._state_with_candidates([])
        result = wf._synthesize_validation_path(state)
        assert result["synthesis"]["routing"] == "validation"

    def test_candidates_preserved(self, wf):
        candidates = [
            {"gene": "BRD4", "source": "cascade_drug_discovery"},
            {"gene": "CDK9", "source": "regnetagents_pagerank"},
        ]
        state = self._state_with_candidates(candidates)
        result = wf._synthesize_validation_path(state)
        assert len(result["synthesis"]["validated_targets"]) == 2

    def test_empty_candidates(self, wf):
        state = self._state_with_candidates([])
        result = wf._synthesize_validation_path(state)
        assert result["synthesis"]["validated_targets"] == []


# ---------------------------------------------------------------------------
# Report formatting — TF path
# ---------------------------------------------------------------------------

class TestFormatTfReport:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def _synthesis(self, **overrides):
        base = {
            "gene": "TP53",
            "cell_type": "epithelial_cell",
            "routing": "tf",
            "gene_role": "master_regulator",
            "cascade_key_findings": ["TP53 knockdown affects 45 downstream genes"],
            "corroborated_targets": [
                {"symbol": "CDKN1A", "source_count": 3, "sources": ["lincs", "depmap", "dorothea"]}
            ],
            "cross_system_hits": [
                {"symbol": "CDKN1A", "source_count": 3, "sources": ["lincs", "depmap", "dorothea"]}
            ],
            "source_agreements": ["CDKN1A: network down + LINCS down"],
            "network_context": {},
            "regnetagents_target_count": 45,
            "errors": {},
        }
        base.update(overrides)
        return base

    def test_header_contains_gene_and_cell_type(self, wf):
        lines = wf._format_tf_report(self._synthesis())
        header = "\n".join(lines)
        assert "TP53" in header
        assert "epithelial_cell" in header

    def test_cross_system_hits_section_present(self, wf):
        lines = wf._format_tf_report(self._synthesis())
        text = "\n".join(lines)
        assert "Cross-System Corroboration" in text
        assert "CDKN1A" in text

    def test_no_cross_system_hits_message(self, wf):
        lines = wf._format_tf_report(self._synthesis(cross_system_hits=[]))
        text = "\n".join(lines)
        assert "No cross-system hits" in text

    def test_key_findings_included(self, wf):
        lines = wf._format_tf_report(self._synthesis())
        text = "\n".join(lines)
        assert "45 downstream genes" in text

    def test_errors_section_when_present(self, wf):
        lines = wf._format_tf_report(self._synthesis(errors={"network": "timeout"}))
        text = "\n".join(lines)
        assert "Partial Data Warnings" in text
        assert "timeout" in text

    def test_returns_list_of_strings(self, wf):
        lines = wf._format_tf_report(self._synthesis())
        assert isinstance(lines, list)
        assert all(isinstance(l, str) for l in lines)


# ---------------------------------------------------------------------------
# Report formatting — effector path
# ---------------------------------------------------------------------------

class TestFormatEffectorReport:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def _synthesis(self, **overrides):
        base = {
            "gene": "APC",
            "cell_type": "epithelial_cell",
            "routing": "effector",
            "tf_partner": "CTNNB1",
            "cascade_key_findings": ["CTNNB1 overexpression activates 78 Wnt targets"],
            "corroborated_targets": [],
            "source_agreements": [],
            "source_disagreements": [],
            "network_context": {},
            "errors": {},
        }
        base.update(overrides)
        return base

    def test_routing_label(self, wf):
        text = "\n".join(wf._format_effector_report(self._synthesis()))
        assert "effector" in text

    def test_tf_partner_present(self, wf):
        text = "\n".join(wf._format_effector_report(self._synthesis()))
        assert "CTNNB1" in text

    def test_key_findings_included(self, wf):
        text = "\n".join(wf._format_effector_report(self._synthesis()))
        assert "Wnt targets" in text

    def test_missing_tf_partner_shows_not_identified(self, wf):
        text = "\n".join(wf._format_effector_report(
            self._synthesis(tf_partner=None)
        ))
        assert "not identified" in text


# ---------------------------------------------------------------------------
# Report formatting — validation path
# ---------------------------------------------------------------------------

class TestFormatValidationReport:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def _synthesis(self, candidates, **overrides):
        base = {
            "gene": "MYC",
            "cell_type": "cd4_t_cells",
            "routing": "validation",
            "validated_targets": candidates,
            "network_context": {},
            "errors": {},
        }
        base.update(overrides)
        return base

    def test_candidate_count_in_header(self, wf):
        candidates = [
            {"gene": "BRD4", "source": "cascade_drug_discovery",
             "key_findings": ["BRD4 knockdown reduces MYC"], "multi_source_genes": []},
        ]
        text = "\n".join(wf._format_validation_report(self._synthesis(candidates)))
        assert "1" in text

    def test_candidate_names_appear(self, wf):
        candidates = [
            {"gene": "BRD4", "source": "cascade_drug_discovery",
             "key_findings": [], "multi_source_genes": []},
            {"gene": "CDK9", "source": "regnetagents_pagerank",
             "key_findings": [], "multi_source_genes": [],
             "pagerank": 0.05, "downstream_targets": 20},
        ]
        text = "\n".join(wf._format_validation_report(self._synthesis(candidates)))
        assert "BRD4" in text
        assert "CDK9" in text

    def test_no_candidates_message(self, wf):
        text = "\n".join(wf._format_validation_report(self._synthesis([])))
        assert "No therapeutic target candidates" in text

    def test_cascade_error_shown(self, wf):
        candidates = [
            {"gene": "BRD4", "source": "cascade_drug_discovery",
             "cascade_error": "timeout", "key_findings": [], "multi_source_genes": []},
        ]
        text = "\n".join(wf._format_validation_report(self._synthesis(candidates)))
        assert "timeout" in text

    def test_not_validated_message_for_candidate_beyond_top3(self, wf):
        """Candidates 4+ should note they were not run through CASCADE validation."""
        candidates = [
            {"gene": g, "source": "cascade_drug_discovery",
             "key_findings": [], "multi_source_genes": []}
            for g in ["BRD4", "CDK9", "CCNT1", "MED1"]
        ]
        text = "\n".join(wf._format_validation_report(self._synthesis(candidates)))
        assert "not run" in text


# ---------------------------------------------------------------------------
# Integration test — skipped unless ORCHESTRA_INTEGRATION_TESTS=1
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not os.getenv("ORCHESTRA_INTEGRATION_TESTS"),
    reason="requires running child servers; set ORCHESTRA_INTEGRATION_TESTS=1"
)
@pytest.mark.asyncio
async def test_run_analysis_returns_state():
    """End-to-end: run_analysis opens MCP connections and returns a state dict."""
    workflow = OrchestraWorkflow()
    result = await workflow.run_analysis(
        gene="TP53",
        cell_type="epithelial_cell",
    )
    assert isinstance(result, dict)
    assert "completed_steps" in result
