"""
Tests for Issue #13: compare_network_contexts (GREmLN vs TCGA network context comparison).

Unit tests run without live child servers.
Integration test is gated on ORCHESTRA_INTEGRATION_TESTS=1.
"""

import os
import pytest
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


def _nc_state(**overrides) -> dict:
    base = {
        "gene": "FOXM1",
        "cell_type": "epithelial_cell",
        "cancer_type": "hnsc",
        "analysis_type": "network_comparison",
        "analysis_depth": "comprehensive",
        "gene_role": None,
        "ensembl_id": None,
        "tf_partner": None,
        "network_analysis": None,
        "pathway_enrichment": None,
        "domain_insights": None,
        "perturbation_result": None,
        "ppi_interactions": None,
        "lincs_effects": None,
        "depmap_essentiality": None,
        "validated_targets": None,
        "causal_chain": None,
        "gene_signature": None,
        "master_regulators": None,
        "cell_types": None,
        "comparison_results": None,
        "network_comparison": None,
        "completed_steps": [],
        "errors": {},
        "final_report": None,
        "synthesis": None,
    }
    base.update(overrides)
    return base


def _context_result(
    conserved=None, pop_only=None, tumor_only=None,
    conserved_fraction=0.5, rewiring="moderate",
    pop_total=10, tumor_total=8,
    cascade_validation=None,
):
    return {
        "gene": "FOXM1",
        "population_averaged_context": "epithelial_cell",
        "tumor_state_context": "tcga_hnsc",
        "regulators": {
            "population_averaged_total": pop_total,
            "tumor_state_total": tumor_total,
            "conserved": conserved or ["E2F1", "MYC", "BRCA1"],
            "conserved_count": len(conserved or ["E2F1", "MYC", "BRCA1"]),
            "conserved_fraction": conserved_fraction,
            "population_averaged_only": pop_only or ["TP53"],
            "tumor_state_only": tumor_only or ["EGFR", "KRAS"],
        },
        "targets": {
            "population_averaged_total": 20,
            "tumor_state_total": 18,
            "conserved_count": 12,
            "conserved_fraction": 0.6,
            "population_averaged_only": ["CDKN1A"],
            "tumor_state_only": ["MMP9"],
        },
        "interpretation": {
            "regulatory_rewiring": rewiring,
            "conserved_fraction_regulators": conserved_fraction,
            "tumor_specific_regulator_count": len(tumor_only or ["EGFR", "KRAS"]),
        },
        "cascade_validation": cascade_validation or {},
    }


def _cascade_result(gene: str, has_lincs=True, has_depmap=False) -> dict:
    findings = []
    if has_lincs:
        findings.append(f"LINCS: {gene} knockdown suppresses downstream targets")
    if has_depmap:
        findings.append(f"DepMap: {gene} is essential in cancer cell lines")
    return {
        "evidence_synthesis": {
            "key_findings": findings,
            "multi_source_genes": [],
        }
    }


# ---------------------------------------------------------------------------
# Routing
# ---------------------------------------------------------------------------

class TestNetworkComparisonRouting:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def test_network_comparison_routes_correctly(self, wf):
        state = _nc_state(analysis_type="network_comparison")
        assert wf._routing_decision(state) == "network_comparison_path"

    def test_network_comparison_overrides_gene_role(self, wf):
        state = _nc_state(analysis_type="network_comparison", gene_role="master_regulator")
        assert wf._routing_decision(state) == "network_comparison_path"

    def test_other_types_unaffected(self, wf):
        state = _nc_state(analysis_type="causal_chain", gene_role="master_regulator")
        assert wf._routing_decision(state) == "tf_path"

    def test_cell_context_comparison_unaffected(self, wf):
        state = _nc_state(analysis_type="cell_context_comparison")
        assert wf._routing_decision(state) == "comparison_path"


# ---------------------------------------------------------------------------
# OrchestraState has new fields
# ---------------------------------------------------------------------------

class TestStateSchemaNetworkComparison:
    def test_cancer_type_field_in_state(self):
        assert "cancer_type" in OrchestraState.__annotations__

    def test_network_comparison_field_in_state(self):
        assert "network_comparison" in OrchestraState.__annotations__


# ---------------------------------------------------------------------------
# _synthesize_network_comparison_path
# ---------------------------------------------------------------------------

class TestSynthesizeNetworkComparisonPath:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def test_routing_is_network_comparison(self, wf):
        state = _nc_state(
            network_comparison=_context_result(),
        )
        result = wf._synthesize_network_comparison_path(state)
        assert result["synthesis"]["routing"] == "network_comparison"

    def test_rewiring_classification_propagated(self, wf):
        state = _nc_state(
            network_comparison=_context_result(rewiring="high", conserved_fraction=0.2),
        )
        result = wf._synthesize_network_comparison_path(state)
        assert result["synthesis"]["rewiring_classification"] == "high"

    def test_conserved_regulators_listed(self, wf):
        state = _nc_state(
            network_comparison=_context_result(conserved=["E2F1", "MYC"]),
        )
        result = wf._synthesize_network_comparison_path(state)
        assert "E2F1" in result["synthesis"]["conserved_regulators"]
        assert "MYC" in result["synthesis"]["conserved_regulators"]

    def test_tumor_only_regulators_listed(self, wf):
        state = _nc_state(
            network_comparison=_context_result(tumor_only=["EGFR", "KRAS"]),
        )
        result = wf._synthesize_network_comparison_path(state)
        assert "EGFR" in result["synthesis"]["tumor_state_only_regulators"]

    def test_cascade_validated_tier_assigned(self, wf):
        cascade_val = {"E2F1": _cascade_result("E2F1", has_lincs=True)}
        state = _nc_state(
            network_comparison=_context_result(
                conserved=["E2F1"], cascade_validation=cascade_val
            ),
        )
        result = wf._synthesize_network_comparison_path(state)
        vc = result["synthesis"]["validated_conserved"]
        assert len(vc) == 1
        assert vc[0]["gene"] == "E2F1"
        assert vc[0]["tier"] == "conserved_cascade_validated"

    def test_not_validated_tier_when_no_cascade_evidence(self, wf):
        cascade_val = {"E2F1": {"evidence_synthesis": {"key_findings": []}}}
        state = _nc_state(
            network_comparison=_context_result(
                conserved=["E2F1"], cascade_validation=cascade_val
            ),
        )
        result = wf._synthesize_network_comparison_path(state)
        vc = result["synthesis"]["validated_conserved"]
        assert vc[0]["tier"] == "conserved_not_validated"

    def test_cascade_error_gives_not_validated(self, wf):
        cascade_val = {"E2F1": {"error": "timeout"}}
        state = _nc_state(
            network_comparison=_context_result(
                conserved=["E2F1"], cascade_validation=cascade_val
            ),
        )
        result = wf._synthesize_network_comparison_path(state)
        vc = result["synthesis"]["validated_conserved"]
        assert vc[0]["tier"] == "conserved_not_validated"
        assert vc[0]["cascade_error"] == "timeout"

    def test_empty_network_comparison_produces_empty_synthesis(self, wf):
        state = _nc_state(network_comparison=None)
        result = wf._synthesize_network_comparison_path(state)
        syn = result["synthesis"]
        assert syn["conserved_regulators"] == []
        assert syn["validated_conserved"] == []

    def test_regnetagents_available_false_when_error(self, wf):
        state = _nc_state(
            network_comparison=None,
            errors={"network_comparison": "compare_network_contexts failed"},
        )
        result = wf._synthesize_network_comparison_path(state)
        assert result["synthesis"]["regnetagents_available"] is False

    def test_top_3_conserved_validated(self, wf):
        conserved = ["A", "B", "C", "D", "E"]
        cascade_val = {g: _cascade_result(g, has_lincs=True) for g in conserved[:3]}
        state = _nc_state(
            network_comparison=_context_result(
                conserved=conserved, cascade_validation=cascade_val
            ),
        )
        result = wf._synthesize_network_comparison_path(state)
        assert len(result["synthesis"]["validated_conserved"]) == 3


# ---------------------------------------------------------------------------
# _format_network_comparison_report
# ---------------------------------------------------------------------------

class TestFormatNetworkComparisonReport:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def _make_synthesis(self, validated_conserved=None, rewiring="moderate",
                        conserved_regs=None, tumor_only=None, pop_only=None,
                        regnetagents_available=True, cascade_available=True,
                        errors=None):
        return {
            "routing": "network_comparison",
            "gene": "FOXM1",
            "cell_type": "epithelial_cell",
            "cancer_type": "hnsc",
            "tumor_state_context": "tcga_hnsc",
            "rewiring_classification": rewiring,
            "reg_conserved_fraction": 0.5,
            "reg_pop_total": 10,
            "reg_tumor_total": 8,
            "conserved_regulators": conserved_regs or ["E2F1", "MYC"],
            "tumor_state_only_regulators": tumor_only or ["EGFR"],
            "population_averaged_only_regulators": pop_only or ["TP53"],
            "validated_conserved": validated_conserved or [],
            "tgt_conserved_count": 12,
            "tgt_conserved_fraction": 0.6,
            "tgt_tumor_only": ["MMP9"],
            "regnetagents_available": regnetagents_available,
            "cascade_available": cascade_available,
            "errors": errors or {},
        }

    def test_header_contains_gene(self, wf):
        syn = self._make_synthesis()
        report = "\n".join(wf._format_network_comparison_report(syn))
        assert "FOXM1" in report

    def test_header_contains_cancer_type(self, wf):
        syn = self._make_synthesis()
        report = "\n".join(wf._format_network_comparison_report(syn))
        assert "hnsc" in report.lower() or "HNSC" in report

    def test_rewiring_classification_shown(self, wf):
        syn = self._make_synthesis(rewiring="high")
        report = "\n".join(wf._format_network_comparison_report(syn))
        assert "HIGH" in report or "high" in report

    def test_conserved_cascade_validated_label(self, wf):
        vc = [{"gene": "E2F1", "tier": "conserved_cascade_validated",
               "cascade_key_findings": ["LINCS: E2F1 knockdown confirmed"], "cascade_error": None}]
        syn = self._make_synthesis(validated_conserved=vc)
        report = "\n".join(wf._format_network_comparison_report(syn))
        assert "E2F1" in report
        assert "CASCADE-validated" in report

    def test_conserved_not_validated_label(self, wf):
        vc = [{"gene": "MYC", "tier": "conserved_not_validated",
               "cascade_key_findings": [], "cascade_error": None}]
        syn = self._make_synthesis(validated_conserved=vc)
        report = "\n".join(wf._format_network_comparison_report(syn))
        assert "MYC" in report
        assert "no CASCADE experimental support" in report or "Conserved (no CASCADE" in report

    def test_tumor_only_regulators_shown(self, wf):
        syn = self._make_synthesis(tumor_only=["EGFR", "KRAS"])
        report = "\n".join(wf._format_network_comparison_report(syn))
        assert "EGFR" in report
        assert "KRAS" in report

    def test_regnetagents_unavailable_warning(self, wf):
        syn = self._make_synthesis(regnetagents_available=False,
                                   errors={"network_comparison": "timeout"})
        report = "\n".join(wf._format_network_comparison_report(syn))
        assert "RegNetAgents unavailable" in report

    def test_cascade_unavailable_warning(self, wf):
        syn = self._make_synthesis(cascade_available=False)
        report = "\n".join(wf._format_network_comparison_report(syn))
        assert "CASCADE unavailable" in report

    def test_target_overlap_shown(self, wf):
        syn = self._make_synthesis()
        report = "\n".join(wf._format_network_comparison_report(syn))
        assert "60.00%" in report or "Conserved targets" in report


# ---------------------------------------------------------------------------
# Graceful degradation
# ---------------------------------------------------------------------------

class TestNetworkComparisonGracefulDegradation:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    @pytest.mark.asyncio
    async def test_missing_cancer_type_sets_error(self, wf):
        state = _nc_state(cancer_type="")
        result = await wf._run_network_comparison_path(state)
        assert "network_comparison" in result.get("errors", {})

    def test_synthesize_with_no_network_comparison(self, wf):
        state = _nc_state(network_comparison=None, errors={"network_comparison": "failed"})
        result = wf._synthesize_network_comparison_path(state)
        syn = result["synthesis"]
        assert syn["routing"] == "network_comparison"
        assert syn["conserved_regulators"] == []
        assert syn["regnetagents_available"] is False


# ---------------------------------------------------------------------------
# Integration test (requires live servers)
# ---------------------------------------------------------------------------

INTEGRATION = pytest.mark.skipif(
    not os.environ.get("ORCHESTRA_INTEGRATION_TESTS"),
    reason="Set ORCHESTRA_INTEGRATION_TESTS=1 to run integration tests",
)


@INTEGRATION
class TestNetworkComparisonIntegration:
    """FOXM1 in epithelial_cell vs TCGA HNSC — end-to-end smoke test.

    HNSC (head/neck squamous) is the closest TCGA proxy for cervical squamous
    carcinoma: shared squamous histology and HPV etiology.
    """

    async def test_foxm1_hnsc_pipeline(self):
        from mcp_client import make_cascade_client, make_regnetagents_client
        from contextlib import AsyncExitStack

        wf = OrchestraWorkflow()
        async with AsyncExitStack() as stack:
            cascade = await stack.enter_async_context(make_cascade_client())
            regnetagents = await stack.enter_async_context(make_regnetagents_client())
            wf._persistent_cascade = cascade
            wf._persistent_regnetagents = regnetagents
            wf._persistent_ready.set()

            result = await wf.run_analysis(
                gene="FOXM1",
                cell_type="epithelial_cell",
                analysis_type="network_comparison",
                cancer_type="hnsc",
            )

        syn = result.get("synthesis") or {}

        assert syn.get("routing") == "network_comparison", (
            f"Wrong routing: {syn.get('routing')}"
        )
        assert "Network Context Comparison" in result.get("final_report", ""), (
            "Report header missing"
        )
        assert result.get("errors") == {}, f"Unexpected errors: {result.get('errors')}"

        # Structural checks on synthesis
        assert "conserved_regulators" in syn
        assert "tumor_state_only_regulators" in syn
        assert "rewiring_classification" in syn
        assert syn["rewiring_classification"] in ("low", "moderate", "high")
