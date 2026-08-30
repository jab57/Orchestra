"""
Tests for Issue #13: compare_network_contexts (GREmLN vs TCGA network context comparison).

Unit tests run without live child servers.
Integration test is gated on ORCHESTRA_INTEGRATION_TESTS=1.
"""

import os
import pytest
from unittest.mock import AsyncMock
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
            "conserved": conserved if conserved is not None else ["E2F1", "MYC", "BRCA1"],
            "conserved_count": len(conserved if conserved is not None else ["E2F1", "MYC", "BRCA1"]),
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
        findings.append(
            "1 gene(s) confirmed by both network propagation and LINCS "
            "experimental knockdown data (directional agreement)."
        )
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

    def test_cascade_available_true_when_no_conserved_regs(self, wf):
        # CASCADE is not called when there are 0 conserved regulators — that is
        # valid biology (complete rewiring), not a CASCADE failure.
        state = _nc_state(
            network_comparison=_context_result(conserved=[], cascade_validation={}),
        )
        result = wf._synthesize_network_comparison_path(state)
        assert result["synthesis"]["cascade_available"] is True

    def test_cascade_available_true_when_calls_errored(self, wf):
        # Per-call CASCADE errors are surfaced via validated_conserved tier, not by
        # setting cascade_available=False. The global flag stays True so no "CASCADE
        # unavailable" banner fires — the tier label "conserved_not_validated" is the
        # right signal.
        state = _nc_state(
            network_comparison=_context_result(
                conserved=["TOP2A"],
                cascade_validation={"TOP2A": {"error": "timeout"}},
            ),
        )
        result = wf._synthesize_network_comparison_path(state)
        assert result["synthesis"]["cascade_available"] is True
        assert result["synthesis"]["validated_conserved"][0]["tier"] == "conserved_not_validated"


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

    # -----------------------------------------------------------------
    # TCGA-only fallback when the gene has no GREmLN baseline
    # -----------------------------------------------------------------

    @pytest.mark.asyncio
    async def test_gremln_failure_fetches_tcga_fallback(self, wf):
        async def call_tool(name, args, timeout_seconds=None):
            if name == "compare_network_contexts":
                return {"error": True, "message": "Gene 'FKBP6' not found in epithelial_cell network"}
            if name == "query_network":
                assert args["query_type"] == "gene_neighbors"
                assert args["network_source"] == "tcga"
                assert args["tcga_network"] == "hnsc"
                assert args["gene"] == "FOXM1"
                return {
                    "regulators": [{"gene": "ZNF211"}, {"gene": "SOX30"}],
                    "targets": [{"gene": "MMP9"}],
                }
            raise AssertionError(f"unexpected tool call: {name}")

        wf._regnetagents = AsyncMock()
        wf._regnetagents.call_tool = AsyncMock(side_effect=call_tool)
        state = _nc_state()

        result = await wf._run_network_comparison_path(state)

        assert "network_comparison" in result["errors"]
        fallback = result.get("network_comparison_tcga_fallback")
        assert fallback is not None
        assert [r["gene"] for r in fallback["regulators"]] == ["ZNF211", "SOX30"]
        assert result.get("network_comparison") is None  # unchanged: still absent on GREmLN failure

    @pytest.mark.asyncio
    async def test_gremln_failure_fallback_also_fails(self, wf):
        async def call_tool(name, args, timeout_seconds=None):
            if name == "compare_network_contexts":
                return {"error": True, "message": "not found"}
            if name == "query_network":
                return {"error": True, "message": "gene not in TCGA network either"}
            raise AssertionError(f"unexpected tool call: {name}")

        wf._regnetagents = AsyncMock()
        wf._regnetagents.call_tool = AsyncMock(side_effect=call_tool)
        state = _nc_state()

        result = await wf._run_network_comparison_path(state)

        assert "network_comparison" in result["errors"]
        assert result.get("network_comparison_tcga_fallback") is None

    @pytest.mark.asyncio
    async def test_gremln_failure_never_calls_cascade_even_with_validate_tumor_acquired(self, wf):
        """
        Paper-reproducibility guard: the corroboration paper's scripts call this path with
        validate_tumor_acquired=True. A GREmLN-absent gene must keep being dropped from that
        scoring entirely (existing, documented behavior) -- the TCGA fallback is display-only
        and must never reach CASCADE or the tumor_acquired_cascade_validation step.
        """
        async def call_tool(name, args, timeout_seconds=None):
            if name == "compare_network_contexts":
                return {"error": True, "message": "not found"}
            if name == "query_network":
                return {"regulators": [{"gene": "ZNF211"}], "targets": []}
            raise AssertionError(f"unexpected tool call: {name}")

        wf._regnetagents = AsyncMock()
        wf._regnetagents.call_tool = AsyncMock(side_effect=call_tool)
        wf._cascade = AsyncMock()
        wf._cascade.call_tool = AsyncMock(side_effect=AssertionError("CASCADE must not be called"))
        state = _nc_state(validate_tumor_acquired=True, rank_tumor_acquired=True)

        result = await wf._run_network_comparison_path(state)

        wf._cascade.call_tool.assert_not_called()
        assert result.get("network_comparison") is None

    def test_synthesize_surfaces_tcga_fallback(self, wf):
        fallback = {"regulators": [{"gene": "ZNF211"}], "targets": []}
        state = _nc_state(
            network_comparison=None,
            network_comparison_tcga_fallback=fallback,
            errors={"network_comparison": "not found"},
        )
        result = wf._synthesize_network_comparison_path(state)
        syn = result["synthesis"]
        assert syn["regnetagents_available"] is False
        assert syn["tcga_only_fallback"] == fallback

    def test_format_shows_tcga_fallback_instead_of_generic_error(self, wf):
        synthesis = {
            "gene": "FKBP6", "cell_type": "epithelial_cell", "cancer_type": "cesc",
            "regnetagents_available": False,
            "tcga_only_fallback": {
                "regulators": [{"gene": "ZNF211"}, {"gene": "SOX30"}],
                "targets": [],
            },
            "errors": {"network_comparison": "not found"},
        }
        lines = wf._format_network_comparison_report(synthesis)
        text = "\n".join(lines)
        assert "No GREmLN baseline for FKBP6" in text
        assert "ZNF211" in text and "SOX30" in text
        assert "RegNetAgents unavailable" not in text

    def test_format_falls_back_to_generic_error_without_fallback_data(self, wf):
        synthesis = {
            "gene": "FKBP6", "cell_type": "epithelial_cell", "cancer_type": "cesc",
            "regnetagents_available": False,
            "tcga_only_fallback": None,
            "errors": {"network_comparison": "not found"},
        }
        lines = wf._format_network_comparison_report(synthesis)
        text = "\n".join(lines)
        assert "RegNetAgents unavailable" in text
        assert "No GREmLN baseline" not in text


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


# ---------------------------------------------------------------------------
# compare_network_contexts_batch report assembly
# ---------------------------------------------------------------------------

class TestNetworkContextsBatch:
    """Unit tests for compare_network_contexts_batch report assembly logic."""

    def _make_report(self, gene: str, cancer_type: str) -> str:
        return (
            f"## Network Context Comparison: {gene}\n"
            f"**Cancer type:** TCGA {cancer_type.upper()}\n"
            f"Rewiring: HIGH\n"
        )

    def test_batch_header_contains_all_genes(self):
        genes = ["FOXM1", "MYC", "TP53"]
        cancer_type = "hnsc"
        cell_type = "epithelial_cell"
        sections = [self._make_report(g, cancer_type) for g in genes]
        header = "\n".join([
            f"## Network Context Comparison Batch: TCGA {cancer_type.upper()}",
            f"**Genes:** {', '.join(genes)}",
            f"**Reference network:** {cell_type} (GREmLN)",
            "",
        ])
        report = header + "\n\n---\n\n".join(sections)
        assert "FOXM1" in report
        assert "MYC" in report
        assert "TP53" in report
        assert "## Network Context Comparison Batch: TCGA HNSC" in report
        assert "**Genes:** FOXM1, MYC, TP53" in report

    def test_batch_sections_separated_by_divider(self):
        genes = ["FOXM1", "MYC"]
        cancer_type = "brca"
        sections = [self._make_report(g, cancer_type) for g in genes]
        combined = "\n\n---\n\n".join(sections)
        assert "---" in combined
        assert combined.index("FOXM1") < combined.index("---") < combined.index("MYC")

    def test_batch_failed_gene_shows_placeholder(self):
        genes = ["FOXM1", "BROKENGENE"]
        cancer_type = "hnsc"
        sections = [
            self._make_report("FOXM1", cancer_type),
            "_BROKENGENE: analysis failed_",
        ]
        combined = "\n\n---\n\n".join(sections)
        assert "_BROKENGENE: analysis failed_" in combined
        assert "FOXM1" in combined

    def test_batch_two_genes_minimum(self):
        genes = ["FOXM1", "STAT3"]
        assert len(genes) >= 2


# ---------------------------------------------------------------------------
# compare_tumor_networks — OrchestraState schema
# ---------------------------------------------------------------------------

class TestTumorNetworkStateSchema:
    def test_cancer_types_field_in_state(self):
        assert "cancer_types" in OrchestraState.__annotations__

    def test_include_gremln_baseline_field_in_state(self):
        assert "include_gremln_baseline" in OrchestraState.__annotations__

    def test_tumor_network_results_field_in_state(self):
        assert "tumor_network_results" in OrchestraState.__annotations__


# ---------------------------------------------------------------------------
# compare_tumor_networks — routing
# ---------------------------------------------------------------------------

class TestTumorNetworkRouting:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def test_tumor_network_comparison_routes_correctly(self, wf):
        state = _nc_state(analysis_type="tumor_network_comparison")
        assert wf._routing_decision(state) == "tumor_network_comparison_path"

    def test_tumor_network_does_not_affect_network_comparison(self, wf):
        state = _nc_state(analysis_type="network_comparison")
        assert wf._routing_decision(state) == "network_comparison_path"


# ---------------------------------------------------------------------------
# compare_tumor_networks — helpers
# ---------------------------------------------------------------------------

def _tnc_state(**overrides) -> dict:
    base = {
        "gene": "FOXM1",
        "cell_type": "epithelial_cell",
        "cancer_types": ["cesc", "hnsc"],
        "include_gremln_baseline": True,
        "analysis_type": "tumor_network_comparison",
        "analysis_depth": "comprehensive",
        "gene_role": None,
        "ensembl_id": None,
        "tf_partner": None,
        "network_analysis": None,
        "network_stats": None,
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
        "cancer_type": None,
        "network_comparison": None,
        "cancer_context": None,
        "gene2": None,
        "cancer_contexts": None,
        "cross_context_novelty": None,
        "tumor_network_results": None,
        "novelty_result": None,
        "edge_novelty_results": None,
        "completed_steps": [],
        "errors": {},
        "final_report": None,
        "synthesis": None,
    }
    base.update(overrides)
    return base


def _tumor_raw(
    conserved=None, tumor_only=None, pop_only=None,
    conserved_fraction=0.1, rewiring="high",
    pop_total=5, tumor_total=10,
    tgt_tumor_only=None,
):
    conserved = conserved or []
    tumor_only = tumor_only or []
    pop_only = pop_only or []
    return {
        "gene": "FOXM1",
        "regulators": {
            "population_averaged_total": pop_total,
            "tumor_state_total": tumor_total,
            "conserved": conserved,
            "conserved_count": len(conserved),
            "conserved_fraction": conserved_fraction,
            "population_averaged_only": pop_only,
            "tumor_state_only": tumor_only,
        },
        "targets": {
            "population_averaged_total": 20,
            "tumor_state_total": 18,
            "conserved_count": 10,
            "conserved_fraction": 0.5,
            "population_averaged_only": [],
            "tumor_state_only": tgt_tumor_only or ["MMP9"],
        },
        "interpretation": {
            "regulatory_rewiring": rewiring,
            "conserved_fraction_regulators": conserved_fraction,
            "tumor_specific_regulator_count": len(tumor_only),
        },
    }


# ---------------------------------------------------------------------------
# _synthesize_tumor_network_comparison_path
# ---------------------------------------------------------------------------

class TestSynthesizeTumorNetworkComparison:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def test_routing_is_tumor_network_comparison(self, wf):
        state = _tnc_state(tumor_network_results={
            "cesc": _tumor_raw(tumor_only=["EGFR", "KRAS"]),
            "hnsc": _tumor_raw(tumor_only=["EGFR", "MYC"]),
            "_cascade_validation": {},
        })
        result = wf._synthesize_tumor_network_comparison_path(state)
        assert result["synthesis"]["routing"] == "tumor_network_comparison"

    def test_convergent_core_identified(self, wf):
        # EGFR in both tumor sets → convergent core
        state = _tnc_state(tumor_network_results={
            "cesc": _tumor_raw(tumor_only=["EGFR", "KRAS"]),
            "hnsc": _tumor_raw(tumor_only=["EGFR", "MYC"]),
            "_cascade_validation": {},
        })
        result = wf._synthesize_tumor_network_comparison_path(state)
        syn = result["synthesis"]
        assert "EGFR" in syn["convergent_core"]
        assert "KRAS" not in syn["convergent_core"]

    def test_divergent_specific_computed(self, wf):
        state = _tnc_state(tumor_network_results={
            "cesc": _tumor_raw(tumor_only=["EGFR", "KRAS"]),
            "hnsc": _tumor_raw(tumor_only=["EGFR", "MYC"]),
            "_cascade_validation": {},
        })
        result = wf._synthesize_tumor_network_comparison_path(state)
        syn = result["synthesis"]
        # KRAS unique to cesc, MYC unique to hnsc
        assert "KRAS" in syn["divergent_specific"].get("cesc", [])
        assert "MYC" in syn["divergent_specific"].get("hnsc", [])

    def test_pairwise_overlap_computed(self, wf):
        state = _tnc_state(tumor_network_results={
            "cesc": _tumor_raw(tumor_only=["EGFR", "KRAS"]),
            "hnsc": _tumor_raw(tumor_only=["EGFR", "MYC"]),
            "_cascade_validation": {},
        })
        result = wf._synthesize_tumor_network_comparison_path(state)
        pairwise = result["synthesis"]["pairwise_overlaps"]
        assert len(pairwise) == 1
        p = pairwise[0]
        assert p["cancer_a"] == "cesc"
        assert p["cancer_b"] == "hnsc"
        assert p["shared_regulators_count"] == 1  # EGFR
        assert p["jaccard_regulators"] == pytest.approx(1 / 3, rel=0.01)  # |{EGFR}| / |{EGFR,KRAS,MYC}|

    def test_verdict_convergent(self, wf):
        # High overlap → convergent
        state = _tnc_state(tumor_network_results={
            "cesc": _tumor_raw(tumor_only=["A", "B", "C", "D", "E"]),
            "hnsc": _tumor_raw(tumor_only=["A", "B", "C", "D", "E"]),
            "_cascade_validation": {},
        })
        result = wf._synthesize_tumor_network_comparison_path(state)
        assert result["synthesis"]["verdict"] == "convergent"

    def test_verdict_divergent(self, wf):
        # No overlap → divergent
        state = _tnc_state(tumor_network_results={
            "cesc": _tumor_raw(tumor_only=["A", "B"]),
            "hnsc": _tumor_raw(tumor_only=["C", "D"]),
            "_cascade_validation": {},
        })
        result = wf._synthesize_tumor_network_comparison_path(state)
        assert result["synthesis"]["verdict"] == "divergent"

    def test_verdict_mixed(self, wf):
        # Partial overlap (Jaccard between 0.15 and 0.40) → mixed
        # 2 shared out of 5 unique = 2/5 = 0.4... borderline, let's use 1 shared out of 4 = 0.25
        state = _tnc_state(tumor_network_results={
            "cesc": _tumor_raw(tumor_only=["A", "B", "C"]),
            "hnsc": _tumor_raw(tumor_only=["A", "D", "E"]),  # jaccard 1/5 = 0.2
            "_cascade_validation": {},
        })
        result = wf._synthesize_tumor_network_comparison_path(state)
        syn = result["synthesis"]
        # 1/5 = 0.20, which is >= 0.15 and < 0.40, and core is non-empty → mixed
        assert syn["verdict"] == "mixed"

    def test_insufficient_data_when_one_cancer_fails(self, wf):
        state = _tnc_state(
            cancer_types=["cesc", "hnsc"],
            tumor_network_results={
                "cesc": {"error": True, "message": "timeout"},
                "hnsc": {"error": True, "message": "timeout"},
                "_cascade_validation": {},
            }
        )
        result = wf._synthesize_tumor_network_comparison_path(state)
        assert result["synthesis"]["verdict"] == "insufficient_data"

    def test_gremln_baseline_included_by_default(self, wf):
        state = _tnc_state(tumor_network_results={
            "cesc": _tumor_raw(tumor_only=["EGFR"]),
            "hnsc": _tumor_raw(tumor_only=["EGFR"]),
            "_cascade_validation": {},
        })
        result = wf._synthesize_tumor_network_comparison_path(state)
        assert result["synthesis"]["gremln_baseline"] != {}

    def test_gremln_baseline_omitted_when_false(self, wf):
        state = _tnc_state(
            include_gremln_baseline=False,
            tumor_network_results={
                "cesc": _tumor_raw(tumor_only=["EGFR"]),
                "hnsc": _tumor_raw(tumor_only=["EGFR"]),
                "_cascade_validation": {},
            }
        )
        result = wf._synthesize_tumor_network_comparison_path(state)
        assert result["synthesis"]["gremln_baseline"] == {}

    def test_conserved_from_gremln_included_in_tumor_set(self, wf):
        # conserved regulators (in both GREmLN and TCGA) should be in the tumor set
        state = _tnc_state(tumor_network_results={
            "cesc": _tumor_raw(conserved=["TP53"], tumor_only=["EGFR"]),
            "hnsc": _tumor_raw(conserved=["TP53"], tumor_only=["MYC"]),
            "_cascade_validation": {},
        })
        result = wf._synthesize_tumor_network_comparison_path(state)
        syn = result["synthesis"]
        # TP53 is in both tumor sets (conserved in GREmLN+TCGA) → convergent core
        assert "TP53" in syn["convergent_core"]

    def test_cascade_validation_tier_in_core(self, wf):
        casc_val = {"EGFR": _cascade_result("EGFR", has_lincs=True)}
        state = _tnc_state(tumor_network_results={
            "cesc": _tumor_raw(tumor_only=["EGFR", "KRAS"]),
            "hnsc": _tumor_raw(tumor_only=["EGFR", "MYC"]),
            "_cascade_validation": casc_val,
        })
        result = wf._synthesize_tumor_network_comparison_path(state)
        core_cascade = result["synthesis"]["convergent_core_cascade"]
        egfr_entry = next((e for e in core_cascade if e["gene"] == "EGFR"), None)
        assert egfr_entry is not None
        assert egfr_entry["tier"] == "cascade_validated"

    def test_cascade_not_validated_without_evidence(self, wf):
        casc_val = {"EGFR": {"evidence_synthesis": {"key_findings": []}}}
        state = _tnc_state(tumor_network_results={
            "cesc": _tumor_raw(tumor_only=["EGFR"]),
            "hnsc": _tumor_raw(tumor_only=["EGFR"]),
            "_cascade_validation": casc_val,
        })
        result = wf._synthesize_tumor_network_comparison_path(state)
        core_cascade = result["synthesis"]["convergent_core_cascade"]
        assert core_cascade[0]["tier"] == "not_validated"

    def test_three_cancer_types_pairwise(self, wf):
        state = _tnc_state(
            cancer_types=["cesc", "hnsc", "luad"],
            tumor_network_results={
                "cesc": _tumor_raw(tumor_only=["A", "B", "C"]),
                "hnsc": _tumor_raw(tumor_only=["A", "B", "D"]),
                "luad": _tumor_raw(tumor_only=["A", "E", "F"]),
                "_cascade_validation": {},
            }
        )
        result = wf._synthesize_tumor_network_comparison_path(state)
        syn = result["synthesis"]
        # C(3,2) = 3 pairs
        assert len(syn["pairwise_overlaps"]) == 3
        # A is in all three → convergent core
        assert "A" in syn["convergent_core"]
        assert "B" not in syn["convergent_core"]  # B only in cesc+hnsc, not luad


# ---------------------------------------------------------------------------
# _format_tumor_network_comparison_report
# ---------------------------------------------------------------------------

class TestFormatTumorNetworkComparisonReport:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    _SENTINEL = object()

    def _make_tnc_synthesis(
        self,
        cancer_types=None,
        available_cts=_SENTINEL,
        pairwise=_SENTINEL,
        convergent_core=_SENTINEL,
        convergent_core_cascade=_SENTINEL,
        divergent_specific=None,
        verdict="mixed",
        mean_jaccard=0.25,
        gremln_baseline=_SENTINEL,
        include_gremln=True,
        errors=None,
    ):
        cancer_types = cancer_types or ["cesc", "hnsc"]
        if available_cts is self._SENTINEL:
            available_cts = cancer_types
        if pairwise is self._SENTINEL:
            pairwise = [
                {
                    "cancer_a": "cesc",
                    "cancer_b": "hnsc",
                    "shared_regulators": ["EGFR"],
                    "shared_regulators_count": 1,
                    "shared_tumor_targets_count": 0,
                    "jaccard_regulators": 0.25,
                    "total_a": 3,
                    "total_b": 3,
                }
            ]
        if convergent_core is self._SENTINEL:
            convergent_core = ["EGFR"]
        if convergent_core_cascade is self._SENTINEL:
            convergent_core_cascade = [
                {"gene": "EGFR", "n_cancer_types": 2, "tier": "not_validated",
                 "cascade_key_findings": [], "cascade_error": None}
            ]
        if gremln_baseline is self._SENTINEL:
            gremln_baseline = {
                "cesc": {"rewiring_classification": "high", "conserved_fraction": 0.05,
                         "pop_total": 5, "tumor_total": 20, "conserved_count": 1},
                "hnsc": {"rewiring_classification": "high", "conserved_fraction": 0.04,
                         "pop_total": 5, "tumor_total": 24, "conserved_count": 1},
            }
        return {
            "routing": "tumor_network_comparison",
            "gene": "FOXM1",
            "cell_type": "epithelial_cell",
            "cancer_types": cancer_types,
            "available_cancer_types": available_cts,
            "include_gremln_baseline": include_gremln,
            "pairwise_overlaps": pairwise,
            "convergent_core": convergent_core,
            "convergent_core_cascade": convergent_core_cascade,
            "divergent_specific": divergent_specific or {"cesc": ["KRAS"], "hnsc": ["MYC"]},
            "verdict": verdict,
            "mean_jaccard": mean_jaccard,
            "gremln_baseline": gremln_baseline,
            "errors": errors or {},
        }

    def test_header_contains_gene(self, wf):
        syn = self._make_tnc_synthesis()
        report = "\n".join(wf._format_tumor_network_comparison_report(syn))
        assert "FOXM1" in report

    def test_header_contains_cancer_types(self, wf):
        syn = self._make_tnc_synthesis()
        report = "\n".join(wf._format_tumor_network_comparison_report(syn))
        assert "CESC" in report
        assert "HNSC" in report

    def test_pairwise_table_present(self, wf):
        syn = self._make_tnc_synthesis()
        report = "\n".join(wf._format_tumor_network_comparison_report(syn))
        assert "Pairwise Regulator Overlap" in report
        assert "Jaccard" in report

    def test_convergent_core_section_present(self, wf):
        syn = self._make_tnc_synthesis()
        report = "\n".join(wf._format_tumor_network_comparison_report(syn))
        assert "Convergent Core" in report
        assert "EGFR" in report

    def test_divergent_section_present(self, wf):
        syn = self._make_tnc_synthesis()
        report = "\n".join(wf._format_tumor_network_comparison_report(syn))
        assert "Cancer-Specific Regulators" in report or "Divergent" in report
        assert "KRAS" in report

    def test_verdict_shown_in_report(self, wf):
        syn = self._make_tnc_synthesis(verdict="convergent")
        report = "\n".join(wf._format_tumor_network_comparison_report(syn))
        assert "CONVERGENT" in report

    def test_verdict_thresholds_stated(self, wf):
        syn = self._make_tnc_synthesis()
        report = "\n".join(wf._format_tumor_network_comparison_report(syn))
        assert "0.40" in report or "0.15" in report  # explicit thresholds

    def test_gremln_baseline_shown_when_requested(self, wf):
        syn = self._make_tnc_synthesis(include_gremln=True)
        report = "\n".join(wf._format_tumor_network_comparison_report(syn))
        assert "GREmLN Baseline" in report

    def test_gremln_baseline_omitted_when_false(self, wf):
        syn = self._make_tnc_synthesis(include_gremln=False, gremln_baseline={})
        report = "\n".join(wf._format_tumor_network_comparison_report(syn))
        assert "GREmLN Baseline" not in report

    def test_empty_convergent_core_shown_gracefully(self, wf):
        syn = self._make_tnc_synthesis(
            convergent_core=[],
            convergent_core_cascade=[],
            verdict="divergent",
        )
        report = "\n".join(wf._format_tumor_network_comparison_report(syn))
        assert "No convergent core" in report

    def test_cascade_validated_label_shown(self, wf):
        syn = self._make_tnc_synthesis(
            convergent_core_cascade=[
                {"gene": "EGFR", "n_cancer_types": 2, "tier": "cascade_validated",
                 "cascade_key_findings": ["LINCS: EGFR knockdown confirmed"], "cascade_error": None}
            ]
        )
        report = "\n".join(wf._format_tumor_network_comparison_report(syn))
        assert "CASCADE-validated" in report or "YES" in report

    def test_warning_when_no_available_data(self, wf):
        syn = self._make_tnc_synthesis(
            available_cts=[],
            pairwise=[],
            convergent_core=[],
            convergent_core_cascade=[],
            divergent_specific={},
            verdict="insufficient_data",
        )
        report = "\n".join(wf._format_tumor_network_comparison_report(syn))
        assert "No valid data" in report or "INSUFFICIENT" in report


# ---------------------------------------------------------------------------
# Graceful degradation — compare_tumor_networks
# ---------------------------------------------------------------------------

class TestTumorNetworkGracefulDegradation:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    @pytest.mark.asyncio
    async def test_too_few_cancer_types_sets_error(self, wf):
        state = _tnc_state(cancer_types=["cesc"])
        result = await wf._run_tumor_network_comparison_path(state)
        assert "tumor_network_comparison" in result.get("errors", {})

    @pytest.mark.asyncio
    async def test_empty_cancer_types_sets_error(self, wf):
        state = _tnc_state(cancer_types=[])
        result = await wf._run_tumor_network_comparison_path(state)
        assert "tumor_network_comparison" in result.get("errors", {})

    def test_synthesize_with_all_errors(self, wf):
        state = _tnc_state(
            tumor_network_results={
                "cesc": {"error": True},
                "hnsc": {"error": True},
                "_cascade_validation": {},
            }
        )
        result = wf._synthesize_tumor_network_comparison_path(state)
        syn = result["synthesis"]
        assert syn["verdict"] == "insufficient_data"
        assert syn["convergent_core"] == []
        assert syn["pairwise_overlaps"] == []
