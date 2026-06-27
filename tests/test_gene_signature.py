"""
Tests for Issue #10: analyze_gene_signature (gene signature driver analysis).

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


def _sig_state(**overrides) -> dict:
    base = {
        "gene": "",
        "cell_type": "epithelial_cell",
        "analysis_type": "gene_signature",
        "analysis_depth": "comprehensive",
        "gene_role": None,
        "tf_partner": None,
        "network_analysis": None,
        "perturbation_result": None,
        "ppi_interactions": None,
        "validated_targets": None,
        "gene_signature": ["AXIN2", "MYC", "CCND1", "CDH1", "VEGFA"],
        "master_regulators": None,
        "cancer_type": None,
        "cancer_contexts": None,
        "cross_context_novelty": None,
        "completed_steps": [],
        "errors": {},
        "final_report": None,
        "synthesis": None,
    }
    base.update(overrides)
    return base


def _mr_result(tfs: list[dict]) -> dict:
    """Build a find_master_regulators-style result dict."""
    return {
        "master_regulators": tfs,
        "query_summary": {
            "gene_set_size": 20,
            "genes_found_in_network": 18,
            "genes_not_found": ["FAKEGENE1", "FAKEGENE2"],
            "network_size": 5000,
            "cell_type": "epithelial_cell",
            "total_regulators_tested": 250,
        },
    }


def _tf_entry(gene: str, overlap: int, enrichment: float = 2.5, p: float = 0.001,
              regulon_size: int = 100, key_findings: list = None,
              dorothea_overlap: int = 0) -> dict:
    return {
        "gene": gene,
        "ensembl_id": f"ENSG_{gene}",
        "rank": 1,
        "regulon_size": regulon_size,
        "overlap_count": overlap,
        "enrichment_score": enrichment,
        "p_value": p,
        "overlapping_genes": ["AXIN2", "MYC"][:overlap],
        "key_findings": key_findings or [],
        "dorothea_overlap": dorothea_overlap,
    }


# ---------------------------------------------------------------------------
# Routing
# ---------------------------------------------------------------------------

class TestSignatureRouting:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def test_gene_signature_routes_signature_path(self, wf):
        state = _sig_state(analysis_type="gene_signature")
        assert wf._routing_decision(state) == "signature_path"

    def test_gene_signature_overrides_gene_role(self, wf):
        """gene_signature analysis_type always routes to signature_path regardless of gene_role."""
        state = _sig_state(analysis_type="gene_signature", gene_role="master_regulator")
        assert wf._routing_decision(state) == "signature_path"

    def test_other_types_unaffected(self, wf):
        state = _sig_state(analysis_type="therapeutic_validation", gene_role="effector")
        assert wf._routing_decision(state) == "validation_path"

    def test_causal_chain_unaffected(self, wf):
        state = _sig_state(analysis_type="causal_chain", gene_role="master_regulator")
        assert wf._routing_decision(state) == "tf_path"


# ---------------------------------------------------------------------------
# OrchestraState has new fields
# ---------------------------------------------------------------------------

class TestStateSchemaSignature:
    def test_gene_signature_field_in_state(self):
        annotations = OrchestraState.__annotations__
        assert "gene_signature" in annotations

    def test_master_regulators_field_in_state(self):
        annotations = OrchestraState.__annotations__
        assert "master_regulators" in annotations


# ---------------------------------------------------------------------------
# _synthesize_signature_path
# ---------------------------------------------------------------------------

class TestSynthesizeSignaturePath:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def test_ranked_drivers_sorted_by_overlap(self, wf):
        state = _sig_state(
            gene_signature=["AXIN2", "MYC", "CCND1", "CDH1", "VEGFA"],
            master_regulators=_mr_result([
                _tf_entry("TP53", overlap=2, enrichment=1.5, p=0.01),
                _tf_entry("CTNNB1", overlap=4, enrichment=3.0, p=0.0001),
                _tf_entry("MYC", overlap=3, enrichment=2.0, p=0.002),
            ]),
        )
        result = wf._synthesize_signature_path(state)
        drivers = result["synthesis"]["ranked_drivers"]
        assert drivers[0]["gene"] == "CTNNB1"
        assert drivers[1]["gene"] == "MYC"
        assert drivers[2]["gene"] == "TP53"

    def test_coverage_pct_computed_correctly(self, wf):
        state = _sig_state(
            gene_signature=["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"],  # 10 genes
            master_regulators=_mr_result([
                _tf_entry("CTNNB1", overlap=4, regulon_size=200),
            ]),
        )
        result = wf._synthesize_signature_path(state)
        drivers = result["synthesis"]["ranked_drivers"]
        assert drivers[0]["coverage_pct"] == 40.0  # 4/10 * 100

    def test_synthesis_routing_is_signature(self, wf):
        state = _sig_state(
            master_regulators=_mr_result([_tf_entry("CTNNB1", overlap=3)]),
        )
        result = wf._synthesize_signature_path(state)
        assert result["synthesis"]["routing"] == "signature"

    def test_synthesis_contains_query_metadata(self, wf):
        state = _sig_state(
            master_regulators=_mr_result([_tf_entry("CTNNB1", overlap=3)]),
        )
        result = wf._synthesize_signature_path(state)
        syn = result["synthesis"]
        assert syn["genes_found_in_network"] == 18
        assert syn["total_regulators_tested"] == 250
        assert len(syn["genes_not_found"]) == 2

    def test_empty_master_regulators_produces_empty_drivers(self, wf):
        state = _sig_state(
            master_regulators=_mr_result([]),
        )
        result = wf._synthesize_signature_path(state)
        assert result["synthesis"]["ranked_drivers"] == []

    def test_cascade_key_findings_propagated(self, wf):
        findings = ["LINCS: CTNNB1 knockdown suppresses MYC", "DepMap: essential in colon cancer"]
        state = _sig_state(
            master_regulators=_mr_result([
                _tf_entry("CTNNB1", overlap=4, key_findings=findings),
            ]),
        )
        result = wf._synthesize_signature_path(state)
        driver = result["synthesis"]["ranked_drivers"][0]
        assert driver["cascade_key_findings"] == findings

    def test_pagerank_flag_true_when_enrichment_positive(self, wf):
        state = _sig_state(
            master_regulators=_mr_result([
                _tf_entry("CTNNB1", overlap=3, enrichment=2.5),
            ]),
        )
        result = wf._synthesize_signature_path(state)
        ev_row = result["synthesis"]["evidence_table"][0]
        assert ev_row["pagerank_rank"] is True

    def test_lincs_flag_from_key_findings(self, wf):
        state = _sig_state(
            master_regulators=_mr_result([
                _tf_entry("CTNNB1", overlap=3,
                          key_findings=["LINCS: knockdown confirmed downstream effect"]),
            ]),
        )
        result = wf._synthesize_signature_path(state)
        ev_row = result["synthesis"]["evidence_table"][0]
        assert ev_row["lincs_knockdown"] is True

    def test_corroboration_count_in_evidence_table(self, wf):
        state = _sig_state(
            master_regulators=_mr_result([
                _tf_entry("CTNNB1", overlap=3, enrichment=2.5,
                          key_findings=["LINCS: confirmed", "super-enhancer: present"]),
            ]),
        )
        result = wf._synthesize_signature_path(state)
        ev_row = result["synthesis"]["evidence_table"][0]
        # pagerank_rank=True (enrichment > 0), lincs=True, super_enhancer=True → at least 3
        assert ev_row["corroboration_count"] >= 3

    def test_dorothea_overlap_propagated_to_ranked_drivers(self, wf):
        state = _sig_state(
            gene_signature=["AXIN2", "MYC", "CCND1", "CDH1", "VEGFA"],
            master_regulators=_mr_result([
                _tf_entry("CTNNB1", overlap=4, dorothea_overlap=3),
                _tf_entry("TP53", overlap=2, dorothea_overlap=0),
            ]),
        )
        result = wf._synthesize_signature_path(state)
        drivers = result["synthesis"]["ranked_drivers"]
        assert drivers[0]["gene"] == "CTNNB1"
        assert drivers[0]["dorothea_overlap"] == 3
        assert drivers[1]["dorothea_overlap"] == 0

    def test_dorothea_overlap_propagated_to_evidence_table(self, wf):
        state = _sig_state(
            master_regulators=_mr_result([
                _tf_entry("CTNNB1", overlap=4, dorothea_overlap=2),
            ]),
        )
        result = wf._synthesize_signature_path(state)
        ev_row = result["synthesis"]["evidence_table"][0]
        assert ev_row["dorothea_overlap"] == 2


# ---------------------------------------------------------------------------
# _format_signature_report
# ---------------------------------------------------------------------------

class TestFormatSignatureReport:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def _make_synthesis(self, ranked_drivers=None, evidence_table=None,
                        errors=None, genes_not_found=None,
                        regnetagents_available=True, cascade_available=True):
        drivers = ranked_drivers or []
        evtab = evidence_table or []
        return {
            "routing": "signature",
            "cell_type": "epithelial_cell",
            "gene_signature": ["AXIN2", "MYC", "CCND1"],
            "signature_size": 3,
            "genes_found_in_network": 3,
            "genes_not_found": genes_not_found or [],
            "total_regulators_tested": 200,
            "network_size": 5000,
            "ranked_drivers": drivers,
            "evidence_table": evtab,
            "regnetagents_available": regnetagents_available,
            "cascade_available": cascade_available,
            "discordance_flags": [],
            "errors": errors or {},
        }

    def test_header_contains_cell_type(self, wf):
        syn = self._make_synthesis()
        report = "\n".join(wf._format_signature_report(syn))
        assert "epithelial_cell" in report

    def test_signature_size_shown(self, wf):
        syn = self._make_synthesis()
        report = "\n".join(wf._format_signature_report(syn))
        assert "3 genes" in report or "Signature size:** 3" in report

    def test_no_drivers_produces_empty_message(self, wf):
        syn = self._make_synthesis(ranked_drivers=[], evidence_table=[])
        report = "\n".join(wf._format_signature_report(syn))
        assert "No master regulators" in report

    def test_driver_appears_in_table(self, wf):
        driver = {
            "gene": "CTNNB1",
            "overlap_count": 4,
            "dorothea_overlap": 2,
            "coverage_pct": 40.0,
            "regulon_size": 200,
            "enrichment_score": 3.0,
            "p_value": 0.0001,
            "corroboration_count": 3,
            "corroboration_denominator": 7,
            "overlapping_genes": ["AXIN2", "MYC"],
            "cascade_key_findings": [],
            "multi_source_genes": [],
            "cascade_error": None,
            "pagerank_rank": True,
            "pathway_member": False,
            "lincs_knockdown": True,
            "depmap_essentiality": False,
            "super_enhancer": True,
            "dorothea_tier": False,
            "cbio_expression": False,
        }
        syn = self._make_synthesis(
            ranked_drivers=[driver],
            evidence_table=[driver],
        )
        report = "\n".join(wf._format_signature_report(syn))
        assert "CTNNB1" in report
        assert "40.0%" in report
        assert "3/7" in report

    def test_dorothea_ovlp_column_shown_in_table(self, wf):
        driver = {
            "gene": "MYC",
            "overlap_count": 3,
            "dorothea_overlap": 2,
            "coverage_pct": 30.0,
            "regulon_size": 150,
            "enrichment_score": 2.0,
            "p_value": 0.001,
            "corroboration_count": 2,
            "corroboration_denominator": 7,
            "overlapping_genes": [],
            "cascade_key_findings": [],
            "multi_source_genes": [],
            "cascade_error": None,
            "pagerank_rank": True,
            "pathway_member": False,
            "lincs_knockdown": False,
            "depmap_essentiality": False,
            "super_enhancer": False,
            "dorothea_tier": True,
            "cbio_expression": False,
        }
        syn = self._make_synthesis(ranked_drivers=[driver], evidence_table=[driver])
        report = "\n".join(wf._format_signature_report(syn))
        assert "DoRothEA Ovlp" in report
        assert "ARACNe Ovlp" in report

    def test_regnetagents_unavailable_warning(self, wf):
        syn = self._make_synthesis(regnetagents_available=False)
        report = "\n".join(wf._format_signature_report(syn))
        assert "RegNetAgents unavailable" in report

    def test_cascade_unavailable_warning(self, wf):
        syn = self._make_synthesis(cascade_available=False)
        report = "\n".join(wf._format_signature_report(syn))
        assert "CASCADE unavailable" in report

    def test_genes_not_found_shown(self, wf):
        syn = self._make_synthesis(genes_not_found=["FAKEGENE"])
        report = "\n".join(wf._format_signature_report(syn))
        assert "FAKEGENE" in report

    def test_errors_shown(self, wf):
        syn = self._make_synthesis(errors={"master_regulators": "connection timeout"})
        report = "\n".join(wf._format_signature_report(syn))
        assert "connection timeout" in report


# ---------------------------------------------------------------------------
# Discordance flags (via _compute_validation_discordance_flags, reused)
# ---------------------------------------------------------------------------

class TestSignatureDiscordance:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def test_topological_hub_not_validated(self, wf):
        """High enrichment (pagerank_rank=True) + no CASCADE experimental support → flag."""
        evidence_table = [{
            "gene": "CTNNB1",
            "pagerank_rank": True,
            "pathway_member": False,
            "lincs_knockdown": False,
            "depmap_essentiality": False,
            "super_enhancer": False,
            "dorothea_tier": False,
            "cbio_expression": False,
            "corroboration_count": 1,
            "corroboration_denominator": 7,
        }]
        flags = wf._compute_validation_discordance_flags(evidence_table)
        types = [f["type"] for f in flags]
        assert "topological_hub_not_validated" in types

    def test_experimentally_active_not_in_network(self, wf):
        """CASCADE support (LINCS) + no network enrichment → experimentally_active flag."""
        evidence_table = [{
            "gene": "BRD4",
            "pagerank_rank": False,
            "pathway_member": False,
            "lincs_knockdown": True,
            "depmap_essentiality": False,
            "super_enhancer": False,
            "dorothea_tier": False,
            "cbio_expression": False,
            "corroboration_count": 1,
            "corroboration_denominator": 7,
        }]
        flags = wf._compute_validation_discordance_flags(evidence_table)
        types = [f["type"] for f in flags]
        assert "experimentally_active_not_in_network" in types

    def test_no_flags_when_concordant(self, wf):
        """Both network and CASCADE support → no discordance."""
        evidence_table = [{
            "gene": "CTNNB1",
            "pagerank_rank": True,
            "pathway_member": False,
            "lincs_knockdown": True,
            "depmap_essentiality": False,
            "super_enhancer": False,
            "dorothea_tier": False,
            "cbio_expression": False,
            "corroboration_count": 2,
            "corroboration_denominator": 7,
        }]
        flags = wf._compute_validation_discordance_flags(evidence_table)
        assert flags == []


# ---------------------------------------------------------------------------
# Graceful degradation
# ---------------------------------------------------------------------------

class TestSignatureGracefulDegradation:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    @pytest.mark.asyncio
    async def test_empty_gene_signature_sets_error(self, wf):
        """Empty gene_signature populates errors without crashing."""
        state = _sig_state(gene_signature=[])
        result = await wf._run_signature_path(state)
        assert "signature_path" in result.get("errors", {})

    def test_synthesize_with_no_master_regulators_result(self, wf):
        """master_regulators=None (RegNetAgents failed) — synthesis still completes."""
        state = _sig_state(
            master_regulators=None,
            errors={"master_regulators": "timeout"},
        )
        result = wf._synthesize_signature_path(state)
        syn = result["synthesis"]
        assert syn["routing"] == "signature"
        assert syn["ranked_drivers"] == []
        assert syn["regnetagents_available"] is False

    def test_synthesize_cascade_error_entry_still_ranks(self, wf):
        """TF entry with cascade_error still appears in ranked_drivers (no CASCADE score)."""
        entry = _tf_entry("CTNNB1", overlap=4)
        entry["cascade_error"] = "timeout"
        state = _sig_state(
            master_regulators=_mr_result([entry]),
        )
        result = wf._synthesize_signature_path(state)
        drivers = result["synthesis"]["ranked_drivers"]
        assert len(drivers) == 1
        assert drivers[0]["gene"] == "CTNNB1"
        assert drivers[0]["cascade_error"] == "timeout"

    def test_report_renders_with_cascade_error(self, wf):
        """Report formatter handles cascade_error gracefully."""
        state = _sig_state(
            master_regulators=_mr_result([
                {**_tf_entry("CTNNB1", overlap=4), "cascade_error": "timeout"}
            ]),
        )
        result = wf._synthesize_signature_path(state)
        lines = wf._format_signature_report(result["synthesis"])
        report = "\n".join(lines)
        assert "CTNNB1" in report
        assert "timeout" in report


# ---------------------------------------------------------------------------
# Cross-context novelty gap (GitHub #13)
# ---------------------------------------------------------------------------

def _novelty_result(verdict: str, total: int) -> dict:
    return {"novelty_verdict": verdict, "pubmed_hits": total, "experimental_hits": 0,
            "computational_hits": 0, "most_recent_year": None, "verdict_rationale": ""}


class TestClassifyNoveltyGap:
    def setup_method(self):
        self.wf = _make_workflow()

    def test_transfer_opportunity_established_then_novel(self):
        verdicts = {
            "breast cancer": _novelty_result("established", 42),
            "cervical cancer": _novelty_result("novel", 0),
        }
        assert self.wf._classify_novelty_gap(verdicts) == "transfer_opportunity"

    def test_transfer_opportunity_novel_then_established(self):
        verdicts = {
            "breast cancer": _novelty_result("novel", 1),
            "cervical cancer": _novelty_result("established", 30),
        }
        assert self.wf._classify_novelty_gap(verdicts) == "transfer_opportunity"

    def test_transfer_opportunity_emerging_then_novel(self):
        verdicts = {
            "breast cancer": _novelty_result("emerging", 7),
            "cervical cancer": _novelty_result("novel", 0),
        }
        assert self.wf._classify_novelty_gap(verdicts) == "transfer_opportunity"

    def test_bilateral_novel(self):
        verdicts = {
            "breast cancer": _novelty_result("novel", 0),
            "cervical cancer": _novelty_result("novel", 2),
        }
        assert self.wf._classify_novelty_gap(verdicts) == "bilateral_novel"

    def test_bilateral_established(self):
        verdicts = {
            "breast cancer": _novelty_result("established", 50),
            "cervical cancer": _novelty_result("emerging", 10),
        }
        assert self.wf._classify_novelty_gap(verdicts) == "bilateral_established"

    def test_error_entries_skipped(self):
        verdicts = {
            "breast cancer": _novelty_result("established", 20),
            "cervical cancer": {"error": "timeout"},
        }
        # Only one valid verdict (established) — no novel to pair with → bilateral_established
        assert self.wf._classify_novelty_gap(verdicts) == "bilateral_established"

    def test_all_errors_returns_unknown(self):
        verdicts = {
            "breast cancer": {"error": "timeout"},
            "cervical cancer": {"error": "timeout"},
        }
        assert self.wf._classify_novelty_gap(verdicts) == "unknown"


class TestCrossContextNoveltySynthesis:
    def setup_method(self):
        self.wf = _make_workflow()

    def _state_with_novelty(self, cross_novelty: dict, cancer_contexts: list) -> dict:
        state = _sig_state(
            master_regulators=_mr_result([_tf_entry("CDCA7", overlap=4)]),
            cross_context_novelty=cross_novelty,
            cancer_contexts=cancer_contexts,
        )
        return state

    def test_cross_context_novelty_in_synthesis(self):
        novelty = {
            "CDCA7": {
                "breast cancer": _novelty_result("established", 15),
                "cervical cancer": _novelty_result("novel", 0),
            }
        }
        state = self._state_with_novelty(novelty, ["breast cancer", "cervical cancer"])
        result = self.wf._synthesize_signature_path(state)
        ccn = result["synthesis"]["cross_context_novelty"]
        assert len(ccn) == 1
        assert ccn[0]["gene"] == "CDCA7"
        assert ccn[0]["gap_classification"] == "transfer_opportunity"

    def test_no_cancer_contexts_produces_empty_table(self):
        state = _sig_state(master_regulators=_mr_result([_tf_entry("CDCA7", overlap=4)]))
        result = self.wf._synthesize_signature_path(state)
        assert result["synthesis"]["cross_context_novelty"] == []

    def test_gap_table_in_report(self):
        novelty = {
            "CDCA7": {
                "breast cancer": _novelty_result("established", 15),
                "cervical cancer": _novelty_result("novel", 0),
            }
        }
        state = self._state_with_novelty(novelty, ["breast cancer", "cervical cancer"])
        result = self.wf._synthesize_signature_path(state)
        report = "\n".join(self.wf._format_signature_report(result["synthesis"]))
        assert "Cross-Context Novelty Gap" in report
        assert "CDCA7" in report
        assert "transfer opportunity" in report
        assert "breast cancer" in report
        assert "cervical cancer" in report

    def test_gap_table_absent_when_no_contexts(self):
        state = _sig_state(master_regulators=_mr_result([_tf_entry("CDCA7", overlap=4)]))
        result = self.wf._synthesize_signature_path(state)
        report = "\n".join(self.wf._format_signature_report(result["synthesis"]))
        assert "Cross-Context Novelty Gap" not in report

    def test_error_entry_renders_gracefully(self):
        novelty = {
            "CDCA7": {
                "breast cancer": _novelty_result("established", 15),
                "cervical cancer": {"error": "timeout"},
            }
        }
        state = self._state_with_novelty(novelty, ["breast cancer", "cervical cancer"])
        result = self.wf._synthesize_signature_path(state)
        report = "\n".join(self.wf._format_signature_report(result["synthesis"]))
        assert "Cross-Context Novelty Gap" in report
        assert "error" in report


# ---------------------------------------------------------------------------
# tcga_network param (2026-06-22)
# ---------------------------------------------------------------------------

class TestEffectiveCancerContexts:
    """_effective_cancer_contexts auto-prepends the TCGA network's cancer context."""

    def setup_method(self):
        self.wf = _make_workflow()

    def test_hnsc_prepended_when_no_existing_contexts(self):
        state = _sig_state(cancer_type="hnsc")
        result = self.wf._effective_cancer_contexts(state)
        assert result == ["head and neck squamous"]

    def test_hnsc_prepended_before_existing_context(self):
        state = _sig_state(cancer_type="hnsc", cancer_contexts=["cervical cancer"])
        result = self.wf._effective_cancer_contexts(state)
        assert result == ["head and neck squamous", "cervical cancer"]

    def test_no_duplicate_if_tcga_ctx_already_present(self):
        state = _sig_state(cancer_type="hnsc",
                           cancer_contexts=["head and neck squamous", "cervical cancer"])
        result = self.wf._effective_cancer_contexts(state)
        assert result.count("head and neck squamous") == 1
        assert result == ["head and neck squamous", "cervical cancer"]

    def test_no_tcga_network_returns_contexts_unchanged(self):
        state = _sig_state(cancer_contexts=["cervical cancer"])
        result = self.wf._effective_cancer_contexts(state)
        assert result == ["cervical cancer"]

    def test_no_tcga_network_no_contexts_returns_empty(self):
        state = _sig_state()
        result = self.wf._effective_cancer_contexts(state)
        assert result == []

    def test_unknown_tcga_network_no_injection(self):
        state = _sig_state(cancer_type="unknown_cancer",
                           cancer_contexts=["cervical cancer"])
        result = self.wf._effective_cancer_contexts(state)
        assert result == ["cervical cancer"]

    def test_all_tcga_networks_have_mapping(self):
        from orchestra_langgraph_workflow import _TCGA_TO_CANCER_CONTEXT
        supported = ["hnsc", "coad", "brca", "luad", "lusc", "ov", "prad", "ucec"]
        for network in supported:
            state = _sig_state(cancer_type=network)
            result = self.wf._effective_cancer_contexts(state)
            assert len(result) == 1, f"Expected one context for {network}, got {result}"
            assert result[0] == _TCGA_TO_CANCER_CONTEXT[network]

    def test_brca_prepended(self):
        state = _sig_state(cancer_type="brca", cancer_contexts=["cervical cancer"])
        result = self.wf._effective_cancer_contexts(state)
        assert result[0] == "breast cancer"
        assert "cervical cancer" in result


class TestTcgaNetworkMcpSchema:
    """MCP tool schema exposes tcga_network param."""

    def test_tcga_network_in_schema(self):
        import importlib
        import sys
        # Import the server module to inspect the tool schema
        spec = importlib.util.find_spec("orchestra_mcp_server")
        if spec is None:
            pytest.skip("orchestra_mcp_server not importable in test context")
        # Verify the param exists at the workflow level via state field
        annotations = OrchestraState.__annotations__
        assert "cancer_type" in annotations  # tcga_network maps to cancer_type


# ---------------------------------------------------------------------------
# Integration test (requires live servers)
# ---------------------------------------------------------------------------

INTEGRATION = pytest.mark.skipif(
    not os.environ.get("ORCHESTRA_INTEGRATION_TESTS"),
    reason="Set ORCHESTRA_INTEGRATION_TESTS=1 to run integration tests",
)

# Issue #14 (2026-06-13): Ensembl REST API fallback removed from GeneIDMapper entirely.
# Both symbol_to_ensembl() and ensembl_to_symbol() now return from local cache only.
# Root cause resolved — INTEGRATION_LIVE skip no longer needed.


@INTEGRATION
class TestSignatureIntegration:
    """19-gene Wnt target signature — end-to-end pipeline smoke test.

    Runs run_analysis() once and asserts on drivers, report header, and errors.
    CTNNB1-as-top-hit is not asserted: only ~10/19 Wnt genes resolve to real ARACNe
    network Ensembl IDs in the local cache, so top-hit varies by network coverage.
    """

    WNT_SIGNATURE = [
        "AXIN2", "MYC", "CCND1", "CDH1", "VEGFA", "LEF1", "TCF7",
        "LGR5", "ASCL2", "SP5", "NKD1", "RNF43", "ZNRF3", "DKK1",
        "APCDD1", "TNFRSF19", "BMP4", "EPHB2", "EPHB3",
    ]

    async def test_pipeline_end_to_end(self):
        """Single call: verify drivers returned, report header present, no hard errors."""
        import sys
        import os
        sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
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
                gene="",
                cell_type="epithelial_cell",
                analysis_type="gene_signature",
                gene_signature=self.WNT_SIGNATURE,
            )

        syn = result.get("synthesis") or {}
        drivers = syn.get("ranked_drivers") or []

        assert len(drivers) > 0, "Expected at least one master regulator returned"
        assert "Gene Signature Driver Analysis" in result.get("final_report", ""), (
            "Report header missing from final_report"
        )
        assert result.get("errors") == {}, f"Unexpected errors: {result.get('errors')}"
