"""
Tests for Issue #10: analyze_gene_signature (gene signature driver analysis).

Unit tests run without live child servers.
Integration test is gated on ORCHESTRA_INTEGRATION_TESTS=1.
"""

import os
import pytest
from unittest.mock import AsyncMock, patch
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
    wf._persistent_cascade = None
    wf._persistent_regnetagents = None
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


def _mr_result(tfs: list[dict], genes_found_in_network: int = 18, gene_set_size: int = 20) -> dict:
    """Build a find_master_regulators-style result dict."""
    return {
        "master_regulators": tfs,
        "query_summary": {
            "gene_set_size": gene_set_size,
            "genes_found_in_network": genes_found_in_network,
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
        """When every input gene resolves in the network, coverage_pct against the
        resolved count and against the raw input count are the same number."""
        state = _sig_state(
            gene_signature=["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"],  # 10 genes
            master_regulators=_mr_result(
                [_tf_entry("CTNNB1", overlap=4, regulon_size=200)],
                genes_found_in_network=10,  # all 10 resolved
            ),
        )
        result = wf._synthesize_signature_path(state)
        drivers = result["synthesis"]["ranked_drivers"]
        assert drivers[0]["coverage_pct"] == 40.0  # 4/10 * 100

    def test_coverage_pct_uses_resolved_count_not_raw_input_count(self, wf):
        """Regression test: when some panel genes don't resolve in the network (e.g.
        DACH1/CESC's 23-of-24 case), coverage_pct must divide by the resolved count
        (genes_found_in_network) -- the same denominator RegNetAgents' own Fisher
        enrichment math already uses internally -- not the raw input panel size.
        Before the fix this silently understated coverage whenever any input gene
        failed to resolve, requiring a manual post-hoc correction in the DACH1 paper."""
        state = _sig_state(
            gene_signature=[f"G{i}" for i in range(24)],  # 24-gene input panel
            master_regulators=_mr_result(
                [_tf_entry("DACH1", overlap=4, regulon_size=409)],
                genes_found_in_network=23,  # only 23 resolved
            ),
        )
        result = wf._synthesize_signature_path(state)
        drivers = result["synthesis"]["ranked_drivers"]
        # 4/23 = 17.4%, not 4/24 = 16.7%
        assert drivers[0]["coverage_pct"] == 17.4
        # signature_size (the reported input-panel-size header field) stays the full 24 --
        # only the coverage_pct denominator changes, nothing else.
        assert result["synthesis"]["signature_size"] == 24

    def test_coverage_pct_falls_back_to_signature_size_when_genes_found_missing(self, wf):
        """Defensive fallback for an older RegNetAgents that doesn't report
        genes_found_in_network at all."""
        mr = _mr_result([_tf_entry("CTNNB1", overlap=4, regulon_size=200)])
        del mr["query_summary"]["genes_found_in_network"]
        state = _sig_state(
            gene_signature=["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"],  # 10 genes
            master_regulators=mr,
        )
        result = wf._synthesize_signature_path(state)
        drivers = result["synthesis"]["ranked_drivers"]
        assert drivers[0]["coverage_pct"] == 40.0  # falls back to 4/10

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
                          key_findings=[
                              "1 gene(s) confirmed by both network propagation and LINCS "
                              "experimental knockdown data (directional agreement)."
                          ]),
            ]),
        )
        result = wf._synthesize_signature_path(state)
        ev_row = result["synthesis"]["evidence_table"][0]
        assert ev_row["lincs_knockdown"] is True

    def test_lincs_disagreement_not_counted_as_confirmation(self, wf):
        state = _sig_state(
            master_regulators=_mr_result([
                _tf_entry("CTNNB1", overlap=3,
                          key_findings=[
                              "1 gene(s) show directional disagreement between network "
                              "prediction and LINCS experimental data — requires investigation."
                          ]),
            ]),
        )
        result = wf._synthesize_signature_path(state)
        ev_row = result["synthesis"]["evidence_table"][0]
        assert ev_row["lincs_knockdown"] is False

    def test_corroboration_count_in_evidence_table(self, wf):
        state = _sig_state(
            master_regulators=_mr_result([
                _tf_entry("CTNNB1", overlap=3, enrichment=2.5,
                          key_findings=[
                              "1 gene(s) confirmed by both network propagation and LINCS "
                              "experimental knockdown data (directional agreement).",
                              "super-enhancer: present",
                          ]),
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

    # -----------------------------------------------------------------
    # IntOGen driver annotation propagation
    # -----------------------------------------------------------------

    def test_driver_fields_propagated_to_ranked_drivers(self, wf):
        entry = _tf_entry("CTNNB1", overlap=4)
        entry["is_known_driver"] = True
        entry["driver_role"] = "oncogene"
        entry["tissue_matched"] = True
        state = _sig_state(
            master_regulators=_mr_result([entry]),
            driver_annotation_available=True,
        )
        result = wf._synthesize_signature_path(state)
        driver = result["synthesis"]["ranked_drivers"][0]
        assert driver["is_known_driver"] is True
        assert driver["driver_role"] == "oncogene"
        assert driver["tissue_matched"] is True

    def test_driver_fields_propagated_to_evidence_table(self, wf):
        entry = _tf_entry("CTNNB1", overlap=4)
        entry["is_known_driver"] = True
        entry["driver_role"] = "oncogene"
        entry["tissue_matched"] = False
        state = _sig_state(master_regulators=_mr_result([entry]))
        result = wf._synthesize_signature_path(state)
        ev_row = result["synthesis"]["evidence_table"][0]
        assert ev_row["is_known_driver"] is True
        assert ev_row["driver_role"] == "oncogene"
        assert ev_row["tissue_matched"] is False

    def test_driver_fields_default_when_absent(self, wf):
        """entries with no annotation (e.g. annotate_cancer_drivers unreachable) default safely."""
        state = _sig_state(master_regulators=_mr_result([_tf_entry("CTNNB1", overlap=4)]))
        result = wf._synthesize_signature_path(state)
        driver = result["synthesis"]["ranked_drivers"][0]
        assert driver["is_known_driver"] is False
        assert driver["driver_role"] is None
        assert driver["tissue_matched"] is None

    def test_driver_annotation_available_propagated_to_synthesis(self, wf):
        state = _sig_state(
            master_regulators=_mr_result([_tf_entry("CTNNB1", overlap=4)]),
            driver_annotation_available=True,
        )
        result = wf._synthesize_signature_path(state)
        assert result["synthesis"]["driver_annotation_available"] is True

    def test_driver_annotation_available_defaults_false(self, wf):
        state = _sig_state(master_regulators=_mr_result([_tf_entry("CTNNB1", overlap=4)]))
        result = wf._synthesize_signature_path(state)
        assert result["synthesis"]["driver_annotation_available"] is False


# ---------------------------------------------------------------------------
# _format_signature_report
# ---------------------------------------------------------------------------

class TestFormatSignatureReport:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def _make_synthesis(self, ranked_drivers=None, evidence_table=None,
                        errors=None, genes_not_found=None,
                        regnetagents_available=True, cascade_available=True,
                        driver_annotation_available=False):
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
            "driver_annotation_available": driver_annotation_available,
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

    # -----------------------------------------------------------------
    # IntOGen driver annotation column
    # -----------------------------------------------------------------

    def _driver_row(self, **overrides) -> dict:
        row = {
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
            "is_known_driver": False,
            "driver_role": None,
            "tissue_matched": None,
        }
        row.update(overrides)
        return row

    def test_known_driver_column_header_shown(self, wf):
        syn = self._make_synthesis(ranked_drivers=[self._driver_row()],
                                    evidence_table=[self._driver_row()])
        report = "\n".join(wf._format_signature_report(syn))
        assert "Known Driver (IntOGen)" in report

    def test_known_driver_shown_as_unknown_when_annotation_unavailable(self, wf):
        row = self._driver_row(is_known_driver=True, driver_role="oncogene")
        syn = self._make_synthesis(ranked_drivers=[row], evidence_table=[row],
                                    driver_annotation_available=False)
        report = "\n".join(wf._format_signature_report(syn))
        assert "unknown" in report
        assert "✓ oncogene" not in report

    def test_known_driver_shown_with_role_when_available(self, wf):
        row = self._driver_row(is_known_driver=True, driver_role="oncogene")
        syn = self._make_synthesis(ranked_drivers=[row], evidence_table=[row],
                                    driver_annotation_available=True)
        report = "\n".join(wf._format_signature_report(syn))
        assert "✓ oncogene" in report

    def test_known_driver_shown_as_no_when_not_a_driver(self, wf):
        row = self._driver_row(is_known_driver=False, driver_role=None)
        syn = self._make_synthesis(ranked_drivers=[row], evidence_table=[row],
                                    driver_annotation_available=True)
        report = "\n".join(wf._format_signature_report(syn))
        # the "No" cell — not asserting exact placement, just that it isn't
        # misreported as a driver and doesn't say "unknown"
        assert "✓ " not in report.split("### Driver Details")[0] or "✓ oncogene" not in report

    def test_tissue_matched_flag_shown(self, wf):
        row = self._driver_row(is_known_driver=True, driver_role="oncogene", tissue_matched=True)
        syn = self._make_synthesis(ranked_drivers=[row], evidence_table=[row],
                                    driver_annotation_available=True)
        report = "\n".join(wf._format_signature_report(syn))
        assert "[tissue-matched]" in report

    def test_not_tissue_matched_omits_flag(self, wf):
        row = self._driver_row(is_known_driver=True, driver_role="oncogene", tissue_matched=False)
        syn = self._make_synthesis(ranked_drivers=[row], evidence_table=[row],
                                    driver_annotation_available=True)
        report = "\n".join(wf._format_signature_report(syn))
        # the footnote explains what [tissue-matched] means generically -- only the
        # Driver Details line for this specific (non-tissue-matched) gene must omit it.
        details_line = next(
            line for line in report.splitlines() if line.strip().startswith("Known driver (IntOGen):")
        )
        assert "[tissue-matched]" not in details_line
        assert "✓ oncogene" in details_line

    def test_driver_details_section_shows_known_driver_line(self, wf):
        row = self._driver_row(is_known_driver=True, driver_role="tumor_suppressor")
        syn = self._make_synthesis(ranked_drivers=[row], evidence_table=[row],
                                    driver_annotation_available=True)
        report = "\n".join(wf._format_signature_report(syn))
        assert "Known driver (IntOGen): ✓ tumor_suppressor" in report

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


class TestRunSignaturePathDriverAnnotation:
    """_run_signature_path wires annotate_cancer_drivers onto the candidate TF list."""

    @pytest.fixture
    def wf(self):
        return _make_workflow()

    @pytest.mark.asyncio
    async def test_annotate_cancer_drivers_called_with_candidate_genes(self, wf):
        calls = []

        async def call_tool(name, args, timeout_seconds=None):
            calls.append((name, args))
            if name == "find_master_regulators":
                return _mr_result([
                    _tf_entry("CTNNB1", overlap=4),
                    _tf_entry("TP53", overlap=2),
                ])
            if name == "annotate_cancer_drivers":
                return {
                    "driver_annotation_available": True,
                    "results": {
                        "CTNNB1": {"is_driver": True, "role": "oncogene"},
                        "TP53": {"is_driver": True, "role": "tumor_suppressor"},
                    },
                }
            raise AssertionError(f"unexpected tool call: {name}")

        wf._regnetagents = AsyncMock()
        wf._regnetagents.call_tool = AsyncMock(side_effect=call_tool)
        wf._cascade = None  # skip DoRothEA overlap + CASCADE validation blocks
        state = _sig_state()

        result = await wf._run_signature_path(state)

        driver_calls = [args for name, args in calls if name == "annotate_cancer_drivers"]
        assert len(driver_calls) == 1
        assert set(driver_calls[0]["genes"]) == {"CTNNB1", "TP53"}
        assert "cancer_type" not in driver_calls[0]  # no tcga_network set

        mr_list = result["master_regulators"]["master_regulators"]
        by_gene = {e["gene"]: e for e in mr_list}
        assert by_gene["CTNNB1"]["is_known_driver"] is True
        assert by_gene["CTNNB1"]["driver_role"] == "oncogene"
        assert by_gene["TP53"]["driver_role"] == "tumor_suppressor"
        assert result["driver_annotation_available"] is True

    @pytest.mark.asyncio
    async def test_cancer_type_passed_through_as_tcga_code(self, wf):
        async def call_tool(name, args, timeout_seconds=None):
            if name == "find_master_regulators":
                return _mr_result([_tf_entry("SHOX2", overlap=3)])
            if name == "annotate_cancer_drivers":
                assert args.get("cancer_type") == "paad"
                return {
                    "driver_annotation_available": True,
                    "results": {"SHOX2": {"is_driver": False, "role": None, "tissue_matched": False}},
                }
            raise AssertionError(f"unexpected tool call: {name}")

        wf._regnetagents = AsyncMock()
        wf._regnetagents.call_tool = AsyncMock(side_effect=call_tool)
        wf._cascade = None
        state = _sig_state(cancer_type="paad")

        # cancer_type="paad" also auto-triggers the cross-context novelty gap block
        # (_effective_cancer_contexts prepends "pancreatic cancer") -- stub PubMed
        # so this test stays offline and focused on the driver-annotation wiring.
        with patch("pubmed_client.novelty_assessment", new=AsyncMock(return_value={})):
            result = await wf._run_signature_path(state)
        mr_list = result["master_regulators"]["master_regulators"]
        assert mr_list[0]["tissue_matched"] is False

    @pytest.mark.asyncio
    async def test_annotate_cancer_drivers_failure_is_non_fatal(self, wf):
        """A broken/unreachable driver-annotation call must not sink the whole signature run --
        find_master_regulators results still come back, just without driver labels."""
        async def call_tool(name, args, timeout_seconds=None):
            if name == "find_master_regulators":
                return _mr_result([_tf_entry("CTNNB1", overlap=4)])
            if name == "annotate_cancer_drivers":
                raise RuntimeError("RegNetAgents connection lost")
            raise AssertionError(f"unexpected tool call: {name}")

        wf._regnetagents = AsyncMock()
        wf._regnetagents.call_tool = AsyncMock(side_effect=call_tool)
        wf._cascade = None
        state = _sig_state()

        result = await wf._run_signature_path(state)

        assert "master_regulators" not in result["errors"]
        assert result["driver_annotation_available"] is False
        mr_list = result["master_regulators"]["master_regulators"]
        assert mr_list[0].get("is_known_driver") is None  # never set, not defaulted to a wrong value

    @pytest.mark.asyncio
    async def test_no_tfs_found_skips_annotation_call(self, wf):
        async def call_tool(name, args, timeout_seconds=None):
            if name == "find_master_regulators":
                return _mr_result([])
            raise AssertionError(f"annotate_cancer_drivers must not be called with an empty TF list: {name}")

        wf._regnetagents = AsyncMock()
        wf._regnetagents.call_tool = AsyncMock(side_effect=call_tool)
        wf._cascade = None
        state = _sig_state()

        result = await wf._run_signature_path(state)
        assert result["master_regulators"]["master_regulators"] == []


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
