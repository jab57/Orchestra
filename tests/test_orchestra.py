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
# _is_in_pathway_enrichment
# ---------------------------------------------------------------------------

class TestIsInPathwayEnrichment:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def test_gene_found_in_pathway_gene_list(self, wf):
        network = {
            "pathway_analysis": {
                "enriched_pathways": [
                    {"name": "Wnt signaling", "genes": ["CTNNB1", "APC", "BRD4"]}
                ]
            }
        }
        assert wf._is_in_pathway_enrichment("BRD4", network) is True

    def test_gene_found_via_dict_symbol_field(self, wf):
        network = {
            "reactome": {
                "pathways": [
                    {"pathway": "MYC targets", "genes": [{"symbol": "BRD4"}]}
                ]
            }
        }
        assert wf._is_in_pathway_enrichment("BRD4", network) is True

    def test_gene_not_found(self, wf):
        network = {
            "pathway_analysis": {
                "enriched_pathways": [
                    {"name": "Wnt", "genes": ["CTNNB1", "APC"]}
                ]
            }
        }
        assert wf._is_in_pathway_enrichment("BRD4", network) is False

    def test_empty_network_returns_false(self, wf):
        assert wf._is_in_pathway_enrichment("BRD4", {}) is False

    def test_case_insensitive_match(self, wf):
        network = {
            "pathway_analysis": {
                "gene_set": ["brd4", "myc"]
            }
        }
        assert wf._is_in_pathway_enrichment("BRD4", network) is True

    def test_empty_gene_returns_false(self, wf):
        assert wf._is_in_pathway_enrichment("", {"pathway_analysis": {"genes": ["BRD4"]}}) is False


# ---------------------------------------------------------------------------
# _score_candidate_evidence
# ---------------------------------------------------------------------------

class TestScoreCandidateEvidence:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def _candidate(self, source="cascade_drug_discovery", pagerank=0.0,
                   key_findings=None):
        return {
            "gene": "BRD4",
            "source": source,
            "pagerank": pagerank,
            "key_findings": key_findings or [],
        }

    def test_pagerank_source_sets_flag(self, wf):
        c = self._candidate(source="regnetagents_pagerank", pagerank=0.05)
        scores = wf._score_candidate_evidence(c, {})
        assert scores["pagerank_rank"] is True

    def test_non_pagerank_source_clears_flag(self, wf):
        c = self._candidate(source="cascade_drug_discovery", pagerank=0.0)
        scores = wf._score_candidate_evidence(c, {})
        assert scores["pagerank_rank"] is False

    def test_lincs_flag_from_key_findings(self, wf):
        c = self._candidate(key_findings=["LINCS knockdown confirms downregulation"])
        scores = wf._score_candidate_evidence(c, {})
        assert scores["lincs_knockdown"] is True

    def test_depmap_essential_sets_flag(self, wf):
        c = self._candidate(key_findings=["DepMap: gene is essential in 42/50 cell lines"])
        scores = wf._score_candidate_evidence(c, {})
        assert scores["depmap_essentiality"] is True

    def test_depmap_not_essential_clears_flag(self, wf):
        c = self._candidate(key_findings=["DepMap: gene is not essential"])
        scores = wf._score_candidate_evidence(c, {})
        assert scores["depmap_essentiality"] is False

    def test_super_enhancer_flag(self, wf):
        c = self._candidate(key_findings=["MYC has super-enhancers in 32 cell types"])
        scores = wf._score_candidate_evidence(c, {})
        assert scores["super_enhancer"] is True

    def test_bet_inhibitor_sets_super_enhancer_flag(self, wf):
        c = self._candidate(key_findings=["Consider BET inhibitor therapy"])
        scores = wf._score_candidate_evidence(c, {})
        assert scores["super_enhancer"] is True

    def test_dorothea_flag(self, wf):
        c = self._candidate(key_findings=["DoRothEA tier A TF with high confidence"])
        scores = wf._score_candidate_evidence(c, {})
        assert scores["dorothea_tier"] is True

    def test_cbio_flag(self, wf):
        c = self._candidate(key_findings=["cBioPortal: high expression in TCGA BRCA"])
        scores = wf._score_candidate_evidence(c, {})
        assert scores["cbio_expression"] is True

    def test_corroboration_count_correct(self, wf):
        c = self._candidate(
            source="regnetagents_pagerank",
            pagerank=0.05,
            key_findings=[
                "LINCS knockdown confirmed",
                "super-enhancer present",
                "DoRothEA tier B",
            ],
        )
        scores = wf._score_candidate_evidence(c, {})
        # pagerank=✓, pathway=- (empty network), lincs=✓, depmap=-, se=✓, dorothea=✓, cbio=-
        assert scores["corroboration_count"] == 4
        assert scores["corroboration_denominator"] == 7

    def test_empty_findings_all_false(self, wf):
        c = self._candidate()
        scores = wf._score_candidate_evidence(c, {})
        assert scores["corroboration_count"] == 0

    def test_pathway_member_from_network(self, wf):
        network = {
            "pathway_analysis": {
                "enriched_pathways": [{"name": "MYC targets", "genes": ["BRD4", "MYC"]}]
            }
        }
        c = self._candidate()
        scores = wf._score_candidate_evidence(c, network)
        assert scores["pathway_member"] is True


# ---------------------------------------------------------------------------
# _synthesize_validation_path — evidence table + graceful degradation
# ---------------------------------------------------------------------------

class TestSynthesizeValidationPathEnhanced:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def _state(self, candidates, network=None, errors=None):
        return _tf_state(
            analysis_type="therapeutic_validation",
            validated_targets=candidates,
            network_analysis=network or {},
            errors=errors or {},
        )

    def test_evidence_table_computed_for_validated_candidates(self, wf):
        candidates = [
            {
                "gene": "BRD4",
                "source": "cascade_drug_discovery",
                "key_findings": ["super-enhancer present", "LINCS confirmed"],
                "multi_source_genes": [],
            }
        ]
        result = wf._synthesize_validation_path(self._state(candidates))
        table = result["synthesis"]["evidence_table"]
        assert len(table) == 1
        assert table[0]["gene"] == "BRD4"
        assert table[0]["lincs_knockdown"] is True
        assert table[0]["super_enhancer"] is True

    def test_unvalidated_candidates_excluded_from_table(self, wf):
        """Candidates beyond top 3 (no key_findings key) are not scored."""
        candidates = [
            {"gene": "BRD4", "source": "cascade_drug_discovery",
             "key_findings": [], "multi_source_genes": []},
            {"gene": "CDK9", "source": "cascade_ppi"},  # no key_findings
        ]
        result = wf._synthesize_validation_path(self._state(candidates))
        table = result["synthesis"]["evidence_table"]
        assert all(row["gene"] != "CDK9" for row in table)

    def test_evidence_table_sorted_by_corroboration(self, wf):
        candidates = [
            {"gene": "LOW", "source": "cascade_ppi",
             "key_findings": [], "multi_source_genes": []},
            {"gene": "HIGH", "source": "regnetagents_pagerank", "pagerank": 0.1,
             "key_findings": ["LINCS confirmed", "super-enhancer present"],
             "multi_source_genes": []},
        ]
        result = wf._synthesize_validation_path(self._state(candidates))
        table = result["synthesis"]["evidence_table"]
        assert table[0]["gene"] == "HIGH"

    def test_regnetagents_available_flag_true_when_network_present(self, wf):
        candidates = [{"gene": "X", "source": "cascade_ppi",
                       "key_findings": [], "multi_source_genes": []}]
        result = wf._synthesize_validation_path(
            self._state(candidates, network={"some": "data"})
        )
        assert result["synthesis"]["regnetagents_available"] is True

    def test_regnetagents_available_flag_false_on_network_error(self, wf):
        candidates = [{"gene": "X", "source": "cascade_ppi",
                       "key_findings": [], "multi_source_genes": []}]
        result = wf._synthesize_validation_path(
            self._state(candidates, errors={"network": "timeout"})
        )
        assert result["synthesis"]["regnetagents_available"] is False

    def test_cascade_available_false_when_all_cascade_errors(self, wf):
        candidates = [
            {"gene": "BRD4", "source": "cascade_drug_discovery",
             "cascade_error": "timeout", "key_findings": [], "multi_source_genes": []}
        ]
        result = wf._synthesize_validation_path(self._state(candidates))
        assert result["synthesis"]["cascade_available"] is False


# ---------------------------------------------------------------------------
# Graceful degradation warnings in report formatters
# ---------------------------------------------------------------------------

class TestGracefulDegradationWarnings:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def _tf_synthesis(self, **overrides):
        base = {
            "gene": "TP53", "cell_type": "epithelial_cell",
            "routing": "tf", "gene_role": "master_regulator",
            "cascade_key_findings": [], "corroborated_targets": [],
            "cross_system_hits": [], "source_agreements": [],
            "network_context": {}, "regnetagents_target_count": 0,
            "regnetagents_available": True, "cascade_available": True,
            "errors": {},
        }
        base.update(overrides)
        return base

    def _effector_synthesis(self, **overrides):
        base = {
            "gene": "APC", "cell_type": "epithelial_cell",
            "routing": "effector", "tf_partner": "CTNNB1",
            "cascade_key_findings": [], "corroborated_targets": [],
            "source_agreements": [], "source_disagreements": [],
            "network_context": {},
            "regnetagents_available": True, "cascade_available": True,
            "errors": {},
        }
        base.update(overrides)
        return base

    def _validation_synthesis(self, **overrides):
        base = {
            "gene": "MYC", "cell_type": "cd4_t_cells",
            "routing": "validation", "validated_targets": [],
            "evidence_table": [], "network_context": {},
            "regnetagents_available": True, "cascade_available": True,
            "errors": {},
        }
        base.update(overrides)
        return base

    def test_tf_report_warns_when_regnetagents_unavailable(self, wf):
        text = "\n".join(wf._format_tf_report(
            self._tf_synthesis(regnetagents_available=False)
        ))
        assert "RegNetAgents unavailable" in text

    def test_tf_report_warns_when_cascade_unavailable(self, wf):
        text = "\n".join(wf._format_tf_report(
            self._tf_synthesis(cascade_available=False)
        ))
        assert "CASCADE unavailable" in text

    def test_tf_report_no_warnings_when_both_available(self, wf):
        text = "\n".join(wf._format_tf_report(self._tf_synthesis()))
        assert "unavailable" not in text

    def test_effector_report_warns_when_cascade_unavailable(self, wf):
        text = "\n".join(wf._format_effector_report(
            self._effector_synthesis(cascade_available=False)
        ))
        assert "CASCADE unavailable" in text

    def test_validation_report_warns_when_regnetagents_unavailable(self, wf):
        text = "\n".join(wf._format_validation_report(
            self._validation_synthesis(regnetagents_available=False)
        ))
        assert "RegNetAgents unavailable" in text

    def test_validation_report_shows_corroboration_table(self, wf):
        synthesis = self._validation_synthesis(
            validated_targets=[
                {"gene": "BRD4", "source": "cascade_drug_discovery",
                 "key_findings": ["super-enhancer present"], "multi_source_genes": []}
            ],
            evidence_table=[{
                "gene": "BRD4",
                "pagerank_rank": False, "pathway_member": False,
                "lincs_knockdown": False, "depmap_essentiality": False,
                "super_enhancer": True, "dorothea_tier": False, "cbio_expression": False,
                "corroboration_count": 1, "corroboration_denominator": 7,
            }],
        )
        text = "\n".join(wf._format_validation_report(synthesis))
        assert "Evidence Corroboration Table" in text
        assert "BRD4" in text
        assert "1/7" in text


# ---------------------------------------------------------------------------
# Discordance flags
# ---------------------------------------------------------------------------

class TestComputeValidationDiscordanceFlags:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def _row(self, gene, pagerank=False, lincs=False, depmap=False,
             se=False, dorothea=False, cbio=False):
        return {
            "gene": gene,
            "pagerank_rank": pagerank,
            "pathway_member": False,
            "lincs_knockdown": lincs,
            "depmap_essentiality": depmap,
            "super_enhancer": se,
            "dorothea_tier": dorothea,
            "cbio_expression": cbio,
            "corroboration_count": sum([pagerank, lincs, depmap, se, dorothea, cbio]),
            "corroboration_denominator": 7,
        }

    def test_brd4_pattern_triggers_experimentally_active_not_in_network(self, wf):
        """BRD4: super-enhancer support but absent from ARACNe network."""
        table = [self._row("BRD4", se=True)]
        flags = wf._compute_validation_discordance_flags(table)
        types = [f["type"] for f in flags]
        assert "experimentally_active_not_in_network" in types
        brd4_flag = next(f for f in flags if f["type"] == "experimentally_active_not_in_network")
        assert brd4_flag["gene"] == "BRD4"

    def test_pagerank_gene_without_experimental_triggers_topological_hub(self, wf):
        """High-PageRank gene with no CASCADE evidence → topological hub not validated."""
        table = [self._row("CDK7", pagerank=True)]
        flags = wf._compute_validation_discordance_flags(table)
        types = [f["type"] for f in flags]
        assert "topological_hub_not_validated" in types
        flag = next(f for f in flags if f["type"] == "topological_hub_not_validated")
        assert flag["gene"] == "CDK7"

    def test_no_discordance_when_both_pagerank_and_experimental(self, wf):
        table = [self._row("MYC", pagerank=True, lincs=True)]
        flags = wf._compute_validation_discordance_flags(table)
        assert flags == []

    def test_no_discordance_when_neither_pagerank_nor_experimental(self, wf):
        table = [self._row("UNKNOWN")]
        flags = wf._compute_validation_discordance_flags(table)
        assert flags == []

    def test_multiple_candidates_can_each_flag(self, wf):
        table = [
            self._row("HUB", pagerank=True),              # topological_hub_not_validated
            self._row("BRD4", se=True),                   # experimentally_active_not_in_network
        ]
        flags = wf._compute_validation_discordance_flags(table)
        types = [f["type"] for f in flags]
        assert "topological_hub_not_validated" in types
        assert "experimentally_active_not_in_network" in types


class TestComputeTfDiscordanceFlags:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def _gene(self, symbol, sources=("lincs",)):
        return {"symbol": symbol, "source_count": len(sources), "sources": list(sources)}

    def test_cascade_only_gene_triggers_experimentally_active_flag(self, wf):
        """Gene in CASCADE multi_source_genes but NOT in rna_targets → discordance."""
        multi = [self._gene("BRD4")]
        rna = set()  # BRD4 absent from network
        flags = wf._compute_tf_discordance_flags(multi, rna, cross_system_hits=[])
        types = [f["type"] for f in flags]
        assert "experimentally_active_not_in_network" in types
        flag = next(f for f in flags if f["type"] == "experimentally_active_not_in_network")
        assert any(g["symbol"] == "BRD4" for g in flag["genes"])

    def test_network_only_gap_when_no_cross_system_hits(self, wf):
        """Large rna_targets, CASCADE has hits, but no overlap → network gap flag."""
        multi = [self._gene("CDKN1A")]
        rna = {"MYC", "BCL2", "MDM2"}  # CDKN1A absent from network
        flags = wf._compute_tf_discordance_flags(multi, rna, cross_system_hits=[])
        types = [f["type"] for f in flags]
        assert "network_topology_without_experimental_support" in types

    def test_no_discordance_when_all_genes_in_network(self, wf):
        """All CASCADE multi_source genes are also in rna_targets → no cascade_only flag."""
        multi = [self._gene("MYC")]
        rna = {"MYC", "BCL2"}
        cross_hits = [{"symbol": "MYC", "source_count": 2, "sources": ["lincs"]}]
        flags = wf._compute_tf_discordance_flags(multi, rna, cross_system_hits=cross_hits)
        assert flags == []

    def test_no_network_gap_when_empty_rna_targets(self, wf):
        """Empty rna_targets (e.g., effector routed wrong) → no network gap flag."""
        multi = [self._gene("CDKN1A")]
        flags = wf._compute_tf_discordance_flags(multi, rna_targets=set(), cross_system_hits=[])
        types = [f["type"] for f in flags]
        assert "network_topology_without_experimental_support" not in types


class TestDiscordanceFlagsInReports:
    @pytest.fixture
    def wf(self):
        return _make_workflow()

    def _tf_synthesis(self, discordance_flags=None, **overrides):
        base = {
            "gene": "TP53", "cell_type": "epithelial_cell",
            "routing": "tf", "gene_role": "master_regulator",
            "cascade_key_findings": [], "corroborated_targets": [],
            "cross_system_hits": [], "source_agreements": [],
            "network_context": {}, "regnetagents_target_count": 0,
            "regnetagents_available": True, "cascade_available": True,
            "discordance_flags": discordance_flags or [],
            "errors": {},
        }
        base.update(overrides)
        return base

    def _validation_synthesis(self, discordance_flags=None, **overrides):
        base = {
            "gene": "MYC", "cell_type": "cd4_t_cells",
            "routing": "validation", "validated_targets": [],
            "evidence_table": [], "network_context": {},
            "regnetagents_available": True, "cascade_available": True,
            "discordance_flags": discordance_flags or [],
            "errors": {},
        }
        base.update(overrides)
        return base

    def test_tf_report_shows_discordance_section(self, wf):
        flags = [{
            "type": "experimentally_active_not_in_network",
            "description": "BRD4 has CASCADE support but absent from network",
            "genes": [{"symbol": "BRD4", "source_count": 1, "sources": ["super_enhancer"]}],
        }]
        text = "\n".join(wf._format_tf_report(self._tf_synthesis(discordance_flags=flags)))
        assert "Notable Discordances" in text
        assert "BRD4" in text

    def test_tf_report_no_discordance_section_when_empty(self, wf):
        text = "\n".join(wf._format_tf_report(self._tf_synthesis()))
        assert "Notable Discordances" not in text

    def test_validation_report_shows_discordance_section(self, wf):
        flags = [{
            "type": "topological_hub_not_validated",
            "gene": "CDK7",
            "description": "CDK7 ranks in network but no CASCADE experimental support",
            "genes": [],
        }]
        synthesis = self._validation_synthesis(
            discordance_flags=flags,
            validated_targets=[
                {"gene": "CDK7", "source": "regnetagents_pagerank",
                 "key_findings": [], "multi_source_genes": [],
                 "pagerank": 0.05, "downstream_targets": 30}
            ],
        )
        text = "\n".join(wf._format_validation_report(synthesis))
        assert "Notable Discordances" in text
        assert "CDK7" in text

    def test_validation_report_no_discordance_section_when_empty(self, wf):
        text = "\n".join(wf._format_validation_report(self._validation_synthesis()))
        assert "Notable Discordances" not in text


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
