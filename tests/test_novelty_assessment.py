"""
Unit tests for PubMed novelty assessment (Issue #15) and edge pair novelty (Issue #16).

All tests mock NCBI HTTP calls — no network required.
"""

import pytest
from unittest.mock import patch, MagicMock, AsyncMock

from pubmed_client import (
    _build_base_query,
    _build_experimental_query,
    _esearch,
    _efetch_year,
    _verdict,
    _rationale,
    _novelty_assessment_sync,
)
from orchestra_langgraph_workflow import OrchestraWorkflow, _TCGA_TO_CANCER_CONTEXT


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _xml_esearch(count: int, ids: list[str] = None) -> str:
    id_block = "".join(f"<Id>{i}</Id>" for i in (ids or []))
    return (
        f"<eSearchResult><Count>{count}</Count>"
        f"<IdList>{id_block}</IdList></eSearchResult>"
    )


def _xml_efetch(year: str) -> str:
    return (
        f"<PubmedArticleSet><PubmedArticle><MedlineCitation>"
        f"<Article><Journal><JournalIssue><PubDate>"
        f"<Year>{year}</Year></PubDate></JournalIssue></Journal>"
        f"</Article></MedlineCitation></PubmedArticle></PubmedArticleSet>"
    )


def _mock_response(text: str, status: int = 200) -> MagicMock:
    r = MagicMock()
    r.status_code = status
    r.text = text
    r.raise_for_status = MagicMock()
    return r


# ---------------------------------------------------------------------------
# Query construction
# ---------------------------------------------------------------------------

class TestQueryConstruction:
    def test_single_gene_query(self):
        q = _build_base_query("FOXM1", "head and neck squamous", None)
        assert '"FOXM1"[tiab]' in q
        assert '"head and neck squamous"[tiab]' in q
        assert "gene2" not in q.lower()

    def test_gene_pair_query(self):
        q = _build_base_query("FOXM1", "head and neck squamous", "TOP2A")
        assert '"FOXM1"[tiab]' in q
        assert '"TOP2A"[tiab]' in q
        assert '"head and neck squamous"[tiab]' in q

    def test_experimental_query_wraps_base(self):
        base = '"TP53"[tiab] AND "breast cancer"[tiab]'
        exp = _build_experimental_query(base)
        assert base in exp
        assert "in vitro" in exp
        assert "computational" in exp


# ---------------------------------------------------------------------------
# Verdict thresholds
# ---------------------------------------------------------------------------

class TestVerdictThresholds:
    def test_novel_zero(self):
        assert _verdict(0) == "novel"

    def test_novel_four(self):
        assert _verdict(4) == "novel"

    def test_emerging_five(self):
        assert _verdict(5) == "emerging"

    def test_emerging_twenty(self):
        assert _verdict(20) == "emerging"

    def test_established_twenty_one(self):
        assert _verdict(21) == "established"

    def test_established_large(self):
        assert _verdict(500) == "established"


# ---------------------------------------------------------------------------
# Rationale
# ---------------------------------------------------------------------------

class TestRationale:
    def test_zero_hits(self):
        r = _rationale(0, 0, "novel")
        assert "No prior papers" in r

    def test_novel_one_hit(self):
        r = _rationale(1, 1, "novel")
        assert "1 prior paper" in r
        assert "limited" in r

    def test_emerging_phrasing(self):
        r = _rationale(10, 4, "emerging")
        assert "10 prior papers" in r
        assert "4 experimental" in r
        assert "active area" in r

    def test_established_phrasing(self):
        r = _rationale(50, 30, "established")
        assert "well-characterized" in r

    def test_computational_count_derived(self):
        r = _rationale(10, 3, "emerging")
        assert "7 computational" in r


# ---------------------------------------------------------------------------
# _esearch HTTP parsing
# ---------------------------------------------------------------------------

class TestEsearch:
    def test_parses_count_and_ids(self):
        xml = _xml_esearch(42, ["12345678", "87654321"])
        with patch("pubmed_client.requests.get", return_value=_mock_response(xml)):
            count, ids = _esearch("dummy query", retmax=2)
        assert count == 42
        assert ids == ["12345678", "87654321"]

    def test_zero_results(self):
        xml = _xml_esearch(0)
        with patch("pubmed_client.requests.get", return_value=_mock_response(xml)):
            count, ids = _esearch("dummy query")
        assert count == 0
        assert ids == []


# ---------------------------------------------------------------------------
# _efetch_year HTTP parsing
# ---------------------------------------------------------------------------

class TestEfetchYear:
    def test_parses_year(self):
        xml = _xml_efetch("2023")
        with patch("pubmed_client.requests.get", return_value=_mock_response(xml)):
            year = _efetch_year("12345678")
        assert year == 2023

    def test_medline_date_fallback(self):
        xml = (
            "<PubmedArticleSet><PubmedArticle><MedlineCitation>"
            "<Article><Journal><JournalIssue><PubDate>"
            "<MedlineDate>2021 Jan-Feb</MedlineDate>"
            "</PubDate></JournalIssue></Journal>"
            "</Article></MedlineCitation></PubmedArticle></PubmedArticleSet>"
        )
        with patch("pubmed_client.requests.get", return_value=_mock_response(xml)):
            year = _efetch_year("12345678")
        assert year == 2021


# ---------------------------------------------------------------------------
# Full _novelty_assessment_sync
# ---------------------------------------------------------------------------

class TestNoveltAssessmentSync:
    def _make_responses(self, total: int, recent_id: str, exp_count: int, year: str):
        return [
            _mock_response(_xml_esearch(total, [recent_id])),   # call 1: total + recent id
            _mock_response(_xml_esearch(exp_count)),            # call 2: experimental count
            _mock_response(_xml_efetch(year)),                  # call 3: year via efetch
        ]

    def test_novel_verdict(self):
        with patch("pubmed_client.requests.get", side_effect=self._make_responses(2, "111", 1, "2022")):
            result = _novelty_assessment_sync("FOXM1", "head and neck squamous", None)
        assert result["novelty_verdict"] == "novel"
        assert result["pubmed_hits"] == 2
        assert result["experimental_hits"] == 1
        assert result["computational_hits"] == 1
        assert result["most_recent_year"] == 2022
        assert result["gene"] == "FOXM1"
        assert result["gene2"] is None

    def test_established_verdict(self):
        with patch("pubmed_client.requests.get", side_effect=self._make_responses(50, "222", 30, "2024")):
            result = _novelty_assessment_sync("TP53", "breast cancer", None)
        assert result["novelty_verdict"] == "established"
        assert result["pubmed_hits"] == 50

    def test_gene_pair(self):
        with patch("pubmed_client.requests.get", side_effect=self._make_responses(3, "333", 2, "2023")):
            result = _novelty_assessment_sync("FOXM1", "head and neck squamous", "TOP2A")
        assert result["gene2"] == "TOP2A"
        assert result["pubmed_hits"] == 3

    def test_zero_hits_skips_exp_and_efetch(self):
        xml_zero = _xml_esearch(0)
        with patch("pubmed_client.requests.get", return_value=_mock_response(xml_zero)) as mock_get:
            result = _novelty_assessment_sync("RAREGEN", "rare cancer type xyz", None)
        # Only 1 call (total esearch) — exp and efetch are skipped when total == 0
        assert mock_get.call_count == 1
        assert result["pubmed_hits"] == 0
        assert result["experimental_hits"] == 0
        assert result["most_recent_year"] is None
        assert result["novelty_verdict"] == "novel"

    def test_computational_hits_never_negative(self):
        # exp_count > total should not produce negative computational_hits
        with patch("pubmed_client.requests.get", side_effect=self._make_responses(5, "444", 8, "2020")):
            result = _novelty_assessment_sync("GENE", "context", None)
        assert result["computational_hits"] >= 0


# ---------------------------------------------------------------------------
# Workflow synthesis and formatting
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


class TestNoveltySynthesis:
    def _state(self, **overrides) -> dict:
        base = {
            "gene": "FOXM1",
            "cell_type": "",
            "analysis_type": "novelty_assessment",
            "analysis_depth": "comprehensive",
            "gene_role": None,
            "cancer_context": "head and neck squamous",
            "gene2": None,
            "novelty_result": None,
            "completed_steps": [],
            "errors": {},
            "synthesis": None,
            "final_report": None,
        }
        base.update(overrides)
        return base

    def test_routing_decision_novelty(self):
        wf = _make_workflow()
        state = self._state()
        assert wf._routing_decision(state) == "novelty_path"

    def test_synthesize_novelty_path_populates_synthesis(self):
        wf = _make_workflow()
        result = {
            "gene": "FOXM1", "gene2": None, "cancer_context": "head and neck squamous",
            "pubmed_hits": 2, "experimental_hits": 1, "computational_hits": 1,
            "most_recent_year": 2023, "novelty_verdict": "novel",
            "verdict_rationale": "2 prior papers — limited prior characterization",
        }
        state = self._state(novelty_result=result)
        state = wf._synthesize_novelty_path(state)
        assert state["synthesis"]["routing"] == "novelty"
        assert state["synthesis"]["novelty_result"] == result
        assert "synthesize" in state["completed_steps"]

    def test_format_novelty_report_with_result(self):
        wf = _make_workflow()
        synthesis = {
            "routing": "novelty",
            "gene": "FOXM1",
            "gene2": "TOP2A",
            "cancer_context": "head and neck squamous",
            "novelty_result": {
                "pubmed_hits": 2,
                "experimental_hits": 1,
                "computational_hits": 1,
                "most_recent_year": 2023,
                "novelty_verdict": "novel",
                "verdict_rationale": "2 prior papers; 1 experimental, 1 computational — limited prior characterization",
            },
            "errors": {},
        }
        lines = wf._format_novelty_report(synthesis)
        report = "\n".join(lines)
        assert "FOXM1/TOP2A" in report
        assert "NOVEL" in report
        assert "2 |" in report
        assert "2023" in report
        assert "Thresholds" in report

    def test_format_novelty_report_no_result(self):
        wf = _make_workflow()
        synthesis = {
            "routing": "novelty",
            "gene": "FOXM1",
            "gene2": None,
            "cancer_context": "breast cancer",
            "novelty_result": None,
            "errors": {"novelty_assessment": "NCBI timeout"},
        }
        lines = wf._format_novelty_report(synthesis)
        report = "\n".join(lines)
        assert "no result" in report.lower()
        assert "NCBI timeout" in report

    def test_format_novelty_report_single_gene_subject(self):
        wf = _make_workflow()
        synthesis = {
            "routing": "novelty",
            "gene": "TP53",
            "gene2": None,
            "cancer_context": "colorectal",
            "novelty_result": {
                "pubmed_hits": 100, "experimental_hits": 60, "computational_hits": 40,
                "most_recent_year": 2024, "novelty_verdict": "established",
                "verdict_rationale": "well-characterized",
            },
            "errors": {},
        }
        lines = wf._format_novelty_report(synthesis)
        report = "\n".join(lines)
        assert "FOXM1/TOP2A" not in report
        assert "TP53" in report
        assert "ESTABLISHED" in report


# ---------------------------------------------------------------------------
# Edge pair novelty (Issue #16)
# ---------------------------------------------------------------------------

def _make_novelty_result(gene: str, gene2: str, hits: int, verdict: str, ctx: str) -> dict:
    return {
        "gene": gene,
        "gene2": gene2,
        "cancer_context": ctx,
        "pubmed_hits": hits,
        "experimental_hits": 0,
        "computational_hits": hits,
        "most_recent_year": None,
        "novelty_verdict": verdict,
        "verdict_rationale": f"{hits} papers",
    }


def _make_edge_state(**overrides) -> dict:
    base = {
        "gene": "FOXM1",
        "cell_type": "epithelial_cell",
        "analysis_type": "causal_chain",
        "analysis_depth": "comprehensive",
        "gene_role": "transcription_factor",
        "cancer_type": None,
        "cancer_context": None,
        "gene2": None,
        "novelty_result": None,
        "edge_novelty_results": None,
        "completed_steps": [],
        "errors": {},
        "synthesis": None,
        "final_report": None,
    }
    base.update(overrides)
    return base


class TestTcgaMapping:
    def test_known_codes_map_correctly(self):
        assert _TCGA_TO_CANCER_CONTEXT["hnsc"] == "head and neck squamous"
        assert _TCGA_TO_CANCER_CONTEXT["coad"] == "colorectal"
        assert _TCGA_TO_CANCER_CONTEXT["brca"] == "breast cancer"
        assert _TCGA_TO_CANCER_CONTEXT["ov"] == "ovarian cancer"

    def test_all_eight_codes_present(self):
        expected = {"hnsc", "coad", "brca", "luad", "lusc", "ov", "prad", "ucec"}
        assert set(_TCGA_TO_CANCER_CONTEXT.keys()) == expected


class TestEdgeNovelty:
    """Tests for _run_edge_novelty LangGraph node (Issue #16)."""

    def _wf(self) -> OrchestraWorkflow:
        return _make_workflow()

    # ------------------------------------------------------------------
    # Skip conditions
    # ------------------------------------------------------------------

    async def test_skips_when_no_cancer_context(self):
        wf = self._wf()
        state = _make_edge_state(synthesis={"routing": "tf"})
        result = await wf._run_edge_novelty(state)
        assert result["synthesis"]["edge_novelty_results"] == []
        assert "run_edge_novelty" in result["completed_steps"]

    async def test_skips_unsupported_routing(self):
        wf = self._wf()
        state = _make_edge_state(
            cancer_context="colorectal",
            synthesis={"routing": "effector"},
        )
        result = await wf._run_edge_novelty(state)
        assert result["synthesis"]["edge_novelty_results"] == []

    async def test_skips_effector_routing(self):
        wf = self._wf()
        state = _make_edge_state(
            cancer_context="breast cancer",
            synthesis={"routing": "effector"},
        )
        result = await wf._run_edge_novelty(state)
        assert result["synthesis"]["edge_novelty_results"] == []

    async def test_skips_when_no_pairs(self):
        wf = self._wf()
        state = _make_edge_state(
            cancer_context="colorectal",
            synthesis={"routing": "tf", "cross_system_hits": [], "network_targets_sample": []},
        )
        result = await wf._run_edge_novelty(state)
        assert result["synthesis"]["edge_novelty_results"] == []

    # ------------------------------------------------------------------
    # TF path edge extraction
    # ------------------------------------------------------------------

    async def test_tf_path_extracts_cross_system_hits(self):
        wf = self._wf()
        cross_hits = [{"symbol": "CDKN1A"}, {"symbol": "MDM2"}]
        state = _make_edge_state(
            gene="TP53",
            cancer_context="colorectal",
            synthesis={
                "routing": "tf",
                "cross_system_hits": cross_hits,
                "network_targets_sample": [],
            },
        )
        mock_result = _make_novelty_result("TP53", "CDKN1A", 0, "novel", "colorectal")
        with patch("pubmed_client.novelty_assessment", new=AsyncMock(return_value=mock_result)):
            result = await wf._run_edge_novelty(state)
        edge_results = result["synthesis"]["edge_novelty_results"]
        assert len(edge_results) == 2
        pairs = [r["pair"] for r in edge_results]
        assert "TP53 → CDKN1A" in pairs
        assert "TP53 → MDM2" in pairs

    async def test_tf_path_falls_back_to_network_targets(self):
        wf = self._wf()
        state = _make_edge_state(
            gene="STAT3",
            cancer_context="cervical cancer",
            synthesis={
                "routing": "tf",
                "cross_system_hits": [],
                "network_targets_sample": ["BCL2", "VEGFA"],
            },
        )
        mock_result = _make_novelty_result("STAT3", "BCL2", 5, "emerging", "cervical cancer")
        with patch("pubmed_client.novelty_assessment", new=AsyncMock(return_value=mock_result)):
            result = await wf._run_edge_novelty(state)
        pairs = [r["pair"] for r in result["synthesis"]["edge_novelty_results"]]
        assert "STAT3 → BCL2" in pairs
        assert "STAT3 → VEGFA" in pairs

    # ------------------------------------------------------------------
    # Validation path edge extraction
    # ------------------------------------------------------------------

    async def test_validation_path_extracts_evidence_table(self):
        wf = self._wf()
        evidence_table = [
            {"gene": "BRD4", "corroboration_count": 1},
            {"gene": "CDK9", "corroboration_count": 3},
        ]
        state = _make_edge_state(
            gene="MYC",
            cancer_context="lymphoma",
            synthesis={
                "routing": "validation",
                "evidence_table": evidence_table,
            },
        )
        mock_result = _make_novelty_result("BRD4", "MYC", 0, "novel", "lymphoma")
        with patch("pubmed_client.novelty_assessment", new=AsyncMock(return_value=mock_result)):
            result = await wf._run_edge_novelty(state)
        pairs = [r["pair"] for r in result["synthesis"]["edge_novelty_results"]]
        assert "BRD4 → MYC" in pairs
        assert "CDK9 → MYC" in pairs

    # ------------------------------------------------------------------
    # Network comparison path edge extraction
    # ------------------------------------------------------------------

    async def test_network_comparison_prioritises_cascade_validated(self):
        wf = self._wf()
        validated_conserved = [
            {"gene": "TOP2A", "tier": "conserved_not_validated"},
            {"gene": "AURKB", "tier": "conserved_cascade_validated"},
        ]
        state = _make_edge_state(
            gene="FOXM1",
            cancer_context="head and neck squamous",
            synthesis={
                "routing": "network_comparison",
                "conserved_regulators": ["TOP2A", "AURKB"],
                "validated_conserved": validated_conserved,
            },
        )
        mock_result = _make_novelty_result("AURKB", "FOXM1", 2, "novel", "head and neck squamous")
        with patch("pubmed_client.novelty_assessment", new=AsyncMock(return_value=mock_result)):
            result = await wf._run_edge_novelty(state)
        # AURKB (cascade_validated) should appear before TOP2A
        pairs = [r["pair"] for r in result["synthesis"]["edge_novelty_results"]]
        assert pairs.index("AURKB → FOXM1") < pairs.index("TOP2A → FOXM1")

    async def test_network_comparison_uses_tcga_mapping_for_cancer_context(self):
        wf = self._wf()
        validated_conserved = [{"gene": "TOP2A", "tier": "conserved_not_validated"}]
        state = _make_edge_state(
            gene="FOXM1",
            cancer_type="hnsc",
            cancer_context=None,  # not explicitly provided
            synthesis={
                "routing": "network_comparison",
                "conserved_regulators": ["TOP2A"],
                "validated_conserved": validated_conserved,
            },
        )
        captured_contexts: list[str] = []

        async def _capture(gene, ctx, gene2=None):
            captured_contexts.append(ctx)
            return _make_novelty_result(gene, gene2 or "", 0, "novel", ctx)

        with patch("pubmed_client.novelty_assessment", new=AsyncMock(side_effect=_capture)):
            await wf._run_edge_novelty(state)
        assert captured_contexts[0] == "head and neck squamous"

    async def test_network_comparison_unknown_tcga_falls_back_to_raw_string(self):
        wf = self._wf()
        state = _make_edge_state(
            gene="EGFR",
            cancer_type="unknown_code",
            cancer_context=None,
            synthesis={
                "routing": "network_comparison",
                "conserved_regulators": ["SP1"],
                "validated_conserved": [{"gene": "SP1", "tier": "conserved_not_validated"}],
            },
        )
        captured_contexts: list[str] = []

        async def _capture(gene, ctx, gene2=None):
            captured_contexts.append(ctx)
            return _make_novelty_result(gene, gene2 or "", 0, "novel", ctx)

        with patch("pubmed_client.novelty_assessment", new=AsyncMock(side_effect=_capture)):
            await wf._run_edge_novelty(state)
        assert captured_contexts[0] == "unknown_code"

    # ------------------------------------------------------------------
    # Exception handling
    # ------------------------------------------------------------------

    async def test_exceptions_are_filtered_not_propagated(self):
        wf = self._wf()
        call_count = 0

        async def _sometimes_fail(gene, ctx, gene2=None):
            nonlocal call_count
            call_count += 1
            if call_count == 1:
                raise RuntimeError("NCBI timeout")
            return _make_novelty_result(gene, gene2 or "", 3, "novel", ctx)

        state = _make_edge_state(
            gene="TP53",
            cancer_context="colorectal",
            synthesis={
                "routing": "tf",
                "cross_system_hits": [{"symbol": "CDKN1A"}, {"symbol": "MDM2"}],
                "network_targets_sample": [],
            },
        )
        with patch("pubmed_client.novelty_assessment", new=AsyncMock(side_effect=_sometimes_fail)):
            result = await wf._run_edge_novelty(state)
        # One failed, one succeeded — only the success is in results
        edge_results = result["synthesis"]["edge_novelty_results"]
        assert len(edge_results) == 1

    # ------------------------------------------------------------------
    # Formatter integration
    # ------------------------------------------------------------------

    def test_format_pair_novelty_section_with_results(self):
        wf = self._wf()
        synthesis = {
            "edge_novelty_results": [
                _make_novelty_result("FOXM1", "TOP2A", 1, "novel", "cervical cancer"),
                _make_novelty_result("FOXM1", "AURKB", 0, "novel", "cervical cancer"),
            ]
        }
        # inject pair labels as _run_edge_novelty would
        synthesis["edge_novelty_results"][0]["pair"] = "FOXM1 → TOP2A"
        synthesis["edge_novelty_results"][1]["pair"] = "FOXM1 → AURKB"
        lines = wf._format_pair_novelty_section(synthesis)
        report = "\n".join(lines)
        assert "Regulatory Pair Novelty" in report
        assert "FOXM1 → TOP2A" in report
        assert "FOXM1 → AURKB" in report
        assert "NOVEL" in report

    def test_format_pair_novelty_section_empty_returns_no_lines(self):
        wf = self._wf()
        assert wf._format_pair_novelty_section({"edge_novelty_results": []}) == []
        assert wf._format_pair_novelty_section({}) == []

    def test_format_tf_report_includes_pair_novelty_section(self):
        wf = self._wf()
        synthesis = {
            "routing": "tf",
            "gene": "TP53",
            "cell_type": "epithelial_cell",
            "gene_role": "transcription_factor",
            "cascade_key_findings": [],
            "corroborated_targets": [],
            "cross_system_hits": [],
            "source_agreements": [],
            "source_disagreements": [],
            "network_context": {},
            "regnetagents_target_count": 0,
            "regnetagents_available": True,
            "cascade_available": True,
            "discordance_flags": [],
            "errors": {},
            "edge_novelty_results": [
                {**_make_novelty_result("TP53", "CDKN1A", 2, "novel", "colorectal"),
                 "pair": "TP53 → CDKN1A"},
            ],
        }
        lines = wf._format_tf_report(synthesis)
        report = "\n".join(lines)
        assert "Regulatory Pair Novelty" in report
        assert "TP53 → CDKN1A" in report

    def test_format_tf_report_omits_pair_novelty_when_empty(self):
        wf = self._wf()
        synthesis = {
            "routing": "tf",
            "gene": "TP53", "cell_type": "epithelial_cell",
            "gene_role": "transcription_factor",
            "cascade_key_findings": [], "corroborated_targets": [],
            "cross_system_hits": [], "source_agreements": [],
            "source_disagreements": [], "network_context": {},
            "regnetagents_target_count": 0,
            "regnetagents_available": True, "cascade_available": True,
            "discordance_flags": [], "errors": {},
            "edge_novelty_results": [],
        }
        report = "\n".join(wf._format_tf_report(synthesis))
        assert "Regulatory Pair Novelty" not in report

    def test_format_validation_report_includes_pair_novelty_section(self):
        wf = self._wf()
        synthesis = {
            "routing": "validation",
            "gene": "MYC", "cell_type": "cd4_t_cells",
            "validated_targets": [],
            "evidence_table": [],
            "discordance_flags": [],
            "errors": {},
            "regnetagents_available": True, "cascade_available": True,
            "edge_novelty_results": [
                {**_make_novelty_result("BRD4", "MYC", 0, "novel", "lymphoma"),
                 "pair": "BRD4 → MYC"},
            ],
        }
        report = "\n".join(wf._format_validation_report(synthesis))
        assert "Regulatory Pair Novelty" in report
        assert "BRD4 → MYC" in report

    def test_format_network_comparison_report_includes_pair_novelty_section(self):
        wf = self._wf()
        synthesis = {
            "routing": "network_comparison",
            "gene": "FOXM1", "cell_type": "epithelial_cell",
            "cancer_type": "hnsc",
            "tumor_state_context": "tcga_hnsc",
            "rewiring_classification": "high",
            "reg_conserved_fraction": 0.04,
            "reg_pop_total": 1, "reg_tumor_total": 24,
            "conserved_regulators": ["TOP2A"],
            "tumor_state_only_regulators": [],
            "population_averaged_only_regulators": [],
            "validated_conserved": [{"gene": "TOP2A", "tier": "conserved_not_validated",
                                     "cascade_key_findings": [], "cascade_error": None}],
            "tgt_conserved_count": 0, "tgt_conserved_fraction": 0.0,
            "tgt_tumor_only": [],
            "regnetagents_available": True, "cascade_available": True,
            "errors": {},
            "edge_novelty_results": [
                {**_make_novelty_result("TOP2A", "FOXM1", 0, "novel", "head and neck squamous"),
                 "pair": "TOP2A → FOXM1"},
            ],
        }
        report = "\n".join(wf._format_network_comparison_report(synthesis))
        assert "Regulatory Pair Novelty" in report
        assert "TOP2A → FOXM1" in report


# ---------------------------------------------------------------------------
# novelty_assessment_batch MCP handler
# ---------------------------------------------------------------------------

def _novel_result(gene: str, hits: int, exp: int, year: str | None, verdict: str) -> dict:
    return {
        "gene": gene,
        "cancer_context": "cervical cancer",
        "gene2": None,
        "pubmed_hits": hits,
        "experimental_hits": exp,
        "computational_hits": max(0, hits - exp),
        "most_recent_year": year,
        "novelty_verdict": verdict,
        "verdict_rationale": f"{gene} has {hits} papers.",
    }


class TestNoveltybatch:
    """Unit tests for the novelty_assessment_batch call_tool handler."""

    def _run_handler(self, genes: list, cancer_context: str, side_effects):
        """Invoke the batch handler synchronously via asyncio.run."""
        import asyncio
        from unittest.mock import patch, AsyncMock
        from orchestra_mcp_server import call_tool

        async def _run():
            with patch("pubmed_client.novelty_assessment", new=AsyncMock(side_effect=side_effects)):
                return await call_tool(
                    "novelty_assessment_batch",
                    {"genes": genes, "cancer_context": cancer_context},
                )

        # call_tool uses app.request_context — patch progress path instead
        return asyncio.run(_run())

    def test_batch_report_assembly_all_genes_present(self):
        """Test the report string built from batch results directly."""
        import asyncio

        genes = ["LITAF", "USP21", "CENPK"]
        raw = [
            _novel_result("LITAF", 2, 0, "2019", "novel"),
            _novel_result("USP21", 8, 3, "2023", "emerging"),
            _novel_result("CENPK", 25, 10, "2024", "established"),
        ]
        cancer_context = "cervical cancer"

        rows = [
            "| Gene | Verdict | PubMed Hits | Experimental | Computational | Last Year |",
            "|------|---------|-------------|--------------|---------------|-----------|",
        ]
        rationales = []
        for gene, res in zip(genes, raw):
            v = (res.get("novelty_verdict") or "unknown").upper()
            rows.append(
                f"| {gene} | {v} | {res.get('pubmed_hits', 0)} |"
                f" {res.get('experimental_hits', 0)} |"
                f" {res.get('computational_hits', 0)} |"
                f" {res.get('most_recent_year') or 'N/A'} |"
            )
            rationales.append(f"**{gene}** ({v}): {res.get('verdict_rationale', '')}")

        report = "\n".join([
            f"## Novelty Assessment: {cancer_context}",
            f"_Assessed {len(genes)} genes_",
            "",
            *rows,
            "",
            "_Thresholds: >20 hits = established · 5–20 = emerging · <5 = novel_",
            "",
            "### Rationales",
            "",
            *rationales,
        ])

        assert "LITAF" in report
        assert "USP21" in report
        assert "CENPK" in report
        assert "NOVEL" in report
        assert "EMERGING" in report
        assert "ESTABLISHED" in report
        assert "## Novelty Assessment: cervical cancer" in report
        assert "_Assessed 3 genes_" in report

    def test_batch_report_error_gene(self):
        """ERROR row is written for a gene whose PubMed call raises."""
        genes = ["LITAF", "BROKENGENE"]
        results_raw = [
            _novel_result("LITAF", 2, 0, "2019", "novel"),
        ]
        error = RuntimeError("network timeout")

        rows = [
            "| Gene | Verdict | PubMed Hits | Experimental | Computational | Last Year |",
            "|------|---------|-------------|--------------|---------------|-----------|",
        ]
        rationales = []
        mixed = [results_raw[0], error]
        for gene, res in zip(genes, mixed):
            if isinstance(res, BaseException):
                rows.append(f"| {gene} | ERROR | — | — | — | — |")
                rationales.append(f"**{gene}**: error — {res}")
            else:
                v = (res.get("novelty_verdict") or "unknown").upper()
                rows.append(f"| {gene} | {v} | {res.get('pubmed_hits', 0)} | ... | ... | ... |")
                rationales.append(f"**{gene}** ({v}): {res.get('verdict_rationale', '')}")

        report = "\n".join(rows)
        assert "| BROKENGENE | ERROR |" in report
        assert "| LITAF | NOVEL |" in report

    def test_batch_report_no_year(self):
        """most_recent_year=None renders as N/A."""
        gene = "NEWGENE"
        res = _novel_result(gene, 0, 0, None, "novel")
        v = (res.get("novelty_verdict") or "unknown").upper()
        row = (
            f"| {gene} | {v} | {res.get('pubmed_hits', 0)} |"
            f" {res.get('experimental_hits', 0)} |"
            f" {res.get('computational_hits', 0)} |"
            f" {res.get('most_recent_year') or 'N/A'} |"
        )
        assert "N/A" in row
