"""
Unit tests for PubMed novelty assessment (Issue #15).

All tests mock NCBI HTTP calls — no network required.
"""

import pytest
from unittest.mock import patch, MagicMock

from pubmed_client import (
    _build_base_query,
    _build_experimental_query,
    _esearch,
    _efetch_year,
    _verdict,
    _rationale,
    _novelty_assessment_sync,
)
from orchestra_langgraph_workflow import OrchestraWorkflow


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
