"""
Unit tests for cbioportal_client.py and the fetch_tcga_methylation_correlation
Orchestra tool handler. All HTTP calls are mocked — no network required.
"""

import math
import pytest
from unittest.mock import patch, MagicMock

from cbioportal_client import (
    _entrez_id,
    _discover_profiles,
    _get_sample_ids,
    _fetch_values,
    _correlation_sync,
    format_correlation_report,
    _TCGA_STUDY_IDS,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mock_get_response(data) -> MagicMock:
    r = MagicMock()
    r.raise_for_status = MagicMock()
    r.json.return_value = data
    return r


def _mock_post_response(data) -> MagicMock:
    r = MagicMock()
    r.raise_for_status = MagicMock()
    r.json.return_value = data
    return r


def _profile(pid: str, alt_type: str) -> dict:
    return {"molecularProfileId": pid, "molecularAlterationType": alt_type}


# ---------------------------------------------------------------------------
# Study ID mapping
# ---------------------------------------------------------------------------

class TestStudyIdMapping:
    def test_all_fourteen_codes_present(self):
        expected = {"blca", "brca", "cesc", "coad", "hnsc", "kirc", "lihc",
                    "luad", "lusc", "ov", "paad", "prad", "stad", "ucec"}
        assert set(_TCGA_STUDY_IDS.keys()) == expected

    def test_study_id_format(self):
        assert _TCGA_STUDY_IDS["cesc"] == "tcga_cesc"
        assert _TCGA_STUDY_IDS["hnsc"] == "tcga_hnsc"


# ---------------------------------------------------------------------------
# _entrez_id
# ---------------------------------------------------------------------------

class TestEntrezId:
    def test_resolves_known_gene(self):
        payload = {"entrezGeneId": 2353, "hugoGeneSymbol": "FOS"}
        with patch("cbioportal_client.requests.get", return_value=_mock_get_response(payload)):
            assert _entrez_id("FOS") == 2353

    def test_returns_none_on_http_error(self):
        with patch("cbioportal_client.requests.get", side_effect=Exception("404")):
            assert _entrez_id("FAKEGENE") is None

    def test_returns_none_when_field_missing(self):
        with patch("cbioportal_client.requests.get", return_value=_mock_get_response({})):
            assert _entrez_id("GENE") is None

    def test_returns_none_when_response_is_list(self):
        # API occasionally returns a list for unknown symbols
        with patch("cbioportal_client.requests.get", return_value=_mock_get_response([])):
            assert _entrez_id("GENE") is None


# ---------------------------------------------------------------------------
# _discover_profiles
# ---------------------------------------------------------------------------

class TestDiscoverProfiles:
    def _profiles_response(self, study: str = "tcga_cesc") -> list[dict]:
        return [
            _profile(f"{study}_rna_seq_v2_mrna", "MRNA_EXPRESSION"),
            _profile(f"{study}_methylation_hm450", "METHYLATION"),
            _profile(f"{study}_mutations", "MUTATION_EXTENDED"),
        ]

    def test_finds_preferred_rna_and_hm450(self):
        with patch("cbioportal_client.requests.get",
                   return_value=_mock_get_response(self._profiles_response())):
            rna, meth = _discover_profiles("tcga_cesc")
        assert rna == "tcga_cesc_rna_seq_v2_mrna"
        assert meth == "tcga_cesc_methylation_hm450"

    def test_falls_back_to_hm27_when_no_hm450(self):
        profiles = [
            _profile("tcga_cesc_rna_seq_v2_mrna", "MRNA_EXPRESSION"),
            _profile("tcga_cesc_methylation_hm27", "METHYLATION"),
        ]
        with patch("cbioportal_client.requests.get",
                   return_value=_mock_get_response(profiles)):
            _, meth = _discover_profiles("tcga_cesc")
        assert meth == "tcga_cesc_methylation_hm27"

    def test_falls_back_to_any_mrna_expression(self):
        profiles = [
            _profile("tcga_cesc_mrna", "MRNA_EXPRESSION"),
            _profile("tcga_cesc_methylation_hm450", "METHYLATION"),
        ]
        with patch("cbioportal_client.requests.get",
                   return_value=_mock_get_response(profiles)):
            rna, _ = _discover_profiles("tcga_cesc")
        assert rna == "tcga_cesc_mrna"

    def test_returns_none_none_on_api_error(self):
        with patch("cbioportal_client.requests.get", side_effect=Exception("connection refused")):
            rna, meth = _discover_profiles("tcga_cesc")
        assert rna is None
        assert meth is None

    def test_returns_none_when_no_methylation_profile(self):
        profiles = [_profile("tcga_cesc_rna_seq_v2_mrna", "MRNA_EXPRESSION")]
        with patch("cbioportal_client.requests.get",
                   return_value=_mock_get_response(profiles)):
            _, meth = _discover_profiles("tcga_cesc")
        assert meth is None


# ---------------------------------------------------------------------------
# _get_sample_ids
# ---------------------------------------------------------------------------

class TestGetSampleIds:
    def test_single_page(self):
        samples = [{"sampleId": f"TCGA-{i:02d}"} for i in range(10)]
        with patch("cbioportal_client.requests.get",
                   return_value=_mock_get_response(samples)):
            ids = _get_sample_ids("tcga_cesc")
        assert ids == [s["sampleId"] for s in samples]

    def test_paginates_until_partial_page(self):
        # Page 0: 500 samples, page 1: 200 samples → stops
        page0 = [{"sampleId": f"S{i}"} for i in range(500)]
        page1 = [{"sampleId": f"T{i}"} for i in range(200)]
        with patch("cbioportal_client.requests.get",
                   side_effect=[_mock_get_response(page0), _mock_get_response(page1)]):
            ids = _get_sample_ids("tcga_cesc")
        assert len(ids) == 700

    def test_empty_study_returns_empty_list(self):
        with patch("cbioportal_client.requests.get", return_value=_mock_get_response([])):
            assert _get_sample_ids("tcga_cesc") == []


# ---------------------------------------------------------------------------
# _fetch_values
# ---------------------------------------------------------------------------

class TestFetchValues:
    def test_parses_float_values(self):
        rows = [
            {"sampleId": "S1", "value": 3.14},
            {"sampleId": "S2", "value": 0.87},
        ]
        with patch("cbioportal_client.requests.post",
                   return_value=_mock_post_response(rows)):
            result = _fetch_values("profile_id", 2353, ["S1", "S2"])
        assert result == {"S1": 3.14, "S2": 0.87}

    def test_skips_null_values(self):
        rows = [{"sampleId": "S1", "value": None}, {"sampleId": "S2", "value": 1.0}]
        with patch("cbioportal_client.requests.post",
                   return_value=_mock_post_response(rows)):
            result = _fetch_values("pid", 123, ["S1", "S2"])
        assert "S1" not in result
        assert result["S2"] == 1.0

    def test_failed_batch_silently_skipped(self):
        with patch("cbioportal_client.requests.post", side_effect=Exception("timeout")):
            result = _fetch_values("pid", 123, ["S1", "S2"])
        assert result == {}


# ---------------------------------------------------------------------------
# _correlation_sync — integration of all steps
# ---------------------------------------------------------------------------

class TestCorrelationSync:
    def _setup_mocks(
        self,
        n_samples: int = 20,
        expr_vals: dict | None = None,
        meth_vals: dict | None = None,
    ):
        """Return a context that patches all HTTP calls for a clean run."""
        samples = [f"S{i}" for i in range(n_samples)]

        if expr_vals is None:
            # Negative correlation: high expression → low methylation
            expr_vals = {s: float(i) for i, s in enumerate(samples)}
        if meth_vals is None:
            meth_vals = {s: float(n_samples - i) / n_samples
                         for i, s in enumerate(samples)}

        profiles = [
            _profile("tcga_cesc_rna_seq_v2_mrna", "MRNA_EXPRESSION"),
            _profile("tcga_cesc_methylation_hm450", "METHYLATION"),
        ]

        entrez_fos = {"entrezGeneId": 2353, "hugoGeneSymbol": "FOS"}
        entrez_cdkn2a = {"entrezGeneId": 1029, "hugoGeneSymbol": "CDKN2A"}

        def _mock_get(url, *args, **kwargs):
            if "/molecular-profiles" in url:
                return _mock_get_response(profiles)
            if "/studies/" in url and "/samples" in url:
                return _mock_get_response([{"sampleId": s} for s in samples])
            if "/genes/FOS" in url:
                return _mock_get_response(entrez_fos)
            if "/genes/CDKN2A" in url:
                return _mock_get_response(entrez_cdkn2a)
            return _mock_get_response({})

        def _mock_post(url, *args, json=None, params=None, **kwargs):
            profile_id = (params or {}).get("molecularProfileId", "")
            body = json or {}
            entrez = (body.get("entrezGeneIds") or [None])[0]
            batch_ids = body.get("sampleIds", [])
            if "rna_seq" in profile_id:
                vals = expr_vals
            else:
                vals = meth_vals
            rows = [{"sampleId": s, "value": vals.get(s)} for s in batch_ids if s in vals]
            return _mock_post_response(rows)

        return _mock_get, _mock_post

    def test_successful_negative_correlation(self):
        _mock_get, _mock_post = self._setup_mocks(n_samples=20)
        with patch("cbioportal_client.requests.get", side_effect=_mock_get), \
             patch("cbioportal_client.requests.post", side_effect=_mock_post):
            result = _correlation_sync("FOS", ["CDKN2A"], "cesc")

        assert "error" not in result
        assert result["regulator"] == "FOS"
        assert result["tcga_network"] == "cesc"
        assert len(result["correlations"]) == 1
        corr = result["correlations"][0]
        assert corr["target_gene"] == "CDKN2A"
        assert corr["rho"] is not None
        assert corr["rho"] < 0  # constructed to be negative
        assert corr["p_value"] is not None
        assert "inverse" in corr["direction"]

    def test_unknown_tcga_code_returns_error(self):
        result = _correlation_sync("FOS", ["CDKN2A"], "xyz_unknown")
        assert "error" in result
        assert "Unknown TCGA code" in result["error"]

    def test_no_rna_profile_returns_error(self):
        profiles = [_profile("tcga_cesc_methylation_hm450", "METHYLATION")]
        with patch("cbioportal_client.requests.get",
                   return_value=_mock_get_response(profiles)):
            result = _correlation_sync("FOS", ["CDKN2A"], "cesc")
        assert "error" in result
        assert "RNA-seq" in result["error"]

    def test_no_methylation_profile_returns_error(self):
        profiles = [_profile("tcga_cesc_rna_seq_v2_mrna", "MRNA_EXPRESSION")]
        with patch("cbioportal_client.requests.get",
                   return_value=_mock_get_response(profiles)):
            result = _correlation_sync("FOS", ["CDKN2A"], "cesc")
        assert "error" in result
        assert "methylation" in result["error"]

    def test_unresolvable_regulator_returns_error(self):
        profiles = [
            _profile("tcga_cesc_rna_seq_v2_mrna", "MRNA_EXPRESSION"),
            _profile("tcga_cesc_methylation_hm450", "METHYLATION"),
        ]
        def _mock_get(url, *args, **kwargs):
            if "/molecular-profiles" in url:
                return _mock_get_response(profiles)
            # Gene resolution fails for everything
            return _mock_get_response({})

        with patch("cbioportal_client.requests.get", side_effect=_mock_get):
            result = _correlation_sync("FAKEGENE", ["CDKN2A"], "cesc")
        assert "error" in result
        assert "FAKEGENE" in result["error"]

    def test_insufficient_samples_flagged_per_target(self):
        _mock_get, _mock_post = self._setup_mocks(n_samples=3)  # below threshold of 5
        with patch("cbioportal_client.requests.get", side_effect=_mock_get), \
             patch("cbioportal_client.requests.post", side_effect=_mock_post):
            result = _correlation_sync("FOS", ["CDKN2A"], "cesc")

        assert "error" not in result
        corr = result["correlations"][0]
        assert corr["rho"] is None
        assert corr["direction"] == "insufficient_data"

    def test_unresolvable_target_logged_per_target(self):
        profiles = [
            _profile("tcga_cesc_rna_seq_v2_mrna", "MRNA_EXPRESSION"),
            _profile("tcga_cesc_methylation_hm450", "METHYLATION"),
        ]
        samples = [f"S{i}" for i in range(20)]
        entrez_fos = {"entrezGeneId": 2353, "hugoGeneSymbol": "FOS"}

        def _mock_get(url, *args, **kwargs):
            if "/molecular-profiles" in url:
                return _mock_get_response(profiles)
            if "/studies/" in url:
                return _mock_get_response([{"sampleId": s} for s in samples])
            if "/genes/FOS" in url:
                return _mock_get_response(entrez_fos)
            return _mock_get_response({})  # FAKEGENE resolves to {}

        def _mock_post(url, *args, json=None, params=None, **kwargs):
            body = json or {}
            batch_ids = body.get("sampleIds", [])
            rows = [{"sampleId": s, "value": float(i)}
                    for i, s in enumerate(batch_ids)]
            return _mock_post_response(rows)

        with patch("cbioportal_client.requests.get", side_effect=_mock_get), \
             patch("cbioportal_client.requests.post", side_effect=_mock_post):
            result = _correlation_sync("FOS", ["FAKEGENE"], "cesc")

        assert "error" not in result
        corr = result["correlations"][0]
        assert corr["target_gene"] == "FAKEGENE"
        assert corr["direction"] == "error"
        assert "Entrez" in corr["note"]


# ---------------------------------------------------------------------------
# format_correlation_report
# ---------------------------------------------------------------------------

class TestFormatCorrelationReport:
    def _result(self, correlations: list[dict]) -> dict:
        return {
            "regulator": "FOS",
            "tcga_network": "cesc",
            "study_id": "tcga_cesc",
            "rna_profile": "tcga_cesc_rna_seq_v2_mrna",
            "methylation_profile": "tcga_cesc_methylation_hm450",
            "n_total_samples": 297,
            "correlations": correlations,
        }

    def test_report_contains_header(self):
        report = format_correlation_report(self._result([]))
        assert "Methylation-Expression Correlation" in report
        assert "FOS" in report
        assert "CESC" in report
        assert "297" in report

    def test_negative_rho_row_formatted(self):
        corr = [{"target_gene": "CDKN2A", "n_samples": 290,
                  "rho": -0.412, "p_value": 1.2e-5,
                  "direction": "inverse (high expression → low methylation)"}]
        report = format_correlation_report(self._result(corr))
        assert "CDKN2A" in report
        assert "-0.412" in report
        assert "inverse" in report

    def test_none_rho_shows_dash(self):
        corr = [{"target_gene": "SOCS1", "n_samples": 2,
                  "rho": None, "p_value": None,
                  "direction": "insufficient_data",
                  "note": "Only 2 samples with matched data (need ≥ 5)"}]
        report = format_correlation_report(self._result(corr))
        assert "SOCS1" in report
        assert "—" in report

    def test_error_result_shows_error_message(self):
        report = format_correlation_report({"error": "Unknown TCGA code 'xyz'"})
        assert "Error" in report
        assert "xyz" in report

    def test_computational_caveat_present(self):
        report = format_correlation_report(self._result([]))
        assert "Computationally derived" in report
        assert "validation" in report
