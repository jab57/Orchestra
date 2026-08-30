"""
Tests for tcga_network scoping in causal_chain_analysis (TF + effector paths) and
validate_therapeutic_targets. Both previously only ever queried the population-averaged
GREmLN network; cancer_type (already an existing OrchestraState field, reused rather
than adding a new one) now threads through to comprehensive_gene_analysis (RegNetAgents)
and comprehensive_perturbation_analysis (CASCADE), both of which already support TCGA
scoping at the child-server level.

Confirms both directions:
  - cancer_type set -> TCGA-scoped params passed to the RegNetAgents/CASCADE calls
    that support it, and NOT passed to the calls that don't (get_protein_interactions,
    therapeutic_target_discovery) -- those tools have no tcga_network/network_source
    parameter in their schemas at all.
  - cancer_type unset -> behavior is byte-identical to before this change (regression
    guard for the default/existing GREmLN-only usage).

get_gene_metadata gained network_source/tcga_network support on the CASCADE side
2026-08-30 (see FUTURE_ROADMAP.md, "Future work: TCGA-aware gene-role classification"
-- now implemented, not future work) -- _classify_gene now threads cancer_type through
to it too, so routing itself (not just the downstream analysis calls) honors the
requested network.
"""
import pytest
from orchestra_langgraph_workflow import OrchestraWorkflow


def _make_workflow() -> OrchestraWorkflow:
    wf = object.__new__(OrchestraWorkflow)
    wf.use_llm = False
    wf.llm_available = False
    wf.llm_client = None
    wf.graph = None
    return wf


def _tf_state(gene="TP53", cell_type="epithelial_cell", cancer_type=None) -> dict:
    return {
        "gene": gene, "cell_type": cell_type, "cancer_type": cancer_type,
        "completed_steps": [], "errors": {},
    }


class _RecordingClient:
    """Records every call_tool(name, params) invocation; returns canned responses."""

    def __init__(self, responses: dict):
        self.calls: list[tuple[str, dict]] = []
        self._responses = responses

    async def call_tool(self, name, params, timeout_seconds=None):
        self.calls.append((name, dict(params)))
        return self._responses.get(name, {})


class TestParamHelpers:
    def test_rna_params_no_cancer_type(self):
        params = OrchestraWorkflow._rna_gene_analysis_params("TP53", "epithelial_cell", None)
        assert params == {"gene": "TP53", "cell_type": "epithelial_cell"}
        assert "tcga_network" not in params

    def test_rna_params_with_cancer_type(self):
        params = OrchestraWorkflow._rna_gene_analysis_params("TP53", "epithelial_cell", "cesc")
        assert params == {"gene": "TP53", "cell_type": "epithelial_cell", "tcga_network": "cesc"}

    def test_cascade_params_no_cancer_type(self):
        params = OrchestraWorkflow._cascade_perturbation_params("TP53", "epithelial_cell", None)
        assert params == {"gene": "TP53", "cell_type": "epithelial_cell"}
        assert "network_source" not in params

    def test_cascade_params_with_cancer_type(self):
        params = OrchestraWorkflow._cascade_perturbation_params("TP53", "epithelial_cell", "cesc")
        assert params == {
            "gene": "TP53", "cell_type": "epithelial_cell",
            "network_source": "tcga", "tcga_network": "cesc",
        }


class TestClassifyGeneTcgaScoping:
    """_classify_gene threads cancer_type through to get_gene_metadata now that CASCADE
    supports network_source/tcga_network on that tool (2026-08-30)."""

    @pytest.mark.asyncio
    async def test_cancer_type_set_scopes_classification_call(self):
        wf = _make_workflow()
        cascade = _RecordingClient({
            "get_gene_metadata": {"gene_type": "master_regulator", "ensembl_id": "ENSG00000180264"},
        })
        wf._cascade = cascade
        state = _tf_state(gene="ADGRD2", cancer_type="cesc")

        result = await wf._classify_gene(state)

        assert result["gene_role"] == "master_regulator"
        assert result["ensembl_id"] == "ENSG00000180264"
        call_params = cascade.calls[0][1]
        assert call_params["network_source"] == "tcga"
        assert call_params["tcga_network"] == "cesc"

    @pytest.mark.asyncio
    async def test_no_cancer_type_is_unchanged(self):
        """Regression guard: default GREmLN-only classification call is byte-identical
        to before this change -- no network_source/tcga_network keys at all."""
        wf = _make_workflow()
        cascade = _RecordingClient({
            "get_gene_metadata": {"gene_type": "effector", "ensembl_id": "X"},
        })
        wf._cascade = cascade
        state = _tf_state(gene="ADGRD2")  # cancer_type=None

        result = await wf._classify_gene(state)

        assert result["gene_role"] == "effector"
        call_params = cascade.calls[0][1]
        assert call_params == {"gene": "ADGRD2", "cell_type": "epithelial_cell"}
        assert "network_source" not in call_params
        assert "tcga_network" not in call_params

    @pytest.mark.asyncio
    async def test_skipped_for_analysis_types_with_no_single_gene(self):
        """gene_signature etc. skip classification entirely -- cancer_type being set
        must not change that."""
        wf = _make_workflow()
        cascade = _RecordingClient({})
        wf._cascade = cascade
        state = _tf_state(cancer_type="cesc")
        state["analysis_type"] = "gene_signature"

        await wf._classify_gene(state)

        assert cascade.calls == []


class TestRunTfPathTcgaScoping:
    @pytest.mark.asyncio
    async def test_cancer_type_threads_to_all_three_calls(self):
        wf = _make_workflow()
        client = _RecordingClient({
            "comprehensive_gene_analysis": {},
            "comprehensive_perturbation_analysis": {},
            "query_network": {"targets": []},
        })
        wf._regnetagents = client
        wf._cascade = client
        state = _tf_state(cancer_type="cesc")

        result = await wf._run_tf_path(state)

        by_name = {name: params for name, params in client.calls}
        assert by_name["comprehensive_gene_analysis"]["tcga_network"] == "cesc"
        assert by_name["comprehensive_perturbation_analysis"]["network_source"] == "tcga"
        assert by_name["comprehensive_perturbation_analysis"]["tcga_network"] == "cesc"
        assert by_name["query_network"]["network_source"] == "tcga"
        assert by_name["query_network"]["tcga_network"] == "cesc"

    @pytest.mark.asyncio
    async def test_no_cancer_type_is_unchanged(self):
        wf = _make_workflow()
        client = _RecordingClient({
            "comprehensive_gene_analysis": {},
            "comprehensive_perturbation_analysis": {},
            "query_network": {"targets": []},
        })
        wf._regnetagents = client
        wf._cascade = client
        state = _tf_state()  # cancer_type=None

        await wf._run_tf_path(state)

        for name, params in client.calls:
            assert "tcga_network" not in params, f"{name} unexpectedly got tcga_network"
            assert "network_source" not in params, f"{name} unexpectedly got network_source"


class TestRunEffectorPathTcgaScoping:
    @pytest.mark.asyncio
    async def test_step2_scoped_step1_not(self):
        """PPI discovery (Step 1) has no TCGA param at all; only Step 2's TF-partner
        analysis (once found) should carry cancer_type."""
        wf = _make_workflow()
        client = _RecordingClient({
            "get_protein_interactions": {"interactions": [{"partner": "CTNNB1", "combined_score": 0.99}]},
            "get_gene_metadata": {"is_transcription_factor": True, "gene_type": "transcription_factor", "num_targets": 50},
            "comprehensive_perturbation_analysis": {},
            "comprehensive_gene_analysis": {},
        })
        wf._cascade = client
        wf._regnetagents = client
        state = _tf_state(gene="APC", cancer_type="cesc")

        await wf._run_effector_path(state)

        by_name = {}
        for name, params in client.calls:
            by_name.setdefault(name, []).append(params)

        # Step 1: PPI call never receives cancer_type-derived params (no such param exists)
        assert "tcga_network" not in by_name["get_protein_interactions"][0]
        # get_gene_metadata (TF-partner classification) also has no TCGA param in its schema
        for params in by_name.get("get_gene_metadata", []):
            assert "tcga_network" not in params

        # Step 2: TF partner's own analysis IS scoped
        assert by_name["comprehensive_perturbation_analysis"][0]["network_source"] == "tcga"
        assert by_name["comprehensive_perturbation_analysis"][0]["tcga_network"] == "cesc"
        assert by_name["comprehensive_gene_analysis"][0]["tcga_network"] == "cesc"


class TestRunValidationPathTcgaScoping:
    @pytest.mark.asyncio
    async def test_scoped_and_unscoped_calls_split_correctly(self):
        wf = _make_workflow()
        client = _RecordingClient({
            "comprehensive_gene_analysis": {"therapeutic_target_prioritization": {"ranked_regulators": [
                {"regulator": "MYC", "centrality_metrics": {"pagerank": 0.1}, "regulator_downstream_targets": 10}
            ]}},
            "therapeutic_target_discovery": {"therapeutic_targets": []},
            "get_protein_interactions": {"interactions": []},
            "comprehensive_perturbation_analysis": {},
        })
        wf._regnetagents = client
        wf._cascade = client
        state = _tf_state(gene="MYC", cancer_type="cesc")

        await wf._run_validation_path(state)

        by_name = {}
        for name, params in client.calls:
            by_name.setdefault(name, []).append(params)

        # Scoped: RegNetAgents network analysis + per-candidate CASCADE validation
        assert by_name["comprehensive_gene_analysis"][0]["tcga_network"] == "cesc"
        assert by_name["comprehensive_perturbation_analysis"][0]["network_source"] == "tcga"
        assert by_name["comprehensive_perturbation_analysis"][0]["tcga_network"] == "cesc"

        # Not scoped: no TCGA param exists in these tools' schemas at all
        assert "tcga_network" not in by_name["therapeutic_target_discovery"][0]
        assert "tcga_network" not in by_name["get_protein_interactions"][0]


class TestReportShowsNetworkLabel:
    """Report headers surface which network was actually used, so a TCGA-scoped
    request is visibly confirmed rather than silently identical-looking to GREmLN."""

    def _wf(self):
        wf = object.__new__(OrchestraWorkflow)
        return wf

    def test_tf_report_shows_tcga_network(self):
        wf = self._wf()
        synthesis = {"gene": "TP53", "cell_type": "epithelial_cell", "cancer_type": "cesc",
                     "gene_role": "transcription_factor", "regnetagents_available": True,
                     "cascade_available": True, "errors": {}}
        lines = wf._format_tf_report(synthesis)
        assert "**Network:** TCGA CESC" in lines

    def test_tf_report_shows_gremln_when_no_cancer_type(self):
        wf = self._wf()
        synthesis = {"gene": "TP53", "cell_type": "epithelial_cell", "cancer_type": None,
                     "gene_role": "transcription_factor", "regnetagents_available": True,
                     "cascade_available": True, "errors": {}}
        lines = wf._format_tf_report(synthesis)
        assert "**Network:** GREmLN epithelial_cell" in lines

    def test_effector_report_shows_tcga_network(self):
        wf = self._wf()
        synthesis = {"gene": "APC", "cell_type": "epithelial_cell", "cancer_type": "cesc",
                     "regnetagents_available": True, "cascade_available": True, "errors": {}}
        lines = wf._format_effector_report(synthesis)
        assert "**Network:** TCGA CESC" in lines

    def test_validation_report_shows_tcga_network(self):
        wf = self._wf()
        synthesis = {"gene": "MYC", "cell_type": "epithelial_cell", "cancer_type": "brca",
                     "validated_targets": [], "regnetagents_available": True,
                     "cascade_available": True, "errors": {}}
        lines = wf._format_validation_report(synthesis)
        assert "**Network:** TCGA BRCA" in lines
