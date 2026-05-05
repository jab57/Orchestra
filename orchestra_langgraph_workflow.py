"""
Orchestra: LangGraph workflow for MCP-over-MCP orchestration.

Issue #2 implemented: effector path (APC→CTNNB1 proof of concept)
Issue #3 implemented: TF path (TP53, BRD4→MYC — parallel RegNetAgents + CASCADE)
Issue #4 implemented: therapeutic target validation (MYC→BRD4 via super-enhancers + PPI)
"""

import asyncio
import logging
import os
from typing import Any, Optional, TypedDict

from dotenv import load_dotenv
from langgraph.graph import END, StateGraph

load_dotenv()

from mcp_client import (
    TIMEOUT_NETWORK,
    TIMEOUT_PERTURBATION,
    TIMEOUT_PPI,
    make_cascade_client,
    make_regnetagents_client,
)

logger = logging.getLogger(__name__)


class OrchestraState(TypedDict):
    # Input
    gene: str
    cell_type: str
    analysis_type: str       # "causal_chain", "therapeutic_validation", "effector_analysis"
    analysis_depth: str      # "basic", "comprehensive"

    # Gene classification (from CASCADE get_gene_metadata)
    gene_role: Optional[str]      # master_regulator, transcription_factor, effector, isolated
    ensembl_id: Optional[str]

    # Effector path: TF partner found via PPI
    tf_partner: Optional[str]

    # RegNetAgents results (via MCP)
    network_analysis: Optional[dict]
    pathway_enrichment: Optional[dict]
    domain_insights: Optional[dict]

    # CASCADE results (via MCP)
    perturbation_result: Optional[dict]
    ppi_interactions: Optional[dict]
    lincs_effects: Optional[dict]
    depmap_essentiality: Optional[dict]

    # Composite results
    validated_targets: Optional[list]
    causal_chain: Optional[dict]
    synthesis: Optional[dict]      # structured evidence table from _synthesize

    # Workflow state
    completed_steps: list
    errors: dict

    # Output
    final_report: Optional[str]


class OrchestraWorkflow:
    """
    LangGraph DAG that orchestrates RegNetAgents and CASCADE via MCP protocol calls.

    MCP clients are opened in run_analysis for the duration of the workflow.

    LLM synthesis is optional (USE_LLM_SYNTHESIS=false by default) and follows the
    same provider pattern as CASCADE and RegNetAgents:
      LLM_PROVIDER=ollama (default) | anthropic
      OLLAMA_MODEL, OLLAMA_HOST, OLLAMA_TEMPERATURE, OLLAMA_MAX_TOKENS, OLLAMA_TIMEOUT
      LLM_MODEL, LLM_API_KEY  (used when LLM_PROVIDER=anthropic)

    When disabled (default), _generate_report returns formatted structured text and
    Claude Desktop handles narrative interpretation. When enabled, a 2-3 sentence
    biological narrative is prepended — useful for standalone script output.
    """

    def __init__(self):
        self.graph = self._build_graph()
        self._cascade = None
        self._regnetagents = None

        # Optional LLM synthesis — disabled by default
        self.use_llm = os.getenv("USE_LLM_SYNTHESIS", "false").lower() == "true"
        self.llm_client: Any = None
        self.llm_available = self._initialize_llm() if self.use_llm else False
        self.ollama_model = os.getenv("OLLAMA_MODEL", "llama3.1:8b")
        self.ollama_temperature = float(os.getenv("OLLAMA_TEMPERATURE", "0.3"))
        self.ollama_max_tokens = int(os.getenv("OLLAMA_MAX_TOKENS", "2000"))

    # ------------------------------------------------------------------
    # LLM provider initialization (mirrors CASCADE / RegNetAgents pattern)
    # ------------------------------------------------------------------

    def _initialize_llm(self) -> bool:
        provider = os.getenv("LLM_PROVIDER", "ollama").lower()
        if provider == "ollama":
            return self._initialize_ollama_provider()
        elif provider == "anthropic":
            return self._initialize_anthropic_provider()
        else:
            logger.error(f"Unknown LLM_PROVIDER: {provider}. Use: ollama | anthropic")
            return False

    def _initialize_ollama_provider(self) -> bool:
        try:
            import ollama
        except ImportError:
            logger.warning("ollama package not installed. Run: pip install ollama")
            return False

        api_key = os.getenv("OLLAMA_API_KEY")
        if api_key:
            logger.info("Using Ollama Cloud (API key detected)")
            self.llm_client = ollama.Client(
                host="https://ollama.com",
                headers={"Authorization": f"Bearer {api_key}"},
            )
        else:
            host = os.getenv("OLLAMA_HOST", "http://localhost:11434")
            logger.info(f"Using local Ollama at {host}")
            self.llm_client = ollama.Client(host=host)

        try:
            models_response = self.llm_client.list()
            available = []
            for m in models_response.models if hasattr(models_response, "models") else models_response:
                name = getattr(m, "model", None) or getattr(m, "name", None) or (m.get("name") if isinstance(m, dict) else None)
                if name:
                    available.append(name)
            model_name = os.getenv("OLLAMA_MODEL", "llama3.1:8b")
            if model_name not in available:
                logger.error(f"Ollama model '{model_name}' not found. Available: {available}")
                logger.error(f"Run: ollama pull {model_name}")
                return False
            logger.info(f"Ollama available, model: {model_name}")
            return True
        except Exception as e:
            logger.warning(f"Ollama not available: {e}")
            return False

    def _initialize_anthropic_provider(self) -> bool:
        try:
            import anthropic
            api_key = os.getenv("LLM_API_KEY")
            if not api_key:
                logger.error("LLM_API_KEY required for anthropic provider")
                return False
            self.llm_client = anthropic.AsyncAnthropic(api_key=api_key)
            model = os.getenv("LLM_MODEL", "claude-haiku-4-5-20251001")
            logger.info(f"Anthropic provider initialized, model: {model}")
            return True
        except ImportError:
            logger.error("anthropic package not installed. Run: pip install anthropic")
            return False
        except Exception as e:
            logger.error(f"Anthropic provider initialization failed: {e}")
            return False

    async def _call_llm(self, prompt: str, system_prompt: Optional[str] = None) -> str:
        provider = os.getenv("LLM_PROVIDER", "ollama").lower()
        if provider == "ollama":
            return await self._call_ollama_provider(prompt, system_prompt)
        elif provider == "anthropic":
            return await self._call_anthropic_provider(prompt, system_prompt)
        raise ValueError(f"Unknown LLM_PROVIDER: {provider}")

    async def _call_ollama_provider(self, prompt: str, system_prompt: Optional[str] = None) -> str:
        timeout = int(os.getenv("OLLAMA_TIMEOUT", "60"))
        messages = []
        if system_prompt:
            messages.append({"role": "system", "content": system_prompt})
        messages.append({"role": "user", "content": prompt})

        for attempt in range(2):
            try:
                response = await asyncio.wait_for(
                    asyncio.to_thread(
                        self.llm_client.chat,
                        model=self.ollama_model,
                        messages=messages,
                        options={
                            "temperature": self.ollama_temperature,
                            "num_predict": self.ollama_max_tokens,
                        },
                    ),
                    timeout=timeout,
                )
                return response.message.content or ""
            except Exception as e:
                if attempt == 1:
                    raise
                logger.warning(f"Ollama call failed (attempt 1): {e}, retrying...")
                await asyncio.sleep(1)
        raise RuntimeError("unreachable")

    async def _call_anthropic_provider(self, prompt: str, system_prompt: Optional[str] = None) -> str:
        model = os.getenv("LLM_MODEL", "claude-haiku-4-5-20251001")
        timeout = int(os.getenv("OLLAMA_TIMEOUT", "60"))

        for attempt in range(2):
            try:
                kwargs = {
                    "model": model,
                    "max_tokens": self.ollama_max_tokens,
                    "messages": [{"role": "user", "content": prompt}],
                }
                if system_prompt:
                    kwargs["system"] = system_prompt
                response = await asyncio.wait_for(
                    self.llm_client.messages.create(**kwargs),
                    timeout=timeout,
                )
                return response.content[0].text
            except Exception as e:
                if attempt == 1:
                    raise
                logger.warning(f"Anthropic call failed (attempt 1): {e}, retrying...")
                await asyncio.sleep(1)
        raise RuntimeError("unreachable")

    # ------------------------------------------------------------------
    # Graph
    # ------------------------------------------------------------------

    def _build_graph(self) -> StateGraph:
        graph = StateGraph(OrchestraState)

        graph.add_node("initialize", self._initialize)
        graph.add_node("classify_gene", self._classify_gene)
        graph.add_node("route_analysis", self._route_analysis)
        graph.add_node("run_tf_path", self._run_tf_path)
        graph.add_node("run_effector_path", self._run_effector_path)
        graph.add_node("run_validation_path", self._run_validation_path)
        graph.add_node("synthesize", self._synthesize)
        graph.add_node("generate_report", self._generate_report)

        graph.set_entry_point("initialize")
        graph.add_edge("initialize", "classify_gene")
        graph.add_edge("classify_gene", "route_analysis")
        graph.add_conditional_edges(
            "route_analysis",
            self._routing_decision,
            {
                "tf_path": "run_tf_path",
                "effector_path": "run_effector_path",
                "validation_path": "run_validation_path",
            },
        )
        graph.add_edge("run_tf_path", "synthesize")
        graph.add_edge("run_effector_path", "synthesize")
        graph.add_edge("run_validation_path", "synthesize")
        graph.add_edge("synthesize", "generate_report")
        graph.add_edge("generate_report", END)

        return graph.compile()

    # ------------------------------------------------------------------
    # Routing
    # ------------------------------------------------------------------

    def _routing_decision(self, state: OrchestraState) -> str:
        if state.get("analysis_type") == "therapeutic_validation":
            return "validation_path"
        role = state.get("gene_role") or "isolated"
        if role in ("master_regulator", "transcription_factor", "minor_regulator"):
            return "tf_path"
        return "effector_path"

    # ------------------------------------------------------------------
    # Nodes
    # ------------------------------------------------------------------

    def _initialize(self, state: OrchestraState) -> OrchestraState:
        state["completed_steps"] = []
        state["errors"] = {}
        return state

    async def _classify_gene(self, state: OrchestraState) -> OrchestraState:
        """Call CASCADE get_gene_metadata to determine gene role and routing."""
        try:
            meta = await self._cascade.call_tool(
                "get_gene_metadata",
                {"gene": state["gene"], "cell_type": state["cell_type"]},
            )
            state["gene_role"] = meta.get("gene_type")
            state["ensembl_id"] = meta.get("ensembl_id")
        except Exception as e:
            state["errors"]["classify_gene"] = str(e)
            state["gene_role"] = "effector"  # safe default
        state["completed_steps"].append("classify_gene")
        return state

    def _route_analysis(self, state: OrchestraState) -> OrchestraState:
        state["completed_steps"].append("route_analysis")
        return state

    async def _run_tf_path(self, state: OrchestraState) -> OrchestraState:
        """
        TF / master regulator path: run RegNetAgents comprehensive network analysis and
        CASCADE perturbation in parallel, then cross-system synthesis in _synthesize.

        This is the core JOSS contribution — genes corroborated by both RegNetAgents
        network topology and CASCADE experimental data score higher than either alone.
        """
        gene = state["gene"]
        cell_type = state["cell_type"]

        rna_result, cascade_result = await asyncio.gather(
            self._regnetagents.call_tool(
                "comprehensive_gene_analysis",
                {"gene": gene, "cell_type": cell_type},
                timeout_seconds=TIMEOUT_NETWORK,
            ),
            self._cascade.call_tool(
                "comprehensive_perturbation_analysis",
                {"gene": gene, "cell_type": cell_type},
                timeout_seconds=TIMEOUT_PERTURBATION,
            ),
            return_exceptions=True,
        )

        if isinstance(rna_result, Exception):
            state["errors"]["network"] = str(rna_result)
        else:
            state["network_analysis"] = rna_result

        if isinstance(cascade_result, Exception):
            state["errors"]["perturbation"] = str(cascade_result)
        else:
            state["perturbation_result"] = cascade_result

        state["completed_steps"].append("run_tf_path")
        return state

    async def _run_effector_path(self, state: OrchestraState) -> OrchestraState:
        """
        Effector path: PPI → TF partner → parallel CASCADE perturbation + RegNetAgents analysis.
        Implements the APC→CTNNB1 use case and generalizes to other scaffold genes.
        """
        gene = state["gene"]
        cell_type = state["cell_type"]

        # Step 1: get PPI partners and find the most influential TF
        try:
            ppi = await self._cascade.call_tool(
                "get_protein_interactions",
                {"gene": gene},
                timeout_seconds=TIMEOUT_PPI,
            )
            state["ppi_interactions"] = ppi
            tf_partner = await self._find_tf_partner(ppi.get("interactions", []), cell_type)
            state["tf_partner"] = tf_partner
        except Exception as e:
            state["errors"]["ppi"] = str(e)
            state["completed_steps"].append("run_effector_path")
            return state

        if not tf_partner:
            state["errors"]["tf_partner"] = f"No TF partner found in PPI for {gene}"
            state["completed_steps"].append("run_effector_path")
            return state

        # Step 2: parallel CASCADE perturbation + RegNetAgents network analysis on TF partner
        cascade_result, rna_result = await asyncio.gather(
            self._cascade.call_tool(
                "comprehensive_perturbation_analysis",
                {"gene": tf_partner, "cell_type": cell_type},
                timeout_seconds=TIMEOUT_PERTURBATION,
            ),
            self._regnetagents.call_tool(
                "comprehensive_gene_analysis",
                {"gene": tf_partner, "cell_type": cell_type},
                timeout_seconds=TIMEOUT_NETWORK,
            ),
            return_exceptions=True,
        )

        if isinstance(cascade_result, Exception):
            state["errors"]["perturbation"] = str(cascade_result)
        else:
            state["perturbation_result"] = cascade_result

        if isinstance(rna_result, Exception):
            state["errors"]["network"] = str(rna_result)
        else:
            state["network_analysis"] = rna_result

        state["completed_steps"].append("run_effector_path")
        return state

    async def _find_tf_partner(
        self, interactions: list[dict], cell_type: str
    ) -> Optional[str]:
        """
        Find the highest-influence TF partner from PPI interactions.
        Checks top 10 partners in parallel; returns the TF with the most downstream targets.
        """
        top = sorted(interactions, key=lambda x: x.get("combined_score", 0), reverse=True)[:10]
        if not top:
            return None

        async def classify(partner: str) -> tuple[str, bool, int]:
            try:
                meta = await self._cascade.call_tool(
                    "get_gene_metadata", {"gene": partner, "cell_type": cell_type}
                )
                is_tf = meta.get("is_transcription_factor", False) or meta.get(
                    "gene_type"
                ) in ("master_regulator", "transcription_factor")
                return partner, is_tf, meta.get("num_targets", 0)
            except Exception:
                return partner, False, 0

        results = await asyncio.gather(*[classify(p["partner"]) for p in top])
        tf_candidates = [(name, n_targets) for name, is_tf, n_targets in results if is_tf]
        if not tf_candidates:
            return None
        return max(tf_candidates, key=lambda x: x[1])[0]

    async def _run_validation_path(self, state: OrchestraState) -> OrchestraState:
        """
        Therapeutic target validation path:
        1. Parallel: RegNetAgents PageRank regulators + CASCADE drug target discovery
        2. Merge unique candidates from both sources
        3. Validate top 3 via CASCADE comprehensive_perturbation_analysis (parallel)
        4. State carries validated_targets list for _synthesize_validation_path

        Two evidence layers cover complementary regulation:
        - RegNetAgents: classical TF regulators via ARACNe/GREmLN network topology
        - CASCADE: epigenetic/drug targets via super-enhancers, PPI, DoRothEA, DepMap
        Targets appearing in both layers have the highest cross-system confidence.
        """
        gene = state["gene"]
        cell_type = state["cell_type"]

        import re

        # Step 1: parallel — RegNetAgents network analysis + CASCADE drug discovery + CASCADE PPI
        # Three calls cover complementary layers: network topology, drug db, protein interactions.
        # The two CASCADE calls serialize within the same subprocess but are fast.
        rna_result, cascade_discovery, ppi_result = await asyncio.gather(
            self._regnetagents.call_tool(
                "comprehensive_gene_analysis",
                {"gene": gene, "cell_type": cell_type},
                timeout_seconds=TIMEOUT_NETWORK,
            ),
            self._cascade.call_tool(
                "therapeutic_target_discovery",
                {"gene": gene, "cell_type": cell_type},
                timeout_seconds=TIMEOUT_NETWORK,
            ),
            self._cascade.call_tool(
                "get_protein_interactions",
                {"gene": gene},
                timeout_seconds=TIMEOUT_PPI,
            ),
            return_exceptions=True,
        )

        if isinstance(rna_result, Exception):
            state["errors"]["network"] = str(rna_result)
            rna_result = None
        else:
            state["network_analysis"] = rna_result

        if isinstance(cascade_discovery, Exception):
            state["errors"]["cascade_discovery"] = str(cascade_discovery)
            cascade_discovery = None

        if isinstance(ppi_result, Exception):
            state["errors"]["ppi"] = str(ppi_result)
            ppi_result = None

        # Step 2: extract candidate target names — three ordered sources
        candidates: list[dict] = []
        seen: set[str] = set()

        def _add_candidate(name: str, info: dict) -> None:
            if name and name not in seen:
                seen.add(name)
                candidates.append({"gene": name, **info})

        # Source A: RegNetAgents ranked upstream regulators by PageRank
        if rna_result:
            ttp = (rna_result or {}).get("therapeutic_target_prioritization") or {}
            for r in (ttp.get("ranked_regulators") or [])[:5]:
                name = r.get("regulator")
                _add_candidate(name, {
                    "source": "regnetagents_pagerank",
                    "pagerank": (r.get("centrality_metrics") or {}).get("pagerank", 0),
                    "downstream_targets": r.get("regulator_downstream_targets", 0),
                })

        # Source B: CASCADE drug suggestions — explicit target field OR gene names in action text
        _ACTION_GENE_RE = re.compile(r'\bConsider\s+([A-Z][A-Z0-9]{1,5})\b')
        if cascade_discovery:
            for suggestion in (cascade_discovery or {}).get("therapeutic_targets") or []:
                target = suggestion.get("target")
                if target and target != "unknown":
                    _add_candidate(target, {
                        "source": "cascade_drug_discovery",
                        "priority": suggestion.get("priority", ""),
                        "reason": suggestion.get("reason", ""),
                    })
                else:
                    # Extract gene name from action text (e.g. "Consider BRD4/BET inhibitors")
                    action = suggestion.get("action", "")
                    m = _ACTION_GENE_RE.search(action)
                    if m:
                        _add_candidate(m.group(1), {
                            "source": "cascade_drug_discovery",
                            "priority": suggestion.get("priority", ""),
                            "reason": suggestion.get("reason", ""),
                        })

        # Source C: STRING PPI top partners (protein-level interactors as drug targets)
        if ppi_result:
            interactions = sorted(
                (ppi_result or {}).get("interactions", []),
                key=lambda x: x.get("combined_score", 0),
                reverse=True,
            )[:3]
            for i in interactions:
                partner = i.get("partner")
                _add_candidate(partner, {
                    "source": "cascade_ppi",
                    "ppi_score": i.get("combined_score", 0),
                })

        # Step 3: validate top 3 unique candidates via CASCADE comprehensive perturbation
        top = candidates[:3]
        if top:
            validation_results = await asyncio.gather(
                *[
                    self._cascade.call_tool(
                        "comprehensive_perturbation_analysis",
                        {"gene": c["gene"], "cell_type": cell_type},
                        timeout_seconds=TIMEOUT_PERTURBATION,
                    )
                    for c in top
                ],
                return_exceptions=True,
            )
            for candidate, result in zip(top, validation_results):
                if isinstance(result, Exception):
                    candidate["cascade_error"] = str(result)
                else:
                    ev = (result or {}).get("evidence_synthesis") or {}
                    candidate["key_findings"] = ev.get("key_findings", [])
                    candidate["multi_source_genes"] = [
                        g["symbol"] for g in (ev.get("multi_source_genes") or [])[:10]
                    ]

        # Fill remaining candidates (not validated) with empty evidence
        for candidate in candidates[3:]:
            candidate.setdefault("key_findings", [])
            candidate.setdefault("multi_source_genes", [])

        state["validated_targets"] = candidates
        state["completed_steps"].append("run_validation_path")
        return state

    # ------------------------------------------------------------------
    # Synthesis
    # ------------------------------------------------------------------

    def _synthesize(self, state: OrchestraState) -> OrchestraState:
        if state.get("analysis_type") == "therapeutic_validation":
            return self._synthesize_validation_path(state)
        role = state.get("gene_role") or "effector"
        if role in ("master_regulator", "transcription_factor", "minor_regulator"):
            return self._synthesize_tf_path(state)
        return self._synthesize_effector_path(state)

    def _synthesize_tf_path(self, state: OrchestraState) -> OrchestraState:
        """
        TF path synthesis: cross-system corroboration between RegNetAgents network
        targets and CASCADE experimentally validated genes.

        A gene that appears in both RegNetAgents downstream targets (network topology)
        AND CASCADE multi_source_genes (experimental data) is a cross-system hit —
        the core evidence Orchestra produces that neither system alone can generate.
        """
        perturbation = state.get("perturbation_result") or {}
        network = state.get("network_analysis") or {}

        evidence_synthesis = perturbation.get("evidence_synthesis", {})
        key_findings = evidence_synthesis.get("key_findings", [])
        multi_source_genes = evidence_synthesis.get("multi_source_genes", [])
        source_agreements = evidence_synthesis.get("source_agreements", [])
        source_disagreements = evidence_synthesis.get("source_disagreements", [])

        # Extract RegNetAgents downstream targets for cross-system overlap
        rna_targets = self._extract_regnetagents_targets(network)

        # Cross-system hits: genes supported by both RegNetAgents network topology
        # AND CASCADE experimental evidence (LINCS, DepMap, etc.)
        cross_system_hits = [
            {
                "symbol": g["symbol"],
                "source_count": g["source_count"],
                "sources": g["sources"],
            }
            for g in multi_source_genes
            if g["symbol"] in rna_targets
        ]

        network_summary = (
            network.get("summary")
            or network.get("network_analysis")
            or network.get("workflow_summary")
            or {}
        )

        state["synthesis"] = {
            "gene": state["gene"],
            "cell_type": state["cell_type"],
            "routing": "tf",
            "gene_role": state.get("gene_role"),
            "cascade_key_findings": key_findings,
            "corroborated_targets": [
                {
                    "symbol": g["symbol"],
                    "source_count": g["source_count"],
                    "sources": g["sources"],
                }
                for g in multi_source_genes[:10]
            ],
            "cross_system_hits": cross_system_hits,
            "source_agreements": source_agreements,
            "source_disagreements": source_disagreements,
            "network_context": network_summary,
            "regnetagents_target_count": len(rna_targets),
            "errors": state.get("errors", {}),
        }
        state["completed_steps"].append("synthesize")
        return state

    def _synthesize_effector_path(self, state: OrchestraState) -> OrchestraState:
        """
        Effector path synthesis: consumes CASCADE's pre-synthesized evidence block
        for the TF partner, plus RegNetAgents network context for pathway framing.
        """
        perturbation = state.get("perturbation_result") or {}
        network = state.get("network_analysis") or {}

        evidence_synthesis = perturbation.get("evidence_synthesis", {})
        key_findings = evidence_synthesis.get("key_findings", [])
        multi_source_genes = evidence_synthesis.get("multi_source_genes", [])
        source_agreements = evidence_synthesis.get("source_agreements", [])
        source_disagreements = evidence_synthesis.get("source_disagreements", [])

        network_summary = (
            network.get("summary")
            or network.get("network_analysis")
            or network.get("workflow_summary")
            or {}
        )

        state["synthesis"] = {
            "gene": state["gene"],
            "cell_type": state["cell_type"],
            "routing": "effector",
            "tf_partner": state.get("tf_partner"),
            "cascade_key_findings": key_findings,
            "corroborated_targets": [
                {
                    "symbol": g["symbol"],
                    "source_count": g["source_count"],
                    "sources": g["sources"],
                }
                for g in multi_source_genes[:10]
            ],
            "source_agreements": source_agreements,
            "source_disagreements": source_disagreements,
            "network_context": network_summary,
            "errors": state.get("errors", {}),
        }
        state["completed_steps"].append("synthesize")
        return state

    def _extract_regnetagents_targets(self, network: dict) -> set:
        """
        Extract downstream target gene symbols from RegNetAgents comprehensive_gene_analysis output.

        Primary path: target_analysis.cascade_targets[].gene_symbol
        Fallback: generic extraction from common field names across tool versions.
        """
        targets = set()

        def _collect(container):
            if isinstance(container, list):
                for item in container:
                    if isinstance(item, str):
                        targets.add(item)
                    elif isinstance(item, dict):
                        sym = (
                            item.get("gene_symbol")
                            or item.get("symbol")
                            or item.get("gene")
                            or item.get("name")
                        )
                        if sym:
                            targets.add(sym)

        # Primary: target_analysis.cascade_targets (confirmed structure from RegNetAgents output)
        target_analysis = network.get("target_analysis") or {}
        _collect(target_analysis.get("cascade_targets"))

        # Fallback: generic field names at top level and in common sub-blocks
        for field in ("targets", "downstream_targets", "regulated_genes", "top_targets", "network_targets"):
            _collect(network.get(field))

        for block_key in ("network_analysis", "summary", "workflow_summary"):
            block = network.get(block_key)
            if isinstance(block, dict):
                for field in ("targets", "downstream_targets", "regulated_genes", "top_targets", "cascade_targets"):
                    _collect(block.get(field))

        return targets

    def _synthesize_validation_path(self, state: OrchestraState) -> OrchestraState:
        """
        Validation path synthesis: builds per-candidate evidence table from
        RegNetAgents network rank + CASCADE perturbation confirmation.
        """
        validated_targets = state.get("validated_targets") or []
        network = state.get("network_analysis") or {}
        network_summary = (
            network.get("summary")
            or network.get("network_analysis")
            or network.get("workflow_summary")
            or {}
        )

        state["synthesis"] = {
            "gene": state["gene"],
            "cell_type": state["cell_type"],
            "routing": "validation",
            "validated_targets": validated_targets,
            "network_context": network_summary,
            "errors": state.get("errors", {}),
        }
        state["completed_steps"].append("synthesize")
        return state

    # ------------------------------------------------------------------
    # Report generation
    # ------------------------------------------------------------------

    async def _generate_report(self, state: OrchestraState) -> OrchestraState:
        """
        Format the structured synthesis as a report.

        Always returns well-formatted structured text — Claude Desktop handles
        the narrative interpretation when running as an MCP server.

        When USE_LLM_SYNTHESIS=true, prepends a 2-3 sentence biological narrative
        generated by the configured LLM provider (LLM_PROVIDER=ollama|anthropic).
        Useful for standalone script runs where no host LLM is present.
        """
        synthesis = state.get("synthesis") or {}

        report_sections = self._format_evidence_report(synthesis)

        routing = synthesis.get("routing", "effector")
        llm_trigger = (
            bool(synthesis.get("tf_partner")) if routing == "effector"
            else routing in ("tf", "validation")
        )

        if self.llm_available and llm_trigger:
            try:
                narrative = await self._call_llm_synthesis(synthesis)
                report_sections.insert(0, narrative + "\n")
            except Exception as e:
                logger.warning(f"LLM synthesis failed, using structured report: {e}")

        state["final_report"] = "\n".join(report_sections)
        state["completed_steps"].append("generate_report")
        return state

    def _format_evidence_report(self, synthesis: dict) -> list[str]:
        routing = synthesis.get("routing", "effector")
        if routing == "tf":
            return self._format_tf_report(synthesis)
        if routing == "validation":
            return self._format_validation_report(synthesis)
        return self._format_effector_report(synthesis)

    def _format_tf_report(self, synthesis: dict) -> list[str]:
        gene = synthesis.get("gene", "unknown")
        cell_type = synthesis.get("cell_type", "unknown")
        gene_role = synthesis.get("gene_role", "transcription_factor")
        key_findings = synthesis.get("cascade_key_findings", [])
        corroborated = synthesis.get("corroborated_targets", [])
        cross_system_hits = synthesis.get("cross_system_hits", [])
        agreements = synthesis.get("source_agreements", [])
        network_ctx = synthesis.get("network_context", {})
        rna_target_count = synthesis.get("regnetagents_target_count", 0)
        errors = synthesis.get("errors", {})

        lines = [
            f"## Orchestra Analysis: {gene} in {cell_type}",
            f"**Routing:** TF / {gene_role}",
            "",
            "### CASCADE Evidence",
        ]

        if key_findings:
            for f in key_findings:
                lines.append(f"- {f}")
        else:
            lines.append("- No key findings available")

        if corroborated:
            lines.append("")
            lines.append("**CASCADE multi-source corroborated downstream genes:**")
            for g in corroborated[:8]:
                lines.append(
                    f"- {g['symbol']}: {g['source_count']} sources "
                    f"({', '.join(g['sources'])})"
                )

        if agreements:
            lines.append("")
            lines.append("**Cross-source agreements (within CASCADE):**")
            for a in agreements[:5]:
                lines.append(f"- {a}")

        lines.append("")
        lines.append("### Cross-System Corroboration")
        lines.append(
            f"RegNetAgents downstream targets: {rna_target_count} genes"
        )

        if cross_system_hits:
            lines.append(
                f"**Cross-system hits** ({len(cross_system_hits)} genes in both "
                "RegNetAgents network topology AND CASCADE experimental data):"
            )
            for g in cross_system_hits[:10]:
                lines.append(
                    f"- {g['symbol']}: CASCADE {g['source_count']} sources "
                    f"({', '.join(g['sources'])}) + RegNetAgents network target"
                )
        else:
            lines.append(
                "No cross-system hits detected — CASCADE and RegNetAgents targets "
                "do not overlap for this gene/cell-type combination."
            )
            if rna_target_count == 0:
                lines.append(
                    "  (RegNetAgents returned 0 extractable targets — "
                    "check network_analysis field structure)"
                )

        if network_ctx:
            lines.append("")
            lines.append("### RegNetAgents Network Context")
            lines.append(str(network_ctx)[:600])

        if errors:
            lines.append("")
            lines.append("### Partial Data Warnings")
            for k, v in errors.items():
                lines.append(f"- {k}: {v}")

        return lines

    def _format_effector_report(self, synthesis: dict) -> list[str]:
        gene = synthesis.get("gene", "unknown")
        cell_type = synthesis.get("cell_type", "unknown")
        tf_partner = synthesis.get("tf_partner") or "not identified"
        key_findings = synthesis.get("cascade_key_findings", [])
        corroborated = synthesis.get("corroborated_targets", [])
        agreements = synthesis.get("source_agreements", [])
        network_ctx = synthesis.get("network_context", {})
        errors = synthesis.get("errors", {})

        lines = [
            f"## Orchestra Analysis: {gene} in {cell_type}",
            f"**Routing:** effector/scaffold",
            f"**TF partner (via PPI):** {tf_partner}",
            "",
            "### CASCADE Evidence",
        ]

        if key_findings:
            for f in key_findings:
                lines.append(f"- {f}")
        else:
            lines.append("- No key findings available")

        if corroborated:
            lines.append("")
            lines.append("**Multi-source corroborated downstream genes:**")
            for g in corroborated[:8]:
                lines.append(
                    f"- {g['symbol']}: {g['source_count']} sources "
                    f"({', '.join(g['sources'])})"
                )

        if agreements:
            lines.append("")
            lines.append("**Cross-source agreements:**")
            for a in agreements[:5]:
                lines.append(f"- {a}")

        if network_ctx:
            lines.append("")
            lines.append("### RegNetAgents Network Context")
            lines.append(str(network_ctx)[:600])

        if errors:
            lines.append("")
            lines.append("### Partial Data Warnings")
            for k, v in errors.items():
                lines.append(f"- {k}: {v}")

        return lines

    # ------------------------------------------------------------------
    # LLM synthesis
    # ------------------------------------------------------------------

    def _format_validation_report(self, synthesis: dict) -> list[str]:
        gene = synthesis.get("gene", "unknown")
        cell_type = synthesis.get("cell_type", "unknown")
        candidates = synthesis.get("validated_targets") or []
        errors = synthesis.get("errors", {})

        lines = [
            f"## Orchestra Therapeutic Target Validation: {gene} in {cell_type}",
            f"**Routing:** therapeutic_validation",
            f"**Candidates evaluated:** {len(candidates)}",
            "",
        ]

        if not candidates:
            lines.append("No therapeutic target candidates identified.")
            if errors:
                lines.append("")
                lines.append("### Errors")
                for k, v in errors.items():
                    lines.append(f"- {k}: {v}")
            return lines

        for i, c in enumerate(candidates, 1):
            name = c.get("gene", "unknown")
            source = c.get("source", "unknown")
            lines.append(f"### Candidate {i}: {name}")

            if source == "regnetagents_pagerank":
                lines.append(
                    f"**Source:** RegNetAgents PageRank  |  "
                    f"pagerank={c.get('pagerank', 0):.4f}  |  "
                    f"downstream_targets={c.get('downstream_targets', 0)}"
                )
            elif source == "cascade_drug_discovery":
                lines.append(
                    f"**Source:** CASCADE drug discovery  |  "
                    f"priority={c.get('priority', '')}  |  "
                    f"reason: {c.get('reason', '')}"
                )
            else:
                lines.append(f"**Source:** {source}")

            findings = c.get("key_findings") or []
            if findings:
                lines.append("**CASCADE perturbation evidence:**")
                for f in findings:
                    lines.append(f"  - {f}")
            elif "cascade_error" in c:
                lines.append(f"**CASCADE validation error:** {c['cascade_error']}")
            else:
                lines.append("**CASCADE validation:** not run (candidate beyond top 3)")

            downstream = c.get("multi_source_genes") or []
            if downstream:
                lines.append(f"**Top downstream genes (CASCADE):** {', '.join(downstream[:8])}")

            lines.append("")

        if errors:
            lines.append("### Partial Data Warnings")
            for k, v in errors.items():
                lines.append(f"- {k}: {v}")

        return lines

    async def _call_llm_synthesis(self, synthesis: dict) -> str:
        routing = synthesis.get("routing", "effector")
        if routing == "tf":
            return await self._call_llm_tf_synthesis(synthesis)
        if routing == "validation":
            return await self._call_llm_validation_synthesis(synthesis)
        return await self._call_llm_effector_synthesis(synthesis)

    async def _call_llm_tf_synthesis(self, synthesis: dict) -> str:
        gene = synthesis.get("gene", "unknown")
        cell_type = synthesis.get("cell_type", "unknown")
        gene_role = synthesis.get("gene_role", "transcription_factor")
        key_findings = synthesis.get("cascade_key_findings", [])
        cross_system_hits = synthesis.get("cross_system_hits", [])[:5]
        corroborated = synthesis.get("corroborated_targets", [])[:5]

        cross_str = (
            "\n".join(
                f"- {g['symbol']}: CASCADE {g['source_count']} sources + RegNetAgents target"
                for g in cross_system_hits
            )
            if cross_system_hits
            else "No cross-system overlap detected."
        )

        prompt = (
            f"Gene: {gene} ({gene_role}) in {cell_type}\n\n"
            "CASCADE key findings:\n"
            + "\n".join(f"- {f}" for f in key_findings)
            + "\n\nCross-system corroborated targets (in both network topology and experimental data):\n"
            + cross_str
            + "\n\nTop CASCADE-only corroborated genes:\n"
            + "\n".join(f"- {g['symbol']}: {g['source_count']} sources" for g in corroborated)
            + f"\n\nWrite a concise 2-3 sentence biological interpretation of {gene}'s "
            f"regulatory role and what the cross-system evidence tells us about its "
            "downstream targets and therapeutic relevance. "
            "Be specific about pathways and what the cross-system agreement implies."
        )
        system_prompt = (
            "You are Orchestra, a bioinformatics analysis system. "
            "Narrate structured evidence concisely. Do not speculate beyond the data provided."
        )
        return await self._call_llm(prompt, system_prompt)

    async def _call_llm_validation_synthesis(self, synthesis: dict) -> str:
        gene = synthesis.get("gene", "unknown")
        cell_type = synthesis.get("cell_type", "unknown")
        candidates = synthesis.get("validated_targets") or []

        candidate_lines = []
        for c in candidates[:5]:
            name = c.get("gene", "unknown")
            source = c.get("source", "")
            findings = c.get("key_findings") or []
            findings_str = "; ".join(findings[:2]) if findings else "not validated"
            candidate_lines.append(f"- {name} [{source}]: {findings_str}")

        prompt = (
            f"Therapeutic target analysis for {gene} in {cell_type}.\n\n"
            "Candidate therapeutic targets identified:\n"
            + "\n".join(candidate_lines)
            + f"\n\nWrite a concise 2-3 sentence summary of the top therapeutic targets "
            f"for inhibiting {gene}, which candidates have the strongest evidence, "
            "and what the CASCADE perturbation data says about each. "
            "Be specific about evidence sources and therapeutic relevance."
        )
        system_prompt = (
            "You are Orchestra, a bioinformatics analysis system. "
            "Narrate structured evidence concisely. Do not speculate beyond the data provided."
        )
        return await self._call_llm(prompt, system_prompt)

    async def _call_llm_effector_synthesis(self, synthesis: dict) -> str:
        gene = synthesis.get("gene", "unknown")
        cell_type = synthesis.get("cell_type", "unknown")
        tf_partner = synthesis.get("tf_partner", "unknown")
        key_findings = synthesis.get("cascade_key_findings", [])
        corroborated = synthesis.get("corroborated_targets", [])[:5]

        prompt = (
            f"Gene queried: {gene} in {cell_type}\n"
            f"Analysis: effector/scaffold — TF partner: {tf_partner}\n\n"
            "CASCADE key findings:\n"
            + "\n".join(f"- {f}" for f in key_findings)
            + "\n\nTop corroborated downstream genes:\n"
            + "\n".join(f"- {g['symbol']}: {g['source_count']} sources" for g in corroborated)
            + f"\n\nWrite a concise 2-3 sentence biological interpretation of how {gene} "
            f"acts through {tf_partner} and what the downstream effects mean. "
            "Be specific about pathways and therapeutic relevance."
        )
        system_prompt = (
            "You are Orchestra, a bioinformatics analysis system. "
            "Narrate structured evidence concisely. Do not speculate beyond the data provided."
        )
        return await self._call_llm(prompt, system_prompt)

    # ------------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------------

    async def run_analysis(
        self,
        gene: str,
        cell_type: str,
        analysis_type: str = "causal_chain",
        analysis_depth: str = "comprehensive",
    ) -> dict:
        """Run a full Orchestra analysis. Opens MCP client connections for the duration."""
        initial_state = OrchestraState(
            gene=gene,
            cell_type=cell_type,
            analysis_type=analysis_type,
            analysis_depth=analysis_depth,
            gene_role=None,
            ensembl_id=None,
            tf_partner=None,
            network_analysis=None,
            pathway_enrichment=None,
            domain_insights=None,
            perturbation_result=None,
            ppi_interactions=None,
            lincs_effects=None,
            depmap_essentiality=None,
            validated_targets=None,
            causal_chain=None,
            synthesis=None,
            completed_steps=[],
            errors={},
            final_report=None,
        )

        async with make_cascade_client() as cascade, make_regnetagents_client() as regnetagents:
            self._cascade = cascade
            self._regnetagents = regnetagents
            try:
                result = await self.graph.ainvoke(initial_state)
            finally:
                self._cascade = None
                self._regnetagents = None

        return result
