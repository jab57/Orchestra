"""
Orchestra: LangGraph workflow for MCP-over-MCP orchestration.

Issue #2 implemented: effector path (APC→CTNNB1 proof of concept)
Issue #3 implemented: TF path (TP53, BRD4→MYC — parallel RegNetAgents + CASCADE)
Issue #4 implemented: therapeutic target validation (MYC→BRD4 via super-enhancers + PPI)
Issue #9 implemented: cross-system discordance reporting (discordance_flags in synthesis layer)
Issue #10 implemented: gene signature driver analysis (analyze_gene_signature — DEG list → ranked TF drivers)
Issue #11 implemented: cross-cell-type context comparison (compare_cell_contexts — 7-source heatmap across N cell types)
"""

import asyncio
import logging
import math
import os
from collections.abc import Callable, Coroutine
from contextvars import ContextVar
from typing import Any, Optional, TypedDict

from dotenv import load_dotenv
from langgraph.graph import END, StateGraph

load_dotenv()

from mcp_client import (
    TIMEOUT_MASTER_REGULATORS,
    TIMEOUT_NETWORK,
    TIMEOUT_NETWORK_COMPARISON,
    TIMEOUT_PERTURBATION,
    TIMEOUT_PPI,
    make_cascade_client,
    make_regnetagents_client,
)

logger = logging.getLogger(__name__)

# Task-local progress callback — inherited by asyncio sub-tasks (gather), safe for concurrent calls.
_progress_cb: ContextVar[Callable[[str], Coroutine] | None] = ContextVar(
    "_progress_cb", default=None
)


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

    # Gene signature analysis (Issue #10)
    gene_signature: Optional[list]    # input DEG list for analyze_gene_signature
    master_regulators: Optional[dict] # RegNetAgents find_master_regulators output

    # Cross-cell-type comparison (Issue #11)
    cell_types: Optional[list]        # input list of cell types for compare_cell_contexts
    comparison_results: Optional[dict] # per-cell-type {network, perturbation} results

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
        # Persistent connections set by MCP server at startup to avoid per-call cold starts
        self._persistent_cascade: Any = None
        self._persistent_regnetagents: Any = None

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

    async def _emit(self, msg: str) -> None:
        cb = _progress_cb.get()
        if cb is not None:
            try:
                await cb(msg)
            except Exception:
                pass

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
        graph.add_node("run_signature_path", self._run_signature_path)
        graph.add_node("run_comparison_path", self._run_comparison_path)
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
                "signature_path": "run_signature_path",
                "comparison_path": "run_comparison_path",
            },
        )
        graph.add_edge("run_tf_path", "synthesize")
        graph.add_edge("run_effector_path", "synthesize")
        graph.add_edge("run_validation_path", "synthesize")
        graph.add_edge("run_signature_path", "synthesize")
        graph.add_edge("run_comparison_path", "synthesize")
        graph.add_edge("synthesize", "generate_report")
        graph.add_edge("generate_report", END)

        return graph.compile()

    # ------------------------------------------------------------------
    # Routing
    # ------------------------------------------------------------------

    def _routing_decision(self, state: OrchestraState) -> str:
        if state.get("analysis_type") == "gene_signature":
            return "signature_path"
        if state.get("analysis_type") == "cell_context_comparison":
            return "comparison_path"
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
        # Signature and comparison analyses have no single cell_type to classify — skip.
        if state.get("analysis_type") in ("gene_signature", "cell_context_comparison"):
            state["completed_steps"].append("classify_gene")
            return state
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

        await self._emit(f"[Orchestra] Running RegNetAgents + CASCADE in parallel for {gene} in {cell_type}...")
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
        await self._emit(f"[Orchestra] Finding TF partner for {gene} via STRING PPI...")
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
        await self._emit(f"[Orchestra] Running parallel CASCADE + RegNetAgents for TF partner {tf_partner}...")
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

    async def _run_signature_path(self, state: OrchestraState) -> OrchestraState:
        """
        Gene signature driver analysis path (Issue #10).

        1. RegNetAgents find_master_regulators — ranks TFs by Fisher enrichment in the
           input gene set; uses ARACNe regulon overlap, not MR inference.
        2. Parallel CASCADE comprehensive_perturbation_analysis on the top 3 TFs —
           adds experimental validation (LINCS, DepMap, super-enhancers, etc.).

        Synthesis in _synthesize_signature_path combines:
          - Signature coverage % (overlap_count / signature_size, from RegNetAgents)
          - 7-source corroboration count (from CASCADE evidence_synthesis)
        to produce a ranked driver table neither system alone can generate.
        """
        gene_signature = state.get("gene_signature") or []
        cell_type = state["cell_type"]

        if not gene_signature:
            state["errors"]["signature_path"] = "gene_signature is empty"
            state["completed_steps"].append("run_signature_path")
            return state

        # Step 1: RegNetAgents master regulator enrichment analysis
        await self._emit(
            f"[Orchestra] Running Fisher enrichment across all regulators "
            f"for {len(gene_signature)}-gene signature in {cell_type} "
            f"(this may take several minutes)..."
        )
        try:
            mr_result = await self._regnetagents.call_tool(
                "find_master_regulators",
                {"gene_set": gene_signature, "cell_type": cell_type, "top_n": 10},
                timeout_seconds=TIMEOUT_MASTER_REGULATORS,
            )
            state["master_regulators"] = mr_result
        except Exception as e:
            state["errors"]["master_regulators"] = str(e)
            state["completed_steps"].append("run_signature_path")
            return state

        # Step 2: parallel CASCADE validation on top 3 enriched TFs
        top_tfs = [
            r["gene"]
            for r in (mr_result.get("master_regulators") or [])[:3]
        ]
        if top_tfs:
            await self._emit(f"[Orchestra] Validating top TFs via CASCADE: {', '.join(top_tfs)}")
            validation_results = await asyncio.gather(
                *[
                    self._cascade.call_tool(
                        "comprehensive_perturbation_analysis",
                        {"gene": tf, "cell_type": cell_type},
                        timeout_seconds=TIMEOUT_PERTURBATION,
                    )
                    for tf in top_tfs
                ],
                return_exceptions=True,
            )
            # Attach CASCADE evidence back onto the master regulator entries
            mr_list = mr_result.get("master_regulators") or []
            for i, (tf, result) in enumerate(zip(top_tfs, validation_results)):
                entry = mr_list[i] if i < len(mr_list) else None
                if entry is None:
                    continue
                if isinstance(result, Exception):
                    entry["cascade_error"] = str(result)
                else:
                    ev = (result or {}).get("evidence_synthesis") or {}
                    entry["key_findings"] = ev.get("key_findings", [])
                    entry["multi_source_genes"] = [
                        g["symbol"] for g in (ev.get("multi_source_genes") or [])[:10]
                    ]

        state["completed_steps"].append("run_signature_path")
        return state

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

        await self._emit(f"[Orchestra] Querying RegNetAgents PageRank + CASCADE drug discovery for {gene} in {cell_type}...")
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
            await self._emit(f"[Orchestra] Validating top candidates via CASCADE: {', '.join(c['gene'] for c in top)}")
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

    async def _run_comparison_path(self, state: OrchestraState) -> OrchestraState:
        """
        Cross-cell-type comparison path (Issue #11).

        For each cell type in cell_types, runs RegNetAgents comprehensive_gene_analysis
        and CASCADE comprehensive_perturbation_analysis in parallel (2N total calls).
        Results are stored per-cell-type in comparison_results for conservation scoring
        in _synthesize_comparison_path.
        """
        gene = state["gene"]
        cell_types = state.get("cell_types") or []

        if not cell_types:
            state["errors"]["comparison"] = "cell_types list is empty"
            state["completed_steps"].append("run_comparison_path")
            return state

        await self._emit(
            f"[Orchestra] Running parallel RegNetAgents + CASCADE for {gene} "
            f"across {len(cell_types)} cell types: {', '.join(cell_types)}"
        )
        async def _analyze_one(ct: str):
            rna, casc = await asyncio.gather(
                self._regnetagents.call_tool(
                    "comprehensive_gene_analysis",
                    {"gene": gene, "cell_type": ct},
                    timeout_seconds=TIMEOUT_NETWORK_COMPARISON,
                ),
                self._cascade.call_tool(
                    "comprehensive_perturbation_analysis",
                    {"gene": gene, "cell_type": ct},
                    timeout_seconds=TIMEOUT_PERTURBATION,
                ),
                return_exceptions=True,
            )
            return (
                ct,
                rna if not isinstance(rna, Exception) else None,
                str(rna) if isinstance(rna, Exception) else None,
                casc if not isinstance(casc, Exception) else None,
                str(casc) if isinstance(casc, Exception) else None,
            )

        all_results = await asyncio.gather(
            *[_analyze_one(ct) for ct in cell_types],
            return_exceptions=True,
        )

        comparison_results: dict = {}
        for i, item in enumerate(all_results):
            ct = cell_types[i] if i < len(cell_types) else f"cell_type_{i}"
            if isinstance(item, Exception):
                comparison_results[ct] = {
                    "network": None, "network_error": str(item),
                    "perturbation": None, "perturbation_error": str(item),
                }
            else:
                ct_name, net, net_err, perturb, perturb_err = item
                comparison_results[ct_name] = {
                    "network": net, "network_error": net_err,
                    "perturbation": perturb, "perturbation_error": perturb_err,
                }

        state["comparison_results"] = comparison_results
        state["completed_steps"].append("run_comparison_path")
        return state

    # ------------------------------------------------------------------
    # Synthesis
    # ------------------------------------------------------------------

    def _synthesize(self, state: OrchestraState) -> OrchestraState:
        if state.get("analysis_type") == "gene_signature":
            return self._synthesize_signature_path(state)
        if state.get("analysis_type") == "cell_context_comparison":
            return self._synthesize_comparison_path(state)
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

        errors = state.get("errors", {})
        discordance_flags = self._compute_tf_discordance_flags(
            multi_source_genes, rna_targets, cross_system_hits
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
            "regnetagents_available": bool(network) and "network" not in errors,
            "cascade_available": bool(perturbation) and "perturbation" not in errors,
            "discordance_flags": discordance_flags,
            "errors": errors,
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

        errors = state.get("errors", {})
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
            "regnetagents_available": bool(network) and "network" not in errors,
            "cascade_available": bool(perturbation) and "perturbation" not in errors,
            "discordance_flags": [],
            "errors": errors,
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

    def _is_in_pathway_enrichment(self, gene: str, network_analysis: dict) -> bool:
        """
        Check if a gene symbol appears in RegNetAgents pathway enrichment gene lists.
        Searches pathway-related sections only to avoid matching the query gene itself.
        """
        if not gene or not network_analysis:
            return False
        gene_upper = gene.upper()

        def _in_gene_list(obj) -> bool:
            if not isinstance(obj, list):
                return False
            for item in obj:
                if isinstance(item, str) and item.upper() == gene_upper:
                    return True
                if isinstance(item, dict):
                    sym = (
                        item.get("gene_symbol") or item.get("symbol")
                        or item.get("gene") or item.get("name") or ""
                    )
                    if isinstance(sym, str) and sym.upper() == gene_upper:
                        return True
            return False

        def _search(obj, depth: int = 0) -> bool:
            if depth > 6:
                return False
            if isinstance(obj, dict):
                for k, v in obj.items():
                    k_lower = k.lower()
                    if k_lower in ("genes", "gene_list", "gene_set", "members",
                                   "pathway_genes", "leading_edge"):
                        if _in_gene_list(v):
                            return True
                    elif ("pathway" in k_lower or "reactome" in k_lower
                          or "enrichment" in k_lower or "gene_set" in k_lower):
                        if _search(v, depth + 1):
                            return True
            elif isinstance(obj, list):
                return any(_search(item, depth + 1) for item in obj)
            return False

        return _search(network_analysis)

    def _score_candidate_evidence(
        self, candidate: dict, network_analysis: dict
    ) -> dict:
        """
        Score a validation candidate against 7 independent evidence sources.

        Returns per-source boolean flags and a total corroboration_count (out of 7).
        The 7 sources span two methodologically independent systems:
          RegNetAgents: PageRank rank, pathway membership
          CASCADE: LINCS knockdown, DepMap essentiality, super-enhancer, DoRothEA, cBioPortal
        """
        # Gather all text signals: key_findings + reason + explanation + any other string fields
        text_parts = list(candidate.get("key_findings") or [])
        for field in ("reason", "explanation", "notes", "description"):
            val = candidate.get(field)
            if isinstance(val, str):
                text_parts.append(val)
        findings_text = "\n".join(text_parts).lower()

        # RegNetAgents sources
        pagerank_hit = (
            candidate.get("source") == "regnetagents_pagerank"
            or float(candidate.get("pagerank") or 0) > 0
        )
        pathway_hit = self._is_in_pathway_enrichment(candidate["gene"], network_analysis)

        # CASCADE sources — inferred from key_findings text
        lincs_hit = "lincs" in findings_text
        # DepMap: "not essential" is a negative signal; any other DepMap mention is positive
        depmap_hit = "depmap" in findings_text and "not essential" not in findings_text
        se_hit = (
            "super-enhancer" in findings_text
            or "super_enhancer" in findings_text
            or "bet inhibitor" in findings_text
        )
        dorothea_hit = "dorothea" in findings_text
        cbio_hit = "cbioportal" in findings_text or "cbio" in findings_text

        flags: dict[str, bool] = {
            "pagerank_rank": pagerank_hit,
            "pathway_member": pathway_hit,
            "lincs_knockdown": lincs_hit,
            "depmap_essentiality": depmap_hit,
            "super_enhancer": se_hit,
            "dorothea_tier": dorothea_hit,
            "cbio_expression": cbio_hit,
        }
        return {
            **flags,
            "corroboration_count": sum(1 for v in flags.values() if v),
            "corroboration_denominator": 7,
        }

    def _compute_tf_discordance_flags(
        self,
        multi_source_genes: list,
        rna_targets: set,
        cross_system_hits: list,
    ) -> list:
        """
        Detect discordances between CASCADE experimental evidence and RegNetAgents
        network topology for the TF analysis path.

        Two patterns:
        - Genes with CASCADE experimental support absent from the regulatory network
          (experimentally active but not a classical TF target in this context).
        - Large regulatory network with no experimentally corroborated targets
          (topological predictions without experimental validation).
        """
        flags = []

        cascade_only = [
            g for g in multi_source_genes if g["symbol"] not in rna_targets
        ][:5]
        if cascade_only:
            flags.append({
                "type": "experimentally_active_not_in_network",
                "description": (
                    "Experimentally supported by CASCADE but absent from the ARACNe/GREmLN "
                    "mRNA regulatory network — regulatory mechanism uncharacterized; "
                    "absence from the mRNA network does not imply non-transcriptional biology"
                ),
                "genes": [
                    {"symbol": g["symbol"], "source_count": g["source_count"],
                     "sources": g["sources"]}
                    for g in cascade_only
                ],
            })

        if rna_targets and multi_source_genes and not cross_system_hits:
            flags.append({
                "type": "network_topology_without_experimental_support",
                "description": (
                    f"RegNetAgents identifies {len(rna_targets)} downstream network targets "
                    "but none overlap with CASCADE multi-source experimental data — "
                    "may reflect cell context differences or stale network edges"
                ),
                "genes": [],
            })

        return flags

    def _compute_validation_discordance_flags(self, evidence_table: list) -> list:
        """
        Detect per-candidate discordances in the validation path evidence table.

        Two patterns:
        - High network rank (PageRank) with no CASCADE experimental support →
          topological hub not experimentally validated.
        - CASCADE experimental support with no network rank →
          experimentally active but absent from ARACNe network (BRD4 pattern).
        """
        flags = []

        for row in evidence_table:
            gene_name = row["gene"]
            experimental_support = any([
                row["lincs_knockdown"],
                row["depmap_essentiality"],
                row["super_enhancer"],
                row["dorothea_tier"],
                row["cbio_expression"],
            ])

            if row["pagerank_rank"] and not experimental_support:
                flags.append({
                    "type": "topological_hub_not_validated",
                    "gene": gene_name,
                    "description": (
                        f"{gene_name} ranks in the regulatory network (PageRank) "
                        "but has no experimental support from CASCADE — "
                        "may be topologically central but not functionally essential "
                        "in this cell type"
                    ),
                    "genes": [],
                })

            if experimental_support and not row["pagerank_rank"]:
                flags.append({
                    "type": "experimentally_active_not_in_network",
                    "gene": gene_name,
                    "description": (
                        f"{gene_name} has CASCADE experimental support but is absent "
                        "from the ARACNe/GREmLN mRNA regulatory network — regulatory "
                        "mechanism uncharacterized; absence from the mRNA network does "
                        "not imply non-transcriptional biology"
                    ),
                    "genes": [],
                })

        return flags

    def _synthesize_signature_path(self, state: OrchestraState) -> OrchestraState:
        """
        Signature path synthesis: combines RegNetAgents Fisher enrichment (network
        coverage of the input signature) with CASCADE experimental corroboration
        (7-source evidence table) for each candidate TF driver.

        Ranking: primary = overlap_count (signature genes in TF regulon),
                 secondary = corroboration_count (CASCADE experimental support).
        Reuses _score_candidate_evidence and _compute_validation_discordance_flags.
        """
        gene_signature = state.get("gene_signature") or []
        mr_result = state.get("master_regulators") or {}
        errors = state.get("errors", {})
        mr_list = mr_result.get("master_regulators") or []
        query_summary = mr_result.get("query_summary") or {}

        signature_size = len(gene_signature)

        # Build ranked_drivers list — one entry per TF candidate
        ranked_drivers = []
        evidence_table = []
        for entry in mr_list:
            overlap_count = entry.get("overlap_count", 0)
            coverage_pct = round(overlap_count / signature_size * 100, 1) if signature_size else 0.0

            # Build a candidate dict compatible with _score_candidate_evidence.
            # enrichment_score > 0 acts as the "network" signal (maps to pagerank_rank flag).
            candidate = {
                "gene": entry["gene"],
                "source": "regnetagents_enrichment",
                "pagerank": entry.get("enrichment_score", 0),  # > 0 → pagerank_rank=True
                "key_findings": entry.get("key_findings", []),
            }
            scores = self._score_candidate_evidence(candidate, {})

            row = {
                "gene": entry["gene"],
                "overlap_count": overlap_count,
                "coverage_pct": coverage_pct,
                "regulon_size": entry.get("regulon_size", 0),
                "enrichment_score": entry.get("enrichment_score", 0),
                "p_value": entry.get("p_value", 1.0),
                "overlapping_genes": entry.get("overlapping_genes", []),
                **scores,
            }
            evidence_table.append(row)

            ranked_drivers.append({
                "gene": entry["gene"],
                "overlap_count": overlap_count,
                "coverage_pct": coverage_pct,
                "enrichment_score": entry.get("enrichment_score", 0),
                "p_value": entry.get("p_value", 1.0),
                "corroboration_count": scores["corroboration_count"],
                "corroboration_denominator": scores["corroboration_denominator"],
                "overlapping_genes": entry.get("overlapping_genes", []),
                "cascade_key_findings": entry.get("key_findings", []),
                "multi_source_genes": entry.get("multi_source_genes", []),
                "cascade_error": entry.get("cascade_error"),
            })

        # Sort by overlap_count desc, corroboration_count desc as tiebreaker
        ranked_drivers.sort(key=lambda x: (-x["overlap_count"], -x["corroboration_count"]))
        evidence_table.sort(key=lambda x: (-x["overlap_count"], -x.get("corroboration_count", 0)))

        regnetagents_available = bool(mr_list) and "master_regulators" not in errors
        cascade_available = any(
            "key_findings" in entry and "cascade_error" not in entry
            for entry in mr_list
        )

        discordance_flags = self._compute_validation_discordance_flags(evidence_table)

        state["synthesis"] = {
            "gene": "",
            "cell_type": state["cell_type"],
            "routing": "signature",
            "gene_signature": gene_signature,
            "signature_size": signature_size,
            "genes_found_in_network": query_summary.get("genes_found_in_network", 0),
            "genes_not_found": query_summary.get("genes_not_found", []),
            "network_size": query_summary.get("network_size", 0),
            "total_regulators_tested": query_summary.get("total_regulators_tested", 0),
            "ranked_drivers": ranked_drivers,
            "evidence_table": evidence_table,
            "regnetagents_available": regnetagents_available,
            "cascade_available": cascade_available,
            "discordance_flags": discordance_flags,
            "errors": errors,
        }
        state["completed_steps"].append("synthesize")
        return state

    def _synthesize_validation_path(self, state: OrchestraState) -> OrchestraState:
        """
        Validation path synthesis: builds per-candidate 7-source corroboration table
        from RegNetAgents network rank + CASCADE perturbation evidence.

        Evidence sources:
          RegNetAgents: PageRank rank, pathway membership (2 sources)
          CASCADE: LINCS knockdown, DepMap essentiality, super-enhancer,
                   DoRothEA TF tier, cBioPortal expression (5 sources)
        Corroboration count = number of independent sources supporting each candidate.
        """
        validated_targets = state.get("validated_targets") or []
        network = state.get("network_analysis") or {}
        errors = state.get("errors", {})
        network_summary = (
            network.get("summary")
            or network.get("network_analysis")
            or network.get("workflow_summary")
            or {}
        )

        # Score each validated candidate (top 3 that received CASCADE perturbation calls)
        evidence_table = []
        for candidate in validated_targets:
            if candidate.get("key_findings") is not None:
                scores = self._score_candidate_evidence(candidate, network)
                evidence_table.append({"gene": candidate["gene"], **scores})

        # Sort by corroboration_count descending so highest-evidence candidate leads
        evidence_table.sort(key=lambda x: x.get("corroboration_count", 0), reverse=True)

        # Graceful degradation flags — used by formatter to warn on partial data
        regnetagents_available = bool(network) and "network" not in errors
        cascade_available = any(
            c.get("key_findings") is not None and "cascade_error" not in c
            for c in validated_targets
        )

        discordance_flags = self._compute_validation_discordance_flags(evidence_table)
        state["synthesis"] = {
            "gene": state["gene"],
            "cell_type": state["cell_type"],
            "routing": "validation",
            "validated_targets": validated_targets,
            "evidence_table": evidence_table,
            "network_context": network_summary,
            "regnetagents_available": regnetagents_available,
            "cascade_available": cascade_available,
            "discordance_flags": discordance_flags,
            "errors": errors,
        }
        state["completed_steps"].append("synthesize")
        return state

    def _synthesize_comparison_path(self, state: OrchestraState) -> OrchestraState:
        """
        Cross-cell-type comparison synthesis (Issue #11).

        Scores each cell type against the 7-source evidence framework, then classifies
        each source by conservation level using the 2/3-majority threshold (min 2):
          conserved        — present in >= ceil(2/3 * N) cell types
          enriched         — present in 2 to threshold-1 cell types
          cell_type_specific — present in exactly 1 cell type
          absent           — present in 0 cell types
        """
        gene = state["gene"]
        cell_types = state.get("cell_types") or []
        comparison_results = state.get("comparison_results") or {}

        _SOURCES = [
            "pagerank_rank", "pathway_member", "lincs_knockdown",
            "depmap_essentiality", "super_enhancer", "dorothea_tier", "cbio_expression",
        ]

        per_cell_type: dict = {}
        for ct in cell_types:
            ct_data = comparison_results.get(ct) or {}
            scores = self._score_cell_type_evidence(
                ct_data.get("network"), ct_data.get("perturbation")
            )
            scores["network_error"] = ct_data.get("network_error")
            scores["perturbation_error"] = ct_data.get("perturbation_error")
            per_cell_type[ct] = scores

        n = len(cell_types)
        threshold = max(2, math.ceil(2 / 3 * n))

        conservation_scores: dict = {}
        for src in _SOURCES:
            count = sum(
                1 for ct in cell_types
                if per_cell_type.get(ct, {}).get(src, False)
            )
            conservation_scores[src] = {
                "count": count,
                "n": n,
                "label": self._classify_conservation(count, n, threshold),
            }

        state["synthesis"] = {
            "routing": "comparison",
            "gene": gene,
            "cell_types": cell_types,
            "per_cell_type": per_cell_type,
            "conservation_scores": conservation_scores,
            "errors": state.get("errors") or {},
        }
        state["completed_steps"].append("synthesize")
        return state

    def _classify_conservation(self, count: int, n: int, threshold: int) -> str:
        if count == 0:
            return "absent"
        if count == 1:
            return "cell_type_specific"
        if count >= threshold:
            return "conserved"
        return "enriched"

    def _score_cell_type_evidence(
        self,
        network_result: Optional[dict],
        perturbation_result: Optional[dict],
    ) -> dict:
        """
        Score regulatory evidence for the query gene in one cell type context.

        RegNetAgents sources: hub status (pagerank_rank) and pathway enrichment.
        CASCADE sources: extracted from evidence_synthesis.key_findings text —
        same keyword signals as _score_candidate_evidence.
        """
        net = network_result or {}
        # comprehensive_gene_analysis returns network_analysis = state['gene_info']
        # which has regulatory_role ("hub_regulator", ...) and pagerank (float)
        net_analysis = net.get("network_analysis") or {}
        regulatory_role = net_analysis.get("regulatory_role") or ""
        pagerank_rank = regulatory_role == "hub_regulator"

        enriched_pathways = (
            (net.get("pathway_enrichment") or {}).get("enriched_pathways")
            or (net.get("pathway_enrichment") or {}).get("pathways")
            or []
        )
        pathway_member = bool(enriched_pathways)

        perturb = perturbation_result or {}
        ev = (perturb.get("evidence_synthesis") or {})
        key_findings = ev.get("key_findings") or []
        findings_text = "\n".join(str(f) for f in key_findings).lower()

        lincs_hit = "lincs" in findings_text
        depmap_hit = "depmap" in findings_text and "not essential" not in findings_text
        se_hit = (
            "super-enhancer" in findings_text
            or "super_enhancer" in findings_text
            or "bet inhibitor" in findings_text
        )
        dorothea_hit = "dorothea" in findings_text
        cbio_hit = "cbioportal" in findings_text or "cbio" in findings_text

        flags: dict = {
            "pagerank_rank": pagerank_rank,
            "pathway_member": pathway_member,
            "lincs_knockdown": lincs_hit,
            "depmap_essentiality": depmap_hit,
            "super_enhancer": se_hit,
            "dorothea_tier": dorothea_hit,
            "cbio_expression": cbio_hit,
        }
        return {
            **flags,
            "corroboration_count": sum(1 for v in flags.values() if v),
            "corroboration_denominator": 7,
            "network_available": network_result is not None,
            "perturbation_available": perturbation_result is not None,
        }

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
            else routing in ("tf", "validation", "signature", "comparison")
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
        if routing == "signature":
            return self._format_signature_report(synthesis)
        if routing == "comparison":
            return self._format_comparison_report(synthesis)
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

        regnetagents_available = synthesis.get("regnetagents_available", True)
        cascade_available = synthesis.get("cascade_available", True)

        lines = [
            f"## Orchestra Analysis: {gene} in {cell_type}",
            f"**Routing:** TF / {gene_role}",
            "",
        ]
        if not regnetagents_available:
            lines.append("> ⚠️ **RegNetAgents unavailable** — network topology and pathway evidence missing; showing CASCADE-only results.")
            lines.append("")
        if not cascade_available:
            lines.append("> ⚠️ **CASCADE unavailable** — experimental validation missing; showing RegNetAgents-only results.")
            lines.append("")
        lines.append("### CASCADE Evidence")

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

        discordance_flags = synthesis.get("discordance_flags", [])
        if discordance_flags:
            lines.append("")
            lines.append("### Notable Discordances")
            for flag in discordance_flags:
                lines.append(f"- **{flag['description']}**")
                for g in flag.get("genes", [])[:5]:
                    sym = g.get("symbol", "?")
                    sc = g.get("source_count", "")
                    srcs = g.get("sources", [])
                    detail = f"{sc} sources ({', '.join(srcs)})" if sc else ""
                    lines.append(f"  - {sym}: {detail}" if detail else f"  - {sym}")

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
        regnetagents_available = synthesis.get("regnetagents_available", True)
        cascade_available = synthesis.get("cascade_available", True)

        lines = [
            f"## Orchestra Analysis: {gene} in {cell_type}",
            f"**Routing:** effector/scaffold",
            f"**TF partner (via PPI):** {tf_partner}",
            "",
        ]
        if not regnetagents_available:
            lines.append("> ⚠️ **RegNetAgents unavailable** — pathway context missing; showing CASCADE-only results.")
            lines.append("")
        if not cascade_available:
            lines.append("> ⚠️ **CASCADE unavailable** — PPI and perturbation data missing.")
            lines.append("")
        lines.append("### CASCADE Evidence")

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

        discordance_flags = synthesis.get("discordance_flags", [])
        if discordance_flags:
            lines.append("")
            lines.append("### Notable Discordances")
            for flag in discordance_flags:
                lines.append(f"- **{flag['description']}**")
                for g in flag.get("genes", [])[:5]:
                    sym = g.get("symbol", "?")
                    sc = g.get("source_count", "")
                    srcs = g.get("sources", [])
                    detail = f"{sc} sources ({', '.join(srcs)})" if sc else ""
                    lines.append(f"  - {sym}: {detail}" if detail else f"  - {sym}")

        if errors:
            lines.append("")
            lines.append("### Partial Data Warnings")
            for k, v in errors.items():
                lines.append(f"- {k}: {v}")

        return lines

    def _format_signature_report(self, synthesis: dict) -> list[str]:
        cell_type = synthesis.get("cell_type", "unknown")
        gene_signature = synthesis.get("gene_signature") or []
        signature_size = synthesis.get("signature_size", len(gene_signature))
        genes_found = synthesis.get("genes_found_in_network", 0)
        genes_not_found = synthesis.get("genes_not_found", [])
        total_tested = synthesis.get("total_regulators_tested", 0)
        ranked_drivers = synthesis.get("ranked_drivers") or []
        evidence_table = synthesis.get("evidence_table") or []
        errors = synthesis.get("errors", {})
        regnetagents_available = synthesis.get("regnetagents_available", True)
        cascade_available = synthesis.get("cascade_available", True)

        lines = [
            f"## Orchestra Gene Signature Driver Analysis — {cell_type}",
            f"**Signature size:** {signature_size} genes  |  "
            f"**Found in network:** {genes_found}  |  "
            f"**Regulators tested:** {total_tested}",
            "",
        ]

        if not regnetagents_available:
            lines.append("> ⚠️ **RegNetAgents unavailable** — master regulator enrichment missing.")
            lines.append("")
        if not cascade_available:
            lines.append("> ⚠️ **CASCADE unavailable** — experimental validation missing; showing network enrichment only.")
            lines.append("")

        if not ranked_drivers:
            lines.append("No master regulators identified for this gene signature.")
            if genes_not_found:
                lines.append(f"Genes not found in network: {', '.join(genes_not_found[:10])}")
            if errors:
                lines.append("")
                lines.append("### Errors")
                for k, v in errors.items():
                    lines.append(f"- {k}: {v}")
            return lines

        # Corroboration + coverage table
        if evidence_table:
            lines.append("### Ranked TF Drivers")
            lines.append("")
            lines.append(
                "| TF Driver | Coverage | Overlap | Enrichment | p-value"
                " | PageRank | Pathway | LINCS | DepMap | SE | DoRothEA | cBio | Score |"
            )
            lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")

            def flag(v: bool) -> str:
                return "✓" if v else "-"

            for row in evidence_table:
                lines.append(
                    f"| {row['gene']} "
                    f"| {row['coverage_pct']}% "
                    f"| {row['overlap_count']}/{row.get('regulon_size', '?')} "
                    f"| {row['enrichment_score']:.2f} "
                    f"| {row['p_value']:.2e} "
                    f"| {flag(row['pagerank_rank'])} "
                    f"| {flag(row['pathway_member'])} "
                    f"| {flag(row['lincs_knockdown'])} "
                    f"| {flag(row['depmap_essentiality'])} "
                    f"| {flag(row['super_enhancer'])} "
                    f"| {flag(row['dorothea_tier'])} "
                    f"| {flag(row['cbio_expression'])} "
                    f"| **{row['corroboration_count']}/{row['corroboration_denominator']}** |"
                )
            lines.append("")
            lines.append(
                "_Coverage = signature genes in TF regulon. "
                "PageRank/Pathway = RegNetAgents enrichment. "
                "LINCS, DepMap, SE, DoRothEA, cBio = CASCADE sources._"
            )
            lines.append("")

        # Per-driver detail for top 3
        lines.append("### Driver Details")
        for i, d in enumerate(ranked_drivers[:3], 1):
            lines.append(f"**{i}. {d['gene']}** — {d['coverage_pct']}% coverage "
                         f"({d['overlap_count']} signature genes), "
                         f"enrichment={d['enrichment_score']:.2f}, "
                         f"p={d['p_value']:.2e}")
            if d.get("overlapping_genes"):
                lines.append(f"   Overlapping: {', '.join(d['overlapping_genes'][:10])}")
            cascade_findings = d.get("cascade_key_findings") or []
            if cascade_findings:
                lines.append("   CASCADE evidence:")
                for f in cascade_findings[:3]:
                    lines.append(f"   - {f}")
            elif d.get("cascade_error"):
                lines.append(f"   CASCADE validation error: {d['cascade_error']}")
            lines.append("")

        if genes_not_found:
            lines.append(f"_Genes not found in network ({len(genes_not_found)}): "
                         f"{', '.join(genes_not_found[:10])}_")
            lines.append("")

        discordance_flags = synthesis.get("discordance_flags", [])
        if discordance_flags:
            lines.append("### Notable Discordances")
            for flag_entry in discordance_flags:
                lines.append(f"- **{flag_entry['description']}**")
            lines.append("")

        if errors:
            lines.append("### Partial Data Warnings")
            for k, v in errors.items():
                lines.append(f"- {k}: {v}")

        return lines

    def _format_comparison_report(self, synthesis: dict) -> list[str]:
        gene = synthesis.get("gene", "unknown")
        cell_types = synthesis.get("cell_types") or []
        per_cell_type = synthesis.get("per_cell_type") or {}
        conservation_scores = synthesis.get("conservation_scores") or {}
        errors = synthesis.get("errors") or {}

        _SOURCE_LABELS = {
            "pagerank_rank": "PageRank hub",
            "pathway_member": "Pathway enrichment",
            "lincs_knockdown": "LINCS knockdown",
            "depmap_essentiality": "DepMap essentiality",
            "super_enhancer": "Super-enhancer",
            "dorothea_tier": "DoRothEA TF",
            "cbio_expression": "cBioPortal",
        }
        _CONSERVATION_SYMBOLS = {
            "conserved": "conserved",
            "enriched": "enriched",
            "cell_type_specific": "cell-type-specific",
            "absent": "absent",
        }

        lines: list[str] = []
        lines.append(f"## Cross-Cell-Type Context Comparison: {gene}")
        lines.append("")
        lines.append(f"**Cell types compared:** {', '.join(cell_types)}")
        lines.append(f"**Conservation threshold:** ≥ ⌈2/3 × {len(cell_types)}⌉ = "
                     f"{max(2, math.ceil(2 / 3 * len(cell_types)))} cell types")
        lines.append("")

        if not cell_types:
            lines.append("> ⚠️ No cell types provided.")
            return lines

        # Evidence heatmap table
        lines.append("### Evidence Heatmap")
        lines.append("")
        header = "| Source | " + " | ".join(cell_types) + " | Conservation |"
        separator = "|---|" + "|".join(["---"] * len(cell_types)) + "|---|"
        lines.append(header)
        lines.append(separator)

        for src, label in _SOURCE_LABELS.items():
            cons = conservation_scores.get(src) or {}
            count = cons.get("count", 0)
            n = cons.get("n", len(cell_types))
            cons_label = _CONSERVATION_SYMBOLS.get(cons.get("label", "absent"), "absent")
            cells = []
            for ct in cell_types:
                score = per_cell_type.get(ct) or {}
                present = score.get(src, False)
                net_avail = score.get("network_available", False)
                perturb_avail = score.get("perturbation_available", False)
                if src in ("pagerank_rank", "pathway_member"):
                    avail = net_avail
                else:
                    avail = perturb_avail
                if not avail:
                    cells.append("N/A")
                elif present:
                    cells.append("✓")
                else:
                    cells.append("—")
            cons_cell = f"{cons_label} ({count}/{n})"
            lines.append(f"| {label} | " + " | ".join(cells) + f" | {cons_cell} |")

        lines.append("")
        lines.append(
            "_PageRank hub, Pathway enrichment = RegNetAgents sources (ARACNe/GREmLN). "
            "LINCS, DepMap, SuperEnhancer, DoRothEA, cBioPortal = CASCADE sources. "
            "Independent of mRNA network: LINCS knockdown, DepMap essentiality, SuperEnhancer._"
        )
        lines.append("")

        # Conservation summary
        lines.append("### Conservation Summary")
        lines.append("")
        by_label: dict = {"conserved": [], "enriched": [], "cell_type_specific": [], "absent": []}
        for src, cons in conservation_scores.items():
            by_label.setdefault(cons.get("label", "absent"), []).append(
                _SOURCE_LABELS.get(src, src)
            )

        if by_label.get("conserved"):
            lines.append(f"**Conserved (≥2/3 cell types):** {', '.join(by_label['conserved'])}")
        if by_label.get("enriched"):
            lines.append(f"**Enriched (majority):** {', '.join(by_label['enriched'])}")
        if by_label.get("cell_type_specific"):
            # Find which cell type for each specific source
            specifics = []
            for src, cons in conservation_scores.items():
                if cons.get("label") == "cell_type_specific":
                    for ct in cell_types:
                        if (per_cell_type.get(ct) or {}).get(src, False):
                            specifics.append(f"{_SOURCE_LABELS.get(src, src)} ({ct} only)")
                            break
            lines.append(f"**Cell-type-specific:** {', '.join(specifics)}")
        if by_label.get("absent"):
            lines.append(f"**Absent across all cell types:** {', '.join(by_label['absent'])}")

        lines.append("")

        # Data availability warnings
        unavailable = [
            ct for ct in cell_types
            if not (per_cell_type.get(ct) or {}).get("network_available", True)
            or not (per_cell_type.get(ct) or {}).get("perturbation_available", True)
        ]
        if unavailable:
            lines.append("### Data Availability")
            lines.append("")
            for ct in unavailable:
                scores = per_cell_type.get(ct) or {}
                if scores.get("network_error"):
                    lines.append(f"> ⚠️ **{ct}**: RegNetAgents unavailable — {scores['network_error']}")
                if scores.get("perturbation_error"):
                    lines.append(f"> ⚠️ **{ct}**: CASCADE unavailable — {scores['perturbation_error']}")
            lines.append("")

        if errors:
            lines.append("### Errors")
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
        evidence_table = synthesis.get("evidence_table") or []
        errors = synthesis.get("errors", {})
        regnetagents_available = synthesis.get("regnetagents_available", True)
        cascade_available = synthesis.get("cascade_available", True)

        lines = [
            f"## Orchestra Therapeutic Target Validation: {gene} in {cell_type}",
            f"**Routing:** therapeutic_validation",
            f"**Candidates evaluated:** {len(candidates)}",
            "",
        ]

        if not regnetagents_available:
            lines.append("> ⚠️ **RegNetAgents unavailable** — PageRank and pathway evidence missing; showing CASCADE-only results.")
            lines.append("")
        if not cascade_available:
            lines.append("> ⚠️ **CASCADE unavailable** — experimental validation missing; showing RegNetAgents-only results.")
            lines.append("")

        if not candidates:
            lines.append("No therapeutic target candidates identified.")
            if errors:
                lines.append("")
                lines.append("### Errors")
                for k, v in errors.items():
                    lines.append(f"- {k}: {v}")
            return lines

        # Corroboration table: structured 7-source evidence summary
        if evidence_table:
            lines.append("### Evidence Corroboration Table")
            lines.append("")
            lines.append(
                "| Candidate | PageRank | Pathway | LINCS | DepMap | SuperEnhancer"
                " | DoRothEA | cBioPortal | Score |"
            )
            lines.append("|---|---|---|---|---|---|---|---|---|")
            for row in evidence_table:
                def flag(v: bool) -> str:
                    return "✓" if v else "-"
                lines.append(
                    f"| {row['gene']} "
                    f"| {flag(row['pagerank_rank'])} "
                    f"| {flag(row['pathway_member'])} "
                    f"| {flag(row['lincs_knockdown'])} "
                    f"| {flag(row['depmap_essentiality'])} "
                    f"| {flag(row['super_enhancer'])} "
                    f"| {flag(row['dorothea_tier'])} "
                    f"| {flag(row['cbio_expression'])} "
                    f"| **{row['corroboration_count']}/{row['corroboration_denominator']}** |"
                )
            lines.append("")
            lines.append(
                "_PageRank, Pathway = RegNetAgents sources (ARACNe/GREmLN mRNA network). "
                "LINCS, DepMap, SuperEnhancer, DoRothEA, cBioPortal = CASCADE sources. "
                "Methodologically independent of mRNA network: LINCS knockdown, DepMap CRISPR, SuperEnhancer. "
                "Partially shared biological priors with mRNA network: STRING PPI, DoRothEA, cBioPortal._"
            )
            lines.append("")

        # Detailed per-candidate breakdown
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

        discordance_flags = synthesis.get("discordance_flags", [])
        if discordance_flags:
            lines.append("### Notable Discordances")
            for flag in discordance_flags:
                lines.append(f"- **{flag['description']}**")
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
        if routing == "signature":
            return await self._call_llm_signature_synthesis(synthesis)
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

    async def _call_llm_signature_synthesis(self, synthesis: dict) -> str:
        cell_type = synthesis.get("cell_type", "unknown")
        signature_size = synthesis.get("signature_size", 0)
        genes_found = synthesis.get("genes_found_in_network", 0)
        ranked_drivers = synthesis.get("ranked_drivers") or []

        driver_lines = []
        for d in ranked_drivers[:5]:
            findings_str = "; ".join((d.get("cascade_key_findings") or [])[:2]) or "no CASCADE data"
            driver_lines.append(
                f"- {d['gene']}: {d['coverage_pct']}% coverage "
                f"({d['overlap_count']} genes), enrichment={d['enrichment_score']:.2f}, "
                f"corroboration={d['corroboration_count']}/7, CASCADE: {findings_str}"
            )

        prompt = (
            f"Gene signature driver analysis in {cell_type}.\n"
            f"Signature: {signature_size} genes ({genes_found} found in network).\n\n"
            "Top TF drivers (ranked by signature coverage + cross-system corroboration):\n"
            + "\n".join(driver_lines)
            + "\n\nWrite a concise 2-3 sentence biological interpretation identifying "
            "which transcription factors most likely drive this gene signature, "
            "which have the strongest cross-system evidence (both network enrichment "
            "and CASCADE experimental support), and what this suggests about the "
            "regulatory mechanism. Be specific."
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
        gene_signature: Optional[list] = None,
        cell_types: Optional[list] = None,
        progress: Optional[Callable] = None,
    ) -> dict:
        """Run a full Orchestra analysis. Opens MCP client connections for the duration."""
        token = _progress_cb.set(progress)
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
            gene_signature=gene_signature,
            master_regulators=None,
            cell_types=cell_types,
            comparison_results=None,
            synthesis=None,
            completed_steps=[],
            errors={},
            final_report=None,
        )

        try:
            if self._persistent_cascade is not None and self._persistent_regnetagents is not None:
                self._cascade = self._persistent_cascade
                self._regnetagents = self._persistent_regnetagents
                try:
                    result = await self.graph.ainvoke(initial_state)
                finally:
                    self._cascade = None
                    self._regnetagents = None
            else:
                async with make_cascade_client() as cascade, make_regnetagents_client() as regnetagents:
                    self._cascade = cascade
                    self._regnetagents = regnetagents
                    try:
                        result = await self.graph.ainvoke(initial_state)
                    finally:
                        self._cascade = None
                        self._regnetagents = None
        finally:
            _progress_cb.reset(token)

        return result
