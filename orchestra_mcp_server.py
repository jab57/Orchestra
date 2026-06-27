"""
Orchestra MCP Server

Exposes Orchestra as an MCP server to Claude Desktop and other MCP clients,
while acting as an MCP client to RegNetAgents and CASCADE child servers.

Ten composite tools:
  causal_chain_analysis             — TF path (parallel RegNetAgents + CASCADE) or
                                      effector path (PPI → TF partner → simulate)
  validate_therapeutic_targets      — PageRank + drug discovery + PPI → 7-source corroboration table
  effector_analysis                 — scaffold/effector routing (APC→CTNNB1 pattern)
  analyze_gene_signature            — DEG list → ranked TF drivers (Fisher enrichment + CASCADE validation)
  compare_cell_contexts             — 7-source evidence heatmap across N cell types (Issue #11)
  compare_network_contexts          — GREmLN vs TCGA regulatory rewiring + CASCADE validation (Issue #13)
  compare_tumor_networks            — tumor-vs-tumor cross-cancer convergence analysis (GitHub #14)
  novelty_assessment                — PubMed hit count + novelty verdict for a gene in a cancer context (Issue #15)
  novelty_assessment_batch          — atomic batch variant: novelty_assessment for N genes in one call
  compare_network_contexts_batch    — atomic batch variant: compare_network_contexts for N genes in one call
  fetch_tcga_methylation_correlation — Spearman correlation between regulator RNA-seq expression and
                                       target gene methylation beta values across a TCGA cohort (cBioPortal)
"""

import asyncio
from contextlib import AsyncExitStack
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent, Prompt, PromptArgument, GetPromptResult, PromptMessage
from orchestra_langgraph_workflow import OrchestraWorkflow
from mcp_client import make_cascade_client, make_regnetagents_client

app = Server("orchestra")
workflow = OrchestraWorkflow()


@app.list_tools()
async def list_tools() -> list[Tool]:
    return [
        Tool(
            name="causal_chain_analysis",
            description=(
                "Full causal chain analysis: classifies gene, runs perturbation "
                "simulation (CASCADE) and regulatory network analysis (RegNetAgents) "
                "in parallel, synthesizes results into an integrated causal report."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {"type": "string", "description": "Gene symbol (e.g. TP53)"},
                    "cell_type": {"type": "string", "description": "Cell type (e.g. epithelial_cell)"},
                    "analysis_depth": {
                        "type": "string",
                        "enum": ["basic", "comprehensive"],
                        "default": "comprehensive"
                    },
                    "cancer_context": {
                        "type": "string",
                        "description": "Plain-text cancer context for PubMed pair-novelty queries (e.g. 'cervical cancer', 'colorectal'). Optional — omit to skip pair novelty.",
                    },
                },
                "required": ["gene", "cell_type"],
            },
        ),
        Tool(
            name="validate_therapeutic_targets",
            description=(
                "Rank upstream regulators by PageRank (RegNetAgents) then validate "
                "top candidates via perturbation simulation and LINCS experimental "
                "knockdown data (CASCADE)."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {"type": "string"},
                    "cell_type": {"type": "string"},
                    "cancer_context": {
                        "type": "string",
                        "description": "Plain-text cancer context for PubMed pair-novelty queries (e.g. 'cervical cancer', 'colorectal'). Optional — omit to skip pair novelty.",
                    },
                },
                "required": ["gene", "cell_type"],
            },
        ),
        Tool(
            name="effector_analysis",
            description=(
                "Automated analysis for effector/scaffold genes with no transcriptional "
                "targets (e.g. APC). Detects effector role, finds TF partners via PPI, "
                "simulates the TF partner, and enriches the cascade against pathways."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {"type": "string"},
                    "cell_type": {"type": "string"},
                },
                "required": ["gene", "cell_type"],
            },
        ),
        Tool(
            name="analyze_gene_signature",
            description=(
                "Identify which transcription factors are most likely driving a gene "
                "signature (e.g. a list of differentially expressed genes). "
                "Uses RegNetAgents Fisher enrichment to rank TFs by regulon overlap with "
                "the input gene set, then validates top candidates via CASCADE perturbation "
                "simulation and experimental data (LINCS, DepMap, super-enhancers). "
                "Returns a ranked driver table with signature coverage % and 7-source "
                "corroboration score. Optional: supply cancer_contexts to add a "
                "Cross-Context Novelty Gap table classifying each driver as "
                "transfer_opportunity, bilateral_novel, or bilateral_established."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "genes": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of gene symbols (e.g. differentially expressed genes)",
                        "minItems": 2,
                    },
                    "cell_type": {
                        "type": "string",
                        "description": "Cell type context for network and perturbation analysis (e.g. epithelial_cell)",
                    },
                    "cancer_contexts": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Optional list of cancer contexts for cross-context novelty gap analysis (e.g. [\"breast cancer\", \"cervical cancer\"]). When supplied, PubMed novelty is queried for the top 5 ranked drivers in each context and a gap table is appended to the report.",
                    },
                    "tcga_network": {
                        "type": "string",
                        "description": "Optional TCGA tumor network to use for master regulator enrichment instead of GREmLN (e.g. 'cesc' for cervical cancer). Supported: blca, brca, cesc, coad, hnsc, kirc, lihc, luad, lusc, ov, paad, prad, stad, ucec. When omitted, uses the population-averaged GREmLN network.",
                    },
                },
                "required": ["genes", "cell_type"],
            },
        ),
        Tool(
            name="compare_cell_contexts",
            description=(
                "Compare a gene's regulatory evidence across multiple cell types. "
                "Runs RegNetAgents network analysis and CASCADE perturbation analysis "
                "for each cell type in parallel, then produces a 7-source evidence heatmap "
                "showing which findings are conserved across cell types vs. cell-type-specific. "
                "Use this to determine tissue specificity of a regulatory relationship before "
                "selecting a therapeutic context."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Gene symbol to compare across cell types (e.g. MYC)",
                    },
                    "cell_types": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": (
                            "List of cell type contexts to compare (e.g. [\"epithelial_cell\", \"cd4_t_cells\", \"nk_cells\"]). "
                            "Available: cd4_t_cells, cd8_t_cells, cd14_monocytes, cd16_monocytes, "
                            "nk_cells, nkt_cells, cd20_b_cells, monocyte-derived_dendritic_cells, "
                            "erythrocytes, epithelial_cell"
                        ),
                        "minItems": 2,
                    },
                },
                "required": ["gene", "cell_types"],
            },
        ),
        Tool(
            name="novelty_assessment",
            description=(
                "Query PubMed for a gene (or gene pair) in a cancer context and return "
                "a structured novelty verdict: established (>20 papers), emerging (5–20), "
                "or novel (<5). Reports total hit count split into experimental vs. computational "
                "papers, and the most recent publication year. Use this to gauge how well-characterized "
                "a finding is before writing it up, or to prioritize results from other Orchestra tools."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Primary gene symbol (e.g. FOXM1)",
                    },
                    "cancer_context": {
                        "type": "string",
                        "description": "Plain-text cancer context for the PubMed query (e.g. 'head and neck squamous', 'breast cancer', 'colorectal')",
                    },
                    "gene2": {
                        "type": "string",
                        "description": "Optional second gene for gene-pair queries (e.g. TOP2A). When provided, the query requires both genes to appear in the abstract.",
                    },
                },
                "required": ["gene", "cancer_context"],
            },
        ),
        Tool(
            name="compare_network_contexts",
            description=(
                "Compare a gene's regulatory wiring between population-averaged (GREmLN ARACNe) "
                "and tumor-state (TCGA ARACNe) networks, then validate conserved regulators via "
                "CASCADE experimental data (LINCS, DepMap, super-enhancers, DoRothEA). "
                "Returns three regulator tiers: conserved + CASCADE-validated (highest confidence "
                "therapeutic candidates), conserved without CASCADE support (regulatory inference), "
                "and tumor-state-only (emerging in cancer, not present in population-averaged wiring). "
                "Rewiring classification (low/moderate/high) quantifies how different the tumor "
                "regulatory program is from population-averaged. "
                "TCGA cancer types available: blca, brca, cesc, coad, hnsc, kirc, lihc, luad, lusc, ov, paad, prad, stad, ucec."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Gene symbol (e.g. FOXM1, STAT3, MYC)",
                    },
                    "cancer_type": {
                        "type": "string",
                        "enum": ["blca", "brca", "cesc", "coad", "hnsc", "kirc", "lihc", "luad", "lusc", "ov", "paad", "prad", "stad", "ucec"],
                        "description": (
                            "TCGA cancer type for the tumor-state network. "
                            "blca=bladder, brca=breast, cesc=cervical, coad=colon, "
                            "hnsc=head/neck squamous, kirc=kidney, lihc=liver, "
                            "luad=lung adenocarcinoma, lusc=lung squamous, ov=ovarian, "
                            "paad=pancreatic, prad=prostate, stad=stomach, ucec=uterine."
                        ),
                    },
                    "cell_type": {
                        "type": "string",
                        "description": (
                            "GREmLN population-averaged cell type for the reference network "
                            "(default: epithelial_cell). "
                            "Available: cd4_t_cells, cd8_t_cells, cd14_monocytes, cd16_monocytes, "
                            "nk_cells, nkt_cells, cd20_b_cells, monocyte-derived_dendritic_cells, "
                            "erythrocytes, epithelial_cell"
                        ),
                        "default": "epithelial_cell",
                    },
                },
                "required": ["gene", "cancer_type"],
            },
        ),
        Tool(
            name="novelty_assessment_batch",
            description=(
                "Run PubMed novelty assessment for a list of genes in a single atomic call. "
                "Returns a summary table (verdict, hit counts, experimental vs. computational "
                "breakdown, most-recent-year) for every gene. "
                "Use this instead of calling novelty_assessment repeatedly when assessing "
                "multiple candidates — prevents silent truncation of the candidate list."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "genes": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Gene symbols to assess (e.g. [\"LITAF\", \"USP21\", \"CENPK\"])",
                        "minItems": 2,
                    },
                    "cancer_context": {
                        "type": "string",
                        "description": "Plain-text cancer context for all queries (e.g. 'cervical cancer')",
                    },
                },
                "required": ["genes", "cancer_context"],
            },
        ),
        Tool(
            name="compare_network_contexts_batch",
            description=(
                "Run GREmLN vs TCGA network context comparison for a list of genes in a "
                "single atomic call. Executes compare_network_contexts sequentially for each "
                "gene and returns combined reports. "
                "Use this instead of calling compare_network_contexts repeatedly when "
                "comparing multiple candidate genes in the same cancer type — prevents "
                "silent truncation of the candidate list."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "genes": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Gene symbols to compare (e.g. [\"FOXM1\", \"MYC\", \"TP53\"])",
                        "minItems": 2,
                    },
                    "cancer_type": {
                        "type": "string",
                        "enum": ["blca", "brca", "cesc", "coad", "hnsc", "kirc", "lihc", "luad", "lusc", "ov", "paad", "prad", "stad", "ucec"],
                        "description": "TCGA cancer type applied to all genes in the batch.",
                    },
                    "cell_type": {
                        "type": "string",
                        "description": "GREmLN reference cell type for all genes (default: epithelial_cell).",
                        "default": "epithelial_cell",
                    },
                },
                "required": ["genes", "cancer_type"],
            },
        ),
        Tool(
            name="compare_tumor_networks",
            description=(
                "Tumor-vs-tumor cross-cancer convergence analysis: compare a gene's TCGA tumor "
                "regulatory network across 2–4 cancer types directly against each other "
                "(not against GREmLN normal). Returns pairwise regulator overlap, convergent "
                "core regulators (present in ALL tested cancer types), cancer-type-specific "
                "divergent regulators, CASCADE validation on convergent core, PubMed pair "
                "novelty for convergent regulators, and a convergent/divergent/mixed verdict "
                "with explicit thresholds. Optionally includes the GREmLN-vs-tumor baseline "
                "rewiring stats for each cancer type as context. "
                "Use Pipeline 9 when the question is whether a gene's rewiring is shared "
                "between specific cancer types (tumor-vs-tumor). Use compare_network_contexts "
                "(Pipeline 8) when the question is how much rewiring occurred vs. normal in "
                "each type separately. These pipelines answer different questions — do not "
                "conflate them."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Gene symbol (e.g. FOXM1, MYC, TP53)",
                    },
                    "cancer_types": {
                        "type": "array",
                        "items": {
                            "type": "string",
                            "enum": ["blca", "brca", "cesc", "coad", "hnsc", "kirc", "lihc", "luad", "lusc", "ov", "paad", "prad", "stad", "ucec"],
                        },
                        "description": "2–4 TCGA cancer type codes to compare (e.g. [\"cesc\", \"hnsc\", \"luad\"]). Limit 4 — each call takes 2–3 minutes.",
                        "minItems": 2,
                        "maxItems": 4,
                    },
                    "cell_type": {
                        "type": "string",
                        "description": "GREmLN reference cell type used for the GREmLN baseline comparison in each cancer type call (default: epithelial_cell).",
                        "default": "epithelial_cell",
                    },
                    "include_gremln_baseline": {
                        "type": "boolean",
                        "description": "Include GREmLN-vs-tumor rewiring stats alongside tumor-vs-tumor results (default: true). Useful context: high GREmLN rewiring + high pairwise tumor overlap = convergent oncogenic rewiring, not baseline variation.",
                        "default": True,
                    },
                },
                "required": ["gene", "cancer_types"],
            },
        ),
        Tool(
            name="fetch_tcga_methylation_correlation",
            description=(
                "Compute Spearman correlation between a regulator's RNA-seq expression and "
                "target gene promoter methylation beta values across a TCGA tumour cohort "
                "(via cBioPortal). Use this to test whether a regulator's activity is "
                "associated with epigenetic silencing of its predicted targets: a negative "
                "ρ (high expression → low methylation) supports the regulatory hypothesis. "
                "Returns a per-target table of ρ, p-value, sample count, and direction."
            ),
            inputSchema={
                "type": "object",
                "properties": {
                    "regulator": {
                        "type": "string",
                        "description": "Gene whose RNA-seq expression to use (e.g. 'FOS')",
                    },
                    "target_genes": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Genes whose methylation beta values to correlate against regulator expression (e.g. ['CDKN2A', 'SOCS1', 'TIMP3'])",
                        "minItems": 1,
                    },
                    "tcga_network": {
                        "type": "string",
                        "enum": ["blca", "brca", "cesc", "coad", "hnsc", "kirc", "lihc", "luad", "lusc", "ov", "paad", "prad", "stad", "ucec"],
                        "description": "TCGA cancer cohort code (e.g. 'cesc' for cervical cancer)",
                    },
                },
                "required": ["regulator", "target_genes", "tcga_network"],
            },
        ),
    ]


@app.list_prompts()
async def list_prompts() -> list[Prompt]:
    return [
        Prompt(
            name="biomarker_discovery",
            description=(
                "Cross-cancer regulatory biomarker discovery: identifies upstream TF drivers "
                "of a silenced gene panel and finds which drivers are established in one cancer "
                "context but novel in another. Use this instead of analyzing genes individually."
            ),
            arguments=[
                PromptArgument(
                    name="cancer_context_1",
                    description="First cancer context for novelty comparison (e.g. 'breast cancer')",
                    required=True,
                ),
                PromptArgument(
                    name="cancer_context_2",
                    description="Second cancer context for novelty comparison (e.g. 'cervical cancer')",
                    required=True,
                ),
            ],
        ),
    ]


@app.get_prompt()
async def get_prompt(name: str, arguments: dict | None) -> GetPromptResult:
    if name == "biomarker_discovery":
        ctx1 = (arguments or {}).get("cancer_context_1", "cancer context 1")
        ctx2 = (arguments or {}).get("cancer_context_2", "cancer context 2")
        return GetPromptResult(
            description="Cross-cancer regulatory biomarker discovery workflow",
            messages=[
                PromptMessage(
                    role="user",
                    content=TextContent(
                        type="text",
                        text=(
                            f"I want to identify upstream transcription factor drivers of a "
                            f"silenced gene panel and find which drivers represent transfer "
                            f"opportunities between {ctx1} and {ctx2}.\n\n"
                            f"Please follow these steps in order:\n\n"
                            f"**Step 1 — Compile the gene panel (do not skip)**\n"
                            f"Before calling any Orchestra tool, compile a literature-curated list "
                            f"of 20–30 genes recurrently silenced by promoter hypermethylation in "
                            f"both {ctx1} and {ctx2}. Present the list with a one-line justification "
                            f"for each gene. Wait for my approval before proceeding.\n\n"
                            f"**Step 2 — Run signature enrichment**\n"
                            f"Call `analyze_gene_signature` with the approved gene list and "
                            f"`cell_type='epithelial_cell'`. Do NOT run `causal_chain_analysis` "
                            f"or `comprehensive_gene_analysis` on individual genes — that bypasses "
                            f"the Fisher enrichment and produces unstatistical results.\n\n"
                            f"**Step 3 — Cross-context novelty**\n"
                            f"For the top 5 ranked TF drivers, call `novelty_assessment_batch` "
                            f"twice — once with cancer_context='{ctx1}' and once with "
                            f"cancer_context='{ctx2}'.\n\n"
                            f"**Step 4 — Synthesize the gap table**\n"
                            f"Present a table classifying each driver as:\n"
                            f"- Transfer opportunity: established in one context (>5 papers), novel (<5) in the other\n"
                            f"- Bilateral novel: novel in both\n"
                            f"- Bilateral established: established in both\n\n"
                            f"Prioritize transfer opportunities for the proposal."
                        ),
                    ),
                ),
            ],
        )
    raise ValueError(f"Unknown prompt: {name}")


@app.call_tool()
async def call_tool(name: str, arguments: dict) -> list[TextContent]:
    # Capture session once at request start so the callback is safe to call
    # from asyncio sub-tasks (gather) without re-entering the ContextVar lookup.
    _session = app.request_context.session
    _request_id = app.request_context.request_id

    async def progress(msg: str) -> None:
        try:
            await _session.send_log_message(
                level="info",
                data=msg,
                logger="orchestra",
                related_request_id=_request_id,
            )
        except Exception:
            pass

    cell_type = arguments.get("cell_type", "")

    if name == "novelty_assessment":
        gene = arguments.get("gene", "")
        cancer_context = arguments.get("cancer_context", "")
        gene2 = arguments.get("gene2") or None
        result = await workflow.run_analysis(
            gene=gene,
            cell_type="",
            analysis_type="novelty_assessment",
            cancer_context=cancer_context,
            gene2=gene2,
            progress=progress,
        )
        subject = f"{gene}/{gene2}" if gene2 else gene
        label = f"{subject} in {cancer_context}"
    elif name == "compare_network_contexts":
        gene = arguments.get("gene", "")
        cancer_type = arguments.get("cancer_type", "")
        cell_type = arguments.get("cell_type", "epithelial_cell")
        result = await workflow.run_analysis(
            gene=gene,
            cell_type=cell_type,
            analysis_type="network_comparison",
            cancer_type=cancer_type,
            progress=progress,
        )
        label = f"{gene}: {cell_type} (GREmLN) vs TCGA {cancer_type.upper()}"
    elif name == "novelty_assessment_batch":
        from pubmed_client import novelty_assessment as _pubmed_novelty
        genes = arguments.get("genes", [])
        cancer_context = arguments.get("cancer_context", "")
        await progress(
            f"[Orchestra] Batch novelty assessment: {len(genes)} genes in {cancer_context}..."
        )
        tasks = [_pubmed_novelty(g, cancer_context) for g in genes]
        results = await asyncio.gather(*tasks, return_exceptions=True)
        rows = [
            "| Gene | Verdict | PubMed Hits | Experimental | Computational | Last Year |",
            "|------|---------|-------------|--------------|---------------|-----------|",
        ]
        rationales = []
        for gene, res in zip(genes, results):
            if isinstance(res, BaseException):
                rows.append(f"| {gene} | ERROR | — | — | — | — |")
                rationales.append(f"**{gene}**: error — {res}")
            else:
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
        label = f"{len(genes)} genes in {cancer_context}"
        return [TextContent(type="text", text=f"Orchestra batch novelty assessment for {label}:\n\n{report}")]
    elif name == "compare_network_contexts_batch":
        genes = arguments.get("genes", [])
        cancer_type = arguments.get("cancer_type", "")
        cell_type = arguments.get("cell_type", "epithelial_cell")
        sections = []
        for gene in genes:
            await progress(
                f"[Orchestra] Comparing network contexts: {gene} "
                f"({cell_type} vs TCGA {cancer_type.upper()})..."
            )
            res = await workflow.run_analysis(
                gene=gene,
                cell_type=cell_type,
                analysis_type="network_comparison",
                cancer_type=cancer_type,
                progress=progress,
            )
            sections.append(res.get("final_report", f"_{gene}: analysis failed_"))
        header = "\n".join([
            f"## Network Context Comparison Batch: TCGA {cancer_type.upper()}",
            f"**Genes:** {', '.join(genes)}",
            f"**Reference network:** {cell_type} (GREmLN)",
            "",
        ])
        report = header + "\n\n---\n\n".join(sections)
        label = f"{', '.join(genes)}: {cell_type} vs TCGA {cancer_type.upper()}"
        return [TextContent(type="text", text=f"Orchestra batch network comparison for {label}:\n\n{report}")]
    elif name == "compare_tumor_networks":
        gene = arguments.get("gene", "")
        cancer_types = arguments.get("cancer_types", [])
        cell_type = arguments.get("cell_type", "epithelial_cell")
        include_gremln_baseline = arguments.get("include_gremln_baseline", True)
        result = await workflow.run_analysis(
            gene=gene,
            cell_type=cell_type,
            analysis_type="tumor_network_comparison",
            cancer_types=cancer_types,
            include_gremln_baseline=include_gremln_baseline,
            progress=progress,
        )
        label = f"{gene}: {', '.join(ct.upper() for ct in cancer_types)} tumor-vs-tumor"
    elif name == "analyze_gene_signature":
        genes = arguments.get("genes", [])
        cancer_contexts = arguments.get("cancer_contexts") or None
        tcga_network = arguments.get("tcga_network") or None
        result = await workflow.run_analysis(
            gene="",
            cell_type=cell_type,
            analysis_type="gene_signature",
            gene_signature=genes,
            cancer_contexts=cancer_contexts,
            cancer_type=tcga_network,
            progress=progress,
        )
        label = f"Gene signature ({len(genes)} genes) in {cell_type}"
    elif name == "compare_cell_contexts":
        gene = arguments.get("gene", "")
        cell_types = arguments.get("cell_types", [])
        result = await workflow.run_analysis(
            gene=gene,
            cell_type="",
            analysis_type="cell_context_comparison",
            cell_types=cell_types,
            progress=progress,
        )
        label = f"{gene} across {', '.join(cell_types)}"
    elif name == "fetch_tcga_methylation_correlation":
        from cbioportal_client import (
            methylation_expression_correlation as _meth_corr,
            format_correlation_report as _fmt_corr,
        )
        regulator = arguments.get("regulator", "")
        target_genes = arguments.get("target_genes", [])
        tcga_network = arguments.get("tcga_network", "")
        await progress(
            f"[Orchestra] Fetching TCGA methylation-expression correlation: "
            f"{regulator} vs {', '.join(target_genes)} in {tcga_network.upper()}..."
        )
        try:
            corr_result = await asyncio.wait_for(
                _meth_corr(regulator, target_genes, tcga_network),
                timeout=120.0,
            )
        except asyncio.TimeoutError:
            corr_result = {
                "error": (
                    "cBioPortal request timed out after 120 s — the public API may be "
                    "under load. Please retry in a moment."
                )
            }
        report = _fmt_corr(corr_result)
        return [TextContent(type="text", text=report)]
    else:
        gene = arguments.get("gene", "")
        depth = arguments.get("analysis_depth", "comprehensive")
        cancer_context = arguments.get("cancer_context") or None
        result = await workflow.run_analysis(
            gene=gene,
            cell_type=cell_type,
            analysis_type=name,
            analysis_depth=depth,
            cancer_context=cancer_context,
            progress=progress,
        )
        label = f"{gene} in {cell_type}"

    return [TextContent(
        type="text",
        text=f"Orchestra analysis for {label}: {result.get('final_report', 'In progress — implementation pending.')}"
    )]


async def main():
    async with AsyncExitStack() as stack:
        # Start stdio_server FIRST so Claude Desktop can connect immediately.
        # Persistent connections are opened in a background task so they never
        # block Orchestra's own MCP handshake.
        async with stdio_server() as (read_stream, write_stream):

            async def _init_persistent() -> None:
                try:
                    cascade = await stack.enter_async_context(make_cascade_client())
                    regnetagents = await stack.enter_async_context(make_regnetagents_client())
                    workflow._persistent_cascade = cascade
                    workflow._persistent_regnetagents = regnetagents
                    workflow._persistent_ready.set()
                except asyncio.CancelledError:
                    pass
                except Exception:
                    workflow._persistent_ready.set()  # unblock waiters even on failure

            init_task = asyncio.create_task(_init_persistent())
            try:
                await app.run(read_stream, write_stream, app.create_initialization_options())
            finally:
                init_task.cancel()

    workflow._persistent_cascade = None
    workflow._persistent_regnetagents = None


if __name__ == "__main__":
    asyncio.run(main())
