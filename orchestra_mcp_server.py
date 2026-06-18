"""
Orchestra MCP Server

Exposes Orchestra as an MCP server to Claude Desktop and other MCP clients,
while acting as an MCP client to RegNetAgents and CASCADE child servers.

Seven composite tools:
  causal_chain_analysis        — TF path (parallel RegNetAgents + CASCADE) or
                                 effector path (PPI → TF partner → simulate)
  validate_therapeutic_targets — PageRank + drug discovery + PPI → 7-source corroboration table
  effector_analysis            — scaffold/effector routing (APC→CTNNB1 pattern)
  analyze_gene_signature       — DEG list → ranked TF drivers (Fisher enrichment + CASCADE validation)
  compare_cell_contexts        — 7-source evidence heatmap across N cell types (Issue #11)
  compare_network_contexts     — GREmLN vs TCGA regulatory rewiring + CASCADE validation (Issue #13)
  novelty_assessment           — PubMed hit count + novelty verdict for a gene in a cancer context (Issue #15)
"""

import asyncio
from contextlib import AsyncExitStack
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent
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
                "corroboration score — a cross-system result neither RegNetAgents nor "
                "CASCADE can produce alone."
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
                "TCGA cancer types available: brca, coad, hnsc, luad, lusc, ov, prad, ucec."
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
                        "enum": ["brca", "coad", "hnsc", "luad", "lusc", "ov", "prad", "ucec"],
                        "description": (
                            "TCGA cancer type for the tumor-state network. "
                            "brca=breast, coad=colon, hnsc=head/neck squamous (HPV-associated, "
                            "closest proxy for cervical), luad=lung adenocarcinoma, "
                            "lusc=lung squamous, ov=ovarian, prad=prostate, ucec=uterine."
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
    ]


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
    elif name == "analyze_gene_signature":
        genes = arguments.get("genes", [])
        result = await workflow.run_analysis(
            gene="",
            cell_type=cell_type,
            analysis_type="gene_signature",
            gene_signature=genes,
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
    else:
        gene = arguments.get("gene", "")
        depth = arguments.get("analysis_depth", "comprehensive")
        result = await workflow.run_analysis(
            gene=gene,
            cell_type=cell_type,
            analysis_type=name,
            analysis_depth=depth,
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
