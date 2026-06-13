"""
Orchestra MCP Server

Exposes Orchestra as an MCP server to Claude Desktop and other MCP clients,
while acting as an MCP client to RegNetAgents and CASCADE child servers.

Five composite tools:
  causal_chain_analysis        — TF path (parallel RegNetAgents + CASCADE) or
                                 effector path (PPI → TF partner → simulate)
  validate_therapeutic_targets — PageRank + drug discovery + PPI → 7-source corroboration table
  effector_analysis            — scaffold/effector routing (APC→CTNNB1 pattern)
  analyze_gene_signature       — DEG list → ranked TF drivers (Fisher enrichment + CASCADE validation)
  compare_cell_contexts        — 7-source evidence heatmap across N cell types (Issue #11)
"""

import asyncio
from contextlib import AsyncExitStack
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent
from orchestra_langgraph_workflow import OrchestraWorkflow
from mcp_client import make_cascade_client, make_regnetagents_client, TIMEOUT_SERVER_WARMUP

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

    if name == "analyze_gene_signature":
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


async def _warmup_regnetagents(client) -> None:
    try:
        await client.call_tool(
            "query_network",
            {"gene": "TP53", "cell_type": "epithelial_cell"},
            timeout_seconds=TIMEOUT_SERVER_WARMUP,
        )
    except Exception:
        pass


async def main():
    async with AsyncExitStack() as stack:
        # Open persistent connections once — eliminates per-call cold starts (~60-90s each)
        cascade = await stack.enter_async_context(make_cascade_client())
        regnetagents = await stack.enter_async_context(make_regnetagents_client())
        workflow._persistent_cascade = cascade
        workflow._persistent_regnetagents = regnetagents

        # Pre-warm RegNetAgents in background so the network cache is ready before
        # the first tool call arrives. Cancelled cleanly when the server exits.
        warmup_task = asyncio.create_task(_warmup_regnetagents(regnetagents))
        stack.callback(warmup_task.cancel)

        async with stdio_server() as (read_stream, write_stream):
            await app.run(read_stream, write_stream, app.create_initialization_options())

    workflow._persistent_cascade = None
    workflow._persistent_regnetagents = None


if __name__ == "__main__":
    asyncio.run(main())
