"""
Orchestra MCP Server

Exposes Orchestra as an MCP server to Claude Desktop and other MCP clients,
while acting as an MCP client to RegNetAgents and CASCADE child servers.

Four composite tools:
  causal_chain_analysis        — TF path (parallel RegNetAgents + CASCADE) or
                                 effector path (PPI → TF partner → simulate)
  validate_therapeutic_targets — PageRank + drug discovery + PPI → 7-source corroboration table
  effector_analysis            — scaffold/effector routing (APC→CTNNB1 pattern)
  analyze_gene_signature       — DEG list → ranked TF drivers (Fisher enrichment + CASCADE validation)
"""

import asyncio
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent
from orchestra_langgraph_workflow import OrchestraWorkflow

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
    ]


@app.call_tool()
async def call_tool(name: str, arguments: dict) -> list[TextContent]:
    cell_type = arguments.get("cell_type", "")

    if name == "analyze_gene_signature":
        genes = arguments.get("genes", [])
        result = await workflow.run_analysis(
            gene="",
            cell_type=cell_type,
            analysis_type="gene_signature",
            gene_signature=genes,
        )
        label = f"Gene signature ({len(genes)} genes) in {cell_type}"
    else:
        gene = arguments.get("gene", "")
        depth = arguments.get("analysis_depth", "comprehensive")
        result = await workflow.run_analysis(
            gene=gene,
            cell_type=cell_type,
            analysis_type=name,
            analysis_depth=depth,
        )
        label = f"{gene} in {cell_type}"

    return [TextContent(
        type="text",
        text=f"Orchestra analysis for {label}: {result.get('final_report', 'In progress — implementation pending.')}"
    )]


async def main():
    async with stdio_server() as (read_stream, write_stream):
        await app.run(read_stream, write_stream, app.create_initialization_options())


if __name__ == "__main__":
    asyncio.run(main())
