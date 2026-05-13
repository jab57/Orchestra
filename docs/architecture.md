# Architecture

## Overview

Orchestra is simultaneously an MCP server (to Claude Desktop or any MCP client) and an MCP client (to RegNetAgents and CASCADE). A LangGraph DAG coordinates the full workflow.

```
Claude Desktop
      │
      ▼
Orchestra MCP Server  (orchestra_mcp_server.py)
      │
      ▼
OrchestraWorkflow  (orchestra_langgraph_workflow.py)
  [LangGraph DAG]
      │                    │
      ▼                    ▼
RegNetAgents MCP      CASCADE MCP
(network analysis,    (perturbation sim,
 pathway enrichment,   PPI, LINCS, DepMap,
 ARACNe/GREmLN)        super-enhancers)
```

## Three-Layer Design

| Layer | Node | What it does |
|---|---|---|
| **Decision** | `classify_gene`, `route_analysis` | Calls CASCADE `get_gene_metadata`; routes to TF, effector, or validation path |
| **Evidence** | `run_tf_path`, `run_effector_path`, `run_validation_path` | Parallel MCP calls to both child servers |
| **Explanation** | `synthesize`, `generate_report` | Cross-system corroboration counting; formatted report |

## LangGraph DAG

```
initialize
    │
classify_gene          ← CASCADE: get_gene_metadata
    │
route_analysis
    │
    ├── tf_path        → run_tf_path         ← parallel: RegNetAgents + CASCADE
    ├── effector_path  → run_effector_path   ← CASCADE PPI → TF partner → parallel
    └── validation_path → run_validation_path ← 3-source candidate extraction → validate
    │
synthesize             ← builds corroboration table
    │
generate_report        ← formats Markdown + optional LLM narrative
    │
END
```

## Routing Logic

`_routing_decision` checks two conditions in order:

1. If `analysis_type == "therapeutic_validation"` → validation path
2. Otherwise, check `gene_role` from CASCADE classification:
   - `master_regulator`, `transcription_factor`, `minor_regulator` → TF path
   - `effector`, `isolated`, or unknown → effector path

## MCP Client Pattern

`MCPClient` (`mcp_client.py`) manages the full lifecycle of one child server connection:

- **Startup:** spawns child server as a subprocess using the child server's own venv Python (`env/Scripts/python.exe`), not Orchestra's. This keeps NumPy and other dependencies isolated.
- **Transport:** MCP stdio — the child server reads/writes JSON-RPC over stdin/stdout.
- **Session:** `mcp.ClientSession` handles `initialize()` and protocol handshake.
- **Teardown:** `AsyncExitStack` ensures the subprocess terminates cleanly on exit or exception.

```python
async with make_cascade_client() as cascade:
    result = await cascade.call_tool(
        "get_gene_metadata",
        {"gene": "APC", "cell_type": "epithelial_cell"},
        timeout_seconds=30.0,
    )
```

`call_tool` prefers `structuredContent` (mcp >1.9) and falls back to JSON-parsing the text content, then returns a plain-text dict as last resort.

Per-tool timeouts enforce latency budgets:
- `TIMEOUT_PERTURBATION = 60s` — comprehensive_perturbation_analysis
- `TIMEOUT_NETWORK = 60s` — comprehensive_gene_analysis, pathway_focused_analysis
- `TIMEOUT_PPI = 15s` — get_protein_interactions

## TF Path

For TFs and master regulators, Orchestra runs two calls in parallel:

```python
rna_result, cascade_result = await asyncio.gather(
    regnetagents.call_tool("comprehensive_gene_analysis", ...),
    cascade.call_tool("comprehensive_perturbation_analysis", ...),
    return_exceptions=True,
)
```

Synthesis (`_synthesize_tf_path`) finds genes that appear in both:
- RegNetAgents `target_analysis.cascade_targets` (network topology targets)
- CASCADE `evidence_synthesis.multi_source_genes` (experimentally corroborated genes)

These **cross-system hits** are the primary Orchestra output: genes independently supported by regulatory network inference and experimental perturbation data.

## Effector Path

Scaffold/effector proteins have no transcriptional targets in the ARACNe network. A direct perturbation query returns empty results. Orchestra handles this automatically:

1. `get_protein_interactions(gene)` — finds STRING PPI partners
2. `_find_tf_partner()` — classifies top 10 PPI partners in parallel; returns the TF with the most downstream targets (heuristic: proxy for network influence)
3. Runs TF-path analysis on the TF partner: `comprehensive_perturbation_analysis(tf_partner)` + `comprehensive_gene_analysis(tf_partner)`

**Note on TF partner selection:** `_find_tf_partner` uses downstream target count as a proxy for network influence. This works reliably for canonical scaffold genes like APC (where CTNNB1 dominates unambiguously) but may require domain knowledge for scaffold proteins with multiple competing high-target TF partners.

## Validation Path

Identifies and ranks therapeutic targets using three independent candidate sources:

1. **RegNetAgents PageRank** (`comprehensive_gene_analysis`) — TF regulators ranked by centrality
2. **CASCADE drug discovery** (`therapeutic_target_discovery`) — super-enhancer, PPI, and drug database targets
3. **CASCADE STRING PPI** (`get_protein_interactions`) — protein-level interactors

Top 3 unique candidates are validated in parallel via `comprehensive_perturbation_analysis`. The synthesis layer scores each candidate against 7 independent evidence sources.

## Synthesis Layer

The corroboration table is Orchestra's core contribution. Before any LLM call, `_score_candidate_evidence` checks each candidate against 7 independent evidence sources:

| Source | System | What it captures |
|---|---|---|
| PageRank rank | RegNetAgents | Topological centrality in ARACNe/GREmLN network |
| Pathway membership | RegNetAgents | Reactome pathway enrichment membership |
| LINCS knockdown | CASCADE | Experimental knockdown → transcriptional response |
| DepMap essentiality | CASCADE | CRISPR fitness dependency in cancer cell lines |
| Super-enhancer | CASCADE | BET inhibitor sensitivity via super-enhancer analysis |
| DoRothEA tier | CASCADE | Curated TF-regulon confidence (literature + ChIP-seq) |
| cBioPortal expression | CASCADE | Primary tumor expression from TCGA/ICGC cohorts |

Corroboration count = number of sources that support the candidate. Sources are checked from key_findings text and candidate metadata, not raw tool outputs.

**Important:** RegNetAgents networks include both healthy and disease/cancer-infiltrating cells (GREmLN corpus, Zhang et al. 2025). CASCADE LINCS/DepMap are cancer cell lines. Agreement between systems is methodologically independent but does not imply matched biological context.

## Graceful Degradation

All evidence calls use `asyncio.gather(..., return_exceptions=True)`. Exceptions are stored in `state["errors"]` and the workflow continues with partial data. Report formatters check `regnetagents_available` and `cascade_available` flags and prepend a warning banner when a system is absent.

## CASCADE Output Contract

CASCADE's `comprehensive_perturbation_analysis` returns a pre-synthesized `evidence_synthesis` block:

```python
cascade_output["evidence_synthesis"]["key_findings"]    # list of strings
cascade_output["evidence_synthesis"]["multi_source_genes"]  # list of {symbol, source_count, sources}
cascade_output["evidence_synthesis"]["source_agreements"]
cascade_output["evidence_synthesis"]["source_disagreements"]
```

Orchestra consumes this block directly. It does not re-implement within-CASCADE cross-source referencing — that is CASCADE's responsibility. Orchestra synthesis is exclusively cross-system: RegNetAgents network topology + CASCADE's already-synthesized experimental evidence.
