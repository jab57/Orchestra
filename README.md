# Orchestra

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/jab57/Orchestra/actions/workflows/test.yml/badge.svg)](https://github.com/jab57/Orchestra/actions/workflows/test.yml)
[![Draft JOSS Paper](https://github.com/jab57/Orchestra/actions/workflows/draft-pdf.yml/badge.svg)](https://github.com/jab57/Orchestra/actions/workflows/draft-pdf.yml)

**MCP Orchestrator for Multi-System Causal Reasoning in Bioinformatics**

Orchestra automates multi-step causal reasoning across gene regulatory networks by composing two specialized MCP servers — [RegNetAgents](https://github.com/jab57/RegNetAgents) and [CASCADE](https://github.com/jab57/CASCADE) — via the Model Context Protocol. A single natural language query triggers an orchestrated LangGraph workflow that combines regulatory network analysis, perturbation simulation, protein interactions, and pathway enrichment into a unified causal report.

Neither system alone can answer: *"What transcription factors drive this gene signature, and what happens downstream if we inhibit the top candidate?"* Orchestra can.

## Quick Start

```bash
# 1. Clone and install
git clone https://github.com/jab57/Orchestra.git
cd Orchestra
python -m venv env
env\Scripts\activate       # Windows
# source env/bin/activate  # macOS/Linux
pip install -r requirements.txt

# 2. Configure paths to child servers
cp .env.example .env
# Edit .env: set REGNETAGENTS_SERVER_PATH and CASCADE_SERVER_PATH

# 3. Run a validation case
python run_validation.py apc
```

> **Prerequisites:** [RegNetAgents](https://github.com/jab57/RegNetAgents) and [CASCADE](https://github.com/jab57/CASCADE) must be installed separately, each with their own virtual environments.

## Architecture

Orchestra acts as both an MCP server (to Claude Desktop or any MCP client) and an MCP client (to RegNetAgents and CASCADE). A LangGraph DAG coordinates three layers:

- **Decision layer** — classifies gene type (TF vs. effector), routes to the appropriate analysis path
- **Evidence layer** — parallel MCP calls to RegNetAgents and CASCADE
- **Explanation layer** — cross-system corroboration counting; optional LLM narrative

```
Claude Desktop
      │
      ▼
Orchestra (MCP Server + LangGraph DAG)
      │                    │
      ▼                    ▼
RegNetAgents           CASCADE
(network analysis,     (perturbation sim,
 pathway enrichment,    PPI, LINCS, DepMap,
 domain agents)         super-enhancers,
                        DoRothEA, cBioPortal)
```

### Key Benefits of LangGraph Architecture

- **Intelligent Routing**: Automatically selects analysis strategy based on gene role (TF, effector, scaffold)
- **Parallel Execution**: RegNetAgents and CASCADE evidence calls run concurrently
- **Cross-System Synthesis**: Builds a corroboration count per candidate across 7 independent evidence sources
- **Graceful Degradation**: Falls back to single-system results if one child server is unavailable
- **Optional LLM Narrative**: Structured evidence always returned; LLM synthesis prepended when enabled

## Composite Tools

Orchestra exposes six tools — analytical capabilities that require both child servers and cannot be replicated by either alone:

### `causal_chain_analysis(gene, cell_type)`

Classifies a gene, runs regulatory network analysis (RegNetAgents) and perturbation simulation (CASCADE) in parallel, and synthesizes results into an integrated causal report. Two routing paths:

- **TF path** (master regulators, transcription factors): parallel RegNetAgents comprehensive analysis + CASCADE perturbation; identifies downstream genes corroborated by both network topology and experimental data
- **Effector path** (scaffold proteins, effectors): PPI → TF partner → simulate TF partner → pathway enrichment

### `validate_therapeutic_targets(gene, cell_type)`

Ranks upstream regulators by PageRank centrality (RegNetAgents), combines with drug target discovery (CASCADE super-enhancers, PPI), and validates top candidates against LINCS experimental knockdown data. Output: 7-source corroboration table.

### `effector_analysis(gene, cell_type)`

Handles scaffold/effector genes (e.g. APC) that have no direct transcriptional targets. Detects effector role, finds TF partners via protein-protein interactions (CASCADE), simulates the TF partner, and enriches the downstream cascade against Reactome pathways (RegNetAgents).

### `analyze_gene_signature(genes, cell_type)`

Identifies which transcription factors are most likely driving a list of differentially expressed genes. RegNetAgents ranks TFs by Fisher's exact test enrichment in the input gene set (ARACNe regulon overlap); CASCADE validates the top candidates with 7-source perturbation evidence. Output: ranked driver table with signature coverage % and cross-system corroboration count.

### `compare_cell_contexts(gene, cell_types)`

Compares a gene's regulatory evidence across multiple cell types. Runs RegNetAgents network analysis and CASCADE perturbation analysis for each cell type in parallel (2N total MCP calls), then produces a 7-source evidence heatmap classifying each source as conserved (≥ 2/3 of cell types), enriched, cell-type-specific, or absent. Use this to determine tissue specificity before selecting a therapeutic context.

### `compare_network_contexts(gene, cancer_type, cell_type="epithelial_cell")`

Compares a gene's regulatory wiring between population-averaged GREmLN ARACNe networks and TCGA tumor-state ARACNe networks. Classifies rewiring as low/moderate/high (Jaccard ≥ 0.6 / 0.3 thresholds), then validates conserved regulators via CASCADE experimental data (LINCS, DepMap, super-enhancers, DoRothEA). Output: tiered regulator list (conserved + CASCADE-validated, conserved without experimental support, tumor-acquired only).

Available TCGA cancer types: `brca`, `coad`, `hnsc`, `luad`, `lusc`, `ov`, `prad`, `ucec`. HNSC (head/neck squamous) is the closest available proxy for HPV-associated cervical squamous carcinoma.

## How It Works

Each tool makes a specific sequence of MCP calls to the child servers.

### `causal_chain_analysis` — TF path

```
1. CASCADE  ← get_gene_metadata(gene, cell_type)
              → gene_role: master_regulator | transcription_factor

2. [parallel]
   CASCADE  ← comprehensive_perturbation_analysis(gene, cell_type)
              → evidence_synthesis.key_findings (LINCS, DepMap, super-enhancer, DoRothEA)
   RegNetAgents ← comprehensive_gene_analysis(gene, cell_type)
              → network_rank, pathway_membership, domain_insights

3. Synthesis: for each candidate target — count independent evidence sources
   - network topology rank        [RegNetAgents]
   - pathway membership (Reactome)[RegNetAgents]
   - LINCS experimental knockdown [CASCADE]
   - DepMap CRISPR essentiality   [CASCADE]
   - super-enhancer / BET status  [CASCADE]
   - DoRothEA TF confidence tier  [CASCADE]
   - cBioPortal tumor expression  [CASCADE]
```

### `causal_chain_analysis` — Effector path

```
1. CASCADE  ← get_gene_metadata(gene, cell_type)
              → gene_role: effector | isolated

2. CASCADE  ← get_protein_interactions(gene)
              → PPI partners with STRING confidence scores

3. Find TF partner: classify top 10 PPI partners in parallel; pick highest-target TF

4. [parallel on TF partner]
   CASCADE  ← comprehensive_perturbation_analysis(tf_partner, cell_type)
   RegNetAgents ← comprehensive_gene_analysis(tf_partner, cell_type)

5. Return: gene → TF partner → downstream cascade → pathway annotation
   (e.g. APC → CTNNB1 → Wnt signaling)
```

This path is unreachable by either child server alone: CASCADE dead-ends on empty perturbation; RegNetAgents has no PPI data to bridge from effector to TF partner.

### `validate_therapeutic_targets`

```
1. [parallel]
   RegNetAgents ← comprehensive_gene_analysis → ranked upstream regulators (PageRank)
   CASCADE      ← therapeutic_target_discovery → super-enhancer + drug db candidates
   CASCADE      ← get_protein_interactions → top STRING PPI partners

2. Merge unique candidates from all three sources

3. [parallel, top 3 candidates]
   CASCADE ← comprehensive_perturbation_analysis(candidate, cell_type)

4. Return: 7-source corroboration table sorted by evidence count
```

### Synthesis Layer

The core of what Orchestra produces is a **corroboration count** per candidate — how many independent evidence sources agree. The LLM narrates this pre-scored table; it does not generate the scores.

| Evidence source | System | What it captures |
|---|---|---|
| PageRank rank | RegNetAgents | Topological centrality in ARACNe/GREmLN network |
| Pathway membership | RegNetAgents | Reactome enrichment membership |
| LINCS knockdown | CASCADE | Experimental CRISPR knockdown response |
| DepMap essentiality | CASCADE | Fitness dependency across cancer cell lines |
| Super-enhancer | CASCADE | BET inhibitor sensitivity via dbSUPER |
| DoRothEA tier | CASCADE | Curated TF-regulon confidence (A–E) |
| cBioPortal expression | CASCADE | Primary tumor expression from TCGA/ICGC |

**Important limitation:** RegNetAgents networks are inferred from heterogeneous cell states (GREmLN corpus, healthy + disease cells). CASCADE LINCS/DepMap data are from cancer cell lines. Agreement between systems is methodologically independent but does not imply a matched biological context — it is a hypothesis generator, not experimental validation.

## Why Orchestra

### APC — the strict-necessity case

APC is a scaffold protein with no transcriptional targets. A perturbation query returns empty results. Orchestra:
1. Detects APC as effector automatically
2. Queries STRING PPI → CTNNB1 (combined_score 0.999)
3. Simulates CTNNB1 overexpression → 2,740 affected genes
4. Enriches cascade against Reactome → Wnt signaling pathway

Neither child server completes this analysis alone.

### BRD4→MYC — the complementarity case

CASCADE's super-enhancer analysis identifies BRD4 as a therapeutic target for MYC (MYC has super-enhancers in 32 cell types → BET inhibitor sensitivity). BRD4 is absent from the ARACNe TF network — expected, because BRD4 acts through chromatin-level co-activation, not direct transcriptional regulation. Orchestra presents both views together: CASCADE epigenetic evidence + RegNetAgents network absence, which is itself informative. BRD4→MYC via super-enhancers is published (Lovén et al. 2013, *Cell*); BET inhibitors are in clinical trials.

## Supported Cell Types

Both child servers support these cell types (use exact strings as `cell_type` parameter):

```
cd4_t_cells          cd8_t_cells
cd14_monocytes       cd16_monocytes
nk_cells             nkt_cells
cd20_b_cells         monocyte-derived_dendritic_cells
erythrocytes         epithelial_cell
```

## Installation

### Required

```bash
git clone https://github.com/jab57/Orchestra.git
cd Orchestra
python -m venv env
env\Scripts\activate       # Windows
# source env/bin/activate  # macOS/Linux
pip install -r requirements.txt
```

### Child Server Prerequisites

Orchestra spawns RegNetAgents and CASCADE as subprocesses. Each must be installed with its own virtual environment:

```
c:\Dev\RegNetAgents\   (or your path)
  ├── env\Scripts\python.exe   ← Orchestra uses this Python
  └── regnetagents_langgraph_mcp_server.py

c:\Dev\CASCADE\
  ├── env\Scripts\python.exe   ← Orchestra uses this Python
  └── cascade_langgraph_mcp_server.py
```

Orchestra uses each child server's own `env/Scripts/python.exe`, not its own Python. This keeps NumPy and other dependencies isolated across projects.

If your child servers are in non-default locations, edit the factory functions in `mcp_client.py`:

```python
def make_cascade_client(cwd: str = r"c:\Dev\CASCADE") -> MCPClient:
def make_regnetagents_client(cwd: str = r"c:\Dev\RegNetAgents") -> MCPClient:
```

### Configure Environment

```bash
cp .env.example .env
```

Key variables:

| Variable | Default | Description |
|---|---|---|
| `REGNETAGENTS_SERVER_PATH` | — | Path to RegNetAgents server script (informational) |
| `CASCADE_SERVER_PATH` | — | Path to CASCADE server script (informational) |
| `USE_LLM_SYNTHESIS` | `false` | Enable LLM biological narrative in reports |
| `LLM_PROVIDER` | `ollama` | `ollama` \| `anthropic` |
| `ORCHESTRA_SSL_NO_VERIFY` | — | Set to `1` on networks with corporate SSL inspection |
| `REGNETAGENTS_TIMEOUT` | `20` | Per-tool timeout (seconds) |
| `CASCADE_TIMEOUT` | `30` | Per-tool timeout (seconds) |

### Verify Installation

```bash
pytest tests/
```

208 unit tests should pass. Integration tests (requiring live child servers) are skipped by default:

```bash
set ORCHESTRA_INTEGRATION_TESTS=1   # Windows
pytest tests/
```

Or run a live validation case:

```bash
python run_validation.py apc
```

### Claude Desktop Configuration

Add to `%APPDATA%\Claude\claude_desktop_config.json` (Windows) or `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS):

```json
{
  "mcpServers": {
    "orchestra": {
      "command": "C:\\Dev\\Orchestra\\env\\Scripts\\python.exe",
      "args": ["C:\\Dev\\Orchestra\\orchestra_mcp_server.py"],
      "env": {
        "PYTHONPATH": "C:\\Dev\\Orchestra"
      }
    }
  }
}
```

Restart Claude Desktop after editing. The six Orchestra tools will appear in the tools list.

## Usage

### Run the Validation Script

```bash
python run_validation.py            # all three biological validation cases
python run_validation.py apc        # APC effector analysis (strict-necessity case)
python run_validation.py tp53       # TP53 causal chain (TF path sanity check)
python run_validation.py brd4       # MYC therapeutic targets (BRD4 complementarity)
```

Results saved to `outputs/validation_<case>_<timestamp>.txt`.

### Use from Python

```python
import asyncio
from orchestra_langgraph_workflow import OrchestraWorkflow

async def main():
    workflow = OrchestraWorkflow()

    # TF path: parallel RegNetAgents + CASCADE, cross-system corroboration
    result = await workflow.run_analysis(
        gene="TP53",
        cell_type="epithelial_cell",
        analysis_type="causal_chain",
    )
    print(result["final_report"])

    # Effector path: APC → CTNNB1 via PPI
    result = await workflow.run_analysis(
        gene="APC",
        cell_type="epithelial_cell",
        analysis_type="causal_chain",
    )

    # Therapeutic target validation: 7-source corroboration table
    result = await workflow.run_analysis(
        gene="MYC",
        cell_type="cd4_t_cells",
        analysis_type="therapeutic_validation",
    )

asyncio.run(main())
```

See `examples/` for focused scripts: `apc_effector_analysis.py`, `tp53_causal_chain.py`, `brd4_target_validation.py`.

### Example Prompts (Claude Desktop)

**Causal chain analysis:**
- "Run Orchestra causal chain analysis on TP53 in epithelial cells"
- "What drives APC loss in colorectal cancer? Use effector analysis."
- "Analyze MYC regulatory network in CD4 T cells"

**Therapeutic target validation:**
- "Find therapeutic targets for MYC in CD4 T cells using Orchestra"
- "What are the druggable upstream regulators of CTNNB1 in epithelial cells?"

### LLM Synthesis Configuration

LLM synthesis is **off by default** (`USE_LLM_SYNTHESIS=false`). Orchestra returns structured text; Claude Desktop handles narrative interpretation. This is the recommended mode for MCP clients.

To enable a 2–3 sentence biological narrative prepended to each report:

**Option A: Ollama local (default, no API key)**

```bash
ollama pull llama3.1:8b
# Then set in .env:
USE_LLM_SYNTHESIS=true
LLM_PROVIDER=ollama
OLLAMA_MODEL=llama3.1:8b
```

**Option B: Anthropic**

```bash
USE_LLM_SYNTHESIS=true
LLM_PROVIDER=anthropic
LLM_API_KEY=sk-ant-...
LLM_MODEL=claude-haiku-4-5-20251001
```

| Variable | Default | Description |
|---|---|---|
| `USE_LLM_SYNTHESIS` | `false` | Enable LLM narrative |
| `LLM_PROVIDER` | `ollama` | `ollama` \| `anthropic` |
| `OLLAMA_MODEL` | `llama3.1:8b` | Ollama model name |
| `OLLAMA_HOST` | `http://localhost:11434` | Ollama server URL |
| `LLM_MODEL` | — | Override model (Anthropic only) |
| `LLM_API_KEY` | — | Required for Anthropic provider |

**Troubleshooting:**
- Ollama not running → `ollama serve` then retry
- Model not found → `ollama pull llama3.1:8b`
- SSL errors (corporate networks) → set `ORCHESTRA_SSL_NO_VERIFY=1` in `.env`

## Project Structure

```
Orchestra/
├── orchestra_mcp_server.py          # MCP server — exposes 6 composite tools to Claude Desktop
├── orchestra_langgraph_workflow.py  # LangGraph DAG, OrchestraState, all analysis paths
├── mcp_client.py                    # MCPClient class, subprocess lifecycle, factory functions
├── run_validation.py                # Standalone validation runner (3 biological cases)
├── generate_figure.py               # Architecture figure generator (JOSS paper)
├── figure_architecture.png/.pdf     # Architecture figure
├── paper.md                         # JOSS paper draft
├── paper.bib                        # JOSS bibliography
├── pyproject.toml                   # pytest configuration
├── requirements.txt                 # Python dependencies (exact versions)
├── .env.example                     # Environment variable template
├── docs/
│   ├── installation.md              # Step-by-step setup with troubleshooting
│   ├── usage.md                     # Tool reference and output interpretation
│   └── architecture.md             # Three-layer design and LangGraph DAG details
├── examples/
│   ├── apc_effector_analysis.py     # APC→CTNNB1 canonical effector use case
│   ├── tp53_causal_chain.py         # TP53 TF path demonstration
│   └── brd4_target_validation.py    # MYC/BRD4 therapeutic target validation
└── tests/
    ├── test_orchestra.py            # Routing, synthesis, report formatting (75+ unit tests)
    ├── test_mcp_client.py           # MCP client lifecycle and tool call tests (17 tests)
    ├── test_graceful_degradation.py # Degradation with mock clients (17 unit tests)
    ├── test_effector_analysis.py    # APC integration test (8 tests; requires child servers)
    ├── test_causal_chain.py         # TP53 integration test (9 tests; requires child servers)
    ├── test_gene_signature.py       # Gene signature path: routing, enrichment, synthesis (30 unit tests + 1 integration)
    └── test_network_comparison.py   # GREmLN vs TCGA network comparison (27 unit tests + 1 integration)
```

## Performance

### Per-Tool Latency Budget

| Tool call | Timeout | Notes |
|---|---|---|
| `comprehensive_perturbation_analysis` | 60s | CASCADE — most variable; depends on LINCS/STRING API |
| `comprehensive_gene_analysis` | 60s | RegNetAgents — network + pathway + domain agents |
| `get_protein_interactions` | 15s | CASCADE — STRING API |
| Other CASCADE calls | 30s | get_gene_metadata, therapeutic_target_discovery |

### End-to-End Workflow

| Workflow | Typical time | Notes |
|---|---|---|
| TF path (`causal_chain_analysis`) | ~30–60s | RegNetAgents + CASCADE run in parallel |
| Effector path (`causal_chain_analysis`) | ~30–60s | PPI lookup + parallel analysis on TF partner |
| `validate_therapeutic_targets` | ~60–120s | 3 parallel candidate extractions + 3 parallel validations |

### Manual vs. Orchestra

| Approach | Steps | Tool calls | Estimated time |
|---|---|---|---|
| Manual (APC use case) | 7 steps (query → dead-end → PPI → reason → simulate → enrich → interpret) | 6+ | ~20 min |
| Orchestra | 1 query | 1 | ~30s |

## Requirements

- Python 3.10+
- [RegNetAgents](https://github.com/jab57/RegNetAgents) — regulatory network analysis
- [CASCADE](https://github.com/jab57/CASCADE) — perturbation simulation and experimental corroboration
- `mcp==1.9.1`, `langgraph==0.2.34`, `python-dotenv==1.0.1`
- `ollama==0.6.1` — required only when `USE_LLM_SYNTHESIS=true` and `LLM_PROVIDER=ollama`
- `anthropic==0.97.0` — required only when `LLM_PROVIDER=anthropic`

## Running Tests

```bash
pip install pytest pytest-cov pytest-asyncio
pytest tests/ -v
```

Unit tests (208) run without live child servers. Integration tests (2) require RegNetAgents and CASCADE:

```bash
set ORCHESTRA_INTEGRATION_TESTS=1
pytest tests/ -v
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for coverage commands and testing guidelines.

## Citation

JOSS submission planned September 2026 after v1.0.0 release. Until then, please cite the GitHub repository.

```bibtex
@software{bird2026orchestra,
  author       = {Bird, Jose},
  title        = {Orchestra: MCP-Level Composition of Bioinformatics Servers for Multi-System Causal Reasoning},
  year         = {2026},
  url          = {https://github.com/jab57/Orchestra},
  note         = {JOSS submission planned September 2026}
}
```

## License

MIT License — see [LICENSE](LICENSE)
