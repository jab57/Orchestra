# Orchestra

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

**MCP Orchestrator for Multi-System Causal Reasoning in Bioinformatics**

Orchestra automates multi-step causal reasoning across gene regulatory networks by composing two specialized MCP servers — [RegNetAgents](https://github.com/jab57/RegNetAgents) and CASCADE — via the Model Context Protocol. A single natural language query triggers an orchestrated workflow that combines regulatory network analysis, perturbation simulation, protein interactions, and pathway enrichment into a unified causal report.

Neither system alone can answer: *"What transcription factors drive this gene signature, and what happens downstream if we inhibit the top candidate?"* Orchestra can.

## Architecture

Orchestra acts as both an MCP server (to Claude Desktop or any MCP client) and an MCP client (to RegNetAgents and CASCADE). A LangGraph DAG coordinates three layers:

- **Decision layer** — classifies gene type, routes to appropriate composite analysis
- **Evidence layer** — parallel MCP calls to RegNetAgents and CASCADE
- **Explanation layer** — synthesis across both systems; optional LLM narrative

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
 domain agents)         super-enhancers)
```

## Tools

Orchestra exposes three composite tools — analytical capabilities that require both child servers and cannot be replicated by either alone:

**`causal_chain_analysis`**
Classifies a gene, runs regulatory network analysis (RegNetAgents) and perturbation simulation (CASCADE) in parallel, and synthesizes results into an integrated causal report. Identifies key upstream regulators and downstream perturbation effects in a single query.

**`validate_therapeutic_targets`**
Ranks upstream regulators by PageRank centrality (RegNetAgents), then validates top candidates against perturbation simulation and LINCS experimental knockdown data (CASCADE). Output: ranked targets with computational and experimental evidence.

**`effector_analysis`**
Handles scaffold/effector genes (e.g. APC) that have no direct transcriptional targets. Detects effector role, finds TF partners via protein-protein interactions (CASCADE), simulates the TF partner, and enriches the downstream cascade against Reactome pathways (RegNetAgents).

## Why Orchestra

Neither RegNetAgents nor CASCADE alone can answer: *"What transcription factors drive this gene signature, and what happens downstream if we inhibit the top candidate?"*

The canonical example is **BRD4→MYC**: RegNetAgents ranks BRD4 as the top upstream regulator of MYC by PageRank; CASCADE confirms via LINCS experimental knockdown data and identifies MYC super-enhancers indicating BET inhibitor sensitivity. Orchestra's synthesis layer connects these into a therapeutic recommendation that neither system produces independently.

A second canonical example is **APC mutation analysis**: APC is a scaffold protein with no transcriptional targets — a perturbation query dead-ends with empty results. Orchestra automatically detects this, queries protein interactions to find CTNNB1 as the key TF partner, simulates CTNNB1 overexpression, and enriches the downstream cascade against Reactome pathways, returning a complete APC→CTNNB1→Wnt causal explanation.

## Status

**Under active development.** Not yet ready for production use.

Target: v1.0.0 September 2026 | JOSS submission September 2026

## Installation

```bash
git clone https://github.com/jab57/Orchestra.git
cd Orchestra
python -m venv env
source env/bin/activate  # Windows: env\Scripts\activate
pip install -r requirements.txt
cp .env.example .env
# Edit .env with your server paths
```

Requires [RegNetAgents](https://github.com/jab57/RegNetAgents) and CASCADE to be installed separately.

## Usage

Coming soon.

## Citation

Coming soon (JOSS submission planned September 2026).

## License

MIT License — see [LICENSE](LICENSE)
