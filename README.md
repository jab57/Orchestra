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

## Status

**Under active development.** Not yet ready for production use.

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

Requires RegNetAgents and CASCADE to be installed separately.

## Usage

Coming soon.

## Citation

Coming soon (JOSS submission planned September 2026).

## License

MIT License — see [LICENSE](LICENSE)
