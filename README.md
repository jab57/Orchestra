# Orchestra

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/jab57/Orchestra/actions/workflows/test.yml/badge.svg)](https://github.com/jab57/Orchestra/actions/workflows/test.yml)
[![Draft JOSS Paper](https://github.com/jab57/Orchestra/actions/workflows/draft-pdf.yml/badge.svg)](https://github.com/jab57/Orchestra/actions/workflows/draft-pdf.yml)

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

## How It Works

Each tool makes a specific sequence of MCP calls to the child servers. The flows below show the actual calls and what data moves between systems.

### `causal_chain_analysis(gene, cell_type)`

```
1. CASCADE  ← get_gene_metadata(gene, cell_type)
              → gene_role: master_regulator | transcription_factor | effector | isolated

2. [parallel]
   CASCADE  ← comprehensive_perturbation_analysis(gene, cell_type)
              → evidence_synthesis.key_findings  (LINCS, DepMap, super-enhancer, DoRothEA)
   RegNetAgents ← comprehensive_gene_analysis(gene, cell_type)
              → network_rank, pathway_membership, domain_insights

3. Synthesis: for each candidate target — count how many independent sources agree
   - network topology rank        [RegNetAgents]
   - pathway membership (Reactome)[RegNetAgents]
   - LINCS experimental knockdown [CASCADE]
   - DepMap CRISPR essentiality   [CASCADE]
   - super-enhancer / BET status  [CASCADE]
   - DoRothEA TF confidence tier  [CASCADE]
   - cBioPortal tumor expression  [CASCADE]

4. Claude API ← structured evidence table → causal narrative
```

Targets corroborated by both systems score higher than those supported by one alone. The LLM narrates the evidence table; it does not generate the scores.

---

### `validate_therapeutic_targets(gene, cell_type)`

```
1. RegNetAgents ← query_network(gene_neighbors, gene, cell_type)
                  → ranked upstream regulators with PageRank scores

2. CASCADE      ← therapeutic_target_discovery(gene, cell_type)
                  → ranked upstream candidates with vulnerability scores

3. Cross-reference: regulators ranked highly by BOTH systems → priority candidates

4. [parallel, for top N candidates]
   CASCADE ← comprehensive_perturbation_analysis(candidate, cell_type)
             → key_findings per candidate (LINCS, DepMap, super-enhancer)

5. Synthesis: per-candidate evidence table (network rank + experimental evidence)

6. Return: ranked target list with corroboration counts and druggability notes
```

Example: CASCADE's super-enhancer analysis identifies BRD4 as a druggable target for MYC — MYC has super-enhancers across 32 cell types indicating BET inhibitor sensitivity. RegNetAgents simultaneously provides MYC's regulatory network context. BRD4 is absent from the ARACNe TF network (expected — it acts through chromatin-level co-activation, not direct mRNA regulation), which is itself informative. Neither system produces this combined epigenetic + network view alone.

---

### `effector_analysis(gene, cell_type)`

For scaffold/effector proteins with no transcriptional targets (e.g. APC, AXIN1):

```
1. CASCADE  ← get_gene_metadata(gene, cell_type)
              → gene_role: effector | isolated  (confirms why perturbation would be empty)

2. CASCADE  ← get_protein_interactions(gene)
              → PPI partners with confidence scores

3. Identify TF partners: filter PPI results for genes classified as TF in any cell type

4. [for top TF partner, e.g. CTNNB1]
   CASCADE  ← comprehensive_perturbation_analysis(tf_partner, overexpression)
              → downstream cascade (genes affected by TF overexpression)

5. RegNetAgents ← pathway_focused_analysis(cascade_gene_list, cell_type)
                  → Reactome pathway enrichment (which pathways are activated?)

6. Synthesis: gene → TF partner → downstream cascade → pathway annotation
   (e.g. APC → CTNNB1 overexpression → Wnt target genes → Wnt signaling pathway)

7. Return: causal explanation with PPI evidence + perturbation cascade + pathway context
```

This path is unreachable by either child server alone: CASCADE dead-ends on empty perturbation results; RegNetAgents has no PPI data to bridge from effector → TF partner.

---

### Synthesis Layer

The core of what Orchestra produces is a **corroboration count** per candidate: how many independent evidence sources — network topology (RegNetAgents) and experimental measurements (CASCADE) — agree on the same target. Agreement between methodologically independent systems is the signal. The LLM's role is interpretation of the pre-scored evidence table, not detection.

An important limitation: RegNetAgents networks are inferred from heterogeneous cell states (healthy + disease), while CASCADE's LINCS and DepMap data come from cancer cell lines. Agreement between systems is meaningful but does not imply a matched biological context — it is a hypothesis generator, not experimental validation.

## Why Orchestra

Neither RegNetAgents nor CASCADE alone can answer: *"What transcription factors drive this gene signature, and what happens downstream if we inhibit the top candidate?"*

The canonical example is **BRD4→MYC**: CASCADE's super-enhancer analysis identifies BRD4 as a therapeutic target for MYC (MYC has super-enhancers in 32 cell types → BET inhibitor sensitivity, confirmed by Lovén et al. 2013). RegNetAgents provides MYC's classical regulatory network context. BRD4 is absent from the ARACNe TF network — a biologically meaningful finding, because BRD4 acts through chromatin-level co-activation rather than direct transcriptional regulation. Orchestra presents both the epigenetic evidence and the network context together, producing a therapeutic recommendation neither system generates alone.

A second canonical example is **APC mutation analysis**: APC is a scaffold protein with no transcriptional targets — a perturbation query dead-ends with empty results. Orchestra automatically detects this, queries protein interactions to find CTNNB1 as the key TF partner, simulates CTNNB1 overexpression, and enriches the downstream cascade against Reactome pathways, returning a complete APC→CTNNB1→Wnt causal explanation.

## Status

**Under active development.**

| Issue | Description | Status |
|---|---|---|
| #1 MCP client | Round-trip calls to both child servers verified | ✓ Done |
| #2 Effector analysis | `effector_analysis` — APC→CTNNB1 end-to-end | ✓ Done |
| #3 Causal chain analysis | `causal_chain_analysis` — TF path (TP53, BRD4→MYC) | ✓ Done |
| #4 Therapeutic target validation | `validate_therapeutic_targets` — PageRank + CASCADE LINCS | ✓ Done |
| #5 Full synthesis layer | LangGraph DAG wired end-to-end with evidence scoring | ✓ Done |
| #6 Biological validation | APC, TP53, BRD4→MYC cases documented with metrics | ✓ Done |
| #7 Test suite | Full JOSS-quality coverage + CI | Partial (74 unit tests) |
| #8 Documentation + examples | Installable by an independent researcher | Pending |

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

## Testing

Orchestra ships a test suite covering routing logic, cross-system synthesis, and report formatting — all runnable without live child servers.

```bash
pytest tests/
```

Integration tests (requiring running RegNetAgents and CASCADE) are skipped by default and can be enabled with:

```bash
ORCHESTRA_INTEGRATION_TESTS=1 pytest tests/
```

CI runs automatically on every push and pull request via GitHub Actions.

## Usage

Coming soon.

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
