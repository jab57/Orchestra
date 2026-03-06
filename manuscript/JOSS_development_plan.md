# Orchestra: JOSS Development Plan

**Created:** 2026-03-06
**Target journal:** Journal of Open Source Software (JOSS)
**Target submission:** September 2026
**Repo public:** March 2026 (starts 6-month JOSS clock)

---

## What Orchestra Is

Orchestra is a standalone MCP server that orchestrates two independent bioinformatics MCP servers — RegNetAgents and CASCADE — via the Model Context Protocol. A single natural language query triggers a coordinated LangGraph workflow that combines regulatory network analysis (RegNetAgents) and perturbation simulation (CASCADE) into an integrated causal report neither system can produce alone.

```
Claude Desktop
      |
      v
Orchestra (MCP server + LangGraph DAG)
      |                    |
      v                    v
RegNetAgents           CASCADE
(network analysis,     (perturbation sim,
 pathway enrichment,    PPI, LINCS, DepMap)
 domain agents)
```

---

## JOSS Contribution

**What this paper is about:**
- First published example of MCP server composition in bioinformatics
- A reusable orchestration pattern: one MCP server acting as client to multiple specialized MCP servers
- Automated causal reasoning across two independent systems via protocol-level composition
- Deterministic evidence layer (tool calls) + LLM interpretation layer (synthesis only), enforced by architecture

**What this paper is NOT about:**
- New biological algorithms or methods
- New databases or data sources
- Benchmarking RegNetAgents or CASCADE independently

**JOSS fit:** JOSS evaluates software quality, documentation, and utility — not biological novelty. Orchestra's contribution is architectural and infrastructural: a composable, protocol-level orchestration layer that makes multi-system causal analysis accessible from a single query.

---

## JOSS Requirements Checklist

| Requirement | Status | Notes |
|---|---|---|
| Open source license (MIT) | DONE | LICENSE in repo root |
| Public GitHub repository | DONE | March 2026 |
| 6 months public history | September 2026 | JOSS minimum |
| README with installation + usage | In progress | |
| API/function documentation | Pending | |
| Test suite | Pending | |
| `paper.md` (<1000 words) | Pending | August 2026 |
| `paper.bib` references | Pending | August 2026 |
| Statement of need | Pending | |
| Example usage | Pending | |
| Zenodo DOI | September 2026 | Tag v1.0.0 |

---

## Core Contribution: What Makes This Novel

**The problem Orchestra solves:**

Neither RegNetAgents nor CASCADE alone can answer: *"What transcription factors drive this gene signature, and what happens downstream if we inhibit the top candidate?"*

- RegNetAgents can rank upstream TFs by PageRank centrality
- CASCADE can simulate perturbation effects and validate via LINCS knockdown data
- Combining them manually requires the user to know which tools to call, in what order, and how to interpret cross-system results

**What Orchestra automates:**

1. Gene classification (TF vs. effector — determines analysis strategy)
2. Parallel evidence gathering across both systems
3. Cross-system synthesis (e.g., RegNetAgents top TF validated by CASCADE LINCS data)
4. LLM-generated causal interpretation over structured evidence

**The APC use case (canonical example):**

APC is a scaffold protein with no transcriptional targets. Manual analysis dead-ends at an empty perturbation result. Orchestra:
1. Classifies APC as effector (no targets)
2. Queries PPI (CASCADE) — finds CTNNB1 as key partner
3. Simulates CTNNB1 overexpression (CASCADE) — gets downstream cascade
4. Enriches cascade against Reactome pathways (RegNetAgents) — Wnt context
5. Returns integrated report explaining APC→CTNNB1→Wnt axis

---

## Implementation Plan

### Phase 1: MCP Client Foundation (April 2026)

**Goal:** Orchestra can call CASCADE and RegNetAgents tools via MCP protocol.

Tasks:
- [ ] Implement `MCPClient` class — spawns child server as subprocess, communicates via stdio MCP transport
- [ ] Connect to CASCADE: call `get_gene_metadata`, `comprehensive_perturbation_analysis`, `get_protein_interactions`
- [ ] Connect to RegNetAgents: call `comprehensive_gene_analysis`, `pathway_focused_analysis`
- [ ] Verify round-trip: Orchestra calls CASCADE tool, gets structured result back
- [ ] Subprocess lifecycle: startup, health check, graceful shutdown
- [ ] Basic error handling: child server timeout, tool call failure → store in state, continue

**Deliverable:** `mcp_client.py` — verified round-trip calls to both child servers

---

### Phase 2: Automated Effector Analysis (April 2026)

**Goal:** Prove the APC→CTNNB1 use case works end-to-end.

Tasks:
- [ ] Implement `_classify_gene` — calls CASCADE `get_gene_metadata` via MCP
- [ ] Implement `_routing_decision` — routes to `effector_path` for scaffold genes
- [ ] Implement `_run_effector_path`:
  - CASCADE: `get_protein_interactions(gene)` → find TF partners
  - CASCADE: `comprehensive_perturbation_analysis(tf_partner, overexpression)`
  - RegNetAgents: `pathway_focused_analysis(affected_genes)`
- [ ] Implement `_synthesize` — merge CASCADE + RegNetAgents outputs
- [ ] Implement `_generate_report` — Claude API call over structured evidence
- [ ] End-to-end test: `APC` in `epithelial_cell` → causal narrative

**Deliverable:** Working `automated_effector_analysis` composite tool

---

### Phase 3: Full Composite Tools (May 2026)

**Goal:** All three Orchestra tools implemented and tested.

#### `causal_chain_analysis` (TF path)
- RegNetAgents `comprehensive_gene_analysis` (network + pathway + domain agents)
- CASCADE `comprehensive_perturbation_analysis` (perturbation + LINCS + DepMap)
- Run in parallel via `asyncio.gather`
- Synthesis: which findings appear in both regulatory network AND experimental data?

Tasks:
- [ ] Implement `_run_tf_path` with parallel MCP calls
- [ ] Implement cross-system synthesis logic (overlap between RegNetAgents targets and CASCADE LINCS hits)
- [ ] Test: `TP53` in `epithelial_cell`
- [ ] Test: `MYC` in `cd4_t_cells`

#### `validate_therapeutic_targets`
- RegNetAgents: rank upstream regulators by PageRank
- CASCADE: validate top N regulators via perturbation + LINCS + DepMap essentiality
- Output: ranked targets with computational + experimental evidence

Tasks:
- [ ] Implement sequential → parallel validation logic
- [ ] Test: `BRD4→MYC` — RegNetAgents ranks BRD4, CASCADE confirms via LINCS
- [ ] Test: colorectal cancer TF panel

**Deliverable:** All three Orchestra tools (`causal_chain_analysis`, `validate_therapeutic_targets`, `effector_analysis`) working end-to-end

---

### Phase 4: LangGraph DAG + Synthesis (June 2026)

**Goal:** Full orchestrated workflow with proper state management and LLM synthesis.

Tasks:
- [ ] Wire all workflow nodes to real MCP calls (replace all TODOs)
- [ ] Implement `asyncio.gather` with `return_exceptions=True` in evidence nodes
- [ ] Implement synthesis node: cross-reference structured outputs from both systems
- [ ] Implement report node: Claude API call (`claude-sonnet-4-6`) with structured prompt
- [ ] Add `errors` dict population — failed MCP calls stored, workflow continues with partial data
- [ ] Add per-tool timeouts (30s perturbation, 15s PPI, 20s network analysis)
- [ ] Graceful degradation: if CASCADE unavailable, return RegNetAgents results + note; vice versa

**Deliverable:** Full LangGraph DAG running real analyses, producing causal reports

---

### Phase 5: Test Suite + Documentation (July 2026)

**Goal:** JOSS-quality software — tests, docs, installable by an independent researcher.

#### Test suite (`tests/`)
- [ ] `test_mcp_client.py` — mock MCP server, verify client round-trip
- [ ] `test_workflow_routing.py` — TF vs effector routing decision unit tests
- [ ] `test_effector_analysis.py` — APC integration test (requires CASCADE + RegNetAgents running)
- [ ] `test_causal_chain.py` — TP53 integration test
- [ ] `test_graceful_degradation.py` — behavior when one child server is down

#### Documentation
- [ ] `README.md` — value proposition, architecture diagram, installation, quick start
- [ ] `docs/installation.md` — step-by-step with troubleshooting
- [ ] `docs/usage.md` — example queries, expected outputs, interpreting reports
- [ ] `docs/architecture.md` — three-layer design, MCP client pattern, LangGraph DAG
- [ ] Inline docstrings on all public methods
- [ ] `.env.example` — all required environment variables documented

#### Examples
- [ ] `examples/apc_effector_analysis.py` — canonical APC→CTNNB1 use case
- [ ] `examples/tp53_causal_chain.py` — TF path demonstration
- [ ] `examples/brd4_target_validation.py` — therapeutic target validation

**Deliverable:** Installable, documented, tested package

---

### Phase 6: Biological Validation Cases (May–June 2026, parallel with Phases 3–4)

**Goal:** 2–3 biological use cases that demonstrate Orchestra's value over manual chaining.

| Use case | Gene | Cell type | What it demonstrates | Paper role |
|---|---|---|---|---|
| Effector routing | APC | epithelial_cell | Automatic scaffold protein handling — Orchestra doesn't dead-end on hard genes | Opens paper (origin story) |
| TF causal chain | TP53 | epithelial_cell | Parallel evidence gathering — sanity check on canonical gene | Middle (establishes baseline) |
| Cross-system therapeutic validation | BRD4→MYC | cd4_t_cells | Full Orchestra value proposition — RegNetAgents ranking + CASCADE LINCS + super-enhancers → actionable therapeutic insight | Closes paper (strongest case) |

**Why BRD4→MYC is the primary demonstration:**

This case requires both systems in a way that is immediately obvious to reviewers:
- RegNetAgents ranks BRD4 as a top upstream regulator of MYC by PageRank centrality (network topology)
- CASCADE confirms: BRD4 knockdown → MYC downregulation in LINCS experimental data (real CRISPR knockdown)
- CASCADE adds: MYC has super-enhancers → BET inhibitor sensitivity (druggability context)
- Synthesis layer connects all three: *target BRD4 with JQ1/OTX015 to suppress MYC in MYC-driven cancers*

Neither system alone produces this conclusion. RegNetAgents gives the ranking but no experimental validation. CASCADE gives the experimental data but no upstream regulator prioritization. The synthesis is Orchestra's contribution.

The BRD4→MYC relationship via super-enhancers is published (Lovén et al. 2013, *Cell*), making it independently verifiable. BET inhibitors (JQ1, OTX015) are in clinical trials — a reviewer who knows cancer biology will immediately recognize the result as meaningful and correct.

For each case, document:
- Manual chaining: tool calls required, steps, estimated time
- Orchestra: single query, time to report
- Cross-system agreement: which findings are corroborated by both systems
- Verification: how the result maps to published biology

---

### Phase 7: Paper + Submission (August–September 2026)

#### paper.md (JOSS format, <1000 words)
- [ ] **Summary** (~150 words): what Orchestra does, one-sentence architecture
- [ ] **Statement of need** (~200 words): why cross-system causal analysis requires orchestration; open with APC problem as the origin story — a scaffold protein that dead-ends manual analysis, motivating automated routing
- [ ] **Architecture** (~200 words): three-layer design, MCP composition pattern, LangGraph DAG
- [ ] **Example usage** (~150 words): APC→CTNNB1 code snippet (opens) + BRD4→MYC result summary (closes) — strongest case last
- [ ] **Performance** (~100 words): latency budget across 3 servers; comparison to manual chaining step count
- [ ] **References**: RegNetAgents (Bird et al. 2026a), CASCADE (Bird et al. 2026b), Lovén et al. 2013 (BRD4→MYC), MCP protocol, LangGraph

#### paper.bib
- [ ] RegNetAgents citation (own JOSS paper, pending)
- [ ] CASCADE citation (own JOSS paper, pending)
- [ ] Anthropic MCP protocol specification
- [ ] LangGraph / LangChain
- [ ] Relevant bioinformatics workflow papers

#### Submission prep
- [ ] Tag `v1.0.0` on GitHub
- [ ] Create Zenodo archive + DOI
- [ ] Verify all JOSS checklist items
- [ ] Submit via JOSS web form

---

## Timeline Summary

| Month | Phase | Key Deliverable |
|---|---|---|
| March 2026 | Repo public | 6-month clock starts; skeleton committed |
| April 2026 | Phase 1: MCP client | Round-trip calls to CASCADE + RegNetAgents verified |
| April 2026 | Phase 2: Effector proof of concept | APC end-to-end working |
| May 2026 | Phase 3: Full composite tools | All three Orchestra tools implemented |
| May–June 2026 | Phase 6: Biological validation | 3 use cases documented with metrics |
| June 2026 | Phase 4: LangGraph DAG | Full workflow, parallel evidence, LLM synthesis |
| July 2026 | Phase 5: Tests + docs | JOSS-quality software |
| August 2026 | Phase 7a: paper.md | Draft complete |
| September 2026 | Phase 7b: Submission | JOSS submission, v1.0.0 tagged, Zenodo DOI |

---

## Submission Strategy

**Recommended order:** RegNetAgents → CASCADE → Orchestra (staggered)

| Order | Paper | Target |
|---|---|---|
| 1 | RegNetAgents | Already in progress |
| 2 | CASCADE | After RegNetAgents accepted |
| 3 | Orchestra | September 2026 — cites both as published |

**Why order matters:**

Submitting all three simultaneously creates a circular dependency problem — each paper cites the others, and reviewers may ask why they are not one paper. By the time Orchestra submits in September 2026, RegNetAgents and CASCADE should be published JOSS papers it can cite directly. That transforms the statement of need from circular ("see our other papers") into concrete and verifiable ("as described in Bird et al. 2026a, 2026b").

---

## Anticipated Reviewer Concerns

### "Is this just glue code?"

The most likely reviewer pushback. The counter-argument must be explicit in `paper.md`:

- The **routing layer** (gene-role classification → analysis strategy selection) is non-trivial domain logic, not a passthrough
- The **synthesis layer** (cross-system evidence integration — identifying findings corroborated by both RegNetAgents network topology and CASCADE experimental data) is the core intellectual contribution
- The **graceful degradation** across subprocess boundaries (partial results if one child server fails) adds real engineering value
- The **MCP composition pattern** itself is a reusable contribution — documented and reproducible by others building multi-server bioinformatics pipelines

### "Why not one paper covering all three tools?"

Each tool solves a distinct problem and targets a distinct audience:
- RegNetAgents: regulatory network analysis (network biology community)
- CASCADE: perturbation simulation + experimental corroboration (systems biology / drug discovery)
- Orchestra: cross-system orchestration (bioinformatics infrastructure / MCP ecosystem)

Orchestra is only viable after the child servers exist. The three-paper structure reflects the actual development sequence, not redundancy.

---

## What This Paper Is NOT Competing With

Orchestra is not a general bioinformatics workflow tool. It does not compete with Nextflow, Snakemake, or Galaxy. It is specifically:
- An MCP composition layer for two specific bioinformatics servers
- A proof-of-concept for the MCP server composition pattern in scientific software
- A concrete solution to the multi-system causal analysis problem in gene regulatory biology

The JOSS paper argues that value, not biological novelty.

---

## Dependencies

| Dependency | Status | Required for |
|---|---|---|
| RegNetAgents v1.0.1 | DONE | Evidence layer (network analysis) |
| CASCADE v1.0.0 | Target March 26, 2026 | Evidence layer (perturbation + PPI) |
| `mcp` Python SDK | Available | MCP client + server implementation |
| `langgraph` | Available | Workflow orchestration |
| `anthropic` SDK | Available | LLM synthesis in report node |

---

## Document History

| Date | Change |
|---|---|
| 2026-03-06 | Initial plan created |
| 2026-03-06 | Added submission strategy (staggered order) and anticipated reviewer concerns |
| 2026-03-06 | Expanded Phase 6 with BRD4→MYC as primary demonstration case; added paper narrative structure |
