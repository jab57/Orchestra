# CASCADE Future Roadmap: Unified Bio-Orchestrator

> **STATUS UPDATE (2025-02)**: LangGraph orchestration has been implemented for CASCADE. Option 3 (MCP-over-MCP Orchestra) is now the preferred future direction for cross-project integration with RegNetAgents — see below.

## Architecture Options

| Aspect | Option 1: Two Separate Servers | Option 2: One Integrated Server | Option 3: MCP-over-MCP Orchestra |
|--------|--------------------------------|---------------------------------|----------------------------------|
| **MCP Servers** | 2 (RegNetAgents + CASCADE) | 1 (Bio-Orchestrator) | 3 (Orchestra + 2 child servers) |
| **Orchestration** | Claude conversation-level | LangGraph StateGraph | LangGraph via MCP protocol |
| **Tool Chaining** | Manual (user/Claude decides) | Automatic (gene-type routing) | Automatic (MCP-composed workflows) |
| **State Management** | None (stateless tools) | Unified state across analyses | Unified state, MCP transport |
| **Complexity** | Simple, modular | More complex, integrated | Moderate (MCP client plumbing) |
| **Maintenance** | Independent updates | Coordinated updates | Independent (protocol boundary) |
| **User Experience** | Multiple tool calls needed | Single query → full analysis | Single query → full analysis |
| **Integration** | N/A | Direct Python imports | MCP protocol (language-agnostic) |
| **Novelty** | Standard | Standard | First MCP composition in bioinformatics |
| **Status** | **In Production** | Postponed | **Preferred future direction** |

---

## Option 1: Two Separate MCP Servers (Current)

```
┌─────────────────────────────────────────────────────────────┐
│                    Claude Desktop / MCP Client              │
│              (orchestrates via conversation)                │
└──────────────────────┬──────────────────────────────────────┘
                       │
           ┌───────────┴───────────┐
           ▼                       ▼
┌─────────────────────┐   ┌─────────────────────┐
│   RegNetAgents      │   │      CASCADE        │
│    MCP Server       │   │    MCP Server       │
│                     │   │                     │
│ - Network analysis  │   │ - Perturbation sim  │
│ - Pathway enrich    │   │ - Gene embeddings   │
│ - Domain insights   │   │ - PPI (STRING)      │
│ - LangGraph workflow│   │ - Smart suggestions │
└─────────────────────┘   └─────────────────────┘

Location: c:\Dev\RegNetAgents   Location: c:\Dev\CASCADE
Status: Frozen (no changes)     Status: Active development
```

### Current Capabilities

**RegNetAgents MCP Server:**
- `comprehensive_gene_analysis` - Full LangGraph workflow
- `multi_gene_analysis` - Parallel gene processing
- `pathway_focused_analysis` - Reactome integration
- Domain-specific agents (cancer, drug, clinical, systems biology)

**CASCADE MCP Server (25 tools):**
- `comprehensive_perturbation_analysis` - Full automated workflow with intelligent routing
- `quick_perturbation` - Fast knockdown/overexpression without full workflow
- `multi_gene_analysis` - Parallel multi-gene processing
- `cross_cell_comparison` - Compare gene across all cell types
- `therapeutic_target_discovery` - Find druggable regulators and partners
- `get_gene_metadata` - Gene type classification
- `get_protein_interactions` - STRING database
- `find_similar_genes` - Embedding-based similarity
- `analyze_network_vulnerability` - Drug target discovery
- LINCS L1000 tools - Experimental knockdown corroboration
- Super-enhancer tools - BRD4/BET inhibitor sensitivity (dbSUPER)
- DoRothEA tools - TF regulon validation (via decoupler)
- Smart suggestions when results are empty

### Current Limitation

When analyzing genes like APC (scaffold proteins with no transcriptional targets), Claude must manually recognize the biological context and chain tools appropriately:

```
APC knockdown → empty results
    ↓ (Claude recognizes APC is effector)
APC PPI → shows CTNNB1 as key partner
    ↓ (Claude reasons about APC→CTNNB1 relationship)
CTNNB1 overexpression → rich cascade data
```

The new `gene_metadata` and `suggestions` fields help Claude make these decisions, but the reasoning still happens at the conversation level.

### Partial Solution: Claude Code Skill (Implemented 2025-01)

The CASCADE skill (`.claude/skills/cascade/SKILL.md`) provides workflow guidance to Claude Code, teaching it:
- When to use knockdown vs PPI analysis
- How to handle effector genes (scaffold proteins)
- Suggested tool chaining (e.g., empty knockdown → check metadata → use PPI)
- Trigger keywords to distinguish from other MCP servers (RegNetAgents)

This reduces manual intervention but doesn't fully automate the orchestration. Option 2 remains the long-term vision for seamless single-tool analysis.

---

## Option 2: One Integrated Server (Future)

```
┌─────────────────────────────────────────────────────────────┐
│                    Claude Desktop / MCP Client              │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│              Unified Bio-Orchestrator MCP Server            │
│                    (NEW PROJECT)                            │
│                                                             │
│  ┌─────────────────────────────────────────────────────┐   │
│  │              LangGraph StateGraph                    │   │
│  │                                                      │   │
│  │  analyze_mutation_implications()                     │   │
│  │    → assess_gene_type                               │   │
│  │    → route_to_appropriate_analysis                  │   │
│  │    → synthesize_results                             │   │
│  └─────────────────────────────────────────────────────┘   │
│                          │                                  │
│            ┌─────────────┴─────────────┐                   │
│            ▼                           ▼                    │
│   ┌─────────────────┐         ┌─────────────────┐          │
│   │  RegNetAgents   │         │    CASCADE       │          │
│   │  (as library)   │         │  (as library)   │          │
│   └─────────────────┘         └─────────────────┘          │
└─────────────────────────────────────────────────────────────┘

Total MCP Servers exposed to Claude: 1
```

### Key Benefits

1. **Intelligent Routing**: Agent automatically decides whether to use network analysis, PPI, or embeddings based on gene type
2. **Unified State**: Single state object tracks all findings across both systems
3. **Biological Reasoning**: LLM-powered decisions about analysis strategy
4. **Seamless UX**: User asks "Analyze APC mutation" → gets complete cascade analysis

### Proposed State Schema

```python
class UnifiedBioState(TypedDict):
    # Input
    gene: str
    cell_type: str
    analysis_type: str  # "mutation", "overexpression", "pathway"

    # Gene Classification
    gene_metadata: dict  # From CASCADE get_gene_metadata
    is_transcription_factor: bool
    known_complexes: list[str]

    # Analysis Results
    network_context: dict      # From RegNetAgents
    perturbation_results: dict # From CASCADE
    ppi_data: dict             # From CASCADE STRING
    pathway_enrichment: dict   # From RegNetAgents
    embedding_similar: list    # From CASCADE (GREmLN embeddings)

    # Routing State
    analyses_completed: list[str]
    next_analysis: str
    suggestions: list[dict]

    # Output
    final_report: str
    confidence_score: float
```

### Proposed Workflow Graph

```
                    ┌─────────────┐
                    │   START     │
                    └──────┬──────┘
                           ▼
                ┌──────────────────────┐
                │  assess_gene_type    │
                │  (calls CASCADE      │
                │   get_gene_metadata) │
                └──────────┬───────────┘
                           ▼
              ┌────────────────────────────┐
              │     route_analysis         │
              │  (conditional edge)        │
              └─────┬──────────┬───────────┘
                    │          │
        ┌───────────┘          └───────────┐
        ▼                                  ▼
┌───────────────┐                 ┌───────────────┐
│ TF Analysis   │                 │ Effector      │
│ Path          │                 │ Analysis Path │
├───────────────┤                 ├───────────────┤
│ 1. knockdown  │                 │ 1. get_ppi    │
│ 2. pathways   │                 │ 2. find TF    │
│ 3. domain     │                 │    partners   │
│    insights   │                 │ 3. analyze TF │
└───────┬───────┘                 │    partners   │
        │                         └───────┬───────┘
        │                                 │
        └─────────────┬───────────────────┘
                      ▼
              ┌───────────────┐
              │  synthesize   │
              │  results      │
              └───────┬───────┘
                      ▼
              ┌───────────────┐
              │     END       │
              └───────────────┘
```

---

## Option 3: Standalone Orchestra Server (MCP-over-MCP) — Preferred

Rather than merging CASCADE and RegNetAgents into a single codebase (Option 2), this approach creates a **standalone MCP server that orchestrates both systems via the MCP protocol itself**. The Orchestra server is an MCP client to CASCADE and RegNetAgents, while being an MCP server to Claude Desktop.

No one in bioinformatics has published MCP server composition. This makes Option 3 both technically cleaner and a stronger contribution for a JOSS paper.

```
┌─────────────────────────────────────────────────────────────┐
│                    Claude Desktop / MCP Client               │
└──────────────────────────┬──────────────────────────────────┘
                           │ MCP protocol
                           ▼
┌─────────────────────────────────────────────────────────────┐
│              Orchestra Server                            │
│              (NEW — c:\Dev\Orchestra)                    │
│                                                              │
│  ┌────────────────────────────────────────────────────────┐  │
│  │            LangGraph StateGraph Workflow                │  │
│  │                                                         │  │
│  │  initialize → classify_gene → route_analysis            │  │
│  │                    ↓                                    │  │
│  │         ┌──────────┼──────────┐                        │  │
│  │         ↓          ↓          ↓                        │  │
│  │  cascade_tasks  regnet_tasks  combined_tasks           │  │
│  │   (parallel)    (parallel)     (parallel)              │  │
│  │         └──────────┼──────────┘                        │  │
│  │                    ↓                                    │  │
│  │            synthesize → generate_report → END          │  │
│  └────────────────────────────────────────────────────────┘  │
│                    │                     │                    │
│         MCP calls ↓                     ↓ MCP calls          │
│  ┌─────────────────────┐     ┌─────────────────────┐        │
│  │   CASCADE            │     │   RegNetAgents       │        │
│  │   MCP Server         │     │   MCP Server         │        │
│  │   (subprocess)       │     │   (subprocess)       │        │
│  │                      │     │                      │        │
│  │ - Perturbation sim   │     │ - Network analysis   │        │
│  │ - Gene embeddings    │     │ - Pathway enrichment │        │
│  │ - PPI (STRING)       │     │ - Domain agents      │        │
│  │ - LINCS data         │     │ - Causal reasoning   │        │
│  │ - Super-enhancers    │     │ - PageRank centrality│        │
│  └─────────────────────┘     └─────────────────────┘        │
└─────────────────────────────────────────────────────────────┘

Total MCP Servers exposed to Claude: 1 (Orchestra)
Child servers managed internally as subprocesses
```

### Why MCP-over-MCP Instead of Direct Imports (Option 2)

| Consideration | Option 2 (Direct Imports) | Option 3 (MCP-over-MCP) |
|---------------|---------------------------|--------------------------|
| **Integration** | Python imports, tight coupling | MCP protocol, clean boundary |
| **RegNetAgents changes** | May need refactoring for import | Zero changes (already an MCP server) |
| **Reusability** | Only usable as Python library | Any MCP client can connect |
| **Testing** | Mock Python objects | Mock MCP responses (protocol-level) |
| **Language lock-in** | Python only | Protocol-agnostic |
| **Novelty for JOSS** | Standard integration | First MCP composition in bioinformatics |
| **Failure isolation** | Shared process, shared crashes | Subprocess isolation |
| **Development overhead** | Import wrangling, dependency conflicts | ~200-300 lines MCP client plumbing |

### Composite Tools (Neither Server Alone)

These tools create genuinely new analytical capabilities by combining both systems:

#### 1. `validate_therapeutic_targets`
RegNetAgents ranks targets via PageRank centrality → CASCADE validates top candidates via perturbation simulation + LINCS experimental knockdown data.
```
RegNetAgents: PageRank → [MYC, TP53, STAT3, ...]
CASCADE per target: perturbation_sim + LINCS + super_enhancer check
Output: ranked targets with experimental + computational validation
```

#### 2. `pathway_perturbation_analysis`
CASCADE simulates gene perturbation → RegNetAgents enriches the affected gene set against Reactome pathways.
```
CASCADE: knockdown(STAT3) → [affected genes]
RegNetAgents: pathway_enrichment([affected genes]) → Reactome pathways
Output: which biological pathways are disrupted by this perturbation
```

#### 3. `bidirectional_causal_analysis`
RegNetAgents identifies upstream regulators of a gene → CASCADE simulates knockdown of top regulators to predict downstream effects.
```
RegNetAgents: upstream_regulators(MYC) → [BRD4, MAX, ...]
CASCADE per regulator: knockdown_simulation → downstream effects
Output: full causal chain (upstream regulators → gene → downstream cascade)
```

#### 4. `automated_effector_analysis`
Solves the APC→CTNNB1 problem (see [Reference: The Original Problem](#reference-the-original-problem)): detect effector gene → PPI → find TF partner → simulate the TF.
```
CASCADE: classify(APC) → effector (no transcriptional targets)
CASCADE: PPI(APC) → CTNNB1 (TF partner)
CASCADE: overexpression(CTNNB1) → downstream cascade
RegNetAgents: pathway_enrichment(cascade) → Wnt pathway context
Output: complete analysis despite APC having no direct targets
```

#### 5. `cross_celltype_causal_map`
Both systems' cross-cell analyses combined: CASCADE's cross-cell perturbation effects + RegNetAgents' cross-cell network topology.
```
CASCADE: perturbation across [cd4, cd8, nk, monocytes, ...]
RegNetAgents: network centrality across same cell types
Output: cell-type-specific vulnerability map with perturbation + topology
```

### Proposed State Schema

```python
class OrchestraState(TypedDict):
    # Input
    gene: str
    cell_type: str
    analysis_type: str  # "therapeutic_validation", "pathway_perturbation", etc.
    analysis_depth: str  # "basic", "comprehensive", "focused"

    # Gene Classification (from CASCADE)
    ensembl_id: str
    gene_symbol: str
    gene_role: str  # master_regulator, transcription_factor, effector, isolated

    # CASCADE Results (via MCP)
    perturbation_result: dict
    ppi_interactions: dict
    similar_genes: list
    lincs_effects: dict
    super_enhancer_status: dict
    vulnerability_analysis: dict

    # RegNetAgents Results (via MCP)
    network_analysis: dict       # comprehensive_gene_analysis output
    pathway_enrichment: dict     # pathway_focused_analysis output
    centrality_scores: dict      # PageRank, betweenness, etc.
    domain_insights: dict        # cancer/drug/clinical agent outputs

    # Composite Results (neither server alone)
    validated_targets: list[dict]
    causal_chain: dict
    pathway_disruptions: list[dict]

    # Routing State
    completed_steps: list[str]
    pending_steps: list[str]

    # Output
    comprehensive_report: str
    therapeutic_suggestions: list[dict]
    confidence_score: float
```

### LangGraph Workflow DAG

```
                    ┌─────────────┐
                    │  initialize  │
                    └──────┬──────┘
                           ▼
                ┌──────────────────────┐
                │   classify_gene      │
                │   (CASCADE:          │
                │    get_gene_metadata)│
                └──────────┬───────────┘
                           ▼
              ┌────────────────────────────┐
              │     route_analysis         │
              │  (gene_role + depth →      │
              │   select composite tools)  │
              └─────┬──────────┬───────────┘
                    │          │
        ┌───────────┘          └───────────┐
        ▼                                  ▼
┌───────────────────┐            ┌───────────────────┐
│ TF / Master Reg   │            │ Effector Path     │
│ Path              │            │                   │
├───────────────────┤            ├───────────────────┤
│ CASCADE:          │            │ CASCADE:          │
│  perturbation_sim │            │  get_ppi          │
│ RegNetAgents:     │            │  find TF partners │
│  pathway_enrich   │            │  simulate TF      │
│  domain_agents    │            │ RegNetAgents:     │
│ Combined:         │            │  pathway_enrich   │
│  validate_targets │            │  domain_agents    │
└───────┬───────────┘            └───────┬───────────┘
        │                                │
        └─────────────┬──────────────────┘
                      ▼
              ┌───────────────────┐
              │   synthesize      │
              │   (merge results, │
              │    score targets) │
              └───────┬───────────┘
                      ▼
              ┌───────────────────┐
              │  generate_report  │
              └───────┬───────────┘
                      ▼
              ┌───────────────────┐
              │       END         │
              └───────────────────┘
```

### Engineering Considerations

- **MCP client plumbing**: ~200-300 lines to spawn child servers as subprocesses and communicate via stdio MCP transport. The `mcp` Python SDK provides `ClientSession` and `StdioServerParameters` for this.
- **Subprocess lifecycle**: Orchestra starts CASCADE and RegNetAgents as child processes on startup, monitors health, restarts on failure.
- **Error propagation**: MCP errors from child servers are caught and surfaced in the Orchestra state (e.g., `cascade_error`, `regnet_error`) rather than crashing the workflow.
- **Timeout handling**: Per-tool timeouts (e.g., 30s for perturbation, 15s for PPI) with graceful degradation — if one child server is slow, the other's results are still returned.
- **Proposed location**: `c:\Dev\Orchestra` (standalone project, no modifications to CASCADE or RegNetAgents).

### Implementation Phases

#### Phase 1: MCP Client Foundation
- [ ] Create `c:\Dev\Orchestra` project structure
- [ ] Implement MCP client that spawns CASCADE as subprocess
- [ ] Implement MCP client that spawns RegNetAgents as subprocess
- [ ] Verify round-trip: Orchestra calls CASCADE tool → gets result
- [ ] Basic health check and subprocess restart logic

#### Phase 2: Composite Tool — `automated_effector_analysis`
- [ ] Implement the APC→CTNNB1 workflow as proof of concept
- [ ] classify gene → PPI → find TF partner → simulate → pathway enrich
- [ ] Validate against the original problem scenario

#### Phase 3: Full Composite Tools
- [ ] `validate_therapeutic_targets` (PageRank → perturbation)
- [ ] `pathway_perturbation_analysis` (simulate → enrich)
- [ ] `bidirectional_causal_analysis` (upstream → simulate)
- [ ] `cross_celltype_causal_map` (both systems' cross-cell)

#### Phase 4: LangGraph Orchestration
- [ ] Implement `OrchestraState` schema
- [ ] Build LangGraph DAG with routing logic
- [ ] Parallel batch execution for independent analyses
- [ ] Synthesis and report generation

#### Phase 5: MCP Server Wrapper & Testing
- [ ] Expose Orchestra as MCP server (server to Claude, client to children)
- [ ] APC mutation end-to-end test
- [ ] TP53, MYC, BRCA1 test cases
- [ ] Performance benchmarks (latency budget across 3 servers)

---

## Prerequisites Before Starting

### Technical Requirements

1. **Python Environment Compatibility**
   - Both projects use Python 3.10+
   - Verify no dependency conflicts between RegNetAgents and CASCADE

2. **Import Strategy**
   - RegNetAgents: Import `RegNetAgentsWorkflow` class directly
   - CASCADE: Import tool functions from `tools/` modules
   - Option 2: Direct Python imports; Option 3: MCP protocol calls

3. **State Management**
   - Design unified state schema that captures both systems' outputs
   - Handle async execution (both systems use asyncio)

### Architectural Decisions Needed (for Option 2 / Option 3)

1. **Where does orchestrator live?**
   - New standalone project (`c:\Dev\Orchestra` for Option 3) - *Recommended*
   - Inside CASCADE as optional module
   - Inside RegNetAgents as extension

2. **How to handle RegNetAgents being frozen?**
   - Import as-is without modification
   - Create wrapper classes if needed
   - Fork if changes become necessary

3. **Deployment model**
   - Single MCP server (orchestrator only)
   - Keep individual servers available for direct access?

### Testing Strategy

1. **Unit tests**: Each routing decision in isolation
2. **Integration tests**: Full workflows (APC → CTNNB1 cascade)
3. **Comparison tests**: Orchestrator results vs manual chaining

---

## Implementation Phases (When Ready)

### Phase 1: Foundation
- [ ] Create new project structure
- [ ] Set up shared state schema
- [ ] Import RegNetAgents and CASCADE as libraries
- [ ] Basic LangGraph workflow skeleton

### Phase 2: Gene Type Routing
- [ ] Implement `assess_gene_type` node
- [ ] Create conditional routing logic
- [ ] TF path: knockdown → pathways
- [ ] Effector path: PPI → partner analysis

### Phase 3: Intelligent Synthesis
- [ ] Result aggregation across systems
- [ ] LLM-powered biological interpretation
- [ ] Confidence scoring

### Phase 4: MCP Server Wrapper
- [ ] Expose orchestrator via MCP
- [ ] Single tool: `analyze_mutation_implications`
- [ ] Streaming progress updates

### Phase 5: Testing & Validation
- [ ] APC mutation test case (the original problem)
- [ ] TP53, MYC, BRCA1 test cases
- [ ] Performance benchmarks

---

## Reference: The Original Problem

This roadmap exists because of the APC mutation analysis scenario:

**What happened:**
1. User ran `quick_perturbation("APC")` → empty results
2. Tool returned no suggestions about why
3. User had to manually recognize APC is a scaffold protein
4. User manually switched to `get_protein_interactions("APC")`
5. User identified CTNNB1 as key partner
6. User reasoned that APC inhibits CTNNB1
7. User ran `quick_perturbation("CTNNB1", perturbation_type="overexpression")` → rich cascade

**What we want:**
1. User runs `analyze_mutation_implications("APC")`
2. Orchestrator detects APC is an effector (no transcriptional targets)
3. Orchestrator automatically queries PPI
4. Orchestrator identifies CTNNB1 as TF partner
5. Orchestrator reasons about APC→CTNNB1 relationship
6. Orchestrator runs CTNNB1 overexpression analysis
7. User receives integrated report with full cascade

---

## Near-Term Enhancements

| Feature | Status | Date | Description |
|---------|--------|------|-------------|
| MCP Resources for metadata discovery | **COMPLETE Issue #5** | 2026-02-24 | 5 browsable resource endpoints (cascade://cell-types, lincs/summary, model/status, network/{cell_type}/summary, gene/{symbol}/{cell_type}); fixes startup timeout via background DoRothEA pre-warm |
| True parallel batch execution | **COMPLETE Issue #4** | 2026-02-17 | asyncio.gather() for 3 batch groups, asyncio.to_thread() for 9 blocking methods, thread-safe model singleton, network/STRING caching; ~25-35s comprehensive analysis |
| DoRothEA TF Regulon Validation | **COMPLETE** | 2026-02 | Curated TF regulons via decoupler (A-E confidence) |
| LLM-Powered Insights | **COMPLETE** | 2025-02-02 | Ollama integration for biological interpretation |
| LangGraph Orchestration | **COMPLETE** | 2025-02-02 | Intelligent workflow routing with parallel execution |
| LINCS L1000 (Harmonizome) | **COMPLETE** | 2025-01-30 | Expression perturbation from CRISPR knockdowns |
| Super-Enhancer Annotations | **COMPLETE** | 2025-01-30 | BRD4 druggability from dbSUPER |
| DepMap CRISPR Essentiality | **SCHEDULED Issue #7** | 2026-03 | Empirical essentiality validation (Chronos scores + CNV/mutation stratification); 8th evidence source; phenotypic validation layer |
| cBioPortal Primary Tumors | **SCHEDULED Issue #8** | 2026-03 | TCGA primary tumor expression + mutation frequency via REST API; closes cell-line vs. primary-tissue gap |
| Raw LINCS Integration | **SCHEDULED Issue #16** | 2026-04 | Full LINCS data replacing Harmonizome; fixes BRD4→MYC gap |

---

### Completed: LLM-Powered Biological Insights (2025-02-02)

Added optional Ollama integration for AI-generated interpretation of perturbation analysis results:

**New Feature:**
- `include_llm_insights=True` parameter in `comprehensive_perturbation_analysis`

**Key Capabilities:**
- **Mechanism Summary**: 2-3 sentence explanation of what the perturbation does mechanistically
- **Therapeutic Implications**: Drug development relevance assessment
- **Pathway Identification**: Automated identification of key affected pathways
- **Confidence Assessment**: High/medium/low confidence with rationale
- **Follow-up Suggestions**: Intelligent recommendations for further analysis
- **Biological Narrative**: 3-4 sentence synthesis suitable for research reports

**Configuration:**
- Auto-detects local vs cloud Ollama (set `OLLAMA_API_KEY` for cloud)
- Configurable model, temperature, and timeout via environment variables
- Graceful fallback: Returns structured data if LLM unavailable
- Default OFF to avoid latency for quick queries

**Files:**
- `cascade_langgraph_workflow.py` - Added `_synthesize_insights` node, `_call_llm_synthesis` method
- `cascade_langgraph_mcp_server.py` - Added `include_llm_insights` parameter to tool schema
- `.env.example` - LLM configuration template

---

### Completed: LangGraph Orchestration (2025-02-02)

Implemented LangGraph-based workflow orchestration within CASCADE MCP server:

**New Files:**
- `cascade_langgraph_workflow.py` - Core LangGraph StateGraph workflow
- `cascade_langgraph_mcp_server.py` - MCP server exposing 25 tools

**Key Features:**
- **Intelligent Routing**: Automatically selects analysis strategy based on gene type (master_regulator, transcription_factor, effector, isolated)
- **Parallel Batch Execution**: Independent analyses run concurrently (batch_core_analysis, batch_external_data, batch_insights)
- **Automatic Synthesis**: Generates comprehensive reports with actionable recommendations
- **Graceful Degradation**: Falls back to network-only if embeddings unavailable

**New Workflow Tools:**
- `comprehensive_perturbation_analysis` - Full automated analysis with intelligent routing
- `quick_perturbation` - Fast knockdown/overexpression without full workflow
- `multi_gene_analysis` - Analyze multiple genes in parallel
- `cross_cell_comparison` - Compare gene across all cell types
- `therapeutic_target_discovery` - Find druggable regulators and partners

**Performance:**
- Basic analysis: ~3s
- Comprehensive analysis: ~8-10s
- Multi-gene parallel (3 genes): ~10s

**Note:** This addresses the "CASCADE-only" orchestration need. The Unified Bio-Orchestrator (Option 2 - cross-project integration with RegNetAgents) remains postponed.

---

### Completed: LINCS L1000 Expression Perturbation (2025-01-30)

Added tools to find regulatory relationships from experimental CRISPR knockdown data:
- `find_expression_regulators(gene)` - what knockdowns affect this gene?
- `get_knockdown_effects(gene)` - what does this knockdown affect?

**Data source**: Harmonizome LINCS L1000 CRISPR Knockout Consensus Signatures

**Limitation**: Harmonizome pre-filters data, removing some validated relationships (e.g., BRD4 → MYC).

### Scheduled: Raw LINCS Integration (Issue #16, April 2026)

**Status**: SCHEDULED — Issue #16, target Apr 20 open / Apr 25 fix, included in v1.0.2 release

**Problem**: BRD4 → MYC (validated BET inhibitor target) not in Harmonizome data. The paper claims LINCS as a corroboration layer; a reviewer familiar with LINCS could flag this gap.

**Solution**: Integrate raw LINCS L1000 data from clue.io (GEO: GSE92742)

**Effort**: Medium (larger dataset, GCTX format parsing)

**Validation**: `find_expression_regulators("MYC")` should return BRD4

**Tasks**:
- [ ] Download raw LINCS L1000 data (level 5 - moderated z-scores)
- [ ] Implement GCTX file parser (HDF5 format via h5py)
- [ ] Create gene-knockdown → affected-genes mapping with directional z-score thresholding
- [ ] Add `find_expression_regulators_raw(gene)` tool (or augment existing tool)
- [ ] Validate: BRD4 knockdown shows MYC downregulation

### Completed: Super-Enhancer Annotations (2025-01-30)

Added tools to identify BRD4/BET inhibitor sensitive genes:
- `check_super_enhancer(gene)` - Check if gene has super-enhancers
- `check_genes_super_enhancers(genes)` - Screen multiple genes

**Data source**: dbSUPER (69K associations, 10K genes, 102 cell types)

**Validation**: MYC has SE in 32 cell types (BRD4-sensitive), TP53 has none

**Therapeutic value**: Enables drug discovery for "undruggable" targets like MYC - if a gene has super-enhancers, BRD4 inhibitors (JQ1, OTX015) may reduce its expression

---

## Document History

| Date | Change |
|------|--------|
| 2026-02-24 | Added Issue #5 complete (MCP Resources) and Issue #4 complete (parallel batch execution) to Near-Term Enhancements table |
| 2026-02-22 | Scheduled Raw LINCS as Issue #15 (April 2026, pre-JOSS); removed Expression Data Fetching (stub deleted Feb 18) |
| 2026-02-24 | Added DepMap CRISPR essentiality as Issue #7 (March 2026); renumbered Issues #7-14 to #8-15 |
| 2026-02-24 | Added cBioPortal primary tumor data as Issue #8; folded CNV/mutation stratification into Issue #7; renumbered Issues #8-15 to #9-16; added static network limitation to paper.md |
| 2026-02-12 | Updated tool names to match current 25-tool MCP server; added DoRothEA to completed features |
| 2025-02-11 | Added Option 3: MCP-over-MCP Orchestra as preferred future direction |
| 2025-02-06 | Renamed project from GREmLN to CASCADE |
| 2025-02-02 | Added LLM-powered biological insights (Ollama integration) |
| 2025-02-02 | Added LangGraph orchestration (COMPLETE), updated architecture diagrams |
| 2025-02-01 | Added status summary table, Expression Data Fetching section, task checklists |
| 2025-01-30 | Added super-enhancer annotations (dbSUPER) for BRD4 druggability |
| 2025-01-30 | Added LINCS L1000 tools, documented limitations, planned raw LINCS integration |
| 2025-01-25 | Initial roadmap created (Jose A. Bird, PhD) |

**Current Status**: Active development - Issues #1–5 complete (DoRothEA, installation verification, test coverage, parallel batch execution, MCP Resources). DepMap CRISPR essentiality (Issue #7) and cBioPortal primary tumor data (Issue #8) scheduled for March 2026. Raw LINCS scheduled as Issue #16 (April 2026, pre-JOSS submission). Option 3 (MCP-over-MCP Orchestra) is the preferred path for post-JOSS cross-project integration.
---

## Orchestra Implementation Plan

> **Added 2026-03-05.** Based on review of ChatGPT recommendations — concept adopted, timeline revised to respect JOSS 6-month public history requirement.

### Goal
Automate multi-step causal reasoning by composing RegNetAgents + CASCADE via the MCP protocol (Option 3). A single natural language query triggers an orchestrated workflow across both systems, producing integrated causal analysis neither tool can provide alone.

### Dependencies
- RegNetAgents v1.0.1 ✅
- CASCADE v1.0.0 (target March 26, 2026)

### Architecture
Three layers:
- **Decision layer** — LangGraph DAG classifies gene, routes analysis, selects composite tools
- **Evidence layer** — parallel MCP calls to CASCADE (perturbation, PPI, LINCS, DepMap) and RegNetAgents (network topology, pathway enrichment, domain agents)
- **Explanation layer** — synthesis across both systems; LLM narrative optional

See Option 3 architecture diagram above for full detail.

### Milestones

| Month | Milestone |
|---|---|
| March 2026 | Create \ repo, make public immediately (starts 6-month clock) |
| April 2026 | Phase 1: MCP client foundation — spawn CASCADE + RegNetAgents as subprocesses, verify round-trip |
| April 2026 | Phase 2: Prototype \ (APC→CTNNB1 proof of concept) |
| May 2026 | Phase 3: Full composite tools (\, \, \) |
| May–June 2026 | Run 2–3 biological use cases: BRD4→MYC, TP53 upstream causal chain, TF ranking for CRC panel |
| June 2026 | Phase 4: LangGraph DAG orchestration + synthesis |
| June 2026 | Collect metrics: speedup vs manual chaining, reproducibility, cross-source validation consistency |
| July 2026 | Phase 5: MCP server wrapper, test suite, documentation, examples (Python script + MCP client) |
| August 2026 | Draft paper.md: intro, architecture, use cases, results, discussion, conclusion |
| August 2026 | Figures: architecture diagram, workflow DAG, validation table |
| September 2026 | Tag v1.0.0, finalize README, Zenodo archive |
| September 2026 | JOSS submission (6-month public history satisfied) |

### Why September JOSS — Not July
JOSS requires 6 months of public development history. If the Orchestra repo goes public in March 2026, the earliest eligible JOSS submission is September 2026. A July submission would be rejected.

### Notes
- No changes to RegNetAgents or CASCADE codebases — Orchestra connects to both via MCP protocol only
- If either child server is unavailable, Orchestra degrades gracefully using whichever server is running
- This is a standalone JOSS paper — novel contribution is MCP server composition, first published example in bioinformatics

