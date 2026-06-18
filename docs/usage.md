# Usage

## Running as a Standalone Script

The quickest way to run Orchestra is through the validation script:

```bash
python run_validation.py            # all three validation cases
python run_validation.py apc        # APC effector analysis only
python run_validation.py tp53       # TP53 causal chain only
python run_validation.py brd4       # MYC therapeutic target validation only
```

Results are saved to `outputs/validation_<case>_<timestamp>.txt`.

See `examples/` for focused single-case scripts you can adapt for your own genes.

## Three Composite Tools

### `causal_chain_analysis(gene, cell_type, cancer_context=None)`

Full causal chain: classifies a gene, runs regulatory network analysis (RegNetAgents) and perturbation simulation (CASCADE) in parallel, and synthesizes the results.

The optional `cancer_context` parameter (plain-text, e.g. `"colorectal"`, `"breast cancer"`) enables embedded pair novelty: after synthesis, Orchestra queries PubMed for each of the top 5 identified regulatory edges concurrently and appends a "Regulatory Pair Novelty" table to the report. Omit to skip pair novelty queries.

**Routing:**
- If the gene is a TF or master regulator → parallel RegNetAgents + CASCADE, cross-system corroboration
- If the gene is an effector/scaffold → PPI to find TF partner, then analyze TF partner

**Example — TP53 in epithelial cells:**

```python
import asyncio
from orchestra_langgraph_workflow import OrchestraWorkflow

async def main():
    workflow = OrchestraWorkflow()
    result = await workflow.run_analysis(
        gene="TP53",
        cell_type="epithelial_cell",
        analysis_type="causal_chain",
    )
    print(result["final_report"])

asyncio.run(main())
```

**Example — APC in epithelial cells (effector path):**

```python
result = await workflow.run_analysis(
    gene="APC",
    cell_type="epithelial_cell",
    analysis_type="causal_chain",
)
```

### `validate_therapeutic_targets(gene, cell_type, cancer_context=None)`

Identifies druggable regulators upstream of a gene using three candidate sources:
1. RegNetAgents: upstream regulators ranked by PageRank centrality
2. CASCADE: drug discovery database (super-enhancers, known vulnerabilities)
3. CASCADE STRING PPI: protein-level interactors

Top candidates are validated via CASCADE comprehensive perturbation analysis. Returns a 7-source corroboration table.

The optional `cancer_context` parameter enables the same embedded pair novelty as `causal_chain_analysis` — top 5 regulator→gene edges are queried against PubMed and a "Regulatory Pair Novelty" table is appended.

**Example — therapeutic targets for MYC in T cells:**

```python
result = await workflow.run_analysis(
    gene="MYC",
    cell_type="cd4_t_cells",
    analysis_type="therapeutic_validation",
)
```

### Supported Cell Types

Both child servers support these cell types:

```
cd4_t_cells          cd8_t_cells
cd14_monocytes       cd16_monocytes
nk_cells             nkt_cells
cd20_b_cells         monocyte-derived_dendritic_cells
erythrocytes         epithelial_cell
```

## Interpreting the Output

### Report Structure

The `final_report` field contains a Markdown-formatted report. The `synthesis` field contains the structured dict that the report was generated from — use it for programmatic access.

```python
result = await workflow.run_analysis(...)
report = result["final_report"]         # Markdown string
synthesis = result["synthesis"]          # structured dict
errors = result["errors"]               # {} on success; partial data otherwise
```

### Corroboration Table (validation path)

The validation path generates a 7-source corroboration table:

| Candidate | PageRank | Pathway | LINCS | DepMap | SuperEnhancer | DoRothEA | cBioPortal | Score |
|---|---|---|---|---|---|---|---|---|
| BRD4 | - | - | - | - | ✓ | - | - | **1/7** |

**PageRank, Pathway** — RegNetAgents sources (regulatory network topology)  
**LINCS, DepMap, SuperEnhancer, DoRothEA, cBioPortal** — CASCADE sources (experimental data)

A candidate scoring ≥3/7 has evidence from multiple independent methods. Agreement between RegNetAgents (network topology, inferred from GREmLN/CellxGene) and CASCADE (LINCS/DepMap, cancer cell lines) is methodologically independent — it is a hypothesis generator, not experimental validation.

**Important limitation:** RegNetAgents networks are inferred from heterogeneous cell states (healthy + disease); CASCADE experimental sources are from cancer cell lines. Cross-system agreement does not imply a matched biological context.

### Cross-System Hits (TF path)

The TF path reports genes found in both RegNetAgents downstream targets (network topology) and CASCADE `multi_source_genes` (experimental data). These are the cross-system corroborated targets — the core output that neither child server produces alone.

```
Cross-system hits (N genes in both RegNetAgents topology AND CASCADE experimental data):
- CDKN1A: CASCADE 3 sources (dorothea, lincs, string) + RegNetAgents network target
```

### Graceful Degradation

If one child server is unavailable, Orchestra continues with the other:

```
⚠️ RegNetAgents unavailable — network topology and pathway evidence missing;
   showing CASCADE-only results.
```

The `errors` dict reports which calls failed:

```python
result["errors"]  # e.g. {"network": "TimeoutError: ..."}
```

## LLM Synthesis (Optional)

By default Orchestra returns structured text and lets Claude Desktop (or your own LLM) handle narrative interpretation. To enable a brief biological narrative prepended to each report:

```env
USE_LLM_SYNTHESIS=true
LLM_PROVIDER=ollama          # ollama (default) | anthropic
OLLAMA_MODEL=llama3.1:8b
```

For Anthropic:

```env
LLM_PROVIDER=anthropic
LLM_MODEL=claude-haiku-4-5-20251001
LLM_API_KEY=your-api-key
```

This follows the same pattern as CASCADE (`USE_LLM_INSIGHTS`) and RegNetAgents (`USE_LLM_AGENTS`).
