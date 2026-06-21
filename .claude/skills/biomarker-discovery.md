Cross-cancer regulatory biomarker discovery using Orchestra MCP tools.

Identify upstream transcription factor drivers of a silenced gene panel and find which drivers represent transfer opportunities between two cancer contexts — established in one, novel in the other.

## Arguments

- $CANCER_CONTEXT_1: First cancer type (e.g. "breast cancer")
- $CANCER_CONTEXT_2: Second cancer type (e.g. "cervical cancer")
- $CELL_TYPE: GREmLN cell type to use (default: epithelial_cell)

## Workflow — follow these steps in order, do not skip or reorder

### Step 1 — Compile the gene panel (STOP and wait for approval)

Before calling any Orchestra tool, compile a literature-curated list of 20–30 genes recurrently silenced by promoter hypermethylation or other epigenetic mechanisms in both $CANCER_CONTEXT_1 and $CANCER_CONTEXT_2.

For each gene provide:
- Gene symbol
- Evidence type (promoter methylation / LOH / transcriptional silencing)
- One sentence of literature justification

Present the list and **wait for the user to approve or modify it before proceeding to Step 2.**

### Step 2 — Run signature enrichment (one call only)

Call `analyze_gene_signature` with:
- `genes`: the approved gene list from Step 1
- `cell_type`: $CELL_TYPE (default epithelial_cell)

**Do NOT call `causal_chain_analysis`, `comprehensive_gene_analysis`, or any per-gene tool on individual genes from the panel.** That bypasses Fisher enrichment and produces unstatistical overlap results.

### Step 3 — Cross-context novelty (two calls)

For the top 5 ranked TF drivers from Step 2, call `novelty_assessment_batch` twice:
1. `genes`: top 5 drivers, `cancer_context`: "$CANCER_CONTEXT_1"
2. `genes`: top 5 drivers, `cancer_context`: "$CANCER_CONTEXT_2"

### Step 4 — Synthesize the gap table

Present a table with one row per driver:

| Driver | Fisher p-value | CASCADE score | $CANCER_CONTEXT_1 | $CANCER_CONTEXT_2 | Gap classification |
|---|---|---|---|---|---|

Gap classification rules:
- **Transfer opportunity**: established (>5 papers) in one context, novel (<5 papers) in the other — highest priority
- **Bilateral novel**: novel in both contexts — hypothesis-generating
- **Bilateral established**: established in both — lower priority for cross-cancer work

### Step 5 — Proposal narrative

For each transfer opportunity, write 2–3 sentences covering:
- What the driver regulates (from Step 2 CASCADE evidence)
- Why the literature gap between the two cancer types is scientifically meaningful
- What experimental validation would look like

Flag all findings as computationally derived, requiring wet-lab confirmation.
