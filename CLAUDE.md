# Orchestra MCP Server — Usage Rules

This project exposes an MCP server (`orchestra_mcp_server.py`) that composes
RegNetAgents and CASCADE for regulatory network analysis. When assisting a user
with Orchestra tools, follow the pipelines and rules below.

## Reference

**Cell types**: cd4_t_cells, cd8_t_cells, cd14_monocytes, cd16_monocytes, nk_cells,
nkt_cells, cd20_b_cells, monocyte-derived_dendritic_cells, erythrocytes, epithelial_cell

**TCGA tumor networks**: blca (bladder), brca (breast), cesc (cervical), coad (colorectal),
hnsc (head/neck SCC), kirc (kidney), lihc (liver), luad (lung adeno), lusc (lung SCC),
ov (ovarian), paad (pancreatic), prad (prostate), stad (stomach), ucec (endometrial)

---

## ⚠️ HARD RULE — pipeline isolation
Each pipeline is self-contained. Do not call tools from another pipeline during
execution. Do not add steps not listed. If additional analysis seems useful,
suggest it to the user after completing the current pipeline — do not run it
automatically.

---

## Pipeline 1a — Compile Gene Panel
**Trigger**: user asks to compile a gene panel, find epigenetically silenced genes,
or says "do not run analysis yet".

Steps:
1. Build a literature-curated list of 20–30 genes recurrently silenced by promoter
   methylation or other epigenetic mechanisms in BOTH cancer contexts specified.
   For each gene: symbol, evidence type, one-sentence justification.
2. Present the list. End your response with:
   "Review this panel and reply 'approved, run analysis' when ready to proceed."
3. Do not call any Orchestra tool. Stop here.

---

## Pipeline 1b — Run Biomarker Discovery
**Trigger**: user says "approved, run analysis" or "proceed" after reviewing a
gene panel in the conversation.

Steps — run exactly these tools in this order, nothing else:
1. Call `analyze_gene_signature` with the panel genes, cell_type, tcga_network
   (if specified), and cancer_contexts set to both cancer contexts being compared.
   Do NOT call per-gene tools on individual panel members.
2. Present the gap table from the built-in cross-context novelty output:
   Driver | p-adj (BH) | CASCADE score | Context 1 papers | Context 2 papers | Classification.
   Use BH-adjusted p as the primary significance criterion, not raw Fisher p.
   Flag any driver marked † (≤ 2 overlapping panel genes) as a provisional finding —
   its fold-enrichment is not robust to gene list variation.
3. For each transfer opportunity or bilateral novel finding write 2–3 sentences:
   what the driver regulates, why the finding is meaningful, what validation would
   look like. Flag all findings as computationally derived.
4. If further analysis of a specific driver seems warranted, suggest it —
   do not run it automatically.

---

## Pipeline 2 — Regulatory Rewiring
**Trigger**: user explicitly asks for regulatory rewiring, how a gene's regulation
changes in a tumor, or which regulators are gained/lost in cancer.

Steps — run exactly these tools in this order, nothing else:
1. Ask for gene and TCGA network if not provided.
2. Call `compare_network_contexts` with gene, cell_type (default: epithelial_cell),
   cancer_type.
3. If the tool returns an error because the gene is not found in the GREmLN network:
   report this to the user, explain that Pipeline 2 requires the gene to be present
   in both the normal tissue (GREmLN) and tumor (TCGA) networks, and suggest
   Pipeline 3 (causal chain analysis) as an alternative.
4. If successful, report: conserved regulators, tumor-acquired regulators, rewiring
   classification (low/moderate/high), CASCADE validation status.
5. If further analysis seems relevant, suggest it — do not run it automatically.

---

## Pipeline 3 — Causal Chain Analysis
**Trigger**: user explicitly asks how gene X drives a phenotype, downstream effects,
causal mechanism, or Pipeline 2 failed due to missing GREmLN coverage.

Steps — run exactly these tools in this order, nothing else:
1. Call `causal_chain_analysis` with gene, cell_type, cancer_context if provided.
2. Summarize causal chain, cross-system hits, discordance flags.
3. If further analysis seems relevant, suggest it — do not run it automatically.

---

## Pipeline 4 — Therapeutic Target Discovery
**Trigger**: user explicitly asks for drug targets, druggable regulators, therapeutic
vulnerabilities.

Steps — run exactly these tools in this order, nothing else:
1. Call `validate_therapeutic_targets` with gene, cell_type, cancer_context if provided.
2. Rank by PageRank + CASCADE validation. Highlight multi-source validated targets.
3. If further analysis seems relevant, suggest it — do not run it automatically.

---

## Pipeline 5 — Gene Signature Analysis
**Trigger**: user provides a gene list and asks which TFs drive it, master regulators
of a pathway, or asks to analyze a specific gene signature.

Steps — run exactly these tools in this order, nothing else:
1. Confirm the gene list with the user before calling anything.
2. Call `analyze_gene_signature` with genes, cell_type, tcga_network (if tumor
   context), cancer_contexts (for novelty gap if provided).
3. Report top TF drivers with Fisher p and BH-adjusted p-values, CASCADE scores,
   and novelty gap if available. Call out any driver marked † (≤ 2 overlapping
   genes) as a provisional finding — do not present its fold-enrichment as reliable.
4. If further analysis seems relevant, suggest it — do not run it automatically.

---

## Pipeline 6 — Cell Context Comparison
**Trigger**: user explicitly asks how a gene behaves across cell types, which cell
types show enriched vs. conserved activity.

Steps — run exactly these tools in this order, nothing else:
1. Call `compare_cell_contexts` with gene and relevant cell types.
2. Present heatmap: conserved / enriched / specific / absent per cell type.
3. If further analysis seems relevant, suggest it — do not run it automatically.

---

## Pipeline 7 — Novelty Assessment
**Trigger**: user explicitly asks how well-studied a gene is in a cancer, or whether
a finding is novel.

Steps — run exactly these tools in this order, nothing else:
1. Call `novelty_assessment` (single) or `novelty_assessment_batch` (multiple) with
   cancer_context.
2. Report PubMed counts, novelty verdict, experimental vs. computational ratio.

---

## Pipeline 8 — Cross-Cancer Rewiring Comparison
**Trigger**: user asks whether a gene's regulatory rewiring is consistent across
cancer types, wants to compare rewiring patterns in multiple tumors, or wants to
identify which cancer types show the most regulatory change for a specific gene.

Steps — run exactly these tools in this order, nothing else:
1. Ask the user to select 2–4 TCGA cancer types (do not run more than 4 — each
   call takes 2–3 minutes and results become difficult to compare beyond this).
   If the user has not specified, suggest biologically related types (e.g. luad +
   lusc for lung, or hnsc + cesc for squamous epithelial) and wait for confirmation
   before proceeding.
2. Call `compare_network_contexts` once per selected cancer type, sequentially.
   Use the same gene and cell_type for each call.
3. Present a summary table:
   Cancer type | Rewiring classification | Conserved regulators (n) |
   Tumor-acquired regulators (n) | CASCADE-validated conserved (n).
4. Identify two cross-cancer patterns separately:
   - Conserved regulators appearing across all tested types are more robust
     baseline candidates than those conserved in only one type.
   - Tumor-acquired regulators appearing in multiple cancer types indicate
     convergent oncogenic rewiring and are higher priority than cancer-type-specific
     ones. Tumor-acquired regulators unique to one cancer type indicate divergent
     rewiring.
5. If further analysis seems relevant, suggest it — do not run it automatically.

---

## Global rules
- Default cell type: epithelial_cell. Ask if unclear.
- For tumor analyses, ask which TCGA network if not specified.
- Never call per-gene tools on members of a gene panel — use `analyze_gene_signature`.
- Use BH-adjusted p (p-adj) as the primary significance criterion whenever
  `analyze_gene_signature` results are discussed. Raw Fisher p may be shown
  alongside but should not be the basis for ranking or conclusions.
- Drivers marked † have ≤ 2 overlapping genes — always flag these explicitly as
  provisional and do not cite their fold-enrichment as evidence of regulatory control.
- Flag all computational findings as requiring wet-lab confirmation.
- After completing any pipeline, suggest relevant follow-up pipelines but do not
  run them unless the user asks.
