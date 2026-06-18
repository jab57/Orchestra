---
title: 'Orchestra: MCP-Level Composition of Bioinformatics Servers for Multi-System Causal Reasoning'
tags:
  - Python
  - gene regulatory networks
  - causal reasoning
  - Model Context Protocol
  - LangGraph
  - therapeutic target discovery
authors:
  - name: Jose Bird
    orcid: 0009-0006-2744-0606
    affiliation: 1
affiliations:
  - name: Bird AI Solutions
    index: 1
date: 1 September 2026  # TODO: update to actual submission date
archive_doi: 10.5281/zenodo.TBD
bibliography: paper.bib
---

# Summary

Orchestra is a Python orchestration layer that composes two specialized bioinformatics MCP servers — RegNetAgents [@bird2026regnetagents] and CASCADE [@bird2026cascade] — via the Model Context Protocol [@mcp], enabling multi-system causal reasoning from a single query. A LangGraph-based workflow [@langgraph] classifies each gene's regulatory role, routes to an appropriate analysis strategy, and executes coordinated tool calls to both child servers in parallel before synthesizing their outputs into a unified evidence report. Orchestra also provides a cross-cell-type comparison path that runs both child servers across multiple biological contexts simultaneously, a tumor-state regulatory rewiring path that compares population-averaged and TCGA ARACNe networks, and a literature novelty assessment path that queries PubMed to place computational findings in the context of published evidence.

Orchestra exposes seven composite tools. `causal_chain_analysis` runs regulatory network analysis (RegNetAgents) and perturbation simulation (CASCADE) in parallel and identifies targets corroborated by both computational network topology and experimental evidence. `validate_therapeutic_targets` ranks upstream regulators by PageRank centrality (RegNetAgents), then validates top candidates against LINCS knockdown data, super-enhancer evidence, and DepMap essentiality (CASCADE). `effector_analysis` handles scaffold proteins with no transcriptional targets by detecting this dead-end condition, finding transcription factor partners via protein-protein interactions, and routing the analysis through the TF partner. `analyze_gene_signature` accepts a user-supplied gene list, applies Fisher's exact test enrichment via RegNetAgents to identify the master regulators most significantly overrepresented in the signature, and validates the top candidates with CASCADE perturbation analysis — bridging expression-based pattern discovery and mechanistic perturbation evidence. `compare_cell_contexts` runs both child servers across multiple user-supplied cell types in parallel and classifies each evidence source as conserved, enriched, cell-type-specific, or absent — surfacing which regulatory relationships are robust across biological contexts and which are context-dependent. `compare_network_contexts` compares a gene's regulatory wiring between population-averaged (GREmLN ARACNe) and tumor-state (TCGA ARACNe) networks, quantifies regulatory rewiring (low/moderate/high based on regulator Jaccard overlap), and validates conserved regulators via CASCADE experimental data — producing a tiered output that distinguishes conserved-and-CASCADE-validated candidates from tumor-acquired regulatory inputs not present in the population-averaged program. `novelty_assessment` queries PubMed via NCBI E-utilities for a gene or gene pair in a plain-text cancer context, returning total hit count split into experimental and computational papers, most-recent publication year, and a structured novelty verdict (novel: <5 papers; emerging: 5–20; established: >20). Unlike the other six tools, `novelty_assessment` does not call RegNetAgents or CASCADE; it adds a literature arm to the synthesis layer that contextualises computational findings within the published record.

# Statement of Need

Compound causal questions in regulatory biology — *"Who drives this gene, and what happens downstream if we inhibit that driver?"* — require combining two methodologically independent evidence types: network topology inferred from transcriptomics and experimental perturbation data from real cellular assays. Neither RegNetAgents nor CASCADE alone can answer these questions; combining them manually requires knowing which tools to call, in which order, and how to interpret cross-system results.

A concrete example is the APC tumor suppressor. APC is a scaffold protein with no transcriptional targets — a perturbation query returns empty results with no guidance. Answering the question manually requires seven steps: recognizing the failure mode, querying protein-protein interactions, identifying CTNNB1 as the key Wnt pathway partner, simulating CTNNB1 perturbation separately, enriching the downstream cascade against pathways, and synthesizing a causal explanation. Orchestra automates all seven steps from a single query.

General bioinformatics workflow managers (Nextflow, Snakemake) orchestrate file-based pipelines and do not compose interactive AI tool servers or perform gene-role-aware routing. No existing tool implements protocol-level composition of independent bioinformatics MCP servers with cross-system evidence synthesis.

# Architecture

Orchestra follows a three-layer design (\autoref{fig:architecture}). The **decision layer** classifies each gene's regulatory role — master regulator, transcription factor, effector, or isolated — via a CASCADE metadata call and routes to one of seven analysis paths. The **evidence layer** executes coordinated MCP tool calls to RegNetAgents and CASCADE, always in parallel where calls are independent, and accumulates results in a typed LangGraph state object. The **explanation layer** synthesizes outputs from both systems: for each candidate target, it counts how many independent sources agree — network topology rank and pathway membership from RegNetAgents; LINCS knockdown, DepMap essentiality, super-enhancer status, DoRothEA TF confidence, STRING PPI, and cBioPortal tumor expression from CASCADE's pre-synthesized evidence block. An optional LLM synthesis node narrates the scored evidence table without altering the structured output.

![Orchestra architecture. Claude Desktop sends a query to the Orchestra MCP server. The decision layer classifies the gene and selects an analysis path. The evidence layer executes parallel MCP calls to RegNetAgents and CASCADE. The explanation layer identifies cross-system hits — targets corroborated by both network topology and experimental evidence — and returns a structured report.\label{fig:architecture}](figure_architecture.png)

# Functionality

Orchestra is installed and launched as follows:

```bash
pip install -r requirements.txt
cp .env.example .env   # configure paths to RegNetAgents and CASCADE
python orchestra_mcp_server.py
```

Once running, any MCP-compatible client can call Orchestra's seven composite tools.

For a transcription factor, `causal_chain_analysis(gene="TP53", cell_type="epithelial_cell")` classifies TP53 as a master regulator, executes RegNetAgents network analysis and CASCADE perturbation simulation in parallel, and reports targets corroborated by both systems — genes appearing in RegNetAgents' downstream network topology AND in CASCADE's multi-source experimental evidence. These cross-system hits carry higher confidence than targets supported by either system alone.

For a therapeutic target query, `validate_therapeutic_targets(gene="MYC", cell_type="cd4_t_cells")` surfaces BRD4 via CASCADE's super-enhancer analysis: MYC has super-enhancers across 32 cell types, indicating BET inhibitor sensitivity [@loven2013]. RegNetAgents simultaneously returns MYC's classical TF regulatory network context. Orchestra's synthesis reports BRD4 with CASCADE epigenetic support and notes its absence from the ARACNe-inferred TF network — a biologically informative finding, because BRD4 acts through chromatin-level co-activation rather than direct mRNA regulation. Neither child server produces this combined view. The BRD4→MYC relationship is experimentally confirmed, and BET inhibitors targeting this axis are in clinical trials [@loven2013].

For an effector gene, `effector_analysis(gene="APC", cell_type="epithelial_cell")` detects the empty perturbation condition automatically, identifies CTNNB1 as the highest-influence TF partner via STRING PPI, runs CASCADE perturbation analysis and RegNetAgents pathway enrichment in parallel on CTNNB1, and returns an integrated APC→CTNNB1→Wnt causal explanation. This path is unreachable by either child server: CASCADE dead-ends on empty results; RegNetAgents has no PPI data to bridge from scaffold protein to TF partner.

For a gene signature, `analyze_gene_signature(genes=["MYC","TP53","CDKN1A",...], cell_type="epithelial_cell")` applies Fisher's exact test enrichment against the RegNetAgents regulatory network to rank master regulators by how significantly they are overrepresented as upstream drivers of the input gene set. Orchestra then passes the top three transcription factors to CASCADE for independent perturbation validation, and reports which regulators are supported by both enrichment statistics and experimental evidence. This path addresses a common discovery scenario — a differentially expressed gene list from an experiment — where the question is not "what does this gene do?" but "which regulator is driving this whole pattern?"

For tissue specificity assessment, `compare_cell_contexts(gene="MYC", cell_types=["epithelial_cell","cd4_t_cells","nk_cells"])` runs both child servers in parallel for each cell type (2N total MCP calls) and classifies each of the seven evidence sources as conserved (present in ≥ ⌈2/3 × N⌉ cell types), enriched, cell-type-specific, or absent. The output is a heatmap-style evidence table that identifies which regulatory relationships are robust across biological contexts and which are context-dependent — a prerequisite for selecting an appropriate experimental model before therapeutic development.

For literature context, `novelty_assessment(gene="EHF", cancer_context="cervical cancer", gene2="STAT3")` queries PubMed for co-occurrence of EHF and STAT3 in cervical cancer abstracts, returning hit count (3 papers, last 2021), experimental vs. computational split, and a novelty verdict of "novel." When chained after `compare_network_contexts`, this call contextualises CASCADE-validated conserved regulators within the published literature — distinguishing computationally well-characterised targets (e.g. TOP2A: 41 papers, zero experimental) from genuinely underexplored axes where prior work is sparse.

For tumor-state regulatory context, `compare_network_contexts(gene="FOXM1", cell_type="epithelial_cell", cancer_type="hnsc")` queries both RegNetAgents GREmLN ARACNe networks (population-averaged) and TCGA ARACNe networks (tumor-state) for the same gene, computes regulator-level Jaccard overlap to classify rewiring as low (≥60% conserved), moderate (30–60%), or high (<30%), and passes conserved regulators to CASCADE for independent experimental validation. The output assigns each conserved regulator a confidence tier — CASCADE-validated (supported by LINCS, DepMap, super-enhancers, DoRothEA, or cBioPortal) or conserved without experimental support — and lists tumor-acquired regulatory inputs not present in the population-averaged wiring. This path addresses a compound question unreachable by either child server: RegNetAgents holds both GREmLN and TCGA networks but cannot validate conserved regulators experimentally; CASCADE has perturbation evidence but does not contain ARACNe network comparisons.

# Results

Seven biological validation cases confirm that Orchestra's routing, evidence coordination, and synthesis layers function correctly (Table 1).

| Tool | Gene / Input | Cell type | Key result |
|---|---|---|---|
| `effector_analysis` | APC | epithelial_cell | CTNNB1 identified (STRING score 0.999); hub regulator, 310 downstream targets |
| `causal_chain_analysis` | TP53 | epithelial_cell | CDKN1A: 3-source corroboration (DoRothEA, LINCS, STRING); 2 cross-system hits |
| `validate_therapeutic_targets` | MYC | cd4_t_cells | BRD4: 1/7 sources (super-enhancer ✓; absent from ARACNe TF network — expected) |
| `analyze_gene_signature` | 20-gene TP53 stress-response signature | epithelial_cell | KLF5: rank 1 (Fisher p=5×10⁻⁶, 5/9 overlap); DoRothEA-A + super-enhancer (12 GI/epithelial cell types) + DepMap GI essentiality |
| `compare_cell_contexts` | GATA3 | cd4_t_cells, cd8_t_cells, epithelial_cell | DoRothEA TF: conserved (3/3); hub regulator + pathway enrichment: enriched (2/3, cd4+epithelial); CD8 T cells show reduced network centrality — consistent with GATA3's role as Th2/CD4 master regulator |
| `compare_network_contexts` | FOXM1 | epithelial_cell (GREmLN, population-averaged) vs TCGA HNSC (tumor-state) | HIGH rewiring (4.2% conserved regulators; 1/24 HNSC regulators shared with GREmLN); TOP2A conserved in both networks; 23 HNSC tumor-acquired regulatory inputs absent from population-averaged wiring |
| `novelty_assessment` | Cervical cancer therapeutic panel (EHF, STAT3, TOP2A, IDO1, FAP, SERPINB3, RSK4, AK6) | HNSC proxy (HPV/squamous) | EHF→STAT3: Novel (3 papers, last 2021); TOP2A: Established (41 papers, 0 experimental — actionable gap); RSK4/AK6: 0 papers (true white space); STAT3: 37% experimental ratio |

Table: Validation results across Orchestra's composite tools (seven tools implemented; seven validated cases shown).

In the APC case, neither child server alone completes the analysis — CASCADE returns empty perturbation results for the scaffold protein; RegNetAgents has no PPI data to bridge from APC to a transcription factor partner. This is the clearest demonstration that Orchestra's coordinated routing is necessary, not merely convenient. In the MYC case, BRD4 scores 1/7 evidence sources: CASCADE's super-enhancer analysis identifies it as a BET inhibitor target while its absence from RegNetAgents' ARACNe-inferred network is a biologically informative finding — BRD4 acts through chromatin co-activation, not direct mRNA regulation, and is not expected in ARACNe networks. Orchestra presents both views simultaneously, producing a more complete picture than either child server alone.

In the cross-cell-type comparison case, `compare_cell_contexts(gene="GATA3", cell_types=["cd4_t_cells","cd8_t_cells","epithelial_cell"])` runs 6 parallel MCP calls (2 per cell type) and applies conservation scoring across 7 evidence sources. DoRothEA TF regulon confidence is conserved (3/3) — GATA3's identity as a transcription factor is confirmed across all three cell types by an experimentally curated source independent of the mRNA network. PageRank hub status and pathway enrichment are enriched (2/3): both are present in cd4_t_cells and epithelial_cell but absent in cd8_t_cells. This differential is biologically meaningful — GATA3 is the canonical master regulator of Th2 CD4+ helper cell differentiation and plays a substantially reduced role in CD8+ cytotoxic lineages [@zhu2010]. Neither child server alone would surface this context-specificity: a single `comprehensive_gene_analysis(GATA3, cd4_t_cells)` call returns "hub regulator" without any indication that the finding does not extend to CD8. Orchestra's conservation scoring makes the distinction explicit. LINCS knockdown, DepMap essentiality, super-enhancer, and cBioPortal sources are absent in all three contexts, reflecting that those evidence sources are anchored to cancer cell lines rather than the normal immune and epithelial contexts in the GREmLN panel.

In the network context comparison case, `compare_network_contexts(gene="FOXM1", cell_type="epithelial_cell", cancer_type="hnsc")` computes Jaccard overlap between GREmLN epithelial_cell regulators and TCGA HNSC tumor-state regulators. The result is HIGH rewiring (Jaccard 0.04): FOXM1 has one upstream regulator in the population-averaged network (TOP2A) and 24 upstream regulators in the HNSC tumor-state network, with only TOP2A conserved across both contexts (4.2% conserved fraction). The 23 tumor-acquired regulatory inputs that appear exclusively in HNSC represent transcriptional programs active in the cancer context but absent from the population-averaged wiring. TOP2A — a type II topoisomerase and known cell-cycle co-regulator with FOXM1 — is identified as the sole conserved regulator. CASCADE corroboration for TOP2A includes STRING protein interaction partners not detected at the mRNA level; LINCS knockdown and DepMap essentiality evidence is absent, correctly assigning TOP2A the `conserved_not_validated` tier rather than `conserved_cascade_validated`. Neither child server alone produces this comparison: RegNetAgents holds both GREmLN and TCGA ARACNe networks but has no perturbation validation; CASCADE has experimental evidence but does not implement multi-network regulatory comparisons.

In the gene signature case, 20 canonical TP53 stress-response target genes are submitted as input. Fisher enrichment against the epithelial_cell ARACNe regulon ranks KLF5 first (p=5×10⁻⁶; 5 of 9 network-resolved genes fall within its regulon). CASCADE independently validates KLF5 via three orthogonal evidence types: DoRothEA confidence A (the highest tier, combining literature, ChIP-seq, and motif evidence), super-enhancer associations in 12 gastrointestinal and epithelial cell types (colon crypt, esophagus), and DepMap essentiality in bowel and GI lineages. The second- and third-ranked regulators by Fisher p-value (SERTAD2, NCOA7) lack DoRothEA validation and have no GI-relevant DepMap essentiality, illustrating how cross-system scoring correctly deprioritizes high-enrichment candidates that are not independently confirmed as transcription factors.

# Limitations

RegNetAgents regulatory networks are inferred from the CellxGene corpus via GREmLN [@zhang2026gremln], which includes both healthy and disease/cancer-infiltrating cells in heterogeneous proportions. CASCADE's experimental sources (LINCS L1000, DepMap) derive from cancer cell lines. These contexts are not matched: corroboration between systems reflects methodological independence, not biological equivalence. Cross-system agreement is a hypothesis generator and should not substitute for experimental validation in a matched biological context.

TF partner selection in the effector path is heuristic — Orchestra selects the PPI partner with the highest downstream target count. This works reliably for canonical scaffold genes like APC, where CTNNB1 dominates unambiguously. For scaffold proteins with competing TF partners of similar network centrality, the heuristic may not select the most biologically relevant partner, and domain knowledge may be required.

Fisher enrichment in the `analyze_gene_signature` path is sensitive to gene list quality and size. Short or low-quality input lists (fewer than ~20 genes) may produce unreliable enrichment statistics, and only the top three candidate regulators are carried forward for CASCADE validation. Users should treat the ranked output as a hypothesis-generation tool rather than a definitive identification of causal drivers.

Conservation classification in `compare_cell_contexts` uses a 2/3-majority threshold (minimum 2 cell types), which is a heuristic choice. With small N (two or three cell types), the conserved and enriched categories collapse — users should interpret conservation labels in the context of how many cell types were compared. The analysis is also limited by the cell types available in the GREmLN corpus (ten types); tumor-specific or rare cell types are not currently supported.

The `compare_network_contexts` tool covers eight TCGA epithelial cancer types (brca, coad, hnsc, luad, lusc, ov, prad, ucec); hematological malignancies, sarcomas, and cervical squamous carcinoma (CESC) are not available. HNSC (head/neck squamous) is the closest available TCGA proxy for cervical squamous carcinoma owing to shared HPV etiology and squamous histology, but it is not a matched context. Rewiring classification thresholds (Jaccard ≥ 0.6 = low, 0.3–0.6 = moderate, < 0.3 = high) are heuristic; users querying genes with small regulon sizes in either network should interpret rewiring classification cautiously, as low absolute regulator counts make Jaccard overlap sensitive to individual gene-network inclusion decisions.

# Software Availability

Orchestra is available at [https://github.com/jab57/Orchestra](https://github.com/jab57/Orchestra) under the MIT license. The repository includes 236 unit tests and 2 integration tests covering workflow routing, cross-system synthesis, effector path TF partner selection, gene signature enrichment, cross-cell-type conservation scoring, GREmLN vs. TCGA network comparison, PubMed novelty assessment (mocked HTTP), and graceful degradation when one child server is unavailable, with continuous integration via GitHub Actions. Bug reports and feature requests are tracked via the GitHub issue tracker.

# AI Usage Disclosure

Development of Orchestra was assisted by Claude Code (Anthropic), an AI coding tool. The AI assistant was used for code generation, refactoring, test writing, and documentation drafting. All AI-generated code and text were reviewed, tested, and validated by the human author. This paper was drafted collaboratively with AI assistance and reviewed for accuracy by the author.

# Acknowledgements

Orchestra builds on RegNetAgents [@bird2026regnetagents] and CASCADE [@bird2026cascade]. We acknowledge the GREmLN development team at the Chan Zuckerberg Initiative AI for the pre-trained gene embeddings and regulatory networks underlying RegNetAgents. Orchestra uses LangGraph for workflow orchestration and the MCP Python SDK for protocol-level server composition.

# References
