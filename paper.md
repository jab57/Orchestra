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
date: 1 September 2026
archive_doi: 10.5281/zenodo.TBD
bibliography: paper.bib
---

# Summary

Orchestra is a Python orchestration layer that composes two specialized bioinformatics MCP servers — RegNetAgents [@bird2026regnetagents] and CASCADE [@bird2026cascade] — via the Model Context Protocol [@mcp], enabling multi-system causal reasoning from a single query. A LangGraph-based workflow [@langgraph] classifies each gene's regulatory role, routes to an appropriate analysis strategy, and executes coordinated tool calls to both child servers in parallel before synthesizing their outputs into a unified evidence report.

Orchestra exposes three composite tools that require both child servers and cannot be replicated by either alone. `causal_chain_analysis` runs regulatory network analysis (RegNetAgents) and perturbation simulation (CASCADE) in parallel and identifies targets corroborated by both computational network topology and experimental evidence. `validate_therapeutic_targets` ranks upstream regulators by PageRank centrality (RegNetAgents), then validates top candidates against LINCS knockdown data, super-enhancer evidence, and DepMap essentiality (CASCADE). `effector_analysis` handles scaffold proteins with no transcriptional targets by detecting this dead-end condition, finding transcription factor partners via protein-protein interactions, and routing the analysis through the TF partner.

# Statement of Need

Compound causal questions in regulatory biology — *"Who drives this gene, and what happens downstream if we inhibit that driver?"* — require combining two methodologically independent evidence types: network topology inferred from transcriptomics and experimental perturbation data from real cellular assays. Neither RegNetAgents nor CASCADE alone can answer these questions; combining them manually requires knowing which tools to call, in which order, and how to interpret cross-system results.

A concrete example is the APC tumor suppressor. APC is a scaffold protein with no transcriptional targets — a perturbation query returns empty results with no guidance. Answering the question manually requires seven steps: recognizing the failure mode, querying protein-protein interactions, identifying CTNNB1 as the key Wnt pathway partner, simulating CTNNB1 perturbation separately, enriching the downstream cascade against pathways, and synthesizing a causal explanation. Orchestra automates all seven steps from a single query.

General bioinformatics workflow managers (Nextflow, Snakemake) orchestrate file-based pipelines and do not compose interactive AI tool servers or perform gene-role-aware routing. No existing tool implements protocol-level composition of independent bioinformatics MCP servers with cross-system evidence synthesis.

# Architecture

Orchestra follows a three-layer design (\autoref{fig:architecture}). The **decision layer** classifies each gene's regulatory role — master regulator, transcription factor, effector, or isolated — via a CASCADE metadata call and routes to one of three analysis paths. The **evidence layer** executes coordinated MCP tool calls to RegNetAgents and CASCADE, always in parallel where calls are independent, and accumulates results in a typed LangGraph state object. The **explanation layer** synthesizes outputs from both systems: for each candidate target, it counts how many independent sources agree — network topology rank and pathway membership from RegNetAgents; LINCS knockdown, DepMap essentiality, super-enhancer status, DoRothEA TF confidence, STRING PPI, and cBioPortal tumor expression from CASCADE's pre-synthesized evidence block. An optional LLM synthesis node narrates the scored evidence table without altering the structured output.

![Orchestra architecture. Claude Desktop sends a query to the Orchestra MCP server. The decision layer classifies the gene and selects an analysis path. The evidence layer executes parallel MCP calls to RegNetAgents and CASCADE. The explanation layer identifies cross-system hits — targets corroborated by both network topology and experimental evidence — and returns a structured report.\label{fig:architecture}](figure_architecture.png)

# Functionality

Orchestra is installed and launched as follows:

```bash
pip install -r requirements.txt
cp .env.example .env   # configure paths to RegNetAgents and CASCADE
python orchestra_mcp_server.py
```

Once running, any MCP-compatible client can call Orchestra's three composite tools.

For a transcription factor, `causal_chain_analysis(gene="TP53", cell_type="epithelial_cell")` classifies TP53 as a master regulator, executes RegNetAgents network analysis and CASCADE perturbation simulation in parallel, and reports targets corroborated by both systems — genes appearing in RegNetAgents' downstream network topology AND in CASCADE's multi-source experimental evidence. These cross-system hits carry higher confidence than targets supported by either system alone.

For a therapeutic target query, `validate_therapeutic_targets(gene="MYC", cell_type="cd4_t_cells")` surfaces BRD4 via CASCADE's super-enhancer analysis: MYC has super-enhancers across 32 cell types, indicating BET inhibitor sensitivity [@loven2013]. RegNetAgents simultaneously returns MYC's classical TF regulatory network context. Orchestra's synthesis reports BRD4 with CASCADE epigenetic support and notes its absence from the ARACNe-inferred TF network — a biologically informative finding, because BRD4 acts through chromatin-level co-activation rather than direct mRNA regulation. Neither child server produces this combined view. The BRD4→MYC relationship is experimentally confirmed, and BET inhibitors targeting this axis are in clinical trials [@loven2013].

For an effector gene, `effector_analysis(gene="APC", cell_type="epithelial_cell")` detects the empty perturbation condition automatically, identifies CTNNB1 as the highest-influence TF partner via STRING PPI, runs CASCADE perturbation analysis and RegNetAgents pathway enrichment in parallel on CTNNB1, and returns an integrated APC→CTNNB1→Wnt causal explanation. This path is unreachable by either child server: CASCADE dead-ends on empty results; RegNetAgents has no PPI data to bridge from scaffold protein to TF partner.

# Limitations

RegNetAgents regulatory networks are inferred from the CellxGene corpus via GREmLN [@zhang2026gremln], which includes both healthy and disease/cancer-infiltrating cells in heterogeneous proportions. CASCADE's experimental sources (LINCS L1000, DepMap) derive from cancer cell lines. These contexts are not matched: corroboration between systems reflects methodological independence, not biological equivalence. Cross-system agreement is a hypothesis generator and should not substitute for experimental validation in a matched biological context.

TF partner selection in the effector path is heuristic — Orchestra selects the PPI partner with the highest downstream target count. This works reliably for canonical scaffold genes like APC, where CTNNB1 dominates unambiguously. For scaffold proteins with competing TF partners of similar network centrality, the heuristic may not select the most biologically relevant partner, and domain knowledge may be required.

# Software Availability

Orchestra is available at [https://github.com/jab57/Orchestra](https://github.com/jab57/Orchestra) under the MIT license. The repository includes automated tests covering workflow routing, cross-system synthesis, effector path TF partner selection, and graceful degradation when one child server is unavailable, with continuous integration via GitHub Actions.

# AI Usage Disclosure

Development of Orchestra was assisted by Claude Code (Anthropic), an AI coding tool. The AI assistant was used for code generation, refactoring, test writing, and documentation drafting. All AI-generated code and text were reviewed, tested, and validated by the human author. This paper was drafted collaboratively with AI assistance and reviewed for accuracy by the author.

# Acknowledgements

Orchestra builds on RegNetAgents [@bird2026regnetagents] and CASCADE [@bird2026cascade]. We acknowledge the GREmLN development team at the Chan Zuckerberg Initiative AI for the pre-trained gene embeddings and regulatory networks underlying RegNetAgents. Orchestra uses LangGraph for workflow orchestration and the MCP Python SDK for protocol-level server composition.

# References
