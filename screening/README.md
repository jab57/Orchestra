# Screening checks for the DACH1/CESC short communication (version 3)

Scripts and outputs behind the enrichment re-checks and the matched-regulator comparisons reported in version 3. Each script ran once; the pass rules were fixed before the runs (see the paper's Methods).

| Script | What it does | Output |
|---|---|---|
| `gate1_stability_and_null.py` | Reproduces the tool's Fisher enrichment for the 24-gene panel, leave-one-panel-gene-out ranks, and the permutation null (2,000 random panels, seed 42) | `gate1_results.json`, `gate1_leave_one_out.csv` |
| `gate3_matched_tf_null_dach1.py` | Compares DACH1's methylation correlations with regulators matched for regulon size and expression (TCGA-CESC) | `gate3_results.json`, `gate3_null_rhos.csv` |
| `gate3_hnsc_replication.py` | Same procedure in TCGA-HNSC | `gate3_hnsc_results.json`, `gate3_hnsc_null_rhos.csv` |
| `gate4_larger_panel.py` | 67-gene panel (the original 24 plus 43 genes); contains the panel | `gate4_larger_panel_results.json` |
| `gate4_run_orchestra_tool.py` | Runs Orchestra's own `analyze_gene_signature` workflow on the 67-gene panel | `gate4_orchestra_tool_output_67genes.json` |

`gate1_stability_and_null.py` is the DACH1/CESC case only; the output files were filtered to that case to match. Run from the repository root, for example `python -I -X utf8 screening/gate1_stability_and_null.py` (about 1 minute); `gate4_larger_panel.py` takes about 3 minutes.

## What is not included
- **Raw data.** The UCSC Xena GDC hub files (`TCGA-CESC.star_tpm.tsv.gz`, `TCGA-HNSC.star_tpm.tsv.gz`, `TCGA-HNSC.methylation450.tsv.gz`) and the Illumina HM450 manifest (GEO platform GPL13534) are not redistributed; the paper's Data availability section gives the SHA-256 values of the three Xena files. The CESC promoter-methylation matrix is built by `revision_checks/build_gene_promoter_matrix.py` and is not included either.
- **The TCGA ARACNe networks** (Giorgi and Alvarez, CC BY-NC-ND 4.0). The scripts read them from a local RegNetAgents install (`C:\Dev\RegNetAgents\models\networks\tcga`); install them with RegNetAgents' `scripts/setup_tcga_networks.py --accept-license`. The output files contain no network edges or regulons, only regulon sizes, overlap counts and test statistics.

## Notes
- File paths in the scripts are the author's local paths (for example `C:\Dev\xena_data`); change them to run elsewhere.
- In `gate3_null_rhos.csv` and `gate3_hnsc_null_rhos.csv` the first column is labelled `Ensembl_ID` but holds gene symbols (one row per regulator; correlations with promoter methylation of each target, regulon size, mean expression). They hold aggregate values only, no per-sample data.
- Software versions are given in the paper's Methods (Orchestra, RegNetAgents, CASCADE `git describe` values).
- All five scripts, as they appear here, were re-run on 2026-10-24 and reproduced the original output files exactly. `gate3_matched_tf_null_dach1.py` needs the CESC promoter matrix built by `revision_checks/build_gene_promoter_matrix.py` (the matrix itself is not included).
