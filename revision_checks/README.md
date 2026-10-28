# DACH1/CESC revision checks (reviewer-style critique, 2026-10-17/18)

Robustness checks for the DACH1/CESC preprint (rs-10386608), reported in version 3. Original
analysis files (`dach1_partial_correlation.py`, `dach1_partial_corr_cesc.csv` at the repo root)
are untouched.

## Run order (from `c:\Dev\Orchestra`)

| Step | Script | Writes |
|---|---|---|
| 8 | `item8_pubmed_recheck.py` | `item8_pubmed_results.tsv` |
| 0 | `fetch_target_expression.py` | `target_expression_cesc.csv` |
| 1+2 | `purity_and_expression_partials.py <absolute.txt>` | `purity_absolute_cesc.csv`, `item2_purity_results.csv`, `item1_expression_results.csv` |
| 1+3 prep | `build_gene_promoter_matrix.py <meth.tsv.gz> <manifest.csv.gz> <dach1_partial_corr_cesc.csv>` | `gene_promoter_meth_cesc_{A,B}*.csv.gz`, `promoter_definitions.json` |
| 1+3 | `probe_level_and_null.py` | `item1_probe_level_results.csv`, `item3_null_summary.csv`, `item3_null_all_genes_A.csv` |
| 5 diag | `diag_pax1_dcc.py` | `cbioportal_extra_meth_cesc.csv` |
| 4+5+6 | `items4_5_6_ci_burden_bimodal.py` | `item4_ci_table.csv`, `item4_sox8_difference_tests.csv`, `item5_burden_pax1.csv`, `item6_dach1_low_vs_high.csv` |

## Downloaded inputs (not stored here; large, public)

- ABSOLUTE purity: GDC `https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5` (PanCanAtlas)
- Probe-level HM450 betas: UCSC Xena GDC hub `TCGA-CESC.methylation450.tsv.gz` (~1 GB, hg38)
- Illumina manifest: GEO GPL13534 `GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz`

## Software versions

Python 3.12.8, SciPy 1.16.3, NumPy 2.3.5, pandas 2.3.3. The scripts do not record versions; these are the versions installed on the analysis machine. On 2026-10-28 the probe-level chain was rerun from the re-downloaded `TCGA-CESC.methylation450.tsv.gz` under these versions (`build_gene_promoter_matrix.py`, `probe_level_and_null.py`, `screening/gate3_matched_tf_null_dach1.py`, `items4_5_6_ci_burden_bimodal.py`). `promoter_definitions.json` (as parsed), `item1_probe_level_results.csv`, `item3_null_summary.csv`, `gate3_results.json`, `gate3_null_rhos.csv` and the item 4 to 6 tables reproduced unchanged; `item3_null_all_genes_A.csv` has identical values with rows in a different order.

## Done 2026-10-18: RASSF1 rerun (`item5_rassf1_rerun.py`, outputs `item5_rassf1_rerun_{orig23,rassf1_24}.json`)

Versions: Orchestra v1.3.11-3-gd090671, RegNetAgents v1.3.0-1-g54b1cf8, CASCADE v1.5.0-2-g281bb03. Settings: epithelial_cell, tcga_network=cesc, contexts cervical cancer + head and neck squamous.
- `orig23` (RASSF1A omitted) reproduces Table 1 / Supp. Table S1 exactly: same 10 drivers, ranks, overlaps, enrichments, Fisher p and BH p-adj.
- `rassf1_24` (RASSF1 added): RASSF1 resolves; same 10 drivers in the same order, same overlap genes; RASSF1 is in no driver's overlap set. Enrichment and p values shift slightly (e.g. DACH1 8.42 -> 8.07, p-adj 2.82e-3 -> 3.07e-3). Regulators tested 616 -> 649. Still significant.

## Not done

- Random-TF-matched null (item 3, second half): recommended skip; note as a limitation.
- Consensus purity (CPE, Aran 2015) as a second purity source: ABSOLUTE showed no effect.

## Bug caught during the work

First probe-level build had an off-by-one in sample columns (recomputed values agreed with cBioPortal at ρ≈0.04).
Fixed in `build_gene_promoter_matrix.py`; all reported probe-level numbers are from the fixed run
(agreement with cBioPortal values 0.67–0.86).

Note: some outputs here (item1_probe_level_results.csv, item3_null_summary.csv, item4_ci_table.csv, item4_sox8_difference_tests.csv, item6_dach1_low_vs_high.csv, gene_promoter_meth_cesc_B*.csv.gz) also contain a second probe definition ("B": TSS1500 + TSS200 + first-exon probes) that is not used in the paper. The paper reports only the TSS200-based definition (TSS1500 as fallback).

Also run but not reported in the paper: item4_ci_table.csv (95% confidence intervals; the paper cites only the SOX8-TFPI2 interval), item5_burden_pax1.csv (PAX1 added to the burden index) and item6_dach1_low_vs_high.csv (DACH1-low vs -high comparison). They are kept as a record.
