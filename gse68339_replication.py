"""
GSE68339 independent-cohort replication of the DACH1/CESC methylation-
correlation finding (Bird 2026, Research Square, doi:10.21203/rs.3.rs-10386608).

Tests whether DACH1 expression inversely correlates with CADM1/ESR1/SLIT2
promoter methylation (originally found in TCGA CESC, n=253) in GSE68339
cohort 2 (n=121, the "Oslo cohort", same-platform matched methylation +
expression). EDNRB is the negative control. Also computes the burden-index
partial Spearman correlation for ESR1, using the same residuals-on-ranks
method as dach1_partial_correlation.py.

Data prerequisite (not committed -- re-fetch on demand):
  - GSE68339-GPL10558_series_matrix.txt.gz (expression, ~21MB)
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68339/matrix/GSE68339-GPL10558_series_matrix.txt.gz
  - GSE68339-GPL13534_series_matrix.txt.gz (methylation, ~1GB)
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68339/matrix/GSE68339-GPL13534_series_matrix.txt.gz
  - GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz (Illumina 450K manifest, ~64MB)
    https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/suppl/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz
  Decompress all three into RAW_DATA_DIR below before running.

Run from c:\\Dev\\Orchestra:
    python gse68339_replication.py

Outputs (committed, small):
  gse68339_cohort2_sample_map.tsv       -- 121 matched patient IDs (expr GSM <-> meth GSM)
  gse68339_probe_selection.json         -- TSS200/TSS1500 probe-to-gene mapping used
  gse68339_cohort2_per_sample.tsv       -- DACH1 expression + 4 target + 7 burden genes' methylation
  gse68339_results_comparison.csv       -- TCGA-vs-GSE68339 comparison table
  gse68339_burden_partial_result.csv    -- ESR1 burden-controlled partial correlation

Verified 2026-09-12: every re-derived statistic reproduces the previously
documented numbers (ESR1 naive rho=-0.358 p=5.6e-05, partial rho=-0.356
p=6.2e-05; CADM1/SLIT2 trend but miss BH significance at n=121; EDNRB null
in both cohorts) -- see FUTURE_ROADMAP.md Open Item 2 for full context.
"""

import csv
import json
import os

import numpy as np
from scipy import stats

RAW_DATA_DIR = os.environ.get("GSE68339_RAW_DIR", ".")
EXPR_MATRIX = os.path.join(RAW_DATA_DIR, "GSE68339-GPL10558_series_matrix.txt")
METH_MATRIX = os.path.join(RAW_DATA_DIR, "GSE68339-GPL13534_series_matrix.txt")
MANIFEST = os.path.join(RAW_DATA_DIR, "GPL13534_manifest.csv")

TARGET_GENES = ["CADM1", "ESR1", "SLIT2", "EDNRB"]
BH_TARGETS = ["CADM1", "ESR1", "SLIT2"]  # EDNRB is the negative control, excluded from correction
BURDEN_GENES = ["RARB", "DAPK1", "CDH1", "MGMT", "SOCS1", "CDKN2A", "MAL"]
DACH1_PROBE = "ILMN_1755741"  # verified: only one of 3 HT-12 V4 DACH1 probes with signal above background


def _read_sample_metadata(matrix_path: str) -> dict:
    """Parse !Sample_* header lines from a GEO series matrix file."""
    meta: dict[str, list[str]] = {}
    with open(matrix_path, encoding="latin-1") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_"):
                parts = line.rstrip("\n").split("\t")
                key = parts[0]
                vals = [v.strip('"') for v in parts[1:]]
                meta.setdefault(key, [])
                meta[key] = meta[key] + vals if key in meta and meta[key] else vals
    return meta


def build_sample_map() -> list[dict]:
    """Identify the 121 GSE68339 cohort-2 patients with matched expression + methylation."""
    expr_meta = _read_sample_metadata(EXPR_MATRIX)
    meth_meta = _read_sample_metadata(METH_MATRIX)

    expr_gsms = expr_meta["!Sample_geo_accession"]
    expr_titles = expr_meta["!Sample_title"]
    expr_by_pid = {}
    for gsm, title in zip(expr_gsms, expr_titles):
        pid = title.split("[")[1].rstrip("]")
        expr_by_pid[pid] = gsm

    meth_gsms = meth_meta["!Sample_geo_accession"]
    meth_titles = meth_meta["!Sample_title"]
    cohort_rows = [v for v in meth_meta if v.startswith("!Sample_characteristics_ch1")]
    # cohort characteristic is one of the repeated characteristics_ch1 rows; find it by content
    meth_cohort = None
    for i in range(0, len(meth_meta.get("!Sample_characteristics_ch1", [])), len(meth_gsms)):
        pass  # metadata already flattened above; re-derive cohort labels directly:
    # Simpler: re-scan file for the specific characteristics line containing "cohort:"
    with open(METH_MATRIX, encoding="latin-1") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_characteristics_ch1") and "cohort:" in line:
                meth_cohort = [v.strip('"').replace("cohort: ", "") for v in line.rstrip("\n").split("\t")[1:]]
                break
    assert meth_cohort is not None, "could not find cohort characteristic row"

    meth_by_pid_cohort2 = {}
    for gsm, title, cohort in zip(meth_gsms, meth_titles, meth_cohort):
        if cohort != "cohort 2":
            continue
        pid = title.split("(")[1].rstrip(")")
        meth_by_pid_cohort2[pid] = gsm

    matched = sorted(set(expr_by_pid) & set(meth_by_pid_cohort2))
    assert len(matched) == 121, f"expected 121 matched patients, got {len(matched)}"

    return [
        {
            "patient_id": pid,
            "expression_gsm_GPL10558": expr_by_pid[pid],
            "methylation_gsm_GPL13534": meth_by_pid_cohort2[pid],
        }
        for pid in matched
    ]


def select_probes() -> dict:
    """For each target/burden gene, select TSS200 probes (TSS1500 fallback if none)."""
    all_genes = set(TARGET_GENES + BURDEN_GENES)
    probes_by_gene_region = {g: {"TSS200": set(), "TSS1500": set()} for g in all_genes}

    with open(MANIFEST, encoding="latin-1") as f:
        for _ in range(7):
            next(f)
        reader = csv.reader(f)
        header = next(reader)
        idx_ilmn = header.index("IlmnID")
        idx_name = header.index("UCSC_RefGene_Name")
        idx_group = header.index("UCSC_RefGene_Group")
        for row in reader:
            if len(row) <= idx_group:
                continue
            ilmn = row[idx_ilmn]
            names = row[idx_name].split(";") if row[idx_name] else []
            groups = row[idx_group].split(";") if row[idx_group] else []
            for gene, region in zip(names, groups):
                if gene in all_genes and region in ("TSS200", "TSS1500"):
                    probes_by_gene_region[gene][region].add(ilmn)

    selected = {}
    for g in TARGET_GENES + BURDEN_GENES:
        t200 = sorted(probes_by_gene_region[g]["TSS200"])
        t1500 = sorted(probes_by_gene_region[g]["TSS1500"])
        if t200:
            selected[g] = {"probes": t200, "region_used": "TSS200"}
        else:
            selected[g] = {"probes": t1500, "region_used": "TSS1500_fallback"}

    out = {
        "target_genes": TARGET_GENES,
        "burden_genes": BURDEN_GENES,
        "selection_rule": (
            "average beta across all distinct probes annotated TSS200 for the gene; "
            "fall back to TSS1500 only if the gene has zero TSS200 probes"
        ),
        "manifest_source": "GPL13534_HumanMethylation450_15017482_v.1.1.csv (GEO platform GPL13534 supplementary file)",
        "selected_probes": selected,
    }
    with open("gse68339_probe_selection.json", "w") as f:
        json.dump(out, f, indent=2)
    return selected


def build_per_sample_table(patients: list[dict], probe_sel: dict) -> list[dict]:
    """Extract DACH1 expression + per-gene averaged methylation for the 121 matched patients."""
    expr_meta = _read_sample_metadata(EXPR_MATRIX)
    expr_gsm_order = expr_meta["!Sample_geo_accession"]
    dach1_vals = None
    with open(EXPR_MATRIX, encoding="latin-1") as f:
        for line in f:
            if line.startswith(f'"{DACH1_PROBE}"'):
                dach1_vals = [float(v) for v in line.rstrip("\n").split("\t")[1:]]
                break
    assert dach1_vals is not None, f"{DACH1_PROBE} not found in expression matrix"
    expr_gsm_to_dach1 = dict(zip(expr_gsm_order, dach1_vals))

    meth_meta = _read_sample_metadata(METH_MATRIX)
    meth_gsm_order = meth_meta["!Sample_geo_accession"]
    meth_gsm_to_colidx = {gsm: i for i, gsm in enumerate(meth_gsm_order)}

    all_probes = sorted({p for g in probe_sel.values() for p in g["probes"]})
    probe_set = set(all_probes)
    probe_values: dict[str, list] = {}
    with open(METH_MATRIX, encoding="latin-1") as f:
        for line in f:
            if not line.startswith('"cg'):
                continue
            pid = line[1 : line.index('"', 1)]
            if pid not in probe_set:
                continue
            vals = line.rstrip("\n").split("\t")[1:]
            probe_values[pid] = [float(v) if v not in ("", "NA", "null") else None for v in vals]
            if len(probe_values) == len(all_probes):
                break

    genes = list(probe_sel.keys())
    rows = []
    for p in patients:
        colidx = meth_gsm_to_colidx[p["methylation_gsm_GPL13534"]]
        row = {
            "patient_id": p["patient_id"],
            "expression_gsm": p["expression_gsm_GPL10558"],
            "methylation_gsm": p["methylation_gsm_GPL13534"],
            "DACH1_expression": expr_gsm_to_dach1[p["expression_gsm_GPL10558"]],
        }
        for g in genes:
            vals = [probe_values[pr][colidx] for pr in probe_sel[g]["probes"]]
            vals_clean = [v for v in vals if v is not None]
            row[g] = sum(vals_clean) / len(vals_clean) if vals_clean else None
        rows.append(row)

    fieldnames = ["patient_id", "expression_gsm", "methylation_gsm", "DACH1_expression"] + genes
    with open("gse68339_cohort2_per_sample.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)
    return rows


def _partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[float, float]:
    """Partial Spearman of x,y controlling for z. Residuals-on-ranks method,
    identical to dach1_partial_correlation.py's _partial_spearman."""
    rx, ry, rz = stats.rankdata(x), stats.rankdata(y), stats.rankdata(z)

    def _resid(a, b):
        lr = stats.linregress(b, a)
        return a - (float(lr[0]) * b + float(lr[1]))

    r, p = stats.pearsonr(_resid(rx, rz), _resid(ry, rz))
    return float(r), float(p)


def run_statistics(rows: list[dict]) -> None:
    dach1 = np.array([r["DACH1_expression"] for r in rows])

    results = {}
    for t in TARGET_GENES:
        vals = np.array([r[t] for r in rows])
        rho, p = stats.spearmanr(dach1, vals)
        results[t] = {"rho": float(rho), "p": float(p), "n": len(vals)}

    pvals_sorted = sorted([(t, results[t]["p"]) for t in BH_TARGETS], key=lambda x: x[1])
    m = len(pvals_sorted)
    padj, prev = {}, 1.0
    for i in range(m - 1, -1, -1):
        t, p = pvals_sorted[i]
        val = min(prev, p * m / (i + 1))
        padj[t] = val
        prev = val

    tcga_documented = {
        "CADM1": {"rho": -0.249, "p": 6.4e-05, "n": 253},
        "ESR1": {"rho": -0.248, "p": 6.6e-05, "n": 253},
        "SLIT2": {"rho": -0.263, "p": 2.4e-05, "n": 253},
        "EDNRB": {"rho": -0.157, "p": 1.2e-02, "n": 253},
    }
    with open("gse68339_results_comparison.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "tcga_rho", "tcga_p", "tcga_n", "gse68339_rho", "gse68339_p", "gse68339_p_adj_bh", "gse68339_n", "replicates"])
        for t in TARGET_GENES:
            r, tcga = results[t], tcga_documented[t]
            adj = padj.get(t)
            replicates = "yes" if (t in padj and padj[t] < 0.05) else ("n/a - negative control" if t == "EDNRB" else "no - trend only")
            w.writerow([t, tcga["rho"], tcga["p"], tcga["n"], round(r["rho"], 4), f"{r['p']:.3e}", f"{adj:.3e}" if adj else "n/a", r["n"], replicates])
            print(f"{t}: rho={r['rho']:+.4f} p={r['p']:.3e}" + (f" p_adj={adj:.3e}" if adj else ""))

    burden = np.array([np.mean([r[g] for g in BURDEN_GENES]) for r in rows])
    esr1 = np.array([r["ESR1"] for r in rows])
    part_r, part_p = _partial_spearman(dach1, esr1, burden)
    naive_r, naive_p = stats.spearmanr(dach1, esr1)
    with open("gse68339_burden_partial_result.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "n", "naive_rho", "naive_p", "partial_rho", "partial_p", "burden_genes"])
        w.writerow(["ESR1", len(rows), round(float(naive_r), 4), f"{naive_p:.3e}", round(part_r, 4), f"{part_p:.3e}", ";".join(BURDEN_GENES)])
    print(f"ESR1 partial (burden-controlled): rho={part_r:+.4f} p={part_p:.3e}")


def main() -> None:
    print("Building sample map...")
    patients = build_sample_map()
    with open("gse68339_cohort2_sample_map.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["patient_id", "expression_gsm_GPL10558", "methylation_gsm_GPL13534"], delimiter="\t")
        w.writeheader()
        w.writerows(patients)
    print(f"  {len(patients)} matched patients")

    print("Selecting probes...")
    probe_sel = select_probes()

    print("Building per-sample table (this reads the full methylation matrix once)...")
    rows = build_per_sample_table(patients, probe_sel)

    print("Running statistics...")
    run_statistics(rows)


if __name__ == "__main__":
    main()
