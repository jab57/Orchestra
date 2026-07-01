"""
Single-source correlation analysis for the DACH1 CESC paper.

Produces all numbers needed for Tables 2 and 3 and the SOX8/SALL2
exploratory comparison section. All analyses use squamous-only TCGA CESC
samples (CANCER_TYPE_DETAILED = 'Cervical Squamous Cell Carcinoma').

Run from c:\\Dev\\Orchestra with the project venv active:
    python dach1_partial_correlation.py

Output sections:
  [TABLE 2]  DACH1 naive Spearman vs all 4 regulon targets
  [TABLE 3]  DACH1 naive vs partial Spearman (3 confirmed targets)
  [SOX8]     SOX8 naive Spearman vs its 3 ARACNe targets
  [SALL2]    SALL2 naive Spearman vs its 3 ARACNe targets

A per-sample CSV is written to dach1_partial_corr_cesc.csv for scatter plots.
"""

import csv
import os
import sys

import numpy as np
from dotenv import load_dotenv
from scipy import stats

load_dotenv()
sys.path.insert(0, os.path.dirname(__file__))

from cbioportal_client import (
    _TCGA_STUDY_IDS,
    _discover_profiles,
    _entrez_ids,
    _fetch_values_timed,
    _get_sample_ids,
    _get_sample_ids_by_clinical_attr,
)

TCGA_NETWORK = "cesc"

# DACH1: all 4 ARACNe-predicted targets (Table 2); partial only for 3 confirmed (Table 3)
DACH1_TARGETS_ALL = ["CADM1", "ESR1", "SLIT2", "EDNRB"]
DACH1_TARGETS_PARTIAL = ["CADM1", "ESR1", "SLIT2"]

# Exploratory comparisons
SOX8_TARGETS  = ["DCC", "FHIT", "TFPI2"]
SALL2_TARGETS = ["DCC", "FHIT", "TFPI2"]

# Methylation burden index genes (circularity-free: none overlap DACH1 targets)
BURDEN_GENES = ["CDKN2A", "RARB", "DAPK1", "CDH1", "MGMT", "MAL", "SOCS1", "PAX1"]

SQUAMOUS_VALUE = "Cervical Squamous Cell Carcinoma"


def _partial_spearman(
    x: np.ndarray, y: np.ndarray, z: np.ndarray
) -> tuple[float, float]:
    """
    Partial Spearman correlation of x and y controlling for z.
    Residuals-on-ranks method: regress rank(x) and rank(y) on rank(z) separately,
    then take Pearson r of the two residual vectors.
    """
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)

    def _resid(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        lr = stats.linregress(b, a)
        return a - (float(lr[0]) * b + float(lr[1]))  # slope=lr[0], intercept=lr[1]

    r, p = stats.pearsonr(_resid(rx, rz), _resid(ry, rz))
    return float(r), float(p)


def _naive_row(n: int, r: float, p: float) -> str:
    direction = "Inverse" if r < -0.2 else ("Concordant" if r > 0.2 else "No clear trend")
    return f"  n={n}  rho={r:+.3f}  p={p:.2e}  [{direction}]"


def _run_naive(
    regulator: str,
    targets: list[str],
    expr: dict[str, float],
    meth_cache: dict[str, dict[str, float]],
) -> dict[str, tuple[int, float, float]]:
    """Return {gene: (n, rho, p)} for each target with sufficient data."""
    results: dict[str, tuple[int, float, float]] = {}
    for gene in targets:
        meth = meth_cache.get(gene, {})
        shared = sorted(set(expr) & set(meth))
        if len(shared) < 5:
            print(f"  {gene}: n={len(shared)} — insufficient data (< 5), skipped")
            continue
        x = np.array([expr[s] for s in shared])
        y = np.array([meth[s] for s in shared])
        rho, p = stats.spearmanr(x, y)
        results[gene] = (len(shared), float(rho), float(p))
        print(f"  {gene}:" + _naive_row(len(shared), float(rho), float(p)))
    return results


def _run_partial(
    regulator: str,
    targets: list[str],
    expr: dict[str, float],
    meth_cache: dict[str, dict[str, float]],
    burden_index: dict[str, float],
) -> None:
    W = 88
    print("=" * W)
    print(f"  {'Target':<8}  {'n':>4}  {'Naive rho':>9}  {'Naive p':>12}  "
          f"{'Partial rho':>11}  {'Partial p':>12}  Result")
    print("=" * W)
    for gene in targets:
        meth = meth_cache.get(gene, {})
        shared = sorted(set(expr) & set(meth) & set(burden_index))
        if len(shared) < 10:
            print(f"  {gene:<8}  n={len(shared)} — skipped")
            continue
        x = np.array([expr[s] for s in shared])
        y = np.array([meth[s] for s in shared])
        z = np.array([burden_index[s] for s in shared])
        naive_r, naive_p = stats.spearmanr(x, y)
        part_r, part_p = _partial_spearman(x, y, z)
        if part_r < -0.2 and part_p < 0.05:
            verdict = "Holds after burden control"
        elif part_p >= 0.05:
            verdict = "Attenuated — possible global methylation confound"
        else:
            verdict = "|rho| < 0.2 after control"
        print(
            f"  {gene:<8}  {len(shared):>4}  {float(naive_r):>+9.3f}  {float(naive_p):>12.2e}  "
            f"{part_r:>+11.3f}  {part_p:>12.2e}  {verdict}"
        )
    print("=" * W)


def main() -> None:
    study_id = _TCGA_STUDY_IDS[TCGA_NETWORK]

    # ── profiles ──────────────────────────────────────────────────────────────
    rna_profile, meth_profile = _discover_profiles(study_id)
    if not rna_profile or not meth_profile:
        sys.exit("ERROR: could not discover molecular profiles")
    print(f"Study:               {study_id}")
    print(f"RNA profile:         {rna_profile}")
    print(f"Methylation profile: {meth_profile}\n")

    # ── squamous-only sample list ─────────────────────────────────────────────
    squamous_ids = set(
        _get_sample_ids_by_clinical_attr(
            study_id,
            attribute_id="CANCER_TYPE_DETAILED",
            allowed_values=[SQUAMOUS_VALUE],
        )
    )
    all_sample_ids = _get_sample_ids(study_id)
    if squamous_ids:
        sample_ids = [s for s in all_sample_ids if s in squamous_ids]
        print(f"Cohort (all histotypes): {len(all_sample_ids)} samples")
        print(f"Squamous-only filter:    {len(sample_ids)} samples\n")
    else:
        sample_ids = all_sample_ids
        print(f"WARNING: histotype filter unavailable — using all {len(sample_ids)} samples\n")

    # ── resolve all gene IDs in one batch ────────────────────────────────────
    all_genes = (
        ["DACH1", "SOX8", "SALL2"]
        + DACH1_TARGETS_ALL
        + SOX8_TARGETS   # DCC, FHIT, TFPI2 — same as SALL2_TARGETS
        + BURDEN_GENES
    )
    all_genes = list(dict.fromkeys(all_genes))  # deduplicate, preserve order
    print(f"Resolving Entrez IDs for {len(all_genes)} genes...")
    entrez_map = _entrez_ids(all_genes)
    missing = [g for g in all_genes if g not in entrez_map]
    if missing:
        print(f"  WARNING: unresolved — {missing}")
    print()

    # ── fetch expression for all three regulators ─────────────────────────────
    expr: dict[str, dict[str, float]] = {}
    for reg in ["DACH1", "SOX8", "SALL2"]:
        if reg not in entrez_map:
            print(f"  {reg}: Entrez ID not resolved — skipping")
            continue
        print(f"Fetching {reg} expression...")
        expr[reg] = _fetch_values_timed(rna_profile, entrez_map[reg], sample_ids)
        print(f"  {len(expr[reg])} values")
    print()

    # ── fetch methylation for all unique target genes ─────────────────────────
    target_genes = list(dict.fromkeys(DACH1_TARGETS_ALL + SOX8_TARGETS))
    meth_cache: dict[str, dict[str, float]] = {}
    print("Fetching target methylation...")
    for gene in target_genes:
        if gene not in entrez_map:
            print(f"  {gene}: Entrez ID not resolved — skipping")
            continue
        vals = _fetch_values_timed(meth_profile, entrez_map[gene], sample_ids)
        if vals:
            meth_cache[gene] = vals
            print(f"  {gene}: {len(vals)} values")
        else:
            print(f"  {gene}: no data returned")
    print()

    # ── fetch burden gene methylation and build index ─────────────────────────
    print("Fetching burden genes...")
    burden_meth: dict[str, dict[str, float]] = {}
    for gene in BURDEN_GENES:
        if gene not in entrez_map:
            print(f"  {gene}: Entrez ID not resolved — skipping")
            continue
        vals = _fetch_values_timed(meth_profile, entrez_map[gene], sample_ids)
        if vals:
            burden_meth[gene] = vals
            print(f"  {gene}: {len(vals)} values")
        else:
            print(f"  {gene}: no data (skipped from burden index)")

    dach1_expr = expr.get("DACH1", {})
    burden_samples = set(dach1_expr.keys())
    for vals in burden_meth.values():
        burden_samples &= set(vals.keys())
    burden_index: dict[str, float] = {
        s: float(np.mean([burden_meth[g][s] for g in burden_meth]))
        for s in burden_samples
    }
    print(f"\nBurden index: {len(burden_meth)} genes ({', '.join(burden_meth.keys())})")
    print(f"Samples with complete burden data: {len(burden_index)}\n")

    # ══════════════════════════════════════════════════════════════════════════
    # [TABLE 2]  DACH1 naive correlations — all 4 ARACNe targets
    # ══════════════════════════════════════════════════════════════════════════
    print("=" * 60)
    print("[TABLE 2]  DACH1 naive Spearman — all 4 targets (squamous-only)")
    print("=" * 60)
    _run_naive("DACH1", DACH1_TARGETS_ALL, dach1_expr, meth_cache)

    # ══════════════════════════════════════════════════════════════════════════
    # [TABLE 3]  DACH1 naive vs partial — 3 confirmed targets
    # ══════════════════════════════════════════════════════════════════════════
    print()
    print("[TABLE 3]  DACH1 naive vs partial Spearman — 3 confirmed targets (squamous-only)")
    _run_partial("DACH1", DACH1_TARGETS_PARTIAL, dach1_expr, meth_cache, burden_index)

    # ══════════════════════════════════════════════════════════════════════════
    # [SOX8]  Exploratory naive correlations
    # ══════════════════════════════════════════════════════════════════════════
    print()
    print("=" * 60)
    print("[SOX8]  Exploratory naive Spearman (squamous-only)")
    print("=" * 60)
    sox8_expr = expr.get("SOX8", {})
    if sox8_expr:
        _run_naive("SOX8", SOX8_TARGETS, sox8_expr, meth_cache)
    else:
        print("  SOX8 expression not available — skipped")

    # ══════════════════════════════════════════════════════════════════════════
    # [SALL2]  Exploratory naive correlations
    # ══════════════════════════════════════════════════════════════════════════
    print()
    print("=" * 60)
    print("[SALL2]  Exploratory naive Spearman (squamous-only)")
    print("=" * 60)
    sall2_expr = expr.get("SALL2", {})
    if sall2_expr:
        _run_naive("SALL2", SALL2_TARGETS, sall2_expr, meth_cache)
    else:
        print("  SALL2 expression not available — skipped")

    # ── CSV output (wide format, one row per sample) ──────────────────────────
    out_path = "dach1_partial_corr_cesc.csv"
    shared_all = sorted(
        set(dach1_expr)
        & set.union(*(set(meth_cache.get(g, {})) for g in DACH1_TARGETS_ALL))
    )
    if shared_all:
        fieldnames = ["sample_id", "DACH1_expr", "burden_index"] + [
            f"{g}_meth" for g in DACH1_TARGETS_ALL if g in meth_cache
        ]
        with open(out_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            for s in shared_all:
                row: dict = {
                    "sample_id": s,
                    "DACH1_expr": dach1_expr.get(s, ""),
                    "burden_index": burden_index.get(s, ""),
                }
                for g in DACH1_TARGETS_ALL:
                    row[f"{g}_meth"] = meth_cache.get(g, {}).get(s, "")
                writer.writerow(row)
        print(f"\nPer-sample data written to {out_path}  ({len(shared_all)} rows)")


if __name__ == "__main__":
    main()
