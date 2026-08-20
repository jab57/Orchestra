"""
Does >=2-source corroboration add predictive power beyond MI edge weight alone --
the network-topology signal shown in Section "MI edge weight itself predicts
candidate quality" to predict OncoKB membership on its own (top-3-by-weight vs
bottom-3-by-weight within each query's capped 10, p=0.0003)?

Companion to scripts/experiment_dorothea_incremental_value.py, same two-check
design, with MI-top-3 (binary) standing in for DoRothEA as the "strongest single
source" being controlled for:

1. Stratified test: among candidates NOT in their query's top 3 by MI weight, does
   corroboration_count>=2 still lift OncoKB rate above baseline? Fisher's exact
   plus a 1,000-draw permutation test.

2. Logistic regression of OncoKB membership on mi_top3 (binary) and
   corroboration_count (0-4), with a likelihood-ratio test against an
   mi_top3-only reduced model.

MI rank is reconstructed from candidate order within each query's
scored_candidates list, not a stored weight value -- confirmed by code reading
(orchestra_langgraph_workflow.py: asyncio.gather preserves input order; dict()
of its results preserves insertion order; scripts/experiment_corroboration_*.py
iterate tumor_acquired_cascade_validation.keys() in that same order when writing
scored_candidates) that this order is exactly descending MI edge weight, the
same convention Section "MI edge weight itself predicts candidate quality" uses
(top-3/bottom-3 by position within the capped 10, not raw weight magnitude).

Inputs: outputs/corroboration_ranked_brca_coad_raw.json, ..._stad_raw.json,
outputs/oncokb_cancer_gene_list.json
Output: outputs/mi_weight_incremental_value_results.json
"""
import json
from pathlib import Path

import numpy as np
import statsmodels.api as sm
from scipy.stats import chi2, fisher_exact
import random

ROOT = Path(__file__).resolve().parent.parent
OUTPUTS = ROOT / "outputs"


def load_focal_dedup(path, oncokb):
    rows = []
    for q in json.load(open(path, encoding="utf-8")):
        if q["category"] != "focal":
            continue
        for i, c in enumerate(q["scored_candidates"]):
            rows.append({
                "gene": c["gene"],
                "mi_top3": int(i < 3),
                "count": c["corroboration_count"],
                "oncokb": int(c["gene"] in oncokb),
            })
    seen = {}
    for r in rows:
        seen.setdefault(r["gene"], r)
    return list(seen.values())


def stratified_test(rows, n_perm=1000, seed=42):
    strat = [r for r in rows if r["mi_top3"] == 0]
    corrob = [r for r in strat if r["count"] >= 2]
    uncorrob = [r for r in strat if r["count"] < 2]
    a = sum(r["oncokb"] for r in corrob)
    b = sum(r["oncokb"] for r in uncorrob)
    if not corrob or not uncorrob:
        return None
    odds, p = fisher_exact([[a, len(corrob) - a], [b, len(uncorrob) - b]], alternative="greater")

    labels = [r["count"] >= 2 for r in strat]
    truth = [r["oncokb"] for r in strat]
    n_hit = sum(labels)
    obs_gap = (sum(t for l, t in zip(labels, truth) if l) / n_hit) - \
              (sum(t for l, t in zip(labels, truth) if not l) / (len(strat) - n_hit))
    rng = random.Random(seed)
    idx = list(range(len(strat)))
    count_ge = 0
    for _ in range(n_perm):
        rng.shuffle(idx)
        hit_idx = set(idx[:n_hit])
        rc = sum(truth[i] for i in hit_idx) / n_hit
        ru = sum(truth[i] for i in range(len(strat)) if i not in hit_idx) / (len(strat) - n_hit)
        if (rc - ru) >= obs_gap:
            count_ge += 1
    return {
        "n_not_top3": len(strat), "n_corroborated": len(corrob), "n_uncorroborated": len(uncorrob),
        "rate_corroborated": a / len(corrob), "rate_uncorroborated": b / len(uncorrob),
        "OR": odds, "fisher_p": p,
        "observed_gap_pp": obs_gap * 100, "permutation_p": count_ge / n_perm,
    }


def logistic_lr_test(rows):
    y = np.array([r["oncokb"] for r in rows])
    mi_top3 = np.array([r["mi_top3"] for r in rows])
    count = np.array([r["count"] for r in rows])

    m_full = sm.Logit(y, sm.add_constant(np.column_stack([mi_top3, count]))).fit(disp=0)
    m_reduced = sm.Logit(y, sm.add_constant(mi_top3)).fit(disp=0)

    lr_stat = 2 * (m_full.llf - m_reduced.llf)
    lr_p = chi2.sf(lr_stat, df=1)

    return {
        "n": len(rows),
        "mi_top3_coef": float(m_full.params[1]), "mi_top3_p": float(m_full.pvalues[1]),
        "count_coef": float(m_full.params[2]), "count_p": float(m_full.pvalues[2]),
        "lr_stat": float(lr_stat), "lr_p": float(lr_p),
    }


def main():
    oncokb = set(x["hugoSymbol"] for x in json.load(open(OUTPUTS / "oncokb_cancer_gene_list.json", encoding="utf-8")))

    results = {}
    for name, path in [
        ("brca_coad", OUTPUTS / "corroboration_ranked_brca_coad_raw.json"),
        ("stad", OUTPUTS / "corroboration_ranked_stad_raw.json"),
    ]:
        rows = load_focal_dedup(path, oncokb)
        stratified = stratified_test(rows)
        logistic = logistic_lr_test(rows)
        results[name] = {"n_unique_focal": len(rows), "stratified": stratified, "logistic_lr": logistic}

        print(f"\n=== {name} (n={len(rows)} unique focal candidates) ===")
        if stratified:
            print(f"  Stratified (not MI-top-3, n={stratified['n_not_top3']}): "
                  f"corrob n={stratified['n_corroborated']} rate={stratified['rate_corroborated']:.1%}, "
                  f"uncorrob n={stratified['n_uncorroborated']} rate={stratified['rate_uncorroborated']:.1%}, "
                  f"OR={stratified['OR']:.2f}, Fisher p={stratified['fisher_p']:.4f}, "
                  f"permutation p={stratified['permutation_p']:.4f}")
        print(f"  Logistic: mi_top3 p={logistic['mi_top3_p']:.4f}, corroboration_count p={logistic['count_p']:.4f}, "
              f"LR test p={logistic['lr_p']:.4f}")

    with open(OUTPUTS / "mi_weight_incremental_value_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {OUTPUTS / 'mi_weight_incremental_value_results.json'}")


if __name__ == "__main__":
    main()
