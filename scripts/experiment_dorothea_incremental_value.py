"""
Does >=2-source corroboration add predictive power beyond DoRothEA alone -- the
single strongest CASCADE evidence source (Section "Robustness: does corroboration
add value beyond the strongest single source?" of the corroboration paper)?

Two checks per panel, on the deduplicated focal candidate set (one row per unique
gene, matching the paper's existing methodology throughout):

1. Stratified test: among candidates where DoRothEA is negative, does corroboration
   from the remaining three sources (>=2 of LINCS, super-enhancer, DepMap) still
   lift OncoKB rate above the DoRothEA-negative baseline? Fisher's exact test plus
   a 1,000-draw permutation test, matching the paper's other permutation checks.

2. Logistic regression of OncoKB membership on DoRothEA (binary) and the count of
   the other three sources (0-3), with a likelihood-ratio test comparing this model
   against a DoRothEA-only reduced model -- the more powerful of the two tests,
   since it uses the full candidate set (not just the DoRothEA-negative subset) and
   treats the other-source count continuously rather than at a single threshold.

Inputs (already produced by scripts/experiment_corroboration_ranked_brca_coad.py
and scripts/experiment_corroboration_ranked_stad.py):
    outputs/corroboration_ranked_brca_coad_raw.json
    outputs/corroboration_ranked_stad_raw.json
    outputs/oncokb_cancer_gene_list.json

Output:
    outputs/dorothea_incremental_value_results.json
"""
import json
import random
from pathlib import Path

import numpy as np
import statsmodels.api as sm
from scipy.stats import chi2, fisher_exact

ROOT = Path(__file__).resolve().parent.parent
OUTPUTS = ROOT / "outputs"


def load_focal_dedup(path, oncokb):
    rows = []
    for q in json.load(open(path, encoding="utf-8")):
        if q["category"] != "focal":
            continue
        for c in q["scored_candidates"]:
            rows.append({
                "gene": c["gene"],
                "dorothea": int(c["dorothea_tf"]),
                "other3": int(c["lincs_knockdown"]) + int(c["super_enhancer"]) + int(c["depmap_essential"]),
                "oncokb": int(c["gene"] in oncokb),
            })
    seen = {}
    for r in rows:
        seen.setdefault(r["gene"], r)
    return list(seen.values())


def stratified_test(rows, n_perm=1000, seed=42):
    strat = [r for r in rows if r["dorothea"] == 0]
    corrob = [r for r in strat if r["other3"] >= 2]
    uncorrob = [r for r in strat if r["other3"] < 2]
    a = sum(r["oncokb"] for r in corrob)
    b = sum(r["oncokb"] for r in uncorrob)
    if not corrob or not uncorrob:
        return None
    odds, p = fisher_exact([[a, len(corrob) - a], [b, len(uncorrob) - b]], alternative="greater")

    labels = [r["other3"] >= 2 for r in strat]
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
        "n_dorothea_negative": len(strat), "n_corroborated": len(corrob), "n_uncorroborated": len(uncorrob),
        "rate_corroborated": a / len(corrob), "rate_uncorroborated": b / len(uncorrob),
        "OR": odds, "fisher_p": p,
        "observed_gap_pp": obs_gap * 100, "permutation_p": count_ge / n_perm,
    }


def logistic_lr_test(rows):
    y = np.array([r["oncokb"] for r in rows])
    dorothea = np.array([r["dorothea"] for r in rows])
    other3 = np.array([r["other3"] for r in rows])

    m_full = sm.Logit(y, sm.add_constant(np.column_stack([dorothea, other3]))).fit(disp=0)
    m_reduced = sm.Logit(y, sm.add_constant(dorothea)).fit(disp=0)

    lr_stat = 2 * (m_full.llf - m_reduced.llf)
    lr_p = chi2.sf(lr_stat, df=1)

    return {
        "n": len(rows),
        "dorothea_coef": float(m_full.params[1]), "dorothea_p": float(m_full.pvalues[1]),
        "other3_coef": float(m_full.params[2]), "other3_p": float(m_full.pvalues[2]),
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
            print(f"  Stratified (DoRothEA-negative, n={stratified['n_dorothea_negative']}): "
                  f"corrob n={stratified['n_corroborated']} rate={stratified['rate_corroborated']:.1%}, "
                  f"uncorrob n={stratified['n_uncorroborated']} rate={stratified['rate_uncorroborated']:.1%}, "
                  f"OR={stratified['OR']:.2f}, Fisher p={stratified['fisher_p']:.4f}, "
                  f"permutation p={stratified['permutation_p']:.4f}")
        print(f"  Logistic: DoRothEA p={logistic['dorothea_p']:.4f}, other3_count p={logistic['other3_p']:.4f}, "
              f"LR test p={logistic['lr_p']:.4f}")

    with open(OUTPUTS / "dorothea_incremental_value_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {OUTPUTS / 'dorothea_incremental_value_results.json'}")


if __name__ == "__main__":
    main()
