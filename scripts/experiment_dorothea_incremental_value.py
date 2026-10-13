"""
Does >=2-source corroboration add predictive power beyond DoRothEA alone -- the
single strongest CASCADE evidence source (Section "Robustness: does corroboration
add value beyond the strongest single source?" of the corroboration paper)?

Two checks per panel, on the deduplicated focal candidate set (one row per unique
gene, matching the paper's existing methodology throughout):

1. Stratified test: among candidates where DoRothEA is negative, does corroboration
   from the remaining three sources (>=2 of LINCS, super-enhancer, DepMap) still
   lift the known-driver rate above the DoRothEA-negative baseline? Fisher's exact test plus
   a 1,000-draw permutation test, matching the paper's other permutation checks.

2. Logistic regression of known-driver membership on DoRothEA (binary) and the count of
   the other three sources (0-3), with a likelihood-ratio test comparing this model
   against a DoRothEA-only reduced model -- the more powerful of the two tests,
   since it uses the full candidate set (not just the DoRothEA-negative subset) and
   treats the other-source count continuously rather than at a single threshold.

Also emits the DoRothEA-alone-vs-known-driver Fisher test quoted in the opening
paragraph of that section (DoRothEA as a lone predictor, on the deduplicated
focal set and on the non-deduplicated negative controls), so every number in
the section is script-reproducible.

Inputs (already produced by scripts/experiment_corroboration_ranked_brca_coad.py
and scripts/experiment_corroboration_ranked_stad.py):
    <data-dir>/corroboration_ranked_brca_coad_raw.json
    <data-dir>/corroboration_ranked_stad_raw.json
(ground truth chosen with --truth; see ground_truth.py)

Output:
    <data-dir>/dorothea_incremental_value_results<suffix>.json
"""
import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np
import statsmodels.api as sm
from scipy.stats import chi2, fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ground_truth import CHOICES, load_truth, output_suffix  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUTPUTS = ROOT / "outputs"


def load_focal_dedup(path, truth):
    rows = []
    for q in json.load(open(path, encoding="utf-8")):
        if q["category"] != "focal":
            continue
        for c in q["scored_candidates"]:
            rows.append({
                "gene": c["gene"],
                "dorothea": int(c["dorothea_tf"]),
                "other3": int(c["lincs_knockdown"]) + int(c["super_enhancer"]) + int(c["depmap_essential"]),
                "truth": int(c["gene"] in truth),
            })
    seen = {}
    for r in rows:
        seen.setdefault(r["gene"], r)
    return list(seen.values())


def load_controls(path, truth):
    """Negative-control candidates, NOT deduplicated -- matching the paper's
    negative-control convention (Section "Statistical framework"): keeping every
    row gives the control test more power, so a null result there is not an
    artifact of a smaller sample."""
    rows = []
    for q in json.load(open(path, encoding="utf-8")):
        if q["category"] not in ("housekeeping", "non_driver"):
            continue
        for c in q["scored_candidates"]:
            rows.append({
                "gene": c["gene"],
                "dorothea": int(c["dorothea_tf"]),
                "truth": int(c["gene"] in truth),
            })
    return rows


def dorothea_alone_test(rows):
    """DoRothEA as a lone predictor of known-driver membership: Fisher's exact,
    one-tailed 'greater', DoRothEA-positive vs DoRothEA-negative candidates.
    Quoted in the opening paragraph of the incremental-value section."""
    pos = [r for r in rows if r["dorothea"] == 1]
    neg = [r for r in rows if r["dorothea"] == 0]
    if not pos or not neg:
        return None
    a = sum(r["truth"] for r in pos)
    b = sum(r["truth"] for r in neg)
    odds, p = fisher_exact([[a, len(pos) - a], [b, len(neg) - b]], alternative="greater")
    return {
        "n_dorothea_positive": len(pos), "n_dorothea_negative": len(neg),
        "rate_dorothea_positive": a / len(pos), "rate_dorothea_negative": b / len(neg),
        "OR": odds, "fisher_p": p,
    }


SOURCES = ("lincs_knockdown", "super_enhancer", "dorothea_tf", "depmap_essential")


def per_source_test(path, truth):
    """Each CASCADE evidence source as a lone predictor of ground-truth membership, on the
    deduplicated focal set: Fisher's exact, one-tailed 'greater', source-positive vs
    source-negative candidates. Descriptive: shows which source is strongest in a panel."""
    seen = {}
    for q in json.load(open(path, encoding="utf-8")):
        if q["category"] != "focal":
            continue
        for c in q["scored_candidates"]:
            seen.setdefault(c["gene"], c)
    out = {}
    for src in SOURCES:
        pos = [c for c in seen.values() if c[src]]
        neg = [c for c in seen.values() if not c[src]]
        if not pos or not neg:
            out[src] = None
            continue
        a = sum(c["gene"] in truth for c in pos)
        b = sum(c["gene"] in truth for c in neg)
        odds, p = fisher_exact([[a, len(pos) - a], [b, len(neg) - b]], alternative="greater")
        out[src] = {"n_positive": len(pos), "n_negative": len(neg),
                    "truth_pos_positive": a, "truth_pos_negative": b,
                    "rate_positive": a / len(pos), "rate_negative": b / len(neg),
                    "OR": odds, "fisher_p": p}
    return out


def stratified_test(rows, n_perm=1000, seed=42):
    strat = [r for r in rows if r["dorothea"] == 0]
    corrob = [r for r in strat if r["other3"] >= 2]
    uncorrob = [r for r in strat if r["other3"] < 2]
    a = sum(r["truth"] for r in corrob)
    b = sum(r["truth"] for r in uncorrob)
    if not corrob or not uncorrob:
        return None
    odds, p = fisher_exact([[a, len(corrob) - a], [b, len(uncorrob) - b]], alternative="greater")

    labels = [r["other3"] >= 2 for r in strat]
    truth = [r["truth"] for r in strat]
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
    y = np.array([r["truth"] for r in rows])
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
    ap = argparse.ArgumentParser(description="Incremental-value test against a chosen ground truth.")
    ap.add_argument("--data-dir", type=Path, default=OUTPUTS,
                    help="folder holding corroboration_ranked_*_raw.json; results are written there too")
    ap.add_argument("--truth", choices=CHOICES, required=True,
                    help="ground-truth gene set. 'oncokb'/'sanger_cgc' are LEGACY (reproduces v1-v3 only).")
    args = ap.parse_args()
    data_dir = args.data_dir
    truth_set = load_truth(args.truth)
    suffix = output_suffix(args.truth)
    print(f"Ground truth: {args.truth} ({len(truth_set)} genes)")

    results = {}
    for name, path in [
        ("brca_coad", data_dir / "corroboration_ranked_brca_coad_raw.json"),
        ("stad", data_dir / "corroboration_ranked_stad_raw.json"),
    ]:
        rows = load_focal_dedup(path, truth_set)
        controls = load_controls(path, truth_set)
        stratified = stratified_test(rows)
        logistic = logistic_lr_test(rows)
        dorothea_alone_focal = dorothea_alone_test(rows)
        dorothea_alone_controls = dorothea_alone_test(controls)
        results[name] = {
            "n_unique_focal": len(rows),
            "stratified": stratified,
            "logistic_lr": logistic,
            "dorothea_alone_focal_dedup": dorothea_alone_focal,
            "dorothea_alone_controls_raw": dorothea_alone_controls,
            "per_source_focal_dedup": per_source_test(path, truth_set),
        }

        print(f"\n=== {name} (n={len(rows)} unique focal candidates) ===")
        if dorothea_alone_focal:
            print(f"  DoRothEA alone (focal, dedup):    "
                  f"DoRothEA+ {dorothea_alone_focal['rate_dorothea_positive']:.1%} vs "
                  f"DoRothEA- {dorothea_alone_focal['rate_dorothea_negative']:.1%}, "
                  f"OR={dorothea_alone_focal['OR']:.2f}, Fisher p={dorothea_alone_focal['fisher_p']:.4f}")
        if dorothea_alone_controls:
            print(f"  DoRothEA alone (controls, raw):   "
                  f"DoRothEA+ {dorothea_alone_controls['rate_dorothea_positive']:.1%} vs "
                  f"DoRothEA- {dorothea_alone_controls['rate_dorothea_negative']:.1%}, "
                  f"OR={dorothea_alone_controls['OR']:.2f}, Fisher p={dorothea_alone_controls['fisher_p']:.4f}")
        if stratified:
            print(f"  Stratified (DoRothEA-negative, n={stratified['n_dorothea_negative']}): "
                  f"corrob n={stratified['n_corroborated']} rate={stratified['rate_corroborated']:.1%}, "
                  f"uncorrob n={stratified['n_uncorroborated']} rate={stratified['rate_uncorroborated']:.1%}, "
                  f"OR={stratified['OR']:.2f}, Fisher p={stratified['fisher_p']:.4f}, "
                  f"permutation p={stratified['permutation_p']:.4f}")
        print(f"  Logistic: DoRothEA p={logistic['dorothea_p']:.4f}, other3_count p={logistic['other3_p']:.4f}, "
              f"LR test p={logistic['lr_p']:.4f}")

    out_path = data_dir / f"dorothea_incremental_value_results{suffix}.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
