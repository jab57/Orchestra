"""
MI-edge-weight ranking check for the corroboration paper (Section "MI edge weight itself
predicts candidate quality, and shows no detected association with corroboration count").

Pools focal candidate-rows from both panels (not deduplicated). Within each query's capped
candidate list, position equals MI-weight rank (1 = highest weight; see
experiment_mi_weight_incremental_value.py for why order encodes rank).

Part A -- does MI weight predict ground-truth membership?
    top 3 (rank <= 3) vs bottom 3 (rank >= 8, i.e. the bottom 3 of a full capped 10; a query
    with fewer than 10 candidates contributes fewer bottom rows): one-tailed Fisher's exact
    ("greater"), and the Spearman correlation between rank and membership across all rows
    (negative = higher weight, more likely a driver).

Part B -- is MI weight associated with corroboration count? (needs no ground truth)
    Spearman correlation of rank vs corroboration_count; mean count for top 3 vs bottom 3;
    one-tailed Fisher's exact on the >=2 rate, top 3 vs bottom 3.

The Part B definitions were validated against the published v3 values (339 rows, rho=-0.049,
p=0.3687, 102 top / 101 bottom rows, means 0.96 vs 0.87, Fisher p=0.37) using no ground truth.

Inputs:  <data-dir>/corroboration_ranked_brca_coad_raw.json, ..._stad_raw.json
Output:  <data-dir>/mi_rank_check_<truth>.json
"""
import argparse
import json
import sys
from pathlib import Path

from scipy.stats import fisher_exact, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ground_truth import load_truth  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUTPUTS = ROOT / "outputs"
CAP = 10
TOP_N = 3


def load_focal_rows(paths):
    """Focal candidate-rows pooled across panels: dicts with gene, rank (1-based), count."""
    rows = []
    for path in paths:
        for q in json.load(open(path, encoding="utf-8")):
            if q["category"] != "focal" or "scored_candidates" not in q:
                continue
            for i, c in enumerate(q["scored_candidates"]):
                if "corroboration_count" not in c:
                    continue
                rows.append({"gene": c["gene"], "rank": i + 1, "count": c["corroboration_count"]})
    return rows


def split_top_bottom(rows, cap=CAP, top_n=TOP_N):
    top = [r for r in rows if r["rank"] <= top_n]
    bottom = [r for r in rows if r["rank"] > cap - top_n]
    return top, bottom


def _fisher_greater(pos_a, n_a, pos_b, n_b):
    odds, p = fisher_exact([[pos_a, n_a - pos_a], [pos_b, n_b - pos_b]], alternative="greater")
    return odds, p


def corroboration_association(rows):
    """Part B: needs no ground truth."""
    top, bottom = split_top_bottom(rows)
    sp = spearmanr([r["rank"] for r in rows], [r["count"] for r in rows])
    rho, p = float(sp.statistic), float(sp.pvalue)
    a = sum(r["count"] >= 2 for r in top)
    b = sum(r["count"] >= 2 for r in bottom)
    _, fp = _fisher_greater(a, len(top), b, len(bottom))
    return {
        "n_rows": len(rows), "n_top": len(top), "n_bottom": len(bottom),
        "spearman_rho": rho, "spearman_p": p,
        "mean_count_top": sum(r["count"] for r in top) / len(top),
        "mean_count_bottom": sum(r["count"] for r in bottom) / len(bottom),
        "ge2_top": a, "ge2_bottom": b, "fisher_p_ge2_greater": float(fp),
    }


def truth_association(rows, truth_set):
    """Part A: MI rank vs ground-truth membership."""
    for r in rows:
        r["truth"] = int(r["gene"] in truth_set)
    top, bottom = split_top_bottom(rows)
    a = sum(r["truth"] for r in top)
    b = sum(r["truth"] for r in bottom)
    odds, p = _fisher_greater(a, len(top), b, len(bottom))
    sp = spearmanr([r["rank"] for r in rows], [r["truth"] for r in rows])
    rho, rp = float(sp.statistic), float(sp.pvalue)
    return {
        "n_rows": len(rows), "n_top": len(top), "n_bottom": len(bottom),
        "truth_top": a, "truth_bottom": b,
        "rate_top": a / len(top), "rate_bottom": b / len(bottom),
        "OR": odds, "fisher_p": float(p),
        "spearman_rho": rho, "spearman_p": rp,
    }


def main():
    ap = argparse.ArgumentParser(description="MI-weight ranking check against a chosen ground truth.")
    ap.add_argument("--data-dir", type=Path, default=OUTPUTS,
                    help="folder holding corroboration_ranked_*_raw.json; results are written there too")
    ap.add_argument("--truth", choices=("intogen", "vogelstein"), required=True,
                    help="ground-truth gene set (no OncoKB option by design)")
    args = ap.parse_args()
    truth_set = load_truth(args.truth)
    print(f"Ground truth: {args.truth} ({len(truth_set)} genes)")

    panel_rows = {
        "brca_coad": load_focal_rows([args.data_dir / "corroboration_ranked_brca_coad_raw.json"]),
        "stad": load_focal_rows([args.data_dir / "corroboration_ranked_stad_raw.json"]),
    }
    rows = panel_rows["brca_coad"] + panel_rows["stad"]
    res = {"truth": args.truth,
           "truth_association": truth_association(rows, truth_set),
           "corroboration_association": corroboration_association(rows),
           # The pooled figures can hide a panel in which rank carries no signal.
           "per_panel": {name: {"truth_association": truth_association(r, truth_set),
                                "corroboration_association": corroboration_association(r)}
                         for name, r in panel_rows.items()}}

    t, c = res["truth_association"], res["corroboration_association"]
    print(f"\nPart A (rank vs {args.truth}): top3 {t['truth_top']}/{t['n_top']} ({t['rate_top']:.1%}) "
          f"vs bottom3 {t['truth_bottom']}/{t['n_bottom']} ({t['rate_bottom']:.1%}); "
          f"OR={t['OR']:.2f}, Fisher p={t['fisher_p']:.4f}; "
          f"Spearman rho={t['spearman_rho']:.3f} (p={t['spearman_p']:.4f}), n={t['n_rows']}")
    print(f"Part B (rank vs corroboration count): rho={c['spearman_rho']:.3f} (p={c['spearman_p']:.4f}); "
          f"mean count top3 {c['mean_count_top']:.2f} vs bottom3 {c['mean_count_bottom']:.2f}; "
          f">=2 rate Fisher p={c['fisher_p_ge2_greater']:.4f}")

    for name, pr in res["per_panel"].items():
        t = pr["truth_association"]
        print(f"  {name}: top3 {t['truth_top']}/{t['n_top']} ({t['rate_top']:.1%}) vs bottom3 "
              f"{t['truth_bottom']}/{t['n_bottom']} ({t['rate_bottom']:.1%}); Fisher p={t['fisher_p']:.4f}; "
              f"Spearman rho={t['spearman_rho']:.3f} (p={t['spearman_p']:.4f}), n={t['n_rows']}")

    out_path = args.data_dir / f"mi_rank_check_{args.truth}.json"
    with open(out_path, "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
