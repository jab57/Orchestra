"""
Threshold-by-group enrichment tables for the corroboration paper (Tables 1 and 2 of
private/manuscript/arxiv/corroboration/): is a candidate regulator with
corroboration_count >= t more likely to be a known cancer driver than one below t,
in focal genes and, separately, in pooled negative controls?

Per panel (BRCA/COAD, STAD) and group (focal, negative), for t in {1, 2, 3}:
one-tailed Fisher's exact test (alternative = "greater") of ground-truth membership,
corroborated vs uncorroborated. A comparison is untestable when no candidate reaches
the threshold (or none falls below it). Benjamini-Hochberg correction is applied across
the testable comparisons within each panel as one family. Focal rows are also reported
deduplicated to one row per unique regulator (negative controls are never deduplicated,
matching the paper's convention).

This is the statistics step the data-collection scripts leave to the paper; it exists so
the tables are reproducible from the candidate JSONs. It has no OncoKB option by design.

Inputs (produced by scripts/experiment_corroboration_ranked_brca_coad.py and
scripts/experiment_corroboration_ranked_stad.py):
    outputs/corroboration_ranked_brca_coad_raw.json
    outputs/corroboration_ranked_stad_raw.json

Output:
    outputs/threshold_tables_<truth>.json
"""
import argparse
import json
import sys
from pathlib import Path

from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from experiment_cgc_crosscheck import dedup, load_panel  # noqa: E402
from ground_truth import load_truth  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUTPUTS = ROOT / "outputs"
THRESHOLDS = (1, 2, 3)
GROUPS = ("focal", "negative")


def fisher_at_threshold(rows, category, truth_set, threshold):
    """One-tailed Fisher test at `threshold`; returns a dict with testable=False if degenerate."""
    sub = [r for r in rows if r["category"] == category]
    corrob = [r for r in sub if r["count"] >= threshold]
    uncorrob = [r for r in sub if r["count"] < threshold]
    out: dict = {"threshold": threshold, "n_corrob": len(corrob), "n_uncorrob": len(uncorrob)}
    if not corrob or not uncorrob:
        out.update({"testable": False, "truth_pos_corrob": None, "truth_pos_uncorrob": None,
                    "rate_corrob": None, "rate_uncorrob": None, "OR": None, "p": None})
        return out
    a = sum(1 for r in corrob if r["gene"] in truth_set)
    c = sum(1 for r in uncorrob if r["gene"] in truth_set)
    odds, p = fisher_exact([[a, len(corrob) - a], [c, len(uncorrob) - c]], alternative="greater")
    out.update({
        "testable": True, "truth_pos_corrob": a, "truth_pos_uncorrob": c,
        "rate_corrob": a / len(corrob), "rate_uncorrob": c / len(uncorrob),
        "OR": odds, "p": p,
    })
    return out


def analyze_panel(rows, truth_set):
    """Raw and focal-deduplicated tables for one panel, with BH across testable raw tests."""
    raw = {g: [fisher_at_threshold(rows, g, truth_set, t) for t in THRESHOLDS] for g in GROUPS}

    testable = [(g, i) for g in GROUPS for i, r in enumerate(raw[g]) if r["testable"]]
    if testable:
        adj = multipletests([raw[g][i]["p"] for g, i in testable], method="fdr_bh")[1]
        for (g, i), q in zip(testable, adj):
            raw[g][i]["p_adj_bh"] = float(q)
    for g in GROUPS:
        for r in raw[g]:
            r.setdefault("p_adj_bh", None)

    focal_dedup = [fisher_at_threshold(dedup(rows), "focal", truth_set, t) for t in THRESHOLDS]
    return {
        "raw": raw,
        "focal_dedup": focal_dedup,
        "n_bh_family": len(testable),
        "n_focal_rows": sum(1 for r in rows if r["category"] == "focal"),
        "n_focal_unique": len(dedup([r for r in rows if r["category"] == "focal"])),
        "n_negative_rows": sum(1 for r in rows if r["category"] == "negative"),
    }


def main():
    ap = argparse.ArgumentParser(description="Threshold-by-group enrichment tables.")
    ap.add_argument("--data-dir", type=Path, default=OUTPUTS,
                    help="folder holding corroboration_ranked_*_raw.json; results are written there too")
    ap.add_argument("--truth", choices=("intogen", "vogelstein"), required=True,
                    help="ground-truth gene set (no OncoKB option by design)")
    args = ap.parse_args()
    data_dir = args.data_dir
    truth_set = load_truth(args.truth)
    print(f"Ground truth: {args.truth} ({len(truth_set)} genes)")

    results = {"truth": args.truth, "n_truth_genes": len(truth_set)}
    for name, fname in [("brca_coad", "corroboration_ranked_brca_coad_raw.json"),
                        ("stad", "corroboration_ranked_stad_raw.json")]:
        rows = load_panel(data_dir / fname)
        res = analyze_panel(rows, truth_set)
        results[name] = res
        print(f"\n=== {name}: {res['n_focal_rows']} focal rows ({res['n_focal_unique']} unique), "
              f"{res['n_negative_rows']} negative-control rows; BH family = {res['n_bh_family']} ===")
        for g in GROUPS:
            for r in res["raw"][g]:
                if not r["testable"]:
                    print(f"  {g:8s} >={r['threshold']}: untestable (corrob n={r['n_corrob']}, "
                          f"uncorrob n={r['n_uncorrob']})")
                    continue
                print(f"  {g:8s} >={r['threshold']}: {r['truth_pos_corrob']}/{r['n_corrob']} vs "
                      f"{r['truth_pos_uncorrob']}/{r['n_uncorrob']}  OR={r['OR']:.2f}  "
                      f"p={r['p']:.4f}  p_adj={r['p_adj_bh']:.4f}")

    out_path = data_dir / f"threshold_tables_{args.truth}.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
