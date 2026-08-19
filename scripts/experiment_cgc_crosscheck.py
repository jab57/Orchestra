"""
Robustness check for the corroboration-paper (manuscript/Orchestra_research_paper.tex):
does the >=2 corroboration-vs-OncoKB pattern (Sections 3.1/3.6 of the paper) replicate
against a second, independently curated cancer-gene ground truth?

Ground truth used: the Sanger Institute's COSMIC Cancer Gene Census (Sondka et al. 2018,
Nature Reviews Cancer). COSMIC's own portal requires an account to download the census
directly; this script instead retrieves Sanger CGC membership via OncoKB's public
`cancerGeneList` API, which republishes the Sanger-curated `sangerCGC` flag for
convenience alongside its own annotations. The underlying gene-list curation is the
Sanger Institute's, independent of OncoKB's own MSK-based curation and of DoRothEA's
TF-focused curation (one of the four CASCADE evidence sources) -- but the API dependency
on OncoKB itself is a residual limitation, noted in the paper's Limitations section.

Inputs (already produced by scripts/experiment_corroboration_ranked_brca_coad.py and
scripts/experiment_corroboration_ranked_stad.py):
    outputs/corroboration_ranked_brca_coad_raw.json
    outputs/corroboration_ranked_stad_raw.json

Output:
    outputs/cgc_crosscheck_results.json
"""

import json
import random
import urllib.request
from pathlib import Path

from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parent.parent
OUTPUTS = ROOT / "outputs"

ONCOKB_LIST_URL = "https://www.oncokb.org/api/v1/utils/cancerGeneList"


def fetch_oncokb_gene_list():
    """Fetch OncoKB's public cancer gene list (includes the Sanger CGC cross-reference flag)."""
    with urllib.request.urlopen(ONCOKB_LIST_URL, timeout=30) as resp:
        return json.load(resp)


def load_panel(path):
    rows = []
    for q in json.load(open(path)):
        category = "focal" if q["category"] == "focal" else "negative"
        for c in q["scored_candidates"]:
            rows.append({"gene": c["gene"], "category": category, "count": c["corroboration_count"]})
    return rows


def dedup(rows):
    seen = {}
    for r in rows:
        seen.setdefault(r["gene"], r)
    return list(seen.values())


def fisher_ge2(rows, category, truth_set):
    sub = [r for r in rows if r["category"] == category]
    corrob = [r for r in sub if r["count"] >= 2]
    uncorrob = [r for r in sub if r["count"] < 2]
    a = sum(1 for r in corrob if r["gene"] in truth_set)
    b = len(corrob) - a
    c = sum(1 for r in uncorrob if r["gene"] in truth_set)
    d = len(uncorrob) - c
    odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return {
        "n_corrob": a + b, "oncokb_pos_corrob": a, "rate_corrob": a / (a + b) if (a + b) else float("nan"),
        "n_uncorrob": c + d, "oncokb_pos_uncorrob": c, "rate_uncorrob": c / (c + d) if (c + d) else float("nan"),
        "OR": odds, "p": p,
    }


def permutation_test(rows_dedup, category, truth_set, n=1000, seed=42):
    sub = [r for r in rows_dedup if r["category"] == category]
    labels = [r["count"] >= 2 for r in sub]
    truth = [r["gene"] in truth_set for r in sub]
    n_corrob = sum(labels)
    obs_rate_c = sum(t for l, t in zip(labels, truth) if l) / n_corrob
    obs_rate_u = sum(t for l, t in zip(labels, truth) if not l) / (len(sub) - n_corrob)
    obs_gap = obs_rate_c - obs_rate_u
    rng = random.Random(seed)
    idx = list(range(len(sub)))
    count_ge = 0
    for _ in range(n):
        rng.shuffle(idx)
        perm_corrob = set(idx[:n_corrob])
        rc = sum(truth[i] for i in perm_corrob) / n_corrob
        ru = sum(truth[i] for i in range(len(sub)) if i not in perm_corrob) / (len(sub) - n_corrob)
        if (rc - ru) >= obs_gap:
            count_ge += 1
    return {"observed_gap_pp": obs_gap * 100, "n_draws": n, "empirical_p": count_ge / n}


def main():
    oncokb_full = fetch_oncokb_gene_list()
    with open(OUTPUTS / "oncokb_cancer_gene_list.json", "w") as f:
        json.dump(oncokb_full, f, indent=2)

    sanger_cgc = sorted(set(x["hugoSymbol"] for x in oncokb_full if x.get("sangerCGC")))
    with open(OUTPUTS / "sanger_cgc_genes.json", "w") as f:
        json.dump(sanger_cgc, f, indent=2)
    sanger_cgc_set = set(sanger_cgc)
    print(f"Sanger CGC gene set: {len(sanger_cgc_set)} genes (via OncoKB cancerGeneList API)")

    brca_coad = load_panel(OUTPUTS / "corroboration_ranked_brca_coad_raw.json")
    stad = load_panel(OUTPUTS / "corroboration_ranked_stad_raw.json")

    results = {}
    for name, rows in [("brca_coad", brca_coad), ("stad", stad)]:
        panel_result = {}
        for category in ["focal", "negative"]:
            raw = fisher_ge2(rows, category, sanger_cgc_set)
            rows_d = dedup(rows) if category == "focal" else rows
            deduped = fisher_ge2(rows_d, category, sanger_cgc_set)
            panel_result[category] = {"raw": raw, "dedup": deduped}
            print(f"{name} {category}: raw OR={raw['OR']:.2f} p={raw['p']:.4f} | "
                  f"dedup OR={deduped['OR']:.2f} p={deduped['p']:.4f}")
        panel_result["permutation_focal_dedup"] = permutation_test(dedup(rows), "focal", sanger_cgc_set)
        print(f"{name} permutation (focal, dedup): {panel_result['permutation_focal_dedup']}")
        results[name] = panel_result

    with open(OUTPUTS / "cgc_crosscheck_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {OUTPUTS / 'cgc_crosscheck_results.json'}")


if __name__ == "__main__":
    main()
