r"""Items 1 (expression part) and 2: purity confounding + target-expression adjustment.

Inputs (all read-only):
  dach1_partial_corr_cesc.csv      original per-sample table (253 squamous CESC)
  target_expression_cesc.csv       from fetch_target_expression.py
  ABSOLUTE table (argv[1])         PanCanAtlas ABSOLUTE calls, GDC 4f277128-f793-4354-a13d-30cc7fe9f6b5
Outputs (new files in this folder): purity_absolute_cesc.csv, item2_purity_results.csv, item1_expression_results.csv

Partial Spearman = residuals-on-ranks, same method as dach1_partial_correlation.py,
generalised to several covariates (OLS on ranks).
Run from c:\Dev\Orchestra:  python -X utf8 <this file> <absolute.txt>
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

HERE = Path(__file__).parent
TARGETS = ["CADM1", "ESR1", "SLIT2", "EDNRB"]

d = pd.read_csv("dach1_partial_corr_cesc.csv").set_index("sample_id")
e = pd.read_csv(HERE / "target_expression_cesc.csv").set_index("sample_id")
a = pd.read_csv(sys.argv[1], sep="\t").set_index("array")
a = a[a["call status"] == "called"]
pur = a.loc[a.index.intersection(d.index), ["purity", "ploidy"]]
pur.to_csv(HERE / "purity_absolute_cesc.csv")
d = d.join(pur).join(e.drop(columns="DACH1").add_suffix("_expr"))
print(f"samples {len(d)}; with ABSOLUTE purity {d.purity.notna().sum()}")
print("purity: median %.2f, IQR %.2f-%.2f, range %.2f-%.2f" % (
    d.purity.median(), *d.purity.quantile([.25, .75]), d.purity.min(), d.purity.max()))

def partial(x, y, Z):
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    if Z is None or len(Z) == 0: return stats.pearsonr(rx, ry)
    A = np.column_stack([np.ones(len(rx))] + [stats.rankdata(z) for z in Z])
    res = lambda v: v - A @ np.linalg.lstsq(A, v, rcond=None)[0]
    r, p = stats.pearsonr(res(rx), res(ry))
    return r, p   # note: p uses n-2 df (as in original); with k covariates, n-2-k is more exact

def fisher_ci(r, n, k=0):
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3 - k)
    return np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)

rows = []
sub = d.dropna(subset=["purity"])
for t in TARGETS:
    m = f"{t}_meth"
    for label, data, covs in [
        ("naive, all 253", d, []),
        ("naive, purity subset", sub, []),
        ("partial | purity", sub, ["purity"]),
        ("partial | burden", sub, ["burden_index"]),
        ("partial | purity + burden", sub, ["purity", "burden_index"]),
    ]:
        x, y = data.DACH1_expr.values, data[m].values
        r, p = partial(x, y, [data[c].values for c in covs])
        lo, hi = fisher_ci(r, len(data), len(covs))
        rows.append(dict(target=t, analysis=label, n=len(data), rho=round(r, 3), ci_lo=round(lo, 3), ci_hi=round(hi, 3), p=f"{p:.2e}"))
res = pd.DataFrame(rows); res.to_csv(HERE / "item2_purity_results.csv", index=False)
print("\n== DACH1 expression vs target methylation, purity adjustment (Fisher CI approximate)")
print(res.pivot(index="analysis", columns="target", values="rho").loc[
    ["naive, all 253", "naive, purity subset", "partial | burden", "partial | purity", "partial | purity + burden"]].to_string())
print(res[res.analysis.str.contains("purity")].to_string(index=False))

print("\n== Purity vs expression / methylation (Spearman, n=%d)" % len(sub))
for c in ["DACH1_expr", "ESR1_expr", "SLIT2_expr", "CADM1_expr", "EDNRB_expr", "burden_index"] + [f"{t}_meth" for t in TARGETS]:
    if c in sub:
        r, p = stats.spearmanr(sub.purity, sub[c]); print(f"  purity vs {c:<13} rho={r:+.3f} p={p:.2e}")

# ---- item 1 (expression part)
print("\n== Item 1: DACH1 vs target EXPRESSION, and DACH1-meth adjusted for target's own expression (all 253)")
r1 = []
for t in TARGETS:
    x, ex, y = d.DACH1_expr.values, d[f"{t}_expr"].values, d[f"{t}_meth"].values
    rho_e, p_e = stats.spearmanr(x, ex)
    rho_m, p_m = stats.spearmanr(x, y)
    rho_a, p_a = partial(x, y, [ex])
    rho_ab, p_ab = partial(x, y, [ex, d.burden_index.values])
    rho_me, p_me = stats.spearmanr(y, ex)
    r1.append(dict(target=t, DACH1_vs_target_expr_rho=round(rho_e, 3), expr_p=f"{p_e:.1e}",
                   naive_meth_rho=round(rho_m, 3), meth_vs_own_expr_rho=round(rho_me, 3),
                   partial_given_own_expr=round(rho_a, 3), partial_p=f"{p_a:.1e}",
                   partial_given_expr_and_burden=round(rho_ab, 3)))
r1 = pd.DataFrame(r1); r1.to_csv(HERE / "item1_expression_results.csv", index=False)
print(r1.to_string(index=False))
