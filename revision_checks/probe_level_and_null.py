r"""Item 1 (expression-independent promoter definitions) and item 3 (genome-wide null for DACH1).

Inputs (this folder): gene_promoter_meth_cesc_{A,B}*.csv.gz, target_expression_cesc.csv, purity_absolute_cesc.csv
                      ../../../../dach1_partial_corr_cesc.csv  (read from cwd = c:\Dev\Orchestra)
Outputs: item1_probe_level_results.csv, item3_null_summary.csv, item3_null_all_genes_A.csv
Run from c:\Dev\Orchestra:  python -I -X utf8 <this file>
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

HERE = Path(__file__).parent
d = pd.read_csv("dach1_partial_corr_cesc.csv").set_index("sample_id")
e = pd.read_csv(HERE / "target_expression_cesc.csv").set_index("sample_id")
pur = pd.read_csv(HERE / "purity_absolute_cesc.csv").set_index("array" if "array" in pd.read_csv(HERE / "purity_absolute_cesc.csv", nrows=1) else pd.read_csv(HERE / "purity_absolute_cesc.csv", nrows=1).columns[0])
d = d.join(pur[["purity"]])
DEF = {"A": "gene_promoter_meth_cesc_A_tss200_fallback1500.csv.gz", "B": "gene_promoter_meth_cesc_B_tss1500_tss200_exon1.csv.gz"}
M = {k: pd.read_csv(HERE / v, index_col=0)[d.index] for k, v in DEF.items()}   # genes x samples, same order
TARGETS = ["CADM1", "ESR1", "SLIT2", "EDNRB"]

def resid(v, Z):
    A = np.column_stack([np.ones(len(v))] + Z); return v - A @ np.linalg.lstsq(A, v, rcond=None)[0]
def pcorr(x, y, covs=()):
    rx, ry = stats.rankdata(x), stats.rankdata(y); Z = [stats.rankdata(c) for c in covs]
    return stats.pearsonr(resid(rx, Z), resid(ry, Z)) if Z else stats.pearsonr(rx, ry)
def ci(r, n, k=0):
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3 - k); return np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)

# ---------------- item 1
x = d.DACH1_expr.values; rows = []
for t in TARGETS:
    cb = d[f"{t}_meth"]                                    # cBioPortal (anticorrelation-selected probe)
    for k, m in M.items():
        y = m.loc[t].values
        ok = ~np.isnan(y); xs, ys = x[ok], y[ok]
        r, p = pcorr(xs, ys); lo, hi = ci(r, ok.sum())
        pp = d.purity.values[ok]; okp = ~np.isnan(pp)
        rp, pp_p = pcorr(xs[okp], ys[okp], [pp[okp]])
        ex = e[t].values[ok]
        re_, pe = pcorr(xs, ys, [ex])
        rows.append(dict(target=t, definition=k, n_probes=None, n=int(ok.sum()), rho=round(r, 3), ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                         p=f"{p:.1e}", rho_purity_adj=round(rp, 3), rho_own_expr_adj=round(re_, 3), p_own_expr_adj=f"{pe:.1e}",
                         meth_vs_own_expr=round(stats.spearmanr(ys, ex)[0], 3),
                         agreement_with_cbioportal_value=round(stats.spearmanr(ys, cb.values[ok])[0], 3)))
    rows.append(dict(target=t, definition="cBioPortal", n=len(d), rho=round(stats.spearmanr(x, cb)[0], 3)))
r1 = pd.DataFrame(rows); r1.to_csv(HERE / "item1_probe_level_results.csv", index=False)
print("== ITEM 1: DACH1 expr vs promoter methylation, three summaries\n", r1.drop(columns=["n_probes"]).to_string(index=False))

# ---------------- item 3: genome-wide null
def null_rhos(m, covs=None, sub=None):
    mm = m if sub is None else m.loc[:, sub]
    keep = mm.notna().all(axis=1) & (mm.std(axis=1) > 0); mm = mm[keep]
    xr = stats.rankdata(x if sub is None else x[sub])
    R = mm.rank(axis=1).values
    if covs is not None:
        Z = [stats.rankdata(c) for c in covs]
        A = np.column_stack([np.ones(R.shape[1])] + Z)
        B = np.linalg.lstsq(A, np.column_stack([xr, R.T]), rcond=None)[0]
        res = np.column_stack([xr, R.T]) - A @ B; xr2, R2 = res[:, 0], res[:, 1:].T
    else: xr2, R2 = xr - xr.mean(), R - R.mean(axis=1, keepdims=True)
    xr2 = xr2 - xr2.mean(); R2 = R2 - R2.mean(axis=1, keepdims=True)
    rho = (R2 @ xr2) / (np.sqrt((R2 ** 2).sum(axis=1)) * np.sqrt((xr2 ** 2).sum()))
    return pd.Series(rho, index=mm.index)

okp = d.purity.notna().values
summ = []
for k, m in M.items():
    for label, kw in (("unadjusted (n=253)", {}), ("purity-adjusted (n=%d)" % okp.sum(), dict(covs=[d.purity.values[okp]], sub=okp))):
        nr = null_rhos(m, **kw)
        if k == "A" and not kw: nr.rename("rho").to_csv(HERE / "item3_null_all_genes_A.csv")
        for t in TARGETS + ["DCC", "FHIT", "TFPI2"]:
            if t not in nr: continue
            v = nr[t]
            summ.append(dict(definition=k, analysis=label, gene=t, rho=round(v, 3), genes_in_null=len(nr),
                             null_median=round(nr.median(), 3), null_sd=round(nr.std(), 3),
                             pct_more_negative=round(100 * (nr <= v).mean(), 2), emp_p_one_sided=round((nr <= v).mean(), 4),
                             z_vs_null=round((v - nr.median()) / nr.std(), 2),
                             frac_null_abs_ge=round((nr.abs() >= abs(v)).mean(), 3)))
s3 = pd.DataFrame(summ); s3.to_csv(HERE / "item3_null_summary.csv", index=False)
print("\n== ITEM 3: DACH1 expression vs promoter methylation of ALL genes (null)\n", s3.to_string(index=False))
