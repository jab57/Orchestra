r"""Item 4 (CIs + SOX8 control), item 5 (burden index with PAX1), item 6 (DACH1 low vs high).
Methylation sources compared: 'cBioPortal' (paper's values; DCC/PAX1 unavailable), 'A' and 'B' (probe-level, expression-independent).
Run from c:\Dev\Orchestra:  python -I -X utf8 <this file>
Outputs: item4_ci_table.csv, item4_sox8_difference_tests.csv, item5_burden_pax1.csv, item6_dach1_low_vs_high.csv"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

HERE = Path(__file__).parent; rng = np.random.default_rng(20261017); NB = 5000
d = pd.read_csv("dach1_partial_corr_cesc.csv").set_index("sample_id")
e = pd.read_csv(HERE / "target_expression_cesc.csv").set_index("sample_id")
x_ = pd.read_csv(HERE / "cbioportal_extra_meth_cesc.csv").set_index("sample_id")
pur = pd.read_csv(HERE / "purity_absolute_cesc.csv", index_col=0)
A = pd.read_csv(HERE / "gene_promoter_meth_cesc_A_tss200_fallback1500.csv.gz", index_col=0)[d.index].T
B = pd.read_csv(HERE / "gene_promoter_meth_cesc_B_tss1500_tss200_exon1.csv.gz", index_col=0)[d.index].T
cb = pd.concat([d[["CADM1_meth", "ESR1_meth", "SLIT2_meth", "EDNRB_meth"]].rename(columns=lambda c: c[:-5]), x_[["FHIT", "TFPI2", "DCC", "PAX1"]]], axis=1)
SRC = {"cBioPortal": cb, "A": A, "B": B}
def sp(a, b): return stats.spearmanr(a, b)[0]
def fisher(r, n): z, se = np.arctanh(r), 1 / np.sqrt(n - 3); return np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)

# ---------- item 4
PAIRS = [("DACH1", t) for t in ["CADM1", "ESR1", "SLIT2", "EDNRB"]] + [(r, t) for r in ["SOX8", "SALL2"] for t in ["DCC", "FHIT", "TFPI2"]]
expr = lambda g: d.DACH1_expr if g == "DACH1" else e[g]
rows, store = [], {}
for src, M in SRC.items():
    for r, t in PAIRS:
        if t not in M or M[t].isna().all(): rows.append(dict(source=src, regulator=r, target=t, n=0, note="no methylation values")); continue
        xv, yv = expr(r).values, M[t].values; ok = ~np.isnan(yv) & ~np.isnan(xv)
        xv, yv = xv[ok], yv[ok]
        rho = sp(xv, yv); lo, hi = fisher(rho, ok.sum())
        bi = rng.integers(0, ok.sum(), size=(NB, ok.sum()))
        bs = np.array([sp(xv[i], yv[i]) for i in bi]); store[(src, r, t)] = (rho, ok, bs)
        rows.append(dict(source=src, regulator=r, target=t, n=int(ok.sum()), rho=round(rho, 3), fisher_lo=round(lo, 3), fisher_hi=round(hi, 3),
                         boot_lo=round(np.percentile(bs, 2.5), 3), boot_hi=round(np.percentile(bs, 97.5), 3), p=f"{stats.spearmanr(xv, yv)[1]:.1e}"))
t4 = pd.DataFrame(rows); t4.to_csv(HERE / "item4_ci_table.csv", index=False)
print("== ITEM 4: rho with 95% CIs\n", t4.to_string(index=False))

# paired bootstrap of difference: rho(DACH1,T) - rho(SOX8,TFPI2), same resampled samples
drows = []
for src, M in SRC.items():
    for sox_t in ["TFPI2", "FHIT"]:
        if sox_t not in M or M[sox_t].isna().all(): continue
        for t in ["CADM1", "ESR1", "SLIT2"]:
            ok = ~(M[t].isna() | M[sox_t].isna()).values
            xd, yd = d.DACH1_expr.values[ok], M[t].values[ok]
            xs, ys = e["SOX8"].values[ok], M[sox_t].values[ok]
            m = ok.sum(); bi = rng.integers(0, m, size=(NB, m))
            diff = np.array([sp(xd[i], yd[i]) - sp(xs[i], ys[i]) for i in bi]); obs = sp(xd, yd) - sp(xs, ys)
            p_two = min(1.0, 2 * min((diff >= 0).mean(), (diff <= 0).mean()))
            drows.append(dict(source=src, comparison=f"DACH1-{t} minus SOX8-{sox_t}", n=int(m), diff=round(obs, 3),
                              boot_lo=round(np.percentile(diff, 2.5), 3), boot_hi=round(np.percentile(diff, 97.5), 3), boot_p_two_sided=round(p_two, 4)))
t4d = pd.DataFrame(drows); t4d.to_csv(HERE / "item4_sox8_difference_tests.csv", index=False)
print("\n== ITEM 4: is SOX8 (negative control) different from the DACH1 targets? (paired bootstrap, negative diff = DACH1 more inverse)\n", t4d.to_string(index=False))

# ---------- item 6
f = d.DACH1_expr; low = f < 10
print(f"\n== ITEM 6: DACH1 RSEM <10: {low.sum()}/{len(f)} = {100*low.mean():.1f}%; median split cutoff {f.median():.1f}")
r6 = []
for src, M in SRC.items():
    for t in ["CADM1", "ESR1", "SLIT2", "EDNRB"]:
        for cut_name, grp in (("RSEM<10", low), ("median split", f < f.median())):
            a_, b_ = M.loc[grp, t].dropna(), M.loc[~grp, t].dropna()
            u, p = stats.mannwhitneyu(a_, b_, alternative="two-sided")
            r6.append(dict(source=src, target=t, split=cut_name, n_low=len(a_), n_high=len(b_), median_beta_low=round(a_.median(), 3),
                           median_beta_high=round(b_.median(), 3), delta_median=round(a_.median() - b_.median(), 3),
                           AUC_low_gt_high=round(u / (len(a_) * len(b_)), 3), p=f"{p:.1e}"))
t6 = pd.DataFrame(r6); t6.to_csv(HERE / "item6_dach1_low_vs_high.csv", index=False)
print(t6.to_string(index=False))

# ---------- item 5: burden index with PAX1 (probe-level A for all genes)
BUR7 = ["CDKN2A", "RARB", "DAPK1", "CDH1", "MGMT", "MAL", "SOCS1"]
idx7, idx8 = A[BUR7].mean(axis=1), A[BUR7 + ["PAX1"]].mean(axis=1)
ok = pur.index.intersection(d.index)
print("\n== ITEM 5: burden index")
print(f"original index (cBioPortal, 7 genes): SD {d.burden_index.std():.3f}, range {d.burden_index.min():.2f}-{d.burden_index.max():.2f}, vs purity rho={sp(d.burden_index[ok], pur.loc[ok,'purity']):+.3f}")
for nm, ix in (("probe-level A, 7 genes", idx7), ("probe-level A, 8 genes + PAX1", idx8)):
    print(f"{nm}: SD {ix.std():.3f}, range {ix.min():.2f}-{ix.max():.2f}, vs purity rho={sp(ix[ok], pur.loc[ok,'purity']):+.3f}, vs original index rho={sp(ix, d.burden_index):+.3f}")
def pc(xv, yv, zv):
    rx, ry, rz = map(stats.rankdata, (xv, yv, zv)); res = lambda v: v - np.polyval(np.polyfit(rz, v, 1), rz)
    return stats.pearsonr(res(rx), res(ry))[0]
r5 = []
for t in ["CADM1", "ESR1", "SLIT2", "EDNRB"]:
    yv = A[t].values; xv = d.DACH1_expr.values
    r5.append(dict(target=t, naive_A=round(sp(xv, yv), 3), partial_7=round(pc(xv, yv, idx7.values), 3), partial_8_with_PAX1=round(pc(xv, yv, idx8.values), 3),
                   partial_orig_index=round(pc(xv, yv, d.burden_index.values), 3)))
t5 = pd.DataFrame(r5); t5.to_csv(HERE / "item5_burden_pax1.csv", index=False); print(t5.to_string(index=False))
