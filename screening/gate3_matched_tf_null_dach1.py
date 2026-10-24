r"""Gate 3 (matched-TF null) for DACH1/CESC.

Question: DACH1 expression correlates negatively with promoter methylation of SLIT2, ESR1, CADM1 (EDNRB = negative
control). Is that more negative than for random regulators of similar regulon size and expression level, tested against
the same targets?

Inputs (raw data not included; SHA-256 of the Xena files are in the paper's Data availability):
  C:\Dev\xena_data\TCGA-CESC.star_tpm.tsv.gz       Xena GDC hub, log2(TPM+1), Ensembl rows, -01A sample IDs
  C:\Dev\xena_data\gencode_v36_gene_map.tsv         GENCODE v36 id -> symbol (the GDC STAR pipeline's annotation)
  C:\Dev\RegNetAgents\models\networks\tcga\cesc\network_index.pkl   regulator -> targets (regulon sizes)
  revision_checks\gene_promoter_meth_cesc_A_tss200_fallback1500.csv.gz   (promoter definition A; built by revision_checks/build_gene_promoter_matrix.py, not included)
  C:\Dev\Orchestra\dach1_partial_corr_cesc.csv      the 253-sample squamous cohort IDs
Outputs (this folder): gate3_results.json, gate3_null_rhos.csv
Run:  python -I -X utf8 screening/gate3_matched_tf_null_dach1.py
Thresholds are fixed here BEFORE looking at results: match on regulon size within x0.5..x2 of DACH1's and mean
expression within +-1 (log2 units); a target counts as "unusual" if its one-sided empirical p <= 0.05.
"""
import json, pickle
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
XENA = Path(r"C:\Dev\xena_data")
NET = Path(r"C:\Dev\RegNetAgents\models\networks\tcga\cesc\network_index.pkl")
REV = HERE.parent / "revision_checks"
COHORT = Path(r"C:\Dev\Orchestra\dach1_partial_corr_cesc.csv")
CAND, TARGETS, CONTROL = "DACH1", ["SLIT2", "ESR1", "CADM1"], "EDNRB"
K_LO, K_HI, EXPR_TOL, ALPHA = 0.5, 2.0, 1.0, 0.05

ids = pd.read_csv(COHORT, usecols=["sample_id"]).sample_id.tolist()
gmap = pd.read_csv(XENA / "gencode_v36_gene_map.tsv", sep="\t", header=None, names=["eid", "sym", "type"])
expr = pd.read_csv(XENA / "TCGA-CESC.star_tpm.tsv.gz", sep="\t", index_col=0)
expr.columns = [c[:15] for c in expr.columns]                      # TCGA-XX-XXXX-01A -> TCGA-XX-XXXX-01
expr = expr.loc[:, ~expr.columns.duplicated()]
found = [s for s in ids if s in expr.columns]
print(f"cohort samples {len(ids)}; found in Xena expression {len(found)}")
expr = expr[found]
sym = gmap.set_index("eid").sym
expr.index = expr.index.map(sym)
expr = expr[expr.index.notna()]
expr = expr.groupby(level=0).max()                                 # collapse duplicate symbols (rare)
meth = pd.read_csv(REV / "gene_promoter_meth_cesc_A_tss200_fallback1500.csv.gz", index_col=0)[found]

d = pickle.load(open(NET, "rb"))
K = {r.upper(): len({t.upper() for t in ts}) for r, ts in d["regulator_targets"].items() if ts}
reg = [r for r in K if r in expr.index]
print(f"regulators in network {len(K)}; with expression {len(reg)}; {CAND} K={K.get(CAND)}")

X = expr.loc[reg]
X = X[X.std(axis=1) > 0]
mean_expr = X.mean(axis=1)
R = X.rank(axis=1).values
R = R - R.mean(axis=1, keepdims=True)
Rn = np.sqrt((R ** 2).sum(axis=1))

def rhos_for(t):
    y = stats.rankdata(meth.loc[t].values.astype(float)); y = y - y.mean()
    return pd.Series((R @ y) / (Rn * np.sqrt((y ** 2).sum())), index=X.index)

rho = {t: rhos_for(t) for t in TARGETS + [CONTROL]}
rho_df = pd.DataFrame(rho); rho_df["K"] = [K[g] for g in rho_df.index]; rho_df["mean_expr"] = mean_expr
rho_df["mean_rho_3"] = rho_df[TARGETS].mean(axis=1)

c = rho_df.loc[CAND]
pool = rho_df[(rho_df.index != CAND) & rho_df.K.between(K_LO * c.K, K_HI * c.K) & ((rho_df.mean_expr - c.mean_expr).abs() <= EXPR_TOL)]
allreg = rho_df[rho_df.index != CAND]
print(f"{CAND}: K={int(c.K)}, mean log2(TPM+1)={c.mean_expr:.2f}; matched pool n={len(pool)}; all regulators n={len(allreg)}")

out = {"samples": len(found), "candidate": CAND, "K": int(c.K), "mean_expr": round(float(c.mean_expr), 3),
       "matched_pool_n": int(len(pool)), "all_regulators_n": int(len(allreg)),
       "match_rule": f"K in [{K_LO}x,{K_HI}x] of candidate's; mean expr within +-{EXPR_TOL}", "targets": {}}
for t in TARGETS + [CONTROL, "mean_rho_3"]:
    v = float(c[t])
    out["targets"][t] = {
        "rho_candidate": round(v, 3),
        "matched_median": round(float(pool[t].median()), 3), "matched_sd": round(float(pool[t].std()), 3),
        "matched_pct_more_negative": round(100 * float((pool[t] <= v).mean()), 2),
        "matched_emp_p_one_sided": round((float((pool[t] <= v).sum()) + 1) / (len(pool) + 1), 4),
        "all_reg_emp_p_one_sided": round((float((allreg[t] <= v).sum()) + 1) / (len(allreg) + 1), 4)}
out["targets_with_matched_p_le_alpha"] = [t for t in TARGETS if out["targets"][t]["matched_emp_p_one_sided"] <= ALPHA]
rho_df.sort_values("mean_rho_3").to_csv(HERE / "gate3_null_rhos.csv")
(HERE / "gate3_results.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print(json.dumps(out, indent=2))
