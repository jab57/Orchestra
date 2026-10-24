r"""DACH1 v3 replication in TCGA-HNSC: methylation correlation + matched-regulator null (rule fixed before this script was run, as described in the paper; single untuned run).

Same as gate3_matched_tf_null_dach1.py but on HNSC: DACH1 expression vs promoter methylation (definition A: mean beta of
TSS200 probes, TSS1500 fallback) of SLIT2, ESR1, CADM1 (EDNRB = control), compared with regulators of similar regulon size
(0.5-2x DACH1's HNSC regulon) and mean expression (+-1 log2 unit). Cohort: primary tumours (-01) present in both files.
Inputs (raw data not included; SHA-256 of the Xena files are in the paper's Data availability), all under C:\Dev\xena_data:
  TCGA-HNSC.star_tpm.tsv.gz, TCGA-HNSC.methylation450.tsv.gz, GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz,
  gencode_v36_gene_map.tsv; network: C:\Dev\RegNetAgents\models\networks\tcga\hnsc\network_index.pkl
Outputs (this folder): gate3_hnsc_results.json, gate3_hnsc_null_rhos.csv
Run:  python -I -X utf8 screening/gate3_hnsc_replication.py
"""
import csv, gzip, json, pickle
from collections import defaultdict
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
XENA = Path(r"C:\Dev\xena_data")
NET = Path(r"C:\Dev\RegNetAgents\models\networks\tcga\hnsc\network_index.pkl")
CAND, TARGETS, CONTROL = "DACH1", ["SLIT2", "ESR1", "CADM1"], "EDNRB"
K_LO, K_HI, EXPR_TOL, ALPHA = 0.5, 2.0, 1.0, 0.05
WANT = set(TARGETS + [CONTROL])

# promoter probes for the four targets (definition A), from the Illumina manifest annotation only
A = defaultdict(lambda: {"TSS200": set(), "TSS1500": set()})
with gzip.open(XENA / "GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz", "rt", encoding="latin-1") as f:
    for line in f:
        if line.startswith("IlmnID"):
            header = next(csv.reader([line])); break
    ix = {c: i for i, c in enumerate(header)}
    for row in csv.reader(f):
        if len(row) <= ix["UCSC_RefGene_Group"] or not row[ix["UCSC_RefGene_Name"]]: continue
        for g, grp in set(zip(row[ix["UCSC_RefGene_Name"]].split(";"), row[ix["UCSC_RefGene_Group"]].split(";"))):
            if g in WANT and grp in ("TSS200", "TSS1500"): A[g][grp].add(row[0])
defA = {g: sorted(v["TSS200"]) if v["TSS200"] else sorted(v["TSS1500"]) for g, v in A.items()}
print({g: len(p) for g, p in defA.items()})
need = {p for v in defA.values() for p in v}

with gzip.open(XENA / "TCGA-HNSC.methylation450.tsv.gz", "rt") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    keep = {}
    for line in f:
        pid, _, rest = line.partition("\t")
        if pid in need: keep[pid] = np.array([float(x) if x not in ("", "NA") else np.nan for x in rest.rstrip("\n").split("\t")], dtype="float64")
mcols = [c[:15] for c in hdr[1:]]
msamp = [(c, i) for i, c in enumerate(mcols) if c.endswith("-01")]
seen = {}
for c, i in msamp: seen.setdefault(c, i)
meth = pd.DataFrame({g: np.nanmean(np.vstack([keep[p] for p in ps if p in keep]), axis=0) for g, ps in defA.items()}).T
meth.columns = mcols
meth = meth.loc[:, ~meth.columns.duplicated()]
meth = meth[[c for c in meth.columns if c.endswith("-01")]]

gmap = pd.read_csv(XENA / "gencode_v36_gene_map.tsv", sep="\t", header=None, names=["eid", "sym", "type"])
expr = pd.read_csv(XENA / "TCGA-HNSC.star_tpm.tsv.gz", sep="\t", index_col=0)
expr.columns = [c[:15] for c in expr.columns]
expr = expr.loc[:, ~expr.columns.duplicated()]
expr = expr[[c for c in expr.columns if c.endswith("-01")]]
expr.index = expr.index.map(gmap.set_index("eid").sym)
expr = expr[expr.index.notna()].groupby(level=0).max()
found = [s for s in expr.columns if s in meth.columns]
expr, meth = expr[found], meth[found]
print(f"expression primary tumours {expr.shape[1]}, methylation {meth.shape[1]}, shared {len(found)}")

d = pickle.load(open(NET, "rb"))
K = {r.upper(): len({t.upper() for t in ts}) for r, ts in d["regulator_targets"].items() if ts}
reg = [r for r in K if r in expr.index]
print(f"regulators in network {len(K)}; with expression {len(reg)}; {CAND} K={K.get(CAND)}")

X = expr.loc[reg]; X = X[X.std(axis=1) > 0]
mean_expr = X.mean(axis=1)
R = X.rank(axis=1).values; R = R - R.mean(axis=1, keepdims=True); Rn = np.sqrt((R ** 2).sum(axis=1))

def rhos_for(t):
    y = meth.loc[t].values.astype(float)
    ok = ~np.isnan(y)
    assert ok.all(), f"{t}: {(~ok).sum()} missing"   # methylation missing -> would need per-sample handling
    y = stats.rankdata(y); y = y - y.mean()
    return pd.Series((R @ y) / (Rn * np.sqrt((y ** 2).sum())), index=X.index)

rho_df = pd.DataFrame({t: rhos_for(t) for t in TARGETS + [CONTROL]})
rho_df["K"] = [K[g] for g in rho_df.index]; rho_df["mean_expr"] = mean_expr
rho_df["mean_rho_3"] = rho_df[TARGETS].mean(axis=1)
c = rho_df.loc[CAND]
pool = rho_df[(rho_df.index != CAND) & rho_df.K.between(K_LO * c.K, K_HI * c.K) & ((rho_df.mean_expr - c.mean_expr).abs() <= EXPR_TOL)]
allreg = rho_df[rho_df.index != CAND]
print(f"{CAND}: K={int(c.K)}, mean={c.mean_expr:.2f}; matched pool n={len(pool)}; all regulators n={len(allreg)}")

out = {"samples": len(found), "candidate": CAND, "K": int(c.K), "mean_expr": round(float(c.mean_expr), 3),
       "matched_pool_n": int(len(pool)), "all_regulators_n": int(len(allreg)), "targets": {}}
for t in TARGETS + [CONTROL, "mean_rho_3"]:
    v = float(c[t])
    out["targets"][t] = {"rho_candidate": round(v, 3), "matched_median": round(float(pool[t].median()), 3),
        "matched_pct_more_negative": round(100 * float((pool[t] <= v).mean()), 2),
        "matched_emp_p_one_sided": round((float((pool[t] <= v).sum()) + 1) / (len(pool) + 1), 4),
        "all_reg_emp_p_one_sided": round((float((allreg[t] <= v).sum()) + 1) / (len(allreg) + 1), 4)}
out["targets_with_matched_p_le_alpha"] = [t for t in TARGETS if out["targets"][t]["matched_emp_p_one_sided"] <= ALPHA]
rho_df.sort_values("mean_rho_3").to_csv(HERE / "gate3_hnsc_null_rhos.csv")
(HERE / "gate3_hnsc_results.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print(json.dumps(out, indent=2))
