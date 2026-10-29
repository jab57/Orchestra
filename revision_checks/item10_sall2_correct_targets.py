r"""Item 10: SALL2 correlations against SALL2's actual overlap genes (CADM1, SCGB3A1, WIF1).
Earlier S3a/Results ran SALL2 against SOX8's genes (DCC, FHIT, TFPI2) by mistake; the tool output (item5_rassf1_rerun_orig23.json) lists
CADM1, SCGB3A1, WIF1. Same method as items4_5_6_ci_burden_bimodal.py (Spearman, Fisher CI, 5000-resample bootstrap, same seed family).
Sources: cBioPortal gene-level methylation (paper's values) and probe-level A (TSS200, TSS1500 fallback).
Run from c:\Dev\Orchestra:  python -X utf8 <this file>   Output: item10_sall2_correct_targets.csv, cbioportal_sall2_targets_meth_cesc.csv"""
import csv, os, sys, warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from dotenv import load_dotenv
load_dotenv(); sys.path.insert(0, os.getcwd())
from cbioportal_client import _TCGA_STUDY_IDS, _discover_profiles, _entrez_ids, _fetch_values_timed
HERE = Path(__file__).parent; rng = np.random.default_rng(20261029); NB = 5000
d = pd.read_csv("dach1_partial_corr_cesc.csv").set_index("sample_id")
e = pd.read_csv(HERE / "target_expression_cesc.csv").set_index("sample_id")
A = pd.read_csv(HERE / "gene_promoter_meth_cesc_A_tss200_fallback1500.csv.gz", index_col=0)[d.index].T
ids = list(d.index); study = _TCGA_STUDY_IDS["cesc"]; _, meth = _discover_profiles(study)
genes = ["SCGB3A1", "WIF1"]; ent = _entrez_ids(genes); print("Entrez:", ent)
fetched = {g: _fetch_values_timed(meth, ent[g], ids) for g in genes if g in ent}
for g, v in fetched.items(): print(f"{g}: cBioPortal methylation for {len(v)}/{len(ids)} samples")
cbx = pd.DataFrame({g: pd.Series(v) for g, v in fetched.items()}).reindex(ids).rename_axis("sample_id")
cbx.to_csv(HERE / "cbioportal_sall2_targets_meth_cesc.csv")
cb = pd.concat([d[["CADM1_meth"]].rename(columns={"CADM1_meth": "CADM1"}), cbx], axis=1)
sp = lambda a, b: stats.spearmanr(a, b)[0]
rows = []
for srcname, M in {"cBioPortal": cb, "A": A}.items():
    for t in ["CADM1", "SCGB3A1", "WIF1"]:
        if t not in M or M[t].isna().all(): rows.append(dict(source=srcname, regulator="SALL2", target=t, n=0, note="no methylation values")); continue
        xv, yv = e["SALL2"].values, M[t].reindex(ids).values; ok = ~np.isnan(xv) & ~np.isnan(yv); xv, yv = xv[ok], yv[ok]
        rho, p = stats.spearmanr(xv, yv); n = ok.sum(); se = 1 / np.sqrt(n - 3)
        bi = rng.integers(0, n, size=(NB, n)); bs = np.array([sp(xv[i], yv[i]) for i in bi])
        rows.append(dict(source=srcname, regulator="SALL2", target=t, n=int(n), rho=round(rho, 3),
                         fisher_lo=round(np.tanh(np.arctanh(rho) - 1.96 * se), 3), fisher_hi=round(np.tanh(np.arctanh(rho) + 1.96 * se), 3),
                         boot_lo=round(np.percentile(bs, 2.5), 3), boot_hi=round(np.percentile(bs, 97.5), 3), p=f"{p:.1e}"))
out = pd.DataFrame(rows); out.to_csv(HERE / "item10_sall2_correct_targets.csv", index=False); print(out.to_string(index=False))
