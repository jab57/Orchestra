"""DACH1 v3: larger-panel rerun of the enrichment screen (2026-10-23). Read-only. Single run, rule fixed before running.

Rule: DACH1 'holds up' if it stays in the top 10 by raw Fisher p AND its BH-adjusted p (over all regulators with >= 1 panel-gene target)
is <= 0.17 (the 24-gene value). Panel: the original 24 + 43 new genes, chosen on methylation-silencing literature without looking at results.
Reuses Network from gate1_stability_and_null.py (reproduces the tool's Fisher p; checked there). Permutation null: 2,000 random panels
of the same resolved size, seed 42, as in gate 1.
Run:  python -I -X utf8 screening/gate4_larger_panel.py     Writes gate4_larger_panel_results.json
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gate1_stability_and_null import Network, panel_indices

HERE = Path(__file__).resolve().parent
ORIG = ("CDKN2A CADM1 MAL DAPK1 RASSF1A CDH1 TIMP3 SOCS1 FHIT PAX1 WIF1 SFRP1 TFPI2 DLEC1 SCGB3A1 RARB APC MGMT ESR1 GSTP1 EDNRB SLIT2 HOXA9 DCC").split()
NEW = ("CCNA1 CDH13 CDKN2B DKK3 EPB41L3 HIC1 JAM3 MLH1 PAX5 RUNX3 SFRP2 SOX1 TWIST1 HS3ST2 ZNF582 SOX17 LMX1A NKX6-1 GATA4 GATA5 RARRES1 "
       "TMEFF2 CCND2 RASSF2 CXCL12 SFRP4 SFRP5 THBS1 KIF1A IRX1 NID2 TCF21 SYNE1 SHISA3 GRIA4 PCDH10 ID4 CHFR SYK TP73 PTEN BRCA1 STK11").split()
assert len(ORIG) == 24 and len(NEW) == 43 and not set(ORIG) & set(NEW)
N_PERM, SEED, CAND, BASE_BH = 2000, 42, "DACH1", 0.17


def bh(p):
    o = np.argsort(p); n = len(p); q = np.empty(n)
    r = p[o] * n / np.arange(1, n + 1)
    q[o] = np.minimum.accumulate(r[::-1])[::-1]
    return np.minimum(q, 1)


net = Network("cesc")
out = {}
for label, genes in (("original_24", ORIG), ("larger_67", ORIG + NEW)):
    idx = panel_indices(net, genes)
    keep, p = net.enrich(idx)
    q = bh(p)
    order = np.argsort(p)
    names = net.regs[keep]
    j = int(np.nonzero(names == CAND)[0][0])
    ov = np.asarray(net.M[idx, :].sum(axis=0)).ravel()[keep]
    top = [{"regulator": str(names[i]), "overlap": int(ov[i]), "raw_p": float(p[i]), "bh_p": float(q[i])} for i in order[:10]]
    rng = np.random.default_rng(SEED)
    best = np.array([net.enrich(rng.choice(net.N, size=len(idx), replace=False))[1].min() for _ in range(N_PERM)])
    out[label] = {"panel_size": len(genes), "resolved": len(idx),
                  "unresolved": [g for g in genes if g.upper() not in net.gidx],
                  "regulators_tested": int(len(keep)),
                  "DACH1": {"overlap": int(ov[j]), "regulon_size": int(net.K[keep][j]), "raw_p": float(p[j]),
                            "rank": 1 + int(np.sum(p < p[j])), "bh_p": float(q[j])},
                  "familywise_p_DACH1": float((np.sum(best <= p[j]) + 1) / (N_PERM + 1)),
                  "top10": top}
big = out["larger_67"]["DACH1"]
out["rule"] = {"top10": big["rank"] <= 10, "bh_not_worse_than_0.17": big["bh_p"] <= BASE_BH,
               "holds_up": bool(big["rank"] <= 10 and big["bh_p"] <= BASE_BH)}
(HERE / "gate4_larger_panel_results.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print(json.dumps(out, indent=1))
