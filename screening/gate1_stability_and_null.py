"""Screen gate 1 (stability) and the permutation null for a nominated regulator.

Public copy: the DACH1/CESC case only.

For each (panel, TCGA network, candidate) this script:
  0. Reproduces RegNetAgents' one-sided Fisher enrichment from the saved TCGA network (hypergeometric
     survival function) and checks it against the values the tool returned for the real panel.
  1. Leave-one-panel-gene-out: removes each panel gene in turn, recomputes enrichment, and records the
     candidate's raw p and its rank by raw p among all regulators. Stable = candidate stays in the top 3
     in >= 80% of reruns.
  2. Permutation null: N_PERM random gene panels of the same size drawn from the network's gene universe
     (seed 42); for each, the best (minimum) raw Fisher p across all regulators with >= 1 overlapping
     target. Family-wise p = fraction of random panels whose best p is <= the candidate's p.

Read-only. Run from anywhere:  python screening/gate1_stability_and_null.py
Writes gate1_results.json and gate1_leave_one_out.csv next to this script.
"""
import csv, json, pickle, sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import hypergeom

NETDIR = Path(r"C:\Dev\RegNetAgents\models\networks\tcga")
HERE = Path(__file__).resolve().parent
N_PERM = 2000
SEED = 42

CASES = {
    "DACH1/cesc": dict(
        net="cesc", candidate="DACH1", best_other="AR",
        panel=["CDKN2A", "CADM1", "MAL", "DAPK1", "CDH1", "TIMP3", "SOCS1", "FHIT", "PAX1", "WIF1", "SFRP1",
               "TFPI2", "DLEC1", "SCGB3A1", "RARB", "APC", "MGMT", "ESR1", "GSTP1", "EDNRB", "SLIT2", "HOXA9", "DCC"],
        tool_p={"AR": 0.000177, "SOX8": 0.000701, "DACH1": 0.00116},
    ),
}


class Network:
    def __init__(self, net):
        d = pickle.load(open(NETDIR / net / "network_index.pkl", "rb"))
        self.genes = sorted(set(g.upper() for g in d["all_genes"]))
        self.gidx = {g: i for i, g in enumerate(self.genes)}
        regs, rows, cols = [], [], []
        for reg, targets in d["regulator_targets"].items():
            ts = {t.upper() for t in targets if t.upper() in self.gidx}
            if not ts:
                continue
            j = len(regs)
            regs.append(reg.upper())
            for t in ts:
                rows.append(self.gidx[t]); cols.append(j)
        self.regs = np.array(regs)
        self.M = sparse.csr_matrix((np.ones(len(rows), dtype=np.int32), (rows, cols)), shape=(len(self.genes), len(regs)))
        self.K = np.asarray(self.M.sum(axis=0)).ravel()  # regulon sizes
        self.N = len(self.genes)

    def enrich(self, panel_idx):
        """Raw one-sided Fisher p for every regulator with >= 1 overlap. Returns (reg_indices, p)."""
        v = np.zeros(self.N, dtype=np.int32); v[panel_idx] = 1
        ov = np.asarray(self.M.T.dot(v)).ravel()
        keep = np.nonzero(ov > 0)[0]
        k = len(panel_idx)
        p = hypergeom.sf(ov[keep] - 1, self.N, self.K[keep], k)
        return keep, p


def panel_indices(net, genes):
    idx = [net.gidx[g.upper()] for g in genes if g.upper() in net.gidx]
    return idx


def rank_of(net, keep, p, name):
    j = np.nonzero(net.regs[keep] == name.upper())[0]
    if len(j) == 0:
        return None, None
    pj = p[j[0]]
    return pj, 1 + int(np.sum(p < pj))


def main():
    results, loo_rows = {}, []
    rng = np.random.default_rng(SEED)
    for label, c in CASES.items():
        net = Network(c["net"])
        idx = panel_indices(net, c["panel"])
        keep, p = net.enrich(idx)
        out = {"network_genes": net.N, "regulators_with_regulon": len(net.regs), "panel_resolved": len(idx),
               "regulators_with_overlap": int(len(keep))}
        # 0. reproduce the tool
        chk = {}
        for g, tp in c["tool_p"].items():
            pj, rk = rank_of(net, keep, p, g)
            chk[g] = {"tool_p": tp, "replica_p": float(pj), "replica_rank": rk,
                      "match": bool(abs(pj - tp) / tp < 0.02)}
        out["reproduces_tool"] = chk
        cand_p, cand_rank = rank_of(net, keep, p, c["candidate"])
        other_p, _ = rank_of(net, keep, p, c["best_other"])
        out["candidate_p"], out["candidate_rank"] = float(cand_p), cand_rank
        # 1. leave-one-out
        ranks, ps = [], []
        for g in c["panel"]:
            sub = [i for i in idx if net.genes[i] != g.upper()]
            kp, pp = net.enrich(sub)
            pj, rk = rank_of(net, kp, pp, c["candidate"])
            ranks.append(rk); ps.append(None if pj is None else float(pj))
            loo_rows.append([label, g, c["candidate"], rk, "" if pj is None else f"{pj:.3g}"])
        ok = [r for r in ranks if r is not None]
        out["leave_one_out"] = {
            "n_reruns": len(ranks),
            "candidate_in_top3": int(sum(1 for r in ok if r <= 3)),
            "fraction_in_top3": sum(1 for r in ok if r <= 3) / len(ranks),
            "worst_rank": max(ok) if ok else None,
            "reruns_where_candidate_has_no_overlap": int(sum(1 for r in ranks if r is None)),
            "passes_80pct_top3": (sum(1 for r in ok if r <= 3) / len(ranks)) >= 0.8,
        }
        # 2. permutation null (best p across regulators, random panels of the same size)
        best = np.empty(N_PERM)
        for t in range(N_PERM):
            rp = rng.choice(net.N, size=len(idx), replace=False)
            _, pp = net.enrich(rp)
            best[t] = pp.min()
        out["null"] = {
            "n_perm": N_PERM,
            "median_best_p_random": float(np.median(best)),
            "p05_best_p_random": float(np.quantile(best, 0.05)),
            "familywise_p_for_candidate_p": float((np.sum(best <= cand_p) + 1) / (N_PERM + 1)),
            "familywise_p_for_best_other_p": float((np.sum(best <= other_p) + 1) / (N_PERM + 1)),
            "best_other": c["best_other"], "best_other_p": float(other_p),
        }
        results[label] = out
        print(label, json.dumps(out, indent=1, default=str))
    (HERE / "gate1_results.json").write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
    with open(HERE / "gate1_leave_one_out.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["case", "panel_gene_removed", "candidate", "candidate_rank_by_raw_p", "candidate_raw_p"])
        w.writerows(loo_rows)


if __name__ == "__main__":
    main()
