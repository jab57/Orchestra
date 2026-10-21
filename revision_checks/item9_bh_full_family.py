"""Item 9: BH adjustment over the full family of regulators tested (DACH1/CESC, v3).

Orchestra versions before the fix corrected the Fisher p-values of only the 10 regulators shown;
the correct family is every regulator tested (616 in the TCGA CESC network). This reruns both
panels reported in the paper through the patched analyze_gene_signature path and tabulates, per
driver: raw Fisher p, the previously reported p-adj (BH over the top 10 only) and the corrected p-adj.

Run from c:\\Dev\\Orchestra:  python revision_checks/item9_bh_full_family.py
Writes item9_bh_full_family.csv next to this script.
"""
import asyncio, csv, io, subprocess, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]  # c:\Dev\Orchestra
sys.path.insert(0, str(REPO))
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
from dotenv import load_dotenv
load_dotenv(REPO / ".env")
from orchestra_langgraph_workflow import OrchestraWorkflow

PANEL_24 = ["CDKN2A", "CADM1", "MAL", "DAPK1", "RASSF1A", "CDH1", "TIMP3", "SOCS1", "FHIT", "PAX1",
            "WIF1", "SFRP1", "TFPI2", "DLEC1", "SCGB3A1", "RARB", "APC", "MGMT", "ESR1", "GSTP1",
            "EDNRB", "SLIT2", "HOXA9", "DCC"]
RUNS = {
    "orig23": [g for g in PANEL_24 if g != "RASSF1A"],
    "rassf1_24": [("RASSF1" if g == "RASSF1A" else g) for g in PANEL_24],
}
CONTEXTS = ["cervical cancer", "head and neck squamous"]


async def run(label):
    wf = OrchestraWorkflow()
    wf._persistent_ready.set()  # no MCP server here: use the per-call client path
    res = await wf.run_analysis(gene="", cell_type="epithelial_cell", analysis_type="gene_signature",
                                gene_signature=RUNS[label], cancer_contexts=CONTEXTS, cancer_type="cesc")
    syn = res.get("synthesis") or {}
    drivers = syn.get("ranked_drivers") or []
    n_tested = syn.get("total_regulators_tested")
    old = OrchestraWorkflow._bh_correct([d["p_value"] for d in drivers])  # pre-fix: BH over the 10 shown
    rows = []
    for d, o in zip(drivers, old):
        rows.append({"panel": label, "gene": d["gene"], "overlap": d["overlap_count"], "regulators_tested": n_tested,
                     "fisher_p": d["p_value"], "p_adj_top10_only_previously_reported": o,
                     "p_adj_all_tested": d["p_value_adj"]})
    return rows


async def main():
    out = []
    for label in RUNS:
        rows = await run(label)
        out += rows
        for r in rows:
            print(label, r["gene"], r["overlap"], f'{r["fisher_p"]:.3g}', f'{r["p_adj_top10_only_previously_reported"]:.3g}',
                  f'{r["p_adj_all_tested"]:.3g}')
    p = Path(__file__).parent / "item9_bh_full_family.csv"
    with p.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print("git:", subprocess.run(["git", "-C", str(REPO), "describe", "--tags", "--always", "--dirty"],
                                 capture_output=True, text=True).stdout.strip())
    print("wrote", p)


if __name__ == "__main__":
    asyncio.run(main())
