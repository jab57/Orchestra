"""Item 13 / point 5: rerun the DACH1 panel enrichment through Orchestra's analyze_gene_signature.

Run 1 = the ORIGINAL 23 network-resolved genes (RASSF1A omitted). Must reproduce Table 1 /
Supplementary Table S1 (same drivers, ranks, p-values) before Run 2 is trusted.
Run 2 = the same panel with RASSF1 in place of RASSF1A (24 genes).
Writes item5_rassf1_rerun_<label>.json next to this script. Read-only against the networks.
"""
import asyncio, json, subprocess, sys, io
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


def git(path, *a):
    return subprocess.run(["git", "-C", path, *a], capture_output=True, text=True).stdout.strip()


async def main(label):
    wf = OrchestraWorkflow()
    wf._persistent_ready.set()  # no MCP server here: use the per-call client path, not the "warming up" stub
    res = await wf.run_analysis(gene="", cell_type="epithelial_cell", analysis_type="gene_signature",
                                gene_signature=RUNS[label], cancer_contexts=CONTEXTS, cancer_type="cesc")
    out = {
        "label": label, "genes": RUNS[label], "cancer_contexts": CONTEXTS, "tcga_network": "cesc",
        "versions": {p: git(p, "describe", "--tags", "--always") for p in
                     [str(REPO), r"C:\Dev\RegNetAgents", r"C:\Dev\CASCADE"]},
        "query_summary": (res.get("master_regulators") or {}).get("query_summary"),
        "ranked_drivers": (res.get("synthesis") or {}).get("ranked_drivers"),
        "errors": res.get("errors"),
    }
    p = Path(__file__).parent / f"item5_rassf1_rerun_{label}.json"
    p.write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
    print("wrote", p)


if __name__ == "__main__":
    asyncio.run(main(sys.argv[1]))
