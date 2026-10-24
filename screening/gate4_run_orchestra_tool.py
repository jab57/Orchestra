r"""Run Orchestra's own analyze_gene_signature workflow on the 67-gene panel (DACH1 v3 larger-panel check, 2026-10-23).

Calls OrchestraWorkflow.run_analysis (the code behind the analyze_gene_signature tool) on the TCGA CESC network. It opens the CASCADE and
RegNetAgents child servers the same way orchestra_mcp_server.main() does, because the workflow only runs once those connections are set.
Writes the full result as JSON; the copy kept with the paper is gate4_orchestra_tool_output_67genes.json.
Run from C:\Dev\Orchestra:  python -I -X utf8 screening/gate4_run_orchestra_tool.py <output.json>
"""
import asyncio, json, sys, io
sys.path.insert(0, r"C:\Dev\Orchestra")
from pathlib import Path
from dotenv import load_dotenv
load_dotenv(Path(r"C:\Dev\Orchestra") / ".env")
from orchestra_langgraph_workflow import OrchestraWorkflow
ORIG = "CDKN2A CADM1 MAL DAPK1 RASSF1A CDH1 TIMP3 SOCS1 FHIT PAX1 WIF1 SFRP1 TFPI2 DLEC1 SCGB3A1 RARB APC MGMT ESR1 GSTP1 EDNRB SLIT2 HOXA9 DCC".split()
NEW = ("CCNA1 CDH13 CDKN2B DKK3 EPB41L3 HIC1 JAM3 MLH1 PAX5 RUNX3 SFRP2 SOX1 TWIST1 HS3ST2 ZNF582 SOX17 LMX1A NKX6-1 GATA4 GATA5 RARRES1 "
       "TMEFF2 CCND2 RASSF2 CXCL12 SFRP4 SFRP5 THBS1 KIF1A IRX1 NID2 TCF21 SYNE1 SHISA3 GRIA4 PCDH10 ID4 CHFR SYK TP73 PTEN BRCA1 STK11").split()
from contextlib import AsyncExitStack
from mcp_client import make_cascade_client, make_regnetagents_client
async def main():
    wf = OrchestraWorkflow()
    async with AsyncExitStack() as stack:
        wf._persistent_cascade = await stack.enter_async_context(make_cascade_client())
        wf._persistent_regnetagents = await stack.enter_async_context(make_regnetagents_client())
        wf._persistent_ready.set()
        print("ready", flush=True)
        r = await wf.run_analysis(gene="", cell_type="epithelial_cell", analysis_type="gene_signature",
                                  gene_signature=ORIG+NEW, cancer_type="cesc")
    Path(sys.argv[1]).write_text(json.dumps(r, indent=1, default=str), encoding="utf-8")
    print("done", flush=True)
asyncio.run(main())
