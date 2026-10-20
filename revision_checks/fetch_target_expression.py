"""Fetch RNA-seq expression for the DACH1 targets (and SOX8/SALL2 + their targets) for the
253 squamous CESC samples already used in dach1_partial_corr_cesc.csv.
Writes target_expression_cesc.csv (new file; original analysis files untouched).
Run from c:\Dev\Orchestra:  python revision_checks/fetch_target_expression.py
"""
import csv, os, sys
from pathlib import Path
from dotenv import load_dotenv
load_dotenv()
sys.path.insert(0, os.getcwd())
from cbioportal_client import _TCGA_STUDY_IDS, _discover_profiles, _entrez_ids, _fetch_values_timed

HERE = Path(__file__).parent
ids = [r["sample_id"] for r in csv.DictReader(open("dach1_partial_corr_cesc.csv"))]
genes = ["DACH1", "CADM1", "ESR1", "SLIT2", "EDNRB", "SOX8", "SALL2", "FHIT", "TFPI2", "DCC"]
study = _TCGA_STUDY_IDS["cesc"]
rna, _ = _discover_profiles(study)
ent = _entrez_ids(genes)
print("study", study, "rna", rna, "unresolved:", [g for g in genes if g not in ent])
cols = {g: _fetch_values_timed(rna, ent[g], ids) for g in genes if g in ent}
for g, v in cols.items(): print(g, len(v))
with open(HERE / "target_expression_cesc.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["sample_id"] + list(cols))
    for s in ids: w.writerow([s] + [cols[g].get(s, "") for g in cols])
