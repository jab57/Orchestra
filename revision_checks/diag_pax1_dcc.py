r"""Item 5 diagnostics: why is PAX1 methylation 'unavailable' and DCC n=0 in the cBioPortal-based analysis?
Checks (a) Entrez resolution + cBioPortal methylation values for the squamous samples, (b) HM450 manifest annotation,
(c) presence in the probe-level Xena matrix (via gene_promoter_meth_cesc_{A,B}). Also writes cBioPortal methylation for
FHIT/TFPI2/DCC/PAX1 to cbioportal_extra_meth_cesc.csv (needed for the SOX8/SALL2 CIs).
Run from c:\Dev\Orchestra: python -X utf8 <this file>"""
import csv, os, sys, warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import pandas as pd
from dotenv import load_dotenv
load_dotenv(); sys.path.insert(0, os.getcwd())
from cbioportal_client import _TCGA_STUDY_IDS, _discover_profiles, _entrez_ids, _fetch_values_timed
HERE = Path(__file__).parent
ids = [r["sample_id"] for r in csv.DictReader(open("dach1_partial_corr_cesc.csv"))]
genes = ["FHIT", "TFPI2", "DCC", "PAX1", "RASSF1", "RASSF1A"]
study = _TCGA_STUDY_IDS["cesc"]; _, meth = _discover_profiles(study)
ent = _entrez_ids(genes); print("Entrez:", ent)
out = {}
for g in genes:
    if g not in ent: print(g, "-> no Entrez ID"); continue
    v = _fetch_values_timed(meth, ent[g], ids); out[g] = v
    print(f"{g}: cBioPortal methylation values for {len(v)}/{len(ids)} samples")
pd.DataFrame({g: pd.Series(v) for g, v in out.items()}).reindex(ids).rename_axis("sample_id").to_csv(HERE / "cbioportal_extra_meth_cesc.csv")
for tag in ("A_tss200_fallback1500", "B_tss1500_tss200_exon1"):
    m = pd.read_csv(HERE / f"gene_promoter_meth_cesc_{tag}.csv.gz", index_col=0)
    for g in genes:
        print(f"probe-level {tag}: {g} ->", "absent" if g not in m.index else f"{m.loc[g].notna().sum()} samples, mean beta {m.loc[g].mean():.2f}")
