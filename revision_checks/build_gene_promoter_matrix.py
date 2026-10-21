r"""Build expression-independent gene-level promoter methylation for TCGA-CESC (items 1, 3, 5).

Inputs (argv): <TCGA-CESC.methylation450.tsv.gz (UCSC Xena GDC hub, hg38 betas)>
               <GPL13534 HM450 manifest .csv.gz (GEO platform GPL13534)> <dach1_partial_corr_cesc.csv>
Output (this folder): gene_promoter_meth_cesc.csv.gz  (rows = genes, cols = the 253 squamous sample IDs)
                      promoter_definitions.json       (probe counts per gene and definition)

Two promoter definitions, both chosen from the Illumina annotation only (never from expression):
  A_tss200_fallback1500 : mean beta of probes annotated TSS200; TSS1500 only if the gene has none
                          (the rule used for the GSE68339 replication in the v3 draft)
  B_tss1500_tss200_exon1: mean beta of all probes annotated TSS1500, TSS200 or 1stExon
                          (the definition proposed in the reviewer-style critique)
A probe counts for a gene if the gene symbol appears in UCSC_RefGene_Name with that group at the same index.
Run:  python -I -X utf8 build_gene_promoter_matrix.py <meth.tsv.gz> <manifest.csv.gz> <dach1_partial_corr_cesc.csv>
"""
import csv, gzip, io, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).parent
meth_path, man_path = sys.argv[1], sys.argv[2]
ids = pd.read_csv(sys.argv[3], usecols=["sample_id"]).sample_id.tolist()

A, B = defaultdict(lambda: {"TSS200": set(), "TSS1500": set()}), defaultdict(set)
with gzip.open(man_path, "rt", encoding="latin-1") as f:
    for line in f:
        if line.startswith("IlmnID"): header = next(csv.reader([line])); break
    ix = {c: i for i, c in enumerate(header)}
    for row in csv.reader(f):
        if len(row) <= ix["UCSC_RefGene_Group"] or not row[ix["UCSC_RefGene_Name"]]: continue
        names, groups = row[ix["UCSC_RefGene_Name"]].split(";"), row[ix["UCSC_RefGene_Group"]].split(";")
        for g, grp in set(zip(names, groups)):
            if grp in ("TSS200", "TSS1500"): A[g][grp].add(row[0])
            if grp in ("TSS200", "TSS1500", "1stExon"): B[g].add(row[0])
defA = {g: (sorted(v["TSS200"]) if v["TSS200"] else sorted(v["TSS1500"])) for g, v in A.items() if v["TSS200"] or v["TSS1500"]}
defB = {g: sorted(v) for g, v in B.items()}
need = set(p for v in defA.values() for p in v) | set(p for v in defB.values() for p in v)
print(f"genes: A={len(defA)} B={len(defB)}; promoter probes needed: {len(need)}")

# stream the (probes x samples) matrix, keep only needed probes
cols = None; keep = {}
with gzip.open(meth_path, "rt") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    sel = [i for i, c in enumerate(hdr) if i > 0 and c[:15] in set(ids)]
    names = [hdr[i][:15] for i in sel]
    seen = {}
    for j, n in enumerate(names): seen.setdefault(n, j)   # first vial per sample
    uniq = sorted(seen.values())
    print(f"matrix samples {len(hdr)-1}; matched to cohort ids: {len(seen)} / {len(ids)}")
    for line in f:
        pid, _, rest = line.partition("\t")
        if pid not in need: continue
        v = rest.rstrip("\n").split("\t")
        keep[pid] = np.array([float(v[sel[j] - 1]) if v[sel[j] - 1] not in ("", "NA") else np.nan for j in uniq], dtype="float32")  # -1: probe id column was split off into pid
samples = [names[j] for j in uniq]
print("probes kept:", len(keep))

def agg(defn):
    out = {}
    for g, probes in defn.items():
        m = [keep[p] for p in probes if p in keep]
        if m: out[g] = np.nanmean(np.vstack(m), axis=0)
    return pd.DataFrame(out, index=samples).T
for tag, defn in (("A_tss200_fallback1500", defA), ("B_tss1500_tss200_exon1", defB)):
    df = agg(defn); df.to_csv(HERE / f"gene_promoter_meth_cesc_{tag}.csv.gz"); print(tag, df.shape)
json.dump({"A_probe_counts": {g: len(p) for g, p in defA.items()}, "B_probe_counts": {g: len(p) for g, p in defB.items()}},
          open(HERE / "promoter_definitions.json", "w"))
