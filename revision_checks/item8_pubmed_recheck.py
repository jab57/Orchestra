"""Item 8: recheck 'zero PubMed records' for DACH1 in cervical / head-and-neck SCC.
Read-only PubMed E-utilities queries; writes item8_pubmed_results.tsv."""
import json, time, urllib.parse, urllib.request, sys
from pathlib import Path

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")

def get(path, **p):
    time.sleep(0.5)
    url = BASE + path + "?" + urllib.parse.urlencode(p)
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.load(r)

QUERIES = {
    "DACH1 + cervical": "DACH1[tiab] AND (cervical[tiab] OR cervix[tiab] OR CESC[tiab])",
    "DACH1 + head/neck terms": "DACH1[tiab] AND (head and neck[tiab] OR HNSC[tiab] OR laryn*[tiab] OR oral[tiab] OR oropharyn*[tiab] OR hypopharyn*[tiab] OR pharyn*[tiab] OR tongue[tiab] OR nasopharyn*[tiab])",
    "DACH1 + squamous": "DACH1[tiab] AND squamous[tiab]",
    "DACH1 + methylation + squamous/HN": "DACH1[tiab] AND methylat*[tiab] AND (squamous[tiab] OR head and neck[tiab] OR oral[tiab] OR laryn*[tiab])",
    "DACH1 (MeSH/all fields) + cervical neoplasms": "DACH1[All Fields] AND (Uterine Cervical Neoplasms[MeSH] OR Head and Neck Neoplasms[MeSH])",
}
rows = []
for name, term in QUERIES.items():
    d = get("esearch.fcgi", db="pubmed", retmode="json", retmax=200, term=term)["esearchresult"]
    ids = d["idlist"]
    summ = {}
    if ids:
        s = get("esummary.fcgi", db="pubmed", retmode="json", id=",".join(ids))["result"]
        summ = {i: s[i] for i in ids}
    print(f"\n== {name}: {d['count']} records\n   {term}")
    for i in ids:
        r = summ[i]
        doi = next((a["value"] for a in r.get("articleids", []) if a["idtype"] == "doi"), "")
        print(f"   {i} | {r.get('pubdate','')[:4]} | {r.get('source','')} | {r['title'][:110]} | {doi}")
        rows.append((name, i, r.get("pubdate", "")[:4], r.get("source", ""), r["title"], doi))
# the specific papers named in the critique
for pmid in ["37243970"]:
    r = get("esummary.fcgi", db="pubmed", retmode="json", id=pmid)["result"].get(pmid)
    print("\nNamed PMID", pmid, "->", r["title"] if r else "not found", "|", r.get("source") if r else "")
(OUT / "item8_pubmed_results.tsv").write_text(
    "query\tpmid\tyear\tjournal\ttitle\tdoi\n" + "\n".join("\t".join(x) for x in rows), encoding="utf-8")
