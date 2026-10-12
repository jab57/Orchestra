"""
Ground-truth gene sets for the corroboration-paper scoring scripts.

    intogen     IntOGen Compendium of Mutational Cancer Driver Genes, release
                2024.09.20 (CC0). Snapshot shipped with RegNetAgents. Primary.
    vogelstein  Vogelstein et al. 2013, Table S2A (mutation-driver genes).
                Secondary. Snapshot lives in the private repo (not open access).
    oncokb      LEGACY. Reads the cached OncoKB list in outputs/. Kept only so the
                v1-v3 analyses remain reproducible; never used for new results.
    sanger_cgc  LEGACY. Reads the cached Sanger CGC list in outputs/. Same status.

The intogen and vogelstein snapshots are checked against pinned SHA-256 hashes
(recorded in private/manuscript/arxiv/corroboration/ground_truth_switch_design.md
and reference_data/README.md) so a changed file cannot be scored silently.
"""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUTPUTS = ROOT / "outputs"

INTOGEN_PATH = Path(r"c:\Dev\RegNetAgents\regnetagents\reference_data\intogen_drivers.tsv")
INTOGEN_SHA256 = "3668daacc6297d332a1cd84234b5390e4052a8ac63665698c7dd05783a79a05a"
INTOGEN_N = 633

VOGELSTEIN_PATH = (ROOT / "private" / "manuscript" / "arxiv" / "corroboration"
                   / "reference_data" / "vogelstein2013_table_s2a_mut_drivers.tsv")
VOGELSTEIN_SHA256 = "8d9aafb60ed47cfc01286a14f5007c5cef1805877e0ddb027014787481facd4b"
VOGELSTEIN_N = 125

LEGACY = ("oncokb", "sanger_cgc")
CHOICES = ("intogen", "vogelstein") + LEGACY


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _check_snapshot(path, expected_sha, label):
    if not Path(path).exists():
        raise FileNotFoundError(f"{label} snapshot not found: {path}")
    got = _sha256(path)
    if got != expected_sha:
        raise ValueError(f"{label} snapshot hash mismatch: expected {expected_sha}, got {got}")


def _read_symbols(path, column):
    with open(path, encoding="utf-8", newline="") as f:
        lines = (ln for ln in f if not ln.startswith("#"))
        return {row[column].strip() for row in csv.DictReader(lines, delimiter="\t")
                if row.get(column, "").strip()}


def load_truth(name):
    """Return the set of gene symbols counted as 'positive' for `name`.

    Symbols are returned exactly as stored (no case folding), matching how the
    original scripts compared candidate symbols against OncoKB.
    """
    if name == "intogen":
        _check_snapshot(INTOGEN_PATH, INTOGEN_SHA256, "IntOGen")
        genes = _read_symbols(INTOGEN_PATH, "symbol")
        if len(genes) != INTOGEN_N:
            raise ValueError(f"IntOGen: expected {INTOGEN_N} genes, got {len(genes)}")
        return genes
    if name == "vogelstein":
        _check_snapshot(VOGELSTEIN_PATH, VOGELSTEIN_SHA256, "Vogelstein")
        genes = _read_symbols(VOGELSTEIN_PATH, "gene")
        if len(genes) != VOGELSTEIN_N:
            raise ValueError(f"Vogelstein: expected {VOGELSTEIN_N} genes, got {len(genes)}")
        return genes
    if name == "oncokb":
        data = json.load(open(OUTPUTS / "oncokb_cancer_gene_list.json", encoding="utf-8"))
        return {x["hugoSymbol"] for x in data}
    if name == "sanger_cgc":
        data = json.load(open(OUTPUTS / "sanger_cgc_genes.json", encoding="utf-8"))
        return set(data)
    raise ValueError(f"unknown ground truth {name!r}; choose from {CHOICES}")


def output_suffix(name, default="oncokb"):
    """'' for the legacy default (so v1-v3 outputs keep their filenames), else '_<name>'."""
    return "" if name == default else f"_{name}"
