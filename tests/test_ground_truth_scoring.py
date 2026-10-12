"""
Tests for the ground-truth-switch scoring code in scripts/ (ground_truth.py,
experiment_threshold_tables.py and the --truth option of the incremental-value and
cross-check scripts).

Synthetic data only: gene names are made up (G001...), nothing here reads the real
corroboration candidate files or scores any real result. The one exception is the
snapshot-integrity test, which only checks that the pinned IntOGen / Vogelstein files
load with their recorded gene counts and hashes.
"""
import importlib
import json
import sys
from pathlib import Path

import pytest
from scipy.stats import fisher_exact

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))

ground_truth = importlib.import_module("ground_truth")
tables = importlib.import_module("experiment_threshold_tables")
cgc = importlib.import_module("experiment_cgc_crosscheck")
dorothea = importlib.import_module("experiment_dorothea_incremental_value")
mi_weight = importlib.import_module("experiment_mi_weight_incremental_value")


def rows(spec):
    """spec: list of (gene, category, count)."""
    return [{"gene": g, "category": c, "count": n} for g, c, n in spec]


def synthetic_panel():
    # focal: 6 corroborated (>=2) of which 4 positive; 10 uncorroborated of which 2 positive
    spec = []
    truth = set()
    for i in range(6):
        spec.append((f"F{i:02d}", "focal", 2 + (i % 2)))
        if i < 4:
            truth.add(f"F{i:02d}")
    for i in range(10):
        spec.append((f"U{i:02d}", "focal", i % 2))  # counts 0/1
        if i < 2:
            truth.add(f"U{i:02d}")
    # negative controls: 3 corroborated (1 positive), 12 uncorroborated (3 positive)
    for i in range(3):
        spec.append((f"N{i:02d}", "negative", 2))
        if i == 0:
            truth.add(f"N{i:02d}")
    for i in range(12):
        spec.append((f"M{i:02d}", "negative", 0))
        if i < 3:
            truth.add(f"M{i:02d}")
    return rows(spec), truth


def test_fisher_matches_scipy_on_known_table():
    r, truth = synthetic_panel()
    out = tables.fisher_at_threshold(r, "focal", truth, 2)
    odds, p = fisher_exact([[4, 2], [2, 8]], alternative="greater")
    assert out["testable"]
    assert (out["n_corrob"], out["n_uncorrob"]) == (6, 10)
    assert (out["truth_pos_corrob"], out["truth_pos_uncorrob"]) == (4, 2)
    assert out["OR"] == pytest.approx(odds)
    assert out["p"] == pytest.approx(p)


def test_untestable_when_no_candidate_reaches_threshold():
    r, truth = synthetic_panel()
    out = tables.fisher_at_threshold(r, "focal", truth, 4)
    assert out["testable"] is False
    assert out["p"] is None and out["OR"] is None


def test_bh_matches_manual_correction_and_skips_untestable():
    r, truth = synthetic_panel()
    res = tables.analyze_panel(r, truth)
    tests = [x for g in tables.GROUPS for x in res["raw"][g] if x["testable"]]
    assert res["n_bh_family"] == len(tests)
    # manual BH
    ps = sorted((x["p"], i) for i, x in enumerate(tests))
    m = len(ps)
    adj = [0.0] * m
    running = 1.0
    for rank in range(m, 0, -1):
        p, i = ps[rank - 1]
        running = min(running, p * m / rank)
        adj[i] = running
    for x, a in zip(tests, adj):
        assert x["p_adj_bh"] == pytest.approx(a)
    for g in tables.GROUPS:
        for x in res["raw"][g]:
            if not x["testable"]:
                assert x["p_adj_bh"] is None


def test_dedup_collapses_focal_but_not_negative():
    r = rows([("A", "focal", 2), ("A", "focal", 2), ("B", "focal", 0),
              ("N1", "negative", 2), ("N1", "negative", 2), ("N2", "negative", 0)])
    res = tables.analyze_panel(r, {"A"})
    assert res["n_focal_rows"] == 3 and res["n_focal_unique"] == 2
    assert res["n_negative_rows"] == 3
    assert res["focal_dedup"][1]["n_corrob"] == 1  # A counted once at >=2


def test_permutation_is_deterministic_for_fixed_seed():
    r, truth = synthetic_panel()
    d = cgc.dedup(r)
    a = cgc.permutation_test(d, "focal", truth, n=200, seed=42)
    b = cgc.permutation_test(d, "focal", truth, n=200, seed=42)
    assert a == b
    assert 0.0 <= a["empirical_p"] <= 1.0


def test_cgc_fisher_uses_truth_keys():
    r, truth = synthetic_panel()
    out = cgc.fisher_ge2(r, "focal", truth)
    assert "truth_pos_corrob" in out and "oncokb_pos_corrob" not in out
    assert out["truth_pos_corrob"] == 4


def test_incremental_scripts_load_with_generic_truth(tmp_path):
    data = [{"gene": "Q", "category": "focal", "scored_candidates": [
        {"gene": f"G{i}", "lincs_knockdown": i % 2 == 0, "super_enhancer": False,
         "dorothea_tf": i % 3 == 0, "depmap_essential": i % 4 == 0,
         "corroboration_count": i % 4} for i in range(6)]},
        {"gene": "H", "category": "housekeeping", "scored_candidates": [
            {"gene": "C1", "lincs_knockdown": False, "super_enhancer": False,
             "dorothea_tf": True, "depmap_essential": False, "corroboration_count": 1}]}]
    f = tmp_path / "panel.json"
    f.write_text(json.dumps(data), encoding="utf-8")
    truth = {"G0", "G3", "C1"}
    d_rows = dorothea.load_focal_dedup(f, truth)
    assert {r["gene"]: r["truth"] for r in d_rows}["G3"] == 1
    assert all("oncokb" not in r for r in d_rows)
    ctrl = dorothea.load_controls(f, truth)
    assert ctrl == [{"gene": "C1", "dorothea": 1, "truth": 1}]
    m_rows = mi_weight.load_focal_dedup(f, truth)
    assert [r["mi_top3"] for r in m_rows] == [1, 1, 1, 0, 0, 0]
    assert dorothea.stratified_test(d_rows) is None or "fisher_p" in dorothea.stratified_test(d_rows)


def test_output_suffix_keeps_legacy_filenames():
    assert ground_truth.output_suffix("oncokb") == ""
    assert ground_truth.output_suffix("intogen") == "_intogen"
    assert ground_truth.output_suffix("sanger_cgc", default="sanger_cgc") == ""


def test_unknown_truth_name_rejected():
    with pytest.raises(ValueError):
        ground_truth.load_truth("not_a_list")


def test_tampered_snapshot_is_rejected(tmp_path, monkeypatch):
    fake = tmp_path / "intogen_drivers.tsv"
    fake.write_text("symbol\trole\nTP53\toncogene\n", encoding="utf-8")
    monkeypatch.setattr(ground_truth, "INTOGEN_PATH", fake)
    with pytest.raises(ValueError, match="hash mismatch"):
        ground_truth.load_truth("intogen")


@pytest.mark.skipif(
    not (ground_truth.INTOGEN_PATH.exists() and ground_truth.VOGELSTEIN_PATH.exists()),
    reason="pinned snapshots live in the local RegNetAgents checkout and the private repo, "
           "so they are absent on CI",
)
def test_pinned_snapshots_load_with_recorded_counts():
    """Integrity only: confirms the pinned files are intact. Scores nothing."""
    intogen = ground_truth.load_truth("intogen")
    vogel = ground_truth.load_truth("vogelstein")
    assert len(intogen) == ground_truth.INTOGEN_N == 633
    assert len(vogel) == ground_truth.VOGELSTEIN_N == 125
    assert "TP53" in intogen and "TP53" in vogel


# ---- experiment_mi_rank_check.py --------------------------------------------------------
rank_check = importlib.import_module("experiment_mi_rank_check")


def _write_panel(tmp_path, name, queries):
    data = [{"gene": q["gene"], "category": q.get("category", "focal"),
             "scored_candidates": [{"gene": g, "corroboration_count": n} for g, n in q["cands"]]}
            for q in queries]
    f = tmp_path / name
    f.write_text(json.dumps(data), encoding="utf-8")
    return f


def test_rank_check_row_and_split_definitions(tmp_path):
    full = [(f"A{i}", i % 3) for i in range(10)]          # a full capped query
    short = [(f"B{i}", 0) for i in range(9)]               # a 9-candidate query
    neg = [(f"C{i}", 4) for i in range(10)]
    f = _write_panel(tmp_path, "p.json", [
        {"gene": "X", "cands": full}, {"gene": "Y", "cands": short},
        {"gene": "Z", "category": "housekeeping", "cands": neg}])
    rows = rank_check.load_focal_rows([f])
    assert len(rows) == 19                                  # negative controls excluded
    top, bottom = rank_check.split_top_bottom(rows)
    assert len(top) == 6                                    # 3 per query
    assert len(bottom) == 5                                 # 3 from the full query + 2 from the 9-long one
    assert {r["rank"] for r in bottom} == {8, 9, 10}


def test_rank_check_truth_association_matches_hand_computation(tmp_path):
    cands = [(f"G{i:02d}", 0) for i in range(10)]
    f = _write_panel(tmp_path, "p.json", [{"gene": "Q", "cands": cands}])
    rows = rank_check.load_focal_rows([f])
    res = rank_check.truth_association(rows, {"G00", "G01", "G02", "G09"})
    assert (res["truth_top"], res["n_top"]) == (3, 3)
    assert (res["truth_bottom"], res["n_bottom"]) == (1, 3)
    _, p = fisher_exact([[3, 0], [1, 2]], alternative="greater")
    assert res["fisher_p"] == pytest.approx(p)
    assert res["spearman_rho"] < 0                          # earlier rank, more likely positive


def test_rank_check_corroboration_part_needs_no_truth(tmp_path):
    cands = [(f"G{i:02d}", i % 4) for i in range(10)]
    f = _write_panel(tmp_path, "p.json", [{"gene": "Q", "cands": cands}])
    out = rank_check.corroboration_association(rank_check.load_focal_rows([f]))
    assert out["n_rows"] == 10 and out["n_top"] == 3 and out["n_bottom"] == 3
    assert out["mean_count_top"] == pytest.approx((0 + 1 + 2) / 3)
    assert out["mean_count_bottom"] == pytest.approx((3 + 0 + 1) / 3)
