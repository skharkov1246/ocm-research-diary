"""Тесты сверки повторного расчёта (compare_results.py) на придуманных JSON (SYNTHETIC)."""
import json

import compare_results as CR


def test_identical_is_ok_and_service_keys_ignored():
    a = {"date": "2026-11-14", "software": "x", "summary": {"n": 3, "p": [1.0, 2.0]}, "runtime_s": 10}
    b = {"date": "2099-01-01", "software": "y", "summary": {"n": 3, "p": [1.0, 2.0]}, "runtime_s": 99}
    st = CR.compare(a, b)
    assert st["ok"] and st["numbers"] == 3 and st["numbers_bad"] == 0


def test_number_out_of_tolerance_and_structure_change_fail():
    ref = {"r": [{"P": 10.0, "eta": 0.3}], "k": 1}
    st = CR.compare({"r": [{"P": 10.0 * (1 + 1e-4), "eta": 0.3}], "k": 1}, ref)
    assert not st["ok"] and st["numbers_bad"] == 1 and st["bad"][0]["path"] == "/r[0]/P"
    assert CR.compare({"r": [{"P": 10.0 * (1 + 1e-8), "eta": 0.3}], "k": 1}, ref)["ok"]  # в допуске 1e-6
    st = CR.compare({"r": [{"P": 10.0}], "k": 1}, ref)
    assert not st["ok"] and st["missing_in_new"] == ["/r[0]/eta"]
    st = CR.compare({"r": [], "k": 1}, ref)
    assert not st["ok"] and st["length_mismatch"]


def test_nan_matches_nan_and_strings_do_not_fail():
    st = CR.compare({"x": float("nan"), "s": "a"}, {"x": float("nan"), "s": "b"})
    assert st["ok"] and st["strings_diff"] == 1


def test_collect_summarizes_jobs(tmp_path):
    d = tmp_path / "val-cycle"
    d.mkdir()
    json.dump(dict(label="balances", ok=True, numbers=5, numbers_bad=0), open(d / "repro_balances.json", "w"))
    (d / "timing.txt").write_text("balances_s 12\n")
    (d / "pytest.txt").write_text("....\n18 passed in 30.0s\n")
    (tmp_path / "_lost").mkdir()
    doc = CR.collect(str(tmp_path), "SYNTHETIC-run")
    assert doc["ok"] and list(doc["jobs"]) == ["val-cycle"]
    assert doc["jobs"]["val-cycle"]["pytest"]["pytest.txt"][-1] == "18 passed in 30.0s"
    assert CR.main(["--collect", str(tmp_path / "nothing")]) == 2
