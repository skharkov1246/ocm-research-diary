"""Проверки обоснования целей камеры (targets_basis_results.json) и файла целей (targets.json).
Запуск: python3 -m pytest -q test_targets_basis.py   (без пересчёта, несколько секунд)

Проверяется: метаданные расчёта, повтор расчётной точки KV-3, физические тренды (NOx растёт с phi и временем,
срыв отступает с ростом времени пребывания, потери давления снижают мощность), формула допуска PF, и что каждое
число обоснования в targets.json совпадает со своим ключом в results-JSON."""
import json
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "targets_basis_results.json")))
T = json.load(open(os.path.join(HERE, "targets.json")))


def get(path, d=R):
    cur = d
    for part in re.findall(r"[^.\[\]]+(?:\[[^\]]+\])?", path):
        name, f, v = re.match(r"([^\[]+)(?:\[([^=]+)=([^\]]+)\])?$", part).groups()
        cur = cur[name]
        if f:
            cur = next(x for x in cur if abs(float(x[f]) - float(v)) < 1e-9)
    return cur


def test_metadata():
    for k in ("method", "software", "date", "inputs", "assumptions", "status"):
        assert R[k], k
    assert R["date"] == "2026-11-14" and T["date"] == "2026-11-14"
    assert all("file" in i and "field" in i for i in R["inputs"])


def test_cycle_repro():
    assert abs(R["cycle_check"]["rel_diff"]) < 1e-3
    i = R["inlet"]
    assert math.isclose(i["FAR"], i["W_fuel_kg_s"] / i["W_air_comb_kg_s"], rel_tol=1e-5)
    assert math.isclose(i["phi_overall"], i["FAR"] / i["FAR_stoich"], rel_tol=1e-5)


def test_nox_trends():
    lit = [r for r in R["crn_grid"] if r["lit"]]
    nox10 = [get("hot[tau_ms=10].NOx_ppm15", r) for r in lit]
    assert nox10 == sorted(nox10), "NOx должен расти с phi"
    for r in lit:
        n = [h["NOx_ppm15"] for h in r["hot"]]
        assert n == sorted(n), "NOx должен расти со временем пребывания"
    u = [get("by_tau[tau_ms=10].NOx_ppm15", x) for x in R["unmixedness"]]
    assert u == sorted(u), "NOx должен расти с неравномерностью смешения"


def test_lbo_and_part_load():
    lbo = sorted(R["psr_lbo"], key=lambda x: x["tau_ms"])
    phis = [x["phi_lbo"] for x in lbo]
    assert phis == sorted(phis, reverse=True), "с ростом времени пребывания срыв отступает в бедную сторону"
    pl = R["part_load_fixed_split"]
    assert pl[0]["load_frac"] == 1.0
    fr = [x["load_frac"] for x in pl]
    assert fr == sorted(fr, reverse=True)


def test_dp_and_pf():
    rows = R["dp_sensitivity"]["rows"]
    p = [x["P_shaft_kW"] for x in rows]
    assert p == sorted(p, reverse=True)
    e = R["exit_uniformity"]
    assert math.isclose(e["pattern_factor_max"], (e["T_ax_max_C"] - e["T41_C"]) / (e["T41_C"] - e["T3_C"]), rel_tol=1e-5)


def test_targets_reference_results():
    n = 0
    for blk in T["targets"] + T["guides"]:
        for c in blk.get("basis_calc", []):
            assert get(c["key"]) == c["value"], c["key"]
            n += 1
    assert n >= 30
    ids = [t["id"] for t in T["targets"]]
    assert len(ids) == len(set(ids))
    for t in T["targets"]:
        for b in t.get("best_in_class", []):
            assert b["status"] in ("открыт", "по выдаче поиска"), b
    pf = T["targets"][ids.index("exit_uniformity")]["target"]["PF_max"]
    assert pf <= R["exit_uniformity"]["pattern_factor_max"] < pf + 1e-3, "цель PF — расчётный предел с округлением вниз"
