"""Проверки пространства сквозных концепций, оценки «КПД на тонну», MAP-Elites и файла задач AWS.

  python3 -m pytest -q test_process.py
"""
import gzip
import json
import math
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import process_space as PS  # noqa: E402
import process_eval as PE  # noqa: E402
import process_elites as PEL  # noqa: E402


# --------------------------------------------------------------------------- пространство
def test_decode_any_genome_is_valid_concept():
    X = PS.genomes(0, 256, "KV-3")
    for x in X:
        C = PS.decode(x, "KV-3")
        cell = PS.descriptors(C)
        assert all(0 <= c < b for c, b in zip(cell, PS.DESC_BINS))
        assert C["turbine"]["t41"] <= max(C["_t41_limit"], 850.0) + 1e-6
        lo, hi = (PS.T.FS().PI_RANGE[(C["compressor"]["arch"], C["compressor"]["n_ax"])] if C["compressor"]["arch"] == "ОЦ"
                  else PS.T.FS().PI_RANGE[C["compressor"]["arch"]])
        assert lo - 1e-6 <= C["compressor"]["pi"] <= hi + 1e-6
        if C["build"]["case"] == "оболочка":
            assert C["build"]["case_make"] == "печать"
        if C["turbine"]["cool"] == "C":
            assert C["turbine"]["type"] == "Р" and C["turbine"]["mat"] == "Si3N4"
    # граничные геномы
    for v in (0.0, 1.0):
        PS.decode(np.full(PS.D, v), "KV-3")


def test_generator_by_index_and_parts():
    a = PS.genomes(0, 40, "KV-3")
    b = np.vstack([PS.genomes(0, 17, "KV-3"), PS.genomes(17, 23, "KV-3")])
    assert np.allclose(a, b)
    assert np.allclose(PS.genome_at(29, "KV-3"), a[29])
    assert not np.allclose(PS.genomes(0, 5, "KV-1"), PS.genomes(0, 5, "KV-3"))
    n, P = 1003, 12
    seen = []
    for k in range(P):
        lo, hi = PS.part_range(n, P, k)
        seen += list(range(lo, hi))
    assert seen == list(range(n))


def test_encode_decode_roundtrip_reference():
    for name, x in PS.reference_seeds("KV-3"):
        C = PS.decode(x, "KV-3")
        x2 = PS.encode(C)
        C2 = PS.decode(x2, "KV-3")
        for blk, key in (("compressor", "arch"), ("turbine", "type"), ("turbine", "cool"), ("shaft", "drive"), ("combustor", "scheme"),
                         ("fuel", "prep"), ("hx", "reheat"), ("compressor", "intercool")):
            assert C[blk][key] == C2[blk][key], (name, blk, key)
        assert C["_make"] == C2["_make"]
        assert abs(C["compressor"]["pi"] / C2["compressor"]["pi"] - 1) < 1e-4
        assert abs(C["turbine"]["t41"] - C2["turbine"]["t41"]) < 0.2


def test_seeds_file_matches_space():
    d = json.load(open(PS.SEEDS_FILE, encoding="utf-8"))
    for kv in ("KV-3", "KV-1", "KV-6"):
        assert len(d["seeds"][kv]) > 0
        for r in d["seeds"][kv]:
            assert len(r["genome"]) == PS.D


# --------------------------------------------------------------------------- оценка
def _seed(kv, name):
    return dict(PS.reference_seeds(kv))[name]


def test_static_rejection_fast_and_explained():
    x = _seed("KV-3", "рекуператор 0,90 + подогрев, π 13,5 (как есть)").copy()
    x[PS.GI["turb"]] = 0.99            # радиальная турбина + подогрев на одном валу - жёсткое правило совместимости
    rec = PE.evaluate(x, "KV-3")
    assert not rec["feasible"] and rec["stage"] == "совместимость"
    assert "reheat_radial" in rec["reject"] and rec["reasons"][0]["text"]
    assert rec["cpu_s"] < 0.2


def test_feasible_record_and_kpd_per_tonne():
    rec = PE.evaluate(_seed("KV-3", "+ паровой риформинг (S/C 2)"), "KV-3", detail=True)
    assert rec["feasible"], rec["reject"]
    assert rec["mass_t"] > 1.0 and rec["mass_lo_t"] < rec["mass_t"] < rec["mass_hi_t"]
    assert math.isclose(rec["kpd_per_t"], rec["eff_el"] / rec["mass_t"], rel_tol=2e-3)
    assert math.isclose(rec["kW_per_kg"], rec["blocks"]["P_shaft_kW"] / (rec["mass_t"] * 1e3), rel_tol=2e-3)
    names = [r["name"] for r in rec["mass"]["rows"]]
    assert any("риформер" in n for n in names) and any("катализатор" in n for n in names)
    assert rec["fit_trl"] == pytest.approx(rec["kpd_per_t"] * rec["k_trl"], rel=1e-3)


def test_mass_bridge_matches_mass_model_example():
    """frame3 в литье в песок: мост блоки -> модель массы против примера модели массы (паспортный путь) - до 5 %"""
    rows = {r["name"]: r for r in PE.selfcheck()}
    cal = json.load(open(PE.MM.OUT_FILE, encoding="utf-8"))
    ex = {e["name"]: e for e in cal["examples"]["items"]}
    m_ex = ex["KV-3 простой; литьё в песок; масло пакета"]["total_kg"] / 1e3
    m = rows["frame3, традиционная"]["mass_t"]
    assert abs(m / m_ex - 1) < 0.05


def test_process_family_order_of_mass():
    rows = {r["name"]: r for r in PE.selfcheck()}
    assert rows["frame3, традиционная"]["mass_t"] > rows["frame3, литьё КНР"]["mass_t"] > rows["frame3, гибрид"]["mass_t"]


# --------------------------------------------------------------------------- архив
def _c(uid, eff, mass, feas=True, viol=0.0):
    return dict(uid=uid, g=[0.5] * PS.D, fit=(eff / mass if feas else None), eff=eff, mass=mass, feas=feas, viol=viol, parts=100, pnp=7,
                trl=9, ktrl=1.0, fit_trl=eff / mass, label=uid)


def test_archive_best_front_prune_and_merge_order():
    pts = [_c(f"u{i:02d}", 25 + i, 2.0 + 0.5 * i) for i in range(10)] + [_c("bad", 20, 9.0), _c("inf", 0, 1, feas=False, viol=0.3)]
    A = PEL.NicheArchive()
    for p in pts:
        A.add(p, (0, 0, 0, 0))
    e = A.cells[(0, 0, 0, 0)]
    assert e["best"]["uid"] == "u00"                     # 25 / 2,0 = 12,5 - наибольший КПД на тонну
    assert len(e["front"]) <= PEL.FRONT_CAP
    assert all(p["uid"] != "bad" for p in e["front"])
    assert e["front"][0]["uid"] == "u00" and e["front"][-1]["uid"] == "u09"
    assert e["infeas"]["uid"] == "inf"
    B1, B2 = PEL.NicheArchive(), PEL.NicheArchive()
    for p in pts[:6]:
        B1.add(p, (1, 0, 0, 0))
    for p in pts[6:]:
        B2.add(p, (1, 0, 0, 0))
    M1, M2 = PEL.NicheArchive(), PEL.NicheArchive()
    M1.merge(B1.cells)
    M1.merge(B2.cells)
    M2.merge(B2.cells)
    M2.merge(B1.cells)
    assert json.dumps(M1.to_list(), sort_keys=True) == json.dumps(M2.to_list(), sort_keys=True)


def test_island_determinism_and_merge_completeness(tmp_path):
    outs = []
    for w in (1, 2):
        d = tmp_path / f"w{w}"
        PEL.run_part("KV-3", 2, 2, 0, 60, 0.3, 30, PS.SEED, w, str(d), log=None)
        outs.append(gzip.open(PEL.part_file(str(d), "KV-3", 0), "rt", encoding="utf-8").read())
    a, b = (json.loads(o)["archive"] for o in outs)
    assert a == b
    # слияние: нет острова 1 -> отказ; разрешено явно -> отчёт «неполное»
    paths = [PEL.part_file(str(tmp_path / "w1"), "KV-3", 0)]
    with pytest.raises(SystemExit):
        PEL.merge_parts(paths, "KV-3", 2, 60, init_frac=0.3, batch=30)
    _, _, rep, _ = PEL.merge_parts(paths, "KV-3", 2, 60, allow_partial=True, init_frac=0.3, batch=30)
    assert not rep["complete"] and rep["missing"] == [1]
    PEL.run_part("KV-3", 2, 2, 1, 60, 0.3, 30, PS.SEED, 1, str(tmp_path / "w1"), log=None)
    paths.append(PEL.part_file(str(tmp_path / "w1"), "KV-3", 1))
    arch, stats, rep, _ = PEL.merge_parts(paths, "KV-3", 2, 60, init_frac=0.3, batch=30)
    assert rep["complete"] and rep["evals_total"] == 120
    # другая конфигурация (оценок острова) - чужая
    _, _, rep2, _ = PEL.merge_parts(paths, "KV-3", 2, 61, allow_partial=True, init_frac=0.3, batch=30)
    assert not rep2["complete"] and len(rep2["config_mismatch_files"]) == 2


# --------------------------------------------------------------------------- задачи AWS
def test_jobs_file_valid_for_runner():
    path = os.path.join(HERE, "jobs_process.json")
    if not os.path.exists(path):
        pytest.skip("jobs_process.json ещё не собран")
    d = json.load(open(path, encoding="utf-8"))
    sys.path.insert(0, os.path.join(HERE, "..", "platform", "aws"))
    import runner
    jobs = runner.load_jobs(path)
    assert len(jobs) >= 10
    for kv in ("KV-3", "KV-1", "KV-6"):
        assert all(f"--kv {kv}" in j.command for j in jobs)
    assert d["plan_check"]["exit_code"] == 0
    assert d["estimate"]["concepts_by_kv"]["KV-3"] >= 1_000_000
    assert d["estimate"]["upper_bound_usd"] < 2000
    parts = sorted(int(j.command.split("--part ")[1].split()[0]) for j in jobs)
    assert parts == list(range(len(jobs)))
