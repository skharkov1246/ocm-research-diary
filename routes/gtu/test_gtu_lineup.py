"""Линейка KV: порядок по мощности на валу, каждое исполнение в одном типоразмере."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gtu_lineup as L  # noqa: E402


def _fam():
    return json.load(open(os.path.join(HERE, "gtu_family_results.json"), encoding="utf-8"))


def test_номера_по_росту_мощности_на_валу():
    lu = L.lineup(_fam())
    assert [x["name"] for x in lu] == [f"KV-{i}" for i in range(1, len(lu) + 1)]
    p = [x["P_shaft_MW_max"] for x in lu]
    assert p == sorted(p)


def test_все_исполнения_разложены_и_по_одному_разу():
    fam = _fam()
    lu = L.lineup(fam)
    закрыто = [s for s in fam["slots"] if s.get(L.SCENARIO)]
    разложено = [(m["model"], m["variant"]) for x in lu for m in x["members"]]
    assert len(разложено) == len(set(разложено)) == len(закрыто)


def test_типоразмер_3_это_группа_taurus_60():
    lu = L.lineup(_fam())
    kv3 = next(x for x in lu if x["name"] == "KV-3")
    assert "Taurus 60" in kv3["models"] and "Centaur 50" in kv3["models"]


def test_несовпадение_числа_исполнений_ловится():
    fam = _fam()
    fam["frames"][L.SCENARIO][0] = dict(fam["frames"][L.SCENARIO][0], n=99)
    try:
        L.lineup(fam)
    except ValueError:
        return
    raise AssertionError("расхождение с результатом не поймано")
