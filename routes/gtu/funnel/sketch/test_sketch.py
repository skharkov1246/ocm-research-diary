"""Тест генератора эскизов (funnel_sketch.py).

Проверки (решение владельца 14.11.2026, эскизы к каждому из 50):
  1. габарит модели совпадает с габаритом JSON (geometry.envelope) в пределах 1 % - у каждого варианта funnel_l1_top50.json и у
     каждой фикстуры, кроме выносной камеры, где несогласованность компоновки уровня 1 должна быть записана в issues (код
     ext_comb_rule), а не спрятана; длина L сверяется отдельно по цепочке узлов, без габарита с обвязкой;
  2. у каждого варианта три картинки (вид сбоку, с торца, аксонометрия) не меньше 1600×1000;
  3. повтор даёт те же размеры (модель и sketch_dims);
  4. двусторонний вход рисует два входа, рекуператор - рекуператор (в корпусе - барабан, сбоку - все модули), кластер - нужное
     число ядер; условные узлы и модули без размеров не выходят за габарит.
Фикстуры - fixtures_l1_layout.json (make_fixtures.py: настоящие концепции уровня 0 через funnel_l1.layout).
Рендер (pyvista) - только если пакет есть: python -m pytest routes/gtu/funnel/sketch/test_sketch.py
"""
import json, os, sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import funnel_sketch as SK   # noqa: E402

# SKETCH_TOP - итог слияния AWS (aws/l1_full/funnel_l1_top50_KV-3_m<M>.json, RUNBOOK_L0_L1.md, раздел 7а); по умолчанию - проба
TOP = os.environ.get("SKETCH_TOP") or os.path.join(os.path.dirname(HERE), "funnel_l1_top50.json")
FIX = os.path.join(HERE, "fixtures_l1_layout.json")


def _entries(path):
    return SK.load_entries(path)[1]


def _all():
    out = [("top", e) for e in _entries(TOP)]
    if os.path.exists(FIX):
        out += [("fix", e) for e in _entries(FIX)]
    return out


def _id(item):
    kind, e = item
    return f"{kind}-{e.get('div_rank') or e.get('fixture_class')}"


ALL = _all()


@pytest.mark.parametrize("item", ALL, ids=[_id(i) for i in ALL])
def test_envelope_within_1pct(item):
    kind, e = item
    M = SK.build_model(e)
    d = SK.dims_of(M)
    envj = e["geometry"]["envelope"]
    codes = {it["code"] for it in M.issues}
    assert "fit_outside" not in codes, M.issues
    # длина - по цепочке узлов (вход ... выхлоп + муфта), без габарита с обвязкой
    nlo, nhi = SK.nodes_bbox(M, include_fit=False)
    assert abs((nhi[0] - nlo[0]) / envj["L"] - 1) <= SK.TOL_ENV
    # габарит с обвязкой, собранный из размеров узлов по правилу компоновки, - с полем envelope
    for k in "LWH":
        assert abs(d["dev_rule_pct"][k]) <= 100 * SK.TOL_ENV, (k, d)
    if kind == "top" or "ext_comb_rule" not in codes:
        for k in "LWH":
            assert abs(d["dev_pct"][k]) <= 100 * SK.TOL_ENV, (k, d, M.issues)
        assert d["ok"]
    else:                       # выносная камера: модель выходит за габарит, причина записана, а не подрезана
        assert not d["ok"]
        assert any(it["code"] == "protrudes" and it["node"].startswith("выносная") for it in M.issues)


@pytest.mark.parametrize("item", ALL, ids=[_id(i) for i in ALL])
def test_repeat_same_dims(item):
    _, e = item
    a, b = SK.dims_of(SK.build_model(e)), SK.dims_of(SK.build_model(e))
    assert a == b
    ra = SK.dims_record(SK.build_model(e))
    rb = SK.dims_record(SK.build_model(e))
    assert json.dumps(ra, sort_keys=True, ensure_ascii=False) == json.dumps(rb, sort_keys=True, ensure_ascii=False)


def test_run_repeat_same_dims_json(tmp_path):
    r1, _ = SK.run(TOP, str(tmp_path / "a"), str(tmp_path / "d1.json"), None, render=False)
    r2, _ = SK.run(TOP, str(tmp_path / "b"), str(tmp_path / "d2.json"), None, render=False)
    v1 = json.load(open(tmp_path / "d1.json", encoding="utf-8"))
    v2 = json.load(open(tmp_path / "d2.json", encoding="utf-8"))
    for v in (v1, v2):
        for r in v["variants"]:
            r.pop("sheets", None)
    assert v1["variants"] == v2["variants"] and v1["scale_px_per_m"] == v2["scale_px_per_m"]
    assert r1["n_drawn"] == len(_entries(TOP)) and not r1["skipped"]


def _fixture(cls):
    for kind, e in ALL:
        if kind == "fix" and e.get("fixture_class") == cls:
            return e
    pytest.skip(f"нет фикстуры {cls}")


def _double_inlet_entries():
    return [e for _, e in ALL if str(e["geometry"]["inlet"]["type"]).startswith("двусторонний")]


def test_double_inlet_two_inlets():
    es = _double_inlet_entries()
    assert es, "нет варианта с двусторонним входом"
    for e in es:
        M = SK.build_model(e)
        n = max(1, int((e["geometry"].get("cores") or {}).get("n", 1) or 1))
        chambers = [p for p in M.parts if p["name"].startswith("входная камера") and p["kind"] == "solid"]
        assert len(chambers) == 2 * n
        letters = {c["letter"] for c in M.conns}
        assert {"B1", "B2"} <= letters
        single = [p for p in M.parts if p["name"] == "входной коллектор"]
        assert not single


def test_single_inlet_one_inlet():
    e = next(e for _, e in ALL if not str(e["geometry"]["inlet"]["type"]).startswith("двусторонний"))
    M = SK.build_model(e)
    n = max(1, int((e["geometry"].get("cores") or {}).get("n", 1) or 1))
    assert len([p for p in M.parts if p["name"] == "входной коллектор"]) == n
    assert not [p for p in M.parts if p["name"].startswith("входная камера")]


def test_recuperator_side_all_modules():
    e = _fixture("recup_side")
    M = SK.build_model(e)
    rec = e["geometry"]["recuperator"]
    cubes = [p for p in M.parts if p["name"].startswith("модуль рекуператора") and p["kind"] == "solid"]
    assert len(cubes) == rec["modules"]
    side = rec["module_side"]
    for p in cubes:
        assert np.allclose(p["hi"] - p["lo"], side, rtol=1e-9)
    no = next(m["no"] for m in M.modules if m["name"] == "модуль рекуператора")
    assert all(p["module"] == no for p in cubes)
    # блок - сбоку от двигателя, на стороне +y (за двигателем на виде сбоку, смотрим с -y): целиком по одну сторону от корпуса
    assert min(p["lo"][1] for p in cubes) >= 0.5 * M.rule["D_core"] - 1e-9
    assert M.env["lo"][1] == pytest.approx(-0.5 * M.rule["D_core"])


def test_recuperator_integral_drum():
    e = _fixture("recup_integral")
    M = SK.build_model(e)
    rec = e["geometry"]["recuperator"]
    drums = [p for p in M.parts if p["name"].startswith("рекуператор в корпусе")]
    assert len(drums) == 1
    p = drums[0]
    assert abs((p["hi"][1] - p["lo"][1]) / rec["D_outer"] - 1) < 1e-3
    assert abs((p["hi"][0] - p["lo"][0]) / rec["L"] - 1) < 1e-9


def test_no_recuperator_without_one():
    for _, e in ALL:
        if e["geometry"].get("recuperator") is None:
            M = SK.build_model(e)
            assert not [p for p in M.parts if "рекуператор" in p["name"]]


@pytest.mark.parametrize("cls", ["cluster2", "cluster4"])
def test_cluster_core_count(cls):
    e = _fixture(cls)
    n = e["geometry"]["cores"]["n"]
    M = SK.build_model(e)
    assert M.cores == n
    assert len([b for b in M.bodies if b["name"].startswith("ядро")]) == n
    assert len([p for p in M.parts if p["name"].startswith("турбина газогенератора")]) == n
    assert len([p for p in M.parts if p["name"].startswith("силовая турбина")]) == 1
    # ядра не накладываются друг на друга: оси разнесены не меньше чем на диаметр корпуса
    cs = [np.array([b["y"], b["z"]]) for b in M.bodies if b["name"].startswith("ядро")]
    dmin = min(np.linalg.norm(a - b) for i, a in enumerate(cs) for b in cs[i + 1:])
    assert dmin >= M.rule["D_case"] - 1e-9


def test_real_top50_clusters_and_inlets():
    for e in _entries(TOP):
        M = SK.build_model(e)
        n = int((e["geometry"].get("cores") or {}).get("n", 1) or 1)
        assert M.cores == n
        assert len([p for p in M.parts if p["name"].startswith("турбина газогенератора")]) == n


def test_conditional_marked_and_listed():
    for _, e in ALL:
        M = SK.build_model(e)
        rec = SK.dims_record(M)
        assert "опоры" in " ".join(rec["conditional_notes"]) or not e["geometry"].get("supports")
        assert rec["conditional_nodes"], "у каждого варианта есть условные узлы (фланцы, опоры)"
        assert set(rec["solid_nodes"]).isdisjoint(rec["conditional_nodes"])
        for m in rec["modules"]:
            if m["name"] == "контроллер магнитных подшипников":
                assert not m["drawn"]          # места и размера нет - не рисуется


# ------------------------------------------------------------------ рендер (нужен pyvista)
def _png_size(path):
    with open(path, "rb") as f:
        head = f.read(24)
    assert head[:8] == b"\x89PNG\r\n\x1a\n"
    return int.from_bytes(head[16:20], "big"), int.from_bytes(head[20:24], "big")


def test_three_pictures_each(tmp_path):
    pytest.importorskip("pyvista")
    pytest.importorskip("matplotlib")
    for path in (TOP, FIX):
        if not os.path.exists(path):
            continue
        out = tmp_path / os.path.basename(path).split(".")[0]
        res, models = SK.run(path, str(out), str(tmp_path / "dims.json"), None, render=True, workers=min(4, os.cpu_count() or 1))
        assert res["rendered"] == 3 * len(models) == 3 * len(_entries(path))
        for M in models:
            paths = SK.sheet_paths(str(out), M.entry)
            assert len(paths) == 3
            for v, p in paths.items():
                assert os.path.exists(p), p
                w, h = _png_size(p)
                assert w >= 1600 and h >= 1000, (p, w, h)


def test_render_independent_of_workers(tmp_path):
    pytest.importorskip("pyvista")
    e = _entries(TOP)[0]
    M = SK.build_model(e)
    s = SK.common_scale([M])
    a = SK.render_sheet((M, "side", s, str(tmp_path / "a.png"), "KV", 1, "x"))
    import multiprocessing as mp
    with mp.get_context("spawn").Pool(1) as pool:
        b = pool.map(SK.render_sheet, [(M, "side", s, str(tmp_path / "b.png"), "KV", 1, "x")])[0]
    assert open(a, "rb").read() == open(b, "rb").read()


def test_cluster_external_combustor_collision_reported():
    """производный случай (не из выборки): кластер 2×2 фикстуры cluster4 с выносной камерой фикстуры cluster2 - камера нижнего
    ядра входит в верхнее ядро; эскиз обязан это записать (issues, код intersect), а не подрезать"""
    import copy
    e4, e2 = copy.deepcopy(_fixture("cluster4")), _fixture("cluster2")
    if e2["geometry"]["combustor"]["type"] != "выносная":
        pytest.skip("у фикстуры cluster2 камера не выносная")
    cb = copy.deepcopy(e2["geometry"]["combustor"])
    g4 = e4["geometry"]
    cb["x0"] = g4["combustor"]["x0"]
    g4["combustor"] = cb
    M = SK.build_model(e4)
    assert any(it["code"] == "intersect" for it in M.issues), M.issues
    assert any(it["code"] == "ext_comb_rule" for it in M.issues)


# ------------------------------------------------------------------ правки по скептикам (15.11.2026)
def _legend_text(M, view="side"):
    return [ln[1] for ln in SK.legend_lines(M, view, 300, "KV")]


def test_each_connection_point_own_letter():
    """масло у передней и задней опоры - E1 и E2 (а не одна E и «2 шт.»), у каждой точки - строка легенды"""
    seen_oil = False
    for _, e in ALL:
        M = SK.build_model(e)
        oil = [cn for cn in e["geometry"]["connections"] if cn["kind"].startswith("масло")]
        letters = [c["letter"] for c in M.conns if c["letter"].startswith("E")]
        if len(oil) > 1:
            seen_oil = True
            assert letters == [f"E{k + 1}" for k in range(len(oil))], letters
            leg = _legend_text(M)
            for k, cn in enumerate(oil):
                ln = next(t for t in leg if t.startswith(f"E{k + 1} — "))
                assert cn["where"] in ln and "шт." not in ln
        # буквы точек разные, кроме повтора B/D на ядрах кластера
        per_core = [c["letter"] for c in M.conns if c["inst"] == 0]
        assert len(per_core) == len(set(per_core)), per_core
    assert seen_oil, "нет варианта с маслом у двух опор"


def test_two_air_points_numbered():
    """производный случай: перепуск на подогрев входа (funnel_l1, inlet_heat) - вторая точка «воздух» получает G2"""
    import copy
    e = copy.deepcopy(_entries(TOP)[0])
    g = e["geometry"]
    g["connections"].append(dict(kind="воздух: перепуск на подогрев входа", where="за компрессором -> вход", x=g["connections"][-1]["x"], n=1))
    M = SK.build_model(e)
    assert [c["letter"] for c in M.conns if c["letter"].startswith("G")] == ["G1", "G2"]
    leg = _legend_text(M)
    assert any(t.startswith("G1 — ") for t in leg) and any(t.startswith("G2 — ") for t in leg)


def test_double_inlet_wheel_between_chambers():
    for e in _double_inlet_entries():
        M = SK.build_model(e)
        n = max(1, int((e["geometry"].get("cores") or {}).get("n", 1) or 1))
        wheels = [p for p in M.parts if p["name"].startswith("колесо 1-й ступени между камерами")]
        assert len(wheels) == n and all(p["kind"] == "cond" for p in wheels)
        for w in wheels:
            ch = sorted([p for p in M.parts if p["name"].startswith("входная камера") and p["inst"] == w["inst"]], key=lambda p: p["lo"][0])
            assert ch[0]["hi"][0] <= w["lo"][0] + 1e-9 and w["hi"][0] <= ch[1]["lo"][0] + 1e-9   # колесо - между камерами
        assert any(it["code"] == "double_inlet_order" for it in M.issues)
        # площадь фланца: у B1 и B2 - каждого, а не сумма
        A = e["geometry"]["inlet"]["flange_area_m2"]
        ln = next(t for t in _legend_text(M) if t.startswith("B1 — "))
        assert f"{SK.fnum(A / (2 * n), 2)} м² каждый" in ln and f"всего {SK.fnum(A, 2)} м²" in ln


def test_recuperator_integral_letters_and_enclosed():
    e = _fixture("recup_integral")
    M = SK.build_model(e)
    drum = next(p for p in M.parts if p["name"].startswith("рекуператор в корпусе"))
    x0, x1, R = drum["lo"][0], drum["hi"][0], 0.5 * (drum["hi"][2] - drum["lo"][2])
    for letter in ("C", "D"):
        cn = next(c for c in M.conns if c["letter"] == letter)
        p = cn["pts"][0]
        assert not (x0 + 1e-6 < p[0] < x1 - 1e-6 and abs(p[2] - M.z_ax) < R - 1e-6), (letter, p)   # точка не внутри барабана
    names = {m["no"]: m["name"] for m in M.modules}
    assert set(M.enclosed) == {no for no, nm in names.items() if nm in ("картридж ротора", "горячая часть", "камера")}
    assert SK.hidden_outlines(M)
    leg = _legend_text(M)
    for no in M.enclosed:
        assert "скрытым контуром" in next(t for t in leg if t.startswith(f"{no} — "))
    assert any(it["code"] == "recup_overlap" for it in M.issues)


def test_cold_end_axial_inlet_flange_not_on_shaft_axis():
    """привод с холодного конца: торец входа занят валом - фланец B сверху коллектора, не диск на оси вала"""
    hit = False
    for _, e in ALL:
        g = e["geometry"]
        if str(g["inlet"]["type"]).startswith("двусторонний") or not str(g["output_shaft"]["output_end"]).startswith("холодный"):
            continue
        hit = True
        M = SK.build_model(e)
        fl = [p for p in M.parts if p["name"].startswith("фланец входа B")]
        assert fl and all(p["prim"][0] == "box" for p in fl)
        for p in fl:
            assert p["lo"][2] >= M.z_ax + 0.5 * g["inlet"]["D"] - 1e-9
    assert hit, "нет варианта: осевой вход и привод с холодного конца"


def test_forced_scale_checked(tmp_path):
    s_ok = SK.run(FIX, str(tmp_path / "a"), None, None, render=False)[0]["scale_px_per_m"]
    with pytest.raises(SK.ScaleOverflow):
        SK.run(FIX, str(tmp_path / "b"), None, None, render=False, px_per_m=s_ok + 100)
    res, models = SK.run(FIX, str(tmp_path / "c"), None, None, render=False, px_per_m=s_ok + 100, allow_overflow=True)
    assert res["scale_overflow"] and any(M.scale_note for M in models)
    res, models = SK.run(FIX, str(tmp_path / "d"), None, None, render=False, px_per_m=s_ok)
    assert not res["scale_overflow"] and not any(M.scale_note for M in models)


def test_axo_scale_label_along_axes():
    M = SK.build_model(_entries(TOP)[0])
    k = SK.axo_axis_scale(M, 380)
    assert all(v < 380 for v in k)
    t = SK.scale_text(M, "axo", 380)
    for v in k:
        assert str(int(round(v))) in t
    assert SK.scale_text(M, "side", 380) == "1 м = 380 пикс."


def test_wheels_drawn_match_geometry():
    for _, e in ALL:
        M = SK.build_model(e)
        rec = SK.dims_record(M)
        for w, st in zip(rec["compressor_wheels_drawn_m"], e["geometry"]["compressor"]):
            keys = ("D",) if st["type"] == "осевая" else ("D1t", "D2", "D_case")
            for k in keys:
                assert w[k] == pytest.approx(st[k], abs=1e-4), (k, w, st)
            assert w["x0"] == pytest.approx(st["x0"], abs=1e-4)


def test_recuperator_integral_labels_visible():
    """рендер: у рекуператора в корпусе C и D видны хоть на одном листе, модули 1-3 - выносками к скрытому контуру на виде сбоку"""
    pytest.importorskip("pyvista")
    M = SK.build_model(_fixture("recup_integral"))
    s = SK.common_scale([M])
    texts = {}
    for v in SK.VIEWS:
        img, depth, P, cam = SK.render_image(M, v, s)
        texts[v] = [lb["text"] for lb in SK.anchors(M, P, cam, depth, v)]
    allt = ", ".join(t for v in texts for t in texts[v])
    for letter in ("A", "B", "C", "D"):
        assert letter in allt.replace(" ", "").split(","), texts
    for no in M.enclosed:
        assert str(no) in texts["side"], texts
