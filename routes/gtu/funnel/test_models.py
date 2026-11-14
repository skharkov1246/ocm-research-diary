"""Проверки узла воронки (уровень 0, часть 1): генератор, разбиение на части, правила допустимости, четыре модели.
Корпуса проверок - придуманные ведомости и модули (не данные установок). Запуск: python3 -m pytest -q test_models.py"""
import copy, json, math, os, sys, tempfile

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import funnel_space as F       # noqa: E402
import make_model as MK        # noqa: E402
import parts_model as PM       # noqa: E402
import plugplay_model as PP    # noqa: E402
import service_model as SM     # noqa: E402


# =========================================================================== генератор и части
def test_concept_repeatable_and_matrix_deterministic():
    a = F.concept(17, 300, "KV-3")
    fp = F.fingerprint(300)
    F._DM.clear()
    assert F.concept(17, 300, "KV-3") == a
    assert F.fingerprint(300) == fp
    assert F.fingerprint(300, seed=F.SEED + 1) != fp


def test_shards_partition_without_gaps_or_overlaps():
    n, S = 1001, 7
    parts = [F.shard_indices(n, k, S) for k in range(S)]
    flat = [i for p in parts for i in p]
    assert sorted(flat) == list(range(n)) and len(flat) == len(set(flat))
    for bad in ((n, S, S), (n, -1, S), (n, 0, 0), (0, 0, 1)):
        with pytest.raises(ValueError):
            F.shard_indices(*bad)


def test_record_independent_of_order():
    n = 64
    fwd = [F.evaluate(i, n, "KV-1") for i in range(8)]
    rev = [F.evaluate(i, n, "KV-1") for i in reversed(range(8))][::-1]
    assert json.dumps(fwd, sort_keys=True, default=float) == json.dumps(rev, sort_keys=True, default=float)


def test_all_archetypes_present_and_templates_respected():
    n = 256
    seen = set()
    for i in range(n):
        c = F.concept(i, n)
        seen.add(c["arch"])
        T = F.ARCH[c["arch"]]["allowed"]
        for ax, vals in T.items():
            if ax in ("mat_cw",) and c["rotor_join"] in ("заодно", "сварка"):
                continue                       # колесо заодно с турбиной - никелевый сплав
            if ax == "rotors" and c["cores"] > 1:
                continue
            if ax == "m_hot" and c["cool"] == "P":
                continue
            if ax in ("npt", "n_ax", "recup_integral"):
                continue
            assert c[ax] in vals, (c["arch"], ax, c[ax], vals)
        if c["cores"] > 1:
            assert c["rotors"] == 2
        if c["rotor_join"] in ("заодно", "сварка"):
            assert c["mat_cw"] == "Ni" and c["m_cw"] == c["m_tw"]
        if c["mat_tw"] == "SC":
            assert c["turb"] == "О"
        if c["mat_tw"] == "Si3N4":
            assert c["turb"] == "Р"
    assert seen == set(F.ARCH_KEYS)


def test_continuous_ranges():
    n = 200
    for i in range(n):
        c = F.concept(i, n)
        lo, hi = F.PI_RANGE[(c["comp"], c["n_ax"])] if c["comp"] == "ОЦ" else F.PI_RANGE[c["comp"]]
        assert lo - 1e-9 <= c["pi"] <= hi + 1e-9
        assert c["pi"] <= max(c["pi_max_U"], lo) + 1e-6
        tl = MK.t41_limit(c["turb"], c["cool"], c["mat_tw"], c["m_tw"])[0]
        assert F.T41_MIN - 1e-9 <= c["t41"] <= max(tl, F.T41_MIN) + 1e-6
        assert 0.35 <= c["split"] <= 0.65 and 0.75 <= c["kUt"] <= 1.0 and 0.45 <= c["phit"] <= 0.80
        assert (0.85 <= c["kUpt"] <= 1.0) if c["rotors"] == 2 else c["kUpt"] == 0.0
        assert (0.75 <= c["eps_rec"] <= 0.92) if c["recup"] == "печатный" else c["eps_rec"] == 0.0


def test_run_shard_file(tmp_path):
    path, doc = F.run_shard("KV-1", 20, 1, 3, str(tmp_path))
    d = json.load(open(path, encoding="utf-8"))
    assert d["indices"] == list(range(1, 20, 3)) == [r["idx"] for r in d["records"]]
    assert d["fingerprint"] == F.fingerprint(20) and d["date"] == "2026-11-14"
    for k in ("method", "software", "inputs", "target"):
        assert d[k]


# =========================================================================== быстрый цикл: сверка двумя способами
def test_core_matches_gtu_cycle():
    for r in F.cross_check():
        assert abs(r["d_eta_pp"]) <= 0.2, r
        assert abs(r["w_core"] - r["w_cycle"]) <= 2.0, r
        assert abs(r["EGT_core"] - r["EGT_cycle"]) <= 2.0, r


# =========================================================================== правила допустимости
@pytest.fixture(scope="module")
def feasible_case():
    n = 256
    for i in range(n):
        r = F.evaluate(i, n, "KV-1")
        if r["feasible"]:
            c = F.concept(i, n, "KV-1"); t = F.target("KV-1")
            p1 = F.quick_perf(c, t); g1 = F.derive(c, t, p1)
            p = F.quick_perf(c, t, g1["d_eta"]); g = F.derive(c, t, p)
            assert all(x["ok"] for x in F.check(c, t, p, g))
            return c, t, p, g
    pytest.skip("в отладочной выборке нет допустимой концепции KV-1")


def _codes(rules):
    return {r["code"].rstrip("0123456789") if not r["code"].startswith("size:") else "size" for r in rules if not r["ok"]}


def test_rule_tip_speed_by_material(feasible_case):
    c, t, p, g = copy.deepcopy(feasible_case)
    lim = MK.u_limit_comp(g["stages"][0]["mat"], c["m_cw"])[0]
    g["stages"][0]["U"] = lim * 1.05
    assert "U_comp" in _codes(F.check(c, t, p, g))


def test_rule_t41_by_cooling_class(feasible_case):
    c, t, p, g = copy.deepcopy(feasible_case)
    c["t41"] = MK.t41_limit(c["turb"], c["cool"], c["mat_tw"], c["m_tw"])[0] + 10.0
    assert "T41_class" in _codes(F.check(c, t, p, g))


def test_rule_soft_tolerance(feasible_case):
    c, t, p, g = copy.deepcopy(feasible_case)
    p["G"] = t["G0"] * 1.21
    p["EGT"] = t["EGT0"] + 101.0
    codes = _codes(F.check(c, t, p, g))
    assert {"soft_G", "soft_EGT"} <= codes


def test_rule_part_size_printer_and_casting(feasible_case):
    c, t, p, g = copy.deepcopy(feasible_case)
    g["parts_mk"] = g["parts_mk"] + [dict(name="ротор целиком", kind="rotor_mono", method="печать", material="IN718",
                                          D=2.5, L=1.0, mass_kg=100.0, rotating=True),
                                     dict(name="колесо турбины", kind="turb_radial", method="литьё", material="K418",
                                          D=0.75, L=0.2, mass_kg=50.0, rotating=True)]
    bad = [r for r in F.check(c, t, p, g) if not r["ok"]]
    assert {r["code"] for r in bad} >= {"size:ротор целиком", "size:колесо турбины"}
    assert F.violation(bad) > 0 and F.violation([]) == 0.0


def test_fixed_speed_shaft_flow_coefficient(feasible_case):
    c, t, p, g = copy.deepcopy(feasible_case)
    k = next(k for k, s in enumerate(g["stages"]) if "phi" in s)
    g["stages"][k]["phi"] = 0.30
    assert "phi" in _codes(F.check(c, t, p, g))


# =========================================================================== parts_model
def _bom(**kw):
    b = copy.deepcopy(PM.REF_BOM)
    b.update(kw)
    return b


def test_parts_reference_scores_five():
    assert abs(PM.count(PM.REF_BOM)["score"] - 5.0) < 1e-9


def test_parts_score_strictly_decreasing_in_neff():
    xs = [PM.score_from_neff(n) for n in range(0, 2000, 7)]
    assert all(a > b for a, b in zip(xs, xs[1:]))
    assert 0.0 < xs[-1] < xs[0] <= 10.0


def test_parts_more_stages_worse():
    one = PM.count(_bom(comp=["Ц"]))
    two = PM.count(_bom(comp=["Ц", "Ц"]))
    assert two["parts"] > one["parts"] and two["score"] < one["score"]


def test_parts_bladed_disc_worse_than_blisk():
    bl = PM.count(_bom(turb=dict(type="О", stages=[dict(blisk=True, zt=55, zn=40, ngv_ring=True)] * 2)))
    bd = PM.count(_bom(turb=dict(type="О", stages=[dict(blisk=False, zt=55, zn=40, ngv_ring=True)] * 2)))
    assert bd["parts"] - bl["parts"] == 2 * 55 and bd["moving"] > bl["moving"] and bd["score"] < bl["score"]


def test_parts_bolts_vs_clamps_and_monolith():
    bolts = PM.count(_bom(conn="болты", case=dict(method="литьё", D=1.0)))
    clamps = PM.count(_bom(conn="быстроразъём", case=dict(method="литьё", D=1.0)))
    assert bolts["fasteners"] > clamps["fasteners"] and bolts["score"] < clamps["score"]
    mono = PM.count(_bom(rotor_join="заодно", turb=dict(type="Р", stages=[dict(blisk=True, zt=13, zn=18, ngv_ring=True)])))
    built = PM.count(_bom(rotor_join="сборный", turb=dict(type="Р", stages=[dict(blisk=True, zt=13, zn=18, ngv_ring=True)])))
    assert mono["by_unit"]["ротор"]["parts"] < built["by_unit"]["ротор"]["parts"] and mono["score"] > built["score"]


def test_parts_cooling_cores_and_bearings():
    u = PM.count(_bom())
    k = PM.count(_bom(ngv_cooled=True, rotor_cooled=True))
    assert k["cooled"] > u["cooled"] and k["score"] < u["score"]
    c1 = PM.count(_bom(rotors=2, pt=dict(stages=[dict(blisk=True, zt=40, zn=30, ngv_ring=True)]), cores=1))
    c4 = PM.count(_bom(rotors=2, pt=dict(stages=[dict(blisk=True, zt=40, zn=30, ngv_ring=True)]), cores=4))
    assert c4["parts"] > c1["parts"] and c4["unique"] < c4["parts"]
    mag = PM.count(_bom(bear="магнит"))
    assert "маслосистема" not in mag["by_unit"] and "маслосистема" in u["by_unit"]


def test_parts_hot_integral_fewer():
    a = PM.count(_bom(hot_integral=False, comb=dict(kind="кольцевая", method="печать", n_inj=14)))
    b = PM.count(_bom(hot_integral=True, comb=dict(kind="кольцевая", method="печать", n_inj=14)))
    assert b["parts"] < a["parts"] and b["score"] > a["score"]


# =========================================================================== service_model
def _mods(mass=20.0, bolts=0, align=False, tools=(), keyed=True, pipes=0):
    j = [dict(type="болты", n_bolts=bolts)] if bolts else [dict(type="хомут")]
    return dict(modules=[dict(name="горячая часть", mass_kg=mass, joints=j, pipes=pipes, pipes_keyed=keyed, cables=8, cables_keyed=keyed,
                              align=align, tools=list(tools), frequent=True, access=[])], amb=False, single_connector=keyed)


def test_service_bounds_and_mass_monotone():
    s = [SM.score(_mods(mass=m))["score"] for m in (10, 25, 100, 500, 1500, 2000, 5000)]
    assert all(a >= b for a, b in zip(s, s[1:])) and s[0] > s[-1]
    assert SM.score(_mods(mass=25))["s_mass"] == 1.0 and SM.score(_mods(mass=2000))["s_mass"] == 0.0
    for v in s:
        assert 0.0 <= v <= 10.0


def test_service_steps_tools_errors():
    base = SM.score(_mods())["score"]
    assert SM.score(_mods(bolts=40))["score"] < base
    assert SM.score(_mods(align=True))["score"] < base
    assert SM.score(_mods(tools=("гидравлический натяжитель",)))["score"] < base
    assert SM.score(_mods(keyed=False, pipes=3))["score"] < SM.score(_mods(keyed=True, pipes=3))["score"]


def test_service_access_adds_steps():
    a = dict(modules=[dict(name="камера", mass_kg=20, joints=[dict(type="хомут")], frequent=True, access=[]),
                      dict(name="горячая часть", mass_kg=20, joints=[dict(type="болты", n_bolts=30)], frequent=False, access=[])])
    b = copy.deepcopy(a); b["modules"][0]["access"] = ["горячая часть"]
    assert SM.score(b)["steps"] > SM.score(a)["steps"] and SM.score(b)["score"] < SM.score(a)["score"]


# =========================================================================== plugplay_model
def _pp(**kw):
    p = dict(n_out=14950.0, n_pkg=14950.0, G=21.4, G0=21.4, EGT=510.0, EGT0=510.0, D_env=1.8, pkg_W=2.6, pkg_H=3.2, rotors=1, rotors0=1,
             cores=1, inlet="одно", bear="масло", recup=False, recup_integral=False, comb="кольцевая", pi=12.0, pi0=12.2, eta=31.5,
             eta0=31.5, module_kg_max=800.0, engine_kg0=None)
    p.update(kw)
    return p


def test_plugplay_soft_tolerance_blocks():
    assert PP.score(_pp(G=21.4 * 1.25))["blocked"] and PP.score(_pp(G=21.4 * 1.25))["score"] == 0.0
    assert PP.score(_pp(EGT=630.0))["blocked"]
    assert PP.score(_pp(n_out=14000.0))["blocked"] == ["вал и муфта"]
    assert PP.score(_pp(D_env=3.0))["blocked"] == ["опоры и рама"]


def test_plugplay_monotone_inside_tolerance():
    s = [PP.score(_pp(G=21.4 * k))["score"] for k in (1.0, 1.05, 1.10, 1.15, 1.20)]
    assert all(a > b for a, b in zip(s, s[1:]))
    e = [PP.score(_pp(EGT=510.0 + d))["score"] for d in (0, 25, 50, 100)]
    assert all(a > b for a, b in zip(e, e[1:]))


def test_plugplay_interfaces():
    base = PP.score(_pp())
    assert base["items"]["САУ"]["state"] == "переходник" and "пуск" in base["no_data"]
    mag = PP.score(_pp(bear="магнит"))
    assert mag["items"]["САУ"]["state"] == "переделка" and mag["items"]["масло"]["state"] == "сохраняется"
    two = PP.score(_pp(rotors=2))
    assert two["items"]["пуск"]["state"] == "переделка" and two["score"] < base["score"]
    clu = PP.score(_pp(cores=4, rotors=2))
    assert clu["score"] < two["score"]
    assert PP.score(_pp(pi=12.2 * 1.3))["items"]["топливо"]["state"] == "переделка"
    assert abs(sum(PP.WEIGHTS.values()) - 1.0) < 1e-12


# =========================================================================== make_model
def _part(**kw):
    p = dict(name="колесо", kind="turb_radial", method="печать", material="IN718", D=0.3, L=0.1, mass_kg=20.0, rotating=True)
    p.update(kw)
    return p


def test_make_printer_fit_smallest_and_none():
    a = MK.assess_part(_part(D=0.3, L=0.1))
    assert a["fit"] and a["route"] == "EOS M 400"
    b = MK.assess_part(_part(D=1.2, L=0.4))
    assert b["fit"] and b["route"] in ("BLT-S1500", "Farsoon FS1521M-U", "Eplus3D EP-M2050")
    c = MK.assess_part(_part(D=2.5, L=0.4))
    assert not c["fit"] and c["pieces"] >= 2
    s_fit = MK.score(dict(parts=[_part(D=0.3)], techs=["традиция"]))["score"]
    s_no = MK.score(dict(parts=[_part(D=2.5)], techs=["традиция"]))["score"]
    assert s_no < s_fit


def test_make_casting_limits_and_quality():
    assert MK.assess_part(_part(method="литьё", material="K418", D=0.65))["fit"]
    assert not MK.assess_part(_part(method="литьё", material="K418", D=0.75))["fit"]
    assert MK.assess_part(_part(method="литьё", material="K418", D=0.65))["proc_q"] < MK.assess_part(_part(method="литьё", material="K418", D=0.45))["proc_q"]
    good = MK.score(dict(parts=[_part(material="IN718")], techs=["традиция"]))["score"]
    poor = MK.score(dict(parts=[_part(material="CM247LC")], techs=["традиция"]))["score"]
    assert poor < good


def test_make_trl_penalty():
    s9 = MK.score(dict(parts=[_part()], techs=["традиция"]))["score"]
    s3 = MK.score(dict(parts=[_part()], techs=["традиция", "rdc"]))
    assert s3["trl_min"] == 3 and s3["score"] < s9


def test_make_limits_match_track():
    assert MK.u_limit_comp("Ti64", "поковка")[0] == 620.0
    assert MK.u_limit_comp("Al", "поковка")[0] == 520.0
    assert MK.u_limit_comp("Ni", "литьё")[0] == 480.0
    assert MK.t41_limit("Р", "U", "K418", "литьё")[0] == 950.0
    assert MK.t41_limit("О", "N", "SC", "литьё")[0] == 1050.0
    assert MK.t41_limit("О", "K", "SC", "литьё")[0] == 1200.0
    assert MK.t41_limit("О", "U", "SC", "литьё")[0] == 1000.0
    assert MK.u_limit_turb("Р", "K418", 950.0)[0] == pytest.approx(465.0)
    for p in MK.PRINTERS:
        assert p["src"].startswith("https://") and p["status"]
