"""Проверки уровня 2 воронки (funnel_l2.py) и файлов задач уровней 2-3 (make_jobs_l23.py). Корпуса придуманы (SYNTHETIC);
CalculiX и SU2 не вызываются, кроме одной пробы disk.py, если ccx установлен."""
import json, math, os, shutil, sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import funnel_l2 as L2          # noqa: E402


def _variant(r, arch="кластер", turb="О"):
    """придуманный вариант отбора уровня 1 (SYNTHETIC): два центробежных колеса, осевая турбина, силовая турбина"""
    st = dict(PRs=3.0, eta_p=0.85, D2=0.40, b2=0.03, beta2b_deg=40.0, D1t=0.28, D1h=0.11, Zfull=11.0, cl=3e-4, M1t=1.1, U2=470.0, phi=0.1)
    disc = dict(R=0.11, sigma_d_MPa=450.0, t_rim=0.02, t_hub=0.04, ratio=2.0)
    return dict(div_rank=r, j=100 + r, idx=1000 + r, arch=arch, niche=f"{arch}|2р|2Ц|{turb}|P|нет|сварка:печать", scheme="SYNTHETIC",
                score=7.0 - 0.01 * r, eta_th_shaft_pct=30.0, n_rpm=23000.0, G=24.0, pi=9.0, t41_C=1030.0,
                concept=dict(cores=2, rotor_join="сварка", turb=turb, mat_tw="CM247LC"),
                meanline=dict(compressor=[st, dict(st, PRs=2.6, D2=0.38, dbl=False)], turbine=dict(stages=[{}, {}], D4=0.3)),
                geometry=dict(compressor=[dict(type="центробежная", blades=11, splitters=11), dict(type="центробежная", blades=11, splitters=0)],
                              combustor=dict(type="кольцевая", D=0.44, L=0.26, injectors=9), power_turbine=dict(n_rpm=14950.0)),
                strength=dict(comp=[dict(stage=1, material_data="in713c", T_in_C=15.0, T_out_C=150.0),
                                    dict(stage=2, material_data="in713c", T_in_C=150.0, T_out_C=300.0)],
                              turb=dict(T_blade_C=780.0, T_rim_C=780.0, disc=disc), pt=dict(T_blade_C=700.0, disc=dict(disc, R=0.15))))


@pytest.fixture()
def cases(tmp_path):
    src = tmp_path / "top50_SYNTHETIC.json"
    src.write_text(json.dumps(dict(kv="KV-3", M=7, code_sha="SYNTHETIC", selected=[_variant(r) for r in range(1, 8)]), ensure_ascii=False))
    out = tmp_path / "cases.json"
    return L2.prepare(str(src), 50, str(out)), str(out)


def test_prepare_stages_and_inlets(cases):
    doc, _ = cases
    v = doc["variants"][0]
    assert [s["id"] for s in v["stages"]] == ["r01_j101_c1", "r01_j101_c2"]
    a1, a2 = v["stages"][0]["args"], v["stages"][1]["args"]
    assert a1["m"] == pytest.approx(12.0) and a1["t01"] == pytest.approx(288.15) and a1["p01"] == pytest.approx(101325.0)
    assert a2["p01"] == pytest.approx(101325.0 * 3.0, rel=1e-6)                       # вход 2-й ступени - за 1-й
    assert a2["t01"] == pytest.approx(288.15 * 3.0 ** (0.4 / (1.4 * 0.85)), rel=1e-4)
    assert a2.get("no_splitters") and not a1.get("no_splitters")
    assert {d["id"].split(":")[1] for d in v["disks"]} == {"comp1", "comp2", "turb_disc", "pt_disc"}
    assert any(p["id"].endswith("turbine_cfd") for p in v["pending"]) and any(p["id"].endswith("comb1d") for p in v["pending"])


def test_axial_disc_profile():
    disc = dict(R=0.12, sigma_d_MPa=400.0, t_rim=0.02)
    c = L2._axial_disc_case("x", disc, 15000.0, 750.0)
    assert c["profile"]["t"][-1] == pytest.approx(0.02)
    assert c["rim"]["p_MPa"] == pytest.approx(400.0)                                  # равнопрочный: σr на ободе = σ_d
    assert all(a >= b for a, b in zip(c["profile"]["t"], c["profile"]["t"][1:]))     # к ободу тоньше
    assert c["T_C"]["T"][-1] == pytest.approx(750.0) and c["T_C"]["T"][0] == pytest.approx(750.0 - L2.T_HUB_DROP)


def test_single_split_by_variant(cases):
    doc, _ = cases
    S = 3
    seen = []
    for k in range(S):
        seen += [v["key"] for v in L2.mine(doc, S, [k])]
    assert sorted(seen) == sorted(v["key"] for v in doc["variants"])                 # каждый вариант - ровно в одной части
    assert L2.parse_shards("0,2", 3) == [0, 2] and L2.parse_shards("all", 3) == [0, 1, 2]
    with pytest.raises(ValueError):
        L2.parse_shards("3", 3)


def test_su2_commands(cases):
    doc, _ = cases
    st = doc["variants"][0]["stages"][0]
    c1, c2 = L2.su2_commands(st, "coarse", "/w", "/o", 64, 400, "plan_sweep_SYNTHETIC.json", 2)
    assert "--case" in c1 and c1[c1.index("--case") + 1] == "custom" and c1[c1.index("--d2") + 1] == str(st["args"]["d2"])
    assert c2[c2.index("--np-per-point") + 1] == "16" and c2[c2.index("--base-strategy") + 1] == "plan_sweep_SYNTHETIC.json"
    assert int(c2[c2.index("--time-limit") + 1]) == 200 * 60                          # срок делится между оставшимися ступенями
    _, c2b = L2.su2_commands(st, "medium", "/w", "/o", 64, 400)
    assert c2b[c2b.index("--base-plan") + 1] == "smoke075" and c2b[c2b.index("--np-per-point") + 1] == "8"


def _part(path, cases_path, k, S, res):
    d = dict(level=2, cases=os.path.basename(cases_path), cases_sha256=L2._sha(cases_path), shard=k, shards=S, results=res)
    path.write_text(json.dumps(d, ensure_ascii=False))


def test_merge_completeness(cases, tmp_path):
    doc, cp = cases
    S = 2
    pdir = tmp_path / "parts"
    pdir.mkdir()
    for k in range(S):
        vs = L2.mine(doc, S, [k])
        fea = {d["id"]: dict(burst_margin_120=1.2, life_local_h=1e6) for v in vs for d in v["disks"]}
        geom = {s["id"]: dict(rc=0) for v in vs for s in v["stages"]}
        _part(pdir / f"funnel_l2_t_s{k}.json", cp, k, S, dict(geom=geom, fea=fea))
    args = L2.argparse.Namespace(cases=cp, inp=str(pdir), what="geom,fea", shards=S, out_dir=str(tmp_path / "m"))
    assert L2.merge(args) == 0
    r = json.load(open(tmp_path / "m" / "funnel_l2_results.json", encoding="utf-8"))
    assert r["completeness"]["ok"] and r["n_pass"] == 7 and len(r["l3_select"]) == 2   # не больше 2 на архетип
    os.remove(pdir / "funnel_l2_t_s1.json")
    assert L2.merge(args) == 2                                                         # потеря части - отказ
    _part(pdir / "funnel_l2_t_s1.json", cp, 1, S, dict(geom={}, fea={}))
    _part(pdir / "funnel_l2_t_s1b.json", cp, 1, S, dict(geom={"x": dict(rc=1)}, fea={}))
    with pytest.raises(SystemExit):                                                    # часть дважды с разным содержанием
        L2.merge(args)


@pytest.mark.skipif(shutil.which("ccx") is None, reason="нет ccx")
def test_axial_disc_runs_in_disk_py(tmp_path):
    c = L2._axial_disc_case("syn_disc", dict(R=0.11, sigma_d_MPa=450.0, t_rim=0.02), 20000.0, 700.0)
    c["mesh"] = dict(nr=20, nz=4)
    r = L2.ST.run_disk_fea(c, str(tmp_path / "d"))
    assert r["equilibrium_ok"] and math.isfinite(r["burst_margin_120"])


def test_jobs_l23_files_exist_and_parse():
    for f in ("jobs_l2.json", "jobs_l3.json"):
        p = os.path.join(HERE, f)
        if not os.path.exists(p):
            pytest.skip("файлы задач ещё не собраны (make_jobs_l23.py)")
        d = json.load(open(p, encoding="utf-8"))
        assert d["jobs"] and all(set(j) <= {"name", "command", "max_min", "container", "container_setup", "disk_gb", "python", "pip"}
                                 for j in d["jobs"])
        assert d["plan_check"]["exit_code"] == 0
        assert "pending_jobs" in d and d["estimate"]["usd"]["total_expected"] > 0


def test_merge_l1_applies_aws_fea(tmp_path):
    """итог CalculiX лучших 300 (jobs_l1_fea.json) отсеивает до отбора 50 вариант с запасом при 120 % < 1,0 (SYNTHETIC)"""
    import merge_l1 as M1
    p = tmp_path / "fea_SYNTHETIC.json"
    p.write_text(json.dumps(dict(results=[dict(key="j1:comp1", result=dict(burst_margin_120=0.95)),
                                          dict(key="j1:comp2", result=dict(burst_margin_120=1.30)),
                                          dict(key="j2:comp1", result=dict(burst_margin_120=1.10)),
                                          dict(key="j3:turb", error="ccx упал")]), ensure_ascii=False))
    recs = [dict(j=1, feasible=True, reject="", reject_detail=""), dict(j=2, feasible=True, reject="", reject_detail=""),
            dict(j=3, feasible=True, reject="", reject_detail="")]
    by, n_bad = M1.apply_fea_results(recs, str(p))
    assert n_bad == 1 and not recs[0]["feasible"] and "fea_burst" in recs[0]["reject"] and "comp1: 0.95" in recs[0]["reject_detail"]
    assert recs[1]["feasible"] and recs[2]["feasible"] and set(by) == {1, 2, 3}
