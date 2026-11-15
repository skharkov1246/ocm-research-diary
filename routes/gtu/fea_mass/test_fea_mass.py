#!/usr/bin/env python3
"""Проверки узла массовой 3D-прочности (fea_mass). Запуск: python3 test_fea_mass.py [--quick] [--json итог.json]

1. ring_centrif   - вращающееся кольцо (сектор 40°, C3D10 из geom3d, циклическая симметрия, закрепление одним узлом):
                    σθ у отверстия и на ободе против аналитики тонкого диска; средние окружные (post3d.robinson) против
                    равновесия полудиска ρω²(b³ - a³)/(3(b - a)).
2. ring_hoopfree  - мёртвый обод без окружной жёсткости (ccx3d.ring_cards): σθ в кольце ~ 0 при вращении.
3. beam_freq      - консольная балка C3D10: первая частота изгиба против Эйлера - Бернулли.
4. spin_campbell  - кольцо-сектор с частотами при 0 и n (преднагружение, NLGEOM): низшая частота каждого ND растёт,
                    Саутвелл (post3d.southwell) по 3 скоростям - невязка < 1 %; пересечение с EO находится.
5. eigen_reader, frd_reader, alias, verdict - разбор и логика на придуманных данных (корпуса - придуманные, правило 18).
6. repeat         - повторяемость: та же деталь (осевой блиск, грубая сетка) дважды -> те же числа.
7. vs_disk        - сравнение с осесимметричным ../platform/fea/disk.py: колесо компрессора варианта 1 funnel_l1_top50.json;
                    профиль ступицы и масса лопаток disk.py - из той же 3D-геометрии; средние окружные напряжения и запас по
                    разрушающей частоте должны совпасть в пределах 10 %.
8. merge          - слияние частей: полнота, задвоение, чужая часть -> отказ (придуманные части).
--quick - без 6 и 7 (долгие).
"""
import gzip
import json
import math
import os
import shutil
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import geom3d as G           # noqa: E402
import ccx3d as C            # noqa: E402
import post3d as Q           # noqa: E402

RHO, NU, E_MOD = 7850.0, 0.3, 2.1e11


class Steel:
    """придуманный материал проверки (сталь, без температурной зависимости)"""
    name, alloy, key = "сталь (проверка)", "steel", "steel"
    rho, nu, T_ref = RHO, NU, 20.0
    E_tab = [(20.0, E_MOD / 1e9), (1000.0, E_MOD / 1e9)]
    alpha_tab = [(20.0, 12.0), (1000.0, 12.0)]
    uts_tab = [[20.0, 600.0], [1000.0, 600.0]]
    ys_tab, ys20, psi, creep = None, None, None, False

    def uts(self, T):
        return 600.0

    def ys(self, T):
        return None

    def E(self, T):
        return E_MOD


def centrif_hoop(r, w, a, b):
    return (3 + NU) / 8 * RHO * w * w * (a * a + b * b + a * a * b * b / r ** 2 - (1 + 3 * NU) / (3 + NU) * r * r)


def _work():
    return tempfile.mkdtemp(prefix="fea3d_test_", dir=os.environ.get("FEA3D_TMP") or None)


def t_ring_centrif(work):
    a, b, thk, Z, rpm = 0.05, 0.20, 0.004, 9, 10000.0
    m = G.ring_sector(a, b, thk, Z, 0.006)
    d = os.path.join(work, "ring")
    os.makedirs(d)
    m.write_inp(os.path.join(d, "mesh.inp"))
    T = np.full(len(m.node_ids), 20.0)
    r = C.run_ccx(d, "static", C.static_deck(Z, C.material_cards(Steel()), rpm, [1.0], m.node_ids, T))
    assert r["ok"], r
    st = C.read_frd_steps(os.path.join(d, "static.frd"))
    ids, S = st[0]["STRESS"]
    Sn = Q.align(m, ids, S[:, :6])
    w = 2 * math.pi * rpm / 60
    rn = np.hypot(m.nodes[:, 0], m.nodes[:, 1])
    hoop = Q.hoop(Sn, m.nodes)
    rows = []
    for rr in (a, b):
        sel = np.abs(rn - rr) < 1e-7
        fe = float(np.median(hoop[sel]))
        an = centrif_hoop(rr, w, a, b)
        rows.append(dict(r=rr, fe_MPa=fe / 1e6, an_MPa=an / 1e6, dev=fe / an - 1))
    rob = Q.robinson(m, Steel(), Sn, T)
    mean_an = RHO * w * w * (b ** 3 - a ** 3) / (3 * (b - a)) / 1e6
    dev_mean = rob["hoop_mean_MPa"] / mean_an - 1
    ok = all(abs(x["dev"]) < 1e-2 for x in rows) and abs(dev_mean) < 5e-3 and m.meta["periodic"]["ok"]
    return dict(ok=ok, rows=rows, hoop_mean_fe_MPa=rob["hoop_mean_MPa"], hoop_mean_an_MPa=mean_an, dev_mean=dev_mean,
                n_nodes=m.meta["n_nodes"], periodic=m.meta["periodic"], bc="один узел отверстия: окружное и осевое (*TRANSFORM, TYPE=C)")


def t_ring_hoopfree(work):
    """кольцо целиком - мёртвый обод: окружные напряжения должны быть малы против изотропного кольца"""
    a, b, thk, Z, rpm = 0.05, 0.20, 0.004, 9, 10000.0
    m = G.ring_sector(a, b, thk, Z, 0.008)
    m.ecls[:] = G.CLASS_RING
    d = os.path.join(work, "ringfree")
    os.makedirs(d)
    m.write_inp(os.path.join(d, "mesh.inp"))
    T = np.full(len(m.node_ids), 20.0)
    cards = C.ring_cards(Steel(), RHO, 1.0, 1e-3)
    r = C.run_ccx(d, "static", C.static_deck(Z, cards, rpm, [1.0], m.node_ids, T))
    assert r["ok"], r
    st = C.read_frd_steps(os.path.join(d, "static.frd"))
    ids, S = st[0]["STRESS"]
    Sn = Q.align(m, ids, S[:, :6])
    hoop = np.abs(Q.hoop(Sn, m.nodes))
    w = 2 * math.pi * rpm / 60
    ref = centrif_hoop(a, w, a, b)
    ratio = float(np.nanmax(hoop) / ref)
    return dict(ok=ratio < 0.02, hoop_max_over_isotropic=ratio,
                note="E_θ = 1e-3 E: окружная связь обода порвана - σθ << изотропного кольца (проверка ориентации CYLINDRICAL)")


def t_beam(work):
    L, w, h = 0.10, 0.02, 0.005
    m = G.beam(L, w, h, 3)
    d = os.path.join(work, "beam")
    os.makedirs(d)
    m.write_inp(os.path.join(d, "mesh.inp"))
    txt = "\n".join(["*HEADING", "beam", "*INCLUDE, INPUT=mesh.inp"] + C.material_cards(Steel()) +
                    ["*BOUNDARY", "NROOT, 1, 3", "*STEP", "*FREQUENCY", "4", "*END STEP"]) + "\n"
    r = C.run_ccx(d, "beam", txt)
    assert r["ok"], r
    f = C.read_eigen_steps(os.path.join(d, "beam.dat"))[0][0]
    I = w * h ** 3 / 12
    f1 = 1.8751 ** 2 / (2 * math.pi) * math.sqrt(E_MOD * I / (RHO * w * h * L ** 4))
    dev = f[0] / f1 - 1
    return dict(ok=abs(dev) < 0.02, f1_fe_Hz=f[0], f1_euler_Hz=f1, dev=dev, n_nodes=m.meta["n_nodes"],
                note="L/h = 20: поправка Тимошенко около 0,3 %; три элемента по толщине")


def t_spin_campbell(work):
    a, b, thk, Z = 0.05, 0.20, 0.004, 12
    rpm = 12000.0
    m = G.ring_sector(a, b, thk, Z, 0.008)
    d = os.path.join(work, "spin")
    os.makedirs(d)
    m.write_inp(os.path.join(d, "mesh.inp"))
    T = np.full(len(m.node_ids), 20.0)
    fr = [0.0, 0.5, 1.0]
    r = C.run_ccx(d, "freq", C.freq_deck(Z, C.material_cards(Steel()), rpm, fr, m.node_ids, T, 8, Z // 2))
    assert r["ok"], r
    eig = C.read_eigen_steps(os.path.join(d, "freq.dat"))
    assert len(eig) == len(fr), len(eig)
    table = list(zip(fr, eig))
    stiff = all(eig[-1][nd][0] > eig[0][nd][0] for nd in eig[0] if nd in eig[-1])
    sw = [Q.southwell(table, rpm, nd, 0) for nd in sorted(eig[0])]
    res = max(x["max_rel_residual"] for x in sw if x)
    op = dict(range=(0.0, 1.25), info_range=None)
    cb = Q.campbell(table, rpm, Z, [dict(eo=2, kind="low", source="EO 2"), dict(eo=3, kind="low", source="EO 3")], op)
    return dict(ok=bool(stiff and res < 0.01 and len(eig[0]) == Z // 2 + 1), spin_stiffening=stiff, southwell_max_residual=res,
                nd=sorted(eig[0]), crossings=cb["crossings"][:4], f_ND2=[e.get(2, [None])[0] for e in eig])


def t_readers():
    """разбор .dat и .frd на придуманных файлах (корпуса - придуманные)"""
    d = _work()
    try:
        dat = os.path.join(d, "x.dat")
        with open(dat, "w") as f:
            f.write("""                        S T E P       1


     E I G E N V A L U E   O U T P U T

  NODAL   MODE NO    EIGENVALUE                      FREQUENCY
 DIAMETER                                    REAL PART            IMAGINARY PART
                                    (RAD/TIME)      (CYCLES/TIME)   (RAD/TIME)

    0          1   0.1000000E+08   0.3162278E+04   0.5032921E+03   0.0000000E+00
    0          2   0.1000000E+08   0.3162278E+04   0.5032921E+03   0.0000000E+00

     P A R T I C I P A T I O N   F A C T O R S

     E I G E N V A L U E   O U T P U T

  NODAL   MODE NO    EIGENVALUE                      FREQUENCY
    1          1   0.4000000E+08   0.6324555E+04   0.1006584E+04   0.0000000E+00
    1          2   0.4000000E+08   0.6324555E+04   0.1006584E+04   0.0000000E+00

                        S T E P       3


     E I G E N V A L U E   O U T P U T

  NODAL   MODE NO    EIGENVALUE                      FREQUENCY
    0          1   0.1210000E+08   0.3478505E+04   0.5536210E+03   0.0000000E+00
""")
        st = C.read_eigen_steps(dat)
        ok1 = len(st) == 2 and st[0] == {0: [503.2921], 1: [1006.584]} and st[1] == {0: [553.621]}
        frd = os.path.join(d, "x.frd")
        rows = []
        for step, val in ((1, 1.0), (1, 2.0), (2, 3.0)):
            rows.append(f"    1PSTEP                         1           1           {step}          ")
            rows.append(" -4  DISP        4    1")
            rows.append(" -1         1" + f"{val:12.5E}{0.0:12.5E}{0.0:12.5E}")
            rows.append(" -3")
        open(frd, "w").write("\n".join(rows) + "\n")
        fs = C.read_frd_steps(frd)
        ok2 = len(fs) == 2 and abs(fs[0]["DISP"][1][0, 0] - 2.0) < 1e-12 and abs(fs[1]["DISP"][1][0, 0] - 3.0) < 1e-12
        al = [Q.alias(e, 57) for e in (40, 57, 80, 1, 29, 28)]
        ok3 = al == [17, 0, 23, 1, 28, 28]
        return dict(ok=bool(ok1 and ok2 and ok3), eigen=ok1, frd_last_increment=ok2, alias=al)
    finally:
        shutil.rmtree(d, ignore_errors=True)


def t_verdict():
    """логика вердикта на придуманной детали"""
    base = dict(burst=dict(margin_120=1.2), static=dict(steps=[dict(yield_ratio_avg_max=0.8, yield_ratio_node_max=0.9)]),
                creep={}, lcf={}, campbell=dict(crossings=[], margin_min=None, free_stator=None))
    v1 = Q.verdict(base)[0]
    b2 = json.loads(json.dumps(base))
    b2["burst"]["margin_120"] = 0.9
    v2 = Q.verdict(b2)
    b3 = json.loads(json.dumps(base))
    b3["campbell"]["crossings"] = [dict(in_range=True, kind="ngv", mode=1, eo=40, nd=17, rpm_frac=0.99, source="СА 40")]
    v3 = Q.verdict(b3)[0]
    b4 = json.loads(json.dumps(base))
    b4["campbell"]["crossings"] = [dict(in_range=True, kind="low", mode=1, eo=4, nd=4, rpm_frac=0.9, source="EO 4")]
    v4 = Q.verdict(b4)[0]
    ok = v1 == "годен" and v2[0] == "отказ" and "120" in v2[1][0] and v3 == "отказ" and v4 == "условно"
    return dict(ok=ok, verdicts=[v1, v2[0], v3, v4])


def t_merge():
    """слияние: полные части - принято; пропуск и чужая часть - отказ (придуманные результаты)"""
    import fea_mass as F
    d = _work()
    try:
        cases = [dict(j=j, idx=100 + j, kv="KV-T") for j in range(5)]
        cp = os.path.join(d, "cases.json")
        json.dump(dict(source=dict(kv="KV-T"), cases=cases), open(cp, "w"))
        sha = F._sha(cp)

        def part(k, S, rows, path_dir):
            os.makedirs(path_dir, exist_ok=True)
            meta = dict(level=F.LEVEL_TAG, kv="KV-T", n=5, shards=S, shard=k, cases_sha256=sha, code_sha="x", opts={}, j=[r["j"] for r in rows])
            with gzip.open(os.path.join(path_dir, F.part_name("KV-T", 5, k, S)), "wt") as fh:
                json.dump(dict(meta, results=rows), fh)
        row = lambda c: dict(j=c["j"], idx=c["idx"], verdict="годен", reasons_hard=[], reasons_soft=[], errors=[], parts=[], summary={})
        good = os.path.join(d, "good")
        for k in range(2):
            part(k, 2, [row(c) for r, c in enumerate(cases) if r % 2 == k], good)
        rows, ref, checks, _ = F.merge(cp, [good])
        ok_full = len(rows) == 5 and checks["indices_ok"]
        miss = os.path.join(d, "miss")
        part(0, 2, [row(c) for r, c in enumerate(cases) if r % 2 == 0], miss)
        try:
            F.merge(cp, [miss])
            ok_miss = False
        except F.MergeError:
            ok_miss = True
        bad = os.path.join(d, "bad")
        part(0, 2, [row(c) for r, c in enumerate(cases) if r % 2 == 0], bad)
        part(1, 2, [row(c) for r, c in enumerate(cases) if r % 2 == 0], bad)       # чужие номера
        try:
            F.merge(cp, [bad])
            ok_bad = False
        except F.MergeError:
            ok_bad = True
        return dict(ok=bool(ok_full and ok_miss and ok_bad), full=ok_full, missing_rejected=ok_miss, foreign_rejected=ok_bad)
    finally:
        shutil.rmtree(d, ignore_errors=True)


def _top50_variant(i=0):
    import fea_mass as F
    import parts3d as P
    sel, _ = F.prepare_from_top50(os.path.join(os.path.dirname(HERE), "funnel", "funnel_l1_top50.json"))
    return sel[i], P


def t_repeat(work):
    """та же деталь дважды: осевая ступень газогенератора варианта 1 (грубая сетка, частоты при 0 и 100 %)"""
    import fea_mass as F
    v, P = _top50_variant(0)
    parts, _ = P.variant_parts(v)
    p = next(x for x in parts if x["kind"] == "axial")
    opts = dict(F.DEFAULTS, level="coarse", freq_fracs=[0.0, 1.0])
    r1 = F.run_part(p, os.path.join(work, "rep1"), opts)
    r2 = F.run_part(p, os.path.join(work, "rep2"), opts)
    h1, h2 = F._content_hash([r1]), F._content_hash([r2])
    return dict(ok=bool(r1["status"] == "ok" and h1 == h2), part=p["id"], status=[r1["status"], r2["status"]], same_hash=h1 == h2,
                burst=[(r1.get("burst") or {}).get("margin"), (r2.get("burst") or {}).get("margin")], error=r1.get("error"))


def _hub_profile(poly, r_lo, r_hi, n=24):
    """толщина и середина ступицы по меридиональному многоугольнику (z, r) на радиусах r_lo..r_hi"""
    rs = np.linspace(r_lo, r_hi, n)
    ts, zs = [], []
    P = np.asarray(poly)
    for r in rs:
        xs = []
        for (z1, r1), (z2, r2) in zip(P, np.roll(P, -1, axis=0)):
            if (r1 - r) * (r2 - r) < 0 or (r1 == r and r2 != r):
                xs.append(z1 + (r - r1) * (z2 - z1) / (r2 - r1))
        xs = sorted(xs)
        L = sum(xs[k + 1] - xs[k] for k in range(0, len(xs) - 1, 2))
        ts.append(max(L, 1e-4))
        zs.append(0.5 * (xs[0] + xs[-1]) if xs else 0.0)
    return rs, np.asarray(ts), np.asarray(zs)


def t_vs_disk(work, level="coarse"):
    """3D-сектор против disk.py на одной геометрии: колесо компрессора 1 варианта 1 funnel_l1_top50.json"""
    import fea_mass as F
    v, P = _top50_variant(0)
    parts, _ = P.variant_parts(v)
    p = next(x for x in parts if x["id"] == "comp1")
    p = dict(p, geom=dict(p["geom"], bore="tie"))         # отверстие под стяжку: у disk.py то же отверстие
    opts = dict(F.DEFAULTS, level=level, freq_fracs=[])
    t0 = time.time()
    mesh = G.wheel_sector(p["geom"], level=level)
    r3 = F.run_part(p, os.path.join(work, "vs3d"), opts)
    if r3["status"] != "ok":
        return dict(ok=False, error=r3.get("error"))
    # профиль disk.py: толщина ступицы по многоугольнику, лопатки - масса по радиусу из 3D-сетки
    poly = mesh.meta["hub_polygon"]
    r_b, r2 = mesh.meta["r_bore"], mesh.meta["r2"]
    rs, ts, zs = _hub_profile(poly, r_b * 1.0001, r2 * 0.9999)
    vol, cen = Q.elem_geometry(mesh)
    rc = np.hypot(cen[:, 0], cen[:, 1])
    bl = mesh.ecls == G.CLASS_BLADE
    edges = np.concatenate([[r_b], 0.5 * (rs[1:] + rs[:-1]), [r2 * 1.05]])
    Vb = np.histogram(rc[bl], bins=edges, weights=vol[bl])[0] * mesh.meta["Z"]
    hb = Vb / (2 * math.pi * rs * np.diff(edges))
    sys.path.insert(0, os.path.join(os.path.dirname(HERE), "platform", "fea"))
    import disk as DK
    mats = F.part_materials(p)
    mat_key = mats[G.CLASS_HUB].key
    ts_ = p["T_spec"]
    case = dict(name="vs3d", material=mat_key, rpm=p["rpm"], profile=dict(r=rs.tolist(), t=ts.tolist(), zc=zs.tolist(), hb=hb.tolist()),
                profile_interp="linear", blade_fill=1.0,
                T_C=dict(r=rs.tolist(), T=(ts_["T_in"] + (ts_["T_out"] - ts_["T_in"]) * ((rs - r_b) / (r2 - r_b)) ** 2).tolist()),
                burst_k=1.0, mesh=dict(nr=80, nz=12))
    rd = DK.run_case(case, workdir=os.path.join(work, "vs_disk"), keep=False)
    hm3, hmd = r3["burst"]["hoop_mean_MPa"], rd["section"]["hoop_mean_ip_MPa"]
    b3, bd = r3["burst"]["margin"], rd["burst"]["margin"]
    d_h, d_b = hm3 / hmd - 1, b3 / bd - 1
    return dict(ok=bool(abs(d_h) < 0.10 and abs(d_b) < 0.10), hoop_mean_3d_MPa=hm3, hoop_mean_disk_MPa=hmd, dev_hoop=d_h,
                burst_3d=b3, burst_disk=bd, dev_burst=d_b, vm_max_disk_MPa=rd["vm_max_node_MPa"],
                vm_hub_avg_3d_MPa=r3["static"]["steps"][0]["regions"]["hub"]["vm_peak_avg_MPa"], level=level,
                n_nodes_3d=r3["mesh"]["n_nodes"], t_s=round(time.time() - t0, 1),
                note="одна геометрия: профиль ступицы disk.py снят с меридионального многоугольника 3D-сектора, масса лопаток по "
                     "радиусу - из элементов лопаток 3D-сетки (blade_fill = 1, hb = dV/(2πr dr)); лопатки в disk.py окружной "
                     "нагрузки не несут - в 3D несут частично")


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--json", default=None)
    ap.add_argument("--only", default=None)
    a = ap.parse_args(argv)
    work = _work()
    tests = [("readers", lambda: t_readers()), ("verdict", lambda: t_verdict()), ("merge", lambda: t_merge()),
             ("ring_centrif", lambda: t_ring_centrif(work)), ("ring_hoopfree", lambda: t_ring_hoopfree(work)),
             ("beam_freq", lambda: t_beam(work)), ("spin_campbell", lambda: t_spin_campbell(work))]
    if not a.quick:
        tests += [("repeat", lambda: t_repeat(work)), ("vs_disk", lambda: t_vs_disk(work))]
    if a.only:
        tests = [t for t in tests if t[0] in a.only.split(",")]
    out, ok_all = {}, True
    try:
        for name, fn in tests:
            t0 = time.time()
            try:
                r = fn()
            except Exception as e:                       # noqa: BLE001
                import traceback
                r = dict(ok=False, error=f"{type(e).__name__}: {e}", trace=traceback.format_exc()[-1200:])
            r["t_s"] = round(time.time() - t0, 1)
            out[name] = r
            ok_all &= bool(r["ok"])
            print(f"{'ok ' if r['ok'] else 'НЕТ'} {name} ({r['t_s']} с)", {k: v for k, v in r.items() if k not in ("trace",)} if not r["ok"] else "")
    finally:
        shutil.rmtree(work, ignore_errors=True)
    res = dict(all_ok=bool(ok_all), ccx=C.ccx_version(), tests=out, date="2026-11-15")
    if a.json:
        with open(a.json, "w", encoding="utf-8") as fh:
            json.dump(res, fh, ensure_ascii=False, indent=1, default=float)
    print("итог:", "все проверки пройдены" if ok_all else "есть непройденные")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
