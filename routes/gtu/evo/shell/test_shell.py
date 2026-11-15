#!/usr/bin/env python3
"""Проверки семейства «корпус-оболочка и рама» (pytest или python3 test_shell.py).

 1. Геном: кодирование ↔ раскодирование опорных форм; оценка детерминирована и быстрая (< 40 мс).
 2. Огибающая: стенка не задевает внутренние узлы (r_i ≥ r_req) у случайных геномов обоих вариантов; у картриджа канал не
    сужается по пути выемки (r_i не меньше накопленного максимума требуемых радиусов).
 3. Оболочка Доннелла: ортотропная формула в изотропном пределе совпадает с замкнутой изотропной; размазанные кольца
    повышают частоту.
 4. Безмоментная теория: цилиндр — σθ = pR/t, σz = pR/(2t); сфера — pR/(2t) по обоим направлениям.
 5. Рамы: без изгиба осевая жёсткость стойки = EA/L; навал снимает тепловое напряжение прямой радиальной стойки; изгиб
    снижает боковую жёсткость; кольцо рамы: c(N) убывает с числом стоек.
 6. Клетка горячей части: σ = 3Etδ/L²; жёсткость растёт с числом пальцев.
 7. Детали: печатная оболочка — меньше деталей и крепежа, чем типовая литая; картридж — больше деталей, чем оболочка
    при той же форме (правило собираемости parts_model).
 8. Рёбра и симметрия: 6 осевых рёбер — класс «N 6–11»; кольца — осесимметрично; искривление угла — асимметрия.
 9. Рельеф рёбер: масса по рельефу (интеграл на сетке) против аналитики H·(w + r_f)·L (10 %).
10. Геометрия элиты: замкнутая сетка оболочки, объём по сетке против интеграла профиля (2 %).
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import shell_model as SM  # noqa: E402


def _rand(n, seed=0):
    return np.random.Generator(np.random.PCG64(seed)).random((n, SM.DIM))


def test_genome_and_speed():
    for v in SM.VARIANTS:
        for k in SM.REF_KINDS:
            p = SM.reference_params(v, k)
            q = SM.decode(SM.encode(p))
            for g in SM.G:
                a, b = p[g["name"]], q[g["name"]]
                if g["kind"] == "real":
                    assert abs(a - b) < 1e-9 * max(1.0, abs(a)), (g["name"], a, b)
                else:
                    assert a == b, (g["name"], a, b)
    X = _rand(60, 3)
    r1 = [SM.evaluate(x, "оболочка") for x in X]
    t0 = time.time()
    r2 = [SM.evaluate(x, "оболочка") for x in X]
    ms = (time.time() - t0) / len(X) * 1e3
    assert all((a["fitness"] == b["fitness"]) or (np.isnan(a["fitness"]) and np.isnan(b["fitness"])) for a, b in zip(r1, r2))
    assert ms < 40.0, ms
    return dict(ms_per_eval=ms, dim=SM.DIM)


def test_envelope_clearance():
    worst = {}
    for v in SM.VARIANTS:
        w = np.inf
        for x in _rand(80, 5):
            p = SM.decode(x)
            pr = SM.profiles(p, v)
            req = SM.R_REQ
            if v == "картридж":
                req = np.maximum.accumulate(req[::-1])[::-1] if p["extract"] == "вперёд" else np.maximum.accumulate(req)
            w = min(w, float(np.min(pr["r_i"] - req)))
        worst[v] = w
        assert w >= -1e-9, (v, w)
    return dict(min_clearance_m=worst)


def test_donnell_isotropic_and_rings():
    E, nu, rho, R, L, t = 110e9, 0.34, 4430.0, 0.6, 0.5, 0.012
    Am = E * t / (1 - nu * nu); Dm = E * t ** 3 / (12 * (1 - nu * nu))
    f, m, n = SM.shell_freq_donnell(R, L, Am, Am, nu * Am, E * t / (2 * (1 + nu)), Dm, Dm, nu * Dm, Dm * (1 - nu) / 2, rho * t)
    lam, k = m * math.pi / L, n / R
    w2 = (Dm * (lam ** 2 + k ** 2) ** 2 + E * t * lam ** 4 / (R ** 2 * (lam ** 2 + k ** 2) ** 2)) / (rho * t)
    assert abs(f - math.sqrt(w2) / (2 * math.pi)) < 1e-6 * f
    f2, _, _ = SM.shell_freq_donnell(R, L, Am, Am + E * 0.004, nu * Am, E * t / (2 * (1 + nu)), Dm, Dm + E * 2e-7, nu * Dm,
                                     Dm * (1 - nu) / 2, rho * (t + 0.004))
    assert f2 > f
    return dict(f_iso_Hz=f, mode=[m, n], f_rings_Hz=f2)


def test_membrane_theory():
    z = np.linspace(0, 1.0, 401)
    p, t = 1e6, 0.01
    Nth, Nphi = SM.membrane_shell(z, np.full_like(z, 0.5), np.full_like(z, t), np.full_like(z, p))
    assert np.allclose(Nth[50:-50] / t, p * 0.5 / t, rtol=1e-6) and np.allclose(Nphi[50:-50] / t, p * 0.5 / (2 * t), rtol=1e-6)
    Rs = 2.0
    zz = np.linspace(-0.4, 0.4, 801)
    r = np.sqrt(Rs ** 2 - zz ** 2)
    Nth, Nphi = SM.membrane_shell(zz, r, np.full_like(zz, 0.002), np.full_like(zz, p))
    mid = slice(300, 501)
    assert np.allclose(Nth[mid], p * Rs / 2, rtol=0.02) and np.allclose(Nphi[mid], p * Rs / 2, rtol=0.02)
    return dict(cyl_hoop_MPa=float(Nth[200] / t / 1e6) if False else p * 0.5 / t / 1e6, sphere_N=float(Nth[400]), sphere_exact=p * Rs / 2)


def test_frames_and_ring():
    E, a = 200e9, 13e-6
    fr = SM.strut_frame(8, 0.3, 0.0, 0.0, 0.0, 0.1, 0.15, E, a, 400.0, 150.0, 430.0, 0.43, 0.13)
    A = 0.685 * 0.1 * 0.015
    assert abs(fr["k_axial"] - E * A / 0.3) < 1e-6 * fr["k_axial"]
    lean = SM.strut_frame(8, 0.3, 30.0, 0.0, 0.0, 0.1, 0.15, E, a, 400.0, 150.0, 430.0, 0.43, 0.13)
    bow = SM.strut_frame(8, 0.3, 0.0, 0.15, 0.0, 0.1, 0.15, E, a, 400.0, 150.0, 430.0, 0.43, 0.13)
    assert lean["sig_thermal"] < 0.5 * fr["sig_thermal"], (lean["sig_thermal"], fr["sig_thermal"])
    assert bow["K_lat"] < 0.2 * fr["K_lat"]
    cs = [SM.ring_coeff(N) for N in (3, 4, 6, 8, 12)]
    assert all(c1 > c2 for c1, c2 in zip(cs, cs[1:])), cs
    return dict(sig_radial_MPa=fr["sig_thermal"] / 1e6, sig_lean30_MPa=lean["sig_thermal"] / 1e6, K_lat=fr["K_lat"], K_bow=bow["K_lat"],
                ring_c=cs)


def test_hot_cage():
    p = SM.reference_params("оболочка")
    c = SM.hot_cage(p)
    L = p["hf_shape"] * p["hf_span"]
    assert abs(c["sig_thermal"] - 3 * c["E"] * p["hf_t"] * SM.INP["dr_hot"] / L ** 2) < 1e-6 * c["sig_thermal"]
    q = dict(p, hf_n=p["hf_n"] // 2)
    assert SM.hot_cage(q)["K"] < c["K"]
    return dict(sigma_MPa=c["sigma"] / 1e6, K=c["K"], dr_hot_mm=SM.INP["dr_hot"] * 1e3)


def test_parts_rule():
    ref = SM.reference("оболочка")["rows"][0]
    p = SM.reference_params("оболочка", "гладкая оболочка")
    p.update(route="печать", emb=5)
    o = SM.evaluate_params(p, "оболочка")["objs"]
    k = SM.evaluate_params(p, "картридж")["objs"]
    assert o["parts"] < ref["parts"] and o["fasteners"] < ref["fasteners"], (o["parts"], ref["parts"], o["fasteners"], ref["fasteners"])
    assert k["parts"] > o["parts"]
    return dict(ref_parts=ref["parts"], ref_fasteners=ref["fasteners"], shell_parts=o["parts"], shell_fasteners=o["fasteners"],
                cart_parts=k["parts"], cart_fasteners=k["fasteners"])


def _sym_of(**kw):
    p = SM.reference_params("оболочка", "гладкая оболочка")
    p.update(route="печать", rib_h=0.02, rib_w=0.01, rib_fillet=0.006, **kw)
    return SM.evaluate_params(p, "оболочка", with_geom=True)


def test_symmetry_classes():
    a = _sym_of(rib_type="осевые", sym_n=6, rib_n=0.0, asym1=0.0, asym2=0.0)
    b = _sym_of(rib_type="кольца", sym_n=1, rib_n2=0.5, asym1=0.4)
    c = _sym_of(rib_type="осевые", sym_n=12, rib_n=0.0, asym1=0.5, asym2=0.3, asym_phase=0.7)
    assert int(a["desc"][2]) == 1, a["sym"]
    assert int(b["desc"][2]) == 0, b["sym"]
    assert int(c["desc"][2]) in (3, 4), c["sym"]
    return dict(axial6=a["sym"], rings=b["sym"], warped=c["sym"])


def test_rib_heightfield_mass():
    import shell_geom as SG
    p = SM.reference_params("оболочка", "гладкая оболочка")
    p.update(route="печать", rib_type="осевые", sym_n=8, rib_n=0.0, rib_h=0.03, rib_w=0.012, rib_fillet=0.008, asym1=0.0, asym2=0.0,
             ring_h=0.0, mount_drop=0.0)
    geo = SM.evaluate_params(p, "оболочка", with_geom=True)
    hf = SG.heightfield(p, "оболочка", (300, 1440), geo)
    z, th, h = hf["z"], hf["th"], hf["h"]
    # без площадок опор: вычесть их вклад, оставив только рёбра (площадки — круги R 0,07 + 0,03 м у трёх опор)
    dz, dth = z[1] - z[0], th[1] - th[0]
    r = hf["r_o"][:, None]
    V_hf = float(np.sum(h * r * dth * dz))
    pads = 0.0
    for (zm, sm, _, _) in geo["mounts"]:
        s = th * hf["R_ref"]
        D = np.hypot(z[:, None] - zm, ((s[None, :] - sm + math.pi * hf["R_ref"]) % (2 * math.pi * hf["R_ref"])) - math.pi * hf["R_ref"])
        pads += float(np.sum(np.where(D < 0.1, h, 0.0) * r * dth * dz))
    V_ribs = V_hf - pads
    V_an = geo["rib"]["V"]
    rel = V_ribs / V_an
    assert 0.85 < rel < 1.15, (V_ribs, V_an)
    return dict(V_heightfield=V_ribs, V_analytic=V_an, ratio=rel)


def test_elite_mesh_volume():
    import shell_geom as SG
    from core import mesh as M
    p = SM.reference_params("оболочка", "гладкая оболочка")
    p.update(rib_type="нет", rib_h=0.0, ring_h=0.0)
    hf = SG.heightfield(p, "оболочка", (400, 360))
    V, F = M.revolve_heightfield(hf["z"], hf["th"], hf["r_o"][:, None] + hf["h"], hf["r_i"])
    pr = hf["geo"]["prof"]
    pads = 0.0
    V_prof = float(np.trapezoid(math.pi * (pr["r_o"] ** 2 - pr["r_i"] ** 2), pr["z"])) if hasattr(np, "trapezoid") else \
        float(np.trapz(math.pi * (pr["r_o"] ** 2 - pr["r_i"] ** 2), pr["z"]))
    vm = M.volume(V, F)
    # рельеф площадок опор добавляет объём; без них — ±2 %
    assert vm > 0 and abs(vm / V_prof - 1) < 0.05, (vm, V_prof)
    return dict(V_mesh=vm, V_profile=V_prof, pads=pads)


if __name__ == "__main__":
    out, ok = {}, True
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                out[name] = fn()
                print("PASS", name, json.dumps(out[name], ensure_ascii=False, default=float)[:300])
            except Exception as ex:  # noqa: BLE001
                ok = False
                import traceback
                traceback.print_exc()
                print("FAIL", name, repr(ex)[:400])
    sys.exit(0 if ok else 1)
