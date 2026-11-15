#!/usr/bin/env python3
"""KV-3: оптимизация топологии диска турбины ступени 1 и колеса компрессора ступени 2 (осесимметрично).

Исходные профили, поля температур, нагрузка лопаток — routes/gtu/frame3/frame3_fea_results.json (литой / кованый
профиль типоразмера № 3, проверенный CalculiX). Варианты (README узла):
  A «живой профиль»  — сплошной металл или пустота в огибающей: диск турбины — в осевом габарите ступицы
                       исходного диска (прямоугольник t_max), колесо — в огибающей исходного профиля (тракт и полость);
  B «полый печатный» — огибающая исходного профиля, обшивка t_skin сплошная, внутри — гироид переменной плотности.
Расчётные случаи: LC1 — 100 % (14 950 об/мин) + нагрев; LC2 — 120 % + нагрев.
Допускаемые: диск турбины IN-713C — LC1 длительная прочность σ_LM(T, 25 000 ч), LC2 — местный σв(T);
колесо Ti-6Al-4V — LC1 k·σ0,2(T)/1,5, LC2 k·σ0,2(T) (k = 0,9 печать + ГИП против поковки — оценка, как
funnel/make_model.u_limit_comp; исходное колесо — k = 1).

Запуск:  python3 disk_topo.py --all [--h-t 1.0e-3 --h-c 1.5e-3 --workers 4]
Выход:   disk_topo_results.json (сводка), поля — <scratchpad>/organic/disk_*.npz
"""
import argparse
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import topo_axi as TA  # noqa: E402
import lattice as LT  # noqa: E402

ROOT = os.path.normpath(os.path.join(HERE, "..", "..", "..", ".."))
F3 = os.path.join(ROOT, "routes", "gtu", "frame3", "frame3_fea_results.json")
OUTDIR = os.environ.get("ORGANIC_OUT", "/tmp/claude-0/-home-user/a62745e4-d2fd-56e0-974f-12746b1a2ae0/scratchpad/organic")
RESULTS = os.path.join(HERE, "disk_topo_results.json")
RPM0 = 14950.0
K_AM = 0.9          # оценка: печать + ГИП против поковки (funnel/make_model.u_limit_comp: 0,9 × 620 м/с)
T_SKIN = 2.0e-3     # обшивка полого диска, м (оценка: ≥ 4 слоёв дорожки LPBF и запас на обработку)
CREEP_H = 25000.0


def frame3():
    with open(F3, encoding="utf-8") as f:
        return json.load(f)


def poly_profile(r, ztop, zbot):
    """многоугольник сечения по профилю: верх слева направо, низ справа налево (ось r = 0 — ребро)"""
    top = np.c_[r, ztop]; bot = np.c_[r[::-1], zbot[::-1]]
    return np.vstack([top, bot])


def inside(poly, x, y):
    from matplotlib.path import Path
    return Path(poly).contains_points(np.c_[x, y], radius=1e-12)


def dist_to_lines(px, py, lines):
    import shapely
    from shapely.geometry import MultiLineString
    ml = MultiLineString([np.asarray(l) for l in lines])
    pts = shapely.points(px, py)
    return shapely.distance(pts, ml)


# =========================================================================== постановки
def setup(part, variant, h, objective="mass"):
    d = frame3()
    P = d["parts"][part]
    meta = dict(part=part, variant=variant, h=h)
    if part == "turb1":
        mat = TA.DK.load_material("in713c")
        g = P["geometry"]; pi = P["profile_info"]
        prof = P["profile"]
        r_live, w = g["r_live_rim_m"], g["w_rim_m"]
        r_b, t_max = pi["r_b_m"], pi["t_max_m"]
        Tv = P["variants"]["cooled_750"]["T_C"]
        T_of_r = lambda r: np.interp(r, Tv["r"], Tv["T"])  # noqa: E731
        p_rim = P["variants"]["cooled_750"]["cases"]["100"]["rim_p_MPa"] * 1e6
        r = np.asarray(prof["r"]); t = np.asarray(prof["t"])
        base_poly = poly_profile(r, t / 2, -t / 2)
        if variant == "A":
            env = np.array([[0, -t_max / 2], [r_b, -t_max / 2], [r_b, -w / 2], [r_live, -w / 2], [r_live, w / 2],
                            [r_b, w / 2], [r_b, t_max / 2], [0, t_max / 2]])
        else:
            env = base_poly
        nr = int(round(r_live / h)); nz = int(round(t_max / h)); nz += nz % 2
        z0, z1 = -t_max / 2, t_max / 2
        r1 = r_live
        meta.update(material="in713c", T_field="cooled_750 (обод 750 °C, ступица на 200 K холоднее — frame3)",
                    r_live=r_live, w_rim=w, r_b=r_b, t_max=t_max, p_rim_MPa=p_rim / 1e6)
    elif part == "comp2":
        mat = TA.DK.load_material("ti64")
        s = P["variants"]["solid"]
        prof = s["profile"]
        r = np.asarray(prof["r"]); t = np.asarray(prof["t"]); zc = np.asarray(prof["zc"]); hb = np.asarray(prof["hb"])
        ztop, zbot = zc + t / 2, zc - t / 2
        base_poly = poly_profile(r, ztop, zbot)
        env = base_poly
        Tv = s["T_C"]
        T_of_r = lambda rr: np.interp(rr, Tv["r"], Tv["T"])  # noqa: E731
        r1 = float(r[-1]); z0, z1 = float(zbot.min()), float(ztop.max())
        nr = int(round(r1 / h)); nz = int(round((z1 - z0) / h))
        meta.update(material="ti64", T_field="frame3 comp2 solid: T(r) 204 → 411 °C", R2=r1,
                    blade_fill=P["blade_fill"], boss_r=P["boss_r_m"], boss_t=P["boss_t_m"])
        p_rim = 0.0
    else:
        raise ValueError(part)
    dr, dz = r1 / nr, (z1 - z0) / nz
    I, J = np.meshgrid(np.arange(nr), np.arange(nz), indexing="ij")
    rc = (I + 0.5) * dr; zcg = z0 + (J + 0.5) * dz
    active = inside(env, rc.ravel(), zcg.ravel()).reshape(nr, nz)
    in_base = inside(base_poly, rc.ravel(), zcg.ravel()).reshape(nr, nz)
    group = np.ones((nr, nz), int)
    if part == "turb1":
        # пассивен только слой под дном пазов (4 мм): переход «полотно — обод» проектируется (при пассивном ободе от r_b
        # оптимизация упирала тонкое полотно в острый входящий угол — особенность напряжений в гладкой сетке ccx)
        passive = (rc >= meta["r_live"] - 4.0e-3)
    else:
        R2 = meta["R2"]
        hub_line = np.c_[r, ztop]
        dhub = dist_to_lines(rc.ravel(), zcg.ravel(), [hub_line]).reshape(nr, nz)
        passive = (dhub <= 3.0e-3) | (rc >= R2 - 10e-3) | ((rc <= meta["boss_r"]) & (zcg <= z0 + 10e-3))
        meta["skin_hub_m"] = 3.0e-3
    if variant == "B":
        # обшивка по всей границе огибающей, кроме оси r = 0
        bl = [env[1:-1]] if part == "comp2" else [base_poly[1:-1]]
        dsk = dist_to_lines(rc.ravel(), zcg.ravel(), bl).reshape(nr, nz)
        passive = passive | (dsk <= T_SKIN)
        group[:] = 2
        meta["t_skin_m"] = T_SKIN
    group[passive] = 0
    model = TA.AxiModel(0.0, r1, z0, z1, nr, nz, active, group, mat, T_of_r)
    model.rpm0 = RPM0
    model.fix(model.nodes_where(lambda rr, zz: rr < 1e-12), 0)
    if part == "turb1":
        ax = model.nodes_where(lambda rr, zz: rr < 1e-12)
        model.fix([ax[np.argmin(np.abs(model.node_z[ax]))]], 1)
        # давление лопаток на грани r_live, |z| ≤ w/2; полная сила — как в disk.py: p·2π·r·w
        F = model.add_edge_traction_r(meta["r_live"], -meta["w_rim"] / 2 - 1e-9, meta["w_rim"] / 2 + 1e-9, p_rim)
        F_ex = p_rim * 2 * math.pi * meta["r_live"] * meta["w_rim"]
        model.f_fix0 *= F_ex / F
        meta["rim_force_N"] = F_ex
    else:
        ax = model.nodes_where(lambda rr, zz: rr < 1e-12)
        model.fix([ax[np.argmin(model.node_z[ax])]], 1)
        # масса лопаток (доля φ объёма над ступицей высотой hb) — радиальные силы на узлах поверхности ступицы
        w2 = (2 * math.pi * RPM0 / 60) ** 2
        phi = meta["blade_fill"]; rho = model.rho
        f_r = lambda rr: phi * rho * w2 * np.interp(rr, r, hb) * 2 * math.pi * rr ** 2  # noqa: E731
        Ftot = 0.0
        for i in range(nr):
            js = np.nonzero(active[i])[0]
            if len(js) == 0:
                continue
            jt = js.max()
            ra, rb_ = i * dr, (i + 1) * dr
            rr = np.linspace(ra, rb_, 7)
            Fi = float(np.trapezoid(f_r(rr), rr))
            if Fi <= 0:
                continue
            n_a = model.nodes_where(lambda x, y: (np.abs(x - ra) < 1e-9) & (np.abs(y - (z0 + (jt + 1) * dz)) < 1e-9))
            n_b = model.nodes_where(lambda x, y: (np.abs(x - rb_) < 1e-9) & (np.abs(y - (z0 + (jt + 1) * dz)) < 1e-9))
            model.add_force(n_a, fr=np.full(len(n_a), Fi / 2)); model.add_force(n_b, fr=np.full(len(n_b), Fi / 2))
            Ftot += Fi
        meta["blade_force_N"] = Ftot
        meta["blade_f_r"] = dict(phi=phi, note="f(r) = φ·ρ·ω²·hb(r)·2π·r² Н на метр радиуса (как плотность столбца disk.py)")
    # допускаемые
    mt = model.mt
    Te = model.Te
    if part == "turb1":
        uts = mt["uts"](Te) * 1e6
        a1 = np.minimum(mt["slm"](Te, CREEP_H) * 1e6, uts); a2 = uts
        k_am = 1.0
        crit = ["LC1: σ_экв ≤ σ_LM(T, 25 000 ч) (типовые свойства, не выше σв)", "LC2 (120 %): σ_экв ≤ σв(T) местно"]
    else:
        k_am = 1.0 if variant == "base" else K_AM
        ys = mt["ys"](Te) * 1e6
        a1 = k_am * ys / 1.5; a2 = k_am * ys
        crit = [f"LC1: σ_экв ≤ {k_am}·σ0,2(T)/1,5", f"LC2 (120 %): σ_экв ≤ {k_am}·σ0,2(T)"]
    cases = [dict(rpm=RPM0, thermal=True, allow=a1), dict(rpm=1.2 * RPM0, thermal=True, allow=a2)]
    meta.update(criteria=crit, k_am=k_am, nr=nr, nz=nz, elements=int(model.ne), dofs=int(model.ndof),
                envelope=env.tolist(), base_poly=base_poly.tolist())
    prm = dict(q=5.0, Emin=1e-6, qs=0.5, P=10.0, move=0.1)
    if variant == "B":
        prm["lattice"] = LT.default_law()
    design = model.group != 0
    in_base_e = in_base[model.ei, model.ej]
    return model, design, cases, prm, meta, in_base_e


def baseline(part, h):
    """исходный профиль в той же сетке: активны элементы внутри исходного профиля, все пассивные"""
    model, design, cases, prm, meta, in_base_e = setup(part, "A", h)
    xb = np.where(in_base_e, 1.0, 0.0)
    xb[model.group == 0] = 1.0
    # пассивный обод/обшивка вне исходного профиля не нужны: оставляем только исходный профиль
    keep = in_base_e
    if part == "comp2":
        # исходное кованое колесо: k = 1
        mt = model.mt; Te = model.Te; ys = mt["ys"](Te) * 1e6
        cases_k1 = [dict(rpm=RPM0, thermal=True, allow=ys / 1.5), dict(rpm=1.2 * RPM0, thermal=True, allow=ys)]
    else:
        cases_k1 = cases
    xb = np.where(keep, 1.0, 0.0)
    model.group = np.where(keep, 0, 1)        # пустота вне профиля — «solid» с x̄ = 0
    sol = model.solve(xb, prm, cases_k1, need_sens=False)
    em = np.where(keep, 1.0, 0.0)
    out = summarize(model, xb, em, sol, cases_k1, keep, meta)
    out["note"] = "исходный профиль frame3 в той же сетке (элементы с центром внутри профиля); тот же расчёт, что и проекты"
    return out, model, xb


def summarize(model, xb, em, sol, cases, solid_mask, meta):
    mt = model.mt
    res = dict(mass_kg=model.mass(em), J_kgm2=model.inertia(em))
    res["cases"] = []
    for k, c in enumerate(sol["cases"]):
        msk = solid_mask
        u = np.where(msk, c["u"], 0.0)
        e = int(np.argmax(u))
        vm = np.where(msk, c["vm"] * (em if False else 1.0), 0.0)
        rec = dict(rpm=cases[k]["rpm"], u_max=float(u.max()), margin=float(1 / max(u.max(), 1e-12)),
                   at=dict(r_mm=float(model.rc[e] * 1e3), z_mm=float(model.zc[e] * 1e3), T_C=float(model.Te[e])),
                   vm_max_MPa=float(np.max(np.where(solid_mask & (model.group != 2), c["vm"], 0.0)) / 1e6))
        sig = c["sig"]
        solid = solid_mask & (model.group != 2)
        if solid.any():
            rec["hoop_mean_solid_MPa"] = float(np.mean(sig[solid, 2]) / 1e6)
        res["cases"].append(rec)
    # запас по разрушению (100 %): ∫σв(T)·S dA / ∫σθ dA по сечению; решётка — σθ макро, σв·S(ρ*)
    c0 = sol["cases"][0]
    Sfac = np.ones(model.ne)
    if (model.group == 2).any():
        lat = LT.default_law()
        rs = lat["rho_min"] + (1 - lat["rho_min"]) * xb
        Sfac = np.where(model.group == 2, lat["S"](rs), 1.0)
        Efac = np.where(model.group == 2, lat["E"](rs), 1.0)
    else:
        Efac = np.ones(model.ne)
    hoop = c0["sig"][:, 2] * Efac
    uts = mt["uts"](model.Te) * 1e6 * Sfac
    m = solid_mask
    res["burst_margin"] = float(math.sqrt(np.sum(uts[m]) / max(np.sum(hoop[m]), 1e-9)))
    res["hoop_mean_MPa"] = float(np.mean(hoop[m]) / 1e6)
    if "slm" in mt:
        lmp = model.mat["_lmp"]
        vm1 = np.where(m, c0["vm"] * Efac, 0.0) / 1e6
        life, ex = lmp.life_h_v(np.where(vm1 > 0, vm1 / np.where(model.group == 2, Sfac, 1.0), 0.0), model.Te)
        life = np.where(m, life, np.inf)
        e = int(np.argmin(life))
        res["creep_life_min_h"] = float(life[e])
        res["creep_life_at"] = dict(r_mm=float(model.rc[e] * 1e3), z_mm=float(model.zc[e] * 1e3), T_C=float(model.Te[e]),
                                    extrapolated=bool(ex[e]))
    return res


# =========================================================================== прогон варианта
def run_spec(spec):
    t0 = time.time()
    part, variant, h = spec["part"], spec["variant"], spec["h"]
    model, design, cases, prm, meta, in_base_e = setup(part, variant, h)
    tgt = spec.get("u_target", [1.0, 1.0])
    x0 = np.ones(model.ne) if variant == "A" else np.full(model.ne, spec.get("x0", 0.5))
    logf = open(os.path.join(OUTDIR, f"log_{spec['name']}.txt"), "w")

    def log(rec):
        logf.write(json.dumps(rec) + "\n"); logf.flush()

    stages = spec.get("stages")
    if stages is None:
        stages = ([dict(beta=1, it=50, move=0.1), dict(beta=2, it=30, move=0.08), dict(beta=4, it=30, move=0.06),
                   dict(beta=8, it=30, move=0.05), dict(beta=16, it=50, move=0.03)]
                  if variant == "A" else [dict(beta=0, it=150, move=0.08)])
    R = spec.get("R_filter", 2.5 * h)
    # цель с запасом 3 % на переход от проецированного к двоичному проекту; затем досборка (до двух раз), если
    # двоичный проект всё же нарушает ограничение: цель уменьшается пропорционально превышению
    tgt_eff = [t * 0.97 for t in tgt]
    opt = TA.optimize(model, x0, design, cases, prm, objective=spec["objective"], R_filter=R, stages=stages, log=log,
                      u_target=tgt_eff)
    hist = list(opt["hist"])
    repairs = []

    def binarize(xb):
        if variant == "A":
            xe = np.where(xb >= 0.5, 1.0, 0.0); xe[model.group == 0] = 1.0
            return xe, xe > 0.5
        xe = xb.copy(); xe[model.group == 0] = 1.0
        return xe, np.ones(model.ne, bool)

    for rep in range(3):
        xe, solid = binarize(opt["xb"])
        sol = model.solve(xe, prm, cases, need_sens=False)
        ub = [float(np.max(np.where(solid, c["u"], 0.0))) for c in sol["cases"]]
        repairs.append(dict(round=rep, u_binary=ub, target_eff=list(tgt_eff)))
        if all(ub[k] <= tgt[k] * 1.005 for k in range(len(tgt))) or rep == 2:
            break
        tgt_eff = [tgt_eff[k] * min(1.0, tgt[k] / ub[k]) * 0.99 for k in range(len(tgt))]
        last = dict(stages[-1]); last["it"] = 40
        opt = TA.optimize(model, opt["x"], design, cases, prm, objective=spec["objective"], R_filter=R, stages=[last],
                          log=log, u_target=tgt_eff)
        hist += opt["hist"]
    logf.close()
    opt["hist"] = hist
    xb = opt["xb"]
    xe, solid = binarize(xb)
    sol = model.solve(xe, prm, cases, need_sens=False)
    eE, em, es, *_ = model.interp(xe, prm)
    out = summarize(model, xe, em, sol, cases, solid, meta)
    out.update(name=spec["name"], part=part, variant=variant, objective=spec["objective"], u_target=tgt,
               repairs=repairs, constraint_ok=bool(all(out["cases"][k]["u_max"] <= tgt[k] * 1.005 for k in range(len(tgt)))),
               iterations=len(opt["hist"]), seconds=round(time.time() - t0, 1),
               grey_final=float(opt["hist"][-1]["grey"]), R_filter_m=R, meta={k: v for k, v in meta.items()
                                                                              if k not in ("envelope", "base_poly")})
    if variant == "B":
        lat = prm["lattice"]
        rs = lat["rho_min"] + (1 - lat["rho_min"]) * xe
        g2 = model.group == 2
        out["lattice"] = dict(rho_min=lat["rho_min"], volume_share_of_envelope=float(np.sum(model.Ve[g2]) / np.sum(model.Ve)),
                              rho_mean=float(np.sum(rs[g2] * model.Ve[g2]) / max(np.sum(model.Ve[g2]), 1e-12)),
                              share_at_min=float(np.mean(rs[g2] <= lat["rho_min"] + 0.02)),
                              share_solid_like=float(np.mean(rs[g2] >= 0.9)))
    np.savez_compressed(os.path.join(OUTDIR, f"disk_{spec['name']}.npz"), xb=xb, xe=xe, ei=model.ei, ej=model.ej,
                        group=model.group, nr=model.nr, nz=model.nz, r1=model.r1, z0=model.z0, z1=model.z1,
                        hist=json.dumps(opt["hist"]), env=np.asarray(meta["envelope"]), base=np.asarray(meta["base_poly"]))
    return out


def specs(h_t, h_c, u_eq):
    """eq — не хуже исходного ни по одному случаю (цель по загрузке = загрузка исходного профиля в той же сетке);
    lim — до предела материала (u ≤ 1)."""
    ut, uc = u_eq["turb1"], u_eq["comp2"]
    S = [dict(name="turb1_A_mass_eq", part="turb1", variant="A", h=h_t, objective="mass", u_target=ut),
         dict(name="turb1_A_mass_lim", part="turb1", variant="A", h=h_t, objective="mass"),
         dict(name="turb1_A_J_eq", part="turb1", variant="A", h=h_t, objective="inertia", u_target=ut),
         dict(name="comp2_A_mass_eq", part="comp2", variant="A", h=h_c, objective="mass", u_target=uc),
         dict(name="comp2_A_mass_lim", part="comp2", variant="A", h=h_c, objective="mass"),
         dict(name="comp2_A_J_eq", part="comp2", variant="A", h=h_c, objective="inertia", u_target=uc),
         dict(name="turb1_B_mass_eq", part="turb1", variant="B", h=h_t, objective="mass", u_target=ut),
         dict(name="comp2_B_mass_eq", part="comp2", variant="B", h=h_c, objective="mass", u_target=uc),
         dict(name="comp2_B_J_eq", part="comp2", variant="B", h=h_c, objective="inertia", u_target=uc)]
    return S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("--h-t", type=float, default=1.0e-3)
    ap.add_argument("--h-c", type=float, default=1.5e-3)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--results", default=RESULTS)
    a = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    t0 = time.time()
    base = {}
    for part, h in (("turb1", a.h_t), ("comp2", a.h_c)):
        b, _, _ = baseline(part, h)
        base[part] = b
        print(part, "baseline", json.dumps(b, default=float)[:400], flush=True)
    u_eq = {p: [c["u_max"] for c in base[p]["cases"]] for p in base}
    S = specs(a.h_t, a.h_c, u_eq)
    if a.only:
        S = [s for s in S if s["name"] in a.only.split(",")]
    with Pool(a.workers) as pool:
        res = pool.map(run_spec, S, chunksize=1)
    prev = {}
    if os.path.exists(a.results) and a.only:
        with open(a.results, encoding="utf-8") as f:
            prev = {r["name"]: r for r in json.load(f).get("runs", [])}
    for r in res:
        prev[r["name"]] = r
    out = dict(title="KV-3: оптимизация топологии диска турбины ступени 1 и колеса компрессора ступени 2",
               method="Свой осесимметричный КЭ (Q4, 2×2 Гаусса) + плотностный метод (RAMP, фильтр плотности, проекция "
                      "Хевисайда с продолжением по β) + MMA; ограничения по напряжению p-нормой (P = 10) с адаптивной нормировкой "
                      "и сопряжёнными производными; вариант B — гироид переменной плотности по гомогенизации lattice.py",
               software=dict(python=sys.version.split()[0], numpy=np.__version__, scipy=__import__("scipy").__version__),
               date="2026-11-14",
               inputs=[dict(file="routes/gtu/frame3/frame3_fea_results.json",
                            field="parts.turb1: profile, profile_info, geometry, variants.cooled_750 (T_C, rim_p_MPa); "
                                  "parts.comp2: variants.solid (profile r/t/zc/hb, T_C), blade_fill, boss_r_m"),
                       dict(file="routes/gtu/platform/fea/materials.json", field="in713c, ti64 (через disk.load_material)"),
                       dict(file="routes/gtu/organic/structure/lattice_results.json", field="law (E*/E, S от ρ*)")],
               assumptions=dict(K_AM=K_AM, T_SKIN_m=T_SKIN, creep_hours=CREEP_H, rpm=[RPM0, 1.2 * RPM0],
                                h_turb_m=a.h_t, h_comp_m=a.h_c,
                                note="k_AM и t_skin — оценка; поля температур — frame3 (оценка облика)"),
               baseline=base, runs=sorted(prev.values(), key=lambda r: r["name"]),
               seconds=round(time.time() - t0, 1))
    with open(a.results, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1, default=float)
    print("->", a.results)


if __name__ == "__main__":
    main()
