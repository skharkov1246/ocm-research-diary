"""Массовая 3D-прочность: критерии по выходу CalculiX - запасы, ресурс, диаграмма Кэмпбелла, вердикт детали.

Статика (100, 115, 120 %):
  * Мизес в узлах (.frd, экстраполяция ccx) - пик по классам (ступица, лопатка, замок); у цельного колеса узлы ближе
    2 r_отв к технологическому отверстию исключаются (концентрация ×2 ложная). У корня лопатки галтели нет - пик в узле
    зависит от сетки, поэтому в критериях - пик, осреднённый по объёму шара радиусом R_AV_H размеров элемента лопатки
    (оценка; сходимость по двум сеткам - режим conv);
  * σ/σ0,2(T) - осреднённый пик и узловой (для сведения);
  * запас по разрушающей частоте - критерий средних окружных напряжений (Робинсон, как disk.py) по ступице сектора:
    n_разр/n = sqrt(k·∫σв(T)/r dV / ∫σθ/r dV) (интеграл по сектору с весом 1/r - площадь меридиональных сечений,
    осреднённая по углу), k = 1 - верхняя оценка;
  * ползучесть - Ларсон - Миллер (materials3d.Material.rupture_h) на 100 %: по осреднённому пику, по номинальному
    напряжению корня пера (осевые ступени), по среднему окружному ступицы при средней и наибольшей температуре; узловой
    минимум - для сведения (перераспределение ползучестью не учтено);
  * МЦУ пуск - останов (0 -> 100 % -> 0): размах осреднённого пика Мизеса; при размахе больше 2σ0,2 - поправка Нейбера
    для идеальной пластичности Δε = Δσ²/(2σ0,2·E); число циклов - materials3d.Material.lcf_cycles (оценка).
Частоты (0 и рабочие частоты с преднагружением): семейства по ND, для каждой формы k - Саутвелл f² = a + b·(n/60)²
(наименьшие квадраты по рассчитанным частотам вращения), пересечения с f = EO·n/60 для ND = h(EO) (h = EO mod Z, отражённое
к Z/2), запас по частоте в рабочих точках.
"""
import math

import numpy as np

CRIT = dict(
    burst_120_min=1.0,          # запас по разрушающей частоте при 120 % (как правило уровня 1 и merge_l1.FEA_BURST_MIN)
    yield_ratio_max=1.0,        # σ_экв/σ0,2(T) при 100 % по осреднённому пику (оценка)
    life_h_min=25000.0,         # ресурс по длительной прочности, ч (strength_l1.RUPTURE_H, оценка)
    lcf_starts=3000,            # требуемое число пусков (оценка: пуск ежедневно 8 лет)
    lcf_factor=3.0,             # запас по числу циклов к оценке универсальных наклонов (оценка)
    k_modes_hard=3,             # пересечение первых трёх форм семейства в рабочем диапазоне - отказ
    k_modes_soft=6,             # формы 4-6 - условно
    margin_min=0.05,            # запас по частоте в рабочих точках (оценка; нормы заказчика не прочитаны)
    pf_max=1e-4,                # вероятность разрушения керамического колеса (strength_l1.PF_MAX)
)
CRIT_BASIS = ("запас при 120 % >= 1 - правило уровня 1; σ/σ0,2 <= 1, 25 000 ч, 3000 пусков × 3, отстройка 5 % и формы 1-3 - оценка "
              "(нормы заказчика, API 616, 14 CFR 33 в сессии не прочитаны)")
R_AV_H = 1.5                    # радиус осреднения пика - в размерах элемента лопатки (оценка)
HARD_EO = ("ngv", "ngv2", "inj", "igv", "igv2", "low")


def elem_geometry(mesh):
    X = mesh.nodes
    E = mesh.elems - 1
    a, b, c, d = (X[E[:, i]] for i in range(4))
    vol = np.abs(np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a)) / 6.0
    cen = X[E].mean(axis=1)
    return vol, cen


def align(mesh, ids, vals):
    """значения .frd по номерам узлов -> массив в порядке mesh.node_ids (отсутствующие - nan)"""
    out = np.full((len(mesh.node_ids), vals.shape[1]), np.nan)
    pos = np.searchsorted(mesh.node_ids, ids)
    ok = (pos < len(mesh.node_ids)) & (mesh.node_ids[np.minimum(pos, len(mesh.node_ids) - 1)] == ids)
    out[pos[ok]] = vals[ok]
    return out


def von_mises(S):
    sxx, syy, szz, sxy, syz, szx = (S[:, k] for k in range(6))
    return np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3 * (sxy ** 2 + syz ** 2 + szx ** 2))


def hoop(S, xyz):
    th = np.arctan2(xyz[:, 1], xyz[:, 0])
    c, s = np.cos(th), np.sin(th)
    return S[:, 0] * s * s + S[:, 1] * c * c - 2 * S[:, 3] * s * c


def s1_max(S):
    M = np.zeros((len(S), 3, 3))
    M[:, 0, 0], M[:, 1, 1], M[:, 2, 2] = S[:, 0], S[:, 1], S[:, 2]
    M[:, 0, 1] = M[:, 1, 0] = S[:, 3]
    M[:, 1, 2] = M[:, 2, 1] = S[:, 4]
    M[:, 0, 2] = M[:, 2, 0] = S[:, 5]
    return np.linalg.eigvalsh(M)[:, -1]


def _avg_at(point, cen, vol, val, mask, R):
    d = np.linalg.norm(cen - point, axis=1)
    sel = mask & (d <= R)
    if not sel.any():
        sel = mask & (d <= d[mask].min() + 1e-12)
    return float(np.sum(val[sel] * vol[sel]) / np.sum(vol[sel])), int(sel.sum())


def static_eval(mesh, mats, T_nodes, S_nodes, f_rpm, excl_r=None):
    """критерии одного шага статики. mats: {класс: Material}; S_nodes (N,6) Па в порядке mesh.node_ids"""
    vol, cen = elem_geometry(mesh)
    rcen = np.hypot(cen[:, 0], cen[:, 1])
    rn = np.hypot(mesh.nodes[:, 0], mesh.nodes[:, 1])
    E = mesh.elems - 1
    vm_n = von_mises(S_nodes) / 1e6
    vm_e = np.nanmean(vm_n[E], axis=1)
    T_e = T_nodes[E].mean(axis=1)
    h = mesh.meta.get("h_blade") or 0.005
    R_av = R_AV_H * h
    node_ok = np.ones(len(rn), bool) if excl_r is None else rn > excl_r
    elem_ok = np.ones(len(E), bool) if excl_r is None else rcen > excl_r
    # класс узла: по элементам
    ncls = np.full(len(rn), -1)
    for c in sorted(set(int(x) for x in mesh.ecls)):
        ncls[np.unique(E[mesh.ecls == c])] = c
    out = dict(rpm_frac=f_rpm, R_av_m=R_av, regions={})
    worst_yr, worst_yr_node = 0.0, 0.0
    for c, mat in mats.items():
        nm = (ncls == c) & node_ok & np.isfinite(vm_n)
        em = (mesh.ecls == c) & elem_ok
        if not nm.any() or not em.any():
            continue
        i = int(np.nanargmax(np.where(nm, vm_n, -np.inf)))
        avg, n_av = _avg_at(mesh.nodes[i], cen, vol, vm_e, em, R_av)
        Ti = float(T_nodes[i])
        ys = mat.ys(Ti)
        # σ/σ0,2 - по узлам (наибольшее отношение) и осреднённый пик в точке наибольшего отношения
        ys_n = np.array([mat.ys(t) or np.nan for t in T_nodes]) if mat.ys20 or mat.ys_tab else None
        reg = dict(material=mat.alloy, vm_peak_node_MPa=float(vm_n[i]), at=dict(r=float(rn[i]), z=float(mesh.nodes[i, 2]), T_C=Ti),
                   vm_peak_avg_MPa=avg, n_avg_elems=n_av)
        if ys_n is not None:
            ratio = np.where(nm, vm_n / ys_n, -np.inf)
            k = int(np.nanargmax(ratio))
            avg_k, _ = _avg_at(mesh.nodes[k], cen, vol, vm_e, em, R_av)
            yr_avg = avg_k / ys_n[k]
            reg.update(yield_ratio_node=float(ratio[k]), yield_ratio_avg=float(yr_avg), ys_MPa_at=float(ys_n[k]),
                       yield_at=dict(r=float(rn[k]), z=float(mesh.nodes[k, 2]), T_C=float(T_nodes[k])), vm_avg_at_yield_MPa=avg_k)
            worst_yr, worst_yr_node = max(worst_yr, yr_avg), max(worst_yr_node, float(ratio[k]))
        reg["ys_basis"] = mat.ys_basis()
        out["regions"][CLS[c]] = reg
    out["yield_ratio_avg_max"] = worst_yr or None
    out["yield_ratio_node_max"] = worst_yr_node or None
    return out, dict(vol=vol, cen=cen, rcen=rcen, vm_e=vm_e, T_e=T_e, vm_n=vm_n, elem_ok=elem_ok)


CLS = {0: "hub", 1: "blade", 2: "att"}


def robinson(mesh, mat, S_nodes, T_nodes, k=1.0, excl_r=None):
    """запас по разрушающей частоте (средние окружные напряжения ступицы) на шаге с напряжениями S_nodes"""
    if not mat.uts_tab:
        return None
    vol, cen = elem_geometry(mesh)
    rc = np.hypot(cen[:, 0], cen[:, 1])
    E = mesh.elems - 1
    m = mesh.ecls == 0
    s_th = np.nanmean(hoop(S_nodes, mesh.nodes)[E], axis=1) / 1e6
    T_e = T_nodes[E].mean(axis=1)
    uts = np.array([mat.uts(t) for t in T_e[m]])
    w = vol[m] / rc[m]
    I_s, I_u, A = float(np.sum(s_th[m] * w)), float(np.sum(uts * w)), float(np.sum(w))
    if I_s <= 0:
        return dict(margin=None, note="средние окружные напряжения не положительны")
    T_mean = float(np.sum(T_e[m] * w) / A)
    return dict(margin=math.sqrt(k * I_u / I_s), hoop_mean_MPa=I_s / A, uts_mean_MPa=I_u / A, T_mean_C=T_mean,
                T_hot_C=float(T_e[m].max()), k=k)


def creep_eval(mesh, mats, ev, T_nodes, part_meta):
    """ресурс по длительной прочности на 100 %: осреднённый пик по классам, корень пера (осевые), минимум по элементам"""
    out = {}
    vol, cen, vm_e, T_e, ok = ev["vol"], ev["cen"], ev["vm_e"], ev["T_e"], ev["elem_ok"]
    R_av = R_AV_H * (mesh.meta.get("h_blade") or 0.005)
    for c, mat in mats.items():
        em = (mesh.ecls == c) & ok
        if not mat.creep or not em.any():
            continue
        lives = np.array([mat.rupture_h(s, t)[0] for s, t in zip(vm_e[em], T_e[em])], dtype=float)
        idx = np.nonzero(em)[0]
        j = int(idx[int(np.argmin(lives))])
        avg, _ = _avg_at(cen[j], cen, vol, vm_e, em, R_av)
        life_avg, ext = mat.rupture_h(avg, float(T_e[j]))
        out[CLS[c]] = dict(life_min_elem_h=float(lives.min()), life_avg_h=float(life_avg), extrapolated=bool(ext),
                           at=dict(r=float(np.hypot(*cen[j, :2])), z=float(cen[j, 2]), T_C=float(T_e[j]), vm_avg_MPa=avg,
                                   vm_elem_MPa=float(vm_e[j])), basis=mat.creep_basis)
    if part_meta.get("kind") == "axial":
        r0, r1 = part_meta["r_root"], part_meta["r_tip"]
        mat = mats.get(1)
        band = (mesh.ecls == 1) & (np.hypot(cen[:, 0], cen[:, 1]) < r0 + 0.08 * (r1 - r0))
        if mat is not None and band.any():
            s = float(np.sum(vm_e[band] * vol[band]) / np.sum(vol[band]))
            Tb = float(np.sum(T_e[band] * vol[band]) / np.sum(vol[band]))
            life, ext = mat.rupture_h(s, Tb)
            out["blade_root"] = dict(vm_mean_MPa=s, T_C=Tb, life_h=None if life is None else float(life), extrapolated=ext,
                                     uts_MPa=mat.uts(Tb), overspeed_ratio=s * 1.2 ** 2 / mat.uts(Tb) if mat.uts(Tb) else None)
    return out


def lcf_eval(mats, static_100, T_by_region=None):
    """МЦУ пуск - останов по осреднённым пикам статики 100 %"""
    out = {}
    for c, mat in mats.items():
        reg = static_100["regions"].get(CLS[c])
        if not reg:
            continue
        ds = reg["vm_peak_avg_MPa"]
        T = reg["at"]["T_C"]
        E = mat.E(T)
        ys = mat.ys(T)
        if ys and ds > 2 * ys:
            de = (ds * 1e6) ** 2 / (2 * ys * 1e6 * E)
            neuber = True
        else:
            de = ds * 1e6 / E
            neuber = False
        N = mat.lcf_cycles(de, T)
        out[CLS[c]] = dict(d_sigma_MPa=ds, d_eps=de, neuber=neuber, T_C=T, cycles=None if N is None else float(N),
                           basis="оценка: универсальные наклоны, ψ принятое (materials3d.PSI)" if N is not None else "нет данных")
    return out


def weibull_eval(mesh, S_nodes, Z, m=15.0, mean_MPa=594.0, V0_mm3=1.0):
    """вероятность разрушения керамического колеса по Вейбуллу: интеграл по объёму всего колеса (Z секторов), σ1 > 0"""
    vol, cen = elem_geometry(mesh)
    E = mesh.elems - 1
    s1 = np.nanmean(s1_max(S_nodes)[E], axis=1) / 1e6
    s0 = mean_MPa / math.gamma(1 + 1 / m)
    I = float(np.sum((vol * 1e9 / V0_mm3) * (np.maximum(s1, 0.0) / s0) ** m)) * Z
    return dict(pf=1 - math.exp(-I), s1_max_MPa=float(np.nanmax(s1)), m=m, sigma0_MPa=s0, V0_mm3=V0_mm3,
                basis="Вейбулл: средняя прочность SN282A 594 МПа при 1200 °C (по выдержке поиска, strength_l1), m = 15 - оценка")


# ----------------------------------------------------------------------------- частоты
def alias(eo, Z):
    h = int(eo) % int(Z)
    return int(Z) - h if h > Z / 2 else h


def southwell(table, rpm, nd, k):
    """наименьшие квадраты f² = a + b·(n/60)² по рассчитанным частотам вращения для формы k семейства nd"""
    xs, ys = [], []
    for frac, fam in table:
        v = fam.get(nd) or []
        if len(v) > k:
            xs.append((frac * rpm / 60.0) ** 2)
            ys.append(v[k] ** 2)
    if len(xs) < 2:
        return None
    X = np.vstack([np.ones(len(xs)), xs]).T
    (a, b), *_ = np.linalg.lstsq(X, np.asarray(ys), rcond=None)
    pred = np.sqrt(np.maximum(a + b * np.asarray(xs), 0.0))
    res = float(np.max(np.abs(pred / np.sqrt(ys) - 1))) if len(xs) > 2 else 0.0
    return dict(a=float(a), b=float(b), max_rel_residual=res, n=len(xs))


def campbell(table, rpm, Z, eo_sources, op, free_stator=None, k_soft=None):
    """table: [(доля частоты, {nd: [f...]}), ...]. Пересечения форм семейства h(EO) с прямыми EO·n/60 до 1,25 n"""
    k_soft = k_soft or CRIT["k_modes_soft"]
    lo, hi = op["range"]
    fits = {}
    for nd in sorted({nd for _, fam in table for nd in fam}):
        for k in range(k_soft):
            sw = southwell(table, rpm, nd, k)
            if sw:
                fits[(nd, k)] = sw

    def f_of(nd, k, frac):
        sw = fits.get((nd, k))
        if not sw:
            return None
        return math.sqrt(max(sw["a"] + sw["b"] * (frac * rpm / 60.0) ** 2, 0.0))

    def crossings(eo):
        h = alias(eo, Z)
        out = []
        for k in range(k_soft):
            sw = fits.get((h, k))
            if not sw:
                continue
            den = eo ** 2 - sw["b"]
            if den <= 0 or sw["a"] <= 0:
                continue
            frac = 60.0 * math.sqrt(sw["a"] / den) / rpm
            if frac <= 1.25:
                out.append(dict(nd=h, mode=k + 1, rpm_frac=frac, f_Hz=eo * frac * rpm / 60.0, in_range=bool(lo <= frac <= hi),
                                in_info_range=bool(op.get("info_range") and op["info_range"][0] <= frac <= op["info_range"][1])))
        return out
    cr = []
    for s in eo_sources:
        for x in crossings(s["eo"]):
            cr.append(dict(x, eo=s["eo"], kind=s["kind"], source=s["source"]))
    # запас по частоте в рабочих точках (границы диапазона и 100 %)
    margins = []
    for frac in sorted({lo, 1.0, hi}):
        for s in eo_sources:
            h = alias(s["eo"], Z)
            fe = s["eo"] * frac * rpm / 60.0
            for k in range(CRIT["k_modes_hard"]):
                f = f_of(h, k, frac)
                if f:
                    margins.append(dict(rpm_frac=frac, eo=s["eo"], kind=s["kind"], nd=h, mode=k + 1, f_Hz=f, f_exc_Hz=fe,
                                        margin=abs(f - fe) / fe))
    mmin = min(margins, key=lambda x: x["margin"]) if margins else None
    free = None
    if free_stator:
        bad = []
        for zd in range(free_stator[0], free_stator[1] + 1):
            hit = [x for e in (zd, 2 * zd) for x in crossings(e) if x["in_range"] and x["mode"] <= CRIT["k_modes_hard"]]
            if hit:
                bad.append(zd)
        free = dict(range=list(free_stator), forbidden=bad, allowed=[z for z in range(free_stator[0], free_stator[1] + 1) if z not in bad])
    tab = [dict(rpm_frac=frac, nd={str(nd): [round(f, 2) for f in v] for nd, v in sorted(fam.items())}) for frac, fam in table]
    spin = None
    if len(table) >= 2:
        f0 = table[0][1]
        fN = table[-1][1]
        spin = {str(nd): (round(fN[nd][0] / f0[nd][0] - 1, 4) if fN.get(nd) and f0.get(nd) else None) for nd in sorted(f0)}
    sw_res = max((v["max_rel_residual"] for v in fits.values()), default=None)
    return dict(table=tab, crossings=cr, margin_min=mmin, free_stator=free, southwell_max_rel_residual=sw_res, spin_change_mode1=spin,
                op_range=[lo, hi], info_range=op.get("info_range"))


# ----------------------------------------------------------------------------- вердикт детали
def verdict(part_res):
    """отказ / условно / годен с причинами (CRIT)"""
    hard, soft = [], []
    b = part_res.get("burst") or {}
    if b.get("margin_120") is not None and b["margin_120"] < CRIT["burst_120_min"]:
        hard.append(f"запас по разрушающей частоте при 120 % {b['margin_120']:.3f} < {CRIT['burst_120_min']}")
    st = part_res.get("static") or {}
    s100 = (st.get("steps") or [{}])[0]
    yr = s100.get("yield_ratio_avg_max")
    if yr is not None and yr > CRIT["yield_ratio_max"]:
        hard.append(f"σ/σ0,2 при 100 % (осреднённый пик) {yr:.3f} > {CRIT['yield_ratio_max']}")
    yrn = s100.get("yield_ratio_node_max")
    if yrn is not None and yrn > CRIT["yield_ratio_max"] and not (yr is not None and yr > CRIT["yield_ratio_max"]):
        soft.append(f"σ/σ0,2 в узле {yrn:.3f} > 1 (осреднённый {yr:.3f}): местная пластичность у острого угла")
    cr = part_res.get("creep") or {}
    for reg, v in cr.items():
        if reg == "blade_root":
            if v.get("life_h") is not None and v["life_h"] < CRIT["life_h_min"]:
                hard.append(f"ползучесть корня пера: {v['life_h']:.0f} ч < {CRIT['life_h_min']:.0f}")
            if v.get("overspeed_ratio") and v["overspeed_ratio"] > 1.0:
                hard.append(f"корень пера при 120 %: σ/σв {v['overspeed_ratio']:.3f} > 1")
            continue
        if v.get("life_avg_h") is not None and v["life_avg_h"] < CRIT["life_h_min"]:
            hard.append(f"ползучесть ({reg}, осреднённый пик): {v['life_avg_h']:.0f} ч < {CRIT['life_h_min']:.0f}")
        elif v.get("life_min_elem_h") is not None and v["life_min_elem_h"] < CRIT["life_h_min"]:
            soft.append(f"ползучесть ({reg}) по элементу {v['life_min_elem_h']:.0f} ч < {CRIT['life_h_min']:.0f} (осреднённый - выше)")
    sm = (part_res.get("burst") or {}).get("section_life_h")
    if sm is not None and sm < CRIT["life_h_min"]:
        hard.append(f"ползучесть ступицы по среднему окружному при наибольшей температуре: {sm:.0f} ч < {CRIT['life_h_min']:.0f}")
    req = CRIT["lcf_starts"] * CRIT["lcf_factor"]
    for reg, v in (part_res.get("lcf") or {}).items():
        if v.get("cycles") is not None and v["cycles"] < req:
            soft.append(f"МЦУ ({reg}, оценка): {v['cycles']:.0f} циклов < {req:.0f}")
    wb = part_res.get("weibull")
    if wb and wb["pf"] > CRIT["pf_max"]:
        hard.append(f"вероятность разрушения керамики {wb['pf']:.2e} > {CRIT['pf_max']:.0e}")
    cb = part_res.get("campbell") or {}
    for x in cb.get("crossings") or []:
        if not x["in_range"]:
            continue
        txt = f"резонанс: {x['source']} (EO {x['eo']}), ND {x['nd']}, форма {x['mode']}, {x['rpm_frac'] * 100:.1f} % n"
        if x["kind"] in HARD_EO and x["mode"] <= CRIT["k_modes_hard"]:
            hard.append(txt)
        else:
            soft.append(txt)
    mm = cb.get("margin_min")
    if mm and mm["margin"] < CRIT["margin_min"] and not any("резонанс" in h for h in hard):
        soft.append(f"отстройка {mm['margin'] * 100:.1f} % < {CRIT['margin_min'] * 100:.0f} %: EO {mm['eo']} ({mm['kind']}), ND {mm['nd']}, "
                    f"форма {mm['mode']} при {mm['rpm_frac'] * 100:.0f} % n")
    fs = cb.get("free_stator")
    if fs and not fs["allowed"]:
        soft.append(f"нет числа лопаток диффузора {fs['range'][0]}-{fs['range'][1]} без резонанса")
    v = "отказ" if hard else ("условно" if soft else "годен")
    return v, hard, soft
