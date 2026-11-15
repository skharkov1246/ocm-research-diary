#!/usr/bin/env python3
"""Встроенный коллектор топливного газа KV-3: дерево каналов в стенке оболочки против кольцевого коллектора с трубками.

Физика: закон Мюррея (Murray 1926) и конструктальный подход (Bejan 1997) — ветвящаяся сеть, у которой при заданном
объёме каналов наименьшая работа прокачки. Для ламинарного течения D0³ = ΣDi³; здесь течение турбулентное (Re ~10⁶),
поэтому показатель ищется численно: минимум перепада давления дерева при заданном объёме с трением по Haaland.
Сравнение при одинаковом объёме каналов: (A) кольцевой коллектор, подвод в одной точке, 14 трубок к форсункам;
(B) симметричное дерево 1 → 2 → 7 (равные длины путей) в стенке оболочки. Неравномерность расхода по форсункам —
решением сети (уравнения узлов, Ньютон), разброс от допусков печати — методом Монте-Карло (±σ диаметра).
Газ — метан (Cantera gri30): плотность и вязкость при давлении коллектора. Расход — frame3 (Wf).
"""
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import minimize, fsolve

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", "..", ".."))
F3P = os.path.join(ROOT, "routes", "gtu", "frame3", "frame3_point_results.json")
RESULTS = os.path.join(HERE, "manifold_results.json")

N_INJ = 14                 # форсунок (funnel/parts_model.REF_BOM: n_inj 14)
DP_INJ_REL = 0.04          # перепад на форсунке / p3 (оценка: «газовая форсунка 3–5 % p3»; DLE — больше)
EPS_ROUGH = 15e-6          # шероховатость печатного канала после обработки, м (оценка; литые/печатные без обработки — 20–50 мкм)
L_PIGTAIL = 0.35           # длина трубки от кольца к форсунке, м (оценка)
TOL_D = 0.02               # σ диаметра канала печати, доля (оценка)


def gas_props(p, T=288.15):
    import cantera as ct
    g = ct.Solution("gri30.yaml")
    g.TPX = T, p, "CH4:1"
    return float(g.density), float(g.viscosity), ct.__version__


def f_haaland(Re, D):
    Re = max(Re, 1.0)
    if Re < 2300:
        return 64.0 / Re
    return (-1.8 * math.log10((EPS_ROUGH / D / 3.7) ** 1.11 + 6.9 / Re)) ** -2


def dp_pipe(m, D, L, rho, mu):
    """перепад, Па, для расхода m (кг/с, знак — направление)"""
    a = abs(m)
    if a < 1e-12:
        return 0.0
    A = math.pi * D * D / 4
    v = a / (rho * A)
    Re = rho * v * D / mu
    return math.copysign(f_haaland(Re, D) * L / D * rho * v * v / 2, m)


def tree_dp(D, L, m_tot, rho, mu, splits=(2, 7)):
    """симметричное дерево: уровень k — n_k каналов диаметра D[k], длины L[k]; перепад от входа до форсунки"""
    dp, m = 0.0, m_tot
    for k, n in enumerate((1,) + tuple(splits)):
        if k > 0:
            m = m / n
        dp += dp_pipe(m, D[k], L[k], rho, mu)
    return dp


def murray_tree(m_tot, rho, mu, L, V):
    """диаметры уровней дерева: минимум перепада при объёме V (численно)"""
    ns = [1, 2, 14]

    def vol(D):
        return sum(n * math.pi * d * d / 4 * l for n, d, l in zip(ns, D, L))

    def obj(x):
        D = np.exp(x)
        return tree_dp(D, L, m_tot, rho, mu) / 1e3

    cons = dict(type="eq", fun=lambda x: (vol(np.exp(x)) - V) / V)
    x0 = np.log([0.03, 0.022, 0.012])
    r = minimize(obj, x0, constraints=[cons], method="SLSQP", options=dict(ftol=1e-12, maxiter=500))
    D = np.exp(r.x)
    # показатель Мюррея: D_parent^k = n·D_child^k  ->  k = ln n / ln(D_p/D_c)
    k01 = math.log(2) / math.log(D[0] / D[1]); k12 = math.log(7) / math.log(D[1] / D[2])
    return D, k01, k12, r


def ring_network(D_ring, R_ring, d_pig, m_tot, rho, mu, p_in, K_inj, dvec=None):
    """кольцо, подвод в узле 0 на угле 0; форсунки — в N_INJ узлах кольца через трубку и сопло.
    Неизвестные: давления узлов кольца (N_INJ). Возврат: расходы форсунок."""
    N = N_INJ
    Ls = 2 * math.pi * R_ring / N
    dvec = np.ones(N) if dvec is None else dvec
    # подвод: узел 0 получает m_tot (подача — от источника, давление p_in не фиксируем: ищем давления узлов)

    def inj_flow(p_node, j):
        # трубка + сопло: p_node − p3 = dp_pipe(m) + K m²  ->  решаем относительно m (монотонно)
        dpt = p_node - p3_
        lo, hi = 0.0, 5 * m_tot / N
        for _ in range(60):
            mm = 0.5 * (lo + hi)
            if dp_pipe(mm, d_pig * dvec[j], L_PIGTAIL, rho, mu) + K_inj * mm * mm < dpt:
                lo = mm
            else:
                hi = mm
        return 0.5 * (lo + hi)

    def seg_flow(pa, pb):
        dpt = pa - pb
        s = 1.0 if dpt >= 0 else -1.0
        lo, hi = 0.0, 2 * m_tot
        for _ in range(60):
            mm = 0.5 * (lo + hi)
            if dp_pipe(mm, D_ring, Ls, rho, mu) < abs(dpt):
                lo = mm
            else:
                hi = mm
        return s * 0.5 * (lo + hi)

    p3_ = p_in

    def resid(x):
        p = p3_ + x * 1e5
        r = np.zeros(N)
        for j in range(N):
            r[j] = -inj_flow(p[j], j)
            r[j] -= seg_flow(p[j], p[(j + 1) % N])
            r[j] -= seg_flow(p[j], p[(j - 1) % N])
        r[0] += m_tot
        return r / m_tot

    dp0 = K_inj * (m_tot / N) ** 2 / 1e5
    x = fsolve(resid, np.full(N, dp0 * 1.05), xtol=1e-12)
    p = p3_ + x * 1e5
    q = np.array([inj_flow(p[j], j) for j in range(N)])
    return q, p


def tree_flows(D, L, m_tot, rho, mu, p3, K_inj, dvec=None):
    """дерево 1 → 2 → 7 с разбросом диаметров (dvec — множители для 1 + 2 + 14 каналов): расходы форсунок"""
    dvec = np.ones(17) if dvec is None else dvec
    D0, D1s, D2s = D[0] * dvec[0], D[1] * dvec[1:3], D[2] * dvec[3:]

    def leaf_m(dp_avail, d):
        lo, hi = 0.0, m_tot
        for _ in range(60):
            mm = 0.5 * (lo + hi)
            if dp_pipe(mm, d, L[2], rho, mu) + K_inj * mm * mm < dp_avail:
                lo = mm
            else:
                hi = mm
        return 0.5 * (lo + hi)

    def branch(dp_avail, d1, d2s):
        # подбор давления в узле разветвления: сумма листьев = расход ветви
        lo, hi = 0.0, dp_avail
        for _ in range(60):
            pn = 0.5 * (lo + hi)
            ml = sum(leaf_m(pn, d) for d in d2s)
            if dp_pipe(ml, d1, L[1], rho, mu) + pn < dp_avail:
                lo = pn
            else:
                hi = pn
        pn = 0.5 * (lo + hi)
        return [leaf_m(pn, d) for d in d2s]

    # давление на входе — подбор так, чтобы сумма = m_tot
    lo, hi = 0.0, 50e5
    for _ in range(50):
        P = 0.5 * (lo + hi)
        # в узле 0 после ствола
        mt = 0.0
        q = []
        # ствол: делим на две ветви
        lo2, hi2 = 0.0, P
        for _ in range(50):
            pn = 0.5 * (lo2 + hi2)
            qa = branch(pn, D1s[0], D2s[:7]); qb = branch(pn, D1s[1], D2s[7:])
            ms = sum(qa) + sum(qb)
            if dp_pipe(ms, D0, L[0], rho, mu) + pn < P:
                lo2 = pn
            else:
                hi2 = pn
        q = qa + qb
        if sum(q) < m_tot:
            lo = P
        else:
            hi = P
    return np.array(q), P


def main():
    with open(F3P, encoding="utf-8") as f:
        fp = json.load(f)
    Wf = float(fp["cycle"]["frame3"]["Wf"])
    p3 = 13.5 * 101325.0
    p_man = p3 * (1 + DP_INJ_REL)
    rho, mu, ctv = gas_props(p_man)
    m_inj = Wf / N_INJ
    K_inj = DP_INJ_REL * p3 / m_inj ** 2               # сопло форсунки: Δp = K·m²
    D2 = max(fp["stage1"]["D2_m"], fp["stage2"]["D2_m"])
    R_ring = 0.55 * 1.9 * D2 / 2                       # радиус коллектора у камеры (оценка: 0,55 радиуса оболочки)
    # (A) кольцо: D_ring и d_pig — типичные (оценка): 1" кольцо и 1/2" трубки
    D_ring, d_pig = 0.025, 0.012
    V_ring = math.pi * D_ring ** 2 / 4 * 2 * math.pi * R_ring + N_INJ * math.pi * d_pig ** 2 / 4 * L_PIGTAIL
    qA, pA = ring_network(D_ring, R_ring, d_pig, Wf, rho, mu, p3, K_inj)
    # (B) дерево той же ёмкости: длины уровней — ствол 0,3 м, ветви — четверть окружности, листья — до форсунок
    L_tree = [0.30, math.pi * R_ring / 2, 0.45 * R_ring]
    D_opt, k01, k12, r = murray_tree(Wf, rho, mu, L_tree, V_ring)
    qB, PB = tree_flows(D_opt, L_tree, Wf, rho, mu, p3, K_inj)
    # Монте-Карло допусков
    rng = np.random.default_rng(2026)
    devA, devB = [], []
    for _ in range(200):
        qa, _ = ring_network(D_ring, R_ring, d_pig, Wf, rho, mu, p3, K_inj, dvec=1 + TOL_D * rng.standard_normal(N_INJ))
        qb, _ = tree_flows(D_opt, L_tree, Wf, rho, mu, p3, K_inj, dvec=1 + TOL_D * rng.standard_normal(17))
        devA.append((qa.max() - qa.min()) / qa.mean()); devB.append((qb.max() - qb.min()) / qb.mean())
    # закон Мюррея ламинарный для сравнения (тот же объём): D0³ = 2D1³, D1³ = 7D2³
    res = dict(
        title="Встроенный коллектор топливного газа KV-3: дерево по закону Мюррея против кольца с трубками",
        method="сеть каналов: трение Haaland (ε = %.0f мкм), сопла форсунок Δp = K·m², уравнения узлов (fsolve / бисекции); "
               "дерево 1→2→7: диаметры — минимум перепада при объёме кольцевого варианта (SLSQP); Монте-Карло 200 выборок "
               "с σ диаметра %.0f %%" % (EPS_ROUGH * 1e6, TOL_D * 100),
        software=dict(python=sys.version.split()[0], numpy=np.__version__, cantera=ctv),
        date="2026-11-14",
        inputs=[dict(file="routes/gtu/frame3/frame3_point_results.json", field="cycle.frame3.Wf (расход топлива), stage*.D2_m"),
                dict(file="routes/gtu/funnel/parts_model.py", field="REF_BOM.comb.n_inj = 14")],
        assumptions=dict(DP_INJ_REL=DP_INJ_REL, EPS_ROUGH_m=EPS_ROUGH, L_PIGTAIL_m=L_PIGTAIL, TOL_D=TOL_D, D_ring_m=D_ring,
                         d_pigtail_m=d_pig, R_ring_m=R_ring, L_tree_m=L_tree, gas="CH4 при 288 K",
                         note="всё — оценка облика: размеры, перепад форсунки, шероховатость и допуски приняты"),
        gas=dict(p_Pa=p_man, rho_kg_m3=rho, mu_Pa_s=mu),
        Wf_kg_s=Wf, m_inj_kg_s=m_inj, volume_m3=V_ring,
        ring=dict(spread_pct=float((qA.max() - qA.min()) / qA.mean() * 100), dp_ring_kPa=float((pA.max() - pA.min()) / 1e3),
                  Re_ring=float(4 * Wf / (math.pi * D_ring * mu)), parts=1 + N_INJ, joints=2 * N_INJ + 1,
                  mc_spread_pct_mean=float(np.mean(devA) * 100), mc_spread_pct_p95=float(np.percentile(devA, 95) * 100)),
        tree=dict(D_m=D_opt.tolist(), murray_k_01=k01, murray_k_12=k12, spread_pct=float((qB.max() - qB.min()) / qB.mean() * 100),
                  dp_tree_kPa=float(tree_dp(D_opt, L_tree, Wf, rho, mu) / 1e3), parts=0, joints=1,
                  mc_spread_pct_mean=float(np.mean(devB) * 100), mc_spread_pct_p95=float(np.percentile(devB, 95) * 100),
                  Re_trunk=float(4 * Wf / (math.pi * D_opt[0] * mu))),
        note="показатель Мюррея k из оптимума (ламинарный закон — 3); турбулентное течение даёт меньший k — оптимум «толще "
             "к листьям», чем по кубическому закону")
    with open(RESULTS, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1, default=float)
    print(json.dumps({k: res[k] for k in ("ring", "tree")}, default=float))


if __name__ == "__main__":
    main()
