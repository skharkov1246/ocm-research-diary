#!/usr/bin/env python3
"""Встроенный коллектор топливного газа KV-3: дерево каналов в стенке оболочки против кольцевого коллектора с трубками.

Физика: закон Мюррея (Murray 1926) и конструктальный подход (Bejan 1997) — ветвящаяся сеть, у которой при заданном
объёме каналов наименьшая работа прокачки. Для ламинарного течения D0³ = ΣDi³; здесь течение турбулентное (Re ~10⁶),
поэтому показатель ищется численно: минимум перепада давления дерева при заданном объёме с трением по Haaland.
Сравнение при одинаковом объёме каналов: (A) кольцевой коллектор, подвод в одной точке, 14 трубок к форсункам;
(B) симметричное дерево 1 → 2 → 7 (равные длины путей) в стенке оболочки. Неравномерность расхода по форсункам —
решением сети (давления узлов, fsolve; Δp = R·m|m|, трение пересчитывается по Re), разброс от допусков печати — методом Монте-Карло (±σ диаметра).
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


def R_pipe(m, D, L, rho, mu):
    """квадратичное сопротивление Δp = R·m|m| при трении для расхода m (Haaland)"""
    A = math.pi * D * D / 4
    Re = max(4 * abs(m) / (math.pi * D * mu), 1.0)
    return f_haaland(Re, D) * L / D / (2 * rho * A * A)


def solve_net(edges, n_nodes, src, m_tot, sinks, K_inj, rho, mu):
    """сеть: edges [(a, b, D, L)], узлы 0..n−1; src — узел подвода m_tot; sinks [(узел, D_трубки, L_трубки)] —
    трубка + сопло к давлению p3 (принято 0). Неизвестные — давления узлов (над p3). Две итерации по трению."""
    m_e = np.full(len(edges), m_tot / max(len(sinks), 1))
    m_s = np.full(len(sinks), m_tot / len(sinks))
    x = np.full(n_nodes, K_inj * (m_tot / len(sinks)) ** 2 * 1.05)
    for _ in range(3):
        Re_ = [R_pipe(m, D, L, rho, mu) for m, (a, b, D, L) in zip(m_e, edges)]
        Rs = [R_pipe(m, d, l, rho, mu) + K_inj for m, (nd, d, l) in zip(m_s, sinks)]

        def flows(p):
            fe = np.array([math.copysign(math.sqrt(abs(p[a] - p[b]) / R), p[a] - p[b]) for (a, b, D, L), R in zip(edges, Re_)])
            fs = np.array([math.sqrt(max(p[nd], 0.0) / R) for (nd, d, l), R in zip(sinks, Rs)])
            return fe, fs

        def resid(p):
            fe, fs = flows(p)
            r = np.zeros(n_nodes)
            for (a, b, D, L), q in zip(edges, fe):
                r[a] -= q; r[b] += q
            for (nd, d, l), q in zip(sinks, fs):
                r[nd] -= q
            r[src] += m_tot
            return r / m_tot
        x = fsolve(resid, x, xtol=1e-13)
        m_e, m_s = flows(x)
    return m_s, x


def ring_network(D_ring, R_ring, d_pig, m_tot, rho, mu, K_inj, dvec=None):
    """кольцо из N_INJ участков, подвод в узле 0; форсунки — в узлах через трубку длины L_PIGTAIL"""
    N = N_INJ
    Ls = 2 * math.pi * R_ring / N
    dvec = np.ones(N) if dvec is None else dvec
    edges = [(j, (j + 1) % N, D_ring, Ls) for j in range(N)]
    sinks = [(j, d_pig * dvec[j], L_PIGTAIL) for j in range(N)]
    return solve_net(edges, N, 0, m_tot, sinks, K_inj, rho, mu)


def tree_flows(D, L, m_tot, rho, mu, K_inj, dvec=None):
    """дерево: узел 0 (вход) —ствол→ узел 1 —2 ветви→ узлы 2, 3 —по 7 листьев (лист = канал + сопло)"""
    dvec = np.ones(17) if dvec is None else dvec
    edges = [(0, 1, D[0] * dvec[0], L[0]), (1, 2, D[1] * dvec[1], L[1]), (1, 3, D[1] * dvec[2], L[1])]
    sinks = [(2 if k < 7 else 3, D[2] * dvec[3 + k], L[2]) for k in range(14)]
    return solve_net(edges, 4, 0, m_tot, sinks, K_inj, rho, mu)


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
    qA, pA = ring_network(D_ring, R_ring, d_pig, Wf, rho, mu, K_inj)
    # (B) дерево той же ёмкости: длины уровней — ствол 0,3 м, ветви — четверть окружности, листья — до форсунок
    L_tree = [0.30, math.pi * R_ring / 2, 0.45 * R_ring]
    D_opt, k01, k12, r = murray_tree(Wf, rho, mu, L_tree, V_ring)
    qB, PB = tree_flows(D_opt, L_tree, Wf, rho, mu, K_inj)
    # Монте-Карло допусков
    rng = np.random.default_rng(2026)
    devA, devB = [], []
    for _ in range(200):
        qa, _ = ring_network(D_ring, R_ring, d_pig, Wf, rho, mu, K_inj, dvec=1 + TOL_D * rng.standard_normal(N_INJ))
        qb, _ = tree_flows(D_opt, L_tree, Wf, rho, mu, K_inj, dvec=1 + TOL_D * rng.standard_normal(17))
        devA.append((qa.max() - qa.min()) / qa.mean()); devB.append((qb.max() - qb.min()) / qb.mean())
    # закон Мюррея ламинарный для сравнения (тот же объём): D0³ = 2D1³, D1³ = 7D2³
    res = dict(
        title="Встроенный коллектор топливного газа KV-3: дерево по закону Мюррея против кольца с трубками",
        method="сеть каналов: трение Haaland (ε = %.0f мкм), сопла форсунок Δp = K·m², давления узлов (fsolve, 3 итерации трения); "
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
