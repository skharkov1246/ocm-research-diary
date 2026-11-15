"""Воронка камеры сгорания KV, уровень 3: низкопорядковая модель термоакустики (сеть акустических элементов).

Метод - одномерная сеть элементов с плоскими волнами (продольные моды), как в работах Dowling & Stow (J. Propulsion
and Power 19, 2003, «Acoustic analysis of gas turbine combustors») и в открытом коде OSCILOS (Imperial College,
Li & Morgans) - оба по памяти и выдержкам поиска, тексты из среды не открывались. Элементы:
  * труба длины L, площади A, с температурой T и числом Маха M среднего течения: p' = A+·e^(-i k+ x) + A-·e^(i k- x),
    k± = ω / (c·(1 ± M)); ρc·u' = A+·e^(-i k+ x) − A-·e^(i k- x);
  * скачок площади (компактный, малое M): p' непрерывно, A·u' непрерывно;
  * компактное пламя на стыке двух труб (малое M, Ранкин - Гюгонио первого порядка): p' непрерывно,
    A2·u2' = A1·u1'·(1 + (T2/T1 − 1)·F(ω)), F(ω) - передаточная функция пламени Q'/Q̄ = F·u1'/ū1;
    модель n - τ: F = n·e^(−iωτ) (на уровне 3 - F(ω) из LES с возбуждением, README);
  * концы: коэффициент отражения R = A-/A+ на входе (по волнам, идущим от конца) и на выходе; закрытый конец R = 1,
    открытый R = −1; запертое сопло турбины - компактное сопло Марбла - Кандела, R = (1 − (γ−1)/2·M)/(1 + (γ−1)/2·M)
    (по памяти, не сверено).
Собственные частоты - нули определителя системы в комплексной плоскости ω (зависимость p' ~ e^(iωt));
рост колебаний при Im ω < 0: инкремент σ = −Im ω, 1/с. Проверка - аналитика (test_comb_acoustic.py).

Запуск: python comb_acoustic.py --selfcheck | --inputs l1_inputs.json --row 0 [--n 1 --tau-ms 2]
"""
import argparse, cmath, json, math, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GAMMA = 1.4
R_AIR = 287.05
DATE = "2026-11-14"


def sound(T, g=GAMMA, R=R_AIR):
    return math.sqrt(g * R * T)


class Network:
    """последовательность труб; между трубами i и i+1 - стык (площадь и, если задано, пламя)"""

    def __init__(self, ducts, R_in, R_out, flame_at=None, n=0.0, tau=0.0, F=None):
        self.d = ducts                    # [dict(L, A, T, M=0, rho=..)]
        self.R_in, self.R_out = R_in, R_out
        self.flame_at, self.n, self.tau, self.F = flame_at, n, tau, F
        for x in self.d:
            x.setdefault("M", 0.0)
            x.setdefault("p", 1.0e5)
            x["c"] = sound(x["T"])
            x["rho"] = x["p"] / (R_AIR * x["T"])

    def ftf(self, w):
        return self.F(w) if self.F else self.n * cmath.exp(-1j * w * self.tau)

    def det(self, w):
        """определитель системы на амплитуды (A+, A-) всех труб"""
        N = len(self.d)
        Mx = np.zeros((2 * N, 2 * N), complex)
        row = 0
        # вход: A-(0) = R_in·A+(0) -> волна, идущая от входа (A+), отражённая (A-)... на входе x = 0: A+ = R_in·A-
        Mx[row, 0], Mx[row, 1] = 1.0, -self.R_in
        row += 1
        for i in range(N - 1):
            a, b = self.d[i], self.d[i + 1]
            kp = w / (a["c"] * (1 + a["M"]))
            km = w / (a["c"] * (1 - a["M"]))
            ep, em = cmath.exp(-1j * kp * a["L"]), cmath.exp(1j * km * a["L"])
            # p' непрерывно
            Mx[row, 2 * i], Mx[row, 2 * i + 1] = ep, em
            Mx[row, 2 * i + 2], Mx[row, 2 * i + 3] = -1.0, -1.0
            row += 1
            # объёмный расход: A_b·u_b = A_a·u_a·(1 + θ·F) на пламени, иначе A_a·u_a
            fac = 1.0
            if self.flame_at == i:
                fac = 1.0 + (b["T"] / a["T"] - 1.0) * self.ftf(w)
            ua = a["A"] / (a["rho"] * a["c"]) * fac
            ub = b["A"] / (b["rho"] * b["c"])
            Mx[row, 2 * i], Mx[row, 2 * i + 1] = ua * ep, -ua * em
            Mx[row, 2 * i + 2], Mx[row, 2 * i + 3] = -ub, ub
            row += 1
        z = self.d[-1]
        kp = w / (z["c"] * (1 + z["M"]))
        km = w / (z["c"] * (1 - z["M"]))
        ep, em = cmath.exp(-1j * kp * z["L"]), cmath.exp(1j * km * z["L"])
        # выход x = L: A-·e^(i k- L) = R_out·A+·e^(-i k+ L)
        Mx[row, 2 * N - 2], Mx[row, 2 * N - 1] = -self.R_out * ep, em
        s = np.linalg.svd(Mx, compute_uv=False)
        return np.linalg.det(Mx), s[-1] / s[0]

    def modes(self, f_max=2000.0, n_scan=4000, f_min=5.0):
        """собственные моды до f_max: минимумы |det| на вещественной оси -> уточнение секущей в комплексной плоскости"""
        fs = np.linspace(f_min, f_max, n_scan)
        vals = np.array([abs(self.det(2 * math.pi * f)[0]) for f in fs])
        cand = [fs[i] for i in range(1, n_scan - 1) if vals[i] < vals[i - 1] and vals[i] <= vals[i + 1]]
        out = []
        for f0 in cand:
            w0, w1 = 2 * math.pi * f0, 2 * math.pi * f0 * (1 + 1e-3) + 1j
            d0, d1 = self.det(w0)[0], self.det(w1)[0]
            ok = False
            for _ in range(80):
                if d1 == d0:
                    break
                w2 = w1 - d1 * (w1 - w0) / (d1 - d0)
                w0, d0, w1 = w1, d1, w2
                d1 = self.det(w1)[0]
                if abs(w1 - w0) < 1e-9 * abs(w1):
                    ok = True
                    break
            if ok and f_min * 2 * math.pi * 0.5 < w1.real < f_max * 2 * math.pi * 1.05:
                if all(abs(w1 - m["omega"]) > 1e-3 * abs(w1) for m in out):
                    out.append(dict(omega=w1, f_Hz=w1.real / (2 * math.pi), growth_1_s=-w1.imag,
                                    cycle_increment=math.exp(-w1.imag * 2 * math.pi / w1.real) if w1.real > 0 else None))
        return sorted(out, key=lambda m: m["f_Hz"])


def marble_candel_R(M, g=GAMMA):
    """отражение запертого компактного сопла по волнам, приходящим сверху по течению (по памяти, не сверено)"""
    a = 0.5 * (g - 1) * M
    return (1 - a) / (1 + a)


def kv_network(e, n=1.0, tau=None, R_in=1.0, transition_L=0.3):
    """сеть камеры варианта уровня 1/0: объём за компрессором (кожух) -> смеситель -> жаровая труба -> переходный
    канал -> запертый сопловой аппарат. Размеры - comb_space.derive; τ по умолчанию - время смешения + половина
    времени зоны пламени (оценка; на уровне 3 - из LES)"""
    c, d, L = e["concept"], e["derived"], e["loads"][0]
    p3, T3 = L["p3_Pa"], L["T3_K"]
    T4 = L["T4_K"]
    T_pz = L["hot_init"]["T_eq"]
    n_can = e["n_can"]
    A_L = d["A_L"]
    A_case = d["A_ref"]
    A_mix = d["A_air_m2"] / n_can
    l_mix = c["l_mix"] * d["d_air_mm"] * 1e-3
    W = L["W_air"] / n_can
    u_mix = W / (L["rho3"] * A_mix)
    rho4 = p3 / (R_AIR * T4)
    M4 = W / (rho4 * A_L) / sound(T4)
    if tau is None:
        tau = 1e-3 * (d["tau_mix_ms"] + 0.5 * c["tau_pz_ms"])
    ducts = [dict(L=0.6 * d["L_mod"], A=A_case, T=T3, p=p3, M=min(0.2, W / (L["rho3"] * A_case) / sound(T3))),
             dict(L=max(l_mix, 1e-3), A=A_mix, T=T3, p=p3, M=min(0.3, u_mix / sound(T3))),
             dict(L=0.4 * d["L_L"], A=A_L, T=T_pz, p=p3, M=min(0.2, M4)),
             dict(L=0.6 * d["L_L"], A=A_L, T=T4, p=p3, M=min(0.2, M4)),
             dict(L=transition_L, A=0.7 * A_L, T=T4, p=p3, M=min(0.3, M4 / 0.7))]
    net = Network(ducts, R_in=R_in, R_out=marble_candel_R(min(0.3, M4 / 0.7)), flame_at=1, n=n, tau=tau)
    return net, dict(tau_s=tau, n=n, ducts=[{k: v for k, v in x.items()} for x in ducts], R_in=R_in,
                     R_out=net.R_out, note="flame_at=1: пламя на стыке смеситель - зона пламени, θ = T_зоны/T3 − 1")


def azimuthal_estimate(e):
    """частоты азимутальных мод кольцевой камеры без расчёта сети: f_m = m·c/(π·D_m) (оценка нижней границы опасной
    полосы; для трубчатых и выносной - нет)"""
    d = e["derived"]
    if not d.get("D_m"):
        return None
    c = sound(e["loads"][0]["T4_K"])
    return [dict(m=m, f_Hz=m * c / (math.pi * d["D_m"])) for m in (1, 2, 3)]


def selfcheck():
    """аналитика: закрытая - открытая и закрытая - закрытая труба, скачок температуры, знак роста трубы Рийке"""
    out = {}
    T, L = 300.0, 1.0
    c = sound(T)
    net = Network([dict(L=L, A=0.01, T=T)], R_in=1.0, R_out=-1.0)
    m = net.modes(f_max=1000)
    out["closed_open"] = dict(computed=[x["f_Hz"] for x in m[:3]], analytic=[(2 * k - 1) * c / (4 * L) for k in (1, 2, 3)])
    net = Network([dict(L=L, A=0.01, T=T)], R_in=1.0, R_out=1.0)
    m = net.modes(f_max=1000)
    out["closed_closed"] = dict(computed=[x["f_Hz"] for x in m[:3]], analytic=[k * c / (2 * L) for k in (1, 2, 3)])
    # две трубы, скачок температуры T1 -> T2, равные площади, закрытая - открытая, без отклика пламени:
    # p1 = cos(k1 x), p2 = sin(k2 (L − x)); на стыке x = b непрерывны p и u: tan(k1 b)·tan(k2 (L − b)) = ρ1c1/(ρ2c2),
    # при одном давлении ρc ~ 1/c, то есть = c2/c1
    T1, T2, b = 300.0, 1200.0, 0.3
    c1, c2 = sound(T1), sound(T2)
    net = Network([dict(L=b, A=0.01, T=T1), dict(L=L - b, A=0.01, T=T2)], R_in=1.0, R_out=-1.0, flame_at=0, n=0.0)
    m = net.modes(f_max=600)
    from scipy.optimize import brentq
    g = lambda f: math.tan(2 * math.pi * f / c1 * b) * math.tan(2 * math.pi * f / c2 * (L - b)) - c2 / c1
    fs = np.linspace(1, 600, 60000)
    roots = []
    for x0, x1 in zip(fs[:-1], fs[1:]):
        try:
            if g(x0) * g(x1) < 0 and abs(g(x0)) < 50 and abs(g(x1)) < 50:
                roots.append(brentq(g, x0, x1))
        except ValueError:
            pass
    out["temperature_jump"] = dict(computed=[x["f_Hz"] for x in m[:3]], analytic=roots[:3])
    # труба Рийке: открытая - открытая, пламя на L/4, n = 1, малая задержка -> первая мода растёт
    net = Network([dict(L=0.25, A=0.01, T=300.0), dict(L=0.75, A=0.01, T=600.0)], R_in=-1.0, R_out=-1.0, flame_at=0,
                  n=1.0, tau=0.2e-3)
    m = net.modes(f_max=400)
    out["rijke_first_mode"] = dict(f_Hz=m[0]["f_Hz"], growth_1_s=m[0]["growth_1_s"], expect="рост > 0 (Рэлей)")
    ok = all(abs(a - b) / b < 2e-3 for k in ("closed_open", "closed_closed", "temperature_jump")
             for a, b in zip(out[k]["computed"], out[k]["analytic"])) and m[0]["growth_1_s"] > 0
    out["all_ok"] = bool(ok)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--selfcheck", action="store_true")
    ap.add_argument("--inputs")
    ap.add_argument("--row", type=int, default=0)
    ap.add_argument("--n", type=float, default=1.0)
    ap.add_argument("--tau-ms", type=float, default=None)
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    if a.selfcheck:
        r = selfcheck()
        print(json.dumps(r, ensure_ascii=False, indent=1))
        return 0 if r["all_ok"] else 1
    e = json.load(open(a.inputs, encoding="utf-8"))["entries"][a.row]
    net, info = kv_network(e, a.n, None if a.tau_ms is None else 1e-3 * a.tau_ms)
    m = net.modes(f_max=1500)
    res = dict(idx=e["idx"], scheme=e["concept"]["scheme"], chamber=e["concept"]["chamber"], network=info,
               modes=[dict(f_Hz=x["f_Hz"], growth_1_s=x["growth_1_s"]) for x in m], azimuthal=azimuthal_estimate(e))
    print(json.dumps(res, ensure_ascii=False, indent=1, default=float))
    if a.out:
        json.dump(res, open(a.out, "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=float)
    return 0


if __name__ == "__main__":
    sys.exit(main())
