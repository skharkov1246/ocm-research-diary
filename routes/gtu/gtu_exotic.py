"""Нетрадиционные схемы: проверка тем же аппаратом свойств (explore.py). Результат - в exotic.json."""
import math, json
import gtu_explore as E
import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))
out = []
# 1. Дисковая (безлопаточная) турбина Теслы вместо лопаточной. КПД турбины 0.40 - верх диапазона опубликованных измерений (принято).
for eta_t in (0.40, 0.55):
    T1, PR, TIT = E.T0, 4.0, 1173.15
    c = E.comp_stage(T1, E.P0, PR, 8.0, "ti")
    f = E.far_for(c["T2"], TIT) / 0.995; b = E.bfrac(f)
    ER = PR * 0.955 / 1.03
    T2s = E.T_s(E.s0(TIT, b) - E.Rg(b) * math.log(ER), b); dhs = E.h(TIT, b) - E.h(T2s, b)
    w = (eta_t * dhs * (1 + f) - c["w"] / 0.99) * 0.98
    out.append(dict(name=f"Одновальный с дисковой турбиной Теслы (КПД турбины {eta_t:.2f})", PR=PR, TIT=900,
                    eff=round(w / (f * E.LHV) * 100, 1), w=round(w / 1e3)))
# 2. Обращённый цикл Брайтона: сгорание при атмосферном давлении, расширение в разрежение, охлаждение, сжатие до атмосферы.
for eps in (0.85, 0.90):
    best = None
    for PRc in (2.5, 3.0, 3.5, 4.0, 5.0):
        TIT = 1223.15
        p_in = E.P0 * (1 - 0.025) * (1 - 0.03)                    # рекуператор (воздух) и камера
        p_low = E.P0 * 1.02 / PRc                                 # вход компрессора
        p_t_out = p_low / ((1 - 0.035) * (1 - 0.02))              # выход турбины: рекуператор (газ) и охладитель
        ER = p_in / p_t_out
        T3 = E.T0
        for _ in range(8):
            f = E.far_for(T3, TIT) / 0.995; b = E.bfrac(f)
            T5 = E.T_s(E.s0(TIT, b) - 0.86 * E.Rg(b) * math.log(ER), b)
            wt = (E.h(TIT, b) - E.h(T5, b)) * (1 + f)
            T3 = E.T0 + eps * (T5 - E.T0)
        # компрессор сжимает охлаждённые до 35 C продукты; свойства продуктов
        Tc1 = 308.15
        Tc2 = E.T_s(E.s0(Tc1, b) + E.Rg(b) * math.log(PRc) / 0.85, b)
        wc = (E.h(Tc2, b) - E.h(Tc1, b)) * (1 + f)
        w = (wt - wc / 0.99) * 0.98
        if w <= 0: continue
        eff = w / (f * E.LHV)
        m = 6e6 / w
        rho = p_low / (E.Rg(b) * Tc1); U2 = math.sqrt((E.h(Tc2, b) - E.h(Tc1, b)) / 0.70)
        D2 = math.sqrt(m / (rho * U2 * 0.085))
        r = dict(name=f"Обращённый цикл Брайтона, рекуператор {eps:.2f}", PR=PRc, TIT=950, eff=round(eff * 100, 1), w=round(w / 1e3),
                 m=round(m, 1), D2=round(D2, 2), rpm=round(60 * U2 / (math.pi * D2)))
        if best is None or r["eff"] > best["eff"]: best = r
    out.append(best)
json.dump(dict(**{"method": "0D-цикл Брайтона: свойства газа Cantera (GRI-Mech 3.0), метан, ISO; КПД узлов и ограничения - принятые, см. README", "software": "cantera " + __import__("cantera").__version__ + ", python3", "date": __import__("datetime").date.today().isoformat()}, res=out), open("gtu_exotic_results.json", "w"), ensure_ascii=False, indent=1)
for r in out: print(r)
