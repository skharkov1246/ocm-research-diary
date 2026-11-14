"""Цифры для применения на добыче: давление топливного газа, дожимной компрессор, режим с обводом рекуператора."""
import math, json
import gtu_explore as E
out = {}
# 1. Давление топливного газа: p_топл = 1,4 * p_за компрессором (принято), подача ПНГ 0,3 и 0,6 МПа абс.
R_CH4, T1, eta_iso = 518.3, 293.15, 0.62           # изотермический КПД поршневого дожимного компрессора (принято)
cases = [("Турбокластер 6 МВт", 4.25, 6.0, 0.337), ("Taurus 60, привод", 12.2, 5.74, 0.3196), ("SGT-400 13 МВт, привод", 16.8, 13.4, 0.362)]
fuel = []
for name, PR, P, eta in cases:
    p_fuel = 1.4 * PR * 0.101325                     # МПа абс
    mf = P / eta / 50.03                             # кг/с метана
    row = dict(name=name, PR=PR, p_fuel=round(p_fuel, 2), mf_kgh=round(mf * 3600))
    for pin in (0.3, 0.6):
        r = max(1.0, p_fuel / pin)
        W = mf * R_CH4 * T1 * math.log(r) / eta_iso / 1e3   # кВт
        row[f"W_{pin}"] = round(W); row[f"share_{pin}"] = round(W / (P * 1e3) * 100, 1)
    fuel.append(row); print(row)
out["fuel"] = fuel
# 2. Режимы одной машины: рекуператор включён / обведён (та же степень сжатия и температура)
base = dict(fam="C", name="tc", N=4, spools=1, cstages=1, comp="al", rotor="two", ggt="radial", shaft="free", pt="ax1", use="gen", drive="direct")
G = 39.8
modes = {}
import cantera as ct
import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))
g = ct.Solution("gri30.yaml")
for eps in (0.85, 0.0):
    c = dict(base, eps=eps); r = E.evaluate(c, 4.25, 950.0)
    P = r["w_net"] * G / 1e6
    Q_fuel = r["f"] * G * E.LHV / 1e6
    T_out = (r["T5"] - (r["T3"] - r["T2"]) * 0.93) if eps > 0 else r["T5"]
    b = r["b"]
    mgas = G * (r["mg"] + r["x_leak"])
    Q_heat = mgas * (E.h(T_out, b) - E.h(393.15, b)) / 1e6      # охлаждение уходящих газов до 120 C
    modes[eps] = dict(P=round(P, 2), eff=round(r["eff"] * 100, 1), eff_el=round(r["eff"] * 0.972 * 100, 1), T_out=round(T_out - 273.15),
                      Q_heat=round(Q_heat, 1), util=round((P * 0.972 + Q_heat) / Q_fuel * 100), fuel_MW=round(Q_fuel, 1),
                      apg_m3h=round(Q_fuel / 40.0 * 3600))     # ПНГ, низшая теплота 40 МДж/м3 (принято)
    print(eps, modes[eps])
out["modes"] = {str(k): v for k, v in modes.items()}
# 3. Объём ПНГ на одну установку в год (8 000 ч)
out["apg_year_mln"] = round(modes[0.85]["apg_m3h"] * 8000 / 1e6, 1)
out["apg_year_mln_bypass"] = round(modes[0.0]["apg_m3h"] * 8000 / 1e6, 1)
print("ПНГ в год, млн м3:", out["apg_year_mln"], "/ при обводе:", out["apg_year_mln_bypass"])
out.update({"method": "0D-цикл Брайтона: свойства газа Cantera (GRI-Mech 3.0), метан, ISO; КПД узлов и ограничения - принятые, см. README", "software": "cantera " + __import__("cantera").__version__ + ", python3", "date": __import__("datetime").date.today().isoformat()}); json.dump(out, open("gtu_upstream_results.json", "w"), ensure_ascii=False, indent=1)
