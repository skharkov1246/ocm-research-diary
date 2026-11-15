"""Влияние соглашения о воде: химическая эксергия выхлопа с водой-паром (как exergy.py) против состояния
с конденсацией избыточной воды при T0 (ограниченное мёртвое состояние: газ насыщен, избыток — жидкость).
Плюс независимая проверка эксергии CH4 по табличным ΔG_f (NIST-JANAF, 298,15 К, по памяти) с поправкой ΔH, ΔS на T0."""
import sys, math, json, os
import numpy as np
sys.path.insert(0, sys.argv[1])
cwd = os.getcwd()
import exergy as X
import cantera as ct
os.chdir(cwd)
f3 = X.load_frame3(); env = X.Env()
kw0 = X.base_kw(f3)
base, ch0 = X.run_simple(kw0, env, f3)
s5 = ch0["st"]["s5"]; W45 = ch0["m"]["W45"]
g = ct.Solution("gri30.yaml"); g.TPY = X.T0, X.P0, s5["Y"]
x = g.X.copy(); iw = g.species_index("H2O")
xw = x[iw]; xs = env.p_sat / X.P0
R = ct.gas_constant / 1e3  # кДж/(кмоль·К)
R_mol = R / 1e3            # кДж/(моль·К)
# на 1 кмоль смеси: конденсирует n так, что (xw - n)/(1 - n) = xs
n = (xw - xs) / (1 - xs)
# ΔG при изотермической конденсации от xw до xs (интеграл RT ln(xs/x(n)) dn), кДж/кмоль смеси
ns = np.linspace(0, n, 20001); xx = (xw - ns) / (1 - ns)
dG = np.trapezoid(R * X.T0 * np.log(xs / xx), ns)
Mmix = g.mean_molecular_weight
E_fuel = base["fuel"]["Ex_kW"]
corr_kW = dG / Mmix * W45   # кДж/кг × кг/с = кВт (отрицательно)
RES = dict(h2o=dict(x_H2O_exhaust=xw, x_sat_T0=xs, condensed_kmol_per_kmol=n, dEx_kJ_per_kg=dG / Mmix,
                      exhaust_correction_kW=corr_kW, pct_fuel_ex=100 * corr_kW / E_fuel,
                      exhaust_pct_reported=base["nodes"]["exhaust"]["pct_fuel_ex"],
                      exhaust_chem_pct=100 * base["nodes"]["exhaust"]["chemical_W"] / E_fuel))
print(json.dumps(RES, ensure_ascii=False, indent=1))
# CH4 по табличным ΔG_f°(298) и ΔH_f°(298), кДж/моль (NIST-JANAF; по памяти): CH4 -50.72/-74.87, CO2 -394.36/-393.51, H2O(г) -228.58/-241.83
dG298 = -394.36 + 2 * -228.58 - (-50.72); dH298 = -393.51 + 2 * -241.83 - (-74.87)
dS = (dH298 - dG298) / 298.15
dGT0 = dH298 - X.T0 * dS     # без ΔCp (поправка ~0,1 кДж/моль при ΔT 10 К)
xe = env.X
ex = -dGT0 + R_mol * X.T0 * math.log(xe["O2"] ** 2 / (xe["CO2"] * xe["H2O"] ** 2))
RES["ch4"] = dict(ex_ch_tables_kJ_mol=ex, ex_ch_exergy_py_kJ_mol=X.fuel_exergy_check(env)["ex_ch_kJ_mol_mu"],
                  note="ΔG_f°, ΔH_f°(298) NIST-JANAF по памяти; без ΔCp")
print(RES["ch4"])
if len(sys.argv) > 2:
    json.dump(RES, open(sys.argv[2], "w"), ensure_ascii=False, indent=1, default=float)
