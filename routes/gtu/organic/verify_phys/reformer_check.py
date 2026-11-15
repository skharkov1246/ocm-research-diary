"""Независимая проверка риформинга cycles.py:
1) равновесие SMR+WGS решателем на эмпирических константах (K1, K2 — не из Cantera) против Cantera (cycles.equil);
2) тепловой баланс риформера в вариантах с паром (энергетическая невязка узла, тепло газа = тепло сырья);
3) активность углерода по всей составной кривой (равновесный газ при каждой T), а не только на выходе;
4) кинетика Xu-Froment при фактической T риформера варианта (а не 650 °C)."""
import json, math, sys
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, sys.argv[1])
import cycles as C

def eq_emp(sc, T, P_bar):
    K1 = math.exp(-26830.0 / T + 30.114); K2 = math.exp(4400.0 / T - 4.036)
    # моли на 1 моль CH4: x — прореагировал по SMR, y — по WGS
    def comp(x, y):
        n = dict(CH4=1 - x, H2O=sc - x - y, CO=x - y, H2=3 * x + y, CO2=y)
        return n, sum(n.values())
    def y_of(x):
        f = lambda y: (lambda n, N: n["CO2"] * n["H2"] - K2 * n["CO"] * n["H2O"])(*comp(x, y))
        return brentq(f, 0.0, min(x, sc - x) * (1 - 1e-12), xtol=1e-15)
    def g(x):
        y = y_of(x); n, N = comp(x, y)
        return n["CO"] * n["H2"] ** 3 * (P_bar / N) ** 2 - K1 * n["CH4"] * n["H2O"]
    x = brentq(g, 1e-6, min(1.0, sc) * (1 - 1e-9), xtol=1e-12)
    return x

rows = []
for P_bar in (5.0, 17.0):
    for sc in (1.5, 2.0, 3.0):
        Yf, h_in, yC = C.feed_steam(sc)
        for T_C in (450, 550, 650, 750, 850):
            T = T_C + 273.15
            Ye, _ = C.equil(Yf, T, P_bar * 1e5)
            Xc = 1 - (Ye[C.I["CH4"]] / C.MW[C.I["CH4"]]) / (yC / C.MW[C.I["CH4"]])
            Xe = eq_emp(sc, T, P_bar)
            rows.append(dict(P_bar=P_bar, sc=sc, T_C=T_C, X_cantera=Xc, X_empirical=Xe, diff_pp=100 * (Xc - Xe)))
print("max |ΔX| п.п.:", max(abs(r["diff_pp"]) for r in rows))
for r in rows:
    if r["sc"] == 2.0: print(r)

# 2–4: варианты с паровым риформингом из результата
res = json.load(open(sys.argv[2]))
out = dict(equilibrium=rows, variants=[])
keys = [r["key"] for r in res["results"] if r["topo"]["tcr"] == "steam" and r.get("opt", {}).get("ok")]
pick = [k for k in keys if k in ("R00-ic-rh-ibc-none-steam", "R90-ic-RH-ibc-none-steam", "R00-ic-RH-ibc-none-steam", "R92-ic-RH-ibc-none-steam")]
for k in pick:
    o = next(r for r in res["results"] if r["key"] == k)
    for tag in ("opt", "as_is"):
        if not o.get(tag) or not o[tag].get("ok"): continue
        rr = C.evaluate(o["topo"], o[tag]["PR"], o[tag]["prm"], detail=True)
        ea = C.exergy_account(rr)
        ref = [x for x in ea["rows"] if x["node"].startswith("Риформер")][0]
        st = rr["_streams"]; fuel = st["fuel"]; info = fuel["info"]
        sc = info["sc"]; P = st["p_fuel"]; Tref = info["T_ref_C"] + 273.15
        pts = C.cold_curve_steam(sc, P, Tref)
        Yf, h_in, yC = C.feed_steam(sc, P)
        Ts = C.T_sat(P)
        aC_curve = []
        for T, q in pts:
            if T <= Ts + 0.02: continue
            Teq = max(T - C.A["approach_eq"], Ts)
            Ye, _ = C.equil(Yf, Teq, P)
            aC_curve.append((T - 273.15, max(C.carbon_activity(Ye, Teq, P))))
        # тепловой баланс: тепло газа (энтальпия) = тепло сырья
        W1 = rr["W_air"]
        hg_in = C.hs(st["st5"])[0] * st["st5"].m; hg_out = C.hs(st["g_in"])[0] * st["g_in"].m
        Q_gas = (hg_in - hg_out) * W1 / 1e3
        Q_cold = rr["mf"] * fuel["q"] * W1 / 1e3
        # кинетика при фактической T
        Teq = max(Tref - C.A["approach_eq"], Ts)
        Xf = np.zeros(len(C.SP)); Xf[C.I["CH4"]], Xf[C.I["H2O"]] = 1 / (1 + sc), sc / (1 + sc)
        Yk = Xf * C.MW / np.dot(Xf, C.MW)
        Ye, _ = C.equil(Yk, Teq, P)
        Xeq = 1 - (Ye[C.I["CH4"]] / C.MW[C.I["CH4"]]) / (Yk[C.I["CH4"]] / C.MW[C.I["CH4"]])
        Wc, Xc, _ = C.xf_pfr(Xf, 1.0 + sc, Teq, P, 20.0)
        W95 = float(np.interp(0.95 * Xeq, Xc, Wc)) if Xc[-1] >= 0.95 * Xeq else None
        n_ch4 = rr["W_fuel_CH4"] / (C.MW[C.I["CH4"]] / 1000.0)
        v = dict(key=k, tag=tag, PR=o[tag]["PR"], T_ref_C=info["T_ref_C"], X_CH4=info["X_CH4"], pinch_K=info["pinch_K"],
                 aC_outlet=max(info["aC"]), aC_curve_max=max(a for _, a in aC_curve), aC_curve_argmax_C=max(aC_curve, key=lambda z: z[1])[0],
                 reformer_energy_residual_kW=ref["energy_residual_kW"], Q_gas_kW=Q_gas, Q_cold_kW=Q_cold,
                 Teq_C=Teq - 273.15, P_bar=P / 1e5, X_eq_at_Teq=Xeq, W95_kg_per_molCH4_s=W95, CH4_mol_s=n_ch4, sc=sc,
                 X_over_eq_at_1kg=float(np.interp(1.0, Wc, Xc) / Xeq), X_over_eq_at_20kg=float(Xc[-1] / Xeq),
                 catalyst_kg_at_Tref=(W95 * n_ch4 if W95 else None), eff_el=o[tag]["eff_el"])
        out["variants"].append(v)
        print(json.dumps(v, ensure_ascii=False, default=float))
json.dump(out, open(sys.argv[3], "w"), ensure_ascii=False, indent=1, default=float)
