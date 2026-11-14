"""Двигатель-заменитель в существующем пакете: мощность, частота выходного вала, расход через вход/выхлоп пакета.
Модель цикла - cycle.py (откалибрована по SGT-400 и Taurus 60). КПД узлов и пределы - принятые (оценка)."""
import json, math
from gtu_cycle import cycle
import gtu_explore as E
import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))

ARCH = {   # компрессор: (политропный КПД, макс. степень сжатия, колёс компрессора)
 "1Ц":   dict(name="одно центробежное колесо", eta_pc=0.855, PRmax=6.75, nc=1),
 "2Ц":   dict(name="два центробежных колеса", eta_pc=0.860, PRmax=16.0, nc=2),
 "3О+Ц": dict(name="три осевых блиска и центробежное колесо", eta_pc=0.875, PRmax=18.0, nc=4),
}
COOL = {   # класс охлаждения: (T41 макс, охлаждение СА, охлаждение ротора, ступеней ТГГ при малой/большой степени сжатия)
 "U": dict(name="без охлаждения", T41=1000, xn=0.0, xr=0.0),
 "N": dict(name="охлаждаемый сопловой аппарат, рабочая лопатка без охлаждения", T41=1050, xn=0.03, xr=0.0),
 "K": dict(name="охлаждаемые сопловой аппарат и рабочая лопатка 1-й ступени", T41=1200, xn=0.055, xr=0.055),
}
SLOT = {
 "T60":  dict(name="Taurus 60", P=5.74, G=21.6, EGT=510, eta=32.0, n_pt=14300, wheels=15),
 "S400": dict(name="SGT-400 13 МВт", P=13.4, G=39.4, Gcap=44.0, EGT=555, eta=36.2, n_pt=9500, wheels=15),
}
res = []
for ak, a in ARCH.items():
    for ck, c in COOL.items():
        best = None
        PRs = [x / 4 for x in range(16, int(a["PRmax"] * 4) + 1)]
        for PR in PRs:
            try:
                r = cycle(PR, c["T41"], eta_pc=a["eta_pc"], eta_pt_gg=0.865, eta_pt_pt=0.885,
                          x_ngv=c["xn"], x_rot=c["xr"], x_leak=0.015, W1=1.0)
            except Exception:
                continue
            if best is None or r["sp_work_kJkg"] > best["sp_work_kJkg"]:
                best = dict(r)
        ng = 1 if best["GG_ER"] < 3.2 else 2
        npt = 1 if best["PT_ER"] < 3.0 else 2
        wheels = a["nc"] + ng + npt
        row = dict(arch=ak, arch_name=a["name"], cool=ck, cool_name=c["name"], PR=best["PR"], T41=c["T41"],
                   w=round(best["sp_work_kJkg"]), eta=round(best["eff"] * 100, 1), EGT=round(best["EGT_C"]),
                   wheels=wheels, ng=ng, npt=npt, T3=round(best["T3_C"]))
        for sk, s in SLOT.items():
            G = s["P"] * 1e3 / best["sp_work_kJkg"]
            cap = s.get("Gcap", s["G"])
            row[sk] = dict(G=round(G, 1), dG=round((G / cap - 1) * 100), P_at_cap=round(best["sp_work_kJkg"] * cap / 1e3, 2),
                           P_pct=round(best["sp_work_kJkg"] * cap / 1e3 / s["P"] * 100), fuel_pct=round((s["eta"] / (best["eff"] * 100) - 1) * 100),
                           dEGT=round(best["EGT_C"] - s["EGT"]))
            # размер колеса компрессора первой ступени и частота газогенератора
            PRs1 = best["PR"] if a["nc"] == 1 else (best["PR"] ** 0.5 if ak == "2Ц" else 1.0)
            if ak != "3О+Ц":
                st = E.comp_stage(E.T0, E.P0, PRs1, G, "ti")
                row[sk].update(D2=round(st["D2"], 2), rpm=round(st["rpm"]), D2dbl=round(st["D2"] / math.sqrt(2), 2))
        res.append(row)
json.dump(dict(**{"method": "0D-цикл Брайтона: свойства газа Cantera (GRI-Mech 3.0), метан, ISO; КПД узлов и ограничения - принятые, см. README", "software": "cantera " + __import__("cantera").__version__ + ", python3", "date": __import__("datetime").date.today().isoformat()}, arch=ARCH, cool=COOL, slot=SLOT, res=res), open("gtu_dropin_results.json", "w"), ensure_ascii=False, indent=1)
for r in res:
    t, s = r["T60"], r["S400"]
    print(f"{r['arch']:5s} {r['cool']} колёс={r['wheels']} π={r['PR']:5.2f} T41={r['T41']} w={r['w']} η={r['eta']} EGT={r['EGT']} | "
          f"T60: G={t['G']} ({t['dG']:+d}%) P@21.6={t['P_at_cap']} ({t['P_pct']}%) топл {t['fuel_pct']:+d}% ΔEGT {t['dEGT']:+d} D2={t.get('D2')}/{t.get('D2dbl')} n={t.get('rpm')} | "
          f"S400: G={s['G']} ({s['dG']:+d}% к 44) P@44={s['P_at_cap']} ({s['P_pct']}%) топл {s['fuel_pct']:+d}% ΔEGT {s['dEGT']:+d} D2={s.get('D2')}/{s.get('D2dbl')} n={s.get('rpm')}")
