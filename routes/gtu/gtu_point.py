import json
from scipy.optimize import brentq
from gtu_cycle import cycle
import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))
s = json.load(open("gtu_study_results.json")); sB = s["SGT400"]
base = dict(x_ngv=0.055, x_rot=0.055, x_leak=0.015)
def cool(T41):
    k = 1 + (T41-1200)/50*0.12
    return dict(x_ngv=0.055*k, x_rot=0.055*k, x_leak=0.015)
E = {"B": (sB["eta_pc"]+0.005, sB["eta_pt_gg"]+0.005, sB["eta_pt_pt"]+0.005),
     "A": (sB["eta_pc"]-0.010, sB["eta_pt_gg"]-0.005, sB["eta_pt_pt"]-0.005)}
W = {"A": 20.5, "B": 41.0}
out = {}
def run(size, PR, T41, tag):
    ec, eg, ep = E[size]
    r = cycle(PR, T41, eta_pc=ec, eta_pt_gg=eg, eta_pt_pt=ep, W1=W[size], **cool(T41))
    r.update(size=size, tag=tag, eta_pc=ec, eta_gg=eg, eta_pt=ep, W1=W[size], eff_el=r["eff"]*0.985*0.972, P_el=r["P_MW"]*0.985*0.972)
    out[f"{size}:{tag}"] = r
    print(f"{size} {tag:10s} PR={PR} T41={T41:.0f} T4={r['T4_C']:.0f} P={r['P_MW']:.2f} МВт КПД={r['eff']*100:.1f} эл={r['eff_el']*100:.1f}% Pэл={r['P_el']:.2f} EGT={r['EGT_C']:.0f} Gвых={r['W_exh']:.1f} T3={r['T3_C']:.0f} T45={r['T45_C']:.0f} ГГ π_т={r['GG_ER']:.2f} СТ π_т={r['PT_ER']:.2f} Lк={r['wc_kJkg']:.0f} топл={r['Wf']*3600:.0f} кг/ч")
    return r
for size in "AB":
    for PR in (18,):
        run(size, PR, 1200, "номинал")
        run(size, PR, 1250, "форс")
        run(size, PR, 1150, "1150")
# режим мощности аналога: какая T41 нужна
for size, P in (("A", 5.74), ("B", 13.4)):
    ec, eg, ep = E[size]
    f = lambda t: cycle(18, t, eta_pc=ec, eta_pt_gg=eg, eta_pt_pt=ep, W1=W[size], **cool(t))["P_MW"] - P
    t = brentq(f, 950, 1300, xtol=0.5)
    run(size, 18, t, "как аналог")
out.update({"method": "0D-цикл Брайтона: свойства газа Cantera (GRI-Mech 3.0), метан, ISO; КПД узлов и ограничения - принятые, см. README", "software": "cantera " + __import__("cantera").__version__ + ", python3", "date": __import__("datetime").date.today().isoformat()}); json.dump(out, open("gtu_point_results.json", "w"), ensure_ascii=False, indent=1)
