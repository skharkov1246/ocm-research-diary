import json
from scipy.optimize import brentq
from gtu_cycle import cycle
import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))

def fitT(PR, W1, P, **kw):
    t = brentq(lambda t: cycle(PR, t, W1=W1, **kw)["P_MW"] - P, 850, 1450, xtol=0.3)
    return cycle(PR, t, W1=W1, **kw)

def calib(name, PR, W1, P, eff, base, cool):
    def res(d):
        kw = dict(eta_pc=base[0]+d, eta_pt_gg=base[1]+d, eta_pt_pt=base[2]+d, **cool)
        return fitT(PR, W1, P, **kw)["eff"] - eff
    d = brentq(res, -0.04, 0.04, xtol=2e-4)
    kw = dict(eta_pc=base[0]+d, eta_pt_gg=base[1]+d, eta_pt_pt=base[2]+d, **cool)
    r = fitT(PR, W1, P, **kw); r.update(name=name, d=d, **kw)
    return r

out = {}
coolB = dict(x_ngv=0.055, x_rot=0.055, x_leak=0.015)
coolA = dict(x_ngv=0.045, x_rot=0.045, x_leak=0.015)
out["SGT400"] = calib("SGT-400 13 МВт, привод", 16.8, 39.4/1.006, 13.4, 0.362, (0.895, 0.875, 0.895), coolB)
out["T60"] = calib("Taurus 60, привод", 12.2, 21.63/1.004, 5.74, 0.3196, (0.875, 0.865, 0.885), coolA)
for k, r in out.items():
    print(f"{r['name']}: eta_pc={r['eta_pc']:.3f} eta_gg={r['eta_pt_gg']:.3f} eta_pt={r['eta_pt_pt']:.3f} T41={r['T41_C']:.0f} T4={r['T4_C']:.0f} P={r['P_MW']:.2f} eff={r['eff']*100:.2f} EGT={r['EGT_C']:.0f} Wexh={r['W_exh']:.1f} T3={r['T3_C']:.0f}")

# --- сетка для нового двигателя ---
sB = out["SGT400"]; sA = out["T60"]
grid = []
variants = {
 # (типоразмер, компоновка): (eta_pc, eta_gg, eta_pt, cooling)
 ("B", "осевой 11 блисков"): (sB["eta_pc"]+0.005, sB["eta_pt_gg"]+0.005, sB["eta_pt_pt"]+0.005, coolB),
 ("B", "осецентробежный 5+1"): (sB["eta_pc"]-0.010, sB["eta_pt_gg"]+0.005, sB["eta_pt_pt"]+0.005, coolB),
 ("A", "осевой 11 блисков"): (sB["eta_pc"]-0.010, sB["eta_pt_gg"]-0.005, sB["eta_pt_pt"]-0.005, coolB),
 ("A", "осецентробежный 5+1"): (sB["eta_pc"]-0.015, sB["eta_pt_gg"]-0.005, sB["eta_pt_pt"]-0.005, coolB),
}
for (size, arch), (ec, eg, ep, cool) in variants.items():
    for PR in (14, 16, 18):
        for T41 in (1150, 1200, 1250):
            c = dict(cool)
            k = 1 + (T41-1200)/50*0.12        # охлаждение растёт на 12 % отн. на каждые 50 C
            c["x_ngv"] *= k; c["x_rot"] *= k
            r = cycle(PR, T41, eta_pc=ec, eta_pt_gg=eg, eta_pt_pt=ep, **c)
            grid.append(dict(size=size, arch=arch, PR=PR, T41=T41, eff=r["eff"], sp=r["sp_work_kJkg"],
                             EGT=r["EGT_C"], T3=r["T3_C"], T4=r["T4_C"], eta_pc=ec))
out["grid"] = grid
out.update({"method": "0D-цикл Брайтона: свойства газа Cantera (GRI-Mech 3.0), метан, ISO; КПД узлов и ограничения - принятые, см. README", "software": "cantera " + __import__("cantera").__version__ + ", python3", "date": __import__("datetime").date.today().isoformat()}); json.dump(out, open("gtu_study_results.json", "w"), ensure_ascii=False, indent=1)
for g in grid:
    print(f"{g['size']} {g['arch']:22s} PR={g['PR']} T41={g['T41']} КПД={g['eff']*100:.1f} уд.раб={g['sp']:.0f} EGT={g['EGT']:.0f} T3={g['T3']:.0f}")
