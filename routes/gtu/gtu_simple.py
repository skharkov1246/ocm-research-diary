import json
from scipy.optimize import brentq
from gtu_cycle import cycle, cycle_r, recup
import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))
out = {}
def show(tag, r, W1):
    r = dict(r); r["tag"] = tag; r["W1"] = W1
    r["eff_el"] = r["eff"]*0.985*0.972; r["P_el"] = r["P_MW"]*0.985*0.972
    out[tag] = r
    print(f"{tag:34s} P={r['P_MW']:.2f} МВт КПД={r['eff']*100:.1f}% (эл {r['eff_el']*100:.1f}%, {r['P_el']:.2f} МВт) уд.раб={r['sp_work_kJkg']:.0f} EGT={r['EGT_C']:.0f} T3={r['T3_C']:.0f} T4={r['T4_C']:.0f} T45={r['T45_C']:.0f} πгг={r['GG_ER']:.2f} πст={r['PT_ER']:.2f} Lк={r['wc_kJkg']:.0f} Gвых={r['W_exh']:.1f}" + (f" T3r={r['T3r_C']:.0f} Tух={r['Tstack_C']:.0f}" if 'T3r_C' in r else ""))
    return r
# --- проверка на прецеденте: OPRA OP16 (π=6.7 принято; паспорт: 1.876 МВт эл., 25.1 %) ---
S = dict(eta_pc=0.845, eta_pt_gg=0.86, eta_pt_pt=0.86, x_ngv=0.02, x_rot=0.02, x_leak=0.01, dp_comb=0.05, p_exh_ratio=1.03)
for t in (950, 1000, 1050):
    show(f"OP16-подобный одновальный T41={t}", cycle(6.7, t, W1=8.8, **S), 8.8)
# --- простой двухвальный: 1 центробежное колесо, 1 ступень ТГГ, 1 ступень СТ ---
W = 31.0
for PR in (5.5, 6.5, 7.5):
    for t in (1000, 1050, 1100):
        show(f"простой π={PR} T41={t}", cycle(PR, t, W1=W, **S), W)
# --- то же с рекуператором ---
for PR in (4.5, 5.5, 6.5):
    for t in (1000, 1050):
        for eps in (0.85,):
            show(f"простой+рекуп π={PR} T41={t} ε={eps}", recup(PR, t, eps, W1=W, **S), W)
# --- промежуточный: 2 центробежных колеса, охлаждаемая турбина 2 ступени, СТ 2 ступени ---
M = dict(eta_pc=0.85, eta_pt_gg=0.865, eta_pt_pt=0.88, x_ngv=0.04, x_rot=0.04, x_leak=0.015, dp_comb=0.045, p_exh_ratio=1.02)
for PR in (11, 13):
    for t in (1100, 1150):
        show(f"2 центробежных π={PR} T41={t}", cycle(PR, t, W1=W, **M), W)
out.update({"method": "0D-цикл Брайтона: свойства газа Cantera (GRI-Mech 3.0), метан, ISO; КПД узлов и ограничения - принятые, см. README", "software": "cantera " + __import__("cantera").__version__ + ", python3", "date": __import__("datetime").date.today().isoformat()}); json.dump(out, open("gtu_simple_results.json", "w"), ensure_ascii=False, indent=1)
