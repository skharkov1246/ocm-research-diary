"""Опорные точки и сверки библиотеки блоков.

  * frame3_concept() - расчётная точка KV-3 (frame3_point_results.json): одновальный, 2Ц (π 4,49 × 3,007, двусторонний вход),
    коэффициенты работы и cm2/U2 ступеней frame3, осевая турбина 3 ступени при U 420 м/с, класс U, T41 1000 °C, π 13,5;
  * check_frame3() - (а) evaluate против funnel_l0.solve той же концепции воронки (та же модель - совпадение до округления);
    (б) против frame3 (КПД на валу, удельная работа, выхлоп) - в пределах калибровки трека ±2 п. п.;
    (в) с КПД узлов frame3 (eta.source = задано) против gtu_cycle.cycle (Cantera напрямую);
  * CYCLES_KEYS / cycles_concept() / check_cycles() - топологии перебора циклов: концепция с допущениями перебора
    (КПД ступеней cycles.eta_stage, турбины A['eta_t'], без охлаждения, потери камеры 4 %) против cycles_results.json;
  * check_recup_gtu_cycle() - рекуператор: цепочка (температурная запись) против gtu_cycle.recup.
"""
import math

from . import _track as T
from .evaluator import evaluate

CYCLES_KEYS = [
    "R00-ic-rh-ibc-none-none",      # простой
    "R90-ic-rh-ibc-none-none",      # рекуператор 0,90
    "R90-ic-RH-ibc-none-none",      # рекуператор 0,90 + подогрев (рекомендация перебора)
    "R90-IC-rh-ibc-none-none",      # рекуператор 0,90 + промохлаждение
    "R00-ic-rh-ibc-none-fh",        # подогрев топлива
]
EXTRA_KEYS = ["R00-ic-rh-IBC-none-none", "R00-ic-rh-ibc-real-none", "R00-ic-rh-ibc-none-steam"]


def frame3_concept(eta_given=False):
    f3 = T.load_json(T.F_FRAME3)
    s1, s2 = f3["stage1"], f3["stage2"]
    pi, pi1 = float(f3["split"]["pi"]), float(f3["split"]["pi1"])
    C = dict(name="KV-3 frame3: одновальный 2Ц, осевая 3 ст., класс U", kv="KV-3",
             compressor=dict(arch="2Ц", inlet="двух", pi=pi, split=math.log(pi1) / math.log(pi), mat="Ti64", make="поковка",
                             mu=[float(s1["work_coeff_lambda"]), float(s2["work_coeff_lambda"])],
                             kcm=0.5 * (float(s1["phi2_cm2_U2"]) + float(s2["phi2_cm2_U2"])), zc=int(s1["blades_full"])),
             turbine=dict(type="О", nt=int(f3["turbine"]["selected_z"]), cool="U", mat="SC", make="литьё", t41=1000.0, lt=1.2, phit=0.5),
             combustor=dict(scheme="LPP", chamber="кольцевая", dp_source="track"),
             build=dict(case="разъём", rotor_join="сборный", bear="масло", conn="болты", comb_make="обработка", case_make="литьё"),
             shaft=dict(drive="прямой"))
    if eta_given:
        cy = f3["cycle"]["frame3"]
        C["eta"] = dict(source="задано", comp=[float(cy["eta_pc"])] * 2, gg=float(cy["eta_t"]))
    return C


def check_frame3():
    FS, L0 = T.FS(), T.L0()
    from .concept import normalize, to_funnel
    f3 = T.load_json(T.F_FRAME3)["cycle"]["frame3"]
    C = frame3_concept()
    r = evaluate(C)
    Cn = normalize(C)
    tgt = T.package_target("KV-3")
    fc = to_funnel(Cn, tgt, nt=Cn["turbine"]["nt"], npt=0)
    sol = L0.solve(fc, tgt)
    p = sol["perf"]
    rg = evaluate(frame3_concept(eta_given=True))
    gc = T.GC().cycle(float(T.load_json(T.F_FRAME3)["split"]["pi"]), 1000.0, eta_pc=float(f3["eta_pc"]), eta_pt_gg=float(f3["eta_t"]),
                      eta_pt_pt=float(f3["eta_t"]), x_ngv=0.0, x_rot=0.0, x_leak=0.015, dp_comb=0.04, dp_duct=0.0)
    return dict(
        evaluate=dict(eff_shaft=r["eff_shaft"], w_kJkg=r["w_kJkg"], G_air=r["G_air"], EGT_C=r["EGT_C"]),
        funnel=dict(eff_shaft=round(p["eta_th"] * 100, 3), w_kJkg=round(p["w"], 3), G_air=round(p["G"], 4), EGT_C=round(p["EGT"], 2)),
        frame3=dict(eff_shaft=round(f3["eff_shaft"] * 100, 3), w_kJkg=f3["w_kJkg"], EGT_C=f3["EGT_C"], P_shaft_kW=f3["P_shaft_kW"],
                    G=25.55, source="frame3/frame3_point_results.json cycle.frame3"),
        eta_given=dict(eff_shaft=rg["eff_shaft"], w_kJkg=rg["w_kJkg"], EGT_C=rg["EGT_C"]),
        gtu_cycle=dict(eff_shaft=round(gc["eff"] * 100, 3), w_kJkg=round(gc["sp_work_kJkg"], 3), EGT_C=round(gc["EGT_C"], 2),
                       call="gtu_cycle.cycle(13.5, 1000, eta_pc, eta_t, eta_t, x_ngv=0, x_rot=0, x_leak=0.015, dp_comb=0.04, dp_duct=0)"),
        d_eff_vs_funnel_pp=round(r["eff_shaft"] - p["eta_th"] * 100, 4),
        d_eff_vs_frame3_pp=round(r["eff_shaft"] - f3["eff_shaft"] * 100, 3),
        d_w_vs_frame3_pct=round((r["w_kJkg"] / f3["w_kJkg"] - 1) * 100, 3),
        d_EGT_vs_frame3_K=round(r["EGT_C"] - f3["EGT_C"], 2),
        d_eff_given_vs_gtu_cycle_pp=round(rg["eff_shaft"] - gc["eff"] * 100, 4),
        d_w_given_vs_gtu_cycle_pct=round((rg["w_kJkg"] / gc["sp_work_kJkg"] - 1) * 100, 3),
        calibration_pp=2.0, record=r)


def _cycles_row(key):
    d = T.load_json(T.F_CYCLES)
    return next(x for x in d["results"] if x["key"] == key)


def cycles_concept(key, tag="opt"):
    """концепция с допущениями перебора циклов для топологии key в точке tag (opt / pr135) из cycles_results.json"""
    CY = T.CY()
    row = _cycles_row(key)
    o = row[tag]
    t = o["topo"]
    PR = float(o["PR"])
    prm = o["prm"]
    if t["ic"]:
        nc, arch = 2, "2Ц"
    else:
        nc = CY.n_stages(PR)
        arch = {1: "1Ц", 2: "2Ц"}[nc]
    etas = [CY.eta_stage(PR ** (1.0 / nc))] * nc
    scheme = "RDC" if t["pgc"] == "real" else "LPP"
    prep = {"none": "нет", "fh": "подогрев", "steam": "паровой риформинг", "egr": "сухой риформинг"}[t["tcr"]]
    C = dict(name=f"перебор циклов: {row['name']} ({key}, {tag})", kv="KV-3",
             compressor=dict(arch=arch, inlet="одно", pi=PR, split=0.5, intercool="воздух" if t["ic"] else "нет", ic_after=0),
             turbine=dict(type="О", nt="auto", cool="U", mat="SC", make="литьё", t41=float(CY.A["T41_C"])),
             combustor=dict(scheme=scheme, chamber="кольцевая", dp_source="track"),
             hx=dict(recup=float(t["rec"]), recup_make="фольга", reheat=bool(t["rh"]), x_rh=float(prm.get("x_rh", 0.5)), ibc=bool(t["ibc"])),
             fuel=dict(prep=prep, sc=float(prm.get("sc", 2.0)), R_egr=float(prm.get("R_egr", 1.0))),
             shaft=dict(drive="прямой"), build=dict(case="картридж"),
             eta=dict(source="задано", comp=etas, gg=float(CY.A["eta_t"])),
             options=dict(rec_model="enthalpy"))
    return C, o


def check_cycles(keys=None, tag="opt"):
    rows = []
    for key in keys or CYCLES_KEYS:
        C, o = cycles_concept(key, tag)
        r = evaluate(C)
        rows.append(dict(key=key, name=C["name"], PR=o["PR"],
                         cycles=dict(eff_shaft=round(o["eff"] * 100, 3), eff_el=round(o["eff_el"] * 100, 3), W_air=o["W_air"],
                                     W_exh=o["W_exh"], EGT_C=o["EGT_C"], T_rec_in_C=o.get("T_rec_in_C"), w_kJkg=o["w_kJkg"]),
                         blocks=dict(eff_shaft=r.get("eff_shaft"), eff_el=r.get("eff_el"), W_air=r.get("G_air"), W_exh=r.get("W_exh"),
                                     EGT_C=r.get("EGT_C"), T_rec_in_C=r.get("T_rec_in_C"), w_kJkg=r.get("w_kJkg")),
                         d_eff_pp=round(r["eff_shaft"] - o["eff"] * 100, 3) if r.get("eff_shaft") is not None else None,
                         d_W_air_pct=round((r["G_air"] / o["W_air"] - 1) * 100, 3) if r.get("G_air") else None,
                         d_EGT_K=round(r["EGT_C"] - o["EGT_C"], 2) if r.get("EGT_C") is not None else None))
    return rows


def check_recup_gtu_cycle(PR=6.5, T41=1000.0, eps=0.85, eta_c=0.86, eta_t=0.875):
    """рекуператор: цепочка (температурная запись, как gtu_cycle.recup) против gtu_cycle.recup при одинаковых КПД"""
    gc = T.GC().recup(PR, T41, eps, eta_pc=eta_c, eta_pt_gg=eta_t, eta_pt_pt=eta_t, x_ngv=0.0, x_rot=0.0, x_leak=0.015, dp_duct=0.0)
    C = dict(name="сверка с gtu_cycle.recup", compressor=dict(arch="2Ц", pi=PR, split=0.5), turbine=dict(t41=T41, nt=2),
             combustor=dict(dp_source="track"), hx=dict(recup=eps), eta=dict(source="задано", comp=[eta_c, eta_c], gg=eta_t),
             options=dict(rec_model="temperature"))
    r = evaluate(C)
    return dict(PR=PR, T41=T41, eps=eps, gtu_cycle=dict(eff=round(gc["eff"] * 100, 3), w=round(gc["sp_work_kJkg"], 3),
                                                        T3r_C=round(gc["T3r_C"], 2), EGT_C=round(gc["EGT_C"], 2)),
                blocks=dict(eff=r["eff_shaft"], w=r["w_kJkg"], T3r_C=r["T3_C"], T_rec_in_C=r["T_rec_in_C"]),
                d_eff_pp=round(r["eff_shaft"] - gc["eff"] * 100, 3), d_w_pct=round((r["w_kJkg"] / gc["sp_work_kJkg"] - 1) * 100, 3),
                d_T3r_K=round(r["T3_C"] - gc["T3r_C"], 2),
                note="gtu_cycle.recup: EGT_C - температура газа за турбиной (на входе рекуператора); T3r - воздух за рекуператором")


# =========================================================================== примеры: каждый вариант блока на точке KV-3 (+ KV-1, KV-6)
REC90_RH = dict(compressor=dict(arch="2Ц", pi=8.2), hx=dict(recup=0.90, recup_make="фольга", reheat=True, x_rh=0.35))
EXAMPLES = {
    "опорная точка frame3": None,                                     # frame3_concept()
    "простой 2Ц π 13,5, LPP": dict(compressor=dict(arch="2Ц", inlet="двух", pi=13.5)),
    "+ подогрев топлива": dict(compressor=dict(arch="2Ц", inlet="двух", pi=13.5), fuel=dict(prep="подогрев")),
    "+ паровой риформинг (S/C 2)": dict(compressor=dict(arch="2Ц", inlet="двух", pi=13.5), fuel=dict(prep="паровой риформинг", sc=2.0)),
    "+ сухой риформинг (R 1)": dict(compressor=dict(arch="2Ц", inlet="двух", pi=13.5), fuel=dict(prep="сухой риформинг", R_egr=1.0)),
    "промохлаждение + рекуператор 0,90": dict(compressor=dict(arch="2Ц", pi=6.0, intercool="воздух"), hx=dict(recup=0.90)),
    "рекуператор 0,90 (фольга) + подогрев, π 8,2": REC90_RH,
    "рекуператор 0,90 (печать) + подогрев, π 8,2": dict(REC90_RH, hx=dict(recup=0.90, recup_make="печать", reheat=True, x_rh=0.35)),
    "рекуператор 0,90 + подогрев, π 13,5 (как есть)": dict(compressor=dict(arch="2Ц", inlet="двух", pi=13.5),
                                                          hx=dict(recup=0.90, reheat=True, x_rh=0.35)),
    "обращённый Брайтон (IBC)": dict(compressor=dict(arch="2Ц", inlet="двух", pi=13.5), hx=dict(ibc=True)),
    "FLOX + рекуператор 0,88": dict(compressor=dict(arch="2Ц", pi=8.0), combustor=dict(scheme="FLOX"), hx=dict(recup=0.88)),
    "TVC (полость - пассивный пилот)": dict(compressor=dict(arch="2Ц", inlet="двух", pi=13.5), combustor=dict(scheme="TVC")),
    "каталитическая + рекуператор 0,85": dict(compressor=dict(arch="2Ц", pi=8.0), combustor=dict(scheme="CAT"), hx=dict(recup=0.85)),
    "детонационная (RDC), класс N": dict(compressor=dict(arch="2Ц", inlet="двух", pi=13.5), combustor=dict(scheme="RDC"),
                                       turbine=dict(cool="N", t41=1050)),
    "волновой ротор, класс N": dict(compressor=dict(arch="2Ц", inlet="двух", pi=13.5), combustor=dict(scheme="ВР"),
                                   turbine=dict(cool="N", t41=1050)),
    "быстроходный через редуктор пакета": dict(compressor=dict(arch="2Ц", inlet="двух", pi=13.5), shaft=dict(drive="быстроходный")),
    "свободная СТ": dict(compressor=dict(arch="2Ц", pi=13.5), shaft=dict(drive="свободная СТ")),
    "свободная СТ + межтурбинная камера, класс K": dict(compressor=dict(arch="2Ц", pi=13.5), shaft=dict(drive="свободная СТ"),
                                                      turbine=dict(cool="K", t41=1150), hx=dict(reheat=True, t_rh=1000)),
    "класс K, T41 1150, π 16": dict(compressor=dict(arch="2Ц", pi=16.0), turbine=dict(cool="K", t41=1150)),
    "печатный класс P, оболочка": dict(compressor=dict(arch="2Ц", pi=13.5),
                                      turbine=dict(cool="P", mat="CM247LC", make="печать", hot_make="печать", hot_mat="Haynes282", t41=1150),
                                      build=dict(case="оболочка")),
    "керамика C, быстроходный": dict(compressor=dict(arch="1Ц", pi=6.0), shaft=dict(drive="быстроходный"),
                                    turbine=dict(type="Р", cool="C", mat="Si3N4", make="керамика", hot_make="керамика", hot_mat="CMC",
                                                 t41=1250)),
    "ОЦ (1 осевая + центробежная), класс K": dict(compressor=dict(arch="ОЦ", n_ax=1, pi=12.0), turbine=dict(cool="K", t41=1150)),
    "кластер 2 ядра + СТ": dict(build=dict(cores=2), compressor=dict(arch="2Ц", pi=13.5), shaft=dict(drive="свободная СТ"),
                                turbine=dict(cool="N", t41=1050)),
    "KV-1: 1Ц радиальная, быстроходный": dict(kv="KV-1", compressor=dict(arch="1Ц", pi=6.0),
                                             turbine=dict(type="Р", mat="K418", make="литьё", t41=950), shaft=dict(drive="быстроходный")),
    "KV-6: ОЦ (2 осевые + центробежная), свободная СТ, класс K": dict(kv="KV-6", compressor=dict(arch="ОЦ", n_ax=2, pi=16.0), shaft=dict(drive="свободная СТ"),
                                                    turbine=dict(cool="K", t41=1200)),
}


def example(name):
    c = EXAMPLES[name]
    return frame3_concept() if c is None else dict(c, name=name)
