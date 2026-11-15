"""Блок 2. Компрессор: архитектуры воронки (1Ц, 2Ц, Д, ОЦ), вход, материал и способ изготовления колёс;
промежуточное охлаждение воздух-воздух (без воды) с вентилятором.

КПД ступеней - модели воронки (funnel_space.eta_comp + поправки зазора и быстроходности funnel_l0.corrections по геометрии
derive), окружные скорости - make_model.u_limit_comp, диапазоны π - funnel_space.PI_RANGE. Промохладитель - параметры
перебора циклов (organic/cycles A: ε 0,80 к температуре наружного воздуха, потери 2 %, вентилятор: расход холодного воздуха
1,5 тепловой нагрузки, перепад 500 Па, КПД 0,65 × 0,95). Всё - оценка.
"""
from . import _track as T
from .concept import stage_types, stages_pi

ARCH = {
    "1Ц": "одно центробежное колесо (до π 7 - прецедент OPRA OP16, выше - трансзвуковое)",
    "2Ц": "два центробежных колеса (точка KV-3 frame3)",
    "Д": "одно диагональное колесо",
    "ОЦ": "осевые ступени (1-3) + центробежное колесо",
}


def ic_params():
    A = T.CY().A
    return dict(eps=float(A["eps_ic"]), dp=float(A["dp_ic"]), Cr=float(A["Cr_ic"]), dp_fan=float(A["dp_fan"]),
                eta_fan=float(A["eta_fan"]), eta_motor=float(A["eta_motor"]))


def ic_after(C):
    c = C["compressor"]
    if c["intercool"] == "нет":
        return ()
    return (int(c["ic_after"]),)


def pi_range(C):
    FS = T.FS()
    c = C["compressor"]
    return FS.PI_RANGE[(c["arch"], c["n_ax"])] if c["arch"] == "ОЦ" else FS.PI_RANGE[c["arch"]]


def static_rules(C):
    """правила совместимости, проверяемые по описанию концепции (код, жёсткое ли, причина)"""
    c = C["compressor"]
    R = []
    n = len(stage_types(C))
    if c["intercool"] != "нет":
        ok = n >= 2 and 0 <= c["ic_after"] < n - 1
        R.append(dict(code="ic_stages", hard=True, ok=ok,
                      text="промохлаждение ставится между ступенями: нужно не меньше двух ступеней сжатия (1Ц и Д - одна ступень)"))
    lo, hi = pi_range(C)
    R.append(dict(code="pi_range", hard=False, ok=lo <= c["pi"] <= hi, val=c["pi"], lim=[lo, hi],
                  text=f"π {c['arch']} в диапазоне воронки {lo:g}-{hi:g} (funnel_space.PI_RANGE, оценка)"))
    if c["mat"] == "Ni" and C["build"]["rotor_join"] not in ("заодно", "сварка"):
        R.append(dict(code="ni_wheel_join", hard=True, ok=False,
                      text="никелевое колесо компрессора - только заодно или сваркой с турбиной (make_model.MAT['Ni'])"))
    if c["inlet"] == "двух" and C["build"]["cores"] > 1:
        R.append(dict(code="dbl_inlet_cluster", hard=True, ok=False, text="кластер ядер - односторонний вход (шаблон воронки)"))
    return R


def trl_items(C):
    out = []
    if C["compressor"]["intercool"] != "нет":
        CY = T.CY()
        out.append(("промохлаждение воздух-воздух", CY.TRL["ic"], CY.TRL_NOTE["ic"]))
    return out


def report(C, ch, g, perf):
    c = C["compressor"]
    st = []
    for k, s in enumerate(ch.get("stages", [])):
        row = dict(k=k + 1, type=stage_types(C)[k], pi=round(s["pi"], 4), eta_p=round(s["eta"], 4), T2_C=round(s["T2"] - 273.15, 1),
                   w_kJkg=round(s["w"] / 1e3, 2))
        if g:
            gs = g["stages"][k]
            row.update(D2_m=round(gs["D"], 4), U2_ms=round(gs["U"], 1), mat=gs["mat"],
                       phi=round(gs["phi"], 4) if "phi" in gs else None, b2_mm=round(gs["b2"] * 1e3, 1) if "b2" in gs else None,
                       htr=round(gs["htr"], 3) if "htr" in gs else None)
        st.append(row)
    out = dict(variant=f"{c['arch']}{c['n_ax'] if c['arch'] == 'ОЦ' else ''}", name=ARCH[c["arch"]], inlet=c["inlet"], mat=c["mat"],
               make=c["make"], pi=c["pi"], stages=st, intercool=c["intercool"])
    if ch.get("ok"):
        out.update(T2_C=round(ch["T2"] - 273.15, 1), p2_bar=round(ch["p2"] / 1e5, 3), wc_kJkg=round(ch["wc"] / 1e3, 2))
        if ch.get("ic"):
            ic = ch["ic"][0]
            out["ic"] = dict(after_stage=ic["after"] + 1, T_in_C=round(ic["T_in"] - 273.15, 1), T_out_C=round(ic["T_out"] - 273.15, 1),
                             Q_kJ_per_kg_air=round(ic["Q"] / 1e3, 2), W_fan_kJ_per_kg_air=round(ic["W_fan"] / 1e3, 3),
                             params=ic_params(), source="organic/cycles: A eps_ic, dp_ic, Cr_ic, dp_fan, eta_fan, eta_motor")
    if g:
        out["n_rpm"] = round(g["n"])
    return out


def stage_list(C, etas):
    """[(π, η_p)] для цепочки"""
    return list(zip(stages_pi(C), etas))
