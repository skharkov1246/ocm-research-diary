"""Блок 4. Турбина газогенератора и силовая турбина: архитектуры воронки (осевая О / радиальная Р, число ступеней),
классы охлаждения U / N / K / P / C и материалы.

Модели воронки: КПД - funnel_space.quick_perf (осевая - eta_axial по нагрузке ψ и потеря с выходной скоростью; радиальная -
по U/C0 и быстроходности, funnel_l0.corrections), пределы - make_model (t41_limit, u_limit_turb, COOL_AIR, NGV_T41),
размеры и правила - funnel_space.derive / check. Число ступеней «auto» - наименьшее, при котором нагрузка в хорошем окне
(concept.PSI_AUTO, NU_AUTO; оценка по frame3). Охлаждение повторной ступени (подогрев) - тем же классом (оценка).
"""
import math

from . import _track as T
from .concept import NU_AUTO, PSI_AUTO

COOL = {
    "U": "без охлаждения (литые неохлаждаемые лопатки; точка KV-3)",
    "N": "охлаждаемый сопловой аппарат, рабочие лопатки без охлаждения",
    "K": "охлаждаемые сопловой аппарат и рабочие лопатки (монокристалл - до 1200 °C)",
    "P": "печатные каналы охлаждения у стенки (сопловой аппарат; колесо - если напечатано)",
    "C": "керамика (Si3N4 колесо, CMC статор) без охлаждения",
}


def cool_fracs(C):
    """доли воздуха охлаждения (сопловой аппарат, ротор) основной и повторной ступени - make_model.COOL_AIR"""
    MK = T.MK()
    t = C["turbine"]
    xn, xr = MK.COOL_AIR[t["cool"]]
    if not MK.rotor_cooled(t["cool"], t["make"]):
        xr = 0.0
    xn2 = xr2 = 0.0
    if C["hx"]["reheat"] and t["cool"] != "U":
        xn2, xr2 = xn, xr
    return xn, xr, xn2, xr2


def t41_limit(C):
    t = C["turbine"]
    return T.MK().t41_limit(t["type"], t["cool"], t["mat"], t["make"])


def static_rules(C):
    t = C["turbine"]
    FS = T.FS()
    R = []
    if t["mat"] not in FS.MAT_TW.get(t["make"], []):
        R.append(dict(code="turb_mat_make", hard=True, ok=False,
                      text=f"сплав {t['mat']} не изготавливается способом «{t['make']}» (funnel_space.MAT_TW)"))
    if t["cool"] == "C" and not (t["type"] == "Р" and t["mat"] == "Si3N4" and t["make"] == "керамика"):
        R.append(dict(code="cool_C", hard=True, ok=False,
                      text="класс C - керамическое радиальное колесо Si3N4 (CGT302; осевых керамических колёс в треке нет, "
                           "make_model.CERAMIC: только turb_radial)"))
    if t["cool"] != "C" and t["make"] == "керамика":
        R.append(dict(code="ceramic_class", hard=True, ok=False, text="керамическое колесо - только класс C"))
    if t["cool"] == "P" and t["hot_make"] != "печать":
        R.append(dict(code="cool_P_print", hard=True, ok=False,
                      text="класс P - печатные каналы в сопловом аппарате: горячая часть печатная (make_model.COOL_AIR_BASIS)"))
    if t["mat"] == "SC" and t["type"] != "О":
        R.append(dict(code="sc_axial", hard=True, ok=False, text="монокристаллические лопатки - осевая турбина (make_model.MAT['SC'])"))
    if t["hot_mat"] not in FS.MAT_HOT.get(t["hot_make"], []):
        R.append(dict(code="hot_mat_make", hard=True, ok=False,
                      text=f"сплав горячей части {t['hot_mat']} не изготавливается способом «{t['hot_make']}» (funnel_space.MAT_HOT)"))
    tl, tb = t41_limit(C)
    R.append(dict(code="T41_class", hard=True, ok=t["t41"] <= tl + 1e-6, val=t["t41"], lim=tl,
                  text=f"T41 по классу охлаждения {t['cool']} и сплаву {t['mat']}: {tb}"))
    if C["hx"]["reheat"]:
        R.append(dict(code="T41_reheat", hard=True, ok=C["hx"]["t_rh"] <= tl + 1e-6, val=C["hx"]["t_rh"], lim=tl,
                      text="газ после межтурбинной камеры - тот же предел класса охлаждения (оценка)"))
        if C["rotors"] == 1 and t["nt"] != "auto" and int(t["nt"]) < 2:
            R.append(dict(code="reheat_stages", hard=True, ok=False,
                          text="подогрев на одном валу делит расширение: нужны не меньше двух ступеней турбины"))
    return R


def nt_auto(C, w_turb, U_max, dhs_fn=None):
    """число ступеней турбины газогенератора «auto»: осевая - наименьшее с ψ = w/(nt·U²) <= PSI_AUTO при пределе U;
    радиальная - наименьшее с U/C0 >= NU_AUTO (C0 по перепаду ступени)"""
    t = C["turbine"]
    if t["nt"] != "auto":
        return int(t["nt"])
    lo = 2 if (C["hx"]["reheat"] and C["rotors"] == 1) else 1
    if t["type"] == "О":
        return max(lo, min(3, math.ceil(w_turb / (PSI_AUTO * U_max ** 2) - 1e-9)))
    for nt in range(lo, 3):
        dhs = dhs_fn(nt) if dhs_fn else w_turb / nt / 0.86
        if U_max / math.sqrt(2 * max(dhs, 1.0)) >= NU_AUTO:
            return nt
    return 2


def npt_auto(C, w_pt, U_pt):
    t = C["turbine"]
    if C["rotors"] != 2:
        return 0
    if t["npt"] != "auto":
        return int(t["npt"])
    return max(1, min(3, math.ceil(w_pt / (PSI_AUTO * U_pt ** 2) - 1e-9)))


def trl_items(C):
    return []


def report(C, ch, g, etas, nt, npt):
    t = C["turbine"]
    tl, tb = t41_limit(C)
    xn, xr, xn2, xr2 = cool_fracs(C)
    out = dict(variant=f"{t['type']}×{nt}" + (f" + СТ×{npt}" if C["rotors"] == 2 else ""), cool=t["cool"], cool_name=COOL[t["cool"]],
               mat=t["mat"], make=t["make"], hot=f"{t['hot_make']} {t['hot_mat']}", t41_C=t["t41"], t41_limit_C=tl, t41_basis=tb,
               x_cool=dict(ngv=xn, rotor=xr, ngv_rh=xn2, rotor_rh=xr2), eta_gg_p=round(etas["gg"], 4),
               eta_pt_p=round(etas["pt"], 4) if etas.get("pt") else None, nt=nt, npt=npt)
    if ch.get("ok"):
        out["sections"] = [dict(name=s["name"], T_in_C=round(s["T_in"] - 273.15, 1), T_out_C=round(s["T_out"] - 273.15, 1),
                                p_in_bar=round(s["p_in"] / 1e5, 3), p_out_bar=round(s["p_out"] / 1e5, 3), ER=round(s["p_in"] / s["p_out"], 4),
                                w_kJkg=round(s["w"] / 1e3, 2)) for s in ch["sections"]]
    if g:
        tg = g["tg"]
        out["geometry"] = dict(D_m=round(tg["D"], 4), Dtip_m=round(tg["Dtip"], 4), U_ms=round(tg["U"], 1),
                               h_last_m=round(tg["h"], 4) if tg.get("h") else None, psi=round(tg["psi"], 3) if tg.get("psi") else None,
                               nu=round(tg["nu"], 3) if tg.get("nu") else None, AN2=round(tg["AN2"], 0) if tg.get("AN2") else None)
        if g["pt"]:
            p = g["pt"]
            out["pt_geometry"] = dict(D_m=round(p["D"], 4), Dtip_m=round(p["Dtip"], 4), U_ms=round(p["U"], 1), h_last_m=round(p["h"], 4),
                                      psi=round(p["psi"], 3), AN2=round(p["AN2"], 0))
    return out
