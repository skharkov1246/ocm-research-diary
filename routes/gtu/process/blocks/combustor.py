"""Блок 3. Камера сгорания: бедная предварительно перемешанная с пассивным пилотом (LPP; решение владельца 15.11.2026 -
один клапан на топливо, частичную нагрузку держит пассивный элемент), FLOX, захваченный вихрь (TVC), каталитическая (CAT),
горение с ростом давления - детонационная (RDC) и волновой ротор (ВР) - со штрафом готовности.

Числа (каждое - из расчёта или источника):
  * потери давления - итог воронки камеры уровня 0 (combustor/l0_aws, 262 144 варианта KV-3): лучший по рейтингу вариант пары
    «схема × тип камеры» (feed_to_engine_funnel.best, dp_source = comb_l0) или наименьшие потери схемы среди допустимых
    (by_scheme.dp_min, comb_l0_min); track - принятые 4 % (кольцевая) / 5 % (выносная), как gtu_cycle и funnel_space._core;
  * горение с ростом давления - модель перебора циклов: эквивалентное располагаемое давление при сгорании при постоянном
    объёме (cycles.pgc_gain, Cantera), реализация 25 % и штраф КПД турбины 2 п. п. (cycles A pgc_f_real, pgc_deta_t);
    волновой ротор - та же модель (своей модели нет - оценка);
  * выбросы и пределы - combustor/targets.json (NOx 9 и CO 15 ppm при 15 % O2; адиабатическая температура зоны пламени
    1750-1850 K; один клапан на топливо) и NOx, CO лучшего варианта пары из уровня 0 камеры; ни один вариант уровня 0 камеры
    не выполняет все жёсткие цели (hard_ok 0) - выбросы показываются как ограничение, которое ещё не закрыто;
  * готовность - comb_space.TRL_SCHEME (STATE_OF_ART.md, оценка), для RDC и ВР - не выше make_model.TRL (3).
"""
import math

from . import _track as T

SCHEME = {
    "LPP": dict(name="бедная предварительно перемешанная с пассивным пилотом (один клапан)", premix=True, pgc=False, mk_trl=None),
    "FLOX": dict(name="FLOX / MILD: струйная стабилизация с рециркуляцией продуктов", premix=True, pgc=False, mk_trl=None),
    "TVC": dict(name="захваченный вихрь (полость с вихрем - пассивный пилот)", premix=True, pgc=False, mk_trl=None),
    "CAT": dict(name="каталитическая (штраф готовности)", premix=False, pgc=False, mk_trl=None),
    "RDC": dict(name="вращающаяся детонация, рост давления (штраф готовности)", premix=False, pgc=True, mk_trl="rdc"),
    "ВР": dict(name="волновой ротор, рост давления (штраф готовности)", premix=False, pgc=True, mk_trl="волновой_ротор"),
}
DP_TRACK = {"кольцевая": 0.04, "выносная": 0.05}      # gtu_cycle.cycle dp_comb = 0,04; funnel_space._core: выносная 5 %
T3_PREMIX_MAX = 650.0     # °C, воздух на входе камеры предварительного смешения - оценка (самовоспламенение и проскок в смесителе;
#                           прецедент - Mercury 50: рекуператор и бедная камера, NOx 5 ppm, паспорт fleet; температура на входе
#                           его камеры не опубликована)
T3_CAT_MIN = 723.0        # K, вход катализатора без предгорелки (combustor/README: бедная предгорелка до 723 K - порог, оценка)


def l0():
    return T.load_json(T.F_COMB_L0)


def targets():
    d = T.load_json(T.F_COMB_TARGETS)
    tg = {t["id"]: t for t in d["targets"]}
    tad = next(x for x in tg["nox_gas"]["derived"] if x["id"] == "flame_zone_Tad")["target"]
    return dict(nox=tg["nox_gas"]["target"]["max"], co=tg["co_gas"]["target"]["max"], tad_min=tad["min"], tad_max=tad["max"],
                valves=tg["fuel_valves"]["target"]["control_valves_per_fuel"], dp_obj=tg["pressure_loss"]["target"]["max"])


def best_pair(C):
    cb = C["combustor"]
    key = f"{cb['scheme']}|{cb['chamber']}"
    return key, l0()["feed_to_engine_funnel"]["best"].get(key)


def dp_comb(C):
    """потери давления камеры (доля p3) и источник"""
    cb = C["combustor"]
    if cb.get("dp") is not None:
        return float(cb["dp"]), "задано в концепции"
    if SCHEME[cb["scheme"]]["pgc"]:
        A = T.CY().A
        return float(A["dp_comb"]), "organic/cycles A['dp_comb'] (подвод к камере с ростом давления)"
    src = cb["dp_source"]
    if src == "track":
        return DP_TRACK[cb["chamber"]], "принятая величина трека (gtu_cycle dp_comb 0,04; funnel_space._core выносная 0,05)"
    if src == "comb_l0_min":
        v = l0()["by_scheme"][cb["scheme"]]["dp_min"]
        return float(v), f"combustor/l0_aws: by_scheme[{cb['scheme']}].dp_min (наименьшие потери среди допустимых, KV-3)"
    key, b = best_pair(C)
    if b is None:
        return DP_TRACK[cb["chamber"]], f"пары {key} нет в уровне 0 камеры - принятая величина трека"
    return float(b["dp_comb"]), f"combustor/l0_aws: feed_to_engine_funnel.best[{key}].dp_comb (лучший по рейтингу, KV-3)"


def pgc_params(C):
    A = T.CY().A
    f = C["combustor"].get("pgc_real")
    return dict(f_real=float(A["pgc_f_real"] if f is None else f), deta_t=float(A["pgc_deta_t"]), dp=float(A["dp_comb"]))


def PRb(C, t41_C):
    """коэффициент полного давления камеры p4/p3: число или функция (T3, p3, mc) для горения с ростом давления"""
    dp, _ = dp_comb(C)
    if not SCHEME[C["combustor"]["scheme"]]["pgc"]:
        return 1.0 - dp
    CY, E = T.CY(), T.E()
    pp = pgc_params(C)
    T41 = t41_C + 273.15
    h_f = CY.hs(CY.St(1.0, CY.T0, CY.P0, CY.Y_CH4))[0]
    cache = {}

    def f(T3, p3, mc):
        key = (round(T3, 1), round(p3, -1), round(mc, 5))
        if key not in cache:
            mf = mc * E.far_for(T3, T41)
            PG, _ = CY.pgc_gain(CY.St(mc, T3, p3, CY.Y_AIR), CY.Y_CH4, h_f, mf, CY.LHV_CH4)
            cache[key] = (1 - pp["dp"]) + pp["f_real"] * (PG - (1 - pp["dp"]))
            cache["PG"] = PG
        return cache[key]
    f.cache = cache
    return f


def trl(C):
    s = C["combustor"]["scheme"]
    t = T.CS().TRL_SCHEME.get("LPP" if s == "LPP" else s)
    basis = f"comb_space.TRL_SCHEME[{s}] (STATE_OF_ART.md, оценка)" if t is not None else ""
    mk = SCHEME[s]["mk_trl"]
    if mk:
        tm = T.MK().TRL[mk]["trl"]
        if t is None or tm < t:
            t, basis = tm, f"make_model.TRL[{mk}]: {T.MK().TRL[mk]['basis']}"
    return t, basis


def rules(C, ch):
    """правила камеры по результату цепочки: (код, жёсткое, ok, величина, предел, причина, группа)"""
    s = C["combustor"]["scheme"]
    sp = SCHEME[s]
    tg = targets()
    R = []
    if not ch.get("ok"):
        return R
    T3C = ch["T3"] - 273.15
    if sp["premix"]:
        R.append(dict(code="comb_T4_flame", hard=True, ok=ch["T4"] <= tg["tad_max"], val=round(ch["T4"], 1), lim=tg["tad_max"],
                      group="выбросы", text="газ за камерой не горячее верхней границы адиабатической температуры зоны пламени "
                                            "(targets.json nox_gas.derived flame_zone_Tad): иначе NOx больше 9 ppm"))
        if s in ("LPP", "TVC"):
            R.append(dict(code="comb_T3_premix", hard=False, ok=T3C <= T3_PREMIX_MAX, val=round(T3C, 1), lim=T3_PREMIX_MAX,
                          group="камера", text="воздух на входе смесителя предварительного смешения (самовоспламенение, проскок; оценка)"))
    if s == "CAT":
        ok = ch["T3"] >= T3_CAT_MIN
        R.append(dict(code="comb_cat_preburner", hard=True, ok=ok, val=round(ch["T3"], 1), lim=T3_CAT_MIN, group="камера",
                      text="катализатору нужен вход не холоднее 723 K: иначе предгорелка со своим клапаном - против цели "
                           "fuel_valves (один клапан на топливо, решение владельца 15.11.2026)"))
    if C["hx"]["reheat"]:
        R.append(dict(code="fuel_valves_reheat", hard=False, ok=False, val=2, lim=tg["valves"], group="камера",
                      text="повторный подогрев - вторая камера: топливо делится между камерами; без второго регулирующего клапана - "
                           "только пассивным делением (дроссель), не проверено"))
    key, b = best_pair(C)
    if b is not None:
        R.append(dict(code="nox_gas", hard=False, ok=b["NOx_100"] <= tg["nox"], val=round(b["NOx_100"], 2), lim=tg["nox"], group="выбросы",
                      text=f"NOx на 100 % лучшего варианта {key} уровня 0 камеры (targets.json nox_gas 9 ppm)"))
        R.append(dict(code="co_gas", hard=False, ok=b["CO_100"] <= tg["co"], val=round(b["CO_100"], 2), lim=tg["co"], group="выбросы",
                      text=f"CO на 100 % лучшего варианта {key} уровня 0 камеры (targets.json co_gas 15 ppm)"))
    return R


def report(C, ch, PRb_obj):
    cb = C["combustor"]
    dp, src = dp_comb(C)
    t, tb = trl(C)
    key, b = best_pair(C)
    sch = l0()["by_scheme"].get(cb["scheme"])
    out = dict(variant=cb["scheme"], name=SCHEME[cb["scheme"]]["name"], chamber=cb["chamber"], dp=round(dp, 4), dp_source=src,
               trl=t, trl_basis=tb)
    if ch.get("ok"):
        out.update(T3_C=round(ch["T3"] - 273.15, 1), T4_C=round(ch["T4"] - 273.15, 1), PRb=round(ch["PRb"], 4))
    if callable(PRb_obj) and "PG" in PRb_obj.cache:
        out["pgc"] = dict(PG_ideal=round(PRb_obj.cache["PG"], 4), **pgc_params(C),
                          source="cycles.pgc_gain (Cantera, сгорание при постоянном объёме) и cycles A pgc_f_real, pgc_deta_t")
    if b is not None:
        out["comb_l0_best"] = dict(key=key, idx=b.get("idx"), dp_comb=b["dp_comb"], D_m=b["D_comb"], L_m=b["L_comb"], mass_kg=b["m_comb"],
                                   parts=b["parts"], NOx_ppm=[b["NOx_100"], b["NOx_50"], b.get("NOx_30")],
                                   CO_ppm=[b["CO_100"], b["CO_50"], b.get("CO_30")], T_metal_max_K=b["T_metal_max_K"],
                                   hard_ok=b["hard_ok"], source="combustor/l0_aws/comb_l0_KV-3_n262144_results.json")
    if sch:
        out["comb_l0_scheme"] = dict(variants=sch["total"], feasible=sch["ok"], hard_ok=sch["hard_ok"], dp_min=sch.get("dp_min"),
                                     parts_min=sch.get("parts_min"), NOx_100_median=sch.get("NOx_100_median"))
    tg = targets()
    if SCHEME[cb["scheme"]]["premix"]:
        out["T4_max_K"] = tg["tad_max"]
        out["T4_max_basis"] = ("газ за камерой не горячее верхней границы адиабатической температуры зоны пламени 1850 K "
                               "(targets.json nox_gas.derived flame_zone_Tad); T41 после охлаждения соплового аппарата - ниже")
    elif cb["scheme"] == "CAT":
        out["T4_max_K"] = None
        out["T4_max_basis"] = ("катализатор до 1173 K, дальше гомогенное догорание (combustor/README, пороги - оценка); предел T41 - "
                               "по турбине")
    else:
        out["T4_max_K"] = None
        out["T4_max_basis"] = "предел T41 - по турбине (класс охлаждения)"
    out["targets"] = tg
    return out


def eta_t_penalty(C):
    return pgc_params(C)["deta_t"] if SCHEME[C["combustor"]["scheme"]]["pgc"] else 0.0
