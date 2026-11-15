"""Блок 6. Вал и стык пакета: прямой привод на частоте пакета (одновальный, вал = вход редуктора или приводной машины) /
свободная силовая турбина на частоте пакета / (только генераторные пакеты) быстроходный одновальный - вал на своей частоте,
выход через существующий корпус редуктора пакета со сменной зубчатой парой.

Частота быстроходного варианта выбирается по коэффициенту расхода первой центробежной ступени: φ ∝ n² при заданном расходе и
окружной скорости (funnel_space.derive), n = n_пакета·√(φ_opt/φ), φ_opt - funnel_space.PHI_OPT. Окно отношения частот
HS_RATIO - оценка (межосевое расстояние существующего корпуса не меняется, меняются числа зубьев); проверить у изготовителя
редуктора. Номинал редуктора по моменту при большей входной частоте не превышается (та же мощность - меньше момент).
Plug&Play стыка вала - состояние модели воронки (plugplay_model) с поправкой: у быстроходного - «переходник».
"""
import math

from . import _track as T

DRIVES = {
    "прямой": "одновальный, вал двигателя = вход редуктора (генераторный) или приводной машины пакета",
    "свободная СТ": "газогенератор на своей частоте + свободная силовая турбина на частоте пакета",
    "быстроходный": "одновальный на своей частоте; выход через существующий корпус редуктора пакета со сменной зубчатой парой "
                    "(только генераторные пакеты)",
}
HS_RATIO = (0.70, 1.60)     # n_двигателя / n_пакета для сменной пары в существующем корпусе - оценка
HS_TRL = 7                  # смена зубчатой пары в существующем корпусе - оценка (практика изготовителей редукторов; не проверено)


def static_rules(C, tgt):
    d = C["shaft"]["drive"]
    R = []
    gen = tgt.get("variant") == "генератор"
    if d == "быстроходный" and not gen:
        R.append(dict(code="hs_generator_only", hard=True, ok=False,
                      text="быстроходный вариант через редуктор пакета - только генераторные пакеты (приводная машина вращается на своей "
                           "переменной частоте)"))
    if d == "прямой" and not gen:
        R.append(dict(code="single_shaft_md", hard=True, ok=False,
                      text="одновальный на приводе компрессора: частота приводной машины меняется с режимом - нужна свободная силовая "
                           "турбина (исходные приводные машины - двухвальные)"))
    if C["build"]["cores"] > 1 and d != "свободная СТ":
        R.append(dict(code="cluster_free_pt", hard=True, ok=False, text="кластер ядер работает на общую свободную силовую турбину"))
    return R


def n_highspeed(n_pkg, phi_at_pkg, typ):
    """частота быстроходного вала по φ первой центробежной ступени: n = n_пакета·√(φ_opt/φ), в окне HS_RATIO"""
    FS = T.FS()
    phi_opt = FS.PHI_OPT["Д" if typ == "Д" else "Ц"]
    r = math.sqrt(phi_opt / max(phi_at_pkg, 1e-6))
    r = min(HS_RATIO[1], max(HS_RATIO[0], r))
    return round(n_pkg * r / 10.0) * 10.0, r


def rules(C, tgt, n_engine):
    R = []
    if C["shaft"]["drive"] == "быстроходный":
        r = n_engine / tgt["n_pkg_true"]
        R.append(dict(code="hs_ratio", hard=True, ok=HS_RATIO[0] - 1e-9 <= r <= HS_RATIO[1] + 1e-9, val=round(r, 4), lim=list(HS_RATIO),
                      text="отношение частот сменной пары в существующем корпусе редуктора (оценка)"))
    return R


def plugplay_fix(C, tgt, items):
    """поправка состояния «вал и муфта» модели воронки для быстроходного варианта (там частота пакета подменена частотой вала)"""
    if C["shaft"]["drive"] == "быстроходный" and "вал и муфта" in items:
        items["вал и муфта"] = dict(state="переходник", value=0.6, weight=items["вал и муфта"]["weight"],
                                    note="сменная зубчатая пара в существующем корпусе редуктора пакета; фланец и муфта - "
                                         "переходник (оценка)")
    return items


def report(C, tgt, g, ch, G, P_shaft_kW):
    d = C["shaft"]["drive"]
    n_eng = g["n"] if g else None
    n_out = tgt["n_pkg_true"] if d != "быстроходный" else n_eng
    out = dict(variant=d, name=DRIVES[d], n_package_rpm=tgt["n_pkg_true"], n_engine_rpm=round(n_eng) if n_eng else None,
               rotors=C["rotors"])
    if n_out:
        w = 2 * math.pi * n_out / 60.0
        out["torque_out_Nm"] = round(P_shaft_kW * 1e3 / w, 1)
    if d == "быстроходный" and n_eng:
        out["gear_ratio_change"] = round(n_eng / tgt["n_pkg_true"], 4)
        out["torque_pkg_input_Nm"] = round(P_shaft_kW * 1e3 / (2 * math.pi * tgt["n_pkg_true"] / 60.0), 1)
    if g and g.get("pt"):
        out["pt"] = dict(D_m=round(g["pt"]["D"], 4), n_rpm=tgt["n_pkg_true"])
    out["bearings"] = C["build"]["bear"]
    out["case"] = C["build"]["case"]
    return out


def trl_items(C):
    if C["shaft"]["drive"] == "быстроходный":
        return [("сменная зубчатая пара в корпусе редуктора пакета", HS_TRL, "оценка: практика изготовителей редукторов, не проверено")]
    return []
