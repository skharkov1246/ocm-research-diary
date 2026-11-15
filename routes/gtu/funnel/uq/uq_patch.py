"""Устойчивый отбор: подстановка одной реализации неопределённостей в модели уровня 1 без правки их файлов. Всё - оценка.

Код воронки (funnel_l1, meanline_l1, strength_l1, offdesign_l1) не меняется: на время одной выборки в памяти процесса заменяются
модульные величины и функции, к которым модели обращаются по имени модуля (ML.cc_shift, ST._uts, OD.TIP_INC_STALL, ...), а после -
восстанавливаются (Perturb - контекст; проверка восстановления - test_uq.py). Порядок внутри выборки:
  1) коэффициенты формы потерь (uq_factors.SHAPE) -> пересчёт привязки уровня КПД meanline_l1 (calibration_summary): уровень остаётся
     привязанным к формулам трека, меняется только зависимость от геометрии вдали от опорных точек;
  2) затем уровень КПД (сдвиги d_eta_c, d_eta_t), допуски изготовления (зазоры, шероховатость, кромка), свойства материалов,
     температура газа, признаки срыва - поверх пересчитанной привязки (на привязку они не влияют).
Шероховатость зависит от способа изготовления детали варианта - set_variant(c) перед расчётом варианта.
"""
import math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
FUNNEL = os.path.dirname(HERE)
for p in (FUNNEL, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)
import meanline_l1 as ML         # noqa: E402
import strength_l1 as ST         # noqa: E402
import offdesign_l1 as OD        # noqa: E402
import uq_factors as UF          # noqa: E402

_ORIG = {}
_S = {}                          # текущие значения факторов и величины варианта (читают обёртки)
_ACTIVE = [False]


def _capture():
    if _ORIG:
        return
    ML.calibration_summary()     # номинальная привязка - в кэше до первой подстановки (её же восстанавливаем)
    _ORIG.update(cc_shift=ML.cc_shift, eta_cap_l0=ML.eta_cap_l0, ax_shift=ML.ax_shift, rad_shift=ML.rad_shift, cc_losses=ML.cc_losses,
                 cc_stage=ML.cc_stage, ax_turbine=ML.ax_turbine, rad_stage=ML.rad_stage, ko_defaults=ML.ko_Y.__defaults__,
                 A=dict(ML.A), K_DUMP=ML.K_DUMP, SH=dict(ML._SH), _uts=ST._uts, _ys=ST._ys, sigma_allow=ST.sigma_allow,
                 rupture_hours=ST.rupture_hours, K_BURST=ST.K_BURST, DT_COOL=ST.DT_COOL, WEIBULL_M=ST.WEIBULL_M,
                 TIP_INC_STALL=OD.TIP_INC_STALL, DINC_STALL=OD.DINC_STALL)
    _S.update(UF.nominal(), kf_c=1.0, kf_t=1.0, d_eta_rough_t=0.0)


# =========================================================================== обёртки (читают _S)
def w_cc_shift(PRs=None):
    return (_ORIG["cc_shift"](PRs) if PRs is not None else _ORIG["cc_shift"]()) + _S["d_eta_c"]


def w_eta_cap_l0(*a, **k):
    return _ORIG["eta_cap_l0"](*a, **k) + _S["d_eta_c"]


def w_ax_shift():
    return _ORIG["ax_shift"]() + _S["d_eta_t"] + _S["d_eta_rough_t"]


def w_rad_shift():
    return _ORIG["rad_shift"]() + _S["d_eta_t"] + _S["d_eta_rough_t"]


def w_cc_losses(g, f, m_side, dbl, r4r2=None, dump=True):
    kc, kf = _S["k_cl_c"], _S["kf_c"]
    g2 = g if kc == 1.0 else dict(g, cl=g["cl"] * kc)
    L = _ORIG["cc_losses"](g2, f, m_side, dbl, r4r2, dump)
    if kf != 1.0:
        L = dict(L)
        for key in ("sf", "dfr", "df"):
            L[key] = L[key] * kf
    return L


def w_cc_stage(*a, **k):
    if "d_pen" in k and _S["k_pen"] != 1.0:
        k["d_pen"] = k["d_pen"] * _S["k_pen"]
    return _ORIG["cc_stage"](*a, **k)


def w_ax_turbine(T0in, p0in, b, m, U, n, w_tot, nst, phi, zn, zt, k_cl=None, shift=None):
    """зазор × k_cl_t: первый проход - расчётный зазор модели (по ступеням), второй - с увеличенным (средний по ступеням)"""
    kc = _S["k_cl_t"]
    if k_cl is None and kc != 1.0:
        r0 = _ORIG["ax_turbine"](T0in, p0in, b, m, U, n, w_tot, nst, phi, zn, zt, k_cl=None, shift=shift)
        k_cl = kc * sum(s["k_cl"] for s in r0["stages"]) / len(r0["stages"])
    return _ORIG["ax_turbine"](T0in, p0in, b, m, U, n, w_tot, nst, phi, zn, zt, k_cl=k_cl, shift=shift)


def w_rad_stage(T0in, p0in, b, m, U4, ER, n, zt, zn, k_cl=None, shift=None, eta0=0.85):
    kc = _S["k_cl_t"]
    if k_cl is None and kc != 1.0:
        r0 = _ORIG["rad_stage"](T0in, p0in, b, m, U4, ER, n, zt, zn, k_cl=None, shift=shift, eta0=eta0)
        k_cl = kc * r0["k_cl"]
    return _ORIG["rad_stage"](T0in, p0in, b, m, U4, ER, n, zt, zn, k_cl=k_cl, shift=shift, eta0=eta0)


def _hot(T):
    return T + _S["dT_gas"] if T > UF.T_HOT_C else T


def w_uts(mat_key, T):
    v, inr = _ORIG["_uts"](mat_key, _hot(T))
    return v * _S["k_uts"], inr


def w_ys(mat_key, T):
    v, inr = _ORIG["_ys"](mat_key, _hot(T))
    return v * _S["k_uts"], inr


def w_sigma_allow(T_C, hours=None, dT_cap=0.0):
    v, e = _ORIG["sigma_allow"](T_C + _S["dT_gas"], hours, dT_cap)
    return v * _S["k_lmp"], e


def w_rupture_hours(sigma_MPa, T_C, dT_cap=0.0):
    return _ORIG["rupture_hours"](sigma_MPa / _S["k_lmp"], T_C + _S["dT_gas"], dT_cap)


# =========================================================================== подстановка и восстановление
def _restore():
    ML.cc_shift, ML.eta_cap_l0, ML.ax_shift, ML.rad_shift = _ORIG["cc_shift"], _ORIG["eta_cap_l0"], _ORIG["ax_shift"], _ORIG["rad_shift"]
    ML.cc_losses, ML.cc_stage, ML.ax_turbine, ML.rad_stage = _ORIG["cc_losses"], _ORIG["cc_stage"], _ORIG["ax_turbine"], _ORIG["rad_stage"]
    ML.ko_Y.__defaults__ = _ORIG["ko_defaults"]
    ML.A.clear(); ML.A.update(_ORIG["A"])
    ML.K_DUMP = _ORIG["K_DUMP"]
    ML._SH.clear(); ML._SH.update(_ORIG["SH"])
    ST._uts, ST._ys, ST.sigma_allow, ST.rupture_hours = _ORIG["_uts"], _ORIG["_ys"], _ORIG["sigma_allow"], _ORIG["rupture_hours"]
    ST.K_BURST, ST.DT_COOL, ST.WEIBULL_M = _ORIG["K_BURST"], _ORIG["DT_COOL"], _ORIG["WEIBULL_M"]
    OD.TIP_INC_STALL, OD.DINC_STALL = _ORIG["TIP_INC_STALL"], _ORIG["DINC_STALL"]
    _S.clear(); _S.update(UF.nominal(), kf_c=1.0, kf_t=1.0, d_eta_rough_t=0.0)


class Perturb:
    """контекст одной реализации неопределённостей: with Perturb(vals): ... (vals - uq_factors.values / nominal)"""

    def __init__(self, vals):
        self.vals = dict(vals)

    def __enter__(self):
        _capture()
        if _ACTIVE[0]:
            raise RuntimeError("Perturb не вкладывается")
        _ACTIVE[0] = True
        v = self.vals
        try:
            _S.update(v, kf_c=1.0, kf_t=1.0, d_eta_rough_t=0.0)
            # 1) форма потерь -> привязка заново (номинальные коэффициенты дают ту же привязку - детерминированно)
            ML.A["zeta_d"], ML.A["f_inc"] = float(v["zeta_d"]), float(v["f_inc"])
            ML.K_DUMP = float(v["k_dump"])
            ML.ko_Y.__defaults__ = (_ORIG["ko_defaults"][0], float(v["tmax_c"]))
            if any(abs(float(v[k]) - float(UF.FACTORS[UF.INDEX[k]][4])) > 1e-15 for k in UF.SHAPE):
                ML._SH.clear()
                ML.calibration_summary()
            # 2) уровень, изготовление, материалы, температура, срыв - поверх привязки
            ML.ko_Y.__defaults__ = (float(v["t_te_mm"]) * 1e-3, float(v["tmax_c"]))
            ML.cc_shift, ML.eta_cap_l0, ML.ax_shift, ML.rad_shift = w_cc_shift, w_eta_cap_l0, w_ax_shift, w_rad_shift
            ML.cc_losses, ML.cc_stage, ML.ax_turbine, ML.rad_stage = w_cc_losses, w_cc_stage, w_ax_turbine, w_rad_stage
            ST._uts, ST._ys, ST.sigma_allow, ST.rupture_hours = w_uts, w_ys, w_sigma_allow, w_rupture_hours
            ST.K_BURST, ST.DT_COOL, ST.WEIBULL_M = float(v["k_burst"]), float(v["dt_cool"]), float(v["weibull_m"])
            OD.TIP_INC_STALL, OD.DINC_STALL = float(v["tip_inc_stall"]), float(v["dinc_stall"])
        except Exception:
            _restore()
            _ACTIVE[0] = False
            raise
        return self

    def __exit__(self, *exc):
        _restore()
        _ACTIVE[0] = False
        return False


def set_variant(c):
    """величины, зависящие от способа изготовления деталей варианта (внутри Perturb): трение компрессора - по колесу компрессора,
    турбин - по колесу турбины газогенератора (у СТ - тот же способ; оценка)"""
    kr = _S["k_rough"]
    _S["kf_c"] = UF.rough_factor(c.get("m_cw"), kr)
    _S["kf_t"] = UF.rough_factor(c.get("m_tw"), kr)
    _S["d_eta_rough_t"] = -UF.RFR_T * (_S["kf_t"] - 1.0)


def state():
    """снимок подменяемого состояния (для проверки восстановления)"""
    _capture()
    return dict(fn=tuple(id(x) for x in (ML.cc_shift, ML.eta_cap_l0, ML.ax_shift, ML.rad_shift, ML.cc_losses, ML.cc_stage, ML.ax_turbine,
                                         ML.rad_stage, ST._uts, ST._ys, ST.sigma_allow, ST.rupture_hours)),
                ko=ML.ko_Y.__defaults__, A=dict(ML.A), K_DUMP=ML.K_DUMP,
                SH={k: (v if not isinstance(v, dict) else None) for k, v in ML._SH.items() if k in ("cc", "ax", "rad")},
                ST=(ST.K_BURST, ST.DT_COOL, ST.WEIBULL_M), OD=(OD.TIP_INC_STALL, OD.DINC_STALL))
