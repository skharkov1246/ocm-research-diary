"""Справочник библиотеки: блоки, варианты, параметры быстрых моделей, готовность и источники - для results-JSON и README.
Числа берутся из модулей трека в момент вызова (ничего не переписано руками)."""
import inspect

from . import _track as T
from . import combustor, compat, compressor, fuel, hx, shaft, turbine
from .concept import DEFAULT, PSI_AUTO, NU_AUTO


def catalog():
    CY, FS, MK = T.CY(), T.FS(), T.MK()
    A = CY.A
    sig = inspect.signature(T.GC().recup).parameters
    l0 = combustor.l0()
    comb_rows = {}
    for s in combustor.SCHEME:
        t, tb = combustor.trl(dict(combustor=dict(scheme=s)))
        sch = l0["by_scheme"].get(s) or {}
        pair = l0["feed_to_engine_funnel"]["best"].get(f"{s}|кольцевая") or {}
        comb_rows[s] = dict(name=combustor.SCHEME[s]["name"], trl=t, trl_basis=tb, premix=combustor.SCHEME[s]["premix"],
                            pressure_gain=combustor.SCHEME[s]["pgc"],
                            dp_best_annular=pair.get("dp_comb"), dp_min_scheme=sch.get("dp_min"), parts_best_annular=pair.get("parts"),
                            NOx_100_best_annular=pair.get("NOx_100"), CO_100_best_annular=pair.get("CO_100"),
                            l0_variants=sch.get("total"), l0_feasible=sch.get("ok"), l0_hard_ok=sch.get("hard_ok"))
    return dict(
        blocks=dict(
            fuel=dict(title="1. Подготовка топлива", options={k: dict(name=v["name"], trl=fuel.trl(k)[0], trl_basis=fuel.trl(k)[1],
                                                                       water=v["water"], note=v["note"]) for k, v in fuel.OPTIONS.items()},
                      params=dict(T_fh_max_C=A["T_fh_max"], dT_app_K=A["dT_app"], dp_fh_gas=A["dp_fh_g"], sc=list(A["sc"]), R_egr=list(A["R_egr"]),
                                  source="organic/cycles A"),
                      model="подогрев - в цепочке; риформинг - приращение cycles.evaluate (Cantera) к цепочке концепции"),
            compressor=dict(title="2. Компрессор", options=compressor.ARCH, pi_range={str(k): v for k, v in FS.PI_RANGE.items()},
                            materials={m: dict(U_forge=MK.u_limit_comp(m, "поковка")[0], U_cast_print=MK.u_limit_comp(m, "литьё")[0])
                                       for m in ("Ti64", "Al", "Ni")},
                            intercooler=dict(compressor.ic_params(), trl=CY.TRL["ic"], trl_basis=CY.TRL_NOTE["ic"], source="organic/cycles A"),
                            model="КПД ступеней - funnel_space.eta_comp + funnel_l0.corrections (зазор, быстроходность по геометрии derive)"),
            combustor=dict(title="3. Камера сгорания", options=comb_rows, targets=combustor.targets(), dp_track=combustor.DP_TRACK,
                           T3_premix_max_C=combustor.T3_PREMIX_MAX, T3_cat_min_K=combustor.T3_CAT_MIN,
                           pgc=dict(f_real=A["pgc_f_real"], deta_t=A["pgc_deta_t"], dp_inlet=A["dp_comb"], model="cycles.pgc_gain"),
                           owner_decision="один клапан на каждое топливо без приводов; частичную нагрузку держит пассивный элемент "
                                          "(README трека, решения 15.11.2026)",
                           source="combustor/l0_aws (262 144 варианта KV-3), combustor/targets.json, comb_space.TRL_SCHEME"),
            turbine=dict(title="4. Турбины", cooling={k: dict(name=v, x_ngv=MK.COOL_AIR[k][0], x_rotor=MK.COOL_AIR[k][1],
                                                               t41_ngv=MK.NGV_T41[k][0], basis=MK.NGV_T41[k][1])
                                                       for k, v in turbine.COOL.items()},
                         materials=dict(blades={k: list(v) for k, v in FS.MAT_TW.items()}, hot={k: list(v) for k, v in FS.MAT_HOT.items()}),
                         nt_auto=dict(psi_axial=PSI_AUTO, nu_radial=NU_AUTO),
                         model="КПД - funnel_space.quick_perf (eta_axial по ψ, U/C0 радиальной); пределы - make_model"),
            hx=dict(title="5. Теплообменники", recup=dict(eps_range=list(hx.RECUP_EPS), dp_air=sig["dp_air"].default, dp_gas=sig["dp_gas"].default,
                                                          T_rec_max_C=FS.T_REC_MAX, techs={k: CY.HX_TECH[v]["name"] for k, v in hx.TECH.items()},
                                                          trl_foil=CY.TRL["rec"], trl_print=MK.TRL["печать_рекуператор"]["trl"]),
                    reheat=dict(dp=A["dp_rh"], x_rh=list(A["x_rh"]), trl=CY.TRL["rh"], trl_basis=CY.TRL_NOTE["rh"]),
                    ibc=dict(hx.ibc_params(), trl=CY.TRL["ibc"], trl_basis=CY.TRL_NOTE["ibc"]),
                    model="цепочка + размеры cycles.size_hx (ламинарные каналы) + детали cycles.parts"),
            shaft=dict(title="6. Вал и стык пакета", options=shaft.DRIVES, hs_ratio=list(shaft.HS_RATIO), hs_trl=shaft.HS_TRL,
                       soft_tolerance=dict(T.PP().SOFT), note="мягкий допуск пакета: воздух и выхлоп до +20 %, выхлоп до +100 °C"),
        ),
        compatibility=compat.TABLE,
        defaults=DEFAULT,
    )
