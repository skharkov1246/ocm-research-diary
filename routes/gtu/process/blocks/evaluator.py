"""Единый вызов evaluate(concept): КПД, мощность, расход и температура выхлопа против допуска пакета, детали, готовность,
размеры для модели массы - по шести блокам «топливо -> вал».

Порядок (уровень облика; всё - оценка, модели не правятся - только вызываются):
  1. совместимость по описанию (compat.check_static);
  2. воронка: концепция воронки (concept.to_funnel) -> funnel_l0.solve (быстрый цикл + геометрия + поправки КПД по геометрии)
     -> КПД ступеней компрессора и турбин, число ступеней «auto»; у быстроходного вала - частота по φ первой ступени и повтор;
  3. цепочка блоков (chain.run) с этими КПД: промохлаждение, рекуператор, подогрев, рост давления, подогрев топлива, IBC;
     химическая рекуперация - приращение из перебора циклов (fuel.reform_delta);
  4. геометрия при расходе цепочки (funnel_space.derive с расходом и температурами цепочки), правила воронки (check),
     детали / обслуживание / Plug&Play / изготовимость (funnel_space.model_inputs + модели воронки), добавочные узлы
     теплообмена (cycles.parts), размеры сердцевин (cycles.size_hx), готовность (наименьший TRL);
  5. допуск пакета: воздух и выхлоп до +20 %, выхлоп до +100 °C (решение владельца 14.11.2026).
eta.source = «задано» - КПД узлов из концепции (калибровка, сверка с перебором циклов и gtu_cycle), геометрия - та же.
"""
import copy
import math
import time

from . import _track as T
from . import chain, combustor, compat, compressor, fuel, hx, shaft, turbine
from .concept import label, normalize, stage_types, to_funnel

RULE_GROUP_SECONDARY = ("выбросы",)      # не отсеивают: ни одна схема уровня 0 камеры не выполняет цели выбросов


# =========================================================================== воронка: КПД, число ступеней, геометрия
def _nominal_etas(fc):
    FS = T.FS()
    E = FS._E()
    typ = FS.stage_types(fc)
    pis = FS.stages_pi(fc)
    mats = FS._stage_mats(fc)
    mus = FS.stage_mus(fc)
    comp = [FS.eta_comp(t, p, m, k, fc["inlet"], mu) for k, (t, p, m, mu) in enumerate(zip(typ, pis, mats, mus))]
    if fc["turb"] == "О":
        gg = E.eta_axial(fc["lt"])
    else:
        gg = 0.86 * max(0.05, 1 - ((fc["lt"] - 0.70) / 0.70) ** 2)
    return dict(comp=comp, gg=gg, pt=FS.ETA_NOM["pt"])


def funnel_pass(C, tgt, pgain=0.0):
    """концепция воронки с числом ступеней «auto», её решение и КПД узлов. Возврат: fc, sol, etas, nt, npt, note"""
    FS, L0 = T.FS(), T.L0()
    t = C["turbine"]
    nt0 = int(t["nt"]) if t["nt"] != "auto" else (3 if t["type"] == "О" else 1)
    npt0 = int(t["npt"]) if t["npt"] != "auto" else 2
    fc = to_funnel(C, tgt, nt=nt0, npt=npt0, pgain=pgain)
    nt, npt = nt0, npt0
    if t["nt"] == "auto" or (C["rotors"] == 2 and t["npt"] == "auto"):
        p = FS.quick_perf(fc, tgt)
        if p["ok"]:
            tt = p["r"]["t"]
            Umax = FS.turb_umax(fc)[0]
            nt = turbine.nt_auto(C, tt["w"], Umax, dhs_fn=lambda k: FS._stage_dhs(tt, k))
            if p["r"]["pt"]:
                npt = turbine.npt_auto(C, p["r"]["pt"]["w"], fc["kUpt"] * FS.U_PT_MAX)
            fc = to_funnel(C, tgt, nt=nt, npt=npt, pgain=pgain)
    sol = L0.solve(fc, tgt)
    note = "воронка: funnel_l0.solve (быстрый цикл quick_perf, геометрия derive, поправки КПД corrections)"
    if sol["g"] is not None and sol["perf"]["ok"]:
        perf = sol["perf"]
        etas = dict(comp=[s["eta"] for s in perf["comp"]], gg=perf["eta_gg"], pt=perf["eta_pt"] or FS.ETA_NOM["pt"])
    else:
        etas = _nominal_etas(fc)
        note = f"воронка не построила решение ({sol.get('fail')}); КПД - номинальные формулы воронки (eta_comp, eta_axial / U/C0)"
    return fc, sol, etas, nt, npt, note


def chain_spec(C, etas, PRb_obj):
    xn, xr, xn2, xr2 = turbine.cool_fracs(C)
    cores = C["build"]["cores"]
    return dict(stages=compressor.stage_list(C, etas["comp"]), ic_after=compressor.ic_after(C), ic=compressor.ic_params(),
                xl=0.015 + (0.01 if cores > 1 else 0.0), xn=xn, xr=xr, xn2=xn2, xr2=xr2,
                rec=hx.recup_params(C) if C["hx"]["recup"] else None, PRb=PRb_obj, t41=float(C["turbine"]["t41"]),
                shaft="single" if C["rotors"] == 1 else "free", eta_gg=etas["gg"] - combustor.eta_t_penalty(C), eta_pt=etas.get("pt"),
                dp_duct=0.025 if cores > 1 else 0.01, reheat=hx.reheat_params(C) if C["hx"]["reheat"] else None,
                fuel=fuel.chain_spec(C), ibc=hx.ibc_params() if C["hx"]["ibc"] else None)


def _etas_given(C, etas_f):
    e = C["eta"]
    n = len(stage_types(C))
    comp = e["comp"] if isinstance(e["comp"], (list, tuple)) else [float(e["comp"])] * n
    return dict(comp=list(comp), gg=float(e["gg"]), pt=float(e.get("pt") or e["gg"]))


def _patched_perf(perf, ch, C, G, EGT, eff, w):
    """быстрый цикл воронки с расходом, температурами и КПД цепочки - для геометрии, правил и моделей воронки"""
    p = dict(perf)
    r = copy.deepcopy(perf["r"])
    last = ch["sections"][-1]
    if C["rotors"] == 1:
        r["t"].update(T_out=last["T_out"], p_out=last["p_out"], m=last["m"])
    elif r.get("pt"):
        r["pt"].update(T_out=last["T_out"], p_out=last["p_out"], m=last["m"], T_in=last["T_in"], p_in=last["p_in"])
    p.update(r=r, G=G, EGT=EGT, eta_th=eff, w=w / 1e3, T2=ch["T2"] - 273.15,
             T5=(ch["T_rec_in"] if ch["T_rec_in"] else ch["T_turb_exit"]) - 273.15)
    return p


# =========================================================================== Plug&Play: поправки блоков к модели воронки
def _pp_set(items, key, state, value, note):
    it = items.get(key)
    if it is None:
        return
    if value < it["value"]:
        items[key] = dict(state=state, value=value, weight=it["weight"], note=note + ("; " + it["note"] if it.get("note") else ""))


def plugplay_blocks(C, tgt, ppr):
    PP = T.PP()
    items = copy.deepcopy(ppr["items"])
    items = shaft.plugplay_fix(C, tgt, items)
    prep = C["fuel"]["prep"]
    if prep == "подогрев":
        _pp_set(items, "топливо", "переходник", 0.6, "подогреватель топлива на выхлопе, горячий газ к форсункам")
    elif prep == "паровой риформинг":
        _pp_set(items, "топливо", "переделка", 0.2, "риформер, вода: водоподготовка, насос, бак - в пакете нет")
    elif prep == "сухой риформинг":
        _pp_set(items, "топливо", "переделка", 0.2, "риформер и компрессор рециркуляции выхлопа")
    if C["compressor"]["intercool"] != "нет" or C["hx"]["ibc"]:
        _pp_set(items, "опоры и рама", "переделка", 0.2, "вынесенные охладители с вентиляторами")
        _pp_set(items, "укрытие", "переделка", 0.2, "место и вентиляция под охладители")
    if C["hx"]["ibc"]:
        _pp_set(items, "выхлоп", "переходник", 0.6, "IBC в выхлопном тракте (разрежение за турбиной IBC)")
    if C["hx"]["reheat"]:
        _pp_set(items, "топливо", "переходник", 0.6, "второй коллектор к межтурбинной камере")
        _pp_set(items, "САУ", "переделка", 0.2, "вторая камера: розжиг, деление топлива, защита")
    blocked = [k for k, v in items.items() if v["state"] == "невозможно"]
    score = 0.0 if blocked else 10.0 * sum(PP.WEIGHTS[k] * v["value"] for k, v in items.items())
    return dict(score=round(score, 3), blocked=blocked, no_data=ppr["no_data"], items=items)


# =========================================================================== главный вызов
def evaluate(concept, detail=False):
    t_start = time.time()
    FS, L0, PM, SM, MK, CY = T.FS(), T.L0(), T.PM(), T.SM(), T.MK(), T.CY()
    C = normalize(concept)
    tgt0 = T.package_target(C["kv"])
    rules = compat.check_static(C, tgt0)
    tgt = tgt0
    fc, sol, etas_f, nt, npt, eta_note = funnel_pass(C, tgt)
    n_hs = None
    if C["shaft"]["drive"] == "быстроходный" and sol["g"] is not None:
        typ = stage_types(C)
        k0 = next(k for k, x in enumerate(typ) if x != "О")
        n_hs, _ = shaft.n_highspeed(tgt0["n_pkg"], sol["g"]["stages"][k0]["phi"], typ[k0])
        tgt = T.package_target(C["kv"], n_out=C["shaft"].get("n_out") or n_hs)
        fc, sol, etas_f, nt, npt, eta_note = funnel_pass(C, tgt)
    given = C["eta"].get("source") == "задано"
    etas = _etas_given(C, etas_f) if given else etas_f
    if given:
        eta_note = "задано в концепции (калибровка / сверка); геометрия - воронка"
    PRb_obj = combustor.PRb(C, float(C["turbine"]["t41"]))
    ch = chain.run(chain_spec(C, etas, PRb_obj))
    if ch.get("ok") and callable(PRb_obj) and not given:          # рост давления: КПД воронки при том же p4/p3
        fc, sol, etas_f, nt, npt, eta_note = funnel_pass(C, tgt, pgain=ch["PRb"] - 1.0)
        etas = etas_f
        ch = chain.run(chain_spec(C, etas, PRb_obj))
    delta = fuel.reform_delta(C)
    out = dict(name=label(C), kv=C["kv"], package=dict(model=tgt0["model"], variant=tgt0["variant"], P_shaft_kW=tgt0["P_shaft_kW"],
                                                         P_kW=tgt0["P_kW"], G0=tgt0["G0"], W_exh0=tgt0["W_exh0"], EGT0=tgt0["EGT0"],
                                                         eta0_pct=tgt0["eta0"], n_pkg=tgt0["n_pkg_true"], pi0=tgt0["pi0"],
                                                         engine_kg0=tgt0["engine_kg0"]),
               eta_source=eta_note, nt=nt, npt=npt if C["rotors"] == 2 else 0)
    if not ch.get("ok"):
        rules.append(dict(code="chain", hard=True, ok=False, group="цикл", text=ch.get("why", "цепочка не замкнулась")))
        return _finish(out, C, rules, t_start, concept_n=C)
    if delta is not None and not delta.get("ok"):
        rules.append(dict(code="reform", hard=True, ok=False, group="цикл", text=delta.get("why")))
        return _finish(out, C, rules, t_start, concept_n=C)
    # ---------------------------------------------------------------- итог цикла (с приращением химической рекуперации)
    A = CY.A
    d_eff = delta["d_eff"] if delta else 0.0
    k_w = delta["k_w"] if delta else 1.0
    eff = ch["eff"] + d_eff
    w = ch["wnet"] * k_w
    EGT = ch["T_stack"] - 273.15 + (delta["d_EGT"] if delta else 0.0)
    m_exh = ch["m_exh"] + (delta["d_mexh"] if delta else 0.0)
    P = tgt0["P_shaft_kW"]
    G = P * 1e3 / w
    W_exh = G * m_exh
    gen = tgt0["variant"] == "генератор"
    k_el = A["eta_gen"] * A["eta_gb"] if gen else 1.0
    out.update(eff_shaft=round(eff * 100, 3), eff_el=round(eff * k_el * 100, 3), k_el=k_el, w_kJkg=round(w / 1e3, 3), G_air=round(G, 4),
               W_exh=round(W_exh, 4), EGT_C=round(EGT, 2), T41_C=C["turbine"]["t41"], PR=C["compressor"]["pi"],
               T3_C=round(ch["T3"] - 273.15, 2), T_rec_in_C=round(ch["T_rec_in"] - 273.15, 2) if ch["T_rec_in"] else None,
               T_turb_exit_C=round(ch["T_turb_exit"] - 273.15, 2), fuel_kg_s=round(G * ch["fuel"], 5),
               P_shaft_kW=P, P_el_kW=round(P * k_el, 1), water_t_h=round(G * delta["water_t_h_per_kgs"], 3) if delta else 0.0)
    # ---------------------------------------------------------------- допуск пакета (стык)
    soft = T.PP().SOFT
    G_cap = tgt0["G0"] * (1 + soft["alpha"])
    W_cap = tgt0["W_exh0"] * (1 + soft["alpha"]) if tgt0["W_exh0"] else None
    EGT_cap = tgt0["EGT0"] + soft["dT"]
    out["interface"] = dict(G_air=round(G, 3), G_cap=round(G_cap, 3), dG_pct=round((G / tgt0["G0"] - 1) * 100, 2),
                            W_exh=round(W_exh, 3), W_exh_cap=round(W_cap, 3) if W_cap else None,
                            dWexh_pct=round((W_exh / tgt0["W_exh0"] - 1) * 100, 2) if tgt0["W_exh0"] else None,
                            EGT_C=round(EGT, 1), EGT_cap_C=EGT_cap, dEGT=round(EGT - tgt0["EGT0"], 1), n_out_rpm=tgt0["n_pkg_true"],
                            rule="мягкий допуск: воздух и выхлоп до +20 %, выхлоп до +100 °C (решение владельца 14.11.2026)")
    rules.append(dict(code="pkg_Wexh", hard=True, ok=(W_cap is None) or W_exh <= W_cap + 1e-9, val=round(W_exh, 4),
                      lim=round(W_cap, 4) if W_cap else None, group="пакет",
                      text="мягкий допуск: расход выхлопа до +20 % паспорта (fleet, генераторное исполнение)"))
    if C["hx"]["recup"]:
        trm = hx.T_rec_max(C)
        rules.append(dict(code="recup_T", hard=True, ok=ch["T_rec_in"] - 273.15 <= trm + 1e-6, val=round(ch["T_rec_in"] - 273.15, 2),
                          lim=trm, group="теплообмен", text="газ на входе в рекуператор (LIM T_rec_max или задано)"))
    rules += combustor.rules(C, ch)
    # ---------------------------------------------------------------- геометрия и правила воронки при расходе цепочки
    g2, perf2 = None, None
    if sol["perf"] is not None and sol["perf"].get("ok") and sol["g"] is not None:
        perf2 = _patched_perf(sol["perf"], ch, C, G, EGT, eff, w)
        try:
            g2 = FS.derive(fc, tgt, perf2)
        except (ValueError, ZeroDivisionError, OverflowError) as e:
            rules.append(dict(code="geometry", hard=True, ok=False, group="воронка", text=f"геометрия не строится: {e}"))
    else:
        rules.append(dict(code="funnel_solve", hard=True, ok=False, group="воронка",
                          text=f"воронка не построила точку: {sol.get('fail')}"))
    if g2 is not None:
        for r in FS.check(fc, tgt, perf2, g2):
            if r["code"] in ("recup_T",):
                continue
            rules.append(dict(r, hard=True, group="воронка"))
        rules += shaft.rules(C, tgt0, g2["n"])
    # ---------------------------------------------------------------- теплообменники, IBC, добавочные узлы
    sz = hx.sizes(C, ch, G)
    ibc_w = hx.ibc_wheel(C, ch, G, g2["n"] if g2 else tgt0["n_pkg"]) if C["hx"]["ibc"] else None
    hard_fail = [r for r in rules if r.get("hard") and not r["ok"] and r.get("group") not in RULE_GROUP_SECONDARY]
    aux = hx.aux_parts(C, sz, ibc_w, limits_all=not any(r["code"] in ("soft_G", "soft_EGT", "pkg_Wexh") for r in hard_fail))
    # ---------------------------------------------------------------- детали, обслуживание, Plug&Play, изготовимость, TRL
    if g2 is not None:
        bom, sv, pp, mk = FS.model_inputs(fc, tgt, perf2, g2)
        if C["hx"]["recup"] and C["hx"]["recup_make"] == "фольга":
            mk = dict(mk, techs=[x for x in mk["techs"] if x != "печать_рекуператор"])
        pr, sr, ppr, mr = PM.count(bom), SM.score(sv), T.PP().score(pp), MK.score(mk)
        ppb = plugplay_blocks(C, tgt0, ppr)
        n_eff_tot = pr["n_eff"] + aux["added"]
        trl_items = [(x, MK.TRL[x]["trl"], MK.TRL[x]["basis"]) for x in mk["techs"]]
        trl_items += compressor.trl_items(C) + hx.trl_items(C) + shaft.trl_items(C)
        tf, tfb = fuel.trl(C["fuel"]["prep"])
        if C["fuel"]["prep"] != "нет":
            trl_items.append((f"топливо: {C['fuel']['prep']}", tf, tfb))
        tc, tcb = combustor.trl(C)
        if tc is not None:
            trl_items.append((f"камера: {C['combustor']['scheme']}", tc, tcb))
        trl_min = min(t for _, t, _ in trl_items)
        out["parts"] = dict(engine=dict(parts=pr["parts"], assemblies=pr["assemblies"], moving=pr["moving"], cooled=pr["cooled"],
                                        fasteners=pr["fasteners"], unique=pr["unique"], n_eff=round(float(pr["n_eff"]), 2)),
                            aux=aux, n_eff_total=round(float(n_eff_tot), 2),
                            s_parts_engine=pr["score"], s_parts_total=round(PM.score_from_neff(n_eff_tot), 3),
                            source="funnel/parts_model.count (двигатель) + organic/cycles.parts (добавочные узлы, по 1 эфф. детали)")
        out["service"] = dict(score=sr["score"], heaviest=sr["heaviest"], steps=sr["steps"], tools=sr["tools"], wrongable=sr["wrongable"])
        out["plugplay"] = dict(score=ppb["score"], blocked=ppb["blocked"], no_data=ppb["no_data"],
                               states={k: v["state"] for k, v in ppb["items"].items()}, pnp_cycles=aux["pnp"],
                               score_funnel_model=ppr["score"], water=bool(aux["water"]))
        out["make"] = dict(score=mr["score"], weeks=mr["weeks"], techs=mk["techs"])
        out["trl"] = dict(min=trl_min, k_trl=round(L0.trl_factor(trl_min), 4),
                          items=[dict(tech=x, trl=t, basis=b) for x, t, b in sorted(trl_items, key=lambda z: z[1])],
                          rule="готовность - наименьший TRL; множитель (0,5 + 0,5·TRL/9) - решение владельца 15.11.2026")
        env = L0.envelope(fc, tgt, g2)
        s_eta = L0.s_eta_of(eff * k_el * 100, tgt0["eta0"])
        s_size = L0.s_size_of(env["V"], env["V_bay"])
        s_em, _ = L0.effmake(s_eta, s_size, mr["score"])
        s = dict(parts=round(PM.score_from_neff(n_eff_tot), 3), plugplay=ppb["score"], service=sr["score"], effmake=round(s_em, 3))
        out["funnel_score"] = dict(score=round(L0.final_score(s, trl_min), 4), weighted=round(L0.weighted(s), 4), parts=s,
                                   note="оценка воронки (веса владельца 0,35 / 0,25 / 0,20 / 0,20 × множитель TRL) - для сравнения "
                                        "с воронкой; КПД на тонну - модель массы (отдельно)")
        out["sizing"] = sizing(C, fc, tgt, g2, env, sz)
    # ---------------------------------------------------------------- блоки
    out["blocks"] = dict(
        fuel=fuel.report(C, ch, delta),
        compressor=compressor.report(C, ch, g2, perf2),
        combustor=combustor.report(C, ch, PRb_obj),
        turbine=turbine.report(C, ch, g2, etas, nt, npt if C["rotors"] == 2 else 0),
        hx=hx.report(C, ch, sz, ibc_w),
        shaft=shaft.report(C, tgt, g2, ch, G, P),
    )
    if detail:
        out["chain"] = _jsonable(ch)
        out["funnel_concept"] = {k: v for k, v in fc.items()}
    return _finish(out, C, rules, t_start, concept_n=C)


def _finish(out, C, rules, t_start, concept_n):
    FS = T.FS()
    bad = [r for r in rules if r.get("hard") and not r["ok"] and r.get("group") not in RULE_GROUP_SECONDARY]
    soft_bad = [r for r in rules if not r["ok"] and (not r.get("hard") or r.get("group") in RULE_GROUP_SECONDARY)]
    out["feasible"] = not bad
    out["reject"] = [r["code"] for r in bad]
    out["warnings"] = [r["code"] for r in soft_bad]
    out["violation"] = FS.violation([dict(r) for r in bad])
    out["emissions_ok"] = all(r["ok"] for r in rules if r.get("group") == "выбросы")
    out["rules"] = [{k: v for k, v in r.items()} for r in rules]
    out["concept"] = concept_n
    out["runtime_s"] = round(time.time() - t_start, 4)
    out["status"] = "оценка уровня облика: КПД узлов, пределы, потери, детали и TRL - принятые величины трека и оценки"
    return out


def sizing(C, fc, tgt, g, env, sz):
    """размеры для модели массы: двигатель (воронка derive / envelope), сердцевины теплообменников (cycles.size_hx), материалы.
    Массы модулей - формулы воронки (оболочки и колёса по подобию gtu_explore.finalize) - не калиброваны по массам машин."""
    MK = T.MK()
    st = g["stages"]
    rho = lambda m: MK.MAT.get(m, {}).get("rho")
    tech = hx.TECH[C["hx"]["recup_make"]] if C["hx"]["recup"] else "lpbf"
    hx_sel = {k: v.get(tech) for k, v in sz.items()}
    return dict(
        engine=dict(n_rpm=round(g["n"]), cores=fc["cores"], Dc=round(g["Dc"], 4), D_case=round(g["D_case"], 4),
                    L_case=round(g["L_case"], 4), L_rotor=round(g["L_rot"], 4), D_hot=round(g["D_hot"], 4), L_hot=round(g["L_hot"], 4),
                    D_comb=round(g["D_comb"], 4), L_comb=round(g["L_comb"], 4), n_inj=g["n_inj"], D_env=round(g["D_env"], 4),
                    compressor=[dict(type=s["type"], D2=round(s["D"], 4), b2=round(s["b2"], 5) if "b2" in s else None, U2=round(s["U"], 1),
                                     mat=s["mat"], rho=rho(s["mat"])) for s in st],
                    turbine=dict(type=fc["turb"], nt=fc["nt"], D=round(g["tg"]["D"], 4), Dtip=round(g["tg"]["Dtip"], 4),
                                 h_last=round(g["tg"]["h"], 4) if g["tg"].get("h") else None, mat=fc["mat_tw"], rho=rho(fc["mat_tw"])),
                    power_turbine=(dict(npt=fc["npt"], D=round(g["pt"]["D"], 4), Dtip=round(g["pt"]["Dtip"], 4), h_last=round(g["pt"]["h"], 4))
                                   if g["pt"] else None),
                    hot=dict(make=fc["m_hot"], mat=fc["mat_hot"], rho=rho(fc["mat_hot"])), case=dict(make=fc["m_case"], arch=fc.get("case_arch")),
                    envelope=dict(L=round(env["L"], 3), W=round(env["W"], 3), H=round(env["H"], 3), V=round(env["V"], 3))),
        hx=hx_sel, hx_tech=tech,
        module_mass_funnel_kg={k: round(v, 1) for k, v in g["mass"].items()},
        mass_funnel_estimate_kg=round(env["mass"]),
        hx_core_mass_kg=round(sum((v or {}).get("mass_kg", 0.0) for v in hx_sel.values()), 1),
        note="массы - формулы воронки (funnel_space.derive, funnel_l0.envelope MASS_K 1,25) и сердцевины cycles.size_hx: не "
             "калиброваны по массам машин - вход модели массы, не итог")


def _jsonable(x):
    if isinstance(x, dict):
        return {k: _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, float):
        return round(x, 6) if math.isfinite(x) else None
    if hasattr(x, "item"):
        return _jsonable(x.item())
    return x
