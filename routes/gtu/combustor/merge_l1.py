"""Воронка камеры сгорания KV, уровень 1: слияние частей, сверка сетью реакторов Cantera, цели, отбор на уровень 2.

Проверки полноты (как merge_l0): у всех частей одна схема, типоразмер, объём n, число частей, отпечаток входов и кода,
итерации и сетка; номера частей 0..S-1 по одному разу; позиции вариантов = range(n) без повторов. Неполное слияние -
complete=false и код выхода 2.

На вариант (README, «Уровень 1»):
  * NOx и CO: на 100 % - выход трубы (после разбавления); на 50 и 30 % - выход ячейки (до разбавления: CO - оценка
    сверху, разбавление CO не добавляет); полнота сгорания на каждой нагрузке;
  * второй способ (навык openfoam-combustor, раздел 7): реактор Cantera GRI-Mech 3.0 постоянного давления от
    равновесного состава горения при phi среза канала (NO обнулён) на время пребывания горячей зоны CFD; отношение
    NOx CFD / Cantera вне 0,5-2 - флаг «разбор» (число не публикуется);
  * потери: диффузор и переходный канал (уровень 0) + max(перепад смесителя уровня 0, потери ячейки CFD) + потери
    трубы CFD + решётка стенки; RTDF выхода трубы (осесимметрия: PF = RTDF - оценка снизу настоящего PF);
  * стенка: температура металла - одномерный баланс в форме comb_crn.wall_temperature, где температура адиабатной
    стенки - газ пристенных ячеек трубы CFD (плёнка разрешена), а не корреляция плёнки;
  * сетка: разница coarse/medium по ячейке 100 % у каждого k-го варианта;
  * цели targets.json, пересчитанные L1: nox_gas, co_gas, load_range, exit_uniformity (по RTDF), pressure_loss,
    wall_life; остальные - статус уровня 0. Рейтинг - веса merge_l0, КПД и NOx - по L1, штраф жёстких целей.
Отбор на уровень 2: 30 вариантов (--k2) - по кругу ниш по рейтингу L1.

Запуск: python merge_l1.py runs/<run_id> --inputs l1_inputs_KV-3.json --out comb_l1_KV-3_results.json
"""
import argparse, glob, json, math, os, sys, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

DATE = "2026-11-14"
CANTERA_BAND = (0.5, 2.0)


def load_parts(paths):
    files = []
    for p in paths:
        files += glob.glob(os.path.join(p, "**", "comb_l1_*.json"), recursive=True) if os.path.isdir(p) else [p]
    files = sorted(set(f for f in files if not f.endswith(".tmp") and "results" not in os.path.basename(f)))
    return [(f, json.load(open(f, encoding="utf-8"))) for f in files]


def check_complete(parts):
    rep = dict(files=len(parts), problems=[])
    if not parts:
        rep["problems"].append("нет частей")
        return False, rep, {}
    keys = ("schema", "kv", "n", "shards", "inputs_fingerprint", "code_hash", "level", "iters", "generate_only")
    ref = {k: parts[0][1].get(k) for k in keys}
    rep["ref"] = ref
    for f, p in parts:
        for k in keys:
            if p.get(k) != ref[k]:
                rep["problems"].append(f"{os.path.basename(f)}: {k}={p.get(k)!r} != {ref[k]!r}")
    n, S = ref["n"], ref["shards"]
    seen = {}
    for f, p in parts:
        seen.setdefault(p["shard"], []).append(f)
    rep["missing_shards"] = sorted(set(range(S)) - set(seen))
    rep["duplicate_shards"] = {k: v for k, v in seen.items() if len(v) > 1}
    rows, used, dup = {}, set(), 0
    for f, p in parts:
        if p["shard"] in used:
            continue
        used.add(p["shard"])
        expect = set(range(p["shard"], n, S))
        got = [r["pos"] for r in p["rows"]]
        if set(got) != expect or len(got) != len(expect):
            rep["problems"].append(f"{os.path.basename(f)}: позиции не совпадают с разбиением i % {S} == {p['shard']}")
        for r in p["rows"]:
            dup += r["pos"] in rows
            rows[r["pos"]] = r
    rep["duplicate_positions"] = dup
    missing = sorted(set(range(n)) - set(rows))
    rep["missing_positions"] = len(missing)
    complete = not rep["problems"] and not rep["missing_shards"] and not missing and dup == 0
    return complete, rep, rows


def _case(row, name):
    c = row["cases"].get(name)
    return c["result"] if c and c.get("run_status") == "ok" and "result" in c else None


def cantera_check(e, cell, li=0):
    """второй способ: реактор постоянного давления (GRI-Mech 3.0) от равновесия HP при phi среза канала с NO = 0 на
    время пребывания горячей зоны CFD; NOx и CO в ppm сухого при 15 % O2"""
    import cantera as ct
    import comb_crn as CR
    L = e["loads"][li]
    phi = cell.get("phi_exit_mean")
    m_w = cell.get("mass_out")
    if not phi or not m_w or not cell.get("mass_hot"):
        return None
    tau = cell["mass_hot"] / m_w
    g = ct.Solution("gri30.yaml")
    g.set_equivalence_ratio(phi, "CH4:1", CR.AIR)
    g.TP = L["T3_K"], L["p3_Pa"]
    g.equilibrate("HP")
    Y = g.Y.copy()
    for s in ("NO", "NO2", "N2O"):
        Y[g.species_index(s)] = 0.0
    g.TPY = g.T, g.P, Y / Y.sum()
    r = ct.IdealGasConstPressureReactor(g, clone=False)
    net = ct.ReactorNet([r])
    net.advance(tau)
    X = {s: float(r.phase.X[r.phase.species_index(s)]) for s in ("NO", "NO2", "CO", "O2", "H2O")}
    dry = 1 - X["H2O"]
    o2 = 100 * X["O2"] / dry
    corr = (20.9 - 15.0) / max(20.9 - o2, 1e-6)
    nox = 1e6 * (X["NO"] + X["NO2"]) / dry * corr
    return dict(phi=phi, tau_s=tau, T=float(r.T), NOx_ppm15=nox, CO_ppm15=1e6 * X["CO"] / dry * corr,
                method="реактор постоянного давления GRI-Mech 3.0 от равновесия HP при phi среза, NO=0, время - "
                       "масса горячей зоны CFD / расход")


def metal_temperature(e, T_aw, T_gas):
    """одномерный баланс стенки в форме comb_crn.wall_temperature, температура адиабатной стенки - из CFD (T_aw),
    газ для излучения - T_gas (средняя по горячей зоне); константы - comb_crn (оценка, по памяти / по выдержке поиска)"""
    import comb_crn as CR
    import comb_space as CS
    from scipy.optimize import brentq
    c, d = e["concept"], e["derived"]
    bc = CR.boundary()
    L0 = bc["loads"][0]
    p, T3, W = L0["p3_Pa"], L0["T3_K"], L0["W_air"]
    N = CS.N_CAN[c["chamber"]]
    mat = "керамика" if c["cooling"] == "керамика" else "металл"
    ew, kw, tw = CR.EPS_W[mat], CR.K_WALL[mat], CS.WALL_T[mat]
    Dh = d["D_h"]
    far = min(L0["W_fuel"] / (c["a_pz"] * W), bc["far_st_CH4"])
    eg = 1 - math.exp(-290 * (p / 1e3) * (far * 0.85 * Dh) ** 0.5 * T_gas ** -1.5)
    K = CR.kin()
    K.GM.TPX = T_gas, p, "CO2:0.06,H2O:0.12,O2:0.09,N2:0.73"
    kg, mug = K.GM.thermal_conductivity, K.GM.viscosity
    mg = c["a_pz"] * W / N
    h1 = 0.020 * kg / Dh ** 0.2 * (mg / (d["A_L"] * mug)) ** 0.8
    ka, mua = L0["k3"], L0["mu3"]
    m_an = W * max(0.2, 1 - c["a_pz"]) / N
    A_an = m_an / (L0["rho3"] * CS.U_ANN)
    D_an = 2 * A_an / (math.pi * (Dh if c["chamber"] in ("выносная", "трубчатые") else 2 * d["D_m"]))
    h2 = 0.020 * ka / D_an ** 0.2 * (m_an / (A_an * mua)) ** 0.8
    Z = ew * CR.EPS_C / (CR.EPS_C + ew * (1 - CR.EPS_C) * min(1.0, Dh / max(d["D_case"], Dh)))

    def outer(t2):
        return Z * CR.SIGMA * (t2 ** 4 - T3 ** 4) + h2 * (t2 - T3)

    def resid(t1):
        t2 = brentq(lambda t: kw / tw * (t1 - t) - outer(t), T3, t1)
        return 0.5 * CR.SIGMA * (1 + ew) * eg * T_gas ** 1.5 * (T_gas ** 2.5 - t1 ** 2.5) + h1 * (T_aw - t1) \
            - kw / tw * (t1 - t2)
    hi = max(T_aw, T_gas) - 1e-3
    t1 = hi if resid(hi) > 0 else brentq(resid, T3 + 1e-3, hi, xtol=0.05)
    return dict(T_metal_K=t1, T_lim_K=CR.T_WALL_MAX[mat], eps_g=eg, h1=h1, h2=h2, T_aw_cfd=T_aw, T_gas=T_gas, mat=mat)


def summarize(row, e, tg):
    """итоги варианта уровня 1 и соответствие целям"""
    c100, c50, c30 = (_case(row, f"cell_{p}_coarse") for p in (100, 50, 30))
    lin = _case(row, "liner_100_coarse")
    med = _case(row, "cell_100_medium")
    s = dict(pos=row["pos"], idx=row["idx"], why=row["why"], scheme=e["concept"]["scheme"],
             chamber=e["concept"]["chamber"], injector=e["concept"]["injector"], passive=e["concept"]["passive"],
             status=row["status"])
    out = lambda r: (r or {}).get("outlet", {}) if r else {}
    for tag, r in (("100", c100), ("50", c50), ("30", c30)):
        o = out(r)
        s[f"cell_{tag}"] = None if r is None else dict(
            status=r["status"], NOx_ppm15=o.get("NOx_ppm15"), CO_ppm15=o.get("CO_ppm15"), eta_c=o.get("eta_c"),
            T_out=o.get("T"), S_exit=r.get("S_exit"), phi_exit=r.get("phi_exit_mean"), recirc=r.get("recirc_ratio_2cal"),
            flashback=r.get("flashback_or_autoignition"), T_max=r.get("T_max"), dp0=r.get("dp0_rel"),
            converged=(r.get("convergence") or {}).get("ok"), mass_imbalance=r.get("mass_imbalance"))
    lo = out(lin)
    s["liner_100"] = None if lin is None else dict(
        status=lin["status"], NOx_ppm15=lo.get("NOx_ppm15"), CO_ppm15=lo.get("CO_ppm15"), eta_c=lo.get("eta_c"),
        T_out=lo.get("T"), RTDF=lo.get("RTDF"), T_aw_max=lin.get("T_aw_max"), dp0=lin.get("dp0_rel"),
        converged=(lin.get("convergence") or {}).get("ok"))
    if med and c100:
        a, b = out(c100), out(med)
        s["grid"] = dict(dNOx_rel=(a.get("NOx_ppm15", 0) - b.get("NOx_ppm15", 0)) / max(b.get("NOx_ppm15", 0), 0.5),
                         dCO_rel=(a.get("CO_ppm15", 0) - b.get("CO_ppm15", 0)) / max(b.get("CO_ppm15", 0), 1.0),
                         dT_rel=(a.get("T", 0) - b.get("T", 0)) / max(b.get("T", 1), 1.0),
                         dS=(c100.get("S_exit") or 0) - (med.get("S_exit") or 0))
    # второй способ и стенка
    s["cantera"] = cantera_check(e, c100) if c100 else None
    if s["cantera"] and c100 and out(c100).get("NOx_ppm15") is not None:
        ratio = out(c100)["NOx_ppm15"] / max(s["cantera"]["NOx_ppm15"], 0.05)
        s["cantera"]["ratio_cfd"] = ratio
        s["cantera"]["flag"] = "разбор" if not CANTERA_BAND[0] <= ratio <= CANTERA_BAND[1] else "ok"
    d = e["derived"]
    dp_cell = max(e["concept"]["dp_inj"], (c100 or {}).get("dp0_rel") or 0.0)
    if lin:
        s["dp_total"] = d["dp_diff"] + d["dp_duct"] + d["dp_series_wall"] + dp_cell + (lin.get("dp0_rel") or 0.0)
        Tg = e["loads"][0]["hot_init"]["T_eq"]
        s["wall"] = metal_temperature(e, lin["T_aw_max"], Tg) if lin.get("T_aw_max") else None
    # цели (пересчитанные уровнем 1)
    t = {x["id"]: x for x in tg["targets"]}
    lim_nox, lim_co = t["nox_gas"]["target"]["max"], t["co_gas"]["target"]["max"]
    nox_lo, co_lo = t["nox_gas"]["also"]["max"], t["co_gas"]["also"]["max"]
    burnt = all(r is not None and (out(r).get("eta_c") or 0) >= 0.99 for r in (c100, c50))
    burnt30 = c30 is not None and (out(c30).get("eta_c") or 0) >= 0.99
    T = {}

    def put(tid, val, lim):
        T[tid] = dict(value=val, limit=lim, ok=bool(val is not None and val <= lim),
                      viol=1.0 if val is None else max(0.0, (val - lim) / abs(lim)))
    put("nox_gas", max(lo.get("NOx_ppm15", 1e9), out(c50).get("NOx_ppm15", 1e9)) if burnt and lin else None, lim_nox)
    put("co_gas", max(lo.get("CO_ppm15", 1e9), out(c50).get("CO_ppm15", 1e9)) if burnt and lin else None, lim_co)
    put("load_range", max(out(c30).get("NOx_ppm15", 1e9) / nox_lo, out(c30).get("CO_ppm15", 1e9) / co_lo)
        if burnt and burnt30 else None, 1.0)
    put("exit_uniformity", lo.get("RTDF") / t["exit_uniformity"]["target"]["RTDF_max"] if lo.get("RTDF") is not None
        else None, 1.0)
    put("pressure_loss", s.get("dp_total"), t["pressure_loss"]["target"]["max"])
    put("wall_life", (s.get("wall") or {}).get("T_metal_K", 1e9) / t["wall_life"]["target"]["T_metal_max_K"]
        if s.get("wall") else None, 1.0)
    s["targets_l1"] = T
    l0fail = set(e["l0"].get("hard_fail", []))
    hard = tg["funnel_interface"]["hard"]
    s["hard_fail"] = sorted([h for h in hard if (h in T and not T[h]["ok"]) or (h not in T and h in l0fail)])
    s["hard_ok"] = not s["hard_fail"]
    return s


def score_l1(s, e):
    """рейтинг: веса merge_l0, КПД (потери) и NOx - по уровню 1; штраф exp(-2·Σ нарушений), нецелевые - по L0 (1/0)"""
    import merge_l0 as M0
    l0 = e["l0"]
    dp = s.get("dp_total") if s.get("dp_total") is not None else l0.get("dp_total", 0.08)
    nox = (s.get("liner_100") or {}).get("NOx_ppm15")
    nox = nox if nox is not None else 30.0
    V_env = math.pi / 4 * l0["D_env"] ** 2 * l0["L_mod"]
    sc = dict(parts=M0._clip(1 - math.log(max(l0["parts"], 3) / 3) / math.log(40 / 3)), plugplay=l0["plugplay"] / 10,
              service=l0["service"] / 10, eff=M0._clip((0.08 - dp) / 0.06),
              size=M0._clip(1 - math.log(max(V_env, 0.02) / 0.02) / math.log(3.0 / 0.02)),
              make=l0["make_score"] / 10, nox=M0._clip(1 - nox / 30.0))
    viol = sum(min(v["viol"], 3.0) for k, v in s["targets_l1"].items())
    viol += sum(1.0 for h in s["hard_fail"] if h not in s["targets_l1"])
    return M0.score(sc, math.exp(-M0.PENALTY_K * viol), M0.WEIGHTS)


def merge(paths, inputs, out, k2=30):
    t0 = time.time()
    import comb_crn as CR
    parts = load_parts(paths)
    complete, rep, rows = check_complete(parts)
    data = json.load(open(inputs, encoding="utf-8"))
    if rep.get("ref", {}).get("inputs_fingerprint") not in (None, data["fingerprint"]):
        rep["problems"].append("отпечаток входов частей не совпадает с файлом входов")
        complete = False
    ents = data["entries"]
    summ = []
    for pos in sorted(rows):
        s = summarize(rows[pos], ents[pos], CR.TARGETS)
        s["score_l1"] = score_l1(s, ents[pos])
        summ.append(s)
    summ.sort(key=lambda x: -x["score_l1"])
    by = {}
    for s in summ:
        by.setdefault(f'{s["scheme"]}|{s["chamber"]}|{s["injector"]}|{s["passive"]}', []).append(s)
    order = sorted(by, key=lambda n: -by[n][0]["score_l1"])
    sel = []
    while len(sel) < k2 and any(by[n] for n in order):
        for n in order:
            if by[n] and len(sel) < k2:
                x = by[n].pop(0)
                if x["status"] != "ошибка" and (x.get("cantera") or {}).get("flag") != "разбор":
                    sel.append(dict(pos=x["pos"], idx=x["idx"], score_l1=x["score_l1"], niche=n))
    feed = {}
    for s in summ:
        key = f'{s["scheme"]}|{s["chamber"]}'
        if key in feed or s.get("dp_total") is None:
            continue
        l0 = ents[s["pos"]]["l0"]
        feed[key] = dict(idx=s["idx"], level="L1", dp_comb=s["dp_total"], NOx_100=(s["liner_100"] or {}).get("NOx_ppm15"),
                         CO_100=(s["liner_100"] or {}).get("CO_ppm15"), NOx_50=(s["cell_50"] or {}).get("NOx_ppm15"),
                         NOx_30=(s["cell_30"] or {}).get("NOx_ppm15"), RTDF=(s["liner_100"] or {}).get("RTDF"),
                         T_metal_max_K=(s.get("wall") or {}).get("T_metal_K"), D_comb=l0["D_env"], L_comb=l0["L_mod"],
                         m_comb=l0["mass_kg"], parts=l0["parts"], n_inj=l0["n_el"], hard_ok=s["hard_ok"])
    flags = sum(1 for s in summ if (s.get("cantera") or {}).get("flag") == "разбор")
    res = dict(title=f"Камера сгорания {data['kv']}, уровень 1: осесимметричный реагирующий RANS (OpenFOAM)",
               method=("слияние частей comb_l1.py с контролем полноты; итоги кейсов comb_post; второй способ - реактор "
                       "Cantera; стенка - баланс comb_crn с T адиабатной стенки из CFD; рейтинг - веса merge_l0 со "
                       "штрафом целей L1; отбор на L2 по кругу ниш"),
               software=f"python {sys.version.split()[0]}, numpy {np.__version__}, cantera {CR.ct.__version__}",
               date=DATE, status="RANS осесимметричный (уровень 1): числа - оценка для отбора на уровень 2",
               complete=complete, completeness=rep,
               inputs=[dict(file=os.path.basename(inputs), fingerprint=data["fingerprint"])] +
                      [dict(file=os.path.basename(f), shard=p["shard"], runtime_s=p["runtime_s"]) for f, p in parts],
               counts=dict(n=rep.get("ref", {}).get("n"), rows=len(rows), hard_ok=sum(s["hard_ok"] for s in summ),
                           cantera_flags=flags,
                           status={k: sum(1 for s in summ if s["status"] == k) for k in ("ok", "частично", "ошибка")}),
               rows=summ, select_l2=sel,
               feed_to_engine_funnel=dict(note="лучший по рейтингу L1 в каждой паре (схема, тип камеры); dp, NOx, CO, RTDF, "
                                               "стенка - L1, габарит, масса, детали - L0", best=feed),
               runtime_s=round(time.time() - t0, 1))
    json.dump(res, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=float)
    return complete, res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("paths", nargs="+")
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--k2", type=int, default=30)
    a = ap.parse_args(argv)
    complete, res = merge(a.paths, a.inputs, a.out, a.k2)
    c = res["counts"]
    print(f"полнота: {complete}; строк {c['rows']} из {c['n']}; все жёсткие цели - {c['hard_ok']}; флагов сверки "
          f"{c['cantera_flags']}; запись {a.out}")
    return 0 if complete else 2


if __name__ == "__main__":
    sys.exit(main())
