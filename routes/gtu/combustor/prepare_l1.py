"""Воронка камеры сгорания KV, уровень 1: отбор вариантов уровня 0 и входной файл кейсов (l1_inputs_*.json).

Идёт в сессии после слияния уровня 0 (нужны Cantera и comb_bc_results.json того же кода). Отбор - K вариантов
(по умолчанию 500) из допустимых (status = ok) по рейтингу уровня 0 (merge_l0: взвешенная оценка × штраф жёстких
целей) с разнообразием: по кругу ниш «схема × тип камеры × форсунка × пассивное деление», внутри ниши - по рейтингу;
сначала члены фронта Парето. При N = 300 ни один вариант не выполнил все жёсткие цели, поэтому отбор идёт по
рейтингу со штрафом, а не по допуску (README, «Уровень 1»).

Для каждого варианта: оси (comb_space.concept по индексу), геометрия и гидравлика (comb_space.derive), по каждой
нагрузке 100/50/30 % - p3, T3, расходы, доля топлива пилота (comb_crn.pilot_share), доля воздуха, встречающая
топливо при самостратификации (comb_crn.strat_air_frac), состав и температура горячего поля для зажигания
(равновесие Cantera GRI-Mech 3.0 при phi зоны, температура не ниже 1600 K - только начальное поле).

Запуск:
  python prepare_l1.py runs/<run_id> --k 500 --out l1_inputs_KV-3.json         # из частей уровня 0
  python prepare_l1.py comb_l0_debug_results.json --k 20 --out /tmp/l1.json    # из слитого файла (top + pareto)
"""
import argparse, hashlib, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import comb_space as CS  # noqa: E402
import comb_crn as CR    # noqa: E402
import merge_l0 as M0    # noqa: E402

SCHEMA = "kv.combustor.l1.inputs/1"
T_IGN_MIN = 1600.0       # температура горячего поля зажигания не ниже, K (начальное поле, не результат)


def niche(c):
    return f'{c["scheme"]}|{c["chamber"]}|{c["injector"]}|{c["passive"]}'


def _rows_from(paths):
    """полные строки уровня 0 из частей (папки/файлы comb_l0_*.json) или компактные из слитого файла"""
    full, compact, meta = [], [], {}
    for p in paths:
        if os.path.isfile(p) and "top" in json.load(open(p, encoding="utf-8")):
            m = json.load(open(p, encoding="utf-8"))
            ref = m["completeness"].get("ref", {})
            meta = dict(n=ref.get("n"), seed=ref.get("seed"), code_hash=ref.get("code_hash"), bc_key=ref.get("bc_key"),
                        source=os.path.basename(p), source_kind="merged")
            seen = set()
            for r in m.get("pareto", []) + m.get("top", []):
                if r["idx"] not in seen:
                    seen.add(r["idx"]); compact.append(r)
            continue
        parts = M0.load_parts([p])
        ok, rep, rows = M0.check_complete(parts)
        if not ok:
            raise SystemExit(f"уровень 0 неполон: {rep['problems'][:3]} пропущено частей {rep['missing_shards'][:5]}")
        ref = rep["ref"]
        meta = dict(n=ref["n"], seed=ref["seed"], code_hash=ref["code_hash"], bc_key=ref["bc_key"],
                    source=p, source_kind="parts")
        full += list(rows.values())
    return full, compact, meta


def select(full, compact, k):
    """отбор k вариантов: фронт Парето, затем по кругу ниш по рейтингу -> [(idx, причина, оценка, кратко)]"""
    M0.HARD = list(CR.TARGETS["funnel_interface"]["hard"])
    cand = []
    if full:
        ok = [r for r in full if r["status"] == "ok"]
        import numpy as np
        S, RAW, PEN = [], [], []
        for r in ok:
            s, raw = M0.objectives(r)
            S.append(s); RAW.append(raw); PEN.append(M0.penalty(r))
        base = [M0.score(s, p, M0.WEIGHTS) for s, p in zip(S, PEN)]
        F = np.array([[(raw[q] if sense == "min" else -raw[q]) for q, sense in M0.PARETO_OBJ] for raw in RAW]) \
            if ok else np.zeros((0, 6))
        pf = set(int(i) for i in M0.pareto(F)) if len(ok) else set()
        for j, r in enumerate(ok):
            cand.append(dict(idx=r["idx"], score=float(base[j]), pareto=j in pf, niche=niche(r["concept"]),
                             l0=M0.compact(r, float(base[j]), PEN[j])))
    else:
        pf_idx = set()
        for r in compact:
            cand.append(dict(idx=r["idx"], score=r.get("score") or 0.0, pareto=False,
                             niche=f'{r["scheme"]}|{r["chamber"]}|{r["injector"]}|{r["passive"]}', l0=r))
        _ = pf_idx
    cand.sort(key=lambda x: -x["score"])
    chosen, why = [], {}
    for x in cand:
        if x["pareto"] and len(chosen) < k:
            chosen.append(x); why[x["idx"]] = "фронт Парето"
    by = {}
    for x in cand:
        if x["idx"] not in why:
            by.setdefault(x["niche"], []).append(x)
    order = sorted(by, key=lambda n: -by[n][0]["score"])
    rnd = 0
    while len(chosen) < k and any(by[n] for n in order):
        for n in order:
            if by[n] and len(chosen) < k:
                x = by[n].pop(0)
                chosen.append(x); why[x["idx"]] = f"ниша {n}, круг {rnd + 1}"
        rnd += 1
    return [(x["idx"], why[x["idx"]], x["score"], x["l0"]) for x in chosen]


def _hot(K, phi, p, T3):
    """состав и температура горячего поля зажигания: равновесие HP метан - воздух при phi зоны (Cantera)"""
    G = K.GM
    G.set_equivalence_ratio(max(phi, 1e-3), "CH4:1", CR.AIR)
    G.TP = T3, p
    G.equilibrate("HP")
    Y = {s: float(G.Y[G.species_index(s)]) for s in ("CH4", "O2", "N2", "AR", "CO2", "H2O", "CO", "H2", "OH", "NO")}
    return dict(T=max(float(G.T), T_IGN_MIN), T_eq=float(G.T), phi=phi, Y=Y)


def entry(i, n, seed, kv, bc, K, why, score, l0):
    c = CS.concept(i, n, kv, seed)
    d = CS.derive(c, bc)
    x100 = CS.fuel_split100(c, bc["loads"][0]["phi_ov"])
    loads = []
    for L in bc["loads"]:
        m_rel = L["W_fuel"] / bc["loads"][0]["W_fuel"]
        if c["scheme"] == "TVC":
            x = CR.pilot_share(x100, m_rel) if c["passive"] == "гидравлика" else x100
        elif c["passive"] == "гидравлика":
            x = CR.pilot_share(x100, m_rel)
        else:
            x = 0.0
        a_eff = CR.strat_air_frac(c, m_rel) if c["passive"] == "струи" else 1.0
        if c["scheme"] == "RQL":
            phi_z = (L["W_fuel"] / (c["a_pz"] * L["W_air"])) / bc["far_st_CH4"]
        elif c["passive"] == "гидравлика" and c["scheme"] in ("LPP", "MM", "FLOX"):
            phi_z = ((1 - x) * L["W_fuel"] / ((1 - c["a_aux"]) * c["a_pz"] * L["W_air"])) / bc["far_st_CH4"]
        else:
            phi_z = ((1 - x) * L["W_fuel"] / (c["a_pz"] * L["W_air"])) / bc["far_st_CH4"]
        loads.append(dict(load=L["load"], p3_Pa=L["p3_Pa"], T3_K=L["T3_K"], W_air=L["W_air"], W_fuel=L["W_fuel"],
                          rho3=L["rho3"], T4_K=L["T4_K"], phi_ov=L["phi_ov"], m_rel=m_rel, x_pilot=x, a_eff=a_eff,
                          phi_zone=phi_z, hot_init=_hot(K, phi_z / max(a_eff, 1e-6), L["p3_Pa"], L["T3_K"])))
    keep = ("d_air_mm", "u_a", "A_air_m2", "n_jet", "d_jet_mm", "p_fuel_ratio", "u_fuel", "J", "holdeman_C", "S_mix",
            "S_tol", "S_unmix", "fuel_split100", "V_L", "a_cool", "a_aux_out", "a_dil", "D_L", "L_L", "H_L", "D_m",
            "D_h", "A_L", "A_wall", "A_ref", "D_case", "D_env", "L_mod", "dp_diff", "dp_duct", "dp_series_wall",
            "dpL_qref", "phi_pz_100", "dp_mix_fric", "tau_mix_ms")
    der = {k: d[k] for k in keep if k in d}
    der["zones"] = d["zones"]
    return dict(idx=i, why=why, l0_score=score, l0=l0, concept=c, derived=der, loads=loads,
                n_can=CS.N_CAN[c["chamber"]], far_st=bc["far_st_CH4"], D2_m=bc["D2_m"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("paths", nargs="+", help="части уровня 0 (папки/файлы) или слитый файл merge_l0")
    ap.add_argument("--k", type=int, default=500)
    ap.add_argument("--kv", default="KV-3")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    full, compact, meta = _rows_from(a.paths)
    if meta.get("code_hash") and meta["code_hash"] != CS.code_hash():
        print(f"[внимание] код уровня 0 частей {meta['code_hash']} != текущий {CS.code_hash()}: варианты по индексу "
              f"пересчитываются текущим кодом", file=sys.stderr)
    sel = select(full, compact, a.k)
    bc = CR.boundary()
    K = CR.kin()
    ents = [entry(i, meta["n"], meta["seed"], a.kv, bc, K, why, sc, l0) for i, why, sc, l0 in sel]
    blob = json.dumps(ents, ensure_ascii=False, sort_keys=True).encode()
    res = dict(schema=SCHEMA, title=f"Входы уровня 1 воронки камеры {a.kv}: {len(ents)} вариантов уровня 0",
               method=("отбор - фронт Парето, затем по кругу ниш (схема × камера × форсунка × деление) по рейтингу "
                       "merge_l0; оси - comb_space.concept, геометрия - comb_space.derive, нагрузки - comb_bc_results.json, "
                       "пилот и самостратификация - comb_crn.pilot_share / strat_air_frac, горячее поле - равновесие HP "
                       "Cantera GRI-Mech 3.0 (T не ниже 1600 K, только начальное поле)"),
               software=f"python {sys.version.split()[0]}, cantera {CR.ct.__version__}", date=CS.DATE,
               status="входы расчёта (не результат)", kv=a.kv, l0=meta, code_hash_now=CS.code_hash(), bc_key=bc["key"],
               k=len(ents), fingerprint=hashlib.sha256(blob).hexdigest()[:16],
               niches=sorted(set(niche(e["concept"]) for e in ents)), entries=ents)
    json.dump(res, open(a.out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"записано {a.out}: {len(ents)} вариантов, ниш {len(res['niches'])}, отпечаток {res['fingerprint']}")


if __name__ == "__main__":
    main()
