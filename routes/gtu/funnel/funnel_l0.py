"""Воронка двигателя-заменителя, уровень 0, часть 2: расчёт концепций для рейтинга - часть выборки (на AWS).

Для каждой концепции (funnel_space.concept: Соболь, зерно SEED, вариант i - только от индекса):
  * цикл: ядро funnel_space на таблицах gtu_explore (сверено с gtu_cycle.cycle до 0,06 п. п.; выборочная сверка с
    gtu_cycle.cycle - в merge_l0.py); КПД ступеней компрессора - формула comp_stage (gtu_explore) с поправками
    funnel_space, плюс поправка за отход удельной быстроходности от оптимума -K·ln²(ns/ns_opt) и за радиальный зазор;
    ns_opt и K центробежной ступени подобраны по двум ступеням frame3 (calib_ns); радиальная турбина - U/C0
    (gtu_explore.turb_radial) и та же форма поправки по быстроходности; осевые турбины - eta_axial(ψ) gtu_explore;
  * размеры колёс, турбин, силовой турбины, рекуператора (funnel_space.derive), габарит L×W×H и масса двигателя (оценка);
  * детали (parts_model), монтаж (service_model), Plug&Play (plugplay_model), изготовимость (make_model);
  * четвёртый критерий владельца «КПД при минимальном размере и изготовимость»: s_eta (КПД на клеммах против исходной
    установки), s_size (доля объёма отсека двигателя), s_make; взвешенная оценка по весам владельца (WEIGHTS);
  * правила допустимости funnel_space.check: недопустимые - с кодом, величиной, пределом (reject, reject_detail).

Запуск части на AWS (несколько частей одним пулом процессов на все ядра машины):
  python funnel_l0.py --kv KV-3 --n 200000 --shards 30 --shard 0,3,6 --workers "$NPROC" --out "$GTU_OUT"
Результат не зависит от числа процессов: запись - чистая функция индекса, порядок - по индексу.
Слияние, рейтинг, ниши и отбор на уровень 1 - merge_l0.py. Всё - оценка (README узла, раздел «Уровень 0»).
"""
import argparse, datetime, gzip, hashlib, json, math, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import funnel_space as FS        # noqa: E402
import make_model as MK          # noqa: E402
import parts_model as PM         # noqa: E402
import plugplay_model as PP      # noqa: E402
import service_model as SM       # noqa: E402

DATE = FS.DATE
SEED = FS.SEED
LEVEL = "L0"
PART_PREFIX = "funnel_l0"

# =========================================================================== веса и шкалы
# решение владельца 14.11.2026: детали 0,35; Plug&Play 0,25; монтаж «хоть уборщиком» 0,20; КПД при минимальном размере и
# изготовимость 0,20. Внутреннее деление четвёртого критерия - оценка (поровну), его чувствительность - в merge_l0.py.
WEIGHTS = {"parts": 0.35, "plugplay": 0.25, "service": 0.20, "effmake": 0.20}
CRIT = list(WEIGHTS)
CRIT_NAMES = {"parts": "минимум отдельных деталей и узлов", "plugplay": "Plug&Play в существующий пакет",
              "service": "монтаж и ремонт «хоть уборщиком»", "effmake": "КПД при минимальном размере и изготовимость"}
CRIT_COL = {"parts": "s_parts", "plugplay": "s_plugplay", "service": "s_service", "effmake": "s_effmake"}
SPLIT_EFFMAKE = {"etasize": 0.5, "make": 0.5}      # оценка
SPLIT_ETASIZE = {"eta": 0.5, "size": 0.5}          # оценка
ETA_SPAN_PP = 5.0          # s_eta = 5 + 5·(КПД - КПД исходной)/5 п. п.: исходная 5, +5 п. п. - 10, -5 п. п. - 0 (оценка)
BAY_L_FRAC = 1.0 / 3.0     # отсек двигателя: min(ширина, высота пакета)² × треть длины пакета (оценка)
MASS_K = 1.25              # масса двигателя = 1,25 × сумма модулей (обвязка, опоры, трубы; оценка)


def check_weights(w):
    """веса критериев: ровно четыре ключа CRIT, неотрицательные, сумма 1 (до 1e-9)"""
    if set(w) != set(CRIT):
        raise ValueError(f"веса: нужны ключи {CRIT}, дано {sorted(w)}")
    if any((not isinstance(v, (int, float))) or v < 0 for v in w.values()):
        raise ValueError("веса неотрицательные числа")
    s = sum(w.values())
    if abs(s - 1.0) > 1e-9:
        raise ValueError(f"сумма весов {s:.12f} != 1")
    return w


def weighted(s, w=None):
    """взвешенная оценка 0-10 по критериям s (баллы 0-10); нормировка критерия - деление на 10 (шкала постоянная)"""
    w = check_weights(dict(w or WEIGHTS))
    return sum(w[k] * s[k] for k in CRIT)


def effmake(s_eta, s_size, s_make, split_em=None, split_es=None):
    em = split_em or SPLIT_EFFMAKE; es = split_es or SPLIT_ETASIZE
    s_es = es["eta"] * s_eta + es["size"] * s_size
    return em["etasize"] * s_es + em["make"] * s_make, s_es


# =========================================================================== удельная быстроходность
NS_OPT_T, K_NS_T = 0.60, 0.15   # радиальная турбина: оптимум ns (ω·√Q_вых/Δh_s^0,75) и коэффициент поправки - оценка
NS_T_BASIS = ("оценка: форма -K·ln²(ns/ns_opt), ns_opt 0,60 и K 0,15 - по памяти о кривых КПД радиальных турбин "
              "(Rohlik, NASA TN D-4384, 1968 - отчёт существует по выдержке поиска ntrs.nasa.gov/citations/19680006474; "
              "значение оптимума в выдержке не найдено, документ не открывался); проверить по SU2 на уровне 2")
_CAL = {}


def calib_ns():
    """ns_opt и K центробежной ступени по двум ступеням frame3 (frame3_point_results.json): у каждой ступени модель
    потерь frame3 даёт оптимум КПД при частоте n_ref и снижение d_eta_used при рабочей частоте n = k_speed·n_ref.
    При тех же расходе и напоре ns ∝ n, поэтому ns_opt = ns(n_ref), а ns/ns_opt = k_speed. K - МНК через ноль."""
    if _CAL:
        return _CAL
    E = FS._E()
    f3 = json.load(open(os.path.join(GTU, "frame3", "frame3_point_results.json"), encoding="utf-8"))
    pts = []
    for key in ("stage1", "stage2"):
        s = f3[key]; b = s["eta_basis"]
        T1 = s["T01_C"] + 273.15; p1 = s["p01_kPa"] * 1e3
        Q = s["m_per_side_kg_s"] / (p1 / (E.RA * T1))
        dhs = s["w_kJkg"] * 1e3 * s["eta_is"]
        ns_ref = (2 * math.pi * b["n_ref_rpm"] / 60) * math.sqrt(Q) / dhs ** 0.75
        pts.append(dict(stage=key, ns_ref=ns_ref, k_speed=b["k_speed"], d_eta=b["d_eta_used"], phi=s["phi"]))
    ns_opt = sum(p["ns_ref"] for p in pts) / len(pts)
    x2 = [math.log(p["k_speed"]) ** 2 for p in pts]
    K = sum(-p["d_eta"] * x for p, x in zip(pts, x2)) / sum(x * x for x in x2)
    for p, x in zip(pts, x2):
        p["fit"] = -K * x
    mu_c = sum(FS.MU_RANGE["Ц"]) / 2; mu_d = sum(FS.MU_RANGE["Д"]) / 2
    ns_d = ns_opt * math.sqrt(FS.PHI_OPT["Д"] / FS.PHI_OPT["Ц"]) * (mu_c / mu_d) ** 0.75
    _CAL.update(ns_opt_c=ns_opt, k_ns_c=K, ns_opt_d=ns_d, ns_opt_t=NS_OPT_T, k_ns_t=K_NS_T,
                points=[{k: (round(v, 5) if isinstance(v, float) else v) for k, v in p.items()} for p in pts],
                basis=dict(c="подбор по frame3 (две ступени): ns_opt = ns при n_ref модели потерь frame3, K - МНК по d_eta_used",
                           d="оценка: ns_opt центробежной × √(φ_opt Д / φ_opt Ц) × (μ Ц / μ Д)^0,75 по funnel_space.PHI_OPT и "
                             "серединам MU_RANGE",
                           t=NS_T_BASIS, axial="осевые ступени компрессора и турбины - без поправки по ns (формулы funnel_space)"))
    return _CAL


def ns_comp(stage_out, Gs, n):
    """ns = ω·√(Gs/ρ1)/Δh_s^0,75 ступени компрессора (на сторону при двустороннем входе)"""
    E = FS._E()
    T1, p1, pis = stage_out["T1"], stage_out["p1"], stage_out["pi"]
    dhs = E.h(E.T_s(E.s0(T1) + E.RA * math.log(pis))) - E.h(T1)
    Q = Gs / (p1 / (E.RA * T1))
    return (2 * math.pi * n / 60) * math.sqrt(Q) / max(dhs, 1.0) ** 0.75


def _clip(x, lo, hi):
    return min(hi, max(lo, x))


def ns_window_c():
    """окно ns/ns_opt центробежной и диагональной ступени, в котором действует поправка: из диапазона коэффициента
    расхода правила phi (funnel_space.PHI_RANGE) при ns ∝ √φ; за окном концепцию отсекает само правило phi"""
    lo, hi = FS.PHI_RANGE["Ц"]
    return math.sqrt(lo / FS.PHI_OPT["Ц"]), math.sqrt(hi / FS.PHI_OPT["Ц"])


def corrections(c, perf, g):
    """поправки КПД ступеней компрессора (зазор comp_stage + быстроходность) и радиальной турбины (быстроходность).
    Быстроходность в поправке ограничена окном правил (ns_window_c, NS_RANGE funnel_space): формула -K·ln² вне подобранной
    области не продолжается, такие концепции отсекают правила phi и turb_ns"""
    cal = calib_ns()
    wlo, whi = ns_window_c()
    dc, nsc = [], []
    for s, so in zip(g["stages"], perf["comp"]):
        if s["type"] == "О":
            dc.append(0.0); nsc.append(None)
            continue
        ns = ns_comp(so, s["Gs"], g["n"])
        opt = cal["ns_opt_c"] if s["type"] == "Ц" else cal["ns_opt_d"]
        cl = max(0.3e-3, 0.0007 * s["D"])
        dc.append(-0.3 * cl / s["b2"] - cal["k_ns_c"] * math.log(_clip(ns / opt, wlo, whi)) ** 2)
        nsc.append(ns)
    dt, nst = 0.0, None
    if c["turb"] == "Р" and g["tg"].get("ns"):
        nst = g["tg"]["ns"]
        dt = -cal["k_ns_t"] * math.log(_clip(nst, *FS.NS_RANGE) / cal["ns_opt_t"]) ** 2
    return dc, dt, nsc, nst


def perf_pass(c, tgt, d_eta_c=None, d_eta_t=0.0):
    """быстрый цикл как funnel_space.quick_perf, плюс поправка d_eta_t к КПД радиальной турбины газогенератора"""
    E = FS._E()
    pis = FS.stages_pi(c); typ = FS.stage_types(c)
    d_eta_c = d_eta_c or [0.0] * len(pis)
    mats = [c["mat_cw"] if not (c["mat_cw"] == "Al" and k > 0) else "Ti64" for k in range(len(pis))]
    mus = FS.stage_mus(c)
    comp = FS._comp([(p, max(0.60, FS.eta_comp(t, p, m, k, c["inlet"], mu) + de))
                     for k, (p, t, m, de, mu) in enumerate(zip(pis, typ, mats, d_eta_c, mus))])
    Ut, Ulim, _ = FS.turb_speed(c)
    Upt = c["kUpt"] * 380.0
    eg, ep = 0.86, 0.88
    r, load = None, {}
    for _ in range(8):
        r = FS._core(c, comp, eg, ep)
        if "fail" in r:
            break
        if c["turb"] == "О":
            psi = r["t"]["w"] / (c["nt"] * Ut ** 2)
            eg_n = E.eta_axial(psi) - FS.leaving_loss(c["phit"], Ut, r["t"]["w"] / eg)
            load = dict(psi=psi)
        else:
            eg_n, nu, dhs = FS._radial_stage(r["t"], c["nt"], Ut)
            eg_n += d_eta_t
            load = dict(nu=nu, dhs=dhs)
        ep_n = ep
        if r["pt"]:
            psi_pt = r["pt"]["w"] / (c["npt"] * Upt ** 2)
            ep_n = E.eta_axial(psi_pt) - FS.leaving_loss(c["phit"], Upt, r["pt"]["w"] / ep)
            load["psi_pt"] = psi_pt
        eg_n, ep_n = max(0.40, eg_n), max(0.40, ep_n)
        if abs(eg_n - eg) < 1e-4 and abs(ep_n - ep) < 1e-4:
            break
        eg, ep = eg_n, ep_n
    if r is None or "fail" in r:
        return dict(ok=False, fail=(r or {}).get("fail", "нет решения"), comp=comp)
    G = tgt["P_shaft_kW"] * 1e3 / r["wnet"]
    return dict(ok=True, comp=comp, r=r, eta_gg=eg, eta_pt=ep if r["pt"] else None, load=load, G=G, Ut=Ut, Ulim_t=Ulim,
                Upt=Upt, eta_th=r["eta"], w=r["wnet"] / 1e3, EGT=r["Tstack"] - 273.15, T5=r["T5"] - 273.15,
                T2=r["T2"] - 273.15)


# =========================================================================== габарит и масса (оценка)
def envelope(c, tgt, g):
    """габарит L×W×H двигателя и масса (оценка): длина ротора + 0,6·D2max (вход и выхлопной коллектор) [+ барабан
    волнового ротора 0,6·D2max] [+ силовая турбина 0,15·D·ступеней + 0,5·Dtip переходного канала]; ширина и высота -
    D_env funnel_space (корпус + 0,5 м обвязки; кластер, СТ, рекуператор в корпусе); выносная камера - над корпусом
    (+1,2·D камеры); вынесенный рекуператор - сбоку рядами кубов модулей + 0,3 м; масса - MASS_K × сумма модулей"""
    Dc, D = g["Dc"], g["D_env"]
    L = g["L_rot"] + 0.6 * Dc
    if c["comb"] == "ВР":
        L += 0.6 * Dc
    if g["pt"]:
        L += 0.15 * g["pt"]["D"] * c["npt"] + 0.5 * g["pt"]["Dtip"]
    W = H = D
    if c["comb"] == "выносная":
        H += 1.2 * g["D_comb"]
    rec = g["rec"]
    if rec and not c["recup_integral"]:
        along = max(1, int(L // max(rec["side"], 1e-6)))
        W += math.ceil(rec["modules"] / along) * rec["side"] + 0.3
    m = g["mass"]
    core = m["rotor"] + m["case"] + m["hot"] + m["comb"]
    mass = MASS_K * (c["cores"] * core + m["pt"] + (rec["mass_kg"] if rec else 0.0))
    V = L * W * H
    V_bay = min(tgt["pkg_W"], tgt["pkg_H"]) ** 2 * tgt["pkg_L"] * BAY_L_FRAC
    return dict(L=L, W=W, H=H, V=V, V_bay=V_bay, mass=mass)


def s_eta_of(eta_el, eta0):
    return min(10.0, max(0.0, 5.0 + 5.0 * (eta_el - eta0) / ETA_SPAN_PP))


def s_size_of(V, V_bay):
    return 10.0 * min(1.0, max(0.0, 1.0 - V / V_bay))


# =========================================================================== ниша (архив MAP-Elites)
NICHE_AXES = ["arch", "rotors", "comp", "turb", "cool", "recup", "rotor_make"]


def rotor_make(c):
    if c["rotor_join"] in ("заодно", "сварка"):
        return f"{c['rotor_join']}:{c['m_tw']}"
    return f"сборный:{c['m_cw']}/{c['m_tw']}"


def niche_of(c):
    """ниша = архетип × число роторов × тип компрессора × тип турбины × класс охлаждения × рекуператор × изготовление ротора"""
    return "|".join([c["arch"], f"{c['rotors']}р", c["comp"], c["turb"], c["cool"], c["recup"], rotor_make(c)])


# =========================================================================== запись концепции
AX_COLS = [a for a in FS.AXES if a != "arch"]
COLUMNS = (["idx", "kv", "arch", "niche"] + [a for a in AX_COLS if a in FS.CAT] + ["pi_max_U"] + [a for a in AX_COLS if a in FS.CONT]
           + ["feasible", "reject", "reject_detail", "violation",
              "eta_th_shaft", "eta_el", "w_kJkg", "G", "dG_pct", "EGT", "dEGT", "T2", "T5", "pi_stages", "eta_c", "ns_c",
              "deta_c", "eta_t", "ns_t", "deta_t", "eta_pt",
              "n_rpm", "D2", "U2", "b2_mm", "phi", "D_turb", "U_turb", "h_turb", "D_pt", "h_pt", "D_case", "D_env",
              "L_env", "W_env", "H_env", "V_env", "mass_kg", "rotor_kg", "rec_modules", "rec_V",
              "s_parts", "s_plugplay", "s_service", "s_make", "s_eta", "s_size", "s_etasize", "s_effmake", "score",
              "n_parts", "n_assemblies", "n_moving", "n_cooled", "n_fasteners", "n_unique", "n_eff",
              "heaviest_kg", "heaviest_module", "steps", "n_tools", "wrongable", "n_modules",
              "pp_blocked", "pp_no_data", "trl_min", "weeks", "techs"])


def _r(x, k=4):
    return None if x is None else (round(float(x), k) if isinstance(x, float) else x)


def _j(vals, k=4):
    """список чисел -> строка «a/b» (для CSV и компактности)"""
    return "/".join("" if v is None else f"{round(v, k):g}" for v in vals)


def _detail_rules(bad, k=3):
    out = []
    for r in bad[:k]:
        v, lim = r.get("val"), r.get("lim")
        if isinstance(v, (int, float)) and lim is not None:
            out.append(f"{r['code']}: {v:g} при пределе {lim if not isinstance(lim, list) else '–'.join(f'{x:g}' for x in lim)} ({r['text']})")
        else:
            out.append(f"{r['code']}: {r.get('text', '')}")
    return "; ".join(out)


def solve(c, tgt):
    """цикл и геометрия с поправками КПД: проход без поправок -> геометрия -> поправки -> проход (два раза).
    Возврат: dict(perf, g, dc, dt, nsc, nst, fail); g = None, если цикл или геометрия не строятся"""
    perf, g, dc, dt, nsc, nst, fail = None, None, None, 0.0, None, None, None
    try:
        p = perf_pass(c, tgt)
        if p["ok"]:
            g = FS.derive(c, tgt, p)
            for _ in range(2):
                dc, dt, nsc, nst = corrections(c, p, g)
                p = perf_pass(c, tgt, dc, dt)
                if not p["ok"]:
                    g = None
                    break
                g = FS.derive(c, tgt, p)
        perf = p
        if not p["ok"]:
            fail = p["fail"]
    except (ValueError, ZeroDivisionError, OverflowError) as e:
        g, fail = None, f"геометрия не строится: {e}"
    return dict(perf=perf, g=g, dc=dc, dt=dt, nsc=nsc, nst=nst, fail=fail)


def evaluate_l0(i, n, kv="KV-3", seed=SEED, detail=False):
    """концепция i выборки объёма n -> плоская запись COLUMNS (dict); detail - модули, интерфейсы, геометрия"""
    c = FS.concept(i, n, kv, seed)
    tgt = FS.target(kv)
    rec = dict(idx=int(i), kv=kv, arch=c["arch"], niche=niche_of(c))
    for a in COLUMNS[4:]:
        if a in c:
            rec[a] = _r(c[a], 6)
    sol = solve(c, tgt)
    perf, g, dc, dt, nsc, nst, fail = (sol[k] for k in ("perf", "g", "dc", "dt", "nsc", "nst", "fail"))
    if g is None:
        bad = [FS._rule("perf", False, text=fail or "нет решения")]
    else:
        bad = [r for r in FS.check(c, tgt, perf, g) if not r["ok"]]
    rec.update(feasible=int(not bad), reject="|".join(r["code"] for r in bad), reject_detail=_detail_rules(bad),
               violation=FS.violation(bad))
    if g is None:
        for k in COLUMNS:
            rec.setdefault(k, None)
        return ({k: rec[k] for k in COLUMNS}, None) if detail else {k: rec[k] for k in COLUMNS}
    env = envelope(c, tgt, g)
    eta_el = perf["eta_th"] * 100 * 0.985 * 0.972
    st = g["stages"]
    rec.update(eta_th_shaft=_r(perf["eta_th"] * 100, 3), eta_el=_r(eta_el, 3), w_kJkg=_r(perf["w"], 2), G=_r(perf["G"], 3),
               dG_pct=_r((perf["G"] / tgt["G0"] - 1) * 100, 2), EGT=_r(perf["EGT"], 1), dEGT=_r(perf["EGT"] - tgt["EGT0"], 1),
               T2=_r(perf["T2"], 1), T5=_r(perf["T5"], 1), pi_stages=_j([s["pi"] for s in st], 3),
               eta_c=_j([s["eta"] for s in st]), ns_c=_j(nsc or [], 3), deta_c=_j(dc or [], 4), eta_t=_r(perf["eta_gg"], 4),
               ns_t=_r(nst, 3), deta_t=_r(dt, 4), eta_pt=_r(perf["eta_pt"], 4),
               n_rpm=_r(float(g["n"]), 0), D2=_j([s["D"] for s in st]), U2=_j([s["U"] for s in st], 1),
               b2_mm=_j([s["b2"] * 1e3 if "b2" in s else None for s in st], 1),
               phi=_j([s.get("phi") for s in st]), D_turb=_r(g["tg"]["D"]), U_turb=_r(g["tg"]["U"], 1),
               h_turb=_r(g["tg"]["h"]), D_pt=_r(g["pt"]["D"]) if g["pt"] else None, h_pt=_r(g["pt"]["h"]) if g["pt"] else None,
               D_case=_r(g["D_case"], 3), D_env=_r(g["D_env"], 3), L_env=_r(env["L"], 3), W_env=_r(env["W"], 3),
               H_env=_r(env["H"], 3), V_env=_r(env["V"], 3), mass_kg=_r(env["mass"], 0), rotor_kg=_r(g["mass"]["rotor"], 1),
               rec_modules=g["rec"]["modules"] if g["rec"] else 0, rec_V=_r(g["rec"]["V_m3"], 3) if g["rec"] else None)
    bom, sv, pp, mk = FS.model_inputs(c, tgt, perf, g)
    pr, sr, ppr, mr = PM.count(bom), SM.score(sv), PP.score(pp), MK.score(mk)
    s_eta = round(s_eta_of(eta_el, tgt["eta0"]), 3)
    s_size = round(s_size_of(env["V"], env["V_bay"]), 3)
    s_em, s_es = effmake(s_eta, s_size, mr["score"])
    s = dict(parts=pr["score"], plugplay=ppr["score"], service=sr["score"], effmake=round(s_em, 3))
    rec.update(s_parts=pr["score"], s_plugplay=ppr["score"], s_service=sr["score"], s_make=mr["score"], s_eta=s_eta, s_size=s_size,
               s_etasize=round(s_es, 3), s_effmake=round(s_em, 3), score=round(weighted(s), 4),
               n_parts=pr["parts"], n_assemblies=pr["assemblies"], n_moving=pr["moving"], n_cooled=pr["cooled"],
               n_fasteners=pr["fasteners"], n_unique=pr["unique"], n_eff=_r(float(pr["n_eff"]), 1),
               heaviest_kg=sr["heaviest"]["mass_kg"], heaviest_module=sr["heaviest"]["module"], steps=sr["steps"],
               n_tools=len(sr["tools"]), wrongable=sr["wrongable"], n_modules=len(sv["modules"]),
               pp_blocked="|".join(ppr["blocked"]), pp_no_data="|".join(ppr["no_data"]), trl_min=mr["trl_min"],
               weeks=mr["weeks"], techs="|".join(mk["techs"]))
    row = {k: rec.get(k) for k in COLUMNS}
    if not detail:
        return row
    geo = geometry(c, tgt, perf, g, env)
    det = dict(concept=c, geometry=geo,
               parts=dict(by_unit={u: {k: v for k, v in d.items() if k != "rule"} for u, d in pr["by_unit"].items()},
                          welds=pr["welds"]),
               service=dict(modules=sr["modules"], frequent=sr["frequent"], tools=sr["tools"],
                            sub=dict(s_mass=sr["s_mass"], s_steps=sr["s_steps"], s_tools=sr["s_tools"], s_err=sr["s_err"],
                                     s_diag=sr["s_diag"])),
               plugplay={k: dict(state=v["state"], note=v.get("note")) for k, v in ppr["items"].items()},
               make=dict(trl_by=mr["trl_by"], s_fit=mr["s_fit"], s_proc=mr["s_proc"], s_trl=mr["s_trl"], s_time=mr["s_time"]),
               rules_failed=bad)
    return row, det


def geometry(c, tgt, perf, g, env):
    """размеры крупных узлов, м (оценка уровня 0; уровень 1 дополняет для эскизов)"""
    st = g["stages"]
    return dict(
        inlet=dict(type="двусторонний" if c["inlet"] == "двух" else "односторонний"),
        compressor=[dict(type=s["type"], D2=_r(s["D"]), b2=_r(s.get("b2")), U2=_r(s["U"], 1), phi=_r(s.get("phi"))) for s in st],
        D_case=_r(g["D_case"], 3), L_case=_r(g["L_case"], 3), L_rotor=_r(g["L_rot"], 3),
        combustor=dict(type=c["comb"], D=_r(g["D_comb"], 3), L=_r(g["L_comb"], 3), n_injectors=g["n_inj"]),
        turbine=dict(type=c["turb"], stages=c["nt"], D=_r(g["tg"]["D"]), Dtip=_r(g["tg"]["Dtip"]), h_last=_r(g["tg"]["h"]),
                     D_hot=_r(g["D_hot"], 3), L_hot=_r(g["L_hot"], 3)),
        power_turbine=(dict(stages=c["npt"], D=_r(g["pt"]["D"]), Dtip=_r(g["pt"]["Dtip"]), h_last=_r(g["pt"]["h"]),
                            n_rpm=tgt["n_pkg"]) if g["pt"] else None),
        recuperator=(dict(modules=g["rec"]["modules"], side=_r(g["rec"]["side"], 3), V=_r(g["rec"]["V_m3"], 3),
                          integral=bool(c["recup_integral"])) if g["rec"] else None),
        cores=c["cores"], shaft_rpm=round(g["n"]), output_rpm=tgt["n_pkg"],
        envelope=dict(L=_r(env["L"], 3), W=_r(env["W"], 3), H=_r(env["H"], 3)), mass_kg_estimate=round(env["mass"]),
        note="оценка уровня 0: формулы derive/envelope (README узла)")


# =========================================================================== части и пул процессов
def code_sha():
    """отпечаток кода расчёта: слияние сверяет его у всех частей"""
    h = hashlib.sha256()
    for f in ("funnel_l0.py", "funnel_space.py", "parts_model.py", "service_model.py", "plugplay_model.py", "make_model.py"):
        h.update(open(os.path.join(HERE, f), "rb").read())
    return h.hexdigest()[:16]


def parse_shards(spec, shards):
    """«3», «0,3,6», «0-9», «all» -> отсортированный список номеров частей"""
    if spec in ("all", "*"):
        return list(range(shards))
    out = set()
    for tok in str(spec).split(","):
        tok = tok.strip()
        if not tok:
            continue
        if "-" in tok:
            a, b = tok.split("-", 1)
            out.update(range(int(a), int(b) + 1))
        else:
            out.add(int(tok))
    bad = [k for k in out if not 0 <= k < shards]
    if bad or not out:
        raise ValueError(f"номера частей {sorted(bad) or spec} вне 0..{shards - 1}")
    return sorted(out)


def _eval_chunk(args):
    idx, n, kv, seed = args
    return [evaluate_l0(i, n, kv, seed) for i in idx]


def evaluate_many(indices, n, kv, seed=SEED, workers=1):
    """записи по индексам; workers процессов; порядок и содержимое не зависят от workers"""
    indices = list(indices)
    if workers <= 1 or len(indices) < 2:
        return _eval_chunk((indices, n, kv, seed))
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    FS._E(); calib_ns(); FS.target(kv)            # таблицы и калибровка - до ветвления процессов
    size = max(1, math.ceil(len(indices) / (workers * 8)))
    chunks = [indices[k:k + size] for k in range(0, len(indices), size)]
    ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else None
    with ProcessPoolExecutor(workers, mp_context=ctx) as ex:
        got = [r for part in ex.map(_eval_chunk, [(ch, n, kv, seed) for ch in chunks]) for r in part]
    return got


def default_workers():
    for k in ("NPROC", "GTU_WORKERS"):
        v = os.environ.get(k)
        if v and v.isdigit() and int(v) > 0:
            return int(v)
    return os.cpu_count() or 1


def method_text():
    return ("воронка, уровень 0: концепции по индексу из последовательности Соболя (funnel_space, зерно SEED); быстрый цикл "
            "на таблицах gtu_explore (Cantera GRI-Mech 3.0, метан, ISO) с КПД ступеней comp_stage, поправками за зазор и "
            "удельную быстроходность (подбор по frame3) и КПД турбин по нагрузке и U/C0; геометрия derive, габарит и масса "
            "(оценка); правила допустимости funnel_space.check; баллы parts/service/plugplay/make, s_eta, s_size; взвешенная "
            "оценка по весам владельца 14.11.2026. Всё - оценка (README узла).")


def part_path(out, kv, n, shard, shards, gz=False):
    return os.path.join(out, f"{PART_PREFIX}_{kv}_n{n}_s{shard:03d}of{shards:03d}.json" + (".gz" if gz else ""))


def run_parts(kv, n, shard_list, shards, out, seed=SEED, workers=1, gz=False):
    """части shard_list из shards одним пулом; по файлу на часть; возврат: пути, сводка времени"""
    t0 = time.time()
    idx_by = {k: FS.shard_indices(n, k, shards) for k in shard_list}
    allidx = sorted(i for v in idx_by.values() for i in v)
    rows = evaluate_many(allidx, n, kv, seed, workers)
    by_idx = {r["idx"]: r for r in rows}
    el = time.time() - t0
    os.makedirs(out, exist_ok=True)
    meta = dict(method=method_text(), software=FS.software(), date=DATE, level=LEVEL, kv=kv, target=FS.target(kv), n=n, seed=seed,
                shards=shards, fingerprint=FS.fingerprint(n, seed), code_sha=code_sha(), calib=calib_ns(),
                weights=WEIGHTS, split_effmake=SPLIT_EFFMAKE, split_etasize=SPLIT_ETASIZE, workers=workers,
                call=dict(shards=shard_list, concepts=len(allidx), elapsed_s=round(el, 2),
                          s_per_concept_core=round(el * workers / max(1, len(allidx)), 5),
                          host_cpus=os.cpu_count(), started=datetime.datetime.fromtimestamp(t0, datetime.timezone.utc).isoformat()))
    paths = []
    for k in shard_list:
        doc = dict(meta, shard=k, columns=COLUMNS, n_rows=len(idx_by[k]), rows=[[by_idx[i][c] for c in COLUMNS] for i in idx_by[k]])
        p = part_path(out, kv, n, k, shards, gz)
        txt = json.dumps(doc, ensure_ascii=False, separators=(",", ":"))
        if gz:
            with gzip.open(p, "wt", encoding="utf-8") as fh:
                fh.write(txt)
        else:
            with open(p, "w", encoding="utf-8") as fh:
                fh.write(txt)
        paths.append(p)
    return paths, meta["call"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--kv", default="KV-3", choices=sorted(FS.KV))
    ap.add_argument("--n", type=int, required=True, help="объём выборки (концепции 0..n-1)")
    ap.add_argument("--shards", type=int, required=True, help="число частей прогона (10-50 на AWS)")
    ap.add_argument("--shard", default="all", help="номер части или список: 3 | 0,3,6 | 0-9 | all")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--out", default=os.environ.get("GTU_OUT", "."))
    ap.add_argument("--workers", type=int, default=0, help="процессов (0 - $NPROC или все ядра)")
    ap.add_argument("--gzip", action="store_true", help="файлы частей .json.gz")
    a = ap.parse_args(argv)
    if a.n < 1 or a.shards < 1:
        ap.error("n и shards >= 1")
    sl = parse_shards(a.shard, a.shards)
    w = a.workers or default_workers()
    paths, call = run_parts(a.kv, a.n, sl, a.shards, a.out, a.seed, w, a.gzip)
    print(f"{a.kv} n={a.n}: части {sl[0]}..{sl[-1]} ({len(sl)} из {a.shards}), {call['concepts']} концепций, {call['elapsed_s']} с, "
          f"{w} процессов, {call['s_per_concept_core']} с на концепцию на ядро; файлов {len(paths)} -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
