"""Воронка двигателя-заменителя, уровень 1: слияние частей, контроль полноты, отбор 50 разных, масштабирование на KV-1/KV-6.

  python merge_l1.py --kv KV-3 --m 4000 --shards 40 --in runs/<run_id> --out-dir l1_aws [--compare-in runs/<run_id_2x>]
  python merge_l1.py --probe          # короткая локальная проба: 12 вариантов отбора уровня 0 (AWS, N = 2000), 2 старта

Слияние: все части funnel_l1<tag>_<KV>_m<M>_s<k>of<S>.json[.gz]; у всех одинаковы вход (sha256 отбора уровня 0), порядок,
отпечаток кода, число стартов; у части k - ровно номера j % S == k; ни один номер 0..M-1 не потерян и не задвоен (копия части
после повтора задачи - с тем же содержанием; время счёта t_s, t_opt_s не сравнивается). Нарушение - отказ (код 2).
Отбор 50 разных (п. 4 приоритета владельца): кандидаты - допустимые на уровне 1 с оценкой не ниже SCORE_MIN; жадно: на
каждом шаге - наибольшее «оценка - LAMBDA·exp(-d/D0)», d - расстояние до ближайшего уже выбранного (дискретные оси - по
Хэммингу с весами, непрерывные - нормированные); не больше NICHE_MAX на нишу. Цена разнообразия - против чистого топ-50.
Масштабирование 50 отобранных на KV-1 и KV-6: та же дискретная схема, непрерывные оси - окна design() и оптимизация.
Пишет funnel_l1_results<sfx>.json, funnel_l1<sfx>.csv (все варианты), funnel_l1_top50<sfx>.json и .csv.

Пересчёт без нового счёта (v3, решение владельца 15.11.2026, п. 1 - множитель готовности TRL):
  python merge_l1.py --rescore --kv KV-3 --m 3983 --shards 40 --in aws/l1_full --out-dir aws/l1_full_rescore_v3 --suffix _rescore_v3 \
      [--fea-results aws/l1_fea/funnel_l1_fea_results.json]
  итоговая оценка каждой записи пересчитывается функцией funnel_l0.final_score по её частным баллам (s_parts, s_plugplay, s_service,
  s_effmake) и TRL_min; отбор 50 разных, цена разнообразия и сходимость - заново; масштабирование на KV-1/KV-6 и CalculiX не
  запускаются (это новый счёт). Записи без частных баллов (точка не построилась) оценки не получают и в отбор не идут.
"""
import argparse, csv, glob, gzip, hashlib, json, math, os, re, sys, tempfile, shutil, time
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import funnel_space as FS        # noqa: E402
import funnel_l0 as L0           # noqa: E402
import funnel_l1 as L1           # noqa: E402
import meanline_l1 as ML         # noqa: E402
import strength_l1 as ST         # noqa: E402
import offdesign_l1 as OD        # noqa: E402

DATE = FS.DATE
TOP_N = 50
NICHE_MAX = 3
SCORE_MIN = 5.0
SCORE_MIN_BASIS = ("5,0 - середина шкалы 0-10, к которой привязаны модели: деталей - опорная ведомость frame3 даёт 5 (parts_model.N_HALF), "
                   "КПД - исходная установка даёт 5 (s_eta); вариант ниже 5 по взвешенной оценке в среднем не лучше опорной конструкции "
                   "по критериям владельца - нести его на уровень 2 незачем. Плюс допустимость на уровне 1 (все правила, включая "
                   "нерасчётные режимы, запас по помпажу, пуск и прочность)")
LAMBDA, D0 = 1.0, 0.15
DIV_BASIS = ("штраф LAMBDA·exp(-d/D0), LAMBDA = 1,0 (единицы оценки), D0 = 0,15: двойник (d -> 0) теряет 1,0 балла, вариант на расстоянии "
             "0,15 - 0,37, на 0,3 - 0,14, на 0,5 - 0,04 (оценка: клоны вытесняются, разные почти не штрафуются)")
DISC_W = {"arch": 3.0, "rotors": 2.0, "comp": 2.0, "turb": 2.0, "cool": 2.0, "recup": 2.0, "rotor_join": 1.0, "m_cw": 1.0, "m_tw": 1.0,
          "cores": 1.0, "nt": 1.0, "npt": 1.0, "comb": 1.0, "bear": 1.0, "inlet": 1.0, "n_ax": 1.0, "mat_tw": 1.0, "mat_cw": 1.0,
          "m_hot": 1.0, "mat_hot": 1.0, "m_comb": 1.0, "m_case": 1.0, "hot_integral": 1.0, "recup_integral": 1.0, "conn": 1.0,
          "cartridge": 1.0, "split_bl": 1.0, "case_arch": 1.0}   # v3: ось «корпус» (вместе с cartridge разница картридж/оболочка - 2,0)
CONT_N = {"t41": (850.0, 1350.0, "lin"), "pi": (2.5, 22.0, "log"), "eps_rec": (0.0, 0.92, "lin"), "eta_th_shaft": (20.0, 45.0, "lin"),
          "V_env": (2.0, 60.0, "log"), "mass_kg": (500.0, 20000.0, "log")}
CONT_W = 4.0
EXTRA = ["div_rank", "pure_rank", "eligible", "select_reason"]


class MergeError(Exception):
    pass


def find_parts(dirs, kv, M, shards, tag=""):
    pat = re.compile(rf"^{L1.PART_PREFIX}{re.escape(tag)}_{re.escape(kv)}_m{M}_s(\d{{3}})of{shards:03d}\.json(\.gz)?$")
    out = []
    for d in dirs:
        for p in glob.glob(os.path.join(d, "**", f"{L1.PART_PREFIX}{tag}_{kv}_m{M}_s*of{shards:03d}.json*"), recursive=True):
            m = pat.match(os.path.basename(p))
            if m:
                out.append((int(m.group(1)), p))
    return sorted(out)


def load_part(p):
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt", encoding="utf-8") as fh:
        return json.load(fh)


TIMING_COLS = ("t_s", "t_opt_s", "t_cpu_s")


def content_hash(d):
    """отпечаток содержания части без времени счёта: повтор части (прерывание спота, другое число процессов) даёт те же строки и
    разборы, но другое время - такой повтор считается той же частью (замечание скептика отбора, 15.11.2026: слияние отказывало)"""
    idx = [i for i, c in enumerate(d["columns"]) if c not in TIMING_COLS]
    rows = [[r[i] for i in idx] for r in d["rows"]]
    det = [{k: v for k, v in (x or {}).items() if k != "t_s"} for x in d.get("details") or []]
    return hashlib.sha256(json.dumps([rows, det], ensure_ascii=False, sort_keys=True).encode()).hexdigest()


def merge_parts(found, kv, M, shards):
    checks = dict(files=len(found), shards_expected=shards, missing=[], duplicate_same=[], duplicate_conflict=[], meta_mismatch=[],
                  index_errors=[], unreadable=[])
    ref, by = None, {}
    keys = ("kv", "M", "shards", "input_sha256", "order_sha", "code_sha", "starts", "maxfev", "de_gens", "level", "top")
    for k, p in found:
        try:
            d = load_part(p)
        except (OSError, EOFError, ValueError) as e:
            checks["unreadable"].append(dict(file=os.path.basename(p), error=type(e).__name__))
            continue
        meta = {x: d.get(x) for x in keys + ("shard",)}
        bad = {}
        if meta["shard"] != k or meta["kv"] != kv or meta["M"] != M or meta["shards"] != shards or meta["level"] != L1.LEVEL:
            bad["id"] = [meta["shard"], meta["kv"], meta["M"], meta["shards"]]
        if ref is None:
            ref = dict(d)
            ref.pop("rows"); ref.pop("details")
        else:
            for x in keys:
                if meta[x] != ref.get(x):
                    bad[x] = [meta[x], ref.get(x)]
            if d["columns"] != ref["columns"]:
                bad["columns"] = "другой состав"
        if bad:
            checks["meta_mismatch"].append(dict(file=os.path.basename(p), diff=bad))
            continue
        exp = FS.shard_indices(M, k, shards)
        got = [r[0] for r in d["rows"]]
        if got != exp:
            checks["index_errors"].append(dict(file=os.path.basename(p), shard=k, lost=sorted(set(exp) - set(got))[:20],
                                               foreign=sorted(set(got) - set(exp))[:20], dup=len(got) - len(set(got))))
            continue
        h = content_hash(d)
        if k in by:
            (checks["duplicate_same"] if by[k]["hash"] == h else checks["duplicate_conflict"]).append(dict(shard=k, file=os.path.basename(p)))
            continue
        by[k] = dict(rows=d["rows"], details=d["details"], hash=h, call=d.get("call"), workers=d.get("workers"))
    checks["missing"] = [k for k in range(shards) if k not in by]
    errs = {x: checks[x] for x in ("missing", "duplicate_conflict", "meta_mismatch", "index_errors") if checks[x]}
    if errs or ref is None:
        raise MergeError(json.dumps(errs or {"files": "частей не найдено"}, ensure_ascii=False))
    rows = sorted((r for v in by.values() for r in v["rows"]), key=lambda r: r[0])
    if [r[0] for r in rows] != list(range(M)):
        raise MergeError(f"номера после слияния не равны 0..{M - 1}")
    det = {}
    for v in by.values():
        for dd in v["details"]:
            if dd is not None:
                det[dd["j"]] = dd
    checks.update(indices_ok=True, n_rows=len(rows), duplicates=0)
    calls, seen = [], set()
    for k, v in sorted(by.items()):
        c = v["call"] or {}
        key = (c.get("started"), tuple(c.get("shards") or []))
        if key not in seen:
            seen.add(key)
            calls.append(dict(c, workers=v["workers"]))
    return ref["columns"], rows, det, ref, checks, calls


# =========================================================================== расстояние и отбор
def _cn(name, v):
    lo, hi, kind = CONT_N[name]
    if v is None:
        return 0.0
    v = float(v)
    if kind == "log":
        return min(1.0, max(0.0, (math.log(max(v, 1e-9)) - math.log(lo)) / (math.log(hi) - math.log(lo))))
    return min(1.0, max(0.0, (v - lo) / (hi - lo)))


def features(rec):
    return dict(disc={a: rec.get(a) for a in DISC_W}, cont={a: _cn(a, rec.get(a)) for a in CONT_N})


def distance(fa, fb):
    """взвешенное расстояние 0-1: дискретные оси - Хэмминг с весами DISC_W, непрерывные - корень среднего квадрата разностей
    нормированных величин CONT_N с общим весом CONT_W"""
    dd = sum(w for a, w in DISC_W.items() if fa["disc"][a] != fb["disc"][a])
    dc = math.sqrt(sum((fa["cont"][a] - fb["cont"][a]) ** 2 for a in CONT_N) / len(CONT_N))
    return (dd + CONT_W * dc) / (sum(DISC_W.values()) + CONT_W)


def diff_axes(ra, rb):
    out = [f"{a}: {ra.get(a)} -> {rb.get(a)}" for a in DISC_W if ra.get(a) != rb.get(a)]
    for a in ("t41", "pi", "eps_rec"):
        if ra.get(a) is not None and rb.get(a) is not None and abs(float(ra[a]) - float(rb[a])) > 0.05 * max(abs(float(ra[a])), 1e-9):
            out.append(f"{a}: {float(ra[a]):.4g} -> {float(rb[a]):.4g}")
    return out


def eligible(recs, score_min=SCORE_MIN):
    return [r for r in recs if r.get("feasible") and r.get("score") is not None and r["score"] >= score_min]


def greedy(recs, n=TOP_N, niche_max=NICHE_MAX, lam=LAMBDA, d0=D0, score_min=SCORE_MIN):
    """жадный отбор разных; детерминирован (равенство - по меньшему j)"""
    cand = sorted(eligible(recs, score_min), key=lambda r: (-r["score"], r["j"]))
    F = {r["j"]: features(r) for r in cand}
    sel, dmin, cnt = [], {r["j"]: 1e9 for r in cand}, Counter()
    while len(sel) < n:
        best, bv = None, None
        for r in cand:
            if r["j"] in {s["j"] for s in sel} or cnt[r["niche"]] >= niche_max:
                continue
            pen = lam * math.exp(-dmin[r["j"]] / d0) if sel else 0.0
            v = r["score"] - pen
            if bv is None or v > bv + 1e-12 or (abs(v - bv) <= 1e-12 and r["j"] < best["j"]):
                best, bv = r, v
        if best is None:
            break
        sel.append(dict(best, _adj=bv, _pen=best["score"] - bv))
        cnt[best["niche"]] += 1
        for r in cand:
            dmin[r["j"]] = min(dmin[r["j"]], distance(F[r["j"]], F[best["j"]]))
    return sel


def pairwise(sel):
    F = [features(r) for r in sel]
    d = [distance(F[a], F[b]) for a in range(len(F)) for b in range(a + 1, len(F))]
    return dict(mean=round(float(np.mean(d)), 4) if d else None, min=round(float(min(d)), 4) if d else None)


def price_of_diversity(recs, sel, n=TOP_N, score_min=SCORE_MIN):
    pure = sorted(eligible(recs, score_min), key=lambda r: (-r["score"], r["j"]))[:n]
    pj, sj = {r["j"] for r in pure}, {r["j"] for r in sel}
    displaced = [r for r in pure if r["j"] not in sj]
    new = [r for r in sel if r["j"] not in pj]
    ps, ss = sorted([r["score"] for r in pure], reverse=True), sorted([r["score"] for r in sel], reverse=True)
    k = min(len(ps), len(ss))
    loss = [ps[i] - ss[i] for i in range(k)]
    min_disp = min([r["score"] for r in displaced], default=None)
    return dict(n_pure=len(pure), n_diverse=len(sel), mean_pure=round(float(np.mean(ps)), 4) if ps else None,
                mean_diverse=round(float(np.mean(ss)), 4) if ss else None,
                mean_loss=round(float(np.mean(ps[:k]) - np.mean(ss[:k])), 4) if k else None,
                n_displaced=len(displaced), n_new=len(new),
                n_new_worse_than_every_displaced=sum(1 for r in new if min_disp is not None and r["score"] < min_disp),
                positions_worse=sum(1 for x in loss if x > 1e-9), loss_max=round(max(loss), 4) if loss else None,
                loss_by_position=[round(x, 4) for x in loss],
                archetypes=dict(pure=len({r["arch"] for r in pure}), diverse=len({r["arch"] for r in sel})),
                niches=dict(pure=len({r["niche"] for r in pure}), diverse=len({r["niche"] for r in sel})),
                pairwise=dict(pure=pairwise(pure), diverse=pairwise(sel)),
                note="«сумма прячет потерю»: кроме средних - сколько мест по порядку хуже (positions_worse) и наибольшая потеря места; "
                     "каждый новый вариант (не из чистого топа) хуже любого вытесненного - это видно по n_new_worse_than_every_displaced")


def convergence(recs, n=TOP_N):
    """отбор на первых M/4, M/2, M вариантах (порядок архетипов и ниш по кругу: первые M/2 - тоже отбор с разнообразием)"""
    M = len(recs)
    full = {r["j"] for r in greedy(recs, n)}
    out, prev = [], None
    for m in sorted({max(1, M // 4), max(1, M // 2), M}):
        sub = [r for r in recs if r["j"] < m]
        sel = greedy(sub, n)
        s = {r["j"] for r in sel}
        row = dict(M=m, eligible=len(eligible(sub)), selected=len(sel), overlap_with_full=len(s & full),
                   mean_score=round(float(np.mean([r["score"] for r in sel])), 4) if sel else None,
                   archetypes=len({r["arch"] for r in sel}), niches=len({r["niche"] for r in sel}))
        if prev is not None:
            row["new_vs_half"] = len(s - prev)
        out.append(row)
        prev = s
    return out


def sensitivity_lambda(recs):
    out = []
    for lam in (0.5, 1.0, 2.0):
        sel = greedy(recs, lam=lam)
        p = price_of_diversity(recs, sel)
        out.append(dict(lam=lam, mean_diverse=p["mean_diverse"], n_new=p["n_new"], archetypes=p["archetypes"]["diverse"],
                        niches=p["niches"]["diverse"], min_pairwise=p["pairwise"]["diverse"]["min"]))
    return out


# =========================================================================== пересчёт оценки без нового счёта (v3)
RESCORE_NEED = ("s_parts", "s_plugplay", "s_service", "s_effmake", "trl_min")
RESCORE_VOID = ("d_score_vs_l0", "start_score", "opt_gain")


def rescore(recs, w=None):
    """итоговая оценка по записям уровня 1 без нового счёта: score = funnel_l0.final_score(частные баллы записи, TRL_min записи) -
    та же функция, что у счёта уровней 0 и 1. Записи без частных баллов - score None (в отбор не идут). Поля, которые по записи не
    восстановить (оценка старта, прирост оптимизации, разница с уровнем 0 - по формуле счёта), обнуляются; прежняя оценка -
    score_prev. Возврат: сводка (сколько пересчитано, у скольких допустимых не хватило баллов, изменение оценки)"""
    n_ok, miss, miss_feas, d = 0, 0, 0, []
    for r in recs:
        r["score_prev"] = r.get("score")
        if all(r.get(k) is not None for k in RESCORE_NEED):
            s = {"parts": r["s_parts"], "plugplay": r["s_plugplay"], "service": r["s_service"], "effmake": r["s_effmake"]}
            r["score"] = round(float(L0.final_score(s, r["trl_min"], w)), 4)
            r["score_w"] = round(float(L0.weighted(s, w)), 4)
            r["k_trl"] = round(float(L0.trl_factor(r["trl_min"])), 4)
            n_ok += 1
            if r["score_prev"] is not None:
                d.append(r["score"] - r["score_prev"])
        else:
            miss += 1
            miss_feas += int(bool(r.get("feasible")))
            r["score"] = None
        for k in RESCORE_VOID:
            r[k] = None
    return dict(rule="score = funnel_l0.final_score(s_parts, s_plugplay, s_service, s_effmake; trl_min) - взвешенная × (0,5 + 0,5·TRL/9)",
                n=len(recs), rescored=n_ok, missing_scores=miss, missing_scores_feasible=miss_feas, voided=list(RESCORE_VOID),
                d_score=dict(n=len(d), p50=round(float(np.median(d)), 4) if d else None, min=round(float(min(d)), 4) if d else None,
                             max=round(float(max(d)), 4) if d else None),
                note="частные баллы записи - по моделям того кода, которым считался уровень 1 (code_sha частей); пересчитывается только "
                     "итоговая оценка и отбор, масштабирование и CalculiX - нет (новый счёт)")


# =========================================================================== описание варианта
def scheme(r):
    st = {"1Ц": "одно центробежное колесо", "2Ц": "два центробежных колеса", "Д": "диагональное колесо",
          "ОЦ": f"{r.get('n_ax')} осевых + центробежное"}[r["comp"]]
    tb = ("радиальная" if r["turb"] == "Р" else "осевая") + f" турбина, {r['nt']} ст."
    parts = [f"{r['rotors']} вал(а)" + (f", ядер {r['cores']}" if r["cores"] > 1 else ""), st + (", двусторонний вход" if r["inlet"] == "двух" else ""),
             tb, f"охлаждение {r['cool']}", f"камера {r['comb']}", "магнитные подшипники" if r["bear"] == "магнит" else "масляные подшипники",
             f"ротор {r['rotor_join']} (колесо компрессора: {r['m_cw']} {r['mat_cw']}; турбины: {r['m_tw']} {r['mat_tw']})"]
    if r["rotors"] == 2:
        parts.insert(3, f"силовая турбина {r['npt']} ст.")
    if r["recup"] == "печатный":
        parts.append(f"печатный рекуператор{' в корпусе' if r['recup_integral'] else ''}, ε {float(r['eps_rec']):.2f}")
    return "; ".join(parts)


def entry(r, det, sel, k, Fsel):
    """описание одного из 50 для эскиза и решения"""
    j = r["j"]
    others = [s for s in sel if s["j"] != j]
    nn = min(others, key=lambda s: distance(Fsel[j], Fsel[s["j"]])) if others else None
    e = dict(div_rank=k + 1, j=j, idx=r["idx"], arch=r["arch"], archetype=FS.ARCH[r["arch"]]["name"], niche=r["niche"], scheme=scheme(r),
             score=r["score"], adj_score=round(r["_adj"], 4), penalty=round(r["_pen"], 4), l0_score=r["l0_score"],
             scores=dict(parts=r["s_parts"], plugplay=r["s_plugplay"], service=r["s_service"], make=r["s_make"], eta=r["s_eta"],
                         size=r["s_size"], effmake=r["s_effmake"]),
             eta_th_shaft_pct=r["eta_th_shaft"], eta_el_pct=r["eta_el"], eta_th_50pct_load_iso=r["eta_th_50_iso"],
             parts_and_units=dict(parts=r["n_parts"], assemblies=r["n_assemblies"]), size_m=dict(L=r["L_env"], W=r["W_env"], H=r["H_env"]),
             V_env_m3=r["V_env"], mass_kg_estimate=r["mass_kg"], trl_min=r["trl_min"], SM_min=r["SM_min"], SM_min_at=r["SM_min_at"],
             n_self_sustain=r["n_self_sustain"], P_starter_kW_estimate=r["P_starter_kW"], burst_min_120=r["burst_min_120"],
             pi=r["pi"], t41_C=r["t41"], n_rpm=r["n_rpm"], G=r["G"], dG_pct=r["dG_pct"], EGT=r["EGT"], dEGT=r["dEGT"])
    if nn:
        e["nearest"] = dict(j=nn["j"], idx=nn["idx"], arch=nn["arch"], distance=round(distance(Fsel[j], Fsel[nn["j"]]), 4),
                            differs=diff_axes(nn, r))
    if det:
        e.update(modules=det.get("modules"), assembly_steps=det.get("service", {}).get("steps"), tools=det.get("service", {}).get("tools"),
                 plugplay=det.get("plugplay"), technologies=det.get("make", {}).get("techs"), trl_by=det.get("make", {}).get("trl_by"),
                 geometry=det.get("geometry"), meanline=det.get("meanline"), strength=det.get("strength"),
                 recuperator=det.get("recuperator"), offdesign_summary=_od_sum(det.get("offdesign")), start=det.get("start"),
                 concept=det.get("concept"), optimization=dict(runs=det.get("optimization", {}).get("runs"),
                                                               evals=det.get("optimization", {}).get("evals")))
    return e


def _od_sum(od):
    if not od:
        return None
    return dict(SM_min=od.get("SM_min"), SM_req=od.get("SM_req"), P_avail_hot_kW=od.get("P_avail_hot_kW"), eta_th_50_iso=od.get("eta_th_50_iso"),
                n_points=od.get("n_points"), n_unsolved=od.get("n_unsolved"), unsolved_note=od.get("unsolved_note"),
                points=[{k: p.get(k) for k in ("Ta_C", "load", "law", "T41_C", "P_kW", "eta_th", "SM", "n", "n_pt", "EGT_C", "derated")}
                        for p in od.get("points", [])], fails=od.get("fails"))


# =========================================================================== масштабирование и CalculiX
def scale(sel, det, kv_list, n0, seed, starts, maxfev, workers):
    jobs = [(r["j"], det[r["j"]]["concept"], kv, n0, seed, starts, maxfev) for r in sel if r["j"] in det for kv in kv_list]
    if workers > 1 and len(jobs) > 1:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor
        L1._prep(n0, seed, "KV-3")
        ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else None
        with ProcessPoolExecutor(workers, mp_context=ctx) as ex:
            got = list(ex.map(L1.scale_variant, jobs))
    else:
        got = [L1.scale_variant(jb) for jb in jobs]
    out = defaultdict(dict)
    for (j, _, kv, *_), (row, d) in zip(jobs, got):
        out[j][kv] = dict(feasible=bool(row.get("feasible")), score=row.get("score"), eta_th_shaft=row.get("eta_th_shaft"),
                          reject=row.get("reject"), reason=(row.get("reject_detail") or None), pi=row.get("pi"), t41=row.get("t41"),
                          n_rpm=row.get("n_rpm"), D2=row.get("D2"), SM_min=row.get("SM_min"), V_env=row.get("V_env"),
                          mass_kg=row.get("mass_kg"), t_s=row.get("t_s"))
    return dict(out)


def fea(sel, det, workers, limit=None):
    """CalculiX (disk.py) для отобранных: колёса компрессора (Ti, Ni) и радиальной турбины (не керамика) - условный профиль
    примера в масштабе D (strength_l1.disk_case)"""
    jobs = []
    for r in sel[:limit]:
        d = det.get(r["j"])
        if not d:
            continue
        for w, stg in zip(d["strength"]["comp"], d["meanline"]["compressor"]):
            if w.get("material_data") in ("ti64", "in713c"):
                solid = d["concept"]["rotor_join"] in ("заодно", "сварка")
                case = ST.disk_case(f"l1_j{r['j']}_c{w['stage']}", "comp", w["material_data"], r["n_rpm"], stg["D2"], _t_in(w), w["T_out_C"], solid)
                jobs.append((r["j"], f"comp{w['stage']}", case))
        t = d["strength"].get("turb")
        if t and d["concept"]["turb"] == "Р" and d["concept"]["mat_tw"] != "Si3N4":
            D4 = d["meanline"]["turbine"]["D4"]
            case = ST.disk_case(f"l1_j{r['j']}_t", "rturb", "in713c", r["n_rpm"], D4, t["T_rim_C"] - 240.0, t["T_rim_C"], True)
            jobs.append((r["j"], "turb", case))
    tmp = tempfile.mkdtemp(prefix="l1_fea_")

    def one(job):
        j, part, case = job
        try:
            return j, part, ST.run_disk_fea(case, os.path.join(tmp, case["name"]))
        except Exception as e:                                         # noqa: BLE001 - отказ ccx записывается, отбор не рушит
            return j, part, dict(error=f"{type(e).__name__}: {e}")
    try:
        if workers > 1 and len(jobs) > 1:
            from concurrent.futures import ThreadPoolExecutor
            with ThreadPoolExecutor(workers) as ex:
                got = list(ex.map(one, jobs))
        else:
            got = [one(j) for j in jobs]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    out = defaultdict(dict)
    for j, part, res in got:
        out[j][part] = res
    return dict(out), len(jobs)


FEA_AWS_TOP = 300


def _t_in(w):
    """температура на входе колеса для поля температур disk.py: вход ступени (до 15.11.2026 у всех ступеней стояло 40 °C - у 2-й
    ступени это занижало температуру ступицы; замечание скептика прочности)"""
    if w.get("T_in_C") is not None:
        return float(w["T_in_C"])
    if w.get("T_mean_C") is not None and w.get("T_out_C") is not None:
        return 2.0 * float(w["T_mean_C"]) - float(w["T_out_C"])
    return 40.0


def fea_cases(recs, det, top=FEA_AWS_TOP, sel=()):
    """вход задачи AWS jobs_l1_fea.json: случаи disk.py для отобранных и далее лучших допустимых (до top вариантов)"""
    order = list(sel) + [r for r in sorted((r for r in recs if r["feasible"]), key=lambda r: (-r["score"], r["j"]))
                         if r["j"] not in {s["j"] for s in sel}]
    cases = []
    for r in order[:top]:
        d = det.get(r["j"])
        if not d:
            continue
        for w, stg in zip(d["strength"]["comp"], d["meanline"]["compressor"]):
            if w.get("material_data") in ("ti64", "in713c"):
                solid = d["concept"]["rotor_join"] in ("заодно", "сварка")
                cases.append(dict(key=f"j{r['j']}:comp{w['stage']}", case=ST.disk_case(f"l1_j{r['j']}_c{w['stage']}", "comp", w["material_data"],
                                                                                      r["n_rpm"], stg["D2"], _t_in(w), w["T_out_C"], solid)))
        t = d["strength"].get("turb")
        if t and d["concept"]["turb"] == "Р" and d["concept"]["mat_tw"] != "Si3N4":
            cases.append(dict(key=f"j{r['j']}:turb", case=ST.disk_case(f"l1_j{r['j']}_t", "rturb", "in713c", r["n_rpm"],
                                                                         d["meanline"]["turbine"]["D4"], t["T_rim_C"] - 240.0, t["T_rim_C"], True)))
    return cases


# =========================================================================== запись
def write_csv(path, cols, rows, extra_cols=(), extra=None):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(list(cols) + list(extra_cols))
        for i, r in enumerate(rows):
            w.writerow(["" if v is None else v for v in r] + (["" if v is None else v for v in extra[i]] if extra else []))


def summary(recs):
    rej = Counter()
    for r in recs:
        for code in (r["reject"] or "").split("|"):
            if code:
                rej[re.sub(r"[_\d]+$", "", code)] += 1
    by = {}
    for a in FS.ARCH_KEYS:
        rs = [r for r in recs if r["arch"] == a]
        ok = [r for r in rs if r["feasible"]]
        d = dict(n=len(rs), feasible=len(ok))
        oks = [r for r in ok if r.get("score") is not None]          # пересчёт: у записи без частных баллов оценки нет
        if oks:
            b = max(oks, key=lambda r: (r["score"], -r["j"]))
            et = sorted(r["eta_th_shaft"] for r in oks if r.get("eta_th_shaft") is not None)
            d.update(best_j=b["j"], best_idx=b["idx"], best_score=b["score"], eta_th=et[len(et) // 2] if et else None)
        elif rs:
            d["first_reject"] = dict(Counter((r["reject"] or "").split("|")[0] for r in rs).most_common(3))
        by[a] = d
    feas = [r for r in recs if r["feasible"]]
    gains = [r["opt_gain"] for r in recs if r.get("opt_gain") is not None and r.get("start_feasible")]
    return dict(n=len(recs), feasible=len(feas), eligible=len(eligible(recs)), score_min=SCORE_MIN,
                reject_reasons=dict(rej.most_common()), by_arch=by,
                l0_feasible_l1_infeasible=sum(1 for r in recs if r["l0_feasible"] and not r["feasible"]),
                l0_infeasible_l1_feasible=sum(1 for r in recs if not r["l0_feasible"] and r["feasible"]),
                score_l1=dict(p50=(float(np.median(sv)) if (sv := [r["score"] for r in feas if r.get("score") is not None]) else None),
                              max=max([r["score"] for r in feas if r.get("score") is not None], default=None)),
                d_score_vs_l0=dict(p50=(float(np.median(dv)) if (dv := [r["d_score_vs_l0"] for r in feas if r.get("d_score_vs_l0") is not None])
                                        else None)),
                optimization=dict(start_feasible=sum(1 for r in recs if r.get("start_feasible")), n_gain_measured=len(gains),
                                  gain_p50=float(np.median(gains)) if gains else None,
                                  gain_negative=sum(1 for g in gains if g < -1e-9),
                                  note="opt_gain - оценка итога минус оценка точки уровня 0 в моделях уровня 1 (только если та точка "
                                       "допустима в расчётной точке; итог может быть ниже, если точка уровня 0 не прошла нерасчётные режимы)"))


def inputs(ref):
    files = ["funnel/funnel_l1.py", "funnel/merge_l1.py", "funnel/meanline_l1.py", "funnel/strength_l1.py", "funnel/offdesign_l1.py",
             "funnel/funnel_l0.py", "funnel/funnel_space.py", "frame3/frame3_point.py", "frame3/frame3_fea_results.json",
             "platform/fea/materials.json", "platform/fea/disk_results.json", "platform/fea/disk.py", "gtu_explore.py", "gtu_family_results.json"]
    out = [dict(file="routes/gtu/" + f, sha256=hashlib.sha256(open(os.path.join(GTU, f), "rb").read()).hexdigest()) for f in files]
    out.append(dict(file="отбор уровня 0", field=f"{ref['input']} sha256 {ref['input_sha256']}, отпечаток выборки {ref.get('fingerprint')}, "
                                                 f"код уровня 0 {ref.get('l0_code_sha')}"))
    return out


FEA_BURST_MIN = 1.0           # запас по разрушающей частоте при 120 % у disk.py (как правило уровня 1 для быстрой прочности)


def apply_fea_results(recs, path, burst_min=FEA_BURST_MIN):
    """итог задачи AWS jobs_l1_fea.json (strength_l1.run_cases -> funnel_l1_fea_results.json) до отбора 50: вариант, у которого
    CalculiX дал запас по разрушению при 120 % ниже burst_min хоть у одного колеса, становится недопустимым с причиной fea_burst
    (быстрая прочность уровня 1 - подобие disk.py; CalculiX - проверка). Возврат: результаты по j и число отсеянных"""
    doc = json.load(open(path, encoding="utf-8"))
    by = defaultdict(dict)
    for it in doc.get("results", []):
        jk, part = it["key"].split(":", 1)
        by[int(jk.lstrip("j"))][part] = it.get("result") or dict(error=it.get("error"))
    n_bad = 0
    for r in recs:
        res = by.get(r["j"])
        if not res or not r["feasible"]:
            continue
        bad = [f"{p_}: {x['burst_margin_120']:.2f}" for p_, x in sorted(res.items())
               if x.get("burst_margin_120") is not None and x["burst_margin_120"] < burst_min]
        if bad:
            r["feasible"] = False
            r["reject"] = "|".join(x for x in ((r.get("reject") or ""), "fea_burst") if x)
            r["reject_detail"] = "; ".join(x for x in ((r.get("reject_detail") or ""), f"CalculiX disk.py: запас по разрушению при 120 % < "
                                                       f"{burst_min} ({', '.join(bad)})") if x)
            n_bad += 1
    return dict(by), n_bad


def run_merge(found, kv, M, shards, out_dir, sfx="", scale_kv=("KV-1", "KV-6"), scale_starts=None, scale_maxfev=None, workers=1,
              fea_top=0, compare=None, extra_results=None, fea_cases_out=None, fea_results=None, rescore_mode=False):
    t0 = time.time()
    cols, rows, det, ref, checks, calls = merge_parts(found, kv, M, shards)
    recs = [dict(zip(cols, r)) for r in rows]
    fea_aws, n_fea_bad = ({}, 0)
    if fea_results:
        fea_aws, n_fea_bad = apply_fea_results(recs, fea_results)
        rows = [[r.get(c) for c in cols] for r in recs]           # CSV - с новой допустимостью и причиной
    rs = None
    if rescore_mode:                                              # пересчёт без нового счёта: оценка - по частным баллам записей
        rs = rescore(recs)
        cols = list(cols) + [c for c in ("score_prev", "score_w", "k_trl") if c not in cols]
        rows = [[r.get(c) for c in cols] for r in recs]
        scale_kv, fea_top, fea_cases_out, compare = (), 0, None, None      # новый счёт и сравнение с другими частями - не здесь
    sel = greedy(recs)
    pod = price_of_diversity(recs, sel)
    Fsel = {r["j"]: features(r) for r in sel}
    pure = sorted(eligible(recs), key=lambda r: (-r["score"], r["j"]))
    prank = {r["j"]: k + 1 for k, r in enumerate(pure)}
    drank = {r["j"]: k + 1 for k, r in enumerate(sel)}
    el = {r["j"] for r in eligible(recs)}
    extra = []
    for r in recs:
        why = None
        if r["j"] in drank:
            why = f"отобран: место {drank[r['j']]} из {len(sel)}"
        elif r["j"] in el:
            why = "допустим и выше порога, вытеснен штрафом за близость или лимитом ниши"
        elif r["feasible"]:
            why = f"допустим, оценка ниже порога {SCORE_MIN}"
        else:
            why = "отсеян на уровне 1: " + (r["reject_detail"] or r["reject"] or "")
        extra.append([drank.get(r["j"]), prank.get(r["j"]), int(r["j"] in el), why])
    scl, t_scale = {}, None
    if scale_kv and sel:
        ts = time.time()
        scl = scale(sel, det, scale_kv, ref["n0"], ref["seed"], scale_starts or ref["starts"], scale_maxfev or ref["maxfev"], workers)
        t_scale = round(time.time() - ts, 1)
    feares, n_fea, t_fea = {j: v for j, v in fea_aws.items() if j in {r["j"] for r in sel}}, 0, None
    if fea_top and sel:
        tf = time.time()
        feares, n_fea = fea(sel, det, workers, fea_top)
        t_fea = round(time.time() - tf, 1)
    top = []
    for k, r in enumerate(sel):
        e = entry(r, det.get(r["j"]), sel, k, Fsel)
        e["scaling"] = scl.get(r["j"])
        if r["j"] in feares:
            e["calculix"] = feares[r["j"]]
        top.append(e)
    os.makedirs(out_dir, exist_ok=True)
    p_csv = os.path.join(out_dir, f"funnel_l1{sfx}.csv")
    write_csv(p_csv, cols, rows, EXTRA, extra)
    top_cols = ["div_rank", "j", "idx", "arch", "niche", "scheme", "score", "adj_score", "penalty", "l0_score", "eta_th_shaft_pct", "eta_el_pct",
                "eta_th_50pct_load_iso", "parts", "assemblies", "n_modules", "assembly_steps", "L", "W", "H", "V_env_m3", "mass_kg_estimate",
                "trl_min", "SM_min", "n_self_sustain", "P_starter_kW_estimate", "burst_min_120", "pi", "t41_C", "n_rpm", "dG_pct", "dEGT",
                "nearest_j", "nearest_distance", "nearest_differs", "KV-1", "KV-6"]
    p_tcsv = os.path.join(out_dir, f"funnel_l1_top50{sfx}.csv")
    with open(p_tcsv, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(top_cols)
        for e in top:
            sc = e.get("scaling") or {}
            fmt = lambda s: ("" if not s else ("масштабируется: оценка %.3f, КПД %.2f %%" % (s["score"], s["eta_th_shaft"]) if s["feasible"]
                                                else "нет: " + (s.get("reason") or s.get("reject") or "")))
            w.writerow([e["div_rank"], e["j"], e["idx"], e["arch"], e["niche"], e["scheme"], e["score"], e["adj_score"], e["penalty"], e["l0_score"],
                        e["eta_th_shaft_pct"], e["eta_el_pct"], e["eta_th_50pct_load_iso"], e["parts_and_units"]["parts"],
                        e["parts_and_units"]["assemblies"], len(e.get("modules") or []), e.get("assembly_steps"), e["size_m"]["L"], e["size_m"]["W"],
                        e["size_m"]["H"], e["V_env_m3"], e["mass_kg_estimate"], e["trl_min"], e["SM_min"], e["n_self_sustain"],
                        e["P_starter_kW_estimate"], e["burst_min_120"], e["pi"], e["t41_C"], e["n_rpm"], e["dG_pct"], e["dEGT"],
                        (e.get("nearest") or {}).get("j"), (e.get("nearest") or {}).get("distance"),
                        "; ".join((e.get("nearest") or {}).get("differs") or []), fmt(sc.get("KV-1")), fmt(sc.get("KV-6"))])
    soft = dict(method=ref["method"], software=ref["software"], date=DATE)
    cal = ML.calibration_summary()
    n_sel = len(sel)
    sel_note = (f"отобрано {n_sel} из {TOP_N}" + ("" if n_sel >= TOP_N else
                f": выше порога {SCORE_MIN} и допустимых на уровне 1 - {len(eligible(recs))} (с лимитом {NICHE_MAX} на нишу); "
                "остальные варианты отсеяны правилами уровня 1 (причины - funnel_l1.csv, reject_detail) или ниже порога"))
    p_top = os.path.join(out_dir, f"funnel_l1_top50{sfx}.json")
    top_doc = dict(title=f"Воронка, уровень 1: {n_sel} разных лучших вариантов - {kv}", **soft, kv=kv, M=M, input=ref["input"],
                   input_sha256=ref["input_sha256"], code_sha=ref["code_sha"], starts=ref["starts"], maxfev=ref["maxfev"],
                   rule=dict(top_n=TOP_N, niche_max=NICHE_MAX, score_min=SCORE_MIN, score_min_basis=SCORE_MIN_BASIS, lam=LAMBDA, d0=D0,
                             basis=DIV_BASIS, distance=dict(disc_weights=DISC_W, cont=CONT_N, cont_weight=CONT_W)),
                   selection_note=sel_note, price_of_diversity=pod, selected=top,
                   note="все числа - оценка уровня 1 (одномерные модели, привязанные к формулам трека; README узла, раздел «Уровень 1»)")
    with open(p_top, "w", encoding="utf-8") as fh:
        json.dump(L1._clean(top_doc), fh, ensure_ascii=False, indent=1)
    if fea_cases_out:
        fc = fea_cases(recs, det, FEA_AWS_TOP, sel)
        with open(fea_cases_out, "w", encoding="utf-8") as fh:
            json.dump(dict(title="уровень 1: случаи disk.py для лучших (задача AWS jobs_l1_fea.json)", date=DATE, kv=kv, M=M,
                           code_sha=ref["code_sha"], top=FEA_AWS_TOP, cases=fc), fh, ensure_ascii=False)
    comp = None
    if compare:
        comp = compare_starts(recs, compare)
    tm = [r["t_s"] for r in recs if r.get("t_s") is not None]
    to = [r["t_opt_s"] for r in recs if r.get("t_opt_s") is not None]
    tcpu = [r["t_cpu_s"] for r in recs if r.get("t_cpu_s") is not None]
    res = dict(title=f"Воронка двигателя-заменителя, уровень 1: углублённый расчёт отбора уровня 0 - {kv}, M = {M}", **soft,
               status=(extra_results or {}).pop("status", "слияние частей уровня 1"), kv=kv, target=ref["target"], M=M,
               input=dict(file=ref["input"], sha256=ref["input_sha256"], fingerprint=ref.get("fingerprint"), l0_code_sha=ref.get("l0_code_sha"),
                          order_sha=ref["order_sha"], problems=ref.get("input_problems")),
               inputs=inputs(ref), code_sha=ref["code_sha"], starts=ref["starts"], maxfev=ref["maxfev"], de_gens=ref.get("de_gens"),
               de_pop=ref.get("de_pop"), merge_checks=checks,
               calibration=cal, strength_facts=ST.facts(),
               rules_l1=dict(M1T_MAX=L1.M1T_MAX, DH_MIN=L1.DH_MIN, D1D2_MAX=L1.D1D2_MAX, ALPHA2_MAX=L1.ALPHA2_MAX, BETA2_MAX=L1.BETA2_MAX,
                             M2T_MAX=L1.M2T_MAX, MW3_MAX=L1.MW3_MAX, R_MIN=L1.R_MIN, M4_MAX=L1.M4_MAX, MW5_MAX=L1.MW5_MAX,
                             T_TI_MAX=L1.T_TI_MAX, SM_MIN=OD.SM_MIN, T_AMB=OD.T_AMB, LOADS=OD.LOADS, K_NU=OD.K_NU,
                             T41_START_C=OD.T41_START_C, recuperator=L1.REC, penalty=L1.PENALTY,
                             HUB_M2_MAX=ML.HUB_M2_MAX, R_HUB_MIN=ML.R_HUB_MIN, AX_MODEL=ML.AX_MODEL, T_TE=ML.T_TE, TMAX_C=ML.TMAX_C,
                             STANITZ_K=ML.STANITZ_K, CC_CAP="КПД ступени не выше eta_cap_l0 (формула уровня 0 + TOL_CC)", TOL_CC=ML.TOL_CC,
                             OD_ISO_DT41=L1.OD_ISO_DT41, OD_ISO_DP=L1.OD_ISO_DP, OD_COVER="все режимы T_AMB × LOADS решены",
                             DINC_STALL=OD.DINC_STALL, TIP_INC_STALL=OD.TIP_INC_STALL, ALPHA2_STALL=OD.ALPHA2_STALL, A1_STEPS=OD.A1_STEPS,
                             T1H_STEPS=OD.T1H_STEPS, CHOKE_X_MIN=OD.CHOKE_X_MIN, VIGV_Z=L1.VIGV_Z, VIGV_ACT_KG=L1.VIGV_ACT_KG,
                             fixes_15_11=("исправления скептиков уровня 1 (README, раздел «Уровень 1»): счёт всех нерешённых режимов и od_cover, "
                                          "строгий od_iso, ВНА / ветвь запирания / подогрев входа с учётом деталей, граница помпажа по "
                                          "признакам срыва, поправки ступени внутри модели, ограничение КПД компрессора, Kacker - Okapuu и "
                                          "свободный вихрь, число ступеней турбины - ось, Stanitz и цепочка радиальных ступеней, повтор "
                                          "части при слиянии, температура входа колеса в CalculiX, смета по процессорному времени")),
               timing=dict(calls=calls, s_per_variant=dict(mean=round(float(np.mean(tm)), 3) if tm else None,
                                                           p50=round(float(np.median(tm)), 3) if tm else None,
                                                           max=round(float(max(tm)), 3) if tm else None, n=len(tm)),
                           cpu_s_per_variant=dict(mean=round(float(np.mean(tcpu)), 3) if tcpu else None,
                                                  p50=round(float(np.median(tcpu)), 3) if tcpu else None,
                                                  max=round(float(max(tcpu)), 3) if tcpu else None, n=len(tcpu),
                                                  note="процессорное время одного процесса на вариант - основа сметы (настенное на "
                                                       "общем контейнере завышено чужой нагрузкой)"),
                           note="время на вариант - в одном процессе, с оптимизацией (starts × maxfev) и нерасчётными режимами",
                           scale_s=t_scale, fea_s=t_fea, merge_s=None),
               timing_detail=dict(t_opt_mean=round(float(np.mean(to)), 3) if to else None,
                                  t_od_mean=round(float(np.mean(tm) - np.mean(to)), 3) if to and tm else None,
                                  opt_share=round(float(np.mean(to) / np.mean(tm)), 3) if to and tm else None,
                                  evals_per_start=(round(float(np.mean(ev)), 1) if (ev := [r["evals"] / r["starts"] for r in recs if r.get("evals")])
                                                   else None),
                                  note="t_opt - оптимизация (starts × maxfev); остальное - проверка кандидатов с нерасчётными режимами и запись"),
               summary=summary(recs),
               selection=dict(file=os.path.basename(p_top), n=n_sel, note=sel_note, archetypes=dict(Counter(r["arch"] for r in sel).most_common()),
                              niches=len({r["niche"] for r in sel}), price_of_diversity=pod, lambda_sensitivity=sensitivity_lambda(recs)),
               convergence=dict(rule="отбор 50 на первых M/4, M/2 и M вариантах в порядке «архетипы по кругу, внутри - ниши по кругу»; "
                                     "сравнение с отбором на всех",
                                rows=convergence(recs)),
               starts_check=comp,
               scaling=_scale_summary(scl), calculix=dict(cases=n_fea, top=fea_top,
                                                         results={str(j): v for j, v in feares.items()} if feares else None,
                                                         aws=(dict(file=os.path.basename(fea_results), variants=len(fea_aws), rejected=n_fea_bad,
                                                                   burst_min=FEA_BURST_MIN) if fea_results else None)),
               files=dict(csv=os.path.basename(p_csv), top50=os.path.basename(p_top), top50_csv=os.path.basename(p_tcsv)))
    if rs:
        res["rescore"] = rs
        res["status"] = "пересчёт отбора по записям уровня 1 без нового счёта (v3: множитель готовности TRL)"
        top_doc["rescore"] = rs
        with open(p_top, "w", encoding="utf-8") as fh:
            json.dump(L1._clean(top_doc), fh, ensure_ascii=False, indent=1)
    if extra_results:
        res.update(extra_results)
    res["timing"]["merge_s"] = round(time.time() - t0, 1)
    p_res = os.path.join(out_dir, f"funnel_l1_results{sfx}.json")
    with open(p_res, "w", encoding="utf-8") as fh:
        json.dump(L1._clean(res), fh, ensure_ascii=False, indent=1)
    return dict(results=p_res, csv=p_csv, top50=p_top, top50_csv=p_tcsv, res=res, sel=sel, recs=recs)


def _scale_summary(scl):
    if not scl:
        return None
    out = {}
    for kv in ("KV-1", "KV-6"):
        v = [d[kv] for d in scl.values() if kv in d]
        if v:
            out[kv] = dict(n=len(v), feasible=sum(1 for x in v if x["feasible"]),
                           reasons=dict(Counter(re.sub(r"[_\d]+$", "", c) for x in v if not x["feasible"]
                                                for c in (x.get("reject") or "").split("|") if c).most_common(8)))
    return out


def compare_starts(recs, other_dirs):
    """те же варианты с удвоенным числом стартов (задача AWS l1-starts2x): изменение оценки и отбора"""
    found = []
    for d in other_dirs:
        found += glob.glob(os.path.join(d, "**", f"{L1.PART_PREFIX}_starts2x_*.json*"), recursive=True)
    rows = {}
    for p in found:
        dd = load_part(p)
        for r in dd["rows"]:
            rec = dict(zip(dd["columns"], r))
            rows[rec["idx"]] = rec
    if not rows:
        return dict(note="частей с удвоенным числом стартов нет")
    base = {r["idx"]: r for r in recs}
    pairs = [(base[i], rows[i]) for i in rows if i in base]
    d = [(b2.get("score") or 0) - (b.get("score") or 0) for b, b2 in pairs]
    return dict(n=len(pairs), mean_gain=round(float(np.mean(d)), 4), improved_gt_0_05=sum(1 for x in d if x > 0.05),
                worse_gt_0_05=sum(1 for x in d if x < -0.05), feasible_base=sum(1 for b, _ in pairs if b["feasible"]),
                feasible_2x=sum(1 for _, b2 in pairs if b2["feasible"]))


# =========================================================================== короткая локальная проба
PROBE = dict(input=os.path.join(HERE, "l0_aws_debug", "funnel_l0_select_KV-3_n2000.json"), top=12, shards=4, starts=2, maxfev=150, workers=4)


def probe(out_dir=HERE, **kw):
    P = dict(PROBE, **kw)
    if P["top"] > L1.LOCAL_MAX_VARIANTS or P["starts"] > L1.LOCAL_MAX_STARTS:
        raise ValueError("локально - не больше 12 вариантов и 2 стартов (решение владельца 14.11.2026)")
    tmp = tempfile.mkdtemp(prefix="funnel_l1_probe_")
    try:
        t0 = time.time()
        paths, meta = L1.run_parts(P["input"], P["top"], list(range(P["shards"])), P["shards"], tmp, P["starts"], P["maxfev"], P["workers"], gz=True)
        t_par = time.time() - t0
        found = find_parts([tmp], meta["kv"], meta["M"], P["shards"])
        extra = dict(status=(f"короткая локальная проба (4 общих ядра): {meta['kv']}, первые {meta['M']} вариантов отбора уровня 0 AWS N = 2000 "
                             f"в порядке архетипов и ниш по кругу, {P['starts']} старта (дифференциальная эволюция, {L1.DE_GENS} поколений + "
                             f"доводка COBYLA {P['maxfev']} вычислений); отладочный (около 150) и полный "
                             "(3000-5000) прогоны - задачи AWS jobs_l1_debug.json и jobs_l1.json"),
                     probe=dict(params={k: (os.path.relpath(v, HERE) if k == "input" else v) for k, v in P.items()},
                                timing=dict(parts_wall_s=round(t_par, 1), s_per_variant_pool=round(t_par * P["workers"] / max(1, meta["M"]), 2),
                                            note="контейнер общий (4 ядра, нагрузка чужими задачами) - время завышено против машины AWS")))
        r = run_merge(found, meta["kv"], meta["M"], P["shards"], out_dir, "", workers=P["workers"], fea_top=TOP_N,
                      scale_starts=P["starts"], scale_maxfev=P["maxfev"], extra_results=extra)
        return r
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--kv", default="KV-3")
    ap.add_argument("--m", type=int, help="число вариантов уровня 1 (M в именах частей)")
    ap.add_argument("--shards", type=int)
    ap.add_argument("--in", dest="inp", nargs="+", default=["."])
    ap.add_argument("--tag", default="")
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--suffix", default="")
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--no-scale", action="store_true")
    ap.add_argument("--scale-starts", type=int, default=None)
    ap.add_argument("--fea-top", type=int, default=TOP_N, help="CalculiX disk.py для первых отобранных (0 - нет)")
    ap.add_argument("--compare-in", nargs="*", default=None, help="папки частей с удвоенным числом стартов (_starts2x)")
    ap.add_argument("--fea-cases-out", default=None, help=f"файл случаев disk.py для лучших {FEA_AWS_TOP} (вход jobs_l1_fea.json)")
    ap.add_argument("--fea-results", default=None, help="итог jobs_l1_fea.json (funnel_l1_fea_results.json): отсев по CalculiX до отбора 50")
    ap.add_argument("--rescore", action="store_true", help="пересчёт оценки (final_score, множитель TRL) и отбора 50 по записям частей без "
                                                           "нового счёта; масштабирование и CalculiX не запускаются")
    a = ap.parse_args(argv)
    if a.probe:
        r = probe()
    else:
        if not (a.m and a.shards):
            ap.error("нужны --m и --shards")
        try:
            r = run_merge(find_parts(a.inp, a.kv, a.m, a.shards, a.tag), a.kv, a.m, a.shards, a.out_dir, a.suffix,
                          () if (a.no_scale or a.rescore) else ("KV-1", "KV-6"), a.scale_starts, None, a.workers or L0.default_workers(),
                          0 if a.rescore else a.fea_top, a.compare_in, fea_cases_out=a.fea_cases_out, fea_results=a.fea_results,
                          rescore_mode=a.rescore)
        except MergeError as e:
            print(f"слияние отклонено: {e}", file=sys.stderr)
            return 2
    s = r["res"]["summary"]; sl = r["res"]["selection"]
    print(f"{r['res']['kv']} M={r['res']['M']}: допустимых на уровне 1 {s['feasible']}, выше порога {s['eligible']}, отобрано {sl['n']} "
          f"(архетипов {len(sl['archetypes'])}, ниш {sl['niches']}); -> {r['results']}, {r['csv']}, {r['top50']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
