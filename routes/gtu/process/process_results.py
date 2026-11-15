#!/usr/bin/env python3
"""Сводка узла «пространство сквозных концепций» -> process_results.json: самопроверка оценки (опорные концепции против
примеров модели массы), обратимость генома (encode/decode опорных концепций), детерминизм островов, мини-прогон KV-3 и
пробы KV-1, KV-6 (сводки process_elites.py mini), наблюдения по нишам (считаются из архивов, не пишутся руками), что не
сделано. Числа README узла - из этого файла.

  python3 process_results.py --run KV-3=<папка мини-прогона> --run KV-1=<папка пробы> --run KV-6=<папка пробы>
Папка прогона - выход process_elites.py mini (summary.json и process_<KV>_part*.json.gz). Полные сводки и архивы -
вне репозитория (scratchpad), здесь - выжимка.
"""
import argparse
import collections
import datetime
import glob
import gzip
import hashlib
import json
import os
import statistics
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import process_space as PS  # noqa: E402
import process_eval as PE  # noqa: E402
import process_elites as PEL  # noqa: E402
import mass_model as MM  # noqa: E402

OUT = os.path.join(HERE, "process_results.json")
DATE = "2026-11-15"


def _sha(p):
    try:
        return hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]
    except OSError:
        return None


def code_sha():
    h = hashlib.sha256()
    for f in ("process_space.py", "process_eval.py", "process_elites.py"):
        h.update(open(os.path.join(HERE, f), "rb").read())
    return h.hexdigest()[:16]


def selfcheck():
    """опорные концепции (frame3 в четырёх семьях, с рекуператором и подогревом) + сверка с примерами модели массы"""
    rows = PE.selfcheck()
    cal = json.load(open(MM.OUT_FILE, encoding="utf-8"))
    ex = {e["name"]: e for e in cal["examples"]["items"]}
    pairs = [("frame3, традиционная", "KV-3 простой; литьё в песок; масло пакета"),
             ("frame3, литьё КНР", "KV-3 простой; точное литьё КНР; масло пакета")]
    cmp_ = []
    for mine, theirs in pairs:
        a = next((r for r in rows if r["name"] == mine), None)
        b = ex.get(theirs)
        if a and b and a.get("mass_t"):
            cmp_.append(dict(concept=mine, mass_t=a["mass_t"], G_air=a["G_air"], mass_model_example=theirs,
                             example_t=round(b["total_kg"] / 1e3, 4), diff_pct=round(100 * (a["mass_t"] * 1e3 / b["total_kg"] - 1), 2),
                             note="пример модели массы - паспортный путь при расходе frame3 25,55 кг/с; здесь - расход цепочки блоков"))
    return dict(rows=rows, vs_mass_model_examples=cmp_,
                note="frame3 в блоках недопустима по допуску пакета (расход +25 %, как в примерах блоков) - масса посчитана для сверки "
                     "(mass_always); печать - сплав колеса Haynes 282 с пределом класса U 850 °C, поэтому расход больше и масса не сравнима")


def roundtrip():
    """encode -> decode опорных концепций: дискретные выборы и непрерывные величины"""
    from blocks import reference as REF
    from blocks.concept import normalize
    out, bad = [], 0
    for name, x in PS.reference_seeds("KV-3"):
        c0 = REF.example(name)
        Cn = normalize(dict(c0, kv="KV-3"))
        C = PS.decode(x, "KV-3")
        same = dict(arch=C["compressor"]["arch"] == Cn["compressor"]["arch"], turb=C["turbine"]["type"] == Cn["turbine"]["type"],
                    cool=C["turbine"]["cool"] == Cn["turbine"]["cool"], drive=C["shaft"]["drive"] == Cn["shaft"]["drive"],
                    scheme=C["combustor"]["scheme"] == Cn["combustor"]["scheme"], prep=C["fuel"]["prep"] == Cn["fuel"]["prep"],
                    recup=bool(C["hx"]["recup"]) == bool(Cn["hx"]["recup"]), reheat=C["hx"]["reheat"] == bool(Cn["hx"]["reheat"]))
        d = dict(pi=abs(C["compressor"]["pi"] / Cn["compressor"]["pi"] - 1), t41=abs(C["turbine"]["t41"] - Cn["turbine"]["t41"]),
                 eps=abs((C["hx"]["recup"] or 0) - (Cn["hx"]["recup"] or 0)))
        ok = all(same.values())
        bad += not ok
        out.append(dict(name=name, discrete_same=ok, pi_rel=round(d["pi"], 5), t41_abs_K=round(d["t41"], 2), eps_abs=round(d["eps"], 5),
                        repairs=C["_repairs"]))
    return dict(rows=out, all_discrete_same=bad == 0,
                note="T41 вне окна класса (предел минус 200 K .. предел) переносится на ближнюю границу; оболочка без печатного корпуса - "
                     "картридж (repairs)")


def determinism():
    """один остров при 1 и 3 процессах пула - один и тот же архив"""
    import tempfile
    outs = []
    t0 = time.time()
    with tempfile.TemporaryDirectory() as td:
        for w in (1, 3):
            _, doc = PEL.run_part("KV-3", 1, 1, 0, 120, 0.3, 40, PS.SEED, w, os.path.join(td, f"w{w}"), log=None)
            outs.append(json.dumps(doc["archive"], sort_keys=True))
    return dict(evals=120, workers=[1, 3], same_archive=outs[0] == outs[1], wall_s=round(time.time() - t0, 1))


def _load_parts(d, kv):
    arch = PEL.NicheArchive()
    sizes, n_isl = [], 0
    for p in sorted(glob.glob(os.path.join(d, f"process_{kv}_part*.json.gz"))):
        sizes.append(os.path.getsize(p))
        with gzip.open(p, "rt", encoding="utf-8") as f:
            doc = json.load(f)
        n_isl += len(doc["islands"])
        arch.merge(PEL.NicheArchive.from_list(doc["archive"]).cells)
    return arch, sizes, n_isl


def observations(arch):
    """наблюдения по нишам: медианы и лучшие по осям (считаются по элитам архива)"""
    fc = arch.feasible_cells()
    elites = [(c, arch.cells[c]["best"]) for c in fc]
    out = {}
    for j, name in enumerate(PS.DESC_NAMES):
        g = collections.defaultdict(list)
        for c, b in elites:
            g[PS.DESC_LABELS[j][c[j]]].append(b)
        out[name] = {lab: dict(n_niches=len(v), fit_max=round(max(x["fit"] for x in v), 3), fit_median=round(statistics.median(x["fit"] for x in v), 3),
                               mass_median_t=round(statistics.median(x["mass"] for x in v), 3), eff_median=round(statistics.median(x["eff"] for x in v), 2),
                               eff_max=round(max(x["eff"] for x in v), 2), kwkg_median=round(statistics.median(x["kwkg"] for x in v), 3),
                               parts_median=round(statistics.median(x["parts"] for x in v)), trl_min=min(x["trl"] for x in v))
                     for lab, v in sorted(g.items(), key=lambda kv_: -max(x["fit"] for x in kv_[1]))}
    # вал: та же ниша по остальным осям - прямой против свободной СТ и быстроходного (парное сравнение элит)
    pair = collections.defaultdict(list)
    for c, b in elites:
        comp, turb, drive = divmod(c[2], len(PS.TURBS) * len(PS.DRIVES))[0], None, None
        a = c[2]
        comp, rest = divmod(a, len(PS.TURBS) * len(PS.DRIVES))
        turb, drive = divmod(rest, len(PS.DRIVES))
        pair[(c[0], c[1], comp, turb, c[3])].append((PS.DRIVES[drive], b))
    ratios = collections.defaultdict(list)
    for k, v in pair.items():
        d = dict(v)
        if "прямой" in d:
            for other in ("свободная СТ", "быстроходный"):
                if other in d:
                    ratios[other].append(dict(mass=d[other]["mass"] / d["прямой"]["mass"], fit=d[other]["fit"] / d["прямой"]["fit"],
                                              eff=d[other]["eff"] - d["прямой"]["eff"]))
    out["вал_парно"] = {o: dict(pairs=len(v), mass_ratio_median=round(statistics.median(x["mass"] for x in v), 3),
                               fit_ratio_median=round(statistics.median(x["fit"] for x in v), 3),
                               d_eff_median_pp=round(statistics.median(x["eff"] for x in v), 2),
                               worse_fit_share=round(sum(x["fit"] < 1 for x in v) / len(v), 3))
                        for o, v in ratios.items() if v}
    out["вал_парно_правило"] = ("пары элит в одной нише по топологии, камере, компрессору, турбине и изготовлению: «свободная СТ» или "
                                "«быстроходный» против «прямой»; отношение массы и КПД на тонну, разность КПД; worse_fit_share - доля пар, "
                                "где КПД на тонну хуже (сумма прячет потерю)")
    return out


def trim_summary(s, n=12):
    keep = ("kv", "evals", "feasible", "feasible_share", "stages", "rejects_top", "cpu_s_per_eval", "niches", "best", "run_wall_s")
    out = {k: s.get(k) for k in keep}
    out["merge"] = {k: s["merge"].get(k) for k in ("status", "complete", "expected_islands", "present", "missing", "evals_total",
                                                   "evals_expected", "config_hash")}
    out["config"] = s["config"]
    out["top"] = s["top"][:n]
    out["top_constrained"] = dict(rule=s["top_constrained"]["rule"], rows=s["top_constrained"]["rows"][:n])
    out["top_trl"] = dict(rule=s["top_trl"]["rule"], rows=s["top_trl"]["rows"][:n])
    out["pareto_global"] = dict(rule=s["pareto_global"]["rule"], rows=s["pareto_global"]["rows"][:25])
    out["best_by_descriptor"] = s["best_by_descriptor"]
    out["detail_top"] = s["detail_top"][:6]
    out["history"] = s["history"]
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", action="append", default=[], help="KV=папка (выход process_elites.py mini)")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args(argv)
    t0 = time.time()
    runs, obs = {}, {}
    for item in a.run:
        kv, d = item.split("=", 1)
        s = json.load(open(os.path.join(d, "summary.json"), encoding="utf-8"))
        arch, sizes, n_isl = _load_parts(d, kv)
        r = trim_summary(s)
        r["part_bytes_per_island"] = int(sum(sizes) / max(n_isl, 1))
        r["role"] = "мини-прогон" if kv == "KV-3" else "проба"
        runs[kv] = r
        obs[kv] = observations(arch)
    res = dict(
        title="Пространство сквозных концепций «топливо -> вал»: MAP-Elites по КПД на тонну (уровень облика)", date=DATE,
        status="черновик: оценка уровня облика; КПД - цепочка блоков (калибровка трека ±2 п. п.), масса - модель по 8 машинам (LOO 22,8 % "
               "промышленная практика); всё, что сверх них, - оценка с пометкой в строках",
        owner_focus=("15.11.2026: подбираем топовые конструкции всего процесса - от поступления топлива до вала; глубже в расчёты не "
                     "погружаемся; концентрируемся на геометрии, которая даст наибольшую эффективность и значительно повысит КПД на тонну"),
        criterion=dict(kpd_per_t="КПД на клеммах, % / масса агрегата-заменителя, т (двигатель + камера + теплообменники + корпус и рама модуля; "
                                 "генератор, редуктор, вход и выхлоп пакета не входят)",
                       power="мощность на валу - стык пакета (funnel_space.target): KV-3 Taurus 60, KV-1 Saturn 20, KV-6 Titan 130 (генераторные)",
                       constraints="допуск пакета: воздух и выхлоп до +20 %, выхлоп до +100 °C; Plug&Play «невозможно» - отказ; детали и "
                                   "Plug&Play - вторичные цели (отбор с ограничениями)"),
        space=dict(D=PS.D, genes=[dict(id=g["id"], kind=g["kind"], opts=g.get("opts"), what=g["what"]) for g in PS.GENES],
                   grid=dict(names=PS.DESC_NAMES, bins=PS.DESC_BINS, n_cells=PS.N_CELLS, note=PS.DESC_NOTE),
                   make=PS.MAKE, ranges=PS.RANGES, space_hash=PS.space_hash()),
        evaluation=dict(extra_mass_basis=PE.A_EVAL["basis"], hybrid_mass_set=PE.HYBRID_SET, practice=PE.PRACTICE),
        elites=dict(front_cap=PEL.FRONT_CAP, ops=PEL.OPS, sigma=PEL.SIGMA, p_feas_parent=PEL.P_FEAS_PARENT, tournament=PEL.TOURNAMENT,
                    pnp_min_constrained=PEL.PNP_MIN, version=PEL.VERSION),
        selfcheck=selfcheck(), roundtrip=roundtrip(), determinism=determinism(),
        runs=runs, observations=obs,
        inputs=dict(blocks_results=_sha(os.path.join(HERE, "blocks", "blocks_results.json")),
                    mass_calibration=_sha(MM.OUT_FILE), code_sha=code_sha()),
        not_done=[
            "остальные точки KV (KV-2, KV-4, KV-5, KV-7 линейки) - в воронке и блоках есть только KV-1, KV-3, KV-6 (funnel_space.KV): "
            "добавить точку в воронку - и задачи получат её без правки кода (флаг --kv)",
            "кластер ядер (cores 2-4) не входит в пространство: модель массы считает один газогенератор",
            "радиальная турбина из двух колёс считается в модели массы одним колесом на весь перепад (тяжелее - консервативно)",
            "межтурбинная камера подогрева - блок камеры по правилу основной при давлении газа перед ней: при низком давлении она в 2-3 раза "
            "тяжелее основной; последовательные камеры на горячем газе компактнее - уточнить на следующем уровне",
            "горение с ростом давления (RDC, ВР): выигрыш - модель перебора циклов с реализацией 25 % и штрафом турбины 2 п. п. (TRL 3); "
            "масса камеры - как у LPP; волновой ротор - плюс масса камеры (барабан не моделирован)",
            "выбросы (NOx, CO) показываются, но не отсеивают: ни один вариант уровня 0 камеры не выполняет все жёсткие цели",
            "рама модуля, газоходы, вентиляторы, катализатор - не калиброваны (оценки в evaluation.extra_mass_basis)",
            "КПД на тонну у авиационной практики (нижняя граница массы) не считается: все концепции - промышленная K",
            "острова внутри эпохи не обмениваются: эпоха 2 - слитый архив как родители (--parents), по желанию",
            "живой прогон AWS не выполнялся (только runner.py plan)",
        ],
        runtime_s=round(time.time() - t0, 1),
    )
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1, default=str)
    for kv, r in runs.items():
        b = r["best"]
        print(f"{kv}: оценок {r['evals']}, допустимых {100 * r['feasible_share']:.1f} %, ниш {r['niches']['feasible']}, "
              f"ЦП {1000 * r['cpu_s_per_eval']['all']:.1f} мс/концепцию; лучший КПД/т {b['fit'] if b else None} ({b['niche'] if b else ''})")
    print("записано", os.path.relpath(a.out, os.path.dirname(os.path.dirname(os.path.dirname(HERE)))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
