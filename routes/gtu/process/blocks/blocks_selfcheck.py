"""Самопроверка библиотеки блоков «от топлива до вала» -> blocks_results.json (рядом) и полные записи примеров в scratchpad.

python3 routes/gtu/process/blocks/blocks_selfcheck.py [--out <файл>] [--full <папка; по умолчанию $GTU_SCRATCH>]
AWS не вызывается; около 15 с на одном ядре. Всё - оценка уровня облика (README узла).
"""
import argparse
import datetime
import hashlib
import json
import os
import platform
import sys
import tempfile
import time

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import blocks  # noqa: E402
from blocks import _track as T  # noqa: E402
from blocks import reference  # noqa: E402
from blocks.catalog import catalog  # noqa: E402

SCRATCH = os.environ.get("GTU_SCRATCH") or os.path.join(tempfile.gettempdir(), "gtu_process")   # полные записи - вне репозитория
CODE = ["__init__.py", "_track.py", "concept.py", "chain.py", "fuel.py", "compressor.py", "combustor.py", "turbine.py", "hx.py",
        "shaft.py", "compat.py", "evaluator.py", "reference.py", "catalog.py"]


def sha(p):
    with open(p, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()[:16]


def code_sha():
    h = hashlib.sha256()
    for fn in CODE:
        with open(os.path.join(HERE, fn), "rb") as f:
            h.update(f.read())
    return h.hexdigest()[:16]


def compact(r):
    """строка примера: главное для отбора и для модели массы"""
    out = dict(name=r["name"], kv=r["kv"], feasible=r["feasible"], reject=r["reject"], warnings=r["warnings"],
               eff_shaft=r.get("eff_shaft"), eff_el=r.get("eff_el"), w_kJkg=r.get("w_kJkg"), G_air=r.get("G_air"), W_exh=r.get("W_exh"),
               EGT_C=r.get("EGT_C"), T41_C=r.get("T41_C"), PR=r.get("PR"), T_rec_in_C=r.get("T_rec_in_C"), water_t_h=r.get("water_t_h"),
               nt=r.get("nt"), npt=r.get("npt"))
    if r.get("interface"):
        i = r["interface"]
        out.update(dG_pct=i["dG_pct"], dWexh_pct=i["dWexh_pct"], dEGT=i["dEGT"], P_at_G_cap_kW=i["P_at_G_cap_kW"])
    if r.get("parts"):
        p = r["parts"]
        out.update(parts_engine=p["engine"]["parts"], assemblies=p["engine"]["assemblies"], aux_added=p["aux"]["added"],
                   n_eff_total=p["n_eff_total"], s_parts=p["s_parts_total"])
    if r.get("trl"):
        out.update(trl_min=r["trl"]["min"], k_trl=r["trl"]["k_trl"],
                   trl_limiting=[x["tech"] for x in r["trl"]["items"] if x["trl"] == r["trl"]["min"]])
    if r.get("plugplay"):
        out.update(plugplay=r["plugplay"]["score"], pnp_cycles=r["plugplay"]["pnp_cycles"], pp_blocked=r["plugplay"]["blocked"])
    if r.get("sizing"):
        s = r["sizing"]
        out.update(n_rpm=s["engine"]["n_rpm"], D_env=s["engine"]["D_env"], envelope=s["engine"]["envelope"],
                   hx_core_mass_kg=s["hx_core_mass_kg"], mass_funnel_estimate_kg=s["mass_funnel_estimate_kg"])
    if r.get("funnel_score"):
        out["funnel_score"] = r["funnel_score"]["score"]
    cb = r.get("blocks", {}).get("combustor", {})
    out["comb_dp"] = cb.get("dp")
    out["runtime_s"] = r["runtime_s"]
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "blocks_results.json"))
    ap.add_argument("--full", default=SCRATCH)
    a = ap.parse_args(argv)
    t0 = time.time()
    f3 = reference.check_frame3()
    rec = f3.pop("record")
    cyc5 = reference.check_cycles(reference.CYCLES_KEYS)
    cyc_x = reference.check_cycles(reference.EXTRA_KEYS)
    recup = reference.check_recup_gtu_cycle()
    ex_full, ex_rows = {}, []
    for name in reference.EXAMPLES:
        r = blocks.evaluate(reference.example(name), detail=False)
        ex_full[name] = r
        ex_rows.append(compact(r))
    import cantera
    import numpy
    import scipy
    eng = T.load_json(os.path.join(T.FLEET, "gtu_fleet_taurus60.json"))
    eng_kg = next(x["engine_kg"] for x in eng["records"] if x.get("engine_kg"))
    m_f3 = rec["sizing"]["mass_funnel_estimate_kg"]
    inputs = [T.F_FAMILY, T.F_FRAME3, T.F_CYCLES, T.F_COMB_TARGETS, T.F_COMB_L0, os.path.join(T.FLEET, "gtu_fleet_taurus60.json"),
              os.path.join(T.GTU, "gtu_cycle.py"), os.path.join(T.GTU, "gtu_explore.py"), os.path.join(T.FUNNEL, "funnel_space.py"),
              os.path.join(T.FUNNEL, "funnel_l0.py"), os.path.join(T.FUNNEL, "make_model.py"), os.path.join(T.FUNNEL, "parts_model.py"),
              os.path.join(T.FUNNEL, "plugplay_model.py"), os.path.join(T.FUNNEL, "service_model.py"), os.path.join(T.CYCLES, "cycles.py"),
              os.path.join(T.COMB, "comb_space.py")]
    out = dict(
        title="Библиотека блоков процесса «от топлива до вала»: варианты, быстрые модели, совместимость, единый вызов evaluate",
        method=("блоки (подготовка топлива, компрессор, камера, турбины, теплообменники, вал и стык пакета) поверх моделей трека без "
                "их правки: КПД узлов, геометрия, правила, детали, Plug&Play, изготовимость - воронка (funnel_space, funnel_l0, "
                "make_model, parts/service/plugplay_model); сквозная 0D-цепочка на таблицах gtu_explore (как funnel_space._core) с "
                "промохлаждением, рекуператором, подогревом, ростом давления, подогревом топлива, IBC; химическая рекуперация - "
                "приращение cycles.evaluate (Cantera); теплообменники - cycles.size_hx и cycles.parts; камера - итог уровня 0 "
                "воронки камеры и targets.json; сверка - gtu_cycle.cycle / recup, frame3, cycles_results.json"),
        software=f"python {platform.python_version()}, numpy {numpy.__version__}, scipy {scipy.__version__}, cantera {cantera.__version__}",
        date=datetime.date.today().isoformat(),
        status="оценка уровня облика; КПД узлов, потери, пределы, детали и TRL - принятые величины трека и оценки; массы не калиброваны",
        code_sha=code_sha(),
        inputs=[dict(file=T.rel(p), sha256_16=sha(p)) for p in inputs],
        assumptions=[
            "свойства газа - таблицы gtu_explore (Cantera GRI-Mech 3.0, метан, ISO), смесь воздуха и стехиометрических продуктов",
            "КПД ступеней компрессора и турбин - модели воронки для базового ядра концепции; добавочные блоки (промохлаждение, "
            "подогрев, IBC) нагрузку ступеней и их КПД не пересчитывают",
            "охлаждение: доли make_model.COOL_AIR; повторная ступень (подогрев) - тем же классом",
            "рекуператор: потери 2,5 % / 3,5 % (gtu_cycle.recup), эффективность по наименьшей теплоёмкости; газ на входе <= 700 °C",
            "камера: потери - лучший по рейтингу вариант пары «схема × тип камеры» уровня 0 камеры (ни один вариант уровня 0 не "
            "выполняет все цели выбросов); track - 4 % (кольцевая) / 5 % (выносная)",
            "горение с ростом давления: 25 % идеального выигрыша (эквивалентное давление при постоянном объёме), турбина -2 п. п.",
            "химическая рекуперация: приращение перебора циклов при тех же π, T41, ε, промохлаждении, подогреве (суперпозиция)",
            "быстроходный вал: частота по φ первой центробежной ступени, отношение частот 0,70-1,60 (сменная пара в корпусе редуктора)",
            "добавочные узлы теплообмена в счёт деталей - по 1 эффективной детали на позицию cycles.parts",
            "массы модулей - формулы воронки; это вход модели массы, не итог",
        ],
        checks=dict(
            kv3_frame3=f3,
            cycles_5=cyc5,
            cycles_extra=cyc_x,
            gtu_cycle_recup=recup,
            tolerances=dict(same_model="до округления (0,005 п. п.)", other_property_impl_pp=0.15, flow_pct=0.6, temp_K=2.0,
                            frame3_calibration_pp=2.0),
        ),
        mass_reference=dict(engine_kg_fleet=eng_kg, record="fleet/gtu_fleet_taurus60.json (MD - единственная масса двигателя в паспортах)",
                            funnel_estimate_frame3_kg=m_f3, ratio=round(eng_kg / m_f3, 3),
                            note="формула массы воронки (funnel_l0.envelope) на опорной точке KV-3 даёт меньше паспортной массы двигателя "
                                 "Taurus 60 в указанное число раз: модель массы должна калиброваться по машинам, а не брать эту оценку"),
        catalog=catalog(),
        examples=ex_rows,
        runtime_s=round(time.time() - t0, 2),
    )
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    os.makedirs(a.full, exist_ok=True)
    full = os.path.join(a.full, "blocks_examples_full.json")
    with open(full, "w", encoding="utf-8") as f:
        json.dump(dict(frame3=rec, examples=ex_full), f, ensure_ascii=False, indent=1)
    print(f"записано: {a.out}; полные записи: {full}; {out['runtime_s']} с")
    print(f"KV-3: evaluate {f3['evaluate']['eff_shaft']} %, воронка {f3['funnel']['eff_shaft']} %, frame3 {f3['frame3']['eff_shaft']} % "
          f"(разница {f3['d_eff_vs_frame3_pp']} п. п.); с КПД frame3 {f3['eta_given']['eff_shaft']} против gtu_cycle "
          f"{f3['gtu_cycle']['eff_shaft']}")
    for r in cyc5 + cyc_x:
        print(f"  {r['key']}: Δη {r['d_eff_pp']} п. п., ΔG {r['d_W_air_pct']} %, ΔEGT {r['d_EGT_K']} K")
    for r in ex_rows:
        print(f"  {r['name'][:52]:52s} η {r['eff_shaft']} G {r['G_air']} EGT {r['EGT_C']} TRL {r.get('trl_min')} "
              f"{'допустим' if r['feasible'] else 'нет: ' + ','.join(r['reject'][:3])}")


if __name__ == "__main__":
    main()
