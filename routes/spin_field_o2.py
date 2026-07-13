#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""routes/spin_field_o2.py — ДЕРЗКАЯ НОВАЯ СТАВКА: электрическое переключение спина
кислорода. Может ли ориентированное поле управлять спиновым состоянием активации
O₂ — гасить инертный триплет и включать реакционный канал — «спиновый транзистор»
оксидации на электроде?

Идея (флип, промышленно не применяется). O₂ инертен, потому что основной терм —
триплет, а большинство органических субстратов синглетны: спиновый запрет держит
мягкую оксидацию под замком (отсюда жёсткие условия и дорогие катализаторы).
Металлоцентр снимает запрет через двухспиновую реактивность (two-state reactivity)
— переход между спиновыми поверхностями у точки пересечения (MECP). Гипотеза:
внешнее ЭЛЕКТРИЧЕСКОЕ поле сдвигает относительные энергии спиновых поверхностей
(разная поляризуемость/перенос заряда Fe→π*(O₂) в разных спинах) и потому двигает
спиновый зазор — а на электроде поле двойного слоя (0.1–1 В/Å) даётся даром. Если
dΔE(спин)/dF велик — поле есть ручка спин-селекции оксидации.

Модель v1 (минимальная): [Fe(NH₃)₄(O₂)] end-on — Fe(II) в N₄-поле (гем-подобное
окружение), O₂ аксиально. Сканируем спины всего комплекса 2S=0,2,4 (синглет/
триплет/квинтет-подобные) в аксиальном поле 0..±0.51 В/Å. Дескрипторы: E(спин,F),
спиновый зазор ΔE(HS−LS)(F) и его наклон dGap/dF; заряд Малликена O₂ (степень
активации: супероксо/пероксо) vs поле. Ловим СМЕНУ основного спина полем.

ОГРАНИЧЕНИЯ (важно): спиновые зазоры на PBE функционал-зависимы и НЕНАДЁЖНЫ по
абсолюту — проба меряет ЧУВСТВИТЕЛЬНОСТЬ зазора к полю (сдвиг), не абсолютный
зазор; кластер без белка/растворителя; замороженная идеализированная геометрия
(сдвиги по полю при фиксированном каркасе честны как чувствительность); def2-SVP;
v1. Гибрид-верификация (PBE0/B3LYP) — отдельным прогоном, как в вулкане.

Запуск: OMP_NUM_THREADS=4 python3 -u routes/spin_field_o2.py
Smoke:  --smoke (F=0, три спина; ~5-10 мин)
Выход:  routes/spin_field_o2_results.json
"""
import os, sys, json, time, argparse, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
from pyscf import gto, dft, lib  # noqa: E402
from pyscf.data.elements import charge as Z  # noqa: E402

METAL = os.environ.get("METAL", "Fe")
XC = os.environ.get("XC", "pbe")
_tag = ("" if METAL == "Fe" else "_" + METAL.lower()) + ("" if XC == "pbe" else "_" + XC)
OUT = os.path.join(HERE, f"spin_field_o2{_tag}_results.json")
EV = 27.211386245988
VA = 51.42206747
FIELDS = [0.0, 0.005, -0.005, 0.010, -0.010]     # а.е. → 0, ±0.257, ±0.514 В/Å
SPINS = [0, 1, 2, 3, 4]                           # 2S всего комплекса (фильтр по чётности)

def say(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

def atomic_save(o):
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w") as f: json.dump(o, f, ensure_ascii=False, indent=1)
    os.replace(tmp, OUT)

def build():
    """[M(NH3)4(O2)] end-on: M в нуле, 4 N в плоскости xy, O2 вдоль +z."""
    fe_n, fe_o, o_o, n_h = 2.10, 1.90, 1.25, 1.02
    atoms = [[METAL, (0.0, 0.0, 0.0)]]
    for k in range(4):
        a = np.pi / 2 * k
        npos = np.array([fe_n * np.cos(a), fe_n * np.sin(a), 0.0])
        atoms.append(["N", tuple(npos)])
        # 3 H на N: один «наружу» от Fe, два вбок с наклоном вверх
        out = npos / np.linalg.norm(npos)
        up = np.array([0, 0, 1.0])
        side = np.cross(out, up)
        for d in (out * 0.9 + up * 0.5, out * 0.6 - up * 0.4 + side * 0.7,
                  out * 0.6 - up * 0.4 - side * 0.7):
            h = npos + d / np.linalg.norm(d) * n_h
            atoms.append(["H", tuple(h)])
    atoms.append(["O", (0.0, 0.0, fe_o)])
    atoms.append(["O", (0.0, 0.0, fe_o + o_o)])
    return atoms

def scf(spin, f_au):
    mol = gto.M(atom=build(), basis="def2-svp", spin=spin, charge=0, verbose=0)
    mf = dft.UKS(mol).density_fit()
    mf.xc = XC; mf.conv_tol = 1e-8; mf.max_cycle = 200; mf.level_shift = 0.3
    if abs(f_au) > 0:
        mol.set_common_orig((0.0, 0.0, 0.0))
        ao = mol.intor_symmetric("int1e_r", comp=3)
        h0 = mf.get_hcore(); E = np.array([0.0, 0.0, f_au])
        mf.get_hcore = lambda *a: h0 + np.einsum("x,xij->ij", E, ao)
    e = mf.kernel()
    if not mf.converged:
        mfn = mf.newton(); mfn.max_cycle = 80
        e = mfn.kernel(mf.mo_coeff, mf.mo_occ)
        mf.converged = bool(mfn.converged); e = float(mfn.e_tot)
    try:
        pop = mf.mulliken_pop(verbose=0)[1]; q_o2 = float(pop[-1] + pop[-2])
    except Exception:
        q_o2 = None
    return float(e), bool(mf.converged), q_o2

def main(a):
    say(f"spin-field O2 probe — xc={XC} threads={lib.num_threads()}")
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    res.setdefault("meta", {
        "what": "[Fe(NH3)4(O2)] end-on; scan whole-complex spin 2S=0,2,4 in axial "
                "field; descriptors E(spin,F), spin gap and dGap/dF, Mulliken q(O2); "
                f"UKS-{XC}/def2-SVP+DF; looks for FIELD-INDUCED spin-ground-state switch",
        "limits": "PBE spin gaps are functional-sensitive & unreliable in absolute "
                  "terms -> probe measures SENSITIVITY of the gap to field, not the "
                  "gap itself; cluster, gas phase, frozen idealized scaffold; v1",
        "fields_au": FIELDS, "spins_2S": SPINS,
    })
    res.setdefault("scan", {})
    # чётность 2S по числу электронов нейтрального [M(NH3)4(O2)]
    nel = Z(METAL) + 4 * (7 + 3) + 2 * 8           # M + 4 NH3 + O2
    spins = [s for s in SPINS if s % 2 == nel % 2] or [nel % 2]
    res["meta"]["metal"] = METAL; res["meta"]["spins_used"] = spins
    fields = [0.0] if a.smoke else FIELDS
    for f in fields:
        fk = f"{f:+.3f}"; blk = res["scan"].setdefault(fk, {})
        for s in spins:
            sk = f"2S={s}"
            if blk.get(sk, {}).get("e") is not None: continue
            t0 = time.time()
            try:
                e, conv, q = scf(s, f)
                blk[sk] = {"e": round(e, 6), "converged": conv, "q_O2": q,
                           "t_s": round(time.time() - t0, 1)}
                say(f"  F={fk} {sk}: E={e:.5f} conv={conv} q(O2)={q if q is None else round(q,3)} ({blk[sk]['t_s']}s)")
            except Exception as exc:
                blk[sk] = {"e": None, "converged": False, "error": f"{type(exc).__name__}: {str(exc)[:70]}"}
                say(f"  F={fk} {sk}: FAILED {exc}")
            atomic_save(res)
    # анализ: основной спин на каждом поле + зазор к ближайшему конкуренту
    per = []
    for fk, blk in sorted(res["scan"].items(), key=lambda kv: float(kv[0])):
        es = {s: blk[f"2S={s}"]["e"] for s in spins
              if blk.get(f"2S={s}", {}).get("e") is not None and blk[f"2S={s}"].get("converged")}
        if len(es) >= 2:
            gs = min(es, key=es.get)
            qo2 = blk.get(f"2S={gs}", {}).get("q_O2")
            per.append((round(float(fk) * VA, 3), gs, es, qo2))
    # пара переключения = спины, бывшие основными хоть на одном поле
    grounds = sorted(set(p[1] for p in per))
    pair = (grounds[0], grounds[-1]) if len(grounds) >= 2 else (
        (grounds[0], sorted(set().union(*[set(p[2]) for p in per]) - set(grounds))[0])
        if per and len(set().union(*[set(p[2]) for p in per])) >= 2 else None)
    rows = []
    for F, gs, es, qo2 in per:
        gap = round((es[pair[1]] - es[pair[0]]) * EV, 4) if pair and pair[0] in es and pair[1] in es else None
        rows.append({"field_V_per_A": F, "ground_2S": gs,
                     f"gap_2S{pair[1]}_minus_2S{pair[0]}_eV" if pair else "gap_eV": gap,
                     "q_O2_ground": qo2})
    res["summary_rows"] = rows
    res["switching_pair_2S"] = list(pair) if pair else None
    res["ground_switches_with_field"] = len(grounds) > 1
    gaps = [r.get(f"gap_2S{pair[1]}_minus_2S{pair[0]}_eV") for r in rows] if pair else []
    gaps = [g for g in gaps if g is not None]
    if pair and len(gaps) >= 3:
        F = np.array([r["field_V_per_A"] for r in rows if r.get(f"gap_2S{pair[1]}_minus_2S{pair[0]}_eV") is not None])
        res["dGap_dF_eV_per_VA"] = round(float(np.polyfit(F, np.array(gaps), 1)[0]), 4)
    # чувствительность заряда O2 к полю
    qs = [(r["field_V_per_A"], r["q_O2_ground"]) for r in rows if r["q_O2_ground"] is not None]
    if len(qs) >= 3:
        res["dqO2_dF_e_per_VA"] = round(float(np.polyfit([x[0] for x in qs], [x[1] for x in qs], 1)[0]), 4)
    res["status"] = "smoke" if a.smoke else ("ok" if len(rows) == len(FIELDS) else "incomplete")
    atomic_save(res)
    say(f"[{METAL}] status={res['status']} switch={res.get('ground_switches_with_field')} "
        f"pair={res.get('switching_pair_2S')} dGap/dF={res.get('dGap_dF_eV_per_VA')} "
        f"dqO2/dF={res.get('dqO2_dF_e_per_VA')} -> {OUT}")
    return 0

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--smoke", action="store_true")
    sys.exit(main(ap.parse_args()))
