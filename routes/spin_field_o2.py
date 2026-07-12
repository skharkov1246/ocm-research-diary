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

OUT = os.path.join(HERE, "spin_field_o2_results.json")
EV = 27.211386245988
VA = 51.42206747
FIELDS = [0.0, 0.005, -0.005, 0.010, -0.010]     # а.е. → 0, ±0.257, ±0.514 В/Å
SPINS = [0, 2, 4]                                 # 2S всего комплекса
XC = os.environ.get("XC", "pbe")

def say(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

def atomic_save(o):
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w") as f: json.dump(o, f, ensure_ascii=False, indent=1)
    os.replace(tmp, OUT)

def build():
    """[Fe(NH3)4(O2)] end-on: Fe в нуле, 4 N в плоскости xy, O2 вдоль +z."""
    fe_n, fe_o, o_o, n_h = 2.10, 1.90, 1.25, 1.02
    atoms = [["Fe", (0.0, 0.0, 0.0)]]
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
    fields = [0.0] if a.smoke else FIELDS
    for f in fields:
        fk = f"{f:+.3f}"; blk = res["scan"].setdefault(fk, {})
        for s in SPINS:
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
    # анализ: основной спин и зазор HS−LS на каждом поле
    rows = []
    for fk, blk in sorted(res["scan"].items(), key=lambda kv: float(kv[0])):
        es = {s: blk.get(f"2S={s}", {}).get("e") for s in SPINS}
        if all(es[s] is not None for s in SPINS):
            gs = min(SPINS, key=lambda s: es[s])
            gap_qt = round((es[4] - es[2]) * EV, 4)      # квинтет − триплет
            gap_st = round((es[2] - es[0]) * EV, 4)      # триплет − синглет
            qo2 = blk.get(f"2S={gs}", {}).get("q_O2")
            rows.append({"field_V_per_A": round(float(fk) * VA, 3), "ground_2S": gs,
                         "gap_Q_minus_T_eV": gap_qt, "gap_T_minus_S_eV": gap_st,
                         "q_O2_ground": qo2})
    res["summary_rows"] = rows
    if len(rows) >= 3:
        F = np.array([r["field_V_per_A"] for r in rows])
        res["dGapQT_dF_eV_per_VA"] = round(float(np.polyfit(F, [r["gap_Q_minus_T_eV"] for r in rows], 1)[0]), 4)
        res["dGapTS_dF_eV_per_VA"] = round(float(np.polyfit(F, [r["gap_T_minus_S_eV"] for r in rows], 1)[0]), 4)
        res["ground_switches_with_field"] = len(set(r["ground_2S"] for r in rows)) > 1
    res["status"] = "smoke" if a.smoke else ("ok" if len(rows) == len(FIELDS) else "incomplete")
    atomic_save(res)
    say(f"status={res['status']} switch={res.get('ground_switches_with_field')} "
        f"dGapQT/dF={res.get('dGapQT_dF_eV_per_VA')} dGapTS/dF={res.get('dGapTS_dF_eV_per_VA')} -> {OUT}")
    return 0

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--smoke", action="store_true")
    sys.exit(main(ap.parse_args()))
