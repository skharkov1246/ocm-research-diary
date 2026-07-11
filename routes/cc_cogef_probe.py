#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""routes/cc_cogef_probe.py — FLIP-эксперимент №2: механохимия рециклинга (COGEF).

Переворот: рвать C–C полиолефина СИЛОЙ (мельница/экструзия), а не теплом
(пиролиз 400°C+). Классика метода — COGEF (COnstrained Geometries simulate
External Force, Beyer 2000): релаксированный скан растяжения связи E(d), из
него — сила разрыва F_max = max dE/dd и барьер гомолиза под внешней тягой F:
E_F(d) = E(d) − F·(d−d₀), ΔE‡(F) = max[E_F] − min[E_F].

Модель v1: н-бутан, центральная C2–C3 (минимальный зигзаг-фрагмент ПЭ; строго говоря —
это НЕ цепь ПЭ: нет натяжения соседних звеньев и конц. эффекты). UKS-PBE/
def2-SVP+DF, констрейн-релаксация geomeTRIC (дистанция C2–C3 фиксируется по
сетке 1.53→3.0 Å), на растянутых точках — broken-symmetry guess (гомолиз →
бирадикал). Дескриптор: ΔE‡(F) и сила F*, при которой барьер падает до
~0.75 эВ (комнатная кинетика мс-масштаба).

Ограничения: PBE на гомолизе неточен (стат. корреляция бирадикала — потом можно
CASSCF(2,2) поверх, наш конвейер умеет); модель — бутан, не ПЭ-цепь; сила
статическая (не импульс мельницы); без энтропии/температуры.

Запуск: OMP_NUM_THREADS=4 python3 -u routes/cc_cogef_probe.py
Smoke:  OMP_NUM_THREADS=1 python3 routes/cc_cogef_probe.py --smoke (~2 мин)
Выход:  routes/cc_cogef_results.json (atomic, incremental, рестарт по точкам)
"""
import os, sys, json, time, argparse, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pyscf import gto, dft, scf, lib  # noqa: E402

OUT = os.path.join(HERE, "cc_cogef_results.json")
EV = 27.211386245988
NN_PER_HA_BOHR = 82.38723498         # 1 Ha/Bohr = 82.387 нН
A2B = 1.8897259886
DGRID = [1.53, 1.65, 1.80, 1.95, 2.10, 2.25, 2.40, 2.60, 2.80, 3.00]  # Å, C2–C3
FGRID_NN = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]              # нН

# н-бутан, зигзаг (стартовая геометрия; релаксируется). Индексы C: 0,1,2,3 (0-based);
# центральная связь — атомы 2 и 3 в 1-based нумерации geomeTRIC = C1(идекс1)-C2(индекс2)?
# Явно: C-атомы идут первыми, центральная связь C[1]-C[2] (0-based) = "2 3" (1-based).
BUTANE = """
C  -1.930  -0.360   0.000
C  -0.560   0.310   0.000
C   0.560  -0.310   0.000
C   1.930   0.360   0.000
H  -2.720   0.395   0.000
H  -2.060  -0.990   0.885
H  -2.060  -0.990  -0.885
H  -0.470   0.955   0.885
H  -0.470   0.955  -0.885
H   0.470  -0.955   0.885
H   0.470  -0.955  -0.885
H   2.720  -0.395   0.000
H   2.060   0.990   0.885
H   2.060   0.990  -0.885
"""

def say(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

def atomic_save(obj):
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1)
    os.replace(tmp, OUT)

def make_mf(mol, bs_guess=False):
    mf = dft.UKS(mol).density_fit()
    mf.xc = "pbe"
    mf.verbose = 0
    mf.conv_tol = 1e-8
    mf.max_cycle = 150
    if bs_guess:
        dm = mf.get_init_guess()
        # ломаем спин-симметрию: альфа на левом C-фрагменте, бета на правом
        dm = (dm[0] * 1.05, dm[1] * 0.95)
        return mf, dm
    return mf, None

def relax_at_d(mol0, d, prev_xyz=None, maxsteps=60, bs=False):
    """Констрейн-релаксация с фикс. d(C2–C3) через geomeTRIC."""
    from pyscf.geomopt import geometric_solver
    atoms = prev_xyz if prev_xyz is not None else [
        [ln.split()[0], tuple(map(float, ln.split()[1:4]))]
        for ln in BUTANE.strip().splitlines()]
    # ставим стартовое d
    a = np.array([list(x[1]) for x in atoms])
    v = a[2] - a[1]; v /= np.linalg.norm(v)
    shift = (d - np.linalg.norm(a[2] - a[1])) / 2
    for i in (2, 3): a[i] += v * shift
    for i in (0, 1): a[i] -= v * shift
    # H следуют за своими C грубо: релаксация поправит
    atoms = [[atoms[k][0], tuple(a[k] if k < 4 else np.array(atoms[k][1]))]
             for k in range(len(atoms))]
    mol = gto.M(atom=atoms, basis="def2-svp", spin=0, verbose=0)
    mf, dm = make_mf(mol, bs_guess=bs)
    if dm is not None: mf.kernel(dm0=dm)
    else: mf.kernel()
    if not mf.converged:
        mfn = mf.newton(); mfn.max_cycle = 60
        mfn.kernel(mf.mo_coeff, mf.mo_occ)
        mf.mo_coeff, mf.mo_occ = mfn.mo_coeff, mfn.mo_occ
        mf.converged = bool(mfn.converged); mf.e_tot = float(mfn.e_tot)
    fd, cpath = tempfile.mkstemp(suffix=".txt")
    with os.fdopen(fd, "w") as fh:
        fh.write(f"$freeze\ndistance 2 3\n")   # 1-based: центральные C
    try:
        mol_eq = geometric_solver.optimize(mf, constraints=cpath, maxsteps=maxsteps)
    finally:
        os.unlink(cpath)
    mf2, dm2 = make_mf(mol_eq, bs_guess=bs)
    if dm2 is not None: e = mf2.kernel(dm0=dm2)
    else: e = mf2.kernel()
    ss = float(mf2.spin_square()[0])
    xyz = [[mol_eq.atom_symbol(i), tuple(mol_eq.atom_coords()[i] / A2B)]
           for i in range(mol_eq.natm)]
    return float(e), bool(mf2.converged), ss, xyz

def cogef_analysis(res):
    pts = sorted((float(k), v["e"]) for k, v in res["scan"].items()
                 if v.get("converged"))
    if len(pts) < 4: return
    D = np.array([p[0] for p in pts]); E = np.array([p[1] for p in pts])
    E = (E - E.min())
    # COGEF-сила: максимум производной (конечные разности)
    dEdd = np.gradient(E, D)                      # Ha/Å
    fmax_nn = float(dEdd.max() * NN_PER_HA_BOHR / A2B)
    rows = []
    for F in FGRID_NN:
        f_au = F / NN_PER_HA_BOHR * A2B           # Ha/Å
        Etil = E - f_au * (D - D[0])
        i_ts = int(np.argmax(Etil)); i_min = int(np.argmin(Etil[:i_ts + 1])) if i_ts > 0 else 0
        bar = float((Etil[i_ts] - Etil[i_min]) * EV)
        rows.append({"F_nN": F, "barrier_eV": round(bar, 3),
                     "d_ts_A": round(float(D[i_ts]), 2)})
    res["cogef"] = {
        "F_max_nN": round(fmax_nn, 2),
        "barrier_vs_force": rows,
        "note": "ΔE‡(F)=max[E−F·d]−min[E−F·d] по релакс-скану; барьер при F=0 "
                "на сетке до 3.0 Å — нижняя оценка (полный гомолиз дальше)",
    }

def main_real(args):
    say(f"COGEF probe (бутан C2–C3) — threads={lib.num_threads()}")
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    res.setdefault("meta", {
        "what": "COGEF: relaxed constrained stretch of the central C–C of "
                "n-butane (PE minimal fragment), UKS-PBE/def2-SVP+DF, geomeTRIC; "
                "BS guess on stretched points; force-modified barrier ΔE‡(F)",
        "grid_A": DGRID, "forces_nN": FGRID_NN,
        "honesty": "butane != PE chain (no chain tension/entanglement, end "
                   "effects); PBE poor for homolysis (CASSCF(2,2) follow-up "
                   "possible); static force != mill impulse; no entropy/T",
        "flip_motivation": "cc-activation: рвать C-C силой (мельница/экструзия) "
                           "вместо пиролиза 400C+ — mechanochemistry recycling",
    })
    res.setdefault("scan", {})
    atomic_save(res)
    prev = None
    for d in DGRID:
        key = f"{d:.2f}"
        if key in res["scan"] and res["scan"][key].get("converged"):
            say(f"  d={key} Å — уже есть (рестарт)"); prev = res["scan"][key].get("xyz"); continue
        bs = d >= 2.0
        t0 = time.time()
        try:
            e, conv, ss, xyz = relax_at_d(None, d, prev_xyz=prev, bs=bs)
            res["scan"][key] = {"d_A": d, "e": e, "converged": conv,
                                "s2": round(ss, 3), "bs_guess": bs,
                                "t_s": round(time.time() - t0, 1), "xyz": xyz}
            prev = xyz
            say(f"  d={key} Å: E={e:.6f} conv={conv} <S²>={ss:.2f} ({res['scan'][key]['t_s']}s)")
        except Exception as exc:
            res["scan"][key] = {"d_A": d, "converged": False,
                                "error": f"{type(exc).__name__}: {str(exc)[:100]}"}
            say(f"  d={key} Å: FAILED {exc}")
        atomic_save(res)
    cogef_analysis(res)
    done = sum(1 for v in res["scan"].values() if v.get("converged"))
    res["status"] = "ok" if done == len(DGRID) else f"partial({done}/{len(DGRID)})"
    atomic_save(res)
    say(f"status={res['status']} -> {OUT}")

def main_smoke(args):
    global DGRID, FGRID_NN, OUT
    OUT = os.path.join(tempfile.mkdtemp(prefix="cogef_smoke_"), "smoke.json")
    DGRID = [1.53, 2.10]; FGRID_NN = [0.0, 2.0]
    res = {"meta": {"mode": "smoke"}, "scan": {}}
    prev = None
    for d in DGRID:
        e, conv, ss, xyz = relax_at_d(None, d, prev_xyz=prev, maxsteps=15, bs=(d >= 2.0))
        res["scan"][f"{d:.2f}"] = {"d_A": d, "e": e, "converged": conv, "s2": round(ss, 3), "xyz": xyz}
        prev = xyz
        say(f"  d={d}: E={e:.6f} conv={conv} <S²>={ss:.2f}")
    ok = all(v["converged"] for v in res["scan"].values()) and \
        res["scan"]["2.10"]["e"] > res["scan"]["1.53"]["e"]
    say(f"SMOKE {'PASS' if ok else 'FAIL'} (растяжение поднимает энергию: "
        f"{(res['scan']['2.10']['e']-res['scan']['1.53']['e'])*EV:.2f} эВ)")
    return 0 if ok else 1

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    sys.exit(main_smoke(a) if a.smoke else main_real(a))
