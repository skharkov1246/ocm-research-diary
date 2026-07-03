#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/oer_scaling_cubane.py — Stage 2 трека oer-scaling: мультиреференсность
оксо-узла Co(IV)=O в усечённом кубан-коре [Co4O4].

Модель (v1, БЕЗ релаксации — single-point на ИДЕАЛИЗИРОВАННОЙ геометрии):
  Кор молекулярного кубана Co4O4(OAc)4(py)4 (JACS 2011, 10.1021/ja203877y),
  усечённый: ацетат -> OH- (сурогат), пиридин -> NH3 (сурогат).
    resting  : [Co4O4(OH)4(NH3)4]^0,  Co(III)4 d6 LS -> синглет базово
               (+ контроль: триплет);
    oxidized : [Co4O4(OH)3(O)(NH3)4]^0 — одна OH- заменена на оксо O^2-
               (формально PCET-снятие H с одной OH: одна Co(IV)=O),
               дублет базово (+ контроль: квартет).
  Геометрия — ИДЕАЛИЗАЦИЯ, не оптимизация: чередующийся куб Co4O4 с ребром
  Co-O = 1.90 A (середина кристаллографического диапазона кубанов ~1.86-1.92 A);
  Co-OH 1.87, Co(IV)=O 1.68, Co-N 1.96, O-H 0.97, N-H 1.02 A — типовые
  литературные значения-мотивы, НЕ наш расчёт. OH ставится транс к кубановому O
  вдоль локальной оси x, NH3 — транс вдоль y, шестая позиция вакантна
  (в реальном комплексе там второй ацетатный O). Итого 32 атома — чуть больше
  оценки "~25" из записи Stage 1 (честно фиксируем).

Протокол:
 1) UKS-PBE/def2-SVP + density fitting; удержание сходимости: level shift 0.25 ->
    при несходимости Newton(SOSCF) -> цикл ВНУТРЕННЕЙ стабильности до
    internal-stable (образец scf_stable из co2_to_fuels_casscf.py, UHF-вариант).
    Записываем <S^2>, gap, Малликен-спины ключевых атомов (локализация дырки).
 2) UNO-мост: натуральные орбитали UHF/UKS (диагонализация суммарной спиновой
    плотности в S-метрике) -> фиктивный ROHF-контейнер орбиталей. N_u уже на
    уровне UNO (базис-инвариантно) — первый индикатор.
 3) AVAS(Co 3d оксо-узла + O 2p оксо/OH + O 2p трёх кубановых O при этом Co),
    порог поднимается, пока пространство не влезет в КАП (14e,12o) — ниже
    честного FCI-лимита. openshell_option=2: все SOMO принудительно в активном.
 4) CASCI (БЕЗ орбитальной оптимизации, честно = "CASCI на UNO-орбиталях"),
    fix_spin. Пространство фиксируется ОДИН раз на геометрию (по опорному
    спин-состоянию: resting -> синглет, oxidized -> дублет), оба спина решаются
    в ОДНОМ пространстве (урок feo_cation.py v2: иначе сравнение нечестно).
    NOON, N_u = sum n_i(2-n_i), избыток N_u - 2S -> качественный вердикт
    "мультиреференсен ли узел Co(IV)=O".

ЧЕСТНО (ограничения v1):
  - single-point на идеализированной геометрии, БЕЗ релаксации — все энергии
    только качественные; спиновые зазоры на неоптимальной геометрии смещены;
  - кластер != электрод: нет растворителя/потенциала; eta и dG-дескриптор
    отсюда НЕ извлекаемы, честен только N_u;
  - E(CASCI) построена на HF-остове поверх UNO(UKS)-орбиталей — НЕ сравнивать
    с E(UKS) и не выдавать за термохимию;
  - формальная метка "Co(IV)=O на одном Co" проверяется Малликен-спинами:
    дырка может делокализоваться по кубу — тогда метка условна (запишем);
  - UKS-синглет может сломать спин-симметрию (<S^2> > 0) — это сигнал
    мультиреференсности уже на DFT-уровне, фиксируем, не прячем.

Файлы: oer_scaling_cubane_results.json — atomic write после КАЖДОЙ стадии,
       рестартопригодно; oer_scaling_cubane_<state>.chk.npz — орбитали/dm
       (не коммитить). --smoke пишет в отдельный *_smoke.json.

Запуск полный:  OMP_NUM_THREADS=4 python3 -u routes/oer_scaling_cubane.py
Дым-тест (сборка геометрии + xyz в JSON + мини-SCF/UNO/AVAS/CASCI на
фрагменте CoO+; <=3 мин, 1 ядро):
                OMP_NUM_THREADS=1 python3 routes/oer_scaling_cubane.py --smoke
"""
import os
import sys
import json
import time
import argparse
import tempfile
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_RESULTS = os.path.join(HERE, "oer_scaling_cubane_results.json")
CHK_TPL = os.path.join(HERE, "oer_scaling_cubane_{state}.chk.npz")

HARTREE_EV = 27.211386245988
XC = "pbe"
BASIS = "def2-svp"
LEVEL_SHIFT = 0.25
MAX_STAB_LOOPS = 4
CONV_TOL_SCF = 1e-8
CAS_MAX_ORB = 12          # кап активного пространства: <= (14e,12o)
CAS_MAX_ELEC = 14
AVAS_THRESHOLDS = (0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

# ── идеализированные длины связей (литературные мотивы, НЕ наш расчёт) ──
R_CO_O_CUBE = 1.90        # ребро куба Co-O (кристаллография кубанов ~1.86-1.92)
R_CO_OH = 1.87            # Co(III)-OH terminal (типовое)
R_CO_OXO = 1.68           # Co(IV)=O (типовой литературный мотив)
R_CO_N = 1.96             # Co(III)-NH3 (типовое)
R_O_H = 0.97
R_N_H = 1.02
ANG_CO_O_H = 109.5        # deg
ANG_CO_N_H = 110.0        # deg

GEOM_HONESTY = ("idealized single-point v1: cube edge Co-O=1.90 A, Co-OH=1.87, "
                "Co(IV)=O=1.68, Co-N=1.96, O-H=0.97, N-H=1.02 A are typical "
                "literature-motif values, NOT optimized and NOT our calculation; "
                "no geometry relaxation anywhere in v1; azimuthal phases of the "
                "NH3 umbrellas and O-H tilts are set by a deterministic greedy "
                "grid search that maximizes the minimal nonbonded contact of "
                "the rigid cage (steric de-clash only, not an optimization)")

STATES = (  # (geom, state, 2S)
    ("resting", "singlet", 0),
    ("resting", "triplet", 2),
    ("oxidized", "doublet", 1),
    ("oxidized", "quartet", 3),
)
CAS_REFERENCE = {"resting": "singlet", "oxidized": "doublet"}


def say(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def atomic_save(obj, path):
    """Атомарная запись JSON: tmp в той же директории + os.replace."""
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path) or ".",
                               prefix=os.path.basename(path) + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(obj, f, ensure_ascii=False, indent=1)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def _norm(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


# ════════════════════════════════════════════════════════════════════════
#  Сборка усечённого кубан-кора (идеализация, см. GEOM_HONESTY)
# ════════════════════════════════════════════════════════════════════════
def _assemble(oxidized, oxo_site, phases):
    """Сборка при заданных азимутальных фазах: phases[i] = (phi_OH, phi_NH3)."""
    d = R_CO_O_CUBE
    co_pos = [np.array(p, float) for p in
              [(0, 0, 0), (d, d, 0), (d, 0, d), (0, d, d)]]
    o_pos = [np.array(p, float) for p in
             [(d, 0, 0), (0, d, 0), (0, 0, d), (d, d, d)]]
    center = np.array([d / 2] * 3)

    atoms, bonds = [], []
    meta = {"co": [], "cubane_o": [], "lig_o": {}, "lig_h": {},
            "n": {}, "oxo_co": None, "oxo_o": None,
            "cubane_o_of_co": {}}
    for p in co_pos:
        meta["co"].append(len(atoms))
        atoms.append(("Co", p))
    for p in o_pos:
        meta["cubane_o"].append(len(atoms))
        atoms.append(("O", p))
    # связи Co-O куба + кубановые соседи каждого Co (отличаются одной координатой)
    for i, ico in enumerate(meta["co"]):
        nbs = [meta["cubane_o"][j] for j, po in enumerate(o_pos)
               if np.count_nonzero(np.abs(po - co_pos[i]) > 1e-6) == 1]
        assert len(nbs) == 3
        meta["cubane_o_of_co"][ico] = nbs
        bonds += [(ico, j) for j in nbs]

    ang_oh = np.deg2rad(ANG_CO_O_H)
    ang_nh = np.deg2rad(ANG_CO_N_H)
    for i, v in enumerate(co_pos):
        ico = meta["co"][i]
        phi_oh, phi_nh = phases[i]
        e_ext = _norm(v - center)                       # внешняя диагональ
        # кубановые соседи вдоль локальных осей x и y
        x_nb = next(po for po in o_pos
                    if abs(po[1] - v[1]) < 1e-6 and abs(po[2] - v[2]) < 1e-6)
        y_nb = next(po for po in o_pos
                    if abs(po[0] - v[0]) < 1e-6 and abs(po[2] - v[2]) < 1e-6)

        # --- OH (или оксо) транс к x-соседу ---
        dir_o = _norm(v - x_nb)
        is_oxo = bool(oxidized and i == oxo_site)
        r = R_CO_OXO if is_oxo else R_CO_OH
        o_xyz = v + r * dir_o
        io = len(atoms)
        atoms.append(("O", o_xyz))
        bonds.append((ico, io))
        meta["lig_o"][ico] = io
        if is_oxo:
            meta["oxo_co"], meta["oxo_o"] = ico, io
        else:
            t0 = _norm(e_ext - np.dot(e_ext, dir_o) * dir_o)
            t1 = np.cross(dir_o, t0)
            t = np.cos(phi_oh) * t0 + np.sin(phi_oh) * t1
            u = np.cos(np.pi - ang_oh) * (-dir_o) + np.sin(np.pi - ang_oh) * t
            ih = len(atoms)
            atoms.append(("H", o_xyz + R_O_H * _norm(u)))
            bonds.append((io, ih))
            meta["lig_h"][ico] = ih

        # --- NH3 транс к y-соседу ---
        a = _norm(v - y_nb)
        n_xyz = v + R_CO_N * a
        inn = len(atoms)
        atoms.append(("N", n_xyz))
        bonds.append((ico, inn))
        meta["n"][ico] = inn
        p_ = _norm(e_ext - np.dot(e_ext, a) * a)
        q_ = np.cross(a, p_)
        th = np.pi - ang_nh                              # угол N-H от оси +a
        for k in range(3):
            ph = 2 * np.pi * k / 3 + phi_nh
            dirh = (np.cos(th) * a
                    + np.sin(th) * (np.cos(ph) * p_ + np.sin(ph) * q_))
            ihh = len(atoms)
            atoms.append(("H", n_xyz + R_N_H * dirh))
            bonds.append((inn, ihh))
    meta["bonds"] = bonds
    return atoms, meta


def _min_nonbonded(atoms, bonds):
    xyz = np.array([p for _, p in atoms])
    bond_set = {tuple(sorted(b)) for b in bonds}
    dm = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    n = len(atoms)
    iu, ju = np.triu_indices(n, k=1)
    mask = np.array([(i, j) not in bond_set for i, j in zip(iu, ju)])
    vals = dm[iu[mask], ju[mask]]
    k = int(np.argmin(vals))
    return float(vals[k]), (int(iu[mask][k]), int(ju[mask][k]))


def build_cubane(oxidized=False, oxo_site=0):
    """[Co4O4(OH)4(NH3)4] (resting) или [Co4O4(OH)3(O)(NH3)4] (oxidized).

    Чередующийся куб: Co на вершинах чётной чётности, O — нечётной.
    На каждом Co: OH транс к кубановому O вдоль локальной оси x,
    NH3 — транс вдоль y; шестая (транс-z) позиция вакантна.
    Азимутальные фазы OH/NH3 подбираются детерминированным greedy-перебором
    (2 прохода, сетка 24 значения), максимизируя минимальный несвязевый
    контакт жёсткого каркаса — только стерическая разводка, НЕ оптимизация.
    Возвращает (atoms, meta): atoms=[(символ, np.array xyz)], meta — индексы.
    """
    phases = [[0.0, np.pi / 6] for _ in range(4)]
    grid = np.linspace(0.0, 2 * np.pi, 24, endpoint=False)
    for _sweep in range(2):
        for i in range(4):
            for slot in (0, 1):
                best_phi, best_d = phases[i][slot], -1.0
                for phi in grid:
                    phases[i][slot] = float(phi)
                    atoms, meta = _assemble(oxidized, oxo_site, phases)
                    dmin, _ = _min_nonbonded(atoms, meta["bonds"])
                    if dmin > best_d + 1e-9:
                        best_d, best_phi = dmin, float(phi)
                phases[i][slot] = best_phi
    return _assemble(oxidized, oxo_site, phases)


def geom_report(atoms, meta):
    """Формула, xyz, минимальная НЕсвязевая дистанция (проверка клэшей)."""
    n = len(atoms)
    dmin, pair = _min_nonbonded(atoms, meta["bonds"])
    formula = "".join(f"{el}{c}" for el, c in
                      sorted(Counter(s for s, _ in atoms).items()))
    lines = [f"{n}", "idealized truncated Co4O4 cubane core (v1, not optimized)"]
    lines += [f"{s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}" for s, p in atoms]
    return {"natoms": n, "formula": formula,
            "min_nonbonded_A": round(dmin, 3),
            "min_nonbonded_pair": [f"{pair[0]}:{atoms[pair[0]][0]}",
                                   f"{pair[1]}:{atoms[pair[1]][0]}"],
            "xyz": "\n".join(lines)}


def make_mol(atoms, spin, charge=0, verbose=0):
    from pyscf import gto
    astr = "; ".join(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}" for s, p in atoms)
    return gto.M(atom=astr, basis=BASIS, charge=charge, spin=spin,
                 symmetry=False, verbose=verbose)


# ════════════════════════════════════════════════════════════════════════
#  UKS: level shift -> newton fallback -> цикл внутренней стабильности
# ════════════════════════════════════════════════════════════════════════
def _internal_stability(mf):
    """-> (mo_i, stable_i); терпимо к разным версиям API pyscf."""
    try:
        out = mf.stability(internal=True, external=False, return_status=True)
    except TypeError:                                     # нет return_status
        out = mf.stability(internal=True, external=False)
        mo_i = out[0] if isinstance(out, (tuple, list)) else out
        stable = (np.asarray(mo_i).shape == np.asarray(mf.mo_coeff).shape
                  and np.allclose(np.asarray(mo_i), np.asarray(mf.mo_coeff)))
        return mo_i, bool(stable)
    if isinstance(out, (tuple, list)) and len(out) == 4:  # (mo_i, mo_e, st_i, st_e)
        return out[0], bool(out[2])
    if isinstance(out, (tuple, list)) and len(out) == 2:
        return out[0], bool(out[1])
    raise RuntimeError(f"unexpected stability() return: {type(out)}")


def _gaps_ev(mf):
    gaps = []
    for s in (0, 1):
        e = np.asarray(mf.mo_energy[s])
        occ = np.asarray(mf.mo_occ[s])
        if (occ > 0).any() and (occ == 0).any():
            gaps.append(float((e[occ == 0].min() - e[occ > 0].max()) * HARTREE_EV))
    return round(min(gaps), 3) if gaps else None


def uks_stable(mol, tag="", verbose=0):
    """UKS(PBE)/DF: level shift -> newton fallback -> internal-stability loop."""
    from pyscf import dft
    t0 = time.time()
    mf = dft.UKS(mol).density_fit()
    mf.xc = XC
    mf.verbose = verbose
    mf.conv_tol = CONV_TOL_SCF
    mf.level_shift = LEVEL_SHIFT
    mf.max_cycle = 200
    mf.kernel()
    used_newton = False
    if not mf.converged:
        say(f"    [{tag}] plain SCF (level shift {LEVEL_SHIFT}) NOT converged "
            f"-> newton fallback")
        used_newton = True
        mf.level_shift = 0.0
        mfn = mf.newton()
        mfn.verbose = verbose
        mfn.max_cycle = 80
        mfn.kernel(mf.mo_coeff, mf.mo_occ)
        mf.mo_coeff, mf.mo_energy, mf.mo_occ = mfn.mo_coeff, mfn.mo_energy, mfn.mo_occ
        mf.e_tot, mf.converged = float(mfn.e_tot), bool(mfn.converged)
    n_rot, stable = 0, False
    for _ in range(MAX_STAB_LOOPS + 1):
        mo_i, stable = _internal_stability(mf)
        if stable:
            break
        if n_rot >= MAX_STAB_LOOPS:
            break
        n_rot += 1
        used_newton = True
        say(f"    [{tag}] SCF internally UNSTABLE -> rotate + re-converge (#{n_rot})")
        mf.level_shift = 0.0
        mfn = mf.newton()
        mfn.verbose = verbose
        mfn.max_cycle = 80
        mfn.kernel(mo_i, mf.mo_occ)
        mf.mo_coeff, mf.mo_energy, mf.mo_occ = mfn.mo_coeff, mfn.mo_energy, mfn.mo_occ
        mf.e_tot, mf.converged = float(mfn.e_tot), bool(mfn.converged)
    try:
        ss, mult = mf.spin_square()
    except Exception:
        ss, mult = None, None
    info = dict(e_scf=round(float(mf.e_tot), 6),
                scf_converged=bool(mf.converged),
                internally_stable=bool(stable),
                stability_rotations=int(n_rot),
                used_newton=bool(used_newton),
                s2=(round(float(ss), 4) if ss is not None else None),
                gap_eV=_gaps_ev(mf),
                t_scf_s=round(time.time() - t0, 1))
    say(f"    [{tag}] E_UKS={mf.e_tot:.6f}  conv={info['scf_converged']} "
        f"stable={info['internally_stable']} rot={n_rot} <S^2>={info['s2']} "
        f"gap={info['gap_eV']} eV ({info['t_scf_s']}s)")
    return mf, info


def spin_pops(mf, meta=None):
    """Малликен-спины по атомам (локализация дырки Co(IV)=O). Может отказать —
    тогда честно пишем причину."""
    try:
        mol = mf.mol
        dm = mf.make_rdm1()
        s1e = mol.intor_symmetric("int1e_ovlp")
        spin_dm = dm[0] - dm[1]
        pop = np.einsum("ij,ji->i", spin_dm, s1e)
        per_atom = np.zeros(mol.natm)
        for i, lab in enumerate(mol.ao_labels(fmt=None)):
            per_atom[lab[0]] += pop[i]
        out = {f"{i}:{mol.atom_symbol(i)}": round(float(x), 3)
               for i, x in enumerate(per_atom) if abs(x) > 0.05}
        return out
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}


# ════════════════════════════════════════════════════════════════════════
#  UNO-мост -> AVAS с капом -> CASCI (fix_spin) в фиксированном пространстве
# ════════════════════════════════════════════════════════════════════════
def uno_reference(mf):
    """UHF/UKS natural orbitals -> фиктивный ROHF-контейнер + N_u(UNO)."""
    from pyscf import scf
    mol = mf.mol
    dm = mf.make_rdm1()
    dm_t = dm[0] + dm[1]
    s1e = mol.intor_symmetric("int1e_ovlp")
    w, v = np.linalg.eigh(s1e)
    s_half = (v * np.sqrt(w)) @ v.T
    s_mhalf = (v / np.sqrt(w)) @ v.T
    n, u = np.linalg.eigh(s_half @ dm_t @ s_half)
    idx = np.argsort(n)[::-1]
    n = np.clip(n[idx], 0.0, 2.0)
    c_no = s_mhalf @ u[:, idx]
    na, nb = mol.nelec
    occ = np.zeros_like(n)
    occ[:nb] = 2.0
    occ[nb:na] = 1.0
    mf_ro = scf.ROHF(mol)
    mf_ro.mo_coeff = c_no
    mf_ro.mo_occ = occ
    mf_ro.mo_energy = -n            # сортировочный суррогат, НЕ энергии
    n_u_uno = float(np.sum(n * (2.0 - n)))
    frontier = [round(float(x), 4) for x in n if 0.02 < x < 1.98]
    return mf_ro, dict(n_u_uno=round(n_u_uno, 4),
                       uno_fractional_occ=frontier[:20])


def avas_capped(mf_ro, labels, tag=""):
    """AVAS с эскалацией порога до попадания в кап (<=CAS_MAX_ELEC e, CAS_MAX_ORB o)."""
    from pyscf.mcscf import avas
    tried = []
    for thr in AVAS_THRESHOLDS:
        ncas, nelec, mo = avas.avas(mf_ro, labels, threshold=thr,
                                    canonicalize=False, openshell_option=2)
        ncas, nelec = int(ncas), int(nelec)
        tried.append({"threshold": thr, "ncas": ncas, "nelecas": nelec})
        if 0 < ncas <= CAS_MAX_ORB and 0 < nelec <= CAS_MAX_ELEC:
            say(f"    [{tag}] AVAS thr={thr}: CAS({nelec}e,{ncas}o) — в капе")
            return dict(ncas=ncas, nelecas=nelec, threshold=thr,
                        labels=labels, escalation=tried), mo
        say(f"    [{tag}] AVAS thr={thr}: CAS({nelec}e,{ncas}o) — вне капа "
            f"(<= {CAS_MAX_ELEC}e,{CAS_MAX_ORB}o), поднимаем порог")
    return dict(error="AVAS did not fit the cap", escalation=tried,
                labels=labels), None


def run_casci(mf_ro, ncas, nelec_act, spin2, mo, tag="", verbose=0):
    """CASCI (без орбитальной оптимизации) в фиксированном пространстве."""
    from pyscf import mcscf, fci
    t0 = time.time()
    na = (nelec_act + spin2) // 2
    nb = nelec_act - na
    if nb < 0 or na > ncas or nb > ncas or (nelec_act - spin2) % 2:
        return {"error": f"bad active split: {nelec_act}e 2S={spin2} in {ncas}o"}
    mc = mcscf.CASCI(mf_ro, ncas, (na, nb))
    mc.verbose = verbose
    s_val = spin2 / 2.0
    try:
        fci.addons.fix_spin_(mc.fcisolver, ss=s_val * (s_val + 1.0))
    except Exception as exc:
        say(f"    [{tag}] fix_spin_ unavailable: {exc!r}")
    mc.fcisolver.max_cycle = 400
    mc.kernel(mo)
    dm1 = mc.fcisolver.make_rdm1(mc.ci, ncas, (na, nb))
    noon = np.clip(np.linalg.eigvalsh(dm1)[::-1], 0.0, 2.0)
    n_u = float(np.sum(noon * (2.0 - noon)))
    try:
        ss = float(mc.fcisolver.spin_square(mc.ci, ncas, (na, nb))[0])
    except Exception:
        ss = None
    out = dict(e_casci=round(float(mc.e_tot), 6),
               casci_converged=bool(mc.converged),
               ncas=ncas, nelecas=nelec_act, na=na, nb=nb,
               spin2=spin2,
               s2=(round(ss, 4) if ss is not None else None),
               noon=[round(float(x), 4) for x in noon],
               n_u=round(n_u, 4),
               n_u_excess=round(n_u - spin2, 4),
               t_casci_s=round(time.time() - t0, 1))
    say(f"    [{tag}] CASCI({nelec_act}e,{ncas}o) 2S={spin2}: "
        f"E={mc.e_tot:.6f} conv={out['casci_converged']} "
        f"N_u={n_u:.3f} (избыток {out['n_u_excess']:+.3f}) "
        f"<S^2>={out['s2']} ({out['t_casci_s']}s)")
    return out


def mr_verdict(n_u_excess):
    """Качественная эвристика (порог — методическое соглашение, не физконстанта):
    избыток N_u сверх формальных 2S."""
    if n_u_excess is None:
        return "нет данных"
    if n_u_excess > 0.4:
        return "выраженно мультиреференсен (избыток N_u > 0.4)"
    if n_u_excess > 0.15:
        return "умеренно мультиреференсен (избыток N_u 0.15-0.4)"
    return "по существу одно-детерминантен (избыток N_u <= 0.15)"


# ════════════════════════════════════════════════════════════════════════
#  Полный прогон
# ════════════════════════════════════════════════════════════════════════
def avas_site_labels(meta, geom_name):
    """Метки AVAS для оксо-узла: Co 3d + O 2p (оксо/OH-O) + O 2p трёх
    кубановых O при этом Co (site-selective, по индексам атомов)."""
    ico = meta["oxo_co"] if geom_name == "oxidized" else meta["co"][0]
    io = meta["oxo_o"] if geom_name == "oxidized" else meta["lig_o"][meta["co"][0]]
    labels = [f"{ico} Co 3d", f"{io} O 2p"]
    labels += [f"{j} O 2p" for j in meta["cubane_o_of_co"][ico]]
    return labels, ico, io


def _load_chk(path):
    if not os.path.exists(path):
        return None
    with np.load(path) as z:
        return {k: z[k] for k in z.files}


def _mf_from_chk(mol, chk, verbose=0):
    from pyscf import dft
    mf = dft.UKS(mol).density_fit()
    mf.xc = XC
    mf.verbose = verbose
    mf.mo_coeff = (chk["mo_a"], chk["mo_b"])
    mf.mo_occ = (chk["occ_a"], chk["occ_b"])
    mf.converged = bool(chk["converged"])
    return mf


def run_geometry(geom_name, res, args):
    atoms, meta = build_cubane(oxidized=(geom_name == "oxidized"))
    gsec = res.setdefault("geometry", {})
    if geom_name not in gsec:
        rep = geom_report(atoms, meta)
        rep["honesty"] = GEOM_HONESTY
        gsec[geom_name] = rep
        atomic_save(res, OUT_RESULTS)
        say(f"  [{geom_name}] geometry: {rep['formula']} ({rep['natoms']} at), "
            f"min nonbonded {rep['min_nonbonded_A']} A "
            f"{rep['min_nonbonded_pair']}")

    entry = res.setdefault(geom_name, {})
    states = [(st, s2) for g, st, s2 in STATES if g == geom_name]
    mfs = {}
    for st, s2 in states:
        stag = f"{geom_name}-{st}"
        sec = entry.setdefault(st, {})
        chk_path = CHK_TPL.format(state=stag)
        chk = _load_chk(chk_path)
        if sec.get("e_scf") is not None and chk is not None:
            say(f"  [{stag}] SCF уже сделан (restart) — из чекпойнта")
            mfs[st] = _mf_from_chk(make_mol(atoms, spin=s2), chk, args.verbose)
            continue
        say(f"  [{stag}] UKS-{XC.upper()}/{BASIS}/DF, 2S={s2} (cold start)")
        mol = make_mol(atoms, spin=s2, verbose=args.verbose)
        if st == states[0][0]:
            say(f"    [{stag}] nao={mol.nao} nelec={mol.nelectron}")
        mf, sinfo = uks_stable(mol, tag=stag, verbose=args.verbose)
        sinfo["mulliken_spin"] = spin_pops(mf, meta)
        sec.update(sinfo)
        atomic_save(res, OUT_RESULTS)
        np.savez_compressed(chk_path,
                            mo_a=np.asarray(mf.mo_coeff[0]),
                            mo_b=np.asarray(mf.mo_coeff[1]),
                            occ_a=np.asarray(mf.mo_occ[0]),
                            occ_b=np.asarray(mf.mo_occ[1]),
                            converged=np.bool_(mf.converged))
        say(f"  [{stag}] checkpoint -> {os.path.basename(chk_path)}")
        mfs[st] = mf

    # UKS спин-зазор на одной (идеализированной!) геометрии
    st0, st1 = states[0][0], states[1][0]
    if entry.get(st0, {}).get("e_scf") is not None and \
            entry.get(st1, {}).get("e_scf") is not None:
        entry["uks_gap_eV"] = {
            f"{st1}_minus_{st0}": round(
                (entry[st1]["e_scf"] - entry[st0]["e_scf"]) * HARTREE_EV, 3),
            "caveat": "идеализированная single-point геометрия — качественно"}
        atomic_save(res, OUT_RESULTS)

    # ---- CASCI-зонд: пространство от ОПОРНОГО состояния, оба спина в нём ----
    ref = CAS_REFERENCE[geom_name]
    csec = entry.setdefault("casci", {})
    need = [st for st, _ in states if csec.get(st, {}).get("n_u") is None
            and "error" not in csec.get(st, {})]
    if need:
        labels, ico, io = avas_site_labels(meta, geom_name)
        say(f"  [{geom_name}] UNO({ref}) -> AVAS{labels} -> CASCI (оба спина "
            f"в ОДНОМ пространстве)")
        mf_ref = mfs[ref]
        mf_ro, uno_info = uno_reference(mf_ref)
        csec["uno"] = dict(reference_state=ref, **uno_info)
        cas, mo = avas_capped(mf_ro, labels, tag=geom_name)
        csec["avas"] = cas
        csec["site_atoms"] = {"co": f"{ico}:Co", "o": f"{io}:O"}
        atomic_save(res, OUT_RESULTS)
        if mo is None:
            say(f"  [{geom_name}] AVAS не влез в кап — CASCI честно пропущен")
        else:
            for st, s2 in states:
                if csec.get(st, {}).get("n_u") is not None:
                    continue
                csec[st] = run_casci(mf_ro, cas["ncas"], cas["nelecas"], s2, mo,
                                     tag=f"{geom_name}-{st}", verbose=args.verbose)
                atomic_save(res, OUT_RESULTS)
    else:
        say(f"  [{geom_name}] CASCI уже сделан (restart)")
    return entry


def main_real(args):
    from pyscf import lib
    say(f"oer-scaling Stage 2: усечённый кубан [Co4O4] — threads={lib.num_threads()}")
    res = json.load(open(OUT_RESULTS)) if os.path.exists(OUT_RESULTS) else {}
    res.setdefault("meta", dict(
        track="oer-scaling", stage=2,
        system=("truncated cubane core: resting [Co4O4(OH)4(NH3)4]^0 (Co(III)4) "
                "and oxidized [Co4O4(OH)3(O)(NH3)4]^0 (formally one Co(IV)=O "
                "after PCET H removal); OH-/NH3 are surrogates for acetate/"
                "pyridine of Co4O4(OAc)4(py)4 (JACS 2011, 10.1021/ja203877y)"),
        method=(f"UKS-{XC.upper()}(level shift {LEVEL_SHIFT} -> newton fallback "
                f"-> internal-stability loop)/{BASIS}/DF; UNO bridge; "
                f"site-AVAS(Co 3d + oxo/OH O 2p + 3 cubane O 2p, threshold "
                f"escalation, cap ({CAS_MAX_ELEC}e,{CAS_MAX_ORB}o), "
                f"openshell_option=2); CASCI fix_spin, active space fixed once "
                f"per geometry on the reference spin state (resting->singlet, "
                f"oxidized->doublet), both spins solved in the SAME space"),
        metric="N_u = sum n_i(2-n_i) over CASCI NOON; excess = N_u - 2S",
        geometry_honesty=GEOM_HONESTY,
        caveats=("v1 single-point, no relaxation; cluster != electrode (no "
                 "solvent/potential) -> no eta/dG here, only N_u; E(CASCI) on "
                 "HF core over UNO(UKS) orbitals — NOT comparable with E(UKS); "
                 "'one Co(IV)=O' is a formal label — hole may delocalize, see "
                 "mulliken_spin; 32 atoms slightly exceed the ~25 estimate of "
                 "the Stage 1 note"),
        started=time.strftime("%Y-%m-%d %H:%M:%S")))
    atomic_save(res, OUT_RESULTS)

    for geom_name in ("resting", "oxidized"):
        if args.geoms and geom_name not in args.geoms:
            continue
        say(f"--- geometry {geom_name} ---")
        run_geometry(geom_name, res, args)

    # ---- вердикт по оксо-узлу ----
    ox = res.get("oxidized", {}).get("casci", {})
    verdict = {}
    for st in ("doublet", "quartet"):
        exc = ox.get(st, {}).get("n_u_excess")
        if exc is not None:
            verdict[st] = {"n_u": ox[st]["n_u"], "n_u_excess": exc,
                           "verdict": mr_verdict(exc)}
    rest = res.get("resting", {}).get("casci", {}).get("singlet", {})
    if rest.get("n_u") is not None:
        verdict["resting_singlet_baseline"] = {
            "n_u": rest["n_u"], "n_u_excess": rest["n_u_excess"],
            "verdict": mr_verdict(rest["n_u_excess"])}
    if verdict:
        verdict["note"] = ("качественный вопрос Stage 2: мультиреференсен ли "
                           "узел Co(IV)=O; пороги 0.15/0.4 — методическая "
                           "эвристика, единая для всех треков дневника")
        res["verdict"] = verdict
        atomic_save(res, OUT_RESULTS)
        say(f"VERDICT: {json.dumps(verdict, ensure_ascii=False)}")
    say(f"saved -> {OUT_RESULTS}; done.")


# ════════════════════════════════════════════════════════════════════════
#  Дым-тест: сборка геометрий + xyz в JSON + мини-SCF/UNO/AVAS/CASCI на CoO+
# ════════════════════════════════════════════════════════════════════════
def main_smoke(args):
    from pyscf import gto, lib
    t0 = time.time()
    out_path = OUT_RESULTS.replace(".json", "_smoke.json")
    say(f"SMOKE: geometry build + xyz->JSON + mini UKS/UNO/AVAS/CASCI on CoO+ "
        f"fragment; threads={lib.num_threads()}; out={out_path}")
    res = {"meta": dict(mode="smoke",
                        note=("plumbing test only: CoO+ fragment spin=4 is an "
                              "arbitrary high-spin choice for robust SCF, NO "
                              "ground-state claim; r(Co-O)=1.68 A idealized"),
                        started=time.strftime("%Y-%m-%d %H:%M:%S"))}
    ok = True

    # 1) обе геометрии + xyz в JSON + проверка клэшей
    res["geometry"] = {}
    for geom_name, exp_atoms in (("resting", 32), ("oxidized", 31)):
        atoms, meta = build_cubane(oxidized=(geom_name == "oxidized"))
        rep = geom_report(atoms, meta)
        rep["honesty"] = GEOM_HONESTY
        labels, ico, io = avas_site_labels(meta, geom_name)
        rep["avas_labels_planned"] = labels
        mol = make_mol(atoms, spin=(0 if geom_name == "resting" else 1))
        rep["nao"] = int(mol.nao)
        rep["nelectron"] = int(mol.nelectron)
        res["geometry"][geom_name] = rep
        atomic_save(res, out_path)
        good = (rep["natoms"] == exp_atoms and rep["min_nonbonded_A"] > 1.5)
        ok &= good
        say(f"  [{geom_name}] {rep['formula']} natoms={rep['natoms']} "
            f"(exp {exp_atoms}) nao={rep['nao']} nelec={rep['nelectron']} "
            f"min_nonbonded={rep['min_nonbonded_A']} A "
            f"{rep['min_nonbonded_pair']} -> {'OK' if good else 'FAIL'}")

    # 2) мини-SCF + полный CASCI-конвейер на фрагменте CoO+ (2 атома)
    say("  [fragment CoO+] UKS + UNO -> AVAS -> CASCI (полный конвейер в миниатюре)")
    frag = gto.M(atom=f"Co 0 0 0; O 0 0 {R_CO_OXO}", basis=BASIS, charge=1,
                 spin=4, symmetry=False, verbose=0)
    mf, sinfo = uks_stable(frag, tag="CoO+", verbose=args.verbose)
    sinfo["mulliken_spin"] = spin_pops(mf)
    res["fragment_coo_plus"] = {"scf": sinfo}
    atomic_save(res, out_path)
    mf_ro, uno_info = uno_reference(mf)
    res["fragment_coo_plus"]["uno"] = uno_info
    cas, mo = avas_capped(mf_ro, ["Co 3d", "O 2p"], tag="CoO+")
    res["fragment_coo_plus"]["avas"] = cas
    atomic_save(res, out_path)
    if mo is None:
        ok = False
    else:
        ci = run_casci(mf_ro, cas["ncas"], cas["nelecas"], 4, mo,
                       tag="CoO+ 2S=4", verbose=args.verbose)
        res["fragment_coo_plus"]["casci"] = ci
        ok &= ("n_u" in ci)
        atomic_save(res, out_path)

    res["smoke_pass"] = bool(ok)
    res["t_total_s"] = round(time.time() - t0, 1)
    atomic_save(res, out_path)
    say(f"SMOKE {'PASS' if ok else 'FAIL'} in {res['t_total_s']:.0f}s -> {out_path}")
    return 0 if ok else 1


def _parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[2])
    p.add_argument("--smoke", action="store_true",
                   help="сборка геометрии + xyz в JSON + мини-SCF/CASCI на CoO+ "
                        "(<=3 мин, 1 ядро); пишет в *_smoke.json")
    p.add_argument("--geoms", default="",
                   help="подмножество геометрий: resting,oxidized (по умолч. обе)")
    p.add_argument("--verbose", type=int, default=0,
                   help="pyscf verbose (default 0)")
    a = p.parse_args(argv)
    a.geoms = [g.strip() for g in a.geoms.split(",") if g.strip()]
    return a


if __name__ == "__main__":
    _args = _parse_args()
    sys.exit(main_smoke(_args) if _args.smoke else main_real(_args))
