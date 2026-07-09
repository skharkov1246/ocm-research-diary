#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/ocm_mnw_selectivity.py — Этап 15 OCM: закрыть двухсубстратный ΔΔE‡
на реальном радикально-кислородном центре Mn=O (Mn–Na₂WO₄/SiO₂, минимальная модель).

Чинит ровно то, на чём упал Этап 14:
  • ЕДИНОЕ курированное активное пространство CAS(9e,10o) для всех 4 структур
    (R/TS × CH₄/C₂H₆): AVAS(Mn 3d + O 2p + 1s переносимого H). Проверено:
    размер стабилен к порогу 0.05–0.2 (метки с C 2p «прыгают» 13o↔15o — не берём).
  • Настоящий HAT-профиль: 2-координатный пин r(C–H)+r(O–H) (geomeTRIC),
    каркас релаксируется; TS = внутренний максимум профиля (краевой — флагуется).
  • Оба субстрата до конца: CH₄ И C₂H₆, один протокол, одна сетка.

Протокол: UKS-PBE0/def2-SVP геометрия и DFT-барьеры; ROHF→AVAS→CASSCF(9e,10o)
(робастный: newton-SCF старт, натурбитали, fix_spin)→SC-NEVPT2 — энергетика.
Секстет (S=5/2) по всему пути (спин-сохраняющий HAT). Модель — диатомный MnO
(секстет ⁶Σ⁺ — известное основное состояние) как минимальный стенд-ин
рабочего Mn=O; НЕ периодический Mn–Na₂WO₄/SiO₂ (оговорка в дневнике).

Дескриптор (конвенция «карты потолка»): ΔΔE‡ = барьер(C₂H₆) − барьер(CH₄);
< 0 — центр доокисляет C₂-продукт быстрее метана (стена BEP), для 70% нужно ≳ +4.

Запуск (стадии независимы, результаты инкрементально в JSON):
  python3 -u routes/ocm_mnw_selectivity.py geom   ch4|c2h6   # скан-профиль
  python3 -u routes/ocm_mnw_selectivity.py polish ch4|c2h6   # цепочка warm-start+стабильность
  python3 -u routes/ocm_mnw_selectivity.py energy ch4|c2h6   # CASSCF+NEVPT2 в R и TS
  python3 -u routes/ocm_mnw_selectivity.py merge             # итоговый results-JSON
"""
import json
import math
import os
import sys
import time

import numpy as np
from pyscf import dft, fci, gto, mcscf, scf
from pyscf.mcscf import avas
from pyscf.mrpt import NEVPT

HARTREE_KCAL = 627.509474
DIR = os.path.dirname(os.path.abspath(__file__))
BASIS = "def2-svp"
MODEL = os.environ.get("OCM_MODEL", "bare")   # bare (MnO) | embedded (Na/W-оболочка)
SPIN = int(os.environ.get("OCM_SPIN", "5"))   # 2S = число SOMO (bare-секстет: 5)
# Этап 17: активный оксо-металл параметризован — Mn (базовый) | Cr | Fe (окно
# Cr/Mn/Fe вакансионного волкана). Свап меняет чётность/основной спин: сначала
# гоняем стадию `spins`, читаем ground_2S, выставляем OCM_SPIN, потом пайплайн.
METAL = os.environ.get("OCM_METAL", "Mn")
XC = "pbe0"
_defpref = "ocm_mnw" if MODEL == "bare" else "ocm_mnw_emb"
if METAL != "Mn":
    _defpref += "_" + METAL.lower()          # ocm_mnw_emb_cr / _fe — не пересекается
PREF = os.environ.get("OCM_PREFIX", _defpref)
# единое активное пространство фиксированного СОСТАВА: 2 закрытые пары + все SOMO
# + 3 виртуали => (ne,no) зависит только от спина (bare-секстет: CAS(9e,10o))
N_DOCC, N_VIR = 2, 3
CAS_TARGET = (2 * N_DOCC + SPIN, N_DOCC + SPIN + N_VIR)
# валентная d-оболочка металла: 3d-ряд по умолчанию; Этап 19 добавляет 4d-ряд
# (Zr/Nb — гомологи Ti/V; def2-SVP несёт для них ECP-28, ecp=BASIS уже передаётся)
_DSHELL = {"Zr": "4d", "Nb": "4d", "Mo": "4d"}
AVAS_LABELS = [f"{METAL} {_DSHELL.get(METAL, '3d')}",
               "O 2p", "2 H 1s"]   # атом 2 (0-based) — переносимый H
# версия геометрического протокола: при несовпадении resume пересчитывает точки
# (frozen-frame сменил floppy-relax — старые точки несовместимы, отбрасываются)
GEOM_VERSION = "embedded-frozen-frame-v2" if MODEL == "embedded" else "bare-v1"

# HAT-сетка: пары (r_CH, r_OH), Å. R-точка: пин только r(O–H)=2.2 (пре-комплекс).
R_POINT = (None, 2.20)
# bare — плотная 8-точечная сетка (Этап 15); embedded — грубая 5-точечная
# (R + подход + два вокруг TS + продукт): точка эмбеддед-кластера в разы дороже
# (R-точка ≈ 8 ч), а дескриптор ΔΔE‡ требует лишь профиль, брекетящий TS.
# ПЕРВЫЕ ДВА индекса грубой сетки совпадают с плотной — уже посчитанные точки
# (в т.ч. дорогая R) переиспользуются resume-логикой, не пересчитываются.
if MODEL == "embedded":
    # первые 3 индекса (R,1,2) совпадают с плотной сеткой — resume-совместимо
    # с уже начатым 8-точечным прогоном (какие точки тот ни успел сохранить)
    PATH = [(1.15, 1.70), (1.20, 1.50), (1.35, 1.17), (1.50, 1.02)]
    # Этап 18 (OCM_GRID=ext): у Ti/V профиль растёт до края 5-точечной сетки
    # (TS не пойман — ts_interior=false). Продолжаем путь в продуктовую сторону
    # АППЕНДОМ (индексы 0-4 не тронуты ⇒ resume переиспользует посчитанные
    # точки, GEOM_VERSION тот же): O–H дожимается к равновесию ~0.96 Å,
    # C–H рвётся дальше. Перевал обязан оказаться интерьерным.
    if os.environ.get("OCM_GRID", "") == "ext":
        PATH += [(1.60, 0.99), (1.75, 0.97), (1.90, 0.96), (2.10, 0.95)]
else:
    PATH = [(1.15, 1.70), (1.20, 1.50), (1.25, 1.35), (1.30, 1.25),
            (1.35, 1.17), (1.40, 1.10), (1.50, 1.02)]


# ---------------------------------------------------------------- геометрии
def start_atoms(sub):
    """Стартовая геометрия Mn=O ... H–R (коллинеарная HAT-ось по z).
    Порядок атомов фиксирован: Mn(0), O(1), H_перенос(2), C(3), ..."""
    zMn, zO = -1.65, 0.0
    rOH0, rCH0 = 2.20, 1.10
    zH = rOH0
    zC = rOH0 + rCH0
    d, th = 1.09, math.radians(108.0)
    atoms = [(METAL, (0.0, 0.0, zMn)), ("O", (0.0, 0.0, zO)),
             ("H", (0.0, 0.0, zH)), ("C", (0.0, 0.0, zC))]
    if sub == "ch4":
        for k in range(3):
            phi = math.radians(120.0 * k)
            atoms.append(("H", (d * math.sin(th) * math.cos(phi),
                                d * math.sin(th) * math.sin(phi),
                                zC - d * math.cos(th))))
    elif sub == "c2h6":
        rCC = 1.53
        zC2 = zC + rCC * math.cos(math.radians(35.0))
        xC2 = rCC * math.sin(math.radians(35.0))
        atoms.append(("C", (xC2, 0.0, zC2)))
        for k in (1, 2):   # два оставшихся H на C1
            phi = math.radians(120.0 * k + 60.0)
            atoms.append(("H", (d * math.sin(th) * math.cos(phi),
                                d * math.sin(th) * math.sin(phi),
                                zC - d * math.cos(th))))
        for k in range(3):  # H на C2
            phi = math.radians(120.0 * k)
            atoms.append(("H", (xC2 + d * math.sin(th) * math.cos(phi),
                                d * math.sin(th) * math.sin(phi),
                                zC2 + d * math.cos(th))))
    else:
        raise ValueError(sub)
    if MODEL == "embedded":
        atoms += environment_atoms()
    return atoms


def environment_atoms():
    """Первая Na/W-оболочка рабочего центра Mn–Na₂WO₄/SiO₂ (минимальный
    нейтральный мотив): [O=Mn(–O–Na)(μ-O)₂W(=O)(OH)]. Формальные степени:
    Mn(IV), W(VI); заряд 0. Стартовые длины — типовые (Mn=O 1.60 уже в оси,
    Mn–O 1.8–1.9, W–O 1.75–1.95, Na–O ~2.2); всё релаксируется в R-точке.
    Позиции относительно Mn (см. start_atoms: Mn в (0,0,zMn), oxo выше по z)."""
    zMn = -1.65
    def rel(dx, dy, dz):
        return (dx, dy, zMn + dz)
    return [
        ("O",  rel(-1.75, 0.0, 0.45)),    # O(–Na), мостик заряда
        ("Na", rel(-3.05, 0.0, 1.65)),
        ("O",  rel(1.10, 1.25, -0.80)),   # μ-O(1)
        ("O",  rel(1.10, -1.25, -0.80)),  # μ-O(2)
        ("W",  rel(2.65, 0.0, -1.45)),
        ("O",  rel(3.55, 0.0, -2.95)),    # W=O
        ("O",  rel(4.10, 0.0, -0.35)),    # W–O(H)
        ("H",  rel(4.35, 0.85, 0.05)),
    ]


def stage_spins():
    """Лестница спинов UKS-PBE0 на стартовой R-геометрии (кластер + CH₄):
    выбираем основное спиновое состояние ПЕРЕД пайплайном, не постулируем.
    Resume: если лестница уже посчитана (есть ground_2S) — пропускаем (иначе
    каждое поколение впустую гоняло бы 4 UKS на 201-AO кластере, часы)."""
    global SPIN
    path = os.path.join(DIR, f"{PREF}_spins.json")
    if os.path.exists(path):
        try:
            done = json.load(open(path))
            if "ground_2S" in done:
                print(f"[spins] resume: ground 2S={done['ground_2S']} on disk "
                      f"(ladder {done.get('rel_kcal')}) — skip recompute", flush=True)
                return
        except Exception:
            pass
    atoms = start_atoms("ch4")
    # чётность электронов без пробной сборки (ECP-остов чётный, заряд 0)
    parity = sum(gto.charge(a[0]) for a in atoms) % 2
    cands = [s for s in range(parity, 8, 2)][:4]   # 4 нижних мультиплетности
    out = {"model": MODEL, "candidates_2S": cands, "ladder": []}
    path = os.path.join(DIR, f"{PREF}_spins.json")
    for s2 in cands:
        t0 = time.time()
        try:
            mol = gto.M(atom=atoms, basis=BASIS, ecp=BASIS, charge=0, spin=s2,
                        verbose=0, max_memory=6000)
            mf = dft.UKS(mol).density_fit()
            mf.xc = XC
            mf.conv_tol = 1e-7
            mfn = scf.newton(mf)
            mfn.kernel()
            out["ladder"].append({"spin_2S": s2, "e_h": float(mfn.e_tot),
                                  "converged": bool(mfn.converged),
                                  "spin_square": round(float(mfn.spin_square()[0]), 3),
                                  "wall_s": round(time.time() - t0, 1)})
        except Exception as e:
            out["ladder"].append({"spin_2S": s2, "error": str(e)[:200]})
        with open(path, "w") as f:
            json.dump(out, f, indent=1)
        print(f"[spins] 2S={s2}: {out['ladder'][-1]}", flush=True)
    ok = [x for x in out["ladder"] if x.get("converged")]
    if ok:
        best = min(ok, key=lambda x: x["e_h"])
        out["ground_2S"] = best["spin_2S"]
        rel = [(x["e_h"] - best["e_h"]) * HARTREE_KCAL for x in ok]
        out["rel_kcal"] = {str(x["spin_2S"]): round(r, 2) for x, r in zip(ok, rel)}
        print(f"[spins] ground 2S={best['spin_2S']}; ladder {out['rel_kcal']}",
              flush=True)
    with open(path, "w") as f:
        json.dump(out, f, indent=1)


def build_mol(atoms):
    return gto.M(atom=atoms, basis=BASIS, ecp=BASIS, charge=0, spin=SPIN,
                 verbose=0, max_memory=6000)


def uks(mol, dm0=None):
    """UKS-PBE0 (RI/DF — кратно быстрее на градиентах), newton-primary
    (устойчив на Mn), DIIS-фоллбэк."""
    mf = dft.UKS(mol).density_fit()
    mf.xc = XC
    mf.conv_tol = 1e-8
    mfn = scf.newton(mf)
    mfn.kernel(dm0=dm0) if dm0 is not None else mfn.kernel()
    if not mfn.converged:
        mf2 = dft.UKS(mol).density_fit()
        mf2.xc = XC
        mf2.level_shift = 0.3
        mf2.max_cycle = 300
        mf2.kernel(dm0=dm0)
        mfn = scf.newton(mf2)
        mfn.kernel(mf2.make_rdm1())
    return mfn


def constrained_opt(atoms, rCH, rOH, tag, dm0=None, maxsteps=200):
    """Релакс-оптимизация с пином r(C–H) и/или r(O–H) (geomeTRIC, $set).
    dm0 — плотность предыдущей точки (тёплый старт: быстрее и держит одну
    электронную ветку вдоль скана — урок polish Этапа 15). Возвращает также
    итоговую плотность для цепочки."""
    from pyscf.geomopt.geometric_solver import optimize
    # имя файла ограничений неймспейсим ПРЕФИКСОМ (=металлом): при параллельном
    # прогоне нескольких металлов в одной директории общий `_constr_ch4_0.txt`
    # вызывал гонку (один процесс os.remove — у другого FileNotFoundError).
    cfile = os.path.join(DIR, f"_constr_{PREF}_{tag}.txt")
    lines = ["$set"]
    if rOH is not None:
        lines.append(f"distance 2 3 {rOH:.4f}")   # O(2)–H(3), 1-based
    if rCH is not None:
        lines.append(f"distance 3 4 {rCH:.4f}")   # H(3)–C(4)
    if MODEL == "embedded":
        # embedded-cluster приближение: первая Na/W-оболочка удерживается решёткой —
        # замораживаем её (последние len(env) атомов), релаксируем только реактивный
        # центр Mn=O + переносимый H + субстрат. Быстрее (30-50 шагов вместо 200+)
        # И физичнее (каркас не «плывёт» в газовой фазе); каркас идентичен во всех
        # точках → согласованный профиль. Оговорка: каркас при идеализированной
        # геометрии, не крист. координатах.
        n = len(atoms)
        n_env = len(environment_atoms())
        lines.append("$freeze")
        lines.append(f"xyz {n - n_env + 1}-{n}")   # 1-based, последние n_env атомов
    with open(cfile, "w") as f:
        f.write("\n".join(lines) + "\n")
    mf = uks(build_mol(atoms), dm0=dm0)
    if MODEL == "embedded":   # больше DOF: щадящие критерии (релакс-скан, не минимум)
        conv = dict(convergence_energy=2e-6, convergence_grms=4.5e-4,
                    convergence_gmax=9e-4, convergence_drms=2.3e-3,
                    convergence_dmax=3.6e-3)
    else:
        conv = dict(convergence_energy=1e-6, convergence_grms=3e-4,
                    convergence_gmax=6e-4, convergence_drms=1.5e-3,
                    convergence_dmax=2.5e-3)
    # assert_convergence=False: при упоре в maxsteps вернуть последнюю геометрию,
    # а не падать (для скан-точки на floppy-каркасе это приемлемо — не стационар)
    mol_eq = optimize(mf, constraints=cfile, maxsteps=maxsteps,
                      assert_convergence=False, **conv)
    if os.path.exists(cfile):
        os.remove(cfile)
    new_atoms = [(mol_eq.atom_symbol(i), tuple(c))
                 for i, c in enumerate(mol_eq.atom_coords(unit="Angstrom"))]
    mf_final = uks(mol_eq, dm0=mf.make_rdm1())
    ss = mf_final.spin_square()[0]
    return (new_atoms, float(mf_final.e_tot), bool(mf_final.converged),
            float(ss), mf_final.make_rdm1())


def stage_geom(sub):
    fr = "frozen frame" if MODEL == "embedded" else "frame relaxed"
    out = {"substrate": sub, "geom_version": GEOM_VERSION,
           "protocol": f"UKS-{XC.upper()}/{BASIS}, sextet, "
           f"2D-pin r(C-H)+r(O-H) (geomeTRIC), {fr}", "points": []}
    path = os.path.join(DIR, f"{PREF}_{sub}_geom.json")
    atoms = start_atoms(sub)
    grid = [R_POINT] + PATH
    if os.path.exists(path):        # резюме после смерти контейнера
        with open(path) as f:
            saved = json.load(f)
        if saved.get("geom_version") != GEOM_VERSION:
            print(f"[{sub}] geom_version mismatch "
                  f"({saved.get('geom_version')} != {GEOM_VERSION}) — recompute all",
                  flush=True)
        else:
            out["points"] = saved.get("points", [])[:len(grid)]
        if out["points"]:
            atoms = [(s, tuple(c)) for s, c in out["points"][-1]["atoms"]]
            print(f"[{sub}] resume: {len(out['points'])}/{len(grid)} points on disk",
                  flush=True)
    dm = None
    for i, (rCH, rOH) in enumerate(grid):
        if i < len(out["points"]):
            continue
        t0 = time.time()
        if MODEL == "embedded":
            ms = 90 if i == 0 else 130   # ограниченный кран на floppy-каркасе
        else:
            ms = 200
        atoms, e, conv, ss, dm = constrained_opt(
            atoms, rCH, rOH, f"{sub}_{i}", dm0=dm, maxsteps=ms)
        out["points"].append({
            "rCH_pin": rCH, "rOH_pin": rOH, "e_h": e, "scf_converged": conv,
            "spin_square": round(ss, 3), "wall_s": round(time.time() - t0, 1),
            "atoms": [[s, [round(x, 6) for x in c]] for s, c in atoms]})
        with open(path, "w") as f:
            json.dump(out, f, indent=1)
        print(f"[{sub}] point {i}/{len(grid)-1} rCH={rCH} rOH={rOH} "
              f"E={e:.6f} conv={conv} <S2>={ss:.2f} ({out['points'][-1]['wall_s']}s)",
              flush=True)
    es = [p["e_h"] for p in out["points"]]
    rel = [(e - es[0]) * HARTREE_KCAL for e in es]
    imax = int(np.argmax(rel[1:]) + 1)
    out["rel_kcal"] = [round(x, 2) for x in rel]
    out["ts_index"] = imax
    out["ts_interior"] = bool(0 < imax < len(grid) - 1)
    out["dft_barrier_kcal"] = round(rel[imax], 2)
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print(f"[{sub}] DFT barrier {rel[imax]:.2f} kcal/mol at point {imax} "
          f"(interior={out['ts_interior']})", flush=True)


def uks_stable(mol, dm0=None, max_hops=3):
    """UKS + внутренняя стабильность: если решение нестабильно — спуск на
    низшее UKS-решение (лечит перескоки состояний вдоль скана)."""
    mf = uks(mol, dm0=dm0)
    for _ in range(max_hops):
        mo_new = mf.stability()[0]
        if mo_new is mf.mo_coeff or (
                isinstance(mo_new, (tuple, list, np.ndarray))
                and np.allclose(np.asarray(mo_new[0]), np.asarray(mf.mo_coeff[0]))):
            break
        dm = mf.make_rdm1(mo_new, mf.mo_occ)
        mf = uks(mol, dm0=dm)
    return mf


def stage_polish(sub):
    """Полировка профиля: по сохранённым геометриям — цепочка warm-start
    (dm предыдущей точки) + стабильность UKS. Снимает перескоки состояний,
    из-за которых профиль «зубчатый» (артефакт холодных SCF-стартов)."""
    path = os.path.join(DIR, f"{PREF}_{sub}_geom.json")
    with open(path) as f:
        g = json.load(f)
    dm = None
    for i, p in enumerate(g["points"]):
        mol = build_mol([(s, tuple(c)) for s, c in p["atoms"]])
        t0 = time.time()
        mf = uks_stable(mol, dm0=dm)
        dm = mf.make_rdm1()
        p["e_h_raw"] = p["e_h"]
        p["e_h"] = float(mf.e_tot)
        p["spin_square"] = round(float(mf.spin_square()[0]), 3)
        p["scf_converged"] = bool(mf.converged)
        p["polished"] = True
        with open(path, "w") as f:
            json.dump(g, f, indent=1)
        print(f"[{sub}] polish {i}: raw {p['e_h_raw']:.6f} -> {p['e_h']:.6f} "
              f"(d={(p['e_h']-p['e_h_raw'])*HARTREE_KCAL:+.2f} kcal) "
              f"<S2>={p['spin_square']} ({time.time()-t0:.0f}s)", flush=True)
    # обратный проход (против гистерезиса): цепочка dm с конца, берём нижнее
    dm = None
    for i in range(len(g["points"]) - 1, -1, -1):
        p = g["points"][i]
        mol = build_mol([(s, tuple(c)) for s, c in p["atoms"]])
        t0 = time.time()
        mf = uks_stable(mol, dm0=dm)
        dm = mf.make_rdm1()
        if float(mf.e_tot) < p["e_h"] - 1e-6:
            p["e_h"] = float(mf.e_tot)
            p["spin_square"] = round(float(mf.spin_square()[0]), 3)
            p["scf_converged"] = bool(mf.converged)
            print(f"[{sub}] polish-rev {i}: lower -> {p['e_h']:.6f} "
                  f"<S2>={p['spin_square']} ({time.time()-t0:.0f}s)", flush=True)
        with open(path, "w") as f:
            json.dump(g, f, indent=1)
    rel, imin, imax, interior = analyze_scan(g)
    g["rel_kcal"] = [round(x, 2) for x in rel]
    g["r_index"], g["ts_index"], g["ts_interior"] = imin, imax, interior
    g["dft_barrier_kcal"] = round(rel[imax] - rel[imin], 2)
    g["polish_done"] = True
    with open(path, "w") as f:
        json.dump(g, f, indent=1)
    print(f"[{sub}] polished profile: {g['rel_kcal']} "
          f"R={imin} TS={imax} barrier={g['dft_barrier_kcal']} kcal/mol", flush=True)


def analyze_scan(g):
    """R = минимум профиля ДО максимума (истинный пре-комплекс: у CH₄ при
    сближении с O есть яма ниже точки 0), TS = внутренний максимум после него."""
    es = [p["e_h"] for p in g["points"]]
    rel = [(e - es[0]) * HARTREE_KCAL for e in es]
    imax = int(np.argmax(rel[1:]) + 1)
    imin = int(np.argmin(rel[:imax]))
    return rel, imin, imax, bool(0 < imax < len(rel) - 1)


# ---------------------------------------------------------------- энергетика
def rohf(mol):
    mf = scf.ROHF(mol)
    mf.conv_tol = 1e-8
    mfn = scf.newton(mf)
    mfn.kernel()
    if not mfn.converged:
        mf2 = scf.ROHF(mol)
        mf2.level_shift = 0.4
        mf2.max_cycle = 300
        mf2.kernel()
        mfn = scf.newton(mf2)
        mfn.kernel(mf2.make_rdm1())
    return mfn


def forced_avas(mf, n_docc=2, n_vir=3):
    """Fixed-AVAS: активное пространство ФИКСИРОВАННОГО состава на любой
    геометрии — все SOMO (секстет: 5) + топ-n_docc закрытых пар + топ-n_vir
    виртуалей по проекции на референсные АО (Mn 3d, O 2p, 1s переносимого H).
    Гарантирует CAS(9e,10o) конструкцией — подбор порога (Этап 14/первый
    прогон Этапа 15) «прыгал» размером между R и TS."""
    mol = mf.mol
    pmol = mol.copy()
    pmol.basis = "minao"
    pmol.build(False, False)
    baslst = pmol.search_ao_label(AVAS_LABELS)
    s2 = pmol.intor_symmetric("int1e_ovlp")[np.ix_(baslst, baslst)]
    s21 = gto.intor_cross("int1e_ovlp", pmol, mol)[baslst]
    sa = s21.T @ np.linalg.solve(s2, s21)     # метрика проекции в АО-базисе

    occ = np.asarray(mf.mo_occ)
    C = np.asarray(mf.mo_coeff)
    idx_d, idx_s, idx_v = (np.where(occ == 2)[0], np.where(occ == 1)[0],
                           np.where(occ == 0)[0])

    def split(idx, n_keep):
        Cb = C[:, idx]
        w, v = np.linalg.eigh(Cb.T @ sa @ Cb)
        order = np.argsort(w)[::-1]           # по убыванию проекции
        Cr = Cb @ v[:, order]
        return Cr[:, :n_keep], Cr[:, n_keep:], w[order][:n_keep]

    act_d, core_d, w_d = split(idx_d, n_docc)
    act_v, rest_v, w_v = split(idx_v, n_vir)
    mo = np.hstack([core_d, act_d, C[:, idx_s], act_v, rest_v])
    ne = 2 * n_docc + len(idx_s)
    no = n_docc + len(idx_s) + n_vir
    assert (ne, no) == CAS_TARGET, (ne, no)
    weights = [round(float(x), 3) for x in list(w_d) + list(w_v)]
    return mo, weights


def cas_nevpt2(atoms, tag):
    """ROHF → AVAS(единое CAS(9e,10o)) → робастный CASSCF → SC-NEVPT2 + NOON."""
    mol = build_mol(atoms)
    mol.max_memory = 12000
    mf = rohf(mol)
    mo, wsel = forced_avas(mf)
    ne, no = CAS_TARGET
    # шаг 1: CASCI + натурбитали — устойчивый старт для CASSCF (рецепт Этапа 14)
    mc0 = mcscf.CASCI(mf, no, ne)
    mc0.fcisolver = fci.direct_spin1.FCI(mol)
    fci.addons.fix_spin_(mc0.fcisolver, shift=0.2, ss=SPIN / 2 * (SPIN / 2 + 1))
    mc0.kernel(mo)
    mo_nat = mc0.cas_natorb()[0]
    # шаг 2: CASSCF с натурбитального старта
    mc = mcscf.CASSCF(mf, no, ne)
    fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=SPIN / 2 * (SPIN / 2 + 1))
    mc.max_cycle_macro = 150
    mc.conv_tol = 1e-7
    mc.kernel(mo_nat)
    dm1 = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    noon = sorted(np.linalg.eigvalsh(dm1), reverse=True)
    e_pt = NEVPT(mc).kernel()          # SC-NEVPT2 корреляционная поправка
    return {
        "tag": tag, "avas_sel_weights": wsel, "cas": list(CAS_TARGET),
        "e_rohf_h": float(mf.e_tot), "rohf_converged": bool(mf.converged),
        "e_casci_h": float(mc0.e_tot),
        "e_casscf_h": float(mc.e_tot), "casscf_converged": bool(mc.converged),
        "e_nevpt2_h": float(mc.e_tot + e_pt), "nevpt2_corr_h": float(e_pt),
        "noon": [round(float(n), 3) for n in noon],
        "n_unpaired_beyond_spin": round(float(
            sum(min(n, 2 - n) for n in noon) - SPIN), 3),
    }


def stage_energy(sub):
    with open(os.path.join(DIR, f"{PREF}_{sub}_geom.json")) as f:
        g = json.load(f)
    rel, iR, iTS, interior = analyze_scan(g)
    out = {"substrate": sub, "cas": list(CAS_TARGET), "avas_labels": AVAS_LABELS,
           "r_index": iR, "ts_index": iTS, "ts_interior": interior,
           "rel_kcal": [round(x, 2) for x in rel],
           "dft_barrier_kcal": round(rel[iTS] - rel[iR], 2)}
    path = os.path.join(DIR, f"{PREF}_{sub}_energy.json")
    for name, idx in (("R", iR), ("TS", iTS)):
        atoms = [(s, tuple(c)) for s, c in g["points"][idx]["atoms"]]
        t0 = time.time()
        res = cas_nevpt2(atoms, f"{sub}_{name}")
        res["wall_s"] = round(time.time() - t0, 1)
        out[name] = res
        with open(path, "w") as f:
            json.dump(out, f, indent=1)
        print(f"[{sub}] {name}: CASSCF {res['e_casscf_h']:.6f} "
              f"(conv={res['casscf_converged']}) NEVPT2 {res['e_nevpt2_h']:.6f} "
              f"NOON={res['noon']} ({res['wall_s']}s)", flush=True)
    for lvl in ("casscf", "nevpt2"):
        out[f"{lvl}_barrier_kcal"] = round(
            (out["TS"][f"e_{lvl}_h"] - out["R"][f"e_{lvl}_h"]) * HARTREE_KCAL, 2)
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print(f"[{sub}] barriers kcal/mol: DFT {out['dft_barrier_kcal']} | "
          f"CASSCF {out['casscf_barrier_kcal']} | NEVPT2 {out['nevpt2_barrier_kcal']}",
          flush=True)


def stage_profile(sub):
    """CASSCF(9e,10o)+NEVPT2 во ВСЕХ точках скана: коррелированный профиль
    сам выбирает свои R и TS — без наследования индексов у DFT-поверхности
    (реперы честные и симметричные для обоих субстратов)."""
    gpath = os.path.join(DIR, f"{PREF}_{sub}_geom.json")
    with open(gpath) as f:
        g = json.load(f)
    path = os.path.join(DIR, f"{PREF}_{sub}_profile.json")
    out = {"substrate": sub, "cas": list(CAS_TARGET), "points": []}
    if os.path.exists(path):
        with open(path) as f:
            out = json.load(f)
    for i, p in enumerate(g["points"]):
        if i < len(out["points"]):
            continue
        atoms = [(s, tuple(c)) for s, c in p["atoms"]]
        t0 = time.time()
        res = cas_nevpt2(atoms, f"{sub}_p{i}")
        res["wall_s"] = round(time.time() - t0, 1)
        out["points"].append(res)
        with open(path, "w") as f:
            json.dump(out, f, indent=1)
        print(f"[{sub}] profile {i}: CASSCF {res['e_casscf_h']:.6f} "
              f"(conv={res['casscf_converged']}) NEVPT2 {res['e_nevpt2_h']:.6f} "
              f"({res['wall_s']}s)", flush=True)
    for lvl in ("casscf", "nevpt2"):
        es = [p[f"e_{lvl}_h"] for p in out["points"]]
        rel = [(e - es[0]) * HARTREE_KCAL for e in es]
        # форвард-барьер: TS = внутренний максимум, R = минимум слева от TS
        # (глобальный минимум может лежать на продуктовой стороне — экзотермика)
        imax = 1 + int(np.argmax(rel[1:-1]))
        imin = int(np.argmin(rel[:imax + 1]))
        out[f"{lvl}_rel_kcal"] = [round(x, 2) for x in rel]
        out[f"{lvl}_r_index"], out[f"{lvl}_ts_index"] = imin, imax
        out[f"{lvl}_ts_interior"] = bool(rel[imax] >= max(rel[0], rel[-1]))
        out[f"{lvl}_barrier_kcal"] = round(rel[imax] - rel[imin], 2)
    out["all_converged"] = all(p["casscf_converged"] for p in out["points"])
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print(f"[{sub}] correlated profile: CASSCF {out['casscf_barrier_kcal']} | "
          f"NEVPT2 {out['nevpt2_barrier_kcal']} kcal/mol "
          f"(R={out['nevpt2_r_index']} TS={out['nevpt2_ts_index']} "
          f"interior={out['nevpt2_ts_interior']} conv={out['all_converged']})",
          flush=True)


def _casscf_own(mf):
    """Собственное решение точки: fixed-AVAS → CASCI-натурбитали → CASSCF."""
    ne, no = CAS_TARGET
    mo, _ = forced_avas(mf)
    mc0 = mcscf.CASCI(mf, no, ne)
    mc0.fcisolver = fci.direct_spin1.FCI(mf.mol)
    fci.addons.fix_spin_(mc0.fcisolver, shift=0.2, ss=SPIN / 2 * (SPIN / 2 + 1))
    mc0.kernel(mo)
    mo_nat = mc0.cas_natorb()[0]
    mc = _casscf_from(mf, mo_nat)
    return mc


def _casscf_from(mf, mo_guess):
    ne, no = CAS_TARGET
    mc = mcscf.CASSCF(mf, no, ne)
    fci.addons.fix_spin_(mc.fcisolver, shift=0.2, ss=SPIN / 2 * (SPIN / 2 + 1))
    mc.max_cycle_macro = 150
    mc.conv_tol = 1e-7
    mc.kernel(mo_guess)
    return mc


def stage_refine(sub):
    """Орбитально-цепной ремонт коррелированного профиля: CASSCF каждой точки
    сеется проекцией сошедшихся орбиталей соседа (вперёд и назад), берём нижнее
    решение. Лечит выбросы типа «CASSCF ушёл в чужой локальный минимум»
    (c2h6 p6: +95 ккал) и зигзаги PT2 от смены референса."""
    gpath = os.path.join(DIR, f"{PREF}_{sub}_geom.json")
    ppath = os.path.join(DIR, f"{PREF}_{sub}_profile.json")
    with open(gpath) as f:
        geom = json.load(f)
    with open(ppath) as f:
        prof = json.load(f)
    n = len(prof["points"])

    def sweep(order, label):
        prev = None      # (mol, mo) последней сошедшейся точки цепочки
        for i in order:
            p = prof["points"][i]
            atoms = [(s, tuple(c)) for s, c in geom["points"][i]["atoms"]]
            mol = build_mol(atoms)
            mol.max_memory = 12000
            mf = rohf(mol)
            t0 = time.time()
            best = None
            if prev is not None:
                mol_p, mo_p = prev
                mo_g = mcscf.project_init_guess(mcscf.CASSCF(mf, *CAS_TARGET[::-1]),
                                                mo_p, prev_mol=mol_p)
                mc = _casscf_from(mf, mo_g)
                if mc.converged:
                    best = mc
            if best is None or best.e_tot > p["e_casscf_h"] + 5e-4:
                mc_own = _casscf_own(mf)
                if best is None or (mc_own.converged
                                    and mc_own.e_tot < best.e_tot):
                    best = mc_own
            if best.e_tot < p["e_casscf_h"] - 5e-4:
                e_pt = NEVPT(best).kernel()
                dm1 = best.fcisolver.make_rdm1(best.ci, best.ncas, best.nelecas)
                noon = sorted(np.linalg.eigvalsh(dm1), reverse=True)
                p["e_casscf_h"] = float(best.e_tot)
                p["casscf_converged"] = bool(best.converged)
                p["e_nevpt2_h"] = float(best.e_tot + e_pt)
                p["nevpt2_corr_h"] = float(e_pt)
                p["noon"] = [round(float(x), 3) for x in noon]
                p["refined"] = label
                print(f"[{sub}] refine-{label} {i}: lower CASSCF "
                      f"{p['e_casscf_h']:.6f} NEVPT2 {p['e_nevpt2_h']:.6f} "
                      f"({time.time()-t0:.0f}s)", flush=True)
            else:
                print(f"[{sub}] refine-{label} {i}: kept ({time.time()-t0:.0f}s)",
                      flush=True)
            prev = (mol, best.mo_coeff)
            with open(ppath, "w") as f:
                json.dump(prof, f, indent=1)

    sweep(range(n), "fwd")
    sweep(range(n - 1, -1, -1), "rev")
    for lvl in ("casscf", "nevpt2"):
        es = [p[f"e_{lvl}_h"] for p in prof["points"]]
        rel = [(e - es[0]) * HARTREE_KCAL for e in es]
        imax = 1 + int(np.argmax(rel[1:-1]))
        imin = int(np.argmin(rel[:imax + 1]))
        prof[f"{lvl}_rel_kcal"] = [round(x, 2) for x in rel]
        prof[f"{lvl}_r_index"], prof[f"{lvl}_ts_index"] = imin, imax
        prof[f"{lvl}_ts_interior"] = bool(rel[imax] >= max(rel[0], rel[-1]))
        prof[f"{lvl}_barrier_kcal"] = round(rel[imax] - rel[imin], 2)
    prof["all_converged"] = all(p["casscf_converged"] for p in prof["points"])
    prof["refine_done"] = True
    with open(ppath, "w") as f:
        json.dump(prof, f, indent=1)
    print(f"[{sub}] refined profile: CASSCF {prof['casscf_barrier_kcal']} | "
          f"NEVPT2 {prof['nevpt2_barrier_kcal']} kcal/mol "
          f"(R={prof['nevpt2_r_index']} TS={prof['nevpt2_ts_index']} "
          f"interior={prof['nevpt2_ts_interior']} conv={prof['all_converged']})",
          flush=True)


def stage_merge2():
    """Итог из refined-ПРОФИЛЕЙ (каждая поверхность — свои R/TS): барьеры и
    ΔΔE‡ по уровням DFT/CASSCF/NEVPT2 + лестница спинов, один самодостаточный
    JSON (его и забираем с AWS)."""
    res = {"model": MODEL, "spin_2S": SPIN, "cas": list(CAS_TARGET),
           "avas_labels": AVAS_LABELS, "substrates": {}}
    sp = os.path.join(DIR, f"{PREF}_spins.json")
    if os.path.exists(sp):
        with open(sp) as f:
            res["spin_ladder"] = json.load(f)
    for sub in ("ch4", "c2h6"):
        with open(os.path.join(DIR, f"{PREF}_{sub}_geom.json")) as f:
            g = json.load(f)
        with open(os.path.join(DIR, f"{PREF}_{sub}_profile.json")) as f:
            pr = json.load(f)
        res["substrates"][sub] = {
            "dft_rel_kcal": g["rel_kcal"], "dft_r_index": g.get("r_index", 0),
            "dft_ts_index": g["ts_index"], "dft_ts_interior": g["ts_interior"],
            "dft_barrier_kcal": g["dft_barrier_kcal"],
            "noon_ts": pr["points"][pr["nevpt2_ts_index"]]["noon"],
            **{k: pr[k] for k in pr if k.startswith(("casscf_", "nevpt2_", "all_"))}}
    res["ddE_kcal"] = {}
    for lvl in ("dft", "casscf", "nevpt2"):
        b = {s: res["substrates"][s][f"{lvl}_barrier_kcal"] for s in ("ch4", "c2h6")}
        res["ddE_kcal"][lvl] = round(b["c2h6"] - b["ch4"], 2)
    res["quantum_shift_kcal"] = round(
        res["ddE_kcal"]["nevpt2"] - res["ddE_kcal"]["dft"], 2)
    out = os.path.join(DIR, f"{PREF}_final.json")
    with open(out, "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "substrates"},
                     ensure_ascii=False, indent=1))
    print("saved", out, flush=True)


def stage_merge():
    res = {"model": "MnO (sextet) minimal radical-oxygen HAT site, stand-in for "
                    "Mn-Na2WO4/SiO2; NOT the periodic catalyst",
           "protocol": {"geometry": f"UKS-{XC.upper()}/{BASIS}, 2D-pin r(C-H)+r(O-H) "
                                    "relaxed scan (geomeTRIC), sextet",
                        "energy": "ROHF -> AVAS(Mn 3d + O 2p + transferring H 1s) "
                                  "-> CASSCF(9e,10o) [fix_spin, natural-orbital start] "
                                  "-> SC-NEVPT2; SAME curated CAS for all 4 structures",
                        "descriptor": "ddE = barrier(C2H6) - barrier(CH4); "
                                      "<0 over-oxidizes C2 (BEP wall); +4 needed for 70%"},
           "substrates": {}}
    ok = True
    for sub in ("ch4", "c2h6"):
        ep = os.path.join(DIR, f"{PREF}_{sub}_energy.json")
        pp = os.path.join(DIR, f"{PREF}_{sub}_profile.json")
        gp = os.path.join(DIR, f"{PREF}_{sub}_geom.json")
        if not (os.path.exists(pp) and os.path.exists(gp)):
            ok = False
            continue
        with open(pp) as f:
            prof = json.load(f)
        with open(gp) as f:
            geom = json.load(f)
        sd = {"dft_barrier_kcal": geom["dft_barrier_kcal"],
              "dft_rel_kcal": geom["rel_kcal"],
              "dft_r_index": geom["r_index"], "dft_ts_index": geom["ts_index"],
              "dft_ts_interior": geom["ts_interior"]}
        # коррелированные барьеры — из ПОЛНОГО профиля (свои R/TS на своей
        # поверхности), а не из энергостадии с DFT-индексами
        for lvl in ("casscf", "nevpt2"):
            for k in ("barrier_kcal", "rel_kcal", "r_index", "ts_index",
                      "ts_interior"):
                sd[f"{lvl}_{k}"] = prof[f"{lvl}_{k}"]
        sd["all_casscf_converged"] = prof["all_converged"]
        sd["noon_ts"] = prof["points"][prof["nevpt2_ts_index"]]["noon"]
        if os.path.exists(ep):
            with open(ep) as f:
                sd["energy_stage_dft_indices"] = json.load(f)
        res["substrates"][sub] = sd
    if ok:
        res["ddE_kcal"] = {}
        for lvl in ("dft", "casscf", "nevpt2"):
            b = {s: res["substrates"][s][f"{lvl}_barrier_kcal"] for s in ("ch4", "c2h6")}
            res["ddE_kcal"][lvl] = round(b["c2h6"] - b["ch4"], 2)
        res["quantum_shift_kcal"] = round(
            res["ddE_kcal"]["nevpt2"] - res["ddE_kcal"]["dft"], 2)
        res["all_casscf_converged"] = all(
            res["substrates"][s]["all_casscf_converged"] for s in ("ch4", "c2h6"))
        res["all_ts_interior"] = all(
            res["substrates"][s][f"{lvl}_ts_interior"]
            for s in ("ch4", "c2h6") for lvl in ("dft", "nevpt2"))
    with open(os.path.join(DIR, f"{PREF}_selectivity_results.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "substrates"}, indent=1))


if __name__ == "__main__":
    stage = sys.argv[1]
    if stage == "geom":
        stage_geom(sys.argv[2])
    elif stage == "polish":
        stage_polish(sys.argv[2])
    elif stage == "profile":
        stage_profile(sys.argv[2])
    elif stage == "refine":
        stage_refine(sys.argv[2])
    elif stage == "spins":
        stage_spins()
    elif stage == "merge2":
        stage_merge2()
    elif stage == "energy":
        stage_energy(sys.argv[2])
    elif stage == "merge":
        stage_merge()
    else:
        raise SystemExit("usage: geom|energy ch4|c2h6, or merge")
