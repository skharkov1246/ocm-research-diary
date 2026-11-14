#!/usr/bin/env python3
"""Самопроверка CalculiX (ccx) для навыка calculix-rotor.

Кольцевой тонкий стальной диск, расчёт в ccx и сравнение с аналитикой или вторым
способом:
  centrif  — сектор 30 град, C3D20R, *CYCLIC SYMMETRY MODEL, *DLOAD CENTRIF:
             sigma_theta против плоского напряжённого состояния;
  thermal  — тот же сектор, *TEMPERATURE по узлам T(r) = dT*((r-a)/(b-a))^2:
             sigma_theta против формулы тонкого диска;
  heatfile — *HEAT TRANSFER, STEADY STATE -> .frd -> *TEMPERATURE, FILE=..., BSTEP=1:
             T(r) против логарифмического профиля и sigma_theta против формулы;
  campbell — защемлённая ступица, *STATIC с CENTRIF (NLGEOM) и *STEP, PERTURBATION +
             *FREQUENCY + *SELECT CYCLIC SYMMETRY MODES при 0, 6000, 12000 об/мин;
             сектор против полной модели 360 град;
  creep    — один C3D8, *CREEP (Norton) и *VISCO: эквивалентная деформация ползучести
             против A*sigma^n*t; в *STATIC ползучести нет; при *PLASTIC с пределом sigma_y
             скорость A*(sigma - sigma_y)^n (модель перенапряжения);
  lcf      — цикл +-1 % деформации, *PLASTIC, HARDENING=KINEMATIC: размах напряжений против
             линейного кинематического упрочнения Прагера (с NLGEOM и без);
  gmsh     — сектор из Gmsh/OCC (C3D10, setPeriodic), *DLOAD CENTRIF: как centrif.

Запуск:  python3 ccx_selfcheck.py [--only centrif,campbell] [--work DIR] [--json out.json]
Каждый ccx — в один поток (OMP_NUM_THREADS=1), до 4 расчётов параллельно.
Код возврата 0 — все выбранные проверки в допуске.
"""
import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

E_MOD, NU, RHO, ALPHA, KCOND = 210e9, 0.3, 7800.0, 12e-6, 20.0
A_IN, B_OUT, THK = 0.05, 0.20, 0.004          # м
NSEC = 12
THETA_S = 2 * math.pi / NSEC
ENV1 = {k: "1" for k in ("OMP_NUM_THREADS", "CCX_NPROC_EQUATION_SOLVER",
                          "CCX_NPROC_RESULTS", "CCX_NPROC_STIFFNESS")}


# ---------------------------------------------------------------- сетка и вход
def ring_mesh(nr, nth, nz, theta, full):
    """Структурированная сетка C3D20R кольцевого сектора (full=True: 360 град)."""
    J = 2 * nth
    ids, nodes = {}, []
    for i in range(2 * nr + 1):
        for j in range(J if full else J + 1):
            for k in range(2 * nz + 1):
                if i % 2 + j % 2 + k % 2 > 1:
                    continue
                r = A_IN + (B_OUT - A_IN) * i / (2 * nr)
                th = theta * j / J
                z = -THK / 2 + THK * k / (2 * nz)
                ids[(i, j, k)] = len(nodes) + 1
                nodes.append(dict(id=len(nodes) + 1, x=r * math.cos(th), y=r * math.sin(th),
                                  z=z, r=r, th=th, i=i, j=j, k=k))

    def n(i, j, k):
        return ids[(i, j % J if full else j, k)]

    elems = []
    for ie in range(nr):
        for je in range(nth):
            for ke in range(nz):
                i0, j0, k0 = 2 * ie, 2 * je, 2 * ke
                q = [(i0, j0), (i0 + 2, j0), (i0 + 2, j0 + 2), (i0, j0 + 2)]
                m = [(i0 + 1, j0), (i0 + 2, j0 + 1), (i0 + 1, j0 + 2), (i0, j0 + 1)]
                conn = [n(a, b, k0) for a, b in q] + [n(a, b, k0 + 2) for a, b in q]
                conn += [n(a, b, k0) for a, b in m] + [n(a, b, k0 + 2) for a, b in m]
                conn += [n(a, b, k0 + 1) for a, b in q]
                elems.append(conn)
    return nodes, elems


def nset(name, ids):
    ids = sorted(set(ids))
    lines = [f"*NSET, NSET={name}"]
    for s in range(0, len(ids), 10):
        lines.append(", ".join(str(v) for v in ids[s:s + 10]))
    return "\n".join(lines)


def model_text(nodes, elems, etype="C3D20R"):
    out = ["*HEADING", "calculix-rotor selfcheck", "*NODE, NSET=NALL"]
    out += [f"{p['id']}, {p['x']:.12e}, {p['y']:.12e}, {p['z']:.12e}" for p in nodes]
    out.append(f"*ELEMENT, TYPE={etype}, ELSET=EALL")
    for e, c in enumerate(elems, 1):
        row = [e] + c
        out.append(", ".join(str(v) for v in row[:16]) + ("," if len(row) > 16 else ""))
        if len(row) > 16:
            out.append(", ".join(str(v) for v in row[16:]))
    out += ["*MATERIAL, NAME=STEEL", "*ELASTIC", f"{E_MOD:.6e}, {NU}", "*DENSITY", f"{RHO}",
            "*EXPANSION, ZERO=0.", f"{ALPHA:.6e}", "*CONDUCTIVITY", f"{KCOND}",
            "*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL"]
    return "\n".join(out)


CYCLIC = ("*SURFACE, NAME=SSLAVE, TYPE=NODE\nNSLAVE\n*SURFACE, NAME=SMASTER, TYPE=NODE\nNMASTER\n"
          "*TIE, NAME=CYC, CYCLIC SYMMETRY\nSSLAVE, SMASTER\n"
          f"*CYCLIC SYMMETRY MODEL, N={NSEC}, NGRAPH=1\n0., 0., 0., 0., 0., 1.")


def run_ccx(work, job, text, timeout=1200):
    d = os.path.join(work, job)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, job + ".inp"), "w") as f:
        f.write(text + "\n")
    t0 = time.time()
    p = subprocess.run(["ccx", "-i", job], cwd=d, capture_output=True, text=True,
                       timeout=timeout, env=dict(os.environ, **ENV1))
    log = p.stdout + p.stderr
    with open(os.path.join(d, "ccx.log"), "w") as f:
        f.write(log)
    if p.returncode != 0 or "*ERROR" in log:
        err = [ln for ln in log.splitlines() if "ERROR" in ln][:3]
        raise RuntimeError(f"{job}: ccx rc={p.returncode} {err}")
    return d, time.time() - t0


# ---------------------------------------------------------------- разбор выхода
def read_frd(path, block):
    """Узловой блок .frd (DISP, STRESS, NDTEMP) последнего шага -> {узел: [значения]}."""
    res, cur = {}, None
    with open(path) as f:
        for ln in f:
            if ln.startswith(" -4"):
                cur = {} if ln.split()[1] == block else None
                if cur is not None:
                    res = cur
            elif cur is not None and ln.startswith(" -1"):
                body = ln[13:].rstrip("\n")
                cur[int(ln[3:13])] = [float(body[s:s + 12]) for s in range(0, len(body) - 11, 12)]
            elif ln.startswith(" -3"):
                cur = None
    return res


def read_eigen(path):
    """Все таблицы EIGENVALUE OUTPUT из .dat -> [(nd, mode, f_Гц)].

    При циклической симметрии ccx пишет отдельную таблицу на каждый узловой диаметр,
    и каждая частота стоит в ней дважды (вещественная форма задачи удвоенного размера).
    """
    rows, on, got = [], False, 0
    with open(path) as f:
        for ln in f:
            if "E I G E N V A L U E   O U T P U T" in ln:
                on, got = True, 0
                continue
            if not on:
                continue
            p = ln.split()
            if len(p) == 6 and p[0].isdigit() and p[1].isdigit():
                rows.append((int(p[0]), int(p[1]), float(p[4])))
                got += 1
            elif len(p) == 5 and p[0].isdigit() and not p[1].isdigit():
                rows.append((None, int(p[0]), float(p[3])))
                got += 1
            elif p and not p[0].isdigit() and got:
                on = False
    return rows


def distinct(vals, rel):
    out = []
    for v in sorted(vals):
        if not out or v > out[-1] * (1 + rel):
            out.append(v)
    return out


def hoop(s, x, y):
    th = math.atan2(y, x)
    c, si = math.cos(th), math.sin(th)
    return s[0] * si * si + s[1] * c * c - 2 * s[3] * si * c


# ---------------------------------------------------------------- аналитика
def centrif_hoop(r, w):
    a, b = A_IN, B_OUT
    return (3 + NU) / 8 * RHO * w * w * (a * a + b * b + a * a * b * b / r ** 2 - (1 + 3 * NU) / (3 + NU) * r * r)


def thermal_hoop(r, T):
    """Тонкий кольцевой диск со свободными кромками, T(r) относительно T_ref."""
    a, b = A_IN, B_OUT
    rr = np.linspace(a, b, 20001)
    Tr = np.array([T(v) for v in rr])
    I = np.concatenate([[0.0], np.cumsum(0.5 * (Tr[1:] * rr[1:] + Tr[:-1] * rr[:-1]) * np.diff(rr))])
    Ir = np.interp(r, rr, I)
    return ALPHA * E_MOD * ((1 + a * a / r ** 2) / (b * b - a * a) * I[-1] + Ir / r ** 2 - T(r))


# ---------------------------------------------------------------- проверки
def sector_static(nr=32, nth=6, nz=2):
    nodes, elems = ring_mesh(nr, nth, nz, THETA_S, full=False)
    J = 2 * nth
    slave = [p["id"] for p in nodes if p["j"] == J]
    master = [p["id"] for p in nodes if p["j"] == 0]
    zfix = [p["id"] for p in nodes if p["j"] == 0 and p["i"] == 0 and p["k"] == nz]
    txt = "\n".join([model_text(nodes, elems), nset("NSLAVE", slave), nset("NMASTER", master),
                     nset("NZFIX", zfix), CYCLIC, "*BOUNDARY", "NMASTER, 2, 2", "NZFIX, 3, 3"])
    return nodes, txt


def hoop_at(work_dir, job, nodes, radii):
    S = read_frd(os.path.join(work_dir, job + ".frd"), "STRESS")
    J = max(p["j"] for p in nodes)
    out = {}
    for r in radii:
        p = min((q for q in nodes if q["j"] == J // 2 and abs(q["z"]) < 1e-12), key=lambda q: abs(q["r"] - r))
        out[r] = hoop(S[p["id"]], p["x"], p["y"])
    return out


def chk_centrif(work):
    w = 10000 * 2 * math.pi / 60
    nodes, txt = sector_static()
    txt += f"\n*STEP\n*STATIC\n*DLOAD\nEALL, CENTRIF, {w * w:.9e}, 0., 0., 0., 0., 0., 1.\n" \
           "*NODE FILE\nU\n*EL FILE\nS\n*END STEP"
    d, dt = run_ccx(work, "centrif", txt)
    h = hoop_at(d, "centrif", nodes, [A_IN, B_OUT])
    rows = [dict(r=r, fe_MPa=h[r] / 1e6, an_MPa=centrif_hoop(r, w) / 1e6) for r in h]
    for x in rows:
        x["dev"] = x["fe_MPa"] / x["an_MPa"] - 1
    return dict(rows=rows, tol=3e-3, ok=bool(all(abs(x["dev"]) < 3e-3 for x in rows)), sec=dt)


def chk_thermal(work):
    dT = 200.0
    T = lambda r: dT * ((r - A_IN) / (B_OUT - A_IN)) ** 2
    nodes, txt = sector_static()
    temps = "\n".join(f"{p['id']}, {T(p['r']):.6f}" for p in nodes)
    txt += "\n*INITIAL CONDITIONS, TYPE=TEMPERATURE\nNALL, 0.\n*STEP\n*STATIC\n*TEMPERATURE\n" + temps + \
           "\n*NODE FILE\nU, NT\n*EL FILE\nS\n*END STEP"
    d, dt = run_ccx(work, "thermal", txt)
    h = hoop_at(d, "thermal", nodes, [A_IN, B_OUT])
    rows = [dict(r=r, fe_MPa=h[r] / 1e6, an_MPa=thermal_hoop(r, T) / 1e6) for r in h]
    for x in rows:
        x["dev"] = x["fe_MPa"] / x["an_MPa"] - 1
    return dict(rows=rows, tol=3e-3, ok=bool(all(abs(x["dev"]) < 3e-3 for x in rows)), sec=dt)


def chk_heatfile(work):
    T1, T2 = 100.0, 500.0
    nodes, elems = ring_mesh(32, 6, 2, THETA_S, full=False)
    bore = [p["id"] for p in nodes if p["i"] == 0]
    rim = [p["id"] for p in nodes if p["i"] == 64]
    heat = "\n".join([model_text(nodes, elems), nset("NBORE", bore), nset("NRIM", rim),
                      "*INITIAL CONDITIONS, TYPE=TEMPERATURE", "NALL, 0.", "*STEP",
                      "*HEAT TRANSFER, STEADY STATE", "1., 1.", "*BOUNDARY",
                      f"NBORE, 11, 11, {T1}", f"NRIM, 11, 11, {T2}", "*NODE FILE", "NT", "*END STEP"])
    d, dt1 = run_ccx(work, "heat", heat)
    NT = read_frd(os.path.join(d, "heat.frd"), "NDTEMP")
    Tlog = lambda r: T1 + (T2 - T1) * math.log(r / A_IN) / math.log(B_OUT / A_IN)
    dev_T = max(abs(NT[p["id"]][0] - Tlog(p["r"])) for p in nodes)
    _, txt = sector_static()
    txt += "\n*INITIAL CONDITIONS, TYPE=TEMPERATURE\nNALL, 0.\n*STEP\n*STATIC\n" \
           "*TEMPERATURE, FILE=heat.frd, BSTEP=1\n*NODE FILE\nU, NT\n*EL FILE\nS\n*END STEP"
    d2 = os.path.join(work, "heatfile")
    os.makedirs(d2, exist_ok=True)
    shutil.copy(os.path.join(d, "heat.frd"), d2)
    d2, dt2 = run_ccx(work, "heatfile", txt)
    h = hoop_at(d2, "heatfile", nodes, [A_IN, B_OUT])
    rows = [dict(r=r, fe_MPa=h[r] / 1e6, an_MPa=thermal_hoop(r, Tlog) / 1e6) for r in h]
    for x in rows:
        x["dev"] = x["fe_MPa"] / x["an_MPa"] - 1
    ok = bool(dev_T < 0.05 and all(abs(x["dev"]) < 3e-3 for x in rows))
    return dict(T_max_abs_dev_K=dev_T, rows=rows, tol=3e-3, ok=ok, sec=dt1 + dt2)


def freq_deck(rpm, full, nr=12, nth=4, nz=2):
    if full:
        nodes, elems = ring_mesh(nr, nth * NSEC, nz, 2 * math.pi, full=True)
        clamp = [p["id"] for p in nodes if p["i"] == 0]
        parts = [model_text(nodes, elems), nset("NCLAMP", clamp)]
    else:
        nodes, elems = ring_mesh(nr, nth, nz, THETA_S, full=False)
        J = 2 * nth
        slave = [p["id"] for p in nodes if p["j"] == J]
        master = [p["id"] for p in nodes if p["j"] == 0]
        clamp = [p["id"] for p in nodes if p["i"] == 0 and p["j"] != J]   # без узлов slave
        parts = [model_text(nodes, elems), nset("NSLAVE", slave), nset("NMASTER", master),
                 nset("NCLAMP", clamp), CYCLIC]
    parts += ["*BOUNDARY", "NCLAMP, 1, 3"]
    w = rpm * 2 * math.pi / 60
    if rpm > 0:
        parts += ["*STEP, NLGEOM", "*STATIC", "*DLOAD",
                  f"EALL, CENTRIF, {w * w:.9e}, 0., 0., 0., 0., 0., 1.", "*END STEP", "*STEP, PERTURBATION"]
    else:
        parts += ["*STEP"]
    parts += ["*FREQUENCY", "40" if full else "12"]   # 12 строк = 6 разных частот на ND
    if not full:
        parts += [f"*SELECT CYCLIC SYMMETRY MODES, NMIN=0, NMAX={NSEC // 2}"]
    parts += ["*END STEP"]
    return "\n".join(parts)


def chk_campbell(work, rpms=(0, 6000, 12000)):
    jobs = [(f"freq_sec_{n}", freq_deck(n, False)) for n in rpms] + \
           [(f"freq_full_{n}", freq_deck(n, True)) for n in (rpms[0], rpms[-1])]
    with ThreadPoolExecutor(4) as ex:
        res = dict(zip([j for j, _ in jobs], ex.map(lambda jt: run_ccx(work, *jt), jobs)))
    table, cmp_rows, ok = {}, [], True
    for n in rpms:
        rows = read_eigen(os.path.join(res[f"freq_sec_{n}"][0], f"freq_sec_{n}.dat"))
        fam = {}
        for nd, _, f in rows:
            fam.setdefault(nd, [])
            if not any(abs(f - g) <= 1e-6 * max(f, 1e-9) for g in fam[nd]):
                fam[nd].append(f)
        table[n] = {nd: sorted(v) for nd, v in sorted(fam.items())}
        if f"freq_full_{n}" in res:
            full = sorted(f for _, _, f in read_eigen(os.path.join(res[f"freq_full_{n}"][0],
                                                                     f"freq_full_{n}.dat")))
            # Семейство гармонического индекса h содержит и формы с ND = N-h, N+h, ...,
            # а кратности в выводе ccx удвоены, поэтому сравниваются множества разных частот.
            cut = min(min(v[-1] for v in fam.values()), full[-1])   # выше — не всё досчитано
            uf = distinct([f for f in full if f < cut * 0.999], 1e-4)
            uc = distinct([f for v in fam.values() for f in v if f < cut * 0.999], 1e-4)
            dev = max([min(abs(x / y - 1) for y in uc) for x in uf] +
                      [min(abs(x / y - 1) for y in uf) for x in uc])
            cmp_rows.append(dict(rpm=n, below_Hz=round(cut, 1), n_full=len(uf), n_sector=len(uc),
                                 max_rel_dev=dev))
            ok = bool(ok and dev < 2e-3 and len(uf) == len(uc))
    cross = []
    for nd in range(2, NSEC // 2 + 1):   # пересечение низшей формы ND с линией EO = ND
        f = [table[n][nd][0] for n in rpms]
        g = [nd * n / 60 - fv for n, fv in zip(rpms, f)]
        for s in range(len(rpms) - 1):
            if g[s] < 0 <= g[s + 1]:
                t = -g[s] / (g[s + 1] - g[s])
                cross.append(dict(nd=nd, eo=nd, rpm=rpms[s] + t * (rpms[s + 1] - rpms[s])))
    stiff = all(table[rpms[-1]][nd][0] > table[0][nd][0] for nd in table[0])
    # Саутвелл f^2 = f0^2 + K (n/60)^2: K по 0 и средней скорости, прогноз на верхнюю
    sw = {}
    if len(rpms) >= 3:
        n1, n2 = rpms[len(rpms) // 2], rpms[-1]
        for nd in table[0]:
            f0, f1, f2 = table[0][nd][0], table[n1][nd][0], table[n2][nd][0]
            K = (f1 ** 2 - f0 ** 2) / (n1 / 60) ** 2
            sw[str(nd)] = dict(K=round(K, 3), pred_dev=math.sqrt(f0 ** 2 + K * (n2 / 60) ** 2) / f2 - 1)
    return dict(table_Hz={str(n): {str(k): [round(x, 3) for x in v[:3]] for k, v in t.items()}
                          for n, t in table.items()},
                sector_vs_full=cmp_rows, crossings=cross, spin_stiffening=stiff, southwell=sw,
                ok=bool(ok and stiff), sec=sum(v[1] for v in res.values()))


CUBE = ["*NODE, NSET=NALL", "1, 0., 0., 0.", "2, 1., 0., 0.", "3, 1., 1., 0.", "4, 0., 1., 0.",
        "5, 0., 0., 1.", "6, 1., 0., 1.", "7, 1., 1., 1.", "8, 0., 1., 1.",
        "*ELEMENT, TYPE=C3D8, ELSET=EALL", "1, 1, 2, 3, 4, 5, 6, 7, 8", "*NSET, NSET=NTOP", "5, 6, 7, 8"]
CUBE_BC = ["*BOUNDARY", "1, 1, 3", "2, 2, 3", "3, 3, 3", "4, 1, 1", "4, 3, 3", "5, 1, 2", "6, 2, 2", "8, 1, 1"]


def read_elprint(path):
    """Блоки *EL PRINT из .dat -> [(время, {'S': [...], 'E': [...], 'P': [...]})], средние по точкам."""
    out, cur, key = [], None, None
    with open(path) as f:
        for ln in f:
            for tag, k in (("stresses (elem", "S"), ("strains (elem", "E"), ("equivalent plastic strain", "P")):
                if tag in ln:
                    t = float(ln.split()[-1])
                    if not out or out[-1][0] != t:
                        out.append((t, {}))
                    cur, key = out[-1][1], k
                    cur[k] = []
                    break
            else:
                p = ln.split()
                if cur is not None and len(p) >= 3 and p[0].isdigit():
                    cur[key].append([float(x) for x in p[2:]])
    return [(t, {k: list(np.mean(v, axis=0)) for k, v in d.items() if v}) for t, d in out]


def chk_creep(work):
    """Norton в *VISCO; *STATIC ползучесть выключает; при *PLASTIC скорость ~ (sigma - sigma_y)^n."""
    # единицы Н, мм, т, МПа: в СИ (A = 1e-46 Па^-5) ccx 2.21 расходился на первом приращении
    A, n, sig, t, E_mpa = 1e-16, 5.0, 100.0, 1000.0, 210000.0
    variants = dict(visco=("*VISCO, CETOL=2.e-5", None, A * sig ** n * t),
                    static=("*STATIC", None, 0.0),
                    yield50=("*VISCO, CETOL=2.e-5", 50.0, A * (sig - 50.0) ** n * t),
                    yield150=("*VISCO, CETOL=2.e-5", 150.0, 0.0))
    rows, ok, sec = {}, True, 0.0
    for name, (proc, sy, expect) in variants.items():
        mat = ["*MATERIAL, NAME=M", "*ELASTIC", f"{E_mpa}, {NU}"]
        if sy is not None:
            mat += ["*PLASTIC", f"{sy}, 0.", f"{sy}, 10."]
        mat += ["*CREEP", f"{A:.3e}, {n}, 0.", "*SOLID SECTION, ELSET=EALL, MATERIAL=M"]
        txt = "\n".join(["*HEADING", "creep"] + CUBE + mat + CUBE_BC +
                        ["*STEP, NLGEOM, INC=10000", proc, f"1., {t}", "*DLOAD", f"1, P2, {-sig}",
                         "*EL PRINT, ELSET=EALL, FREQUENCY=10000", "S, PEEQ", "*END STEP"])
        d, dt = run_ccx(work, "creep_" + name, txt)
        sec += dt
        last = read_elprint(os.path.join(d, f"creep_{name}.dat"))[-1][1]
        peeq = last["P"][0]
        good = bool(abs(peeq - expect) <= 1e-2 * max(expect, 1e-12) + 1e-9)
        rows[name] = dict(peeq=peeq, expected=expect, szz_MPa=last["S"][2], ok=good)
        ok &= good
    return dict(rows=rows, ok=ok, sec=sec)


def chk_lcf(work):
    """Цикл +-1 % деформации, линейное кинематическое упрочнение: Прагер против ccx."""
    E_mpa, sy, H, eps = 210000.0, 300.0, 1000.0, 0.01
    s_an = (sy + H * eps) / (1 + H / E_mpa)          # петля симметрична: +-s_an
    rows, sec = {}, 0.0
    for name, nl in (("nlgeom", ", NLGEOM"), ("linear", "")):
        steps = []
        for k, u in enumerate((eps, -eps, eps)):
            steps += [f"*STEP, INC=1000{nl}", "*STATIC", "0.02, 1.", "*BOUNDARY", f"NTOP, 3, 3, {u}"]
            if k == 0:   # запрос вывода действует и в следующих шагах
                steps += ["*EL PRINT, ELSET=EALL", "S, E, PEEQ"]
            steps += ["*END STEP"]
        txt = "\n".join(["*HEADING", "lcf"] + CUBE +
                        ["*MATERIAL, NAME=M", "*ELASTIC", f"{E_mpa}, {NU}", "*PLASTIC, HARDENING=KINEMATIC",
                         f"{sy}, 0.", f"{sy + H * 0.1}, 0.1", "*SOLID SECTION, ELSET=EALL, MATERIAL=M"] +
                        CUBE_BC + steps)
        d, dt = run_ccx(work, "lcf_" + name, txt)
        sec += dt
        ends = [v for tt, v in read_elprint(os.path.join(d, f"lcf_{name}.dat")) if abs(tt - round(tt)) < 1e-9]
        s = [v["S"][2] for v in ends]
        rows[name] = dict(s33_MPa=s, s_an_MPa=s_an, max_rel_dev=max(abs(abs(x) / s_an - 1) for x in s))
    ok = bool(rows["nlgeom"]["max_rel_dev"] < 5e-3)
    return dict(rows=rows, ok=ok, note="ok — только по NLGEOM; linear выводится для сведения", sec=sec)


def chk_gmsh(work, h=0.004):
    import gmsh
    w = 10000 * 2 * math.pi / 60
    d = os.path.join(work, "gmsh")
    os.makedirs(d, exist_ok=True)
    gmsh.initialize(interruptible=False)   # в рабочем потоке обработчик сигналов ставить нельзя
    gmsh.option.setNumber("General.Terminal", 0)
    occ = gmsh.model.occ
    o = occ.addCylinder(0, 0, -THK / 2, 0, 0, THK, B_OUT, angle=THETA_S)
    i = occ.addCylinder(0, 0, -THK / 2, 0, 0, THK, A_IN, angle=THETA_S)
    vol, _ = occ.cut([(3, o)], [(3, i)])
    occ.synchronize()
    faces = [f for _, f in gmsh.model.getBoundary(vol, oriented=False)]
    cen = {f: occ.getCenterOfMass(2, f) for f in faces}
    fm = [f for f in faces if abs(cen[f][1]) < 1e-12 and abs(cen[f][2]) < 1e-9][0]
    fs = [f for f in faces if abs(math.atan2(cen[f][1], cen[f][0]) - THETA_S) < 1e-6][0]
    c, s = math.cos(THETA_S), math.sin(THETA_S)
    gmsh.model.mesh.setPeriodic(2, [fs], [fm], [c, -s, 0, 0, s, c, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    gmsh.model.addPhysicalGroup(3, [vol[0][1]], name="EVOL")   # физический — только объём
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.model.mesh.generate(3)
    inp = os.path.join(d, "mesh.inp")
    gmsh.write(inp)
    nm = gmsh.model.mesh.getNodes(2, fm, includeBoundary=True)[0]
    ns = gmsh.model.mesh.getNodes(2, fs, includeBoundary=True)[0]
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    gmsh.finalize()
    xyz = np.asarray(xyz).reshape(-1, 3)
    pos = dict(zip((int(t) for t in tags), xyz))
    nm, ns = sorted(int(t) for t in nm), sorted(int(t) for t in ns)
    zfix = min(nm, key=lambda t: np.hypot(pos[t][0] - A_IN, pos[t][2]))
    mesh = open(inp).read()
    mat = ["*MATERIAL, NAME=STEEL", "*ELASTIC", f"{E_MOD:.6e}, {NU}", "*DENSITY", f"{RHO}",
           "*SOLID SECTION, ELSET=EVOL, MATERIAL=STEEL",
           nset("NSLAVE", ns), nset("NMASTER", nm), nset("NZFIX", [zfix]), CYCLIC,
           "*BOUNDARY", "NMASTER, 2, 2", "NZFIX, 3, 3", "*STEP", "*STATIC", "*DLOAD",
           f"EVOL, CENTRIF, {w * w:.9e}, 0., 0., 0., 0., 0., 1.", "*NODE FILE", "U", "*EL FILE", "S", "*END STEP"]
    _, dt = run_ccx(work, "gmsh", mesh + "\n" + "\n".join(mat))
    S = read_frd(os.path.join(d, "gmsh.frd"), "STRESS")
    rows = []
    for r in (A_IN, B_OUT):
        cand = [t for t in S if abs(np.hypot(pos[t][0], pos[t][1]) - r) < 1e-7]
        vals = [hoop(S[t], pos[t][0], pos[t][1]) for t in cand]
        fe = float(np.median(vals))
        rows.append(dict(r=r, n_nodes=len(cand), fe_MPa=fe / 1e6, an_MPa=centrif_hoop(r, w) / 1e6,
                         dev=fe / centrif_hoop(r, w) - 1))
    return dict(nodes_master=len(nm), nodes_slave=len(ns), rows=rows, tol=1e-2,
                ok=bool(len(nm) == len(ns) and all(abs(x["dev"]) < 1e-2 for x in rows)), sec=dt)


CHECKS = dict(centrif=chk_centrif, thermal=chk_thermal, heatfile=chk_heatfile,
              campbell=chk_campbell, creep=chk_creep, lcf=chk_lcf, gmsh=chk_gmsh)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--only", default=",".join(CHECKS))
    ap.add_argument("--work", help="папка расчётов (по умолчанию временная, удаляется)")
    ap.add_argument("--json", help="записать результаты в JSON")
    a = ap.parse_args(argv)
    if not shutil.which("ccx"):
        sys.exit("ccx не найден: apt-get install -y calculix-ccx")
    work = a.work or tempfile.mkdtemp(prefix="ccx_selfcheck_")
    names = [s for s in a.only.split(",") if s]
    ver = subprocess.run(["ccx", "-v"], capture_output=True, text=True).stdout.strip()
    out = dict(ccx=ver, checks={})
    with ThreadPoolExecutor(2) as ex:   # campbell сам запускает до 4 ccx
        futs = {n: ex.submit(CHECKS[n], work) for n in names}
        for n in names:
            try:
                out["checks"][n] = futs[n].result()
            except Exception as e:  # noqa
                out["checks"][n] = dict(ok=False, error=str(e))
    ok = bool(all(v.get("ok") for v in out["checks"].values()))
    out["all_ok"] = ok
    print(json.dumps(out, ensure_ascii=False, indent=1, default=float))
    if a.json:
        with open(a.json, "w") as f:
            json.dump(out, f, ensure_ascii=False, indent=1, default=float)
    if not a.work:
        shutil.rmtree(work, ignore_errors=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
