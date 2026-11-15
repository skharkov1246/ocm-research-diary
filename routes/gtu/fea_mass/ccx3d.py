"""Массовая 3D-прочность: входы CalculiX (ccx 2.21) для сектора с циклической симметрией, запуск и разбор выхода.

Два входа на деталь (сетка - общий файл mesh.inp через *INCLUDE):
  * static.inp - статика: три шага *STATIC с NLGEOM при 100, 115 и 120 % частоты вращения (*DLOAD CENTRIF) и поле
    температур в узлах; закрепление по правилу 3 навыка calculix-rotor: отверстие свободно, один узел отверстия вне
    периодических граней - окружное и осевое в цилиндрической системе (*TRANSFORM, TYPE=C); узлы slave не закрепляются;
  * freq.inp - частоты: шаг *FREQUENCY без нагрузки (0 об/мин, холодный), затем для каждой частоты вращения пара
    «*STATIC с NLGEOM (CENTRIF + температуры) -> *STEP, PERTURBATION + *FREQUENCY» - частоты с преднагружением во
    вращающейся системе, без Кориолиса; узловые диаметры 0..⌊N/2⌋ (*SELECT CYCLIC SYMMETRY MODES, NMAX задаётся явно -
    без него ccx 2.21 считает только ND 0). Опора - узлы отверстия (посадка на вал), без узлов slave.
Единицы - СИ (м, кг, с, Па, °C). ccx - в один поток (навык calculix-rotor, раздел 2); параллельно - разные детали.
"""
import math
import os
import subprocess
import time

import numpy as np

ENV1 = dict(OMP_NUM_THREADS="1", CCX_NPROC_EQUATION_SOLVER="1", CCX_NPROC_STIFFNESS="1", CCX_NPROC_RESULTS="1")


def ccx_version():
    try:
        p = subprocess.run(["ccx", "-v"], capture_output=True, text=True, timeout=30)
        for ln in (p.stdout + p.stderr).splitlines():
            if "Version" in ln:
                return ln.strip()
    except (OSError, subprocess.SubprocessError):
        return None
    return None


# ----------------------------------------------------------------------------- карточки
def material_cards(mat, name="MAT", elset="EALL"):
    """*MATERIAL: E(T), ν; средний КТР от 20 °C (*EXPANSION, ZERO=20); плотность. mat - materials3d.Material"""
    out = [f"*MATERIAL, NAME={name}", "*ELASTIC"]
    for T, E in mat.E_tab:
        out.append(f"{E * 1e9:.6e}, {mat.nu:.4f}, {T:.2f}")
    out.append(f"*EXPANSION, ZERO={mat.T_ref:.1f}")
    for T, a in mat.alpha_tab:
        out.append(f"{a * 1e-6:.6e}, {T:.2f}")
    out += ["*DENSITY", f"{mat.rho:.3f}", f"*SOLID SECTION, ELSET={elset}, MATERIAL={name}"]
    return out


def ring_cards(mat, rho, fill, e_theta, name="MRING", elset="ERING"):
    """мёртвый обод (хвостовики + выступы диска между пазами): ортотропный материал в цилиндрической системе (1 - радиус,
    2 - окружное, 3 - ось): E_r = E_z = fill·E, E_θ = e_theta·E (паз рвёт окружную связь), сдвиг rθ и θz - 0,1·G (оценка);
    ν23 ограничен условием положительной определённости (|ν23| < sqrt(E2/E3))"""
    out = ["*ORIENTATION, NAME=CYL, SYSTEM=CYLINDRICAL", "0., 0., 0., 0., 0., 1.", f"*MATERIAL, NAME={name}",
           "*ELASTIC, TYPE=ENGINEERING CONSTANTS"]
    nu = mat.nu
    for T, E in mat.E_tab:
        E = E * 1e9
        E1, E2, E3 = fill * E, e_theta * E, fill * E
        G = E / (2 * (1 + nu))
        n23 = 0.3 * math.sqrt(E2 / E3)
        out.append(f"{E1:.6e}, {E2:.6e}, {E3:.6e}, {nu:.4f}, {nu:.4f}, {n23:.6f}, {0.1 * G:.6e}, {fill * G:.6e},")
        out.append(f"{0.1 * G:.6e}, {T:.2f}")
    out.append(f"*EXPANSION, ZERO={mat.T_ref:.1f}")
    for T, a in mat.alpha_tab:
        out.append(f"{a * 1e-6:.6e}, {T:.2f}")
    out += ["*DENSITY", f"{rho * fill:.3f}", f"*SOLID SECTION, ELSET={elset}, MATERIAL={name}, ORIENTATION=CYL"]
    return out


def cyclic_cards(N):
    return ["*SURFACE, NAME=SSLAVE, TYPE=NODE", "NSLAVE", "*SURFACE, NAME=SMASTER, TYPE=NODE", "NMASTER",
            "*TIE, NAME=CYC, CYCLIC SYMMETRY", "SSLAVE, SMASTER", f"*CYCLIC SYMMETRY MODEL, N={int(N)}, NGRAPH=1",
            "0., 0., 0., 0., 0., 1."]


def temperature_cards(node_ids, T):
    return ["*TEMPERATURE"] + [f"{int(n)}, {t:.3f}" for n, t in zip(node_ids, T)]


def omega2(rpm):
    w = 2 * math.pi * rpm / 60.0
    return w * w


def static_deck(N, mat_cards, rpm, fractions, node_ids, T):
    """статика: шаги NLGEOM при rpm·f для f в fractions; вывод U и S в .frd по каждому шагу"""
    L = ["*HEADING", "fea_mass static", "*INCLUDE, INPUT=mesh.inp"] + mat_cards + cyclic_cards(N)
    L += ["*TRANSFORM, NSET=NFIX, TYPE=C", "0., 0., 0., 0., 0., 1.", "*BOUNDARY", "NFIX, 2, 3",
          "*TIME POINTS, NAME=TEND", "1.", "*INITIAL CONDITIONS, TYPE=TEMPERATURE", "NALL, 20."]
    for k, f in enumerate(fractions):
        L += ["*STEP, NLGEOM, INC=200", "*STATIC", "0.125, 1." if k == 0 else "0.5, 1.", "*DLOAD", f"EALL, CENTRIF, {omega2(rpm * f):.9e}, 0., 0., 0., 0., 0., 1."]
        L += temperature_cards(node_ids, T)
        L += ["*NODE FILE, TIME POINTS=TEND", "U", "*EL FILE, TIME POINTS=TEND", "S", "*END STEP"]   # только конец шага
    return "\n".join(L) + "\n"


def freq_deck(N, mat_cards, rpm, fractions, node_ids, T, nmodes, nd_max, cold_zero=False):
    """частоты с преднагружением: для каждой доли f (0 допускается) - *STATIC с NLGEOM (CENTRIF при f > 0 и поле
    температур режима), затем *STEP, PERTURBATION + *FREQUENCY. cold_zero - ещё шаг без нагрузки и нагрева в начале
    (частоты в покое на стенде, для сведения)"""
    sel = f"*SELECT CYCLIC SYMMETRY MODES, NMIN=0, NMAX={int(nd_max)}"
    L = ["*HEADING", "fea_mass frequency", "*INCLUDE, INPUT=mesh.inp"] + mat_cards + cyclic_cards(N)
    L += ["*BOUNDARY", "NBORE, 1, 3", "*INITIAL CONDITIONS, TYPE=TEMPERATURE", "NALL, 20."]
    if cold_zero:
        L += ["*STEP", "*FREQUENCY", str(int(nmodes)), sel, "*TEMPERATURE", "NALL, 20.", "*END STEP"]
    for k, f in enumerate(fractions):
        L += ["*STEP, NLGEOM, INC=200", "*STATIC", "0.125, 1." if k == 0 else "0.5, 1."]
        if f > 0:
            L += ["*DLOAD", f"EALL, CENTRIF, {omega2(rpm * f):.9e}, 0., 0., 0., 0., 0., 1."]
        L += temperature_cards(node_ids, T)
        L += ["*END STEP", "*STEP, PERTURBATION", "*FREQUENCY", str(int(nmodes)), sel, "*END STEP"]
    return "\n".join(L) + "\n"


# ----------------------------------------------------------------------------- запуск
def run_ccx(workdir, job, text, timeout=3600):
    os.makedirs(workdir, exist_ok=True)
    with open(os.path.join(workdir, job + ".inp"), "w") as f:
        f.write(text)
    t0, tc0 = time.time(), _children_cpu()
    try:
        p = subprocess.run(["ccx", "-i", job], cwd=workdir, capture_output=True, text=True, timeout=timeout,
                           env=dict(os.environ, **ENV1))
        log, rc = p.stdout + p.stderr, p.returncode
    except subprocess.TimeoutExpired as e:
        log, rc = (e.stdout or "") if isinstance(e.stdout, str) else "", "timeout"
    with open(os.path.join(workdir, job + ".log"), "w") as f:
        f.write(log[-200000:])
    err = [ln.strip() for ln in log.splitlines() if "*ERROR" in ln][:5]
    ok = rc == 0 and not err
    return dict(ok=ok, rc=rc, errors=err, seconds=round(time.time() - t0, 2), cpu_s=round(_children_cpu() - tc0, 2),
                warnings=sum(1 for ln in log.splitlines() if "*WARNING" in ln))


def _children_cpu():
    t = os.times()
    return t.children_user + t.children_system


# ----------------------------------------------------------------------------- разбор выхода
def read_frd_steps(path, blocks=("DISP", "STRESS")):
    """узловые блоки .frd по шагам: [ {блок: (ids (n,), values (n, k))}, ... ] в порядке номеров шагов; если в файле
    несколько приращений шага, берётся последнее (номер шага - из строки 1PSTEP)"""
    by_step, step, name, ids, vals = {}, None, None, [], []

    def flush():
        if name is not None and step is not None:
            by_step.setdefault(step, {})[name] = (np.asarray(ids, dtype=np.int64), np.asarray(vals, dtype=float))
    with open(path) as f:
        for ln in f:
            if ln.startswith("    1PSTEP"):
                p = ln.split()
                try:
                    step = int(p[3])
                except (IndexError, ValueError):
                    step = None
            elif ln.startswith(" -4"):
                flush()
                name = ln.split()[1]
                if name not in blocks:
                    name = None
                ids, vals = [], []
            elif name is not None and ln.startswith(" -1"):
                body = ln[13:].rstrip("\n")
                ids.append(int(ln[3:13]))
                vals.append([float(body[s:s + 12]) for s in range(0, len(body) - 11, 12)])
            elif ln.startswith(" -3"):
                flush()
                name = None
    flush()
    return [by_step[k] for k in sorted(by_step)]


def read_eigen_steps(path):
    """таблицы собственных частот .dat по шагам: [ {ND: [f, ...]} ] в порядке шагов с выводом частот. Строка таблицы ccx 2.21
    при циклике: «ND, номер формы, собственное число, ω, f [Гц], мнимая часть»; без циклики - без ND. Частоты в таблице ND
    идут парами (вещественная форма удвоенной задачи) - здесь уникальные (отн. допуск 1e-6)"""
    steps, cur, on, step_no, cur_no = [], None, False, None, None
    with open(path) as f:
        for ln in f:
            if "S T E P" in ln:
                try:
                    step_no = int(ln.split()[-1])
                except ValueError:
                    step_no = None
                on = False
                continue
            if "E I G E N V A L U E   O U T P U T" in ln:
                if cur is None or step_no != cur_no:          # таблицы всех ND одного шага - вместе
                    cur, cur_no = {}, step_no
                    steps.append(cur)
                on = True
                continue
            if not on:
                continue
            if "P A R T I C I P A T I O N" in ln or "E F F E C T I V E" in ln or "T O T A L" in ln:
                on = False
                continue
            p = ln.split()
            try:
                if len(p) == 6 and p[0].isdigit() and p[1].isdigit():
                    nd, fhz = int(p[0]), float(p[4])
                elif len(p) == 5 and p[0].isdigit():
                    nd, fhz = 0, float(p[3])
                else:
                    continue
            except ValueError:
                continue
            lst = cur.setdefault(nd, [])
            if not any(abs(fhz - g) <= 1e-6 * max(abs(fhz), 1e-9) for g in lst):
                lst.append(fhz)
    return [{nd: sorted(v) for nd, v in st.items()} for st in steps if st]


def von_mises(S):
    sxx, syy, szz, sxy, syz, szx = (S[:, k] for k in range(6))
    return np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3 * (sxy ** 2 + syz ** 2 + szx ** 2))


def hoop(S, xyz):
    th = np.arctan2(xyz[:, 1], xyz[:, 0])
    c, s = np.cos(th), np.sin(th)
    return S[:, 0] * s * s + S[:, 1] * c * c - 2 * S[:, 3] * s * c


def principal_max(S):
    out = np.empty(len(S))
    for i, (a, b, c, d, e, f) in enumerate(S):
        out[i] = np.linalg.eigvalsh(np.array([[a, d, f], [d, b, e], [f, e, c]]))[-1]
    return out
