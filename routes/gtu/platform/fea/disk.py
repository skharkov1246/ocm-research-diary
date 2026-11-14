#!/usr/bin/env python3
"""Прочность и ресурс диска / колеса ротора: осесимметричная модель в CalculiX.

Профиль диска задаётся таблицей: радиус r, осевая толщина t, (необязательно)
положение средней линии zc и осевая протяжённость лопаток hb. Нагрузка — вращение
(CENTRIF), поле температур T(r), растягивающая нагрузка на ободе от лопаток.
Лопатки, отлитые заодно с диском (колесо компрессора, радиальная турбина), не несут
окружной нагрузки: их масса «размазывается» в плотность столбца элементов под ними
(rho_eff = rho·(1 + phi·hb/t), phi — доля объёма, занятая лопатками).

Элементы CAX8 (8-узловые осесимметричные; ось симметрии — ось y модели CalculiX,
x — радиус). Материал с зависимостью E, alpha от температуры — materials.json.

Выход:
  * эквивалентные (Мизес) напряжения: в узлах (.frd, экстраполяция) и в точках
    интегрирования (.dat);
  * средние окружные напряжения по сечению — двумя способами (точки интегрирования
    .dat и узловое поле .frd) и сверка с равновесием полудиска
    ∫σθ dA = ω²∫ρ r² dA + p·r_o·t_o;
  * запас по разрушающей частоте вращения по средним окружным напряжениям:
    n_разр/n = sqrt(k·∫σв(T) dA / ∫σθ dA);
  * ресурс по ползучести (длительная прочность) через параметр Ларсона — Миллера
    для материала, у которого в materials.json есть блок lmp.

Запуск:
  python3 disk.py --examples            примеры (колёса схемы 130 и пакета Taurus 60)
  python3 disk.py case.json             свой расчёт; формат case — см. README.md
"""
import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MATERIALS = os.path.join(HERE, "materials.json")
G3 = (-math.sqrt(0.6), 0.0, math.sqrt(0.6))
W3 = (5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0)
# натуральные координаты узлов CAX8: углы 1-4 против часовой, затем середины сторон
XI_N = np.array([-1, 1, 1, -1, 0, 1, 0, -1], float)
ETA_N = np.array([-1, -1, 1, 1, -1, 0, 1, 0], float)


# ---------------------------------------------------------------- утилиты
def ccx_path():
    p = shutil.which("ccx")
    if not p:
        raise RuntimeError("ccx не найден: apt-get install -y calculix-ccx")
    return p


def ccx_version():
    try:
        out = subprocess.run([ccx_path(), "-v"], capture_output=True, text=True, timeout=30).stdout
        for line in out.splitlines():
            if "Version" in line:
                return line.strip()
    except Exception as e:  # noqa: BLE001
        return f"не определена ({e})"
    return "не определена"


def tab_interp(tab, x):
    """Линейная интерполяция по таблице [[x, y], ...]; за краем — крайнее значение.
    Возвращает (y, вне_таблицы)."""
    a = np.asarray(tab, float)
    if a.ndim != 2 or a.shape[0] == 0:
        raise ValueError("пустая таблица")
    xs, ys = a[:, 0], a[:, 1]
    out = bool(x < xs[0] - 1e-9 or x > xs[-1] + 1e-9)
    return float(np.interp(x, xs, ys)), out


def shape8(xi, eta):
    """Функции формы 8-узлового четырёхугольника и производные по xi, eta."""
    N = np.empty(8); dxi = np.empty(8); deta = np.empty(8)
    for k in range(8):
        a, b = XI_N[k], ETA_N[k]
        if a != 0 and b != 0:
            N[k] = 0.25 * (1 + xi * a) * (1 + eta * b) * (xi * a + eta * b - 1)
            dxi[k] = 0.25 * a * (1 + eta * b) * (2 * xi * a + eta * b)
            deta[k] = 0.25 * b * (1 + xi * a) * (xi * a + 2 * eta * b)
        elif a == 0:
            N[k] = 0.5 * (1 - xi * xi) * (1 + eta * b)
            dxi[k] = -xi * (1 + eta * b)
            deta[k] = 0.5 * b * (1 - xi * xi)
        else:
            N[k] = 0.5 * (1 + xi * a) * (1 - eta * eta)
            dxi[k] = 0.5 * a * (1 - eta * eta)
            deta[k] = -eta * (1 + xi * a)
    return N, dxi, deta


def von_mises(s11, s22, s33, s12, s13=0.0, s23=0.0):
    return np.sqrt(0.5 * ((s11 - s22) ** 2 + (s22 - s33) ** 2 + (s33 - s11) ** 2)
                   + 3.0 * (s12 ** 2 + s13 ** 2 + s23 ** 2))


# ---------------------------------------------------------------- материал
class LMP:
    """Длительная прочность по параметру Ларсона — Миллера P = T_K·(C + lg t_ч).
    lg σ = полином 2-й степени от P по точкам данных; вне диапазона данных —
    касательное продолжение (результат помечается как экстраполяция)."""

    def __init__(self, block):
        self.C = float(block.get("C", 20.0))
        d = np.array([[r[0], r[1], r[2]] for r in block["data"]], float)
        self.P = (d[:, 0] + 273.15) * (self.C + np.log10(d[:, 1]))
        self.y = np.log10(d[:, 2])
        self.P0 = float(self.P.mean())
        x = (self.P - self.P0) / 1000.0
        deg = 2 if len(x) >= 4 else 1
        self.coef = np.polyfit(x, self.y, deg)
        self.Pmin, self.Pmax = float(self.P.min()), float(self.P.max())
        dp = np.polyder(self.coef)
        for Pe in (self.Pmin, self.Pmax):
            if np.polyval(dp, (Pe - self.P0) / 1000.0) >= 0:
                raise ValueError("подгонка Ларсона — Миллера немонотонна на краю данных")
        self.check = block.get("check", [])
        self.spec_min = block.get("spec_min", [])

    def log_sigma_v(self, P):
        """lg σ(P), векторно; вне [Pmin, Pmax] — касательная к полиному на краю."""
        P = np.asarray(P, float)
        x = (P - self.P0) / 1000.0
        dp = np.polyder(self.coef)
        xl, xh = (self.Pmin - self.P0) / 1000.0, (self.Pmax - self.P0) / 1000.0
        y = np.polyval(self.coef, x)
        y = np.where(P < self.Pmin, np.polyval(self.coef, xl) + np.polyval(dp, xl) * (x - xl), y)
        y = np.where(P > self.Pmax, np.polyval(self.coef, xh) + np.polyval(dp, xh) * (x - xh), y)
        return y

    def log_sigma(self, P):
        return float(self.log_sigma_v(P))

    def sigma(self, T_C, t_h):
        P = (T_C + 273.15) * (self.C + math.log10(t_h))
        return 10.0 ** self.log_sigma(P)

    def P_of_sigma_v(self, s_MPa):
        """Обращение монотонной lg σ(P) бисекцией (векторно)."""
        ls = np.log10(np.asarray(s_MPa, float))
        lo = np.full(ls.shape, self.Pmin - 30000.0)
        hi = np.full(ls.shape, self.Pmax + 30000.0)
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            up = self.log_sigma_v(mid) > ls
            lo = np.where(up, mid, lo)
            hi = np.where(up, hi, mid)
        return 0.5 * (lo + hi)

    def P_of_sigma(self, s_MPa):
        return float(self.P_of_sigma_v(s_MPa))

    def life_h_v(self, s_MPa, T_C):
        """Время до разрушения, ч (векторно); σ <= 0 -> inf. Второй массив — экстраполяция."""
        s = np.asarray(s_MPa, float)
        T = np.asarray(T_C, float) * np.ones_like(s)
        pos = s > 0
        P = self.P_of_sigma_v(np.where(pos, s, 1.0))
        e = P / (T + 273.15) - self.C
        life = np.where(pos & (e < 300), 10.0 ** np.minimum(e, 300), np.inf)
        extrap = pos & ((P < self.Pmin) | (P > self.Pmax))
        return life, extrap

    def life_h(self, s_MPa, T_C):
        """Время до разрушения, ч; (inf, False) при σ <= 0. Второй элемент — экстраполяция."""
        life, ex = self.life_h_v([s_MPa], [T_C])
        return float(life[0]), bool(ex[0])


def load_material(spec):
    """spec — ключ из materials.json или словарь того же вида."""
    if isinstance(spec, str):
        with open(MATERIALS, encoding="utf-8") as f:
            lib = json.load(f)
        if spec not in lib:
            raise KeyError(f"материал {spec!r} не найден в materials.json")
        m = dict(lib[spec]); m["key"] = spec
    else:
        m = dict(spec); m.setdefault("key", m.get("name", "custom"))
    m["_lmp"] = LMP(m["lmp"]) if m.get("lmp") else None
    return m


# ---------------------------------------------------------------- профиль и сетка
def read_profile_csv(path):
    cols = {"r": [], "t": [], "zc": [], "hb": []}
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(line for line in f if not line.lstrip().startswith("#")):
            for k in cols:
                if k in row and row[k] not in (None, ""):
                    cols[k].append(float(row[k]))
    return {k: v for k, v in cols.items() if v}


class Profile:
    """interp = "linear" — толщина и средняя линия линейны между точками таблицы (изломы в
    точках таблицы — настоящие углы); "pchip" — гладкий монотонный сплайн через точки таблицы
    (без ложных концентраторов в изломах). Протяжённость лопаток hb — всегда линейно."""

    def __init__(self, p, interp="linear"):
        self.interp = interp
        self.r = np.asarray(p["r"], float)
        n = len(self.r)
        if n < 2 or np.any(np.diff(self.r) <= 0):
            raise ValueError("радиусы профиля должны строго возрастать, точек не меньше двух")
        self.t = np.asarray(p["t"], float)
        self.zc = np.asarray(p.get("zc", [0.0] * n), float)
        self.hb = np.asarray(p.get("hb", [0.0] * n), float)
        for name, a in (("t", self.t), ("zc", self.zc), ("hb", self.hb)):
            if len(a) != n:
                raise ValueError(f"колонка {name}: {len(a)} значений при {n} радиусах")
        if np.any(self.t <= 0):
            raise ValueError("толщина должна быть положительной")
        if self.r[0] < 0:
            raise ValueError("радиус не может быть отрицательным")
        if interp == "pchip":
            from scipy.interpolate import PchipInterpolator
            self._t, self._zc = PchipInterpolator(self.r, self.t), PchipInterpolator(self.r, self.zc)
        elif interp != "linear":
            raise ValueError("profile_interp: linear или pchip")

    def at(self, r):
        hb = float(np.interp(r, self.r, self.hb))
        if self.interp == "pchip":
            rr = min(max(r, self.r[0]), self.r[-1])
            return float(self._t(rr)), float(self._zc(rr)), hb
        return float(np.interp(r, self.r, self.t)), float(np.interp(r, self.r, self.zc)), hb


def radial_grid(prof, nr):
    """Угловые радиусы узлов: точки излома профиля — узлы сетки; nr элементов
    распределены по участкам пропорционально длине (не меньше одного на участок)."""
    L = prof.r[-1] - prof.r[0]
    seg = np.diff(prof.r)
    n = np.maximum(1, np.round(nr * seg / L).astype(int))
    rs = [prof.r[0]]
    for k in range(len(seg)):
        rs.extend(np.linspace(prof.r[k], prof.r[k + 1], n[k] + 1)[1:])
    return np.array(rs)


def build_mesh(prof, nr, nz):
    """Структурированная сетка CAX8. Узел (i, j) тонкой решётки (2·NR+1)×(2·nz+1);
    узлы с нечётными i и j (центры элементов) не создаются."""
    if nz % 2:
        nz += 1                       # средняя плоскость — угловой ряд узлов
    rc = radial_grid(prof, nr)
    NR = len(rc) - 1
    rf = np.empty(2 * NR + 1)
    rf[0::2] = rc
    rf[1::2] = 0.5 * (rc[:-1] + rc[1:])
    nj = 2 * nz + 1
    nid = lambda i, j: i * nj + j + 1  # noqa: E731
    nodes = {}
    for i, r in enumerate(rf):
        t, zc, _ = prof.at(r)
        for j in range(nj):
            if i % 2 and j % 2:
                continue
            nodes[nid(i, j)] = (float(r), zc + t * (j / (2 * nz) - 0.5), i, j)
    elems = []                         # (eid, [8 узлов], столбец)
    eid = 0
    for ie in range(NR):
        for je in range(nz):
            i0, j0 = 2 * ie, 2 * je
            conn = [nid(i0, j0), nid(i0 + 2, j0), nid(i0 + 2, j0 + 2), nid(i0, j0 + 2),
                    nid(i0 + 1, j0), nid(i0 + 2, j0 + 1), nid(i0 + 1, j0 + 2), nid(i0, j0 + 1)]
            eid += 1
            elems.append((eid, conn, ie))
    return dict(nodes=nodes, elems=elems, rc=rc, NR=NR, nz=nz, nid=nid)


def column_density_factor(prof, rc, phi):
    """1 + phi·∫hb·r dr / ∫t·r dr по каждому столбцу (составная формула Симпсона; при
    линейном профиле точна — подынтегральные функции квадратичные). Масса, размазанная в
    плотность, входит и в КЭ-модель, и в проверку равновесия — одинаково."""
    f = np.ones(len(rc) - 1)
    if phi <= 0:
        return f
    for k in range(len(rc) - 1):
        rs = np.linspace(rc[k], rc[k + 1], 5)          # составная формула Симпсона, 4 отрезка
        w = np.array([1, 4, 2, 4, 1], float)
        num = w @ np.array([prof.at(r)[2] * r for r in rs])
        den = w @ np.array([prof.at(r)[0] * r for r in rs])
        f[k] = 1.0 + phi * num / den
    return f


def temperature_field(case, r):
    T = case.get("T_C", 20.0)
    if isinstance(T, (int, float)):
        return float(T)
    return float(np.interp(r, T["r"], T["T"]))


# ---------------------------------------------------------------- .inp
def write_inp(path, case, mat, prof, mesh, dens_fac, omega, p_rim):
    nodes, elems = mesh["nodes"], mesh["elems"]
    NR, nz, nid = mesh["NR"], mesh["nz"], mesh["nid"]
    Tref = float(case.get("T_ref_C", 20.0))
    L = [f"*HEADING\ndisk.py {case['name']}\n*NODE"]
    for k in sorted(nodes):
        r, z = nodes[k][:2]
        L.append(f"{k}, {r:.9e}, {z:.9e}")
    cols = sorted(set(e[2] for e in elems))
    uniq = len(set(np.round(dens_fac, 12))) > 1
    for c in cols:
        L.append(f"*ELEMENT, TYPE=CAX8, ELSET=C{c + 1}")
        for eid, conn, col in elems:
            if col == c:
                L.append(f"{eid}, " + ", ".join(str(n) for n in conn))
    L.append("*ELSET, ELSET=EALL")
    L.extend(f"C{c + 1}," for c in cols)
    L.append("*NSET, NSET=NALL")
    ids = sorted(nodes)
    for k in range(0, len(ids), 16):
        L.append(", ".join(str(n) for n in ids[k:k + 16]) + ",")
    axis = [k for k, v in nodes.items() if v[2] == 0] if prof.r[0] <= 1e-12 else []
    if axis:
        L.append("*NSET, NSET=NAXIS")
        for k in range(0, len(axis), 16):
            L.append(", ".join(str(n) for n in axis[k:k + 16]) + ",")
    L.append(f"*NSET, NSET=NFIX\n{nid(0, nz)},")
    rho0 = float(mat["density"])
    E_tab = [(e * 1e9, T) for T, e in mat["E_GPa"]]
    a_tab = [(a * 1e-6, T) for T, a in mat["alpha_mean_1e6"]]
    groups = [(c, dens_fac[c]) for c in cols] if uniq else [(None, 1.0)]
    for c, fac in groups:
        name = f"M{c + 1}" if c is not None else "M1"
        L.append(f"*MATERIAL, NAME={name}\n*ELASTIC")
        L.extend(f"{E:.6e}, {mat['nu']}, {T}" for E, T in E_tab)
        L.append(f"*EXPANSION, ZERO={Tref}")
        L.extend(f"{a:.6e}, {T}" for a, T in a_tab)
        L.append(f"*DENSITY\n{rho0 * fac:.6f}")
        L.append(f"*SOLID SECTION, ELSET={'C%d' % (c + 1) if c is not None else 'EALL'}, MATERIAL={name}")
    L.append("*BOUNDARY")
    if axis:
        L.append("NAXIS, 1, 1, 0.")
    L.append("NFIX, 2, 2, 0.")
    L.append(f"*INITIAL CONDITIONS, TYPE=TEMPERATURE\nNALL, {Tref}")
    L.append("*STEP\n*STATIC\n*DLOAD")
    L.append(f"EALL, CENTRIF, {omega ** 2:.9e}, 0., 0., 0., 0., 1., 0.")
    if p_rim:
        for eid, conn, col in elems:
            if col == NR - 1:
                L.append(f"{eid}, P2, {-p_rim:.6e}")   # растяжение — отрицательное давление
    L.append("*TEMPERATURE")
    for k in ids:
        L.append(f"{k}, {temperature_field(case, nodes[k][0]):.4f}")
    L.append("*EL PRINT, ELSET=EALL\nS")
    L.append("*NODE FILE\nU, NT\n*EL FILE\nS\n*END STEP")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")


# ---------------------------------------------------------------- разбор
def parse_frd(path):
    """Блоки узловых результатов .frd: {'DISP': {n: [..]}, 'STRESS': ..., 'NDTEMP': ...}."""
    out, cur = {}, None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(" -4"):
                cur = line.split()[1]
                out[cur] = {}
            elif line.startswith(" -3"):
                cur = None
            elif cur and line.startswith(" -1"):
                n = int(line[3:13])
                vals, pos = [], 13
                while pos + 12 <= len(line.rstrip("\n")) + 0:
                    vals.append(float(line[pos:pos + 12]))
                    pos += 12
                out[cur][n] = vals
    return out


def parse_dat_stress(path):
    """Напряжения в точках интегрирования: {(elem, ip): (sxx, syy, szz, sxy, sxz, syz)}."""
    out, on = {}, False
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if "stresses (elem, integ.pnt." in line:
                on = True
                continue
            if on:
                s = line.split()
                if not s:
                    if out:
                        on = False
                    continue
                if len(s) != 8:
                    on = False
                    continue
                out[(int(s[0]), int(s[1]))] = tuple(float(v) for v in s[2:])
    return out


# ---------------------------------------------------------------- расчёт
def run_case(case, workdir=None, keep=True):
    t0 = time.time()
    mat = load_material(case["material"])
    pdat = case["profile"] if "profile" in case else read_profile_csv(case["profile_csv"])
    prof = Profile(pdat, case.get("profile_interp", "linear"))
    m = case.get("mesh", {})
    mesh = build_mesh(prof, int(m.get("nr", 60)), int(m.get("nz", 8)))
    nodes, elems = mesh["nodes"], mesh["elems"]
    NR, nz, nid = mesh["NR"], mesh["nz"], mesh["nid"]
    omega = 2 * math.pi * float(case["rpm"]) / 60.0
    r_o = prof.r[-1]
    t_o = prof.at(r_o)[0]
    rim = case.get("rim") or {}
    if "p_MPa" in rim:
        p_rim = float(rim["p_MPa"]) * 1e6
    elif "mass_kg" in rim:
        p_rim = float(rim["mass_kg"]) * omega ** 2 * float(rim["r_cg"]) / (2 * math.pi * r_o * t_o)
    else:
        p_rim = 0.0
    dens_fac = column_density_factor(prof, mesh["rc"], float(case.get("blade_fill", 0.0)))

    tag = f"{case['name']}_{NR}x{nz}"
    workdir = workdir or os.path.join(HERE, "runs", tag)
    os.makedirs(workdir, exist_ok=True)
    job = "disk"
    write_inp(os.path.join(workdir, job + ".inp"), case, mat, prof, mesh, dens_fac, omega, p_rim)
    # Решатель и разбор результатов — в один поток: однажды многопоточный ccx 2.21 (SPOOLES,
    # 4 потока) молча выдал решение с невязкой равновесия 1,6 %; в 90+ повторных прогонах того
    # же входа это не повторилось. Модели малы (секунды): параллелим по расчётам, не внутри ccx.
    nth = str(case.get("threads", 1))
    env = dict(os.environ, OMP_NUM_THREADS=nth, CCX_NPROC_EQUATION_SOLVER=nth,
               CCX_NPROC_RESULTS=nth, CCX_NPROC_STIFFNESS=nth, NUMBER_OF_PROCESSORS=nth)
    frd_p, dat_p = os.path.join(workdir, job + ".frd"), os.path.join(workdir, job + ".dat")
    attempts = []
    for attempt in range(2):
        for ext in (".frd", ".dat"):
            if os.path.exists(os.path.join(workdir, job + ext)):
                os.remove(os.path.join(workdir, job + ext))
        pr = subprocess.run([ccx_path(), "-i", job], cwd=workdir, capture_output=True, text=True,
                            timeout=int(case.get("timeout_s", 1200)), env=env)
        log = pr.stdout + pr.stderr
        with open(os.path.join(workdir, "ccx.log"), "w", encoding="utf-8") as f:
            f.write(log)
        if pr.returncode != 0 or "*ERROR" in log or not os.path.exists(frd_p):
            raise RuntimeError(f"ccx завершился с ошибкой (код {pr.returncode}); журнал: {workdir}/ccx.log")
        res = postprocess(case, mat, prof, mesh, parse_frd(frd_p), parse_dat_stress(dat_p),
                          omega, p_rim, dens_fac, r_o, t_o)
        attempts.append(res["quality"]["gate_equilibrium"])
        if res["quality"]["passed"]:
            break
        env.update(OMP_NUM_THREADS="1", CCX_NPROC_EQUATION_SOLVER="1", CCX_NPROC_RESULTS="1",
                   CCX_NPROC_STIFFNESS="1", NUMBER_OF_PROCESSORS="1")
    res["quality"]["attempts"] = attempts
    if not res["quality"]["passed"]:
        res["flags"].append("ПРОВЕРКА РАВНОВЕСИЯ НЕ ПРОЙДЕНА — результату не доверять")
    res["run"] = dict(dir=os.path.relpath(workdir, HERE) if keep else None,
                      seconds=round(time.time() - t0, 2), elements=len(elems), nodes=len(nodes),
                      ccx_threads=int(nth))
    if not keep:
        shutil.rmtree(workdir, ignore_errors=True)
    return res


def postprocess(case, mat, prof, mesh, frd, ips, omega, p_rim, dens_fac, r_o, t_o):
    nodes, elems, NR, nz, nid = mesh["nodes"], mesh["elems"], mesh["NR"], mesh["nz"], mesh["nid"]
    S, U = frd["STRESS"], frd["DISP"]
    TN = {k: v[0] for k, v in frd.get("NDTEMP", {}).items()}
    lmp = mat["_lmp"]
    flags = set()

    # --- узловые величины
    nd = {}
    for k, (r, z, i, j) in nodes.items():
        sr, sz, st, trz = (S[k][0], S[k][1], S[k][2], S[k][3])
        T = TN.get(k, temperature_field(case, r))
        c, d = 0.5 * (sr + sz), math.hypot(0.5 * (sr - sz), trz)
        nd[k] = dict(r=r, z=z, i=i, j=j, sr=sr, sz=sz, st=st, trz=trz,
                     vm=float(von_mises(sr, sz, st, trz)), s1=max(st, c + d), T=T,
                     ur=U[k][0], uz=U[k][1])

    # --- интегралы по сечению: точки интегрирования средней плоскости развёрнутого C3D20
    A = F_ip = F_nd = F_eq = IT = Iu = F_abs = 0.0
    vm_ip_max, vm_ip_loc = -1.0, None
    ip_rec = []                                    # (r, z, T, σ_Мизес) для ресурса по точкам
    for eid, conn, col in elems:
        xy = np.array([[nodes[n][0], nodes[n][1]] for n in conn])
        st_n = np.array([nd[n]["st"] for n in conn])
        T_n = np.array([nd[n]["T"] for n in conn])
        rho = float(mat["density"]) * dens_fac[col]
        for jq in range(3):
            for iq in range(3):
                N, dxi, deta = shape8(G3[iq], G3[jq])
                J = np.array([[dxi @ xy[:, 0], dxi @ xy[:, 1]], [deta @ xy[:, 0], deta @ xy[:, 1]]])
                dA = W3[iq] * W3[jq] * abs(np.linalg.det(J))
                r = float(N @ xy[:, 0]); z = float(N @ xy[:, 1]); T = float(N @ T_n)
                s = ips[(eid, 9 + 3 * jq + iq + 1)]
                vm = float(von_mises(*s))
                ip_rec.append((r, z, T, vm))
                if vm > vm_ip_max:
                    vm_ip_max, vm_ip_loc = vm, (r, z, T)
                su, out = tab_interp(mat["uts_MPa"], T)
                if out:
                    flags.add("σв: температура вне таблицы — взято крайнее значение")
                A += dA
                F_ip += s[2] * dA
                F_abs += abs(s[2]) * dA
                F_nd += float(N @ st_n) * dA
                F_eq += rho * omega ** 2 * r * r * dA
                IT += T * dA
                Iu += su * 1e6 * dA
    F_eq += p_rim * r_o * t_o
    # ворота качества: невязка равновесия полудиска и расхождение двух способов осреднения,
    # отнесённые к ∫|σθ|dA (определено и при чисто температурной нагрузке, где F_eq = 0)
    g_eq = abs(F_ip - F_eq) / F_abs if F_abs > 0 else 0.0
    g_nd = abs(F_nd - F_ip) / F_abs if F_abs > 0 else 0.0
    quality = dict(gate_equilibrium=g_eq, gate_nodal_vs_ip=g_nd, limits=[1e-3, 5e-3],
                   passed=bool(g_eq < 1e-3 and g_nd < 5e-3))
    s_avg_ip, s_avg_nd = F_ip / A, F_nd / A
    T_avg = IT / A
    k_b = float(case.get("burst_k", 1.0))
    burst = math.sqrt(k_b * Iu / F_ip) if F_ip > 0 else math.inf

    # --- максимумы и характерные точки
    kmax = max(nd, key=lambda k: nd[k]["vm"])
    kst = max(nd, key=lambda k: nd[k]["st"])
    kbore = nid(0, nz)
    krim = nid(2 * NR, nz)
    mid = [nd[nid(i, nz)] for i in range(2 * NR + 1)]

    # --- отношение к пределу текучести там, где температура в таблице
    ys = mat.get("ys_MPa")
    y_ratio, y_loc, n_out = None, None, 0
    if ys:
        T_hi = ys[-1][0]
        for k, v in nd.items():
            if v["T"] > T_hi + 1e-9:
                n_out += 1
                continue
            sy, _ = tab_interp(ys, v["T"])
            q = v["vm"] / (sy * 1e6)
            if y_ratio is None or q > y_ratio:
                y_ratio, y_loc = q, (v["r"], v["z"], v["T"])
        if n_out:
            flags.add(f"σ0,2: {n_out} узлов горячее таблицы предела текучести — в отношении не учтены")

    # --- ползучесть
    creep = None
    if lmp is not None:
        keys = list(nd)
        lives, exs = lmp.life_h_v([nd[k]["vm"] / 1e6 for k in keys], [nd[k]["T"] for k in keys])
        ib = int(np.argmin(lives))
        kb = keys[ib]
        lives1, exs1 = lmp.life_h_v([nd[k]["s1"] / 1e6 for k in keys], [nd[k]["T"] for k in keys])
        i1 = int(np.argmin(lives1))
        k1 = keys[i1]
        ipa = np.array(ip_rec)
        lives_ip, exs_ip = lmp.life_h_v(ipa[:, 3] / 1e6, ipa[:, 2])
        jb = int(np.argmin(lives_ip))
        life_sec, ex_sec = lmp.life_h(s_avg_ip / 1e6, T_avg)
        T_max = max(v["T"] for v in nd.values())
        life_hot, ex_hot = lmp.life_h(s_avg_ip / 1e6, T_max)
        creep = dict(
            method="Ларсон — Миллер, C=%.0f; lg σ — полином по P; данные: materials.json" % lmp.C,
            local_min_h=_fin(float(lives[ib])), local_extrapolated=bool(exs[ib]),
            local_at=dict(r_mm=nd[kb]["r"] * 1e3, z_mm=nd[kb]["z"] * 1e3,
                          T_C=nd[kb]["T"], vm_MPa=nd[kb]["vm"] / 1e6),
            local_s1_min_h=_fin(float(lives1[i1])), local_s1_extrapolated=bool(exs1[i1]),
            local_s1_at=dict(r_mm=nd[k1]["r"] * 1e3, z_mm=nd[k1]["z"] * 1e3, T_C=nd[k1]["T"],
                             s1_MPa=nd[k1]["s1"] / 1e6),
            local_ip_min_h=_fin(float(lives_ip[jb])), local_ip_extrapolated=bool(exs_ip[jb]),
            local_ip_at=dict(r_mm=ipa[jb, 0] * 1e3, z_mm=ipa[jb, 1] * 1e3, T_C=ipa[jb, 2],
                             vm_MPa=ipa[jb, 3] / 1e6),
            section_mean_h=_fin(life_sec), section_mean_extrapolated=ex_sec,
            section_hot_h=_fin(life_hot), section_hot_extrapolated=ex_hot,
            note=("local — по Мизесу в узле при его температуре (упругий пик, без перераспределения "
                  "ползучестью — консервативно); local_s1 — по наибольшему главному напряжению, только "
                  "растяжение (сжатие на горячем ободе — термическое, релаксирует ползучестью); "
                  "local_ip — по Мизесу в точках интегрирования; "
                  "section_mean — средние окружные при средней температуре сечения; section_hot — "
                  "средние окружные при наибольшей температуре"))

    rel = lambda a, b: (a - b) / b if b else float("nan")  # noqa: E731
    res = dict(
        name=case["name"], material=mat.get("name", mat["key"]), rpm=float(case["rpm"]),
        omega=omega, r_in_mm=prof.r[0] * 1e3, r_out_mm=r_o * 1e3, U_tip=omega * r_o,
        rim_p_MPa=p_rim / 1e6, blade_fill=float(case.get("blade_fill", 0.0)),
        mesh=dict(nr=NR, nz=nz),
        vm_max_node_MPa=nd[kmax]["vm"] / 1e6,
        vm_max_at=dict(r_mm=nd[kmax]["r"] * 1e3, z_mm=nd[kmax]["z"] * 1e3, T_C=nd[kmax]["T"]),
        vm_max_ip_MPa=vm_ip_max / 1e6,
        vm_max_ip_at=dict(r_mm=vm_ip_loc[0] * 1e3, z_mm=vm_ip_loc[1] * 1e3, T_C=vm_ip_loc[2]),
        hoop_max_MPa=nd[kst]["st"] / 1e6,
        hoop_max_at=dict(r_mm=nd[kst]["r"] * 1e3, z_mm=nd[kst]["z"] * 1e3),
        bore_mid=dict(r_mm=nd[kbore]["r"] * 1e3, hoop_MPa=nd[kbore]["st"] / 1e6,
                      radial_MPa=nd[kbore]["sr"] / 1e6, vm_MPa=nd[kbore]["vm"] / 1e6,
                      T_C=nd[kbore]["T"]),
        rim_mid=dict(hoop_MPa=nd[krim]["st"] / 1e6, radial_MPa=nd[krim]["sr"] / 1e6,
                     ur_mm=nd[krim]["ur"] * 1e3, T_C=nd[krim]["T"]),
        section=dict(area_mm2=A * 1e6, T_mean_C=T_avg,
                     hoop_mean_ip_MPa=s_avg_ip / 1e6, hoop_mean_nodal_MPa=s_avg_nd / 1e6,
                     hoop_mean_equilibrium_MPa=F_eq / A / 1e6,
                     err_ip_vs_equilibrium=rel(F_ip, F_eq), err_nodal_vs_ip=rel(F_nd, F_ip),
                     uts_mean_MPa=Iu / A / 1e6),
        burst=dict(k=k_b, margin=burst, rpm_burst=burst * float(case["rpm"]),
                   method="n_разр/n = sqrt(k·∫σв(T)dA / ∫σθ dA), σθ — точки интегрирования"),
        yield_ratio_max=y_ratio,
        yield_ratio_at=None if y_loc is None else dict(r_mm=y_loc[0] * 1e3, z_mm=y_loc[1] * 1e3, T_C=y_loc[2]),
        creep=creep,
        quality=quality,
        midplane=dict(r_mm=[v["r"] * 1e3 for v in mid], hoop_MPa=[v["st"] / 1e6 for v in mid],
                      radial_MPa=[v["sr"] / 1e6 for v in mid], vm_MPa=[v["vm"] / 1e6 for v in mid],
                      T_C=[v["T"] for v in mid]),
        flags=sorted(flags),
    )
    return res


def _fin(x):
    return None if not math.isfinite(x) else x


# ---------------------------------------------------------------- примеры
def _quad_T(r_in, r_out, T_in, T_out, n=11):
    rs = np.linspace(r_in, r_out, n)
    return dict(r=rs.tolist(), T=(T_in + (T_out - T_in) * ((rs - r_in) / (r_out - r_in)) ** 2).tolist())


def example_cases():
    """Условные профили (НЕ конструкция): размеры колёс — из gtu_explore / gtu_dropin."""
    R = 0.4348 / 2      # схема 130 (и её титановый вариант 139): D2 0,4348 м, 22 205 об/мин, b2 23,4 мм
    ctrl = dict(r=[0.020, 0.050, 0.060, 0.100, 0.130, 0.150, 0.190, R],
                t=[0.110, 0.100, 0.095, 0.065, 0.045, 0.035, 0.020, 0.012],
                hb=[0.0, 0.0, 0.080, 0.075, 0.060, 0.050, 0.035, 0.0234])
    comp = dict(
        name="comp130_ti", material="ti64", rpm=22205.0, profile=ctrl, profile_interp="pchip",
        blade_fill=0.06, T_C=_quad_T(0.020, R, 40.0, 190.0), burst_k=1.0,
        note=("Колесо компрессора схемы 130 (gtu_explore_results.json: id 130 — алюминиевое, id 139 — "
              "титановое, размеры одинаковы). Профиль ступицы, протяжённость и доля лопаток, поле "
              "температур — условные. Отверстие под стяжку 40 мм. T2 = 192 °C по comp_stage."))
    k = 0.6171 / 0.4348   # пакет Taurus 60, одно колесо, ПК 6,75: D2 0,617 м, 19 119 об/мин, b2 28,2 мм
    Rt = 0.6171 / 2
    ctrl_t60 = dict(r=[x * k for x in ctrl["r"][:-1]] + [Rt], t=[x * k for x in ctrl["t"]],
                    hb=[x * k for x in ctrl["hb"][:-1]] + [0.0282])
    comp_t60 = dict(
        name="comp_T60_1C_ti", material="ti64", rpm=19119.0, profile=ctrl_t60, profile_interp="pchip",
        blade_fill=0.06, T_C=_quad_T(0.020 * k, Rt, 50.0, 275.0), burst_k=1.0,
        note=("Колесо пакета Taurus 60, вариант «1Ц, без охлаждения» (gtu_dropin_results.json): "
              "окружная скорость 618 м/с — на принятом пределе для титана 620 м/с. Профиль — "
              "условный, подобие колеса схемы 130 с коэффициентом 1,419 (отверстие 57 мм). "
              "T2 = 277 °C по comp_stage."))
    ctrl_t60s = dict(ctrl_t60, r=[0.0] + ctrl_t60["r"][1:])
    comp_t60s = dict(
        comp_t60, name="comp_T60_1C_ti_solid", profile=ctrl_t60s,
        T_C=_quad_T(0.0, Rt, 50.0, 275.0),
        note="То же колесо Taurus 60 без центрального отверстия (вал-цапфы крепятся к ступице).")
    Rt3 = 0.397 / 2      # литая радиальная турбина ядра схемы 130: Dt 0,397 м, тот же вал
    ctrl_tb = dict(r=[0.0, 0.035, 0.040, 0.080, 0.120, 0.140, 0.160, Rt3],
                   t=[0.090, 0.087, 0.085, 0.060, 0.035, 0.026, 0.016, 0.008],
                   hb=[0.0, 0.0, 0.070, 0.065, 0.055, 0.050, 0.035, 0.024])
    turb = dict(
        name="turb130_in713c", material="in713c", rpm=22205.0, profile=ctrl_tb, profile_interp="pchip",
        blade_fill=0.07, T_C=_quad_T(0.0, Rt3, 620.0, 860.0), burst_k=1.0,
        note=("Литое колесо радиальной турбины ядра схемы 130 (TIT 950 °C, Dt 0,397 м). Колесо без "
              "отверстия (вал приваривается). Профиль, лопатки и поле температур — условные: обод "
              "860 °C ≈ относительная температура торможения T0 − U²/(2cp) при U = 462 м/с, центр 620 °C."))
    return [comp, comp_t60, comp_t60s, turb]


def _run_job(job):
    key, case = job
    return key, run_case(case)


def run_examples(out_json, workers=4):
    from concurrent.futures import ProcessPoolExecutor
    ver = ccx_version()
    cases = example_cases()
    jobs = []
    for c in cases:
        jobs.append(((c["name"], "base"), dict(c, mesh=dict(nr=60, nz=8))))
        jobs.append(((c["name"], "fine"), dict(c, mesh=dict(nr=120, nz=16))))
        if c["material"] == "in713c":
            m = load_material("in713c")
            m["alpha_mean_1e6"] = [[T, a * 1.15] for T, a in m["alpha_mean_1e6"]]
            m.pop("_lmp")
            m["key"] = "in713c_alpha_x1.15"
            jobs.append(((c["name"], "alpha115"),
                         dict(c, name=c["name"] + "_alpha115", material=m, mesh=dict(nr=120, nz=16))))
    with ProcessPoolExecutor(max_workers=workers) as ex:      # 4 ядра — 4 расчёта одновременно
        out = dict(ex.map(_run_job, jobs))
    results = []
    for c in cases:
        base, fine = out[(c["name"], "base")], out[(c["name"], "fine")]
        conv = dict(
            vm_max_node=fine["vm_max_node_MPa"] / base["vm_max_node_MPa"] - 1,
            hoop_mean=fine["section"]["hoop_mean_ip_MPa"] / base["section"]["hoop_mean_ip_MPa"] - 1,
            bore_hoop=fine["bore_mid"]["hoop_MPa"] / base["bore_mid"]["hoop_MPa"] - 1,
            burst_margin=fine["burst"]["margin"] / base["burst"]["margin"] - 1)
        if fine["creep"] and base["creep"]:
            for key in ("local_min_h", "local_s1_min_h", "local_ip_min_h", "section_mean_h"):
                if base["creep"][key] and fine["creep"][key]:
                    conv["creep_" + key] = fine["creep"][key] / base["creep"][key] - 1
        entry = dict(case={k: v for k, v in c.items() if k not in ("profile", "T_C")},
                     profile=c["profile"], T_C=c["T_C"], result=fine,
                     result_60x8=dict(vm_max_node_MPa=base["vm_max_node_MPa"],
                                      bore_hoop_MPa=base["bore_mid"]["hoop_MPa"],
                                      hoop_mean_MPa=base["section"]["hoop_mean_ip_MPa"],
                                      burst_margin=base["burst"]["margin"],
                                      quality_passed=base["quality"]["passed"]),
                     convergence_fine_vs_60x8=conv)
        if (c["name"], "alpha115") in out:
            sens = out[(c["name"], "alpha115")]
            entry["sensitivity_alpha_x1.15"] = dict(
                vm_max_node_MPa=sens["vm_max_node_MPa"], vm_max_at=sens["vm_max_at"],
                bore_hoop_MPa=sens["bore_mid"]["hoop_MPa"],
                creep_local_min_h=sens["creep"]["local_min_h"] if sens["creep"] else None,
                creep_local_ip_min_h=sens["creep"]["local_ip_min_h"] if sens["creep"] else None,
                creep_local_s1_min_h=sens["creep"]["local_s1_min_h"] if sens["creep"] else None,
                creep_local_s1_at=sens["creep"]["local_s1_at"] if sens["creep"] else None,
                creep_local_at=sens["creep"]["local_at"] if sens["creep"] else None,
                burst_margin=sens["burst"]["margin"], quality_passed=sens["quality"]["passed"])
        results.append(entry)
        r = fine
        cr = r["creep"]["local_min_h"] if r["creep"] else None
        print(f"{c['name']}: U={r['U_tip']:.0f} м/с, σМизес max={r['vm_max_node_MPa']:.0f} МПа "
              f"(r={r['vm_max_at']['r_mm']:.0f} мм), σθ ср={r['section']['hoop_mean_ip_MPa']:.0f} МПа, "
              f"запас по разрушению {r['burst']['margin']:.3f}, ползучесть "
              f"{'—' if cr is None else f'{cr:.0f} ч'}, равновесие {'да' if r['quality']['passed'] else 'НЕТ'}")
    lmp = load_material("in713c")["_lmp"]
    lmp_tab = [dict(T_C=T, life_h=t, sigma_MPa=round(lmp.sigma(T, t), 1))
               for T in (750, 800, 850, 900, 950) for t in (1e3, 1e4, 3e4, 1e5)]
    doc = dict(
        method=("Осесимметричная КЭ-модель диска (CAX8, CalculiX): вращение + поле температур + "
                "масса лопаток в плотности столбца; средние окружные напряжения, запас по разрушающей "
                "частоте вращения, длительная прочность по Ларсону — Миллеру"),
        software=f"CalculiX ccx ({ver}), python3, numpy {np.__version__}",
        date=time.strftime("%Y-%m-%d"),
        caveat="ПРОФИЛИ КОЛЁС И ПОЛЯ ТЕМПЕРАТУР УСЛОВНЫЕ; свойства материалов — см. materials.json (статусы)",
        lmp_in713c_typical=lmp_tab,
        cases=results)
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, indent=1, default=_json_default)
    try:
        plot_examples(results, os.path.join(HERE, "disk_examples.png"))
    except ImportError:
        print("matplotlib не установлен — рисунок не построен")
    return doc


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    raise TypeError(type(o))


def plot_examples(results, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # верх — профиль (равный масштаб осей), низ — напряжения и температура по средней плоскости;
    # вторая ось — температура: это нагрузка, порождающая термические напряжения (подписана цветом)
    fig, axes = plt.subplots(2, len(results), figsize=(5 * len(results), 8.5), squeeze=False)
    for col, e in enumerate(results):
        r = e["result"]; p = e["profile"]; mp = r["midplane"]
        ax = axes[0][col]
        pr = Profile(p, e["case"].get("profile_interp", "linear"))
        rg = np.linspace(pr.r[0], pr.r[-1], 400)
        tz = np.array([pr.at(x) for x in rg]) * 1e3
        rr = rg * 1e3
        ax.fill_between(rr, tz[:, 1] - tz[:, 0] / 2, tz[:, 1] + tz[:, 0] / 2, color="0.55", label="диск")
        ax.fill_between(rr, tz[:, 1] + tz[:, 0] / 2, tz[:, 1] + tz[:, 0] / 2 + tz[:, 2], facecolor="0.88",
                        hatch="///", edgecolor="0.35", lw=0.5, label="лопатки (только масса)")
        ax.plot(pr.r * 1e3, pr.zc * 1e3 + pr.t * 1e3 / 2, "ko", ms=3, label="точки таблицы")
        ax.set_title(e["case"]["name"] + " — профиль условный", fontsize=9)
        ax.set_xlabel("r, мм"); ax.set_ylabel("z, мм"); ax.set_aspect("equal"); ax.legend(fontsize=7)
        ax = axes[1][col]
        # цвет + тип линии: кривые различимы и без цвета
        ax.plot(mp["r_mm"], mp["hoop_MPa"], c="#1f4e79", ls="-", lw=1.6, label="σθ окружное")
        ax.plot(mp["r_mm"], mp["radial_MPa"], c="#2e7d32", ls="--", lw=1.4, label="σr радиальное")
        ax.plot(mp["r_mm"], mp["vm_MPa"], c="#6a1b9a", ls="-.", lw=1.4, label="σ по Мизесу")
        ax.axhline(r["section"]["hoop_mean_ip_MPa"], ls=(0, (1, 1)), c="k", lw=1.0,
                   label=f"σθ среднее {r['section']['hoop_mean_ip_MPa']:.0f} МПа")
        ax.set_xlabel("r, мм (средняя плоскость)"); ax.set_ylabel("напряжение, МПа"); ax.grid(alpha=0.3)
        ax2 = ax.twinx()
        ax2.plot(mp["r_mm"], mp["T_C"], c="#b71c1c", ls=":", lw=1.4, label="T, °C (правая ось)")
        ax2.set_ylabel("температура, °C", color="#b71c1c")
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2,
                  frameon=False)
        ax.set_title(f"n = {r['rpm']:.0f} об/мин, U = {r['U_tip']:.0f} м/с, запас по разрушению "
                     f"{r['burst']['margin']:.2f}", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", nargs="?", help="JSON-файл расчёта")
    ap.add_argument("--examples", action="store_true", help="примеры колёс схемы 130 и Taurus 60")
    ap.add_argument("--out", default=None, help="куда записать результат JSON")
    a = ap.parse_args()
    if a.examples:
        run_examples(a.out or os.path.join(HERE, "disk_results.json"))
    elif a.case:
        with open(a.case, encoding="utf-8") as f:
            case = json.load(f)
        res = run_case(case)
        out = a.out or os.path.splitext(a.case)[0] + "_result.json"
        with open(out, "w", encoding="utf-8") as f:
            json.dump(res, f, ensure_ascii=False, indent=1, default=_json_default)
        print(json.dumps({k: res[k] for k in ("vm_max_node_MPa", "section", "burst", "creep")},
                         ensure_ascii=False, indent=1, default=_json_default))
        print("записано:", out)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
