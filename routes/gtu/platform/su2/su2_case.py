#!/usr/bin/env python3
"""Узел SU2 расчётной платформы: межлопаточный канал колеса по сетке и расчётной точке.

По сетке одного канала (.su2 или .su2.gz из routes/gtu/platform/geom/impeller.py и JSON сетки рядом) и
расчётной точке (results-JSON трека, например routes/gtu/frame3/frame3_point_results.json:stage1) собирает
конфигурацию сжимаемого RANS SU2 v8.6.0: SST, вращающаяся система координат, периодичность поворотом,
условия Giles — полные параметры на входе и статическое давление на выходе. Запускает SU2_CFD локально
(mpirun -np N, общий лимит времени, мягкий останов с записью рестарта) и разбирает историю: сходимость
невязок, расход на входе и выходе (дисбаланс), полная степень сжатия и изоэнтропический КПД по средним по
расходу p0 и T0 (MARKER_ANALYZE_AVERAGE= MASSFLUX), сравнение с расчётной точкой.

Команды:
  run       расчёт: разгон частоты и давления (две стадии первого порядка), затем точки по противодавлению
            (--points ДОЛЯ:ИТЕРАЦИЙ,...; каждая следующая — с решения предыдущей) — одна точка или линия частоты
  parse     заново разобрать папку прогона (case_summary.json)
  aws-plan  пакет задач ../aws/runner.py для характеристики (AWS не вызывается) и оценка машино-часов и
            долларов по замеру времени пробного прогона

Примеры:
  python3 su2_case.py run --geom-results ../../frame3/frame3_geom_results.json --stage stage1 --level coarse \\
          --point ../../frame3/frame3_point_results.json:stage1 --points 0.75:600 --np 4 --time-limit 2400 \\
          --shroud slip --vaneless slip --out <папка> --record su2_trial_results.json
  python3 su2_case.py run --mesh passage.su2 --mesh-json impeller.json --point point.json:stage1 --dry
  python3 su2_case.py parse <папка прогона>
  python3 su2_case.py aws-plan --trial-results su2_trial_results.json --out jobs_frame3.json

Процессы SU2 (поиск SU2_CFD, mpirun, лимит времени и мягкий останов, разбор журнала, деление маркеров
стенок за колесом) — функции run_case.py этой папки, проверенные дымовыми прогонами (README).
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import platform
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import run_case as RC  # noqa: E402

DATE = "2026-11-14"
MARKERS = ("inlet", "outlet", "hub", "shroud", "blade")

# Численная схема по умолчанию — та же, что прошла дымовую проверку run_case.py (core130 coarse, README):
# Roe, первый порядок на разгоне и MUSCL + Venkatakrishnan на точке, неявная схема, адаптивный CFL 1-20.
NUMERICS = dict(
    turb="SST", init_mach=0.4, aoa=90.0, tu=0.03, mut_ratio=10.0, gamma=1.4, R=287.058,
    mu_ref=1.716e-5, mu_t_ref=273.15, sutherland=110.4,
    relax_in=0.3, relax_out=0.3, nspan=10, scheme="ROE", limiter="VENKATAKRISHNAN", venkat=0.05,
    entropy_fix=0.05, cfl=5.0, cfl_adapt="YES", cfl_adapt_param="0.5, 1.1, 1.0, 20.0, 0.1", lin_prec="ILU",
    lin_iter=10, lin_error=1e-4, res_min=-10, wrt_freq=500, ramp_freq=5, output_files="RESTART, PARAVIEW",
)

# Критерии сходимости точки — навык su2-turbo, раздел 6 (приняты в навыке, не из документации SU2)
CONV = dict(orders=4.0, mass_tol=0.005, window=500, stab_tol=0.001)


# ----------------------------------------------------------------------------------------------
# вход: точка и сетка
# ----------------------------------------------------------------------------------------------
def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path: str) -> str:
    p = os.path.abspath(path)
    return os.path.relpath(p, REPO) if p.startswith(REPO + os.sep) else p


POINT_FIELDS = ("PR", "p01_kPa", "T01_C", "T02_C", "m_total_kg_s", "n_rpm", "D2_m", "c2_m_s",
                "alpha2_flow_from_radial_deg", "blades_full")


def load_point(spec: str) -> dict:
    """Расчётная точка ступени: 'путь.json:ключ' (поля — как stage1/stage2 в frame3_point_results.json).
    Возвращает величины в СИ и словарь «величина ← поле файла» для results-JSON."""
    path, sep, key = spec.rpartition(":")
    if not sep or not path or os.sep in key:
        path, key = spec, ""
    data = json.load(open(path, encoding="utf-8"))
    s = data[key] if key else data
    miss = [k for k in POINT_FIELDS if k not in s]
    if miss:
        raise SystemExit(f"в точке {spec} нет полей {miss}")
    m_total = float(s["m_total_kg_s"])
    double = s.get("inlet") == "двусторонний"
    n_sides = 2 if double else 1
    m_side = float(s.get("m_per_side_kg_s", m_total / n_sides))
    if abs(m_side * n_sides - m_total) > 1e-9 * m_total:
        raise SystemExit(f"точка {spec}: m_per_side × сторон ({m_side} × {n_sides}) ≠ m_total ({m_total})")
    T01, T02 = float(s["T01_C"]) + 273.15, float(s["T02_C"]) + 273.15
    a2 = math.radians(float(s["alpha2_flow_from_radial_deg"]))
    return dict(
        file=rel(path), key=key, sha256=sha256(path), label=s.get("label", key),
        PR=float(s["PR"]), p01=float(s["p01_kPa"]) * 1000.0, T01=T01, T02=T02,
        m_total=m_total, m_side=m_side, n_sides=n_sides, double_sided=double,
        rpm=float(s["n_rpm"]), omega=float(s["n_rpm"]) * 2 * math.pi / 60.0, r2=float(s["D2_m"]) / 2,
        U2=s.get("U2_m_s"), c_theta2=float(s["c2_m_s"]) * math.sin(a2), Z=int(s["blades_full"]),
        eta_is=s.get("eta_is"), eta_p=s.get("eta_p"),
        fields={"PR": "PR", "p01": "p01_kPa (·1000 Па)", "T01": "T01_C (+273,15 К)", "T02": "T02_C (+273,15 К)",
                "m_total": "m_total_kg_s", "m_side": "m_per_side_kg_s" if "m_per_side_kg_s" in s else "m_total_kg_s / сторон",
                "n_sides": "inlet («двусторонний» → 2)", "rpm": "n_rpm", "r2": "D2_m / 2",
                "c_theta2": "c2_m_s · sin(alpha2_flow_from_radial_deg)", "Z": "blades_full",
                "eta_is": "eta_is", "eta_p": "eta_p"},
    )


def resolve_mesh(a) -> dict:
    """Пути сетки и JSON сетки: из сводки узла геометрии (--geom-results, --stage, --level) или прямо
    (--mesh, --mesh-json). Сверяет SHA-256 файла сетки с записанной в сводке, если она есть."""
    out = {"source": None, "sha256_expected": None}
    if a.geom_results:
        g = json.load(open(a.geom_results, encoding="utf-8"))
        try:
            L = g["stages"][a.stage]["levels"][a.level]
        except KeyError:
            raise SystemExit(f"в {a.geom_results} нет stages.{a.stage}.levels.{a.level}") from None
        f = L["files"]
        ps = f["passage_su2"]
        stored = ps.get("stored_as")
        if os.path.isfile(ps["path"]):
            out.update(mesh=ps["path"], sha256_expected=ps.get("sha256"))
        elif stored and os.path.isfile(stored["path"]):
            out.update(mesh=stored["path"], sha256_expected=stored.get("sha256"))
        else:
            raise SystemExit(f"файла сетки нет: {ps['path']} (и сжатой копии) — постройте её узлом геометрии "
                             f"(routes/gtu/frame3/frame3_geom.py --levels {a.level})")
        out.update(mesh_json=f["impeller_json"]["path"], mesh_json_sha256_expected=f["impeller_json"].get("sha256"),
                   source=f"{rel(a.geom_results)} → stages.{a.stage}.levels.{a.level}.files",
                   geom_check=L.get("check_su2", {}), n_nodes_geom=L.get("n_nodes"), n_cells_geom=L.get("n_cells"))
    else:
        if not a.mesh:
            raise SystemExit("нужна сетка: --geom-results с --stage/--level или --mesh с --mesh-json")
        js = a.mesh_json
        if not js:
            d = os.path.dirname(os.path.abspath(a.mesh))
            cands = sorted(x for x in os.listdir(d) if x.startswith("impeller_") and x.endswith(".json"))
            if len(cands) != 1:
                raise SystemExit(f"JSON сетки не найден однозначно рядом с {a.mesh}: {cands} — укажите --mesh-json")
            js = os.path.join(d, cands[0])
        out.update(mesh=a.mesh, mesh_json=js, source="--mesh/--mesh-json")
    for k in ("mesh", "mesh_json"):
        out[k] = os.path.abspath(out[k])
        if not os.path.isfile(out[k]):
            raise SystemExit(f"нет файла {out[k]}")
    if out.get("sha256_expected") and not a.no_sha:
        got = sha256(out["mesh"])
        out["sha256"] = got
        if got != out["sha256_expected"]:
            raise SystemExit(f"SHA-256 сетки {out['mesh']} не совпал со сводкой геометрии ({got[:12]}… ≠ "
                             f"{out['sha256_expected'][:12]}…)")
    return out


def stage_mesh(src: str, run_dir: str) -> str:
    """Сетка в папке прогона под коротким именем mesh.su2: .gz распаковывается, иначе — ссылка (копия, если
    ссылки нельзя). Короткие имена — потому что SU2 v8.6.0 падал на длинном пути в SOLUTION_FILENAME."""
    dst = os.path.join(run_dir, "mesh.su2")
    if os.path.lexists(dst):
        os.remove(dst)
    if src.endswith(".gz"):
        with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo, 1 << 20)
    else:
        try:
            os.symlink(src, dst)
        except OSError:
            shutil.copyfile(src, dst)
    return dst


def consistency(pt: dict, mj: dict, hdr: dict) -> dict:
    """Сверка точки и сетки: частота, полные параметры входа, расход канала, число каналов, маркеры."""
    s2, d1 = mj["su2"], mj["d1"]
    per = s2["periodic"]
    omega_mesh = s2["rotation_rate_rad_s"]
    n_pass = int(round(360.0 / per["angle_deg"]))
    c = {
        "omega_mesh_rad_s": omega_mesh[2], "omega_point_rad_s": pt["omega"],
        "omega_rel_diff": omega_mesh[2] / pt["omega"] - 1.0,
        "rotation_axis_z": abs(omega_mesh[0]) == 0 and abs(omega_mesh[1]) == 0 and list(per["axis"]) == [0, 0, 1],
        "p01_rel_diff": d1["p01"] / pt["p01"] - 1.0, "T01_rel_diff": d1["T01"] / pt["T01"] - 1.0,
        "n_passages": n_pass, "Z_point": pt["Z"], "n_passages_eq_Z": n_pass == pt["Z"],
        "mesh_double_sided": bool(d1.get("double_sided")), "point_double_sided": pt["double_sided"],
        "m_passage_point_kg_s": pt["m_side"] / n_pass,
        "m_passage_mesh_json_kg_s": s2.get("mass_flow_passage_kg_s"),
        "markers_present": sorted(hdr.get("markers", {})),
    }
    need = list(MARKERS) + [per["master"], per["slave"]]
    c["markers_missing"] = [m for m in need if m not in hdr.get("markers", {})]
    c["ok"] = bool(abs(c["omega_rel_diff"]) < 1e-6 and c["rotation_axis_z"] and abs(c["p01_rel_diff"]) < 1e-9
                   and abs(c["T01_rel_diff"]) < 1e-9 and c["n_passages_eq_Z"] and not c["markers_missing"]
                   and c["mesh_double_sided"] == c["point_double_sided"])
    return c


def estimate_p_out(pt: dict, frac: float, r_out: float, b_out: float, gamma: float = 1.4, R: float = 287.058) -> float:
    """Статическое давление на выходе области (r_out, ширина b_out) по расчётной точке ступени — оценка для
    граничного условия, а не результат. p02 = p01·(1 + η_эф·(T02/T01 − 1))^(γ/(γ−1)), η_эф — из PR и T02 точки;
    при доле частоты frac работа ∝ frac², расход ∝ frac, закрутка ∝ frac; Cθ = Cθ2·r2/r_out (сохранение
    момента в безлопаточном участке), Cm — из неразрывности расхода одной стороны колеса (m_side).
    Идеальный газ с постоянным γ, как в конфигурации SU2. PR точки — ступени (с диффузором), поэтому для
    выхода колеса на r_out оценка занижена на потери диффузора."""
    cp, k = gamma * R / (gamma - 1), gamma / (gamma - 1)
    T01, p01 = pt["T01"], pt["p01"]
    tau = pt["T02"] / T01
    eta_eff = (pt["PR"] ** (1 / k) - 1) / (tau - 1)
    T02 = T01 * (1 + (tau - 1) * frac ** 2)
    p02 = p01 * (1 + eta_eff * (T02 / T01 - 1)) ** k
    ct = pt["c_theta2"] * frac * pt["r2"] / r_out
    m = pt["m_side"] * frac
    A = 2 * math.pi * r_out * b_out
    cm = 0.0
    for _ in range(200):
        T = T02 - (ct ** 2 + cm ** 2) / (2 * cp)
        if T <= 0:
            raise ValueError("оценка p_вых: дозвукового решения на r_out нет (запирание)")
        p = p02 * (T / T02) ** k
        cm_new = m / (p / (R * T) * A)
        if abs(cm_new - cm) < 1e-7:
            if cm_new ** 2 + ct ** 2 >= gamma * R * T:
                raise ValueError("оценка p_вых: на r_out сверхзвук (M ≥ 1)")
            return p
        cm = cm_new
    raise ValueError("оценка p_вых не сошлась за 200 итераций (близко к запиранию)")


# ----------------------------------------------------------------------------------------------
# план стадий
# ----------------------------------------------------------------------------------------------
def parse_points(text: str) -> list:
    """'0.75:600,1.0:3000' → [(0.75, 600), (1.0, 3000)]: доля оценки давления на выходе и итераций точки."""
    pts = []
    for item in (text or "").split(","):
        item = item.strip()
        if not item:
            continue
        try:
            f, n = item.split(":")
            pts.append((float(f), int(n)))
        except ValueError:
            raise SystemExit(f"--points: '{item}' — ожидается ДОЛЯ:ИТЕРАЦИЙ, например 0.75:600") from None
    if not pts:
        raise SystemExit("--points: нужна хотя бы одна точка")
    for f, n in pts:
        if not (0.2 <= f <= 2.0) or n < 1:
            raise SystemExit(f"--points: доля {f} вне [0,2; 2,0] или итераций {n} < 1")
    return pts


def make_plan(speed: float, points: list, startup: bool = True, iter_scale: float = 1.0, point_pramp: int = 100) -> list:
    """Стадии: разгон (s1: до 0,6 частоты линии, давление 0,95 оценки при этой частоте; s2: до частоты линии и
    давления первой точки; первый порядок, давление растёт вдвое дольше частоты — README узла, диагностика), затем
    точки по противодавлению вторым порядком; давление новой точки — рампой point_pramp итераций от прежнего.
    Длины рамп и итераций масштабируются iter_scale (кроме итераций точек, они заданы явно)."""
    sc = lambda n: max(1, int(round(n * iter_scale)))  # noqa: E731
    plan = []
    if startup:
        plan += [dict(name="s1_ramp", frac=0.6 * speed, pfrac=0.95, ramp=sc(80), pramp=sc(160), iter=sc(170), muscl=False),
                 dict(name="s2_ramp", frac=speed, pfrac=points[0][0], ramp=sc(80), pramp=sc(160), iter=sc(200), muscl=False)]
    for i, (pf, n) in enumerate(points):
        pr = 0 if (i == 0 and startup) else point_pramp
        if pr >= n:
            pr = max(0, n - 1)
        plan.append(dict(name=f"p{i + 1}_{pf:.3f}".replace(".", "_"), frac=speed, pfrac=pf, ramp=0, pramp=pr,
                         iter=n, muscl=True))
    return RC.validate_plan(plan)


# ----------------------------------------------------------------------------------------------
# конфигурация
# ----------------------------------------------------------------------------------------------
def g10(x: float) -> str:
    return f"{x:.10g}"


def walls(shroud: str, vaneless: str) -> dict:
    """Стенки: hub и blade — адиабатные с прилипанием, вращаются с колесом (MARKER_HEATFLUX); покрышка
    stationary — неподвижна (MARKER_HEATFLUX + MARKER_SHROUD), rotating — вращается, slip — без трения
    (MARKER_EULER); стенки за колесом (r ≥ r2) — как в сетке (as-is) или отдельными маркерами *_vaneless."""
    heat, still, euler = ["hub", "blade"], [], []
    if shroud == "stationary":
        heat.append("shroud"); still.append("shroud")
    elif shroud == "rotating":
        heat.append("shroud")
    elif shroud == "slip":
        euler.append("shroud")
    else:
        raise ValueError(f"shroud = {shroud}")
    if vaneless != "as-is":
        vl = ["hub_vaneless", "shroud_vaneless"]
        if vaneless == "stationary":
            heat += vl; still += vl
        elif vaneless == "slip":
            euler += vl
        else:
            raise ValueError(f"vaneless = {vaneless}")
    return dict(heat=heat, still=still, euler=euler)


def build_config(c: dict) -> list:
    """Конфигурация SU2 одной стадии — список (ключ, значение). c: словарь параметров (stage_params).
    Имена ключей сверяются с регистрацией в Common/src/CConfig.cpp v8.6.0 (check_keys), значения — SU2_CFD -d."""
    n = dict(NUMERICS)
    n.update(c.get("numerics", {}))
    gamma = float(n["gamma"])
    kv = []

    def add(k, v):
        kv.append((k, str(v)))

    add("SOLVER", "RANS")
    add("KIND_TURB_MODEL", n["turb"])
    if n["turb"] == "SA":
        add("SA_OPTIONS", "NONE")
    add("MATH_PROBLEM", "DIRECT")
    add("RESTART_SOL", "YES" if c.get("restart") else "NO")
    add("SOLUTION_FILENAME", c.get("restart") or "solution.dat")
    # начальное поле: число Маха вдоль +z (AOA= 90), статические p, T — из p01, T01 при этом числе Маха
    m0 = float(n["init_mach"])
    fs_t = c["T01"] / (1 + 0.5 * (gamma - 1) * m0 ** 2)
    fs_p = c["p01"] * (fs_t / c["T01"]) ** (gamma / (gamma - 1))
    add("MACH_NUMBER", g10(m0))
    add("AOA", g10(float(n["aoa"])))
    add("SIDESLIP_ANGLE", "0.0")
    add("FREESTREAM_OPTION", "TEMPERATURE_FS")
    add("FREESTREAM_PRESSURE", g10(fs_p))
    add("FREESTREAM_TEMPERATURE", g10(fs_t))
    add("FREESTREAM_TURBULENCEINTENSITY", g10(float(n["tu"])))
    add("FREESTREAM_TURB2LAMVISCRATIO", g10(float(n["mut_ratio"])))
    add("INIT_OPTION", "TD_CONDITIONS")
    add("REF_DIMENSIONALIZATION", "DIMENSIONAL")
    # газ
    add("FLUID_MODEL", "IDEAL_GAS")
    add("GAMMA_VALUE", g10(gamma))
    add("GAS_CONSTANT", g10(float(n["R"])))
    add("VISCOSITY_MODEL", "SUTHERLAND")
    add("MU_REF", g10(float(n["mu_ref"])))
    add("MU_T_REF", g10(float(n["mu_t_ref"])))
    add("SUTHERLAND_CONSTANT", g10(float(n["sutherland"])))
    add("CONDUCTIVITY_MODEL", "CONSTANT_PRANDTL")
    # вращающаяся система (ось z; рампа частоты меняет только z-компоненту — код CConfig::SetRate)
    add("GRID_MOVEMENT", "ROTATING_FRAME")
    add("MOTION_ORIGIN", "0.0 0.0 0.0")
    add("ROTATION_RATE", f"0.0 0.0 {g10(c['omega'])}")
    add("MACH_MOTION", "0.35")
    add("RAMP_MOTION_FRAME", "YES" if c.get("ramp_rot") else "NO")
    add("RAMP_MOTION_FRAME_COEFF", f"({g10(c.get('omega_ini', 0.0))}, {int(n['ramp_freq'])}, {max(int(c.get('ramp_rot', 0)), 1)})")
    # стенки
    w = walls(c["shroud"], c["vaneless"])
    add("MARKER_HEATFLUX", "(" + ", ".join(f"{m}, 0.0" for m in w["heat"]) + ")")
    if w["still"]:
        add("MARKER_SHROUD", "(" + ", ".join(w["still"]) + ")")
    if w["euler"]:
        add("MARKER_EULER", "(" + ", ".join(w["euler"]) + ")")
    # периодичность: второй маркер = первый, повёрнутый на +angle вокруг +z (узел в узел, CGeometry::MatchPeriodic)
    add("MARKER_PERIODIC", f"({c['per_master']}, {c['per_slave']}, 0.0, 0.0, 0.0, 0.0, 0.0, {g10(c['angle'])}, 0.0, 0.0, 0.0)")
    # вход и выход: Giles на турбомашинных маркерах; вход — p0, T0, направление (1, 0, 0) в осях решётки (по
    # нормали, без закрутки), выход — статическое давление (STATIC_PRESSURE_1D)
    add("MARKER_TURBOMACHINERY", "(inlet, outlet)")
    add("TURBOMACHINERY_KIND", "AXIAL_CENTRIFUGAL")
    add("TURBO_PERF_KIND", "COMPRESSOR")
    add("MARKER_GILES", f"(inlet, TOTAL_CONDITIONS_PT, {g10(c['p01'])}, {g10(c['T01'])}, 1.0, 0.0, 0.0, "
                        f"{g10(float(n['relax_in']))}, 0.0, outlet, STATIC_PRESSURE_1D, {g10(c['p_out'])}, 0.0, 0.0, "
                        f"0.0, 0.0, {g10(float(n['relax_out']))}, 0.0)")
    add("RAMP_OUTLET", "YES" if c.get("ramp_p") else "NO")
    add("RAMP_OUTLET_COEFF", f"({g10(c.get('p_ini', c['p_out']))}, {int(n['ramp_freq'])}, {max(int(c.get('ramp_p', 0)), 1)})")
    add("AVERAGE_PROCESS_KIND", "MIXEDOUT")
    add("PERFORMANCE_AVERAGE_PROCESS_KIND", "MIXEDOUT")
    add("SPANWISE_KIND", "EQUISPACED")
    add("NUM_SPANWISE_SECTIONS", int(n["nspan"]))
    # средние по расходу на входе и выходе — столбцы Avg_*(inlet/outlet) истории
    add("MARKER_ANALYZE", "(inlet, outlet)")
    add("MARKER_ANALYZE_AVERAGE", "MASSFLUX")
    # численная схема
    add("NUM_METHOD_GRAD", "WEIGHTED_LEAST_SQUARES")
    add("CONV_NUM_METHOD_FLOW", n["scheme"])
    add("MUSCL_FLOW", "YES" if c.get("muscl") else "NO")
    add("SLOPE_LIMITER_FLOW", n["limiter"])
    add("VENKAT_LIMITER_COEFF", g10(float(n["venkat"])))
    add("ENTROPY_FIX_COEFF", g10(float(n["entropy_fix"])))
    add("CONV_NUM_METHOD_TURB", "SCALAR_UPWIND")
    add("MUSCL_TURB", "NO")
    add("TIME_DISCRE_FLOW", "EULER_IMPLICIT")
    add("TIME_DISCRE_TURB", "EULER_IMPLICIT")
    add("CFL_NUMBER", g10(float(n["cfl"])))
    add("CFL_ADAPT", n["cfl_adapt"])
    add("CFL_ADAPT_PARAM", f"({n['cfl_adapt_param']})")
    add("LINEAR_SOLVER", "FGMRES")
    add("LINEAR_SOLVER_PREC", n["lin_prec"])
    add("LINEAR_SOLVER_ERROR", g10(float(n["lin_error"])))
    add("LINEAR_SOLVER_ITER", int(n["lin_iter"]))
    add("MGLEVEL", 0)  # Giles — только на самой мелкой сетке (код)
    # останов
    add("ITER", int(c["iter"]))
    add("CONV_FIELD", "RMS_DENSITY")
    add("CONV_RESIDUAL_MINVAL", g10(float(n["res_min"])))
    add("CONV_STARTITER", 10)
    # вывод
    add("SCREEN_OUTPUT", "(INNER_ITER, WALL_TIME, RMS_DENSITY, RMS_ENERGY, RMS_TKE, SURFACE_MASSFLOW, SURFACE_TOTAL_PRESSURE)")
    add("HISTORY_OUTPUT", "(ITER, WALL_TIME, RMS_RES, FLOW_COEFF, FLOW_COEFF_SURF, TURBO_PERF)")
    add("SCREEN_WRT_FREQ_INNER", 10)
    add("VOLUME_OUTPUT", "(COORDINATES, SOLUTION, PRIMITIVE)")
    add("OUTPUT_FILES", f"({n['output_files']})")
    add("OUTPUT_WRT_FREQ", int(n["wrt_freq"]))
    add("TABULAR_FORMAT", "CSV")
    add("CONV_FILENAME", f"history_{c['name']}")
    add("RESTART_FILENAME", f"restart_{c['name']}.dat")
    add("VOLUME_FILENAME", f"vol_{c['name']}")
    add("MESH_FILENAME", c["mesh"])
    add("MESH_FORMAT", "SU2")
    for line in c.get("extra", []):
        k, _, v = line.partition("=")
        kv.append((k.strip(), v.strip()))
    keys = [k for k, _ in kv]
    dup = sorted({k for k in keys if keys.count(k) > 1})
    if dup:
        raise ValueError(f"ключи повторены: {dup}")
    return kv


def cfg_text(kv: list, header: str = "") -> str:
    lines = [f"% {ln}" for ln in header.splitlines()] if header else []
    lines += [f"{k}= {v}" for k, v in kv]
    return "\n".join(lines) + "\n"


# ----------------------------------------------------------------------------------------------
# разбор истории
# ----------------------------------------------------------------------------------------------
COL = dict(m_in="Avg_Massflow(inlet)", m_out="Avg_Massflow(outlet)", pt_in="Avg_TotalPress(inlet)",
           pt_out="Avg_TotalPress(outlet)", tt_in="Avg_TotalTemp(inlet)", tt_out="Avg_TotalTemp(outlet)")


def _fin(x) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(x)


def perf_row(row: dict, gamma: float) -> dict:
    """π, T0-отношение и изоэнтропический КПД полн./полн. по средним по расходу p0, T0 одной строки истории:
    η = (π^((γ−1)/γ) − 1)/(T02/T01 − 1) — идеальный газ с постоянным γ (как IDEAL_GAS конфигурации)."""
    pt1, pt2, tt1, tt2 = (row.get(COL[k]) for k in ("pt_in", "pt_out", "tt_in", "tt_out"))
    if not all(_fin(x) for x in (pt1, pt2, tt1, tt2)) or pt1 <= 0 or tt1 <= 0:
        return {}
    pr, tr = pt2 / pt1, tt2 / tt1
    eta = (pr ** ((gamma - 1) / gamma) - 1) / (tr - 1) if tr > 1 + 1e-12 else None
    return {"PR": pr, "TR": tr, "eta": eta}


def read_rows(path: str) -> list:
    head, d = RC.read_history(path)
    n = len(d.get("Inner_Iter", next(iter(d.values()), [])))
    return [{h: d[h][i] for h in head} for i in range(n)]


def rel_change(vals: list):
    v = [x for x in vals if _fin(x)]
    if len(v) < len(vals) or not v or v[-1] == 0:
        return None
    return (max(v) - min(v)) / abs(v[-1])


def analyse_stage(rows: list, gamma: float, n_passages: int, n_sides: int, run_max: dict, conv: dict) -> dict:
    """Разбор одной стадии: невязки (начало, максимум, конец; падение от максимума всего прогона до этой стадии
    включительно), расход и дисбаланс, π и η на последней строке, их изменение в окне conv['window'] итераций
    стадии и вывод о сходимости по критериям conv."""
    if not rows:
        return {"rows": 0}
    out = {"rows": len(rows)}
    res = {}
    for h in rows[0]:
        if not h.startswith("rms["):
            continue
        v = [r[h] for r in rows]
        fin = [x for x in v if math.isfinite(x)]
        if not fin:
            res[h] = {"nan": True}
            continue
        run_max[h] = max(run_max.get(h, -1e30), max(fin))
        res[h] = {"first": v[0], "max_stage": max(fin), "last": v[-1], "max_run": run_max[h],
                  "drop_from_run_max": run_max[h] - v[-1] if math.isfinite(v[-1]) else None,
                  "nan": len(fin) < len(v)}
    out["residuals"] = res
    last = rows[-1]
    m_in, m_out = last.get(COL["m_in"]), last.get(COL["m_out"])
    perf = {"m_in_passage": m_in, "m_out_passage": m_out}
    if _fin(m_in) and _fin(m_out) and m_in != 0:
        perf["mass_imbalance"] = abs(m_in + m_out) / abs(m_in)
        perf["mass_sign_ok"] = bool(m_in > 0 > m_out)
        perf["m_side"] = m_in * n_passages
        perf["m_wheel"] = m_in * n_passages * n_sides
    pr = perf_row(last, gamma)
    if pr:
        perf.update(PR_tt=pr["PR"], TR_tt=pr["TR"], eta_tt_isentropic=pr["eta"],
                    pt_in=last[COL["pt_in"]], pt_out=last[COL["pt_out"]], tt_in=last[COL["tt_in"]], tt_out=last[COL["tt_out"]])
        if _fin(m_in):
            perf["m_wheel_corrected"] = m_in * n_passages * n_sides * math.sqrt(last[COL["tt_in"]] / 288.15) / (last[COL["pt_in"]] / 101325.0)
    for col, key in (("PRTT", "PRTT_su2"), ("TotTotEff", "TotTotEff_su2"), ("PRTS", "PRTS_su2"),
                     ("TotStaticEff", "TotStaticEff_su2")):
        if col in last:
            perf[key] = last[col]
    out["performance"] = perf
    win = int(conv["window"])
    seg = rows[-win:] if len(rows) >= win else None
    stab = {"window": win, "rows_in_stage": len(rows)}
    if seg:
        ps = [perf_row(r, gamma) for r in seg]
        stab["PR_rel_change"] = rel_change([p.get("PR") for p in ps]) if all(ps) else None
        stab["eta_rel_change"] = rel_change([p.get("eta") if p and p.get("eta") is not None else float("nan") for p in ps])
        stab["m_in_rel_change"] = rel_change([r.get(COL["m_in"]) for r in seg])
    out["stability"] = stab
    for tcol in ("Time(sec)", "Cur_Time"):
        if tcol in last:
            ts = [r[tcol] for r in rows if _fin(r[tcol])]
            # Time(sec) у SU2 v8 — время итерации (WALL_TIME); медиана по стадии устойчива к выбросам загрузки
            if ts:
                srt = sorted(ts)
                out["sec_per_iter_median"] = srt[len(srt) // 2]
                out["sec_per_iter_p10"] = srt[len(srt) // 10]  # нижний дециль: ближе к времени без чужой нагрузки
                out["sec_per_iter_mean"] = sum(ts) / len(ts)
            break
    out["convergence"] = verdict(out, conv)
    return out


def verdict(st: dict, conv: dict) -> dict:
    """Сошлась ли точка: нет NaN, прямое течение, rms[Rho] упала от максимума прогона ≥ conv['orders'] порядков,
    дисбаланс расхода ≤ conv['mass_tol'], π и η за последние conv['window'] итераций стадии меняются
    меньше conv['stab_tol'] (навык su2-turbo, раздел 6)."""
    rho = st.get("residuals", {}).get("rms[Rho]", {})
    perf = st.get("performance", {})
    stab = st.get("stability", {})
    drop = rho.get("drop_from_run_max")
    v = {
        "no_nan": not any(r.get("nan") for r in st.get("residuals", {}).values()),
        "forward_flow": bool(perf.get("mass_sign_ok")),
        "rms_rho_drop_orders": drop,
        "residual_drop_ok": bool(_fin(drop) and drop >= conv["orders"]),
        "mass_imbalance": perf.get("mass_imbalance"),
        "mass_ok": bool(_fin(perf.get("mass_imbalance")) and perf["mass_imbalance"] <= conv["mass_tol"]),
        "PR_window_rel_change": stab.get("PR_rel_change"),
        "eta_window_rel_change": stab.get("eta_rel_change"),
        "stability_ok": bool(_fin(stab.get("PR_rel_change")) and _fin(stab.get("eta_rel_change"))
                             and stab["PR_rel_change"] < conv["stab_tol"] and stab["eta_rel_change"] < conv["stab_tol"]),
        "criteria": dict(conv),
    }
    if not stab.get("PR_rel_change") and st.get("rows", 0) < conv["window"]:
        v["stability_note"] = f"в стадии {st.get('rows', 0)} итераций < окна {conv['window']} — устойчивость не оценивается"
    v["converged"] = bool(v["no_nan"] and v["forward_flow"] and v["residual_drop_ok"] and v["mass_ok"] and v["stability_ok"])
    return v


def compare_with_point(perf: dict, pt: dict) -> dict:
    """Числа SU2 против расчётной точки — сравнение, не проверка: точка — одномерная модель ступени (с
    диффузором), SU2 — колесо с безлопаточным участком до r = 1,3 r2."""
    o = {"note": "сравнение, не проверка: точка — одномерная модель ступени с диффузором (frame3_point.py), "
                 "SU2 — колесо с безлопаточным участком до r_out; одномерную модель под SU2 не подгонять"}
    if _fin(perf.get("m_wheel")):
        o.update(m_wheel_su2=perf["m_wheel"], m_wheel_point=pt["m_total"], m_ratio=perf["m_wheel"] / pt["m_total"])
    if _fin(perf.get("PR_tt")):
        o.update(PR_su2=perf["PR_tt"], PR_point_stage=pt["PR"], PR_ratio=perf["PR_tt"] / pt["PR"],
                 TR_su2=perf["TR_tt"], TR_point=pt["T02"] / pt["T01"])
    if _fin(perf.get("eta_tt_isentropic")):
        o.update(eta_su2=perf["eta_tt_isentropic"], eta_point_is_stage=pt.get("eta_is"),
                 d_eta_points=(perf["eta_tt_isentropic"] - pt["eta_is"]) if _fin(pt.get("eta_is")) else None)
    return o


# ----------------------------------------------------------------------------------------------
# запуск
# ----------------------------------------------------------------------------------------------
def tool_versions(su2: str | None, src: str | None) -> dict:
    v = {"python": platform.python_version(), "os": platform.platform()}
    if su2:
        v["su2"] = RC.su2_version(su2)
        v["su2_exe"] = su2
    if src and os.path.isdir(os.path.join(src, ".git")):
        try:
            v["su2_src_commit"] = subprocess.run(["git", "-C", src, "rev-parse", "HEAD"], capture_output=True,
                                                 text=True, timeout=30).stdout.strip()
            v["su2_src_tag"] = subprocess.run(["git", "-C", src, "describe", "--tags", "--exact-match"],
                                              capture_output=True, text=True, timeout=30).stdout.strip() or None
        except Exception:  # noqa: BLE001
            pass
    mpirun = shutil.which("mpirun")
    if mpirun:
        try:
            out = subprocess.run([mpirun, "--version"], capture_output=True, text=True, timeout=30).stdout
            v["mpirun"] = out.splitlines()[0].strip() if out else None
        except Exception:  # noqa: BLE001
            pass
    return v


def stage_params(base: dict, st: dict, prev: dict, pt: dict, p_target: float) -> dict:
    c = dict(base)
    c.update(name=st["name"], iter=st["iter"], muscl=bool(st.get("muscl")), restart=prev["restart"],
             omega=pt["omega"] * st["frac"], p_out=p_target)
    if st.get("ramp", 0) > 0 and abs(st["frac"] - prev["frac"]) > 1e-9:
        c.update(ramp_rot=st["ramp"], omega_ini=pt["omega"] * prev["frac"])
    pramp = int(st.get("pramp", st.get("ramp", 0)))
    if pramp > 0 and abs(p_target - prev["p"]) > 1.0:
        c.update(ramp_p=pramp, p_ini=prev["p"])
    return c


def run(a) -> int:
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(sig, RC._kill_child)
    conv = dict(CONV, orders=a.conv_orders, mass_tol=a.mass_tol, window=a.window, stab_tol=a.stab_tol)
    pt = load_point(a.point)
    mesh = resolve_mesh(a)
    mj = json.load(open(mesh["mesh_json"], encoding="utf-8"))
    out = os.path.abspath(a.out or os.path.join(HERE, "runs", f"case_{pt['key'] or 'point'}"))
    os.makedirs(os.path.join(out, "TURBOMACHINERY"), exist_ok=True)  # профили по высоте (док.: папку SU2 не создаёт)
    t_mesh = time.monotonic()
    mesh_local = stage_mesh(mesh["mesh"], out)
    hdr = RC.read_mesh_header(mesh_local)
    chk = consistency(pt, mj, hdr)
    if not chk["ok"] and not a.force:
        raise SystemExit(f"точка и сетка не согласованы (--force — считать всё равно): "
                         f"{json.dumps(chk, ensure_ascii=False)}")
    split = None
    if a.vaneless != "as-is":
        split_path = os.path.join(out, "mesh_split.su2")
        split = RC.split_vaneless(mesh_local, split_path, float(mj["d1"]["r2"]))
        mesh_local = split_path
    t_mesh = time.monotonic() - t_mesh
    per = mj["su2"]["periodic"]
    n_pass = chk["n_passages"]
    gamma, R = float(NUMERICS["gamma"]), float(NUMERICS["R"])
    r_out = float(mj.get("contours", {}).get("r_out", 1.3 * float(mj["d1"]["r2"])))
    b_out = float(mj["d1"]["b2"])
    points = parse_points(a.points)
    startup = not a.restart_from
    plan = make_plan(a.speed, points, startup=startup, iter_scale=a.iter_scale, point_pramp=a.point_pramp)
    try:
        p_est_design = estimate_p_out(pt, 1.0, r_out, b_out, gamma, R)
        p_est = {round(st["frac"], 6): estimate_p_out(pt, st["frac"], r_out, b_out, gamma, R) for st in plan}
    except ValueError as e:
        raise SystemExit(f"{e}: оценка давления на выходе для этой точки не годится — задайте --p-out-design")
    if a.p_out_design:  # своё давление на выходе при расчётной частоте: оценки всех частот — в той же пропорции
        k = a.p_out_design / p_est_design
        p_est = {f: p * k for f, p in p_est.items()}
        p_est_design = a.p_out_design

    numerics = {}
    for kv in a.set:
        k, _, v = kv.partition("=")
        numerics[k.strip()] = v.strip()
    if a.output_files:
        numerics["output_files"] = a.output_files
    base = dict(T01=pt["T01"], p01=pt["p01"], per_master=per["master"], per_slave=per["slave"],
                angle=float(per["angle_deg"]), shroud=a.shroud, vaneless=a.vaneless, mesh=os.path.basename(mesh_local),
                numerics=numerics, extra=list(a.extra))
    su2 = RC.find_su2(a.su2_bin)
    src = a.su2_src or os.environ.get("SU2_HOME")
    cconfig = RC.find_cconfig(src)
    summ = {
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "software": tool_versions(su2, src), "np": a.np, "mpi_args": a.mpi_args,
        "point": {k: v for k, v in pt.items()},
        "mesh": {"source": mesh["source"], "file": mesh["mesh"], "sha256": mesh.get("sha256"),
                 "sha256_expected": mesh.get("sha256_expected"), "json": mesh["mesh_json"],
                 "header": hdr, "staged_as": os.path.basename(mesh_local), "vaneless_split_faces": split,
                 "t_prepare_s": round(t_mesh, 1), "n_nodes_geom": mesh.get("n_nodes_geom")},
        "consistency": chk,
        "setup": {"n_passages": n_pass, "n_sides": pt["n_sides"], "r_out": r_out, "b_out": b_out, "gamma": gamma,
                  "R": R, "speed_frac": a.speed, "rpm": pt["rpm"] * a.speed, "p_out_estimate_design": p_est_design,
                  "p_out_estimate_by_frac": {str(f): p for f, p in p_est.items()},
                  "p_out_design_given": a.p_out_design, "shroud": a.shroud, "vaneless_walls": a.vaneless,
                  "points": points, "time_limit_s": a.time_limit, "restart_from": a.restart_from,
                  "restart_state": a.restart_state, "numerics": {**NUMERICS, **numerics}, "extra": list(a.extra)},
        "convergence_criteria": conv, "stages": [],
    }
    deadline = time.monotonic() + a.time_limit
    m0 = float(summ["setup"]["numerics"]["init_mach"])
    fs_t = pt["T01"] / (1 + 0.5 * (gamma - 1) * m0 ** 2)
    fs_p = pt["p01"] * (fs_t / pt["T01"]) ** (gamma / (gamma - 1))
    # рампа давления стартует со статического давления начального поля, а не с p01 (README узла, диагностика)
    prev = {"frac": 0.0, "p": fs_p, "restart": None}
    if a.restart_from:
        srcf = os.path.abspath(a.restart_from)
        if not os.path.isfile(srcf):
            raise SystemExit(f"нет файла рестарта {srcf}")
        shutil.copyfile(srcf, os.path.join(out, "restart_in.dat"))
        prev["restart"] = "restart_in.dat"
        if a.restart_state:
            try:
                f0, p0 = (float(x) for x in a.restart_state.split(","))
            except ValueError:
                raise SystemExit("--restart-state: ДОЛЯ,ДАВЛЕНИЕ, например 1.0,253000") from None
            prev.update(frac=f0, p=p0)
        else:
            prev.update(frac=plan[0]["frac"], p=plan[0]["pfrac"] * p_est[round(plan[0]["frac"], 6)])
    key_check = None
    for st in plan:
        p_target = st["pfrac"] * p_est[round(st["frac"], 6)]
        c = stage_params(base, st, prev, pt, p_target)
        kv = build_config(c)
        hdr_txt = (f"SU2 v8.6.0: канал колеса, стадия {st['name']}; собрано su2_case.py {summ['created_utc']}\n"
                   f"точка: {pt['file']}:{pt['key']}; сетка: {rel(mesh['mesh'])}")
        txt = cfg_text(kv, hdr_txt)
        cfg_name = f"case_{st['name']}.cfg"
        with open(os.path.join(out, cfg_name), "w", encoding="utf-8") as fh:
            fh.write(txt)
        if cconfig and key_check is None:
            unknown, nknown = RC.check_keys(txt, cconfig)
            key_check = {"cconfig": cconfig, "registered_keys": nknown, "unknown": unknown}
            if unknown:
                raise SystemExit(f"ключей нет в {cconfig}: {unknown}")
        rec = {"name": st["name"], "cfg": cfg_name, "frac": st["frac"], "rpm": pt["rpm"] * st["frac"],
               "pfrac": st["pfrac"], "p_out": p_target, "p_out_ramp_from": c.get("p_ini"),
               "omega_ramp_from": c.get("omega_ini"), "ramp_rot": c.get("ramp_rot", 0), "ramp_p": c.get("ramp_p", 0),
               "iter": st["iter"], "muscl": bool(st.get("muscl")), "restart_from": prev["restart"],
               "history_file": f"history_{st['name']}.csv", "log_file": f"su2_{st['name']}.log"}
        summ["stages"].append(rec)
        nxt = {"frac": st["frac"], "p": p_target, "restart": f"restart_{st['name']}.dat"}
        if su2 is None:
            rec["status"] = "no_su2"
            prev = nxt
            continue
        if a.dry:
            if prev["restart"]:
                rec["status"] = "cfg_only"  # SU2_CFD -d читал бы рестарт, а его в сухом прогоне нет
            else:
                lp = os.path.join(out, f"dry_{st['name']}.log")
                rc, status, dt = RC.run_limited([su2, "-d", cfg_name], out, lp, deadline, env=RC.su2_env(1))
                ok = "Exit Success" in open(lp, encoding="utf-8", errors="replace").read()
                rec.update(status="dry_ok" if ok else "dry_failed", seconds=round(dt, 1))
            prev = nxt
            continue
        load0 = os.getloadavg() if hasattr(os, "getloadavg") else None
        rc, status, dt = RC.run_limited(RC.mpi_cmd(a.np, su2, cfg_name, a.mpi_args.split()), out,
                                        os.path.join(out, rec["log_file"]), deadline, env=RC.su2_env(a.np))
        # средняя загрузка машины (1, 5, 15 мин) до и после стадии: на общей машине время итерации завышено
        rec.update(returncode=rc, status=status, seconds=round(dt, 1), loadavg_start=load0,
                   loadavg_end=os.getloadavg() if hasattr(os, "getloadavg") else None)
        rec["log"] = RC.parse_log(os.path.join(out, rec["log_file"]))
        print(f"[{st['name']}] {status} за {dt:.0f} с", flush=True)
        if status != "finished" or RC._CHILD["signal"]:
            break
        if not os.path.isfile(os.path.join(out, nxt["restart"])):
            rec["status"] = "no_restart"
            break
        prev = nxt
    summ["key_check"] = key_check or {"note": "исходник SU2 не найден — ключи не сверялись (--su2-src или $SU2_HOME)"}
    if su2 is None:
        summ["blocker"] = "SU2_CFD не найден: --su2-bin, $SU2_RUN или PATH (сборка — install_su2.sh)"
    analyse_run(summ, out)
    write_json(os.path.join(out, "case_summary.json"), summ)
    print_summary(summ)
    print(f"папка прогона: {out}")
    if a.record and not a.dry and su2:
        record(a.record, summ, out)
    if RC._CHILD["signal"]:
        return 128 + RC._CHILD["signal"]
    return 0 if (a.dry or su2) else 2


def analyse_run(summ: dict, out: str) -> None:
    """Разбор всех стадий прогона по их history_*.csv; итог — последняя стадия с историей."""
    s = summ["setup"]
    run_max = {}
    final = None
    for rec in summ["stages"]:
        hp = os.path.join(out, rec["history_file"])
        if not os.path.isfile(hp):
            continue
        rows = read_rows(hp)
        rec["history"] = analyse_stage(rows, s["gamma"], s["n_passages"], s["n_sides"], run_max,
                                       summ["convergence_criteria"])
        if rec["history"].get("rows"):
            rec["vs_point"] = compare_with_point(rec["history"]["performance"], summ["point"])
            final = rec
        rec["log"] = RC.parse_log(os.path.join(out, rec["log_file"]))
    if final:
        h = final["history"]
        summ["final"] = {"stage": final["name"], "status": final.get("status"), "pfrac": final["pfrac"],
                         "p_out": final["p_out"], "rows": h["rows"], "converged": h["convergence"]["converged"],
                         "convergence": h["convergence"], "performance": h["performance"], "vs_point": final["vs_point"]}
    summ["iterations_total"] = sum(r.get("history", {}).get("rows", 0) for r in summ["stages"])
    summ["seconds_total"] = round(sum(r.get("seconds", 0) or 0 for r in summ["stages"]), 1)
    ts = [r["history"].get("sec_per_iter_median") for r in summ["stages"] if r.get("history", {}).get("sec_per_iter_median")]
    if ts:
        rows = [r["history"]["rows"] for r in summ["stages"] if r.get("history", {}).get("sec_per_iter_median")]
        summ["sec_per_iter_weighted_median"] = sum(t * n for t, n in zip(ts, rows)) / sum(rows)


def write_json(path: str, obj) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, indent=1, default=float)


def print_summary(s: dict) -> None:
    sw = s.get("software", {})
    print(f"SU2: {sw.get('su2')}  np={s.get('np')}")
    hdr = s["mesh"]["header"]
    st = s["setup"]
    print(f"сетка: узлов {hdr.get('NPOIN')}, ячеек {hdr.get('NELEM')}; каналов {st['n_passages']} × сторон {st['n_sides']}; "
          f"{st['rpm']:.0f} об/мин; оценка p_вых при расчётной частоте {st['p_out_estimate_design']:.0f} Па")
    for r in s["stages"]:
        line = f"  {r['name']}: {r.get('status')} {r.get('seconds', '')} с, частота {r['frac']:.3f}, p_вых {r['p_out']:.0f} Па"
        h = r.get("history", {})
        if h.get("rows"):
            rho = h["residuals"].get("rms[Rho]", {})
            pf = h["performance"]
            line += f"; итераций {h['rows']}, rms[Rho] {rho.get('first', float('nan')):.2f} → {rho.get('last', float('nan')):.2f}"
            if _fin(pf.get("mass_imbalance")):
                line += f"; m_вх {pf['m_in_passage']:+.4f}, m_вых {pf['m_out_passage']:+.4f} кг/с, дисбаланс {pf['mass_imbalance']:.2%}"
            if _fin(pf.get("PR_tt")):
                e = pf.get("eta_tt_isentropic")
                line += f"; π {pf['PR_tt']:.3f}" + (f", η {e:.3f}" if _fin(e) else "")
        lg = r.get("log", {})
        if lg.get("periodic_matched_points"):
            line += f"; периодичность {lg['periodic_matched_points']} точек, {lg['periodic_max_distance']:.1e}, плохих {lg['bad_match']}"
        if lg.get("error"):
            line += f"; ошибка: {lg['error']}"
        print(line)
    if s.get("final"):
        f = s["final"]
        print(f"итог ({f['stage']}): сошлось — {f['converged']}; сравнение с точкой: "
              + json.dumps({k: v for k, v in f["vs_point"].items() if k != "note"}, ensure_ascii=False, default=float))


def compact_run(summ: dict, out: str) -> dict:
    """Краткий итог прогона для results-JSON (без больших файлов)."""
    keep = ("name", "status", "seconds", "frac", "rpm", "pfrac", "p_out", "p_out_ramp_from", "ramp_rot", "ramp_p",
            "iter", "muscl", "restart_from", "loadavg_start", "loadavg_end")
    return {
        "created_utc": summ["created_utc"], "run_dir": out, "np": summ["np"], "mpi_args": summ["mpi_args"],
        "mesh": {k: summ["mesh"][k] for k in ("source", "file", "sha256", "sha256_expected", "staged_as",
                                               "vaneless_split_faces")}
                | {"NPOIN": summ["mesh"]["header"].get("NPOIN"), "NELEM": summ["mesh"]["header"].get("NELEM"),
                   "markers": summ["mesh"]["header"].get("markers")},
        "consistency": summ["consistency"], "setup": summ["setup"], "key_check": summ.get("key_check"),
        "stages": [{k: r.get(k) for k in keep}
                   | {"rows": r.get("history", {}).get("rows"),
                      "sec_per_iter_median": r.get("history", {}).get("sec_per_iter_median"),
                      "sec_per_iter_p10": r.get("history", {}).get("sec_per_iter_p10"),
                      "residuals": r.get("history", {}).get("residuals"),
                      "performance": r.get("history", {}).get("performance"),
                      "stability": r.get("history", {}).get("stability"),
                      "convergence": r.get("history", {}).get("convergence"),
                      "vs_point": r.get("vs_point"),
                      "log": {k: r.get("log", {}).get(k) for k in ("periodic_matched_points", "periodic_max_distance",
                                                                   "bad_match", "nan", "exit_success", "error")}}
                   for r in summ["stages"]],
        "final": summ.get("final"), "iterations_total": summ.get("iterations_total"),
        "seconds_total": summ.get("seconds_total"), "sec_per_iter_weighted_median": summ.get("sec_per_iter_weighted_median"),
    }


def record(path: str, summ: dict, out: str) -> None:
    """Дописать прогон в results-JSON (поля method, software, date, inputs заполняются при первом прогоне)."""
    try:
        res = json.load(open(path, encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        res = {}
    pt = summ["point"]
    res.setdefault("title", f"SU2: межлопаточный канал колеса — {pt['label']} ({pt['file']}:{pt['key']})")
    res.setdefault("method", (
        "su2_case.py: сжимаемый RANS SU2 v8.6.0 (SST V2003m), вращающаяся система координат (ROTATING_FRAME, ось z), "
        "один канал (основная лопатка + полулопатка) с периодичностью поворотом узел в узел, Giles: полные p0, T0 на "
        "входе по нормали, статическое давление на выходе r_out = 1,3 r2 (оценка по точке × доля); разгон частоты и "
        "давления рампами первым порядком, затем точка вторым порядком (MUSCL); π и η полн./полн. — по средним по "
        "расходу p0 и T0 входа и выхода (MARKER_ANALYZE_AVERAGE= MASSFLUX), идеальный газ γ = 1,4"))
    res["software"] = summ["software"]
    res["date"] = DATE
    res["inputs"] = [
        {"file": pt["file"], "key": pt["key"], "sha256": pt["sha256"], "fields": pt["fields"]},
        {"file": rel(summ["mesh"]["file"]), "sha256": summ["mesh"]["sha256"] or summ["mesh"]["sha256_expected"],
         "source": summ["mesh"]["source"], "note": "сетка канала узла геометрии (не в репозитории — scratchpad)"},
        {"file": rel(summ["mesh"]["json"]), "fields": "su2.periodic, su2.rotation_rate_rad_s, d1.p01, d1.T01, d1.r2, "
                                                      "d1.b2, contours.r_out"},
    ]
    res.setdefault("runs", []).append(compact_run(summ, out))
    write_json(path, res)
    print(f"записано в {path}")


def parse_cmd(a) -> int:
    path = os.path.join(a.run_dir, "case_summary.json")
    summ = json.load(open(path, encoding="utf-8"))
    analyse_run(summ, a.run_dir)
    write_json(path, summ)
    print_summary(summ)
    if a.record:
        record(a.record, summ, os.path.abspath(a.run_dir))
    return 0


# ----------------------------------------------------------------------------------------------
# план прогонов на AWS (только файл задач и оценка; AWS не вызывается)
# ----------------------------------------------------------------------------------------------
CONTAINER_SETUP = [
    "apt-get update -q",
    "apt-get install -y -q --no-install-recommends build-essential python3 python3-venv python3-pip git pkg-config "
    "ca-certificates openmpi-bin libopenmpi-dev libgl1 libglu1-mesa libxrender1 libxcursor1 libxfixes3 libxft2 "
    "libfontconfig1 libxinerama1",
    "python3 -m venv /opt/venv",
    "/opt/venv/bin/pip install -q -r routes/gtu/platform/docker/requirements.txt",
    "bash routes/gtu/platform/su2/install_su2.sh /opt/su2 /opt/SU2-src",
    # тот же коммит тега v8.6.0, что собран и проверен в сессии (git describe — v8.6.0)
    "test \"$(git -C /opt/SU2-src rev-parse HEAD)\" = 7048fd35f56d4bf52e586429d6746a36a10255f8",
]
# OpenMPI 4.1.6 в контейнере docker без особых прав: без OMPI_MCA_osc=^ucx SU2_CFD под mpirun падает в MPI_Win_create,
# без btl_vader_single_copy_mechanism=none нет прямого копирования между процессами (../docker/Dockerfile, там же проверка)
CONTAINER_ENV = "export OMPI_MCA_osc=^ucx OMPI_MCA_btl_vader_single_copy_mechanism=none"
RUNNER_UPLOAD_MIN = 8  # runner.py оставляет 8 мин срока на выгрузку (README узла aws, «Как идёт прогон», п. 6)


def line_points(pfracs: list, first_iter: int, next_iter: int) -> str:
    return ",".join(f"{pf:g}:{first_iter if i == 0 else next_iter}" for i, pf in enumerate(pfracs))


def aws_plan(a) -> int:
    """Пакет задач runner.py: одна задача = одна линия частоты одной ступени (точки по противодавлению идут
    подряд с рестарта предыдущей — навык su2-turbo, раздел 8). Оценка времени: t = t_ref × (узлы/узлы_ref) ×
    (процессы_ref/процессы) / эффективность — та же формула, что batch/cost_estimate.toml."""
    tr = json.load(open(a.trial_results, encoding="utf-8"))
    run0 = tr["runs"][-1]
    kcol = f"sec_per_iter_{a.t_ref_kind}"
    stages_t = [s for s in run0["stages"] if s.get(kcol) and s.get("rows")]
    if a.t_ref:
        t_ref, t_src = a.t_ref, "--t-ref (задано)"
    elif stages_t:
        t_ref = sum(s[kcol] * s["rows"] for s in stages_t) / sum(s["rows"] for s in stages_t)
        t_src = f"{kcol} (Time(sec) = ITER_TIME истории) по стадиям пробного прогона, взвешенное по итерациям"
    else:
        raise SystemExit("в пробном прогоне нет замера времени итерации (--t-ref)")
    n_ref, np_ref = run0["mesh"]["NPOIN"], run0["np"]
    g = json.load(open(a.geom_results, encoding="utf-8"))
    speeds = [float(x) for x in a.speeds.split(",")]
    pfracs = [float(x) for x in a.pfracs.split(",")]
    startup_iter = int(round((170 + 200) * a.iter_scale))
    iters_line = startup_iter + a.first_iter + (len(pfracs) - 1) * a.next_iter
    jobs, est = [], []
    for stage in a.stages.split(","):
        nodes = g["stages"][stage]["levels"][a.level]["n_nodes"]
        t_iter = t_ref * (nodes / n_ref) * (np_ref / a.procs) / a.parallel_eff
        compute_min = iters_line * t_iter / 60.0
        su2_limit_s = int(math.ceil(compute_min * a.time_margin * 60.0))
        max_min = int(math.ceil(a.setup_min + a.mesh_min + su2_limit_s / 60.0 + RUNNER_UPLOAD_MIN + a.slack_min))
        if max_min > 720:
            raise SystemExit(f"{stage}: линия не помещается в 720 мин ({max_min}) — уменьшите итерации или делите линию")
        for sp in speeds:
            name = f"su2-f3-{stage}-n{int(round(sp * 100)):03d}"
            pts = line_points(pfracs, a.first_iter, a.next_iter)
            cmd = (
                "set -e\n"
                "export PATH=/opt/venv/bin:/opt/su2/bin:$PATH SU2_RUN=/opt/su2/bin SU2_HOME=/opt/SU2-src\n"
                f"{CONTAINER_ENV}\n"
                "git rev-parse HEAD > \"$GTU_OUT/commit.txt\"\n"
                "rc=0\n"
                f"python3 routes/gtu/frame3/frame3_geom.py --out \"$GTU_OUT/geom\" --levels {a.level} --stages {stage} "
                f"--gzip-levels --threads \"$NPROC\" --results \"$GTU_OUT/frame3_geom_results.json\"\n"
                f"python3 routes/gtu/platform/su2/su2_case.py run --geom-results \"$GTU_OUT/frame3_geom_results.json\" "
                f"--stage {stage} --level {a.level} --point routes/gtu/frame3/frame3_point_results.json:{stage} "
                f"--speed {sp:g} --points {pts} --iter-scale {a.iter_scale:g} --shroud {a.shroud} --vaneless {a.vaneless} "
                f"--np \"$NPROC\" --mpi-args \"\" --time-limit {su2_limit_s} --output-files RESTART "
                f"--out \"$GTU_OUT/run\" --record \"$GTU_OUT/su2_line.json\" || rc=$?\n"
                "# сетки (сотни МБ) не выгружаются: их заново строит frame3_geom.py по коммиту\n"
                "rm -f \"$GTU_OUT/run/mesh.su2\" \"$GTU_OUT/run/mesh_split.su2\" \"$GTU_OUT\"/geom/*/passage_custom_*.su2\n"
                "exit \"$rc\""
            )
            jobs.append({"name": name, "container": "ubuntu:24.04", "container_setup": list(CONTAINER_SETUP),
                         "command": cmd, "max_min": max_min})
        hours = max_min / 60.0
        est.append({"stage": stage, "level": a.level, "nodes": nodes, "jobs": len(speeds),
                    "t_iter_s_estimate": t_iter, "iterations_per_line": iters_line, "compute_min_per_line": compute_min,
                    "su2_time_limit_s": su2_limit_s, "max_min": max_min,
                    "expected_hours_per_line": (a.setup_min + a.mesh_min + compute_min + a.upload_min) / 60,
                    "cap_hours_per_line": hours})
    doc = {"defaults": {"instance_type": a.instance_type, "fallback_types": ["c7i.16xlarge"], "image": "al2023",
                        "disk_gb": a.disk_gb},
           "jobs": jobs}
    write_json(a.out, doc)
    exp_h = sum(e["expected_hours_per_line"] * e["jobs"] for e in est)
    cap_h = sum(e["cap_hours_per_line"] * e["jobs"] for e in est)
    disk_ipv4_h = a.disk_gb * 0.08 / 730 + 0.005
    runner_bound = sum((a.price * 1.2 + disk_ipv4_h) * (e["max_min"] + 10) / 60 * e["jobs"] for e in est)
    summary = {
        "jobs_file": rel(a.out), "n_jobs": len(jobs), "speeds": speeds, "pfracs": pfracs,
        "iterations": {"startup": startup_iter, "first_point": a.first_iter, "next_point": a.next_iter,
                       "per_line": iters_line},
        "t_ref_s_per_iter": t_ref, "t_ref_source": t_src,
        "t_ref_by_kind": {k: (sum(s[f"sec_per_iter_{k}"] * s["rows"] for s in run0["stages"] if s.get(f"sec_per_iter_{k}") and s.get("rows"))
                              / max(1, sum(s["rows"] for s in run0["stages"] if s.get(f"sec_per_iter_{k}") and s.get("rows"))))
                          for k in ("median", "p10")},
        "n_ref_nodes": n_ref, "np_ref": np_ref, "procs_aws": a.procs, "parallel_eff": a.parallel_eff,
        "parallel_eff_note": "EXAMPLE: принято (как batch/cost_estimate.toml), не замер — ускорение 4 → 64 процесса не мерилось",
        "setup_min": a.setup_min, "mesh_min": a.mesh_min, "upload_min": a.upload_min, "time_margin": a.time_margin,
        "max_min_formula": "setup_min + mesh_min + su2_time_limit/60 + 8 (выгрузка runner.py) + slack_min",
        "slack_min": a.slack_min,
        "assumption_notes": {"setup_min": "оценка: apt + pip + сборка SU2 install_su2.sh на 64 ядрах (на 4 ядрах здесь "
                                          "сборка шла ~8–10 мин, см. software)",
                             "mesh_min": "оценка: frame3_geom.py medium одной ступени с STEP (здесь 115–181 с на 4 ядрах)"},
        "price_usd_h": a.price, "price_note": "спот c7a.16xlarge ~0,71 долл./ч — допущение задания (цена 1c 14.11.2026 — 0,706)",
        "per_stage": est,
        "machine_hours_expected": exp_h, "usd_expected": exp_h * a.price,
        "machine_hours_cap": cap_h, "usd_cap_at_price": cap_h * a.price,
        "usd_runner_budget_bound": runner_bound,
        "runner_bound_note": "верхняя граница предохранителя runner.py: (MaxPrice = цена × 1,2 + диск + IPv4) × (max_min + "
                             "grace 10 мин) по всем задачам",
        "ceiling_usd": a.ceiling, "within_ceiling": bool(exp_h * a.price <= a.ceiling and runner_bound <= a.ceiling),
        "launch": "не запускалось: AWS здесь не вызывается; проверка файла — runner.py plan (без AWS)",
    }
    if not summary["within_ceiling"]:
        print(f"ВНИМАНИЕ: оценка выше потолка {a.ceiling} долл.", file=sys.stderr)
    tr["aws_plan"] = summary
    write_json(a.trial_results, tr)
    print(json.dumps({k: summary[k] for k in ("n_jobs", "machine_hours_expected", "usd_expected", "machine_hours_cap",
                                              "usd_runner_budget_bound", "within_ceiling")}, ensure_ascii=False, default=float))
    return 0


# ----------------------------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="расчёт")
    r.add_argument("--point", required=True, help="расчётная точка: путь.json:ключ (например frame3_point_results.json:stage1)")
    r.add_argument("--geom-results", default=None, help="сводка узла геометрии (frame3_geom_results.json)")
    r.add_argument("--stage", default=None, help="ступень в сводке геометрии (stage1, stage2)")
    r.add_argument("--level", default="coarse", help="уровень сетки в сводке (coarse, medium)")
    r.add_argument("--mesh", default=None, help="сетка .su2 или .su2.gz (вместо --geom-results)")
    r.add_argument("--mesh-json", default=None, help="JSON сетки impeller.py (по умолчанию impeller_*.json рядом)")
    r.add_argument("--no-sha", action="store_true", help="не сверять SHA-256 сетки со сводкой геометрии")
    r.add_argument("--force", action="store_true", help="считать, даже если точка и сетка не согласованы")
    r.add_argument("--speed", type=float, default=1.0, help="доля расчётной частоты линии")
    r.add_argument("--points", default="1.0:1000", help="ДОЛЯ:ИТЕРАЦИЙ,... — доли оценки давления на выходе")
    r.add_argument("--point-pramp", type=int, default=100, help="рампа давления к следующей точке, итераций")
    r.add_argument("--iter-scale", type=float, default=1.0, help="множитель итераций и рамп разгона")
    r.add_argument("--p-out-design", type=float, default=None, help="своё давление на выходе при расчётной частоте, Па "
                   "(вместо оценки по точке; при других частотах — в той же пропорции к оценке)")
    r.add_argument("--shroud", default="slip", choices=["stationary", "rotating", "slip"],
                   help="покрышка: slip (без трения, по умолчанию — так прошла дымовая проверка), stationary, rotating")
    r.add_argument("--vaneless", default="slip", choices=["as-is", "stationary", "slip"],
                   help="стенки за колесом r ≥ r2 (отдельные маркеры в копии сетки)")
    r.add_argument("--set", action="append", default=[], metavar="ИМЯ=ЗНАЧЕНИЕ", help=f"численный параметр: {', '.join(NUMERICS)}")
    r.add_argument("--extra", action="append", default=[], metavar="'KEY= VAL'", help="строка SU2 вне набора")
    r.add_argument("--output-files", default=None, help="OUTPUT_FILES, например RESTART (без полей ParaView)")
    r.add_argument("--np", type=int, default=4, help="процессов MPI (1 — без mpirun)")
    r.add_argument("--mpi-args", default="--mca mpi_yield_when_idle 1", help="ключи mpirun (общая машина — уступать ядро)")
    r.add_argument("--time-limit", type=float, default=1200.0, help="общий лимит счёта, с")
    r.add_argument("--su2-bin", default=None)
    r.add_argument("--su2-src", default=None, help="исходник SU2 (сверка ключей с CConfig.cpp, коммит)")
    r.add_argument("--out", default=None, help="папка прогона")
    r.add_argument("--restart-from", default=None, help="продолжить с restart_*.dat (без стадий разгона)")
    r.add_argument("--restart-state", default=None, metavar="ДОЛЯ,ДАВЛЕНИЕ", help="частота и давление рестарта")
    r.add_argument("--dry", action="store_true", help="только конфигурации и SU2_CFD -d")
    r.add_argument("--record", default=None, metavar="RESULTS_JSON", help="дописать итог в results-JSON")
    r.add_argument("--conv-orders", type=float, default=CONV["orders"])
    r.add_argument("--mass-tol", type=float, default=CONV["mass_tol"])
    r.add_argument("--window", type=int, default=CONV["window"])
    r.add_argument("--stab-tol", type=float, default=CONV["stab_tol"])
    p = sub.add_parser("parse", help="разобрать папку прогона")
    p.add_argument("run_dir")
    p.add_argument("--record", default=None, metavar="RESULTS_JSON", help="дописать итог в results-JSON")
    w = sub.add_parser("aws-plan", help="пакет задач runner.py и оценка стоимости (AWS не вызывается)")
    w.add_argument("--trial-results", default=os.path.join(HERE, "su2_trial_results.json"))
    w.add_argument("--geom-results", default=os.path.join(REPO, "routes", "gtu", "frame3", "frame3_geom_results.json"))
    w.add_argument("--out", default=os.path.join(HERE, "jobs_frame3.json"))
    w.add_argument("--stages", default="stage1,stage2")
    w.add_argument("--level", default="medium")
    w.add_argument("--speeds", default="0.7,0.8,0.9,1.0,1.05")
    w.add_argument("--pfracs", default="0.55,0.62,0.68,0.73,0.78,0.83,0.88,0.93")
    w.add_argument("--first-iter", type=int, default=3000)
    w.add_argument("--next-iter", type=int, default=2000)
    w.add_argument("--iter-scale", type=float, default=2.0, help="множитель итераций разгона на облаке")
    w.add_argument("--shroud", default="slip")
    w.add_argument("--vaneless", default="slip")
    w.add_argument("--t-ref-kind", default="median", choices=["median", "p10"],
                   help="какое время итерации пробного прогона брать: медиану (с чужой нагрузкой, с запасом) или нижний дециль")
    w.add_argument("--t-ref", type=float, default=None, help="своё время итерации эталона, с (вместо замера)")
    w.add_argument("--procs", type=int, default=64)
    w.add_argument("--parallel-eff", type=float, default=0.6)
    w.add_argument("--setup-min", type=float, default=15.0)
    w.add_argument("--mesh-min", type=float, default=10.0)
    w.add_argument("--upload-min", type=float, default=3.0, help="ожидаемое время выгрузки, мин (в смете)")
    w.add_argument("--slack-min", type=float, default=5.0, help="запас срока max_min сверх оценки, мин")
    w.add_argument("--time-margin", type=float, default=1.5, help="запас лимита счёта к оценке времени")
    w.add_argument("--price", type=float, default=0.71)
    w.add_argument("--instance-type", default="c7a.16xlarge")
    w.add_argument("--disk-gb", type=int, default=60)
    w.add_argument("--ceiling", type=float, default=300.0)
    return ap


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    if a.cmd == "run":
        if a.restart_state and not a.restart_from:
            raise SystemExit("--restart-state имеет смысл только с --restart-from")
        return run(a)
    if a.cmd == "parse":
        return parse_cmd(a)
    return aws_plan(a)


if __name__ == "__main__":
    sys.exit(main())
