#!/usr/bin/env python3
"""Замер скорости CalculiX на машине: трёхмерное колесо (C3D10) — время и память статики и частот (узел bench).

Задача-замер для сметы уровней 2–3 воронки: сколько стоит один трёхмерный расчёт прочности колеса в ccx на
машине раннера. Это замер стоимости, а не расчёт прочности: колесо — твёрдотельная модель impeller.py (STEP), сетка
Gmsh из 10-узловых тетраэдров (C3D10, промежуточные узлы на прямых рёбрах — Mesh.SecondOrderLinear, чтобы грубая
сетка не давала вывернутых элементов), материал Ti-6Al-4V упругий (E 114 ГПа, ν 0,34, ρ 4430 кг/м³ — принятые,
как в platform/fea/materials.json по смыслу), закрепление — все узлы у отверстия (r ≤ r_min + 12 мм × D2/0,435),
нагрузка — вращение с номинальной частотой колеса (*DLOAD CENTRIF). Решатель — SPOOLES пакета calculix-ccx
(другого в сборке Ubuntu нет; PARDISO/PaStiX не подключены — README узла bench).

Что делает `run`:
  1. STEP колеса: impeller.py --case <case> --no-mesh (секунды);
  2. для каждой целевой сетки (--targets, узлов): шаг сетки по пробной сетке h0 (узлы ~ h⁻³), сетка Gmsh (HXT, все
     ядра), запись .inp (узлы, C3D10, ELSET SOLID);
  3. статика при числах потоков --threads (OMP_NUM_THREADS = CCX_NPROC_* = T): время, процессорное время, пиковая
     память процесса (os.wait4), max|U| из .frd — сверка, что расчёт тот же при любом T;
  4. упаковка (--pack T): nproc/T одновременных копий той же статики по T потоков — сколько расчётов в час даёт
     машина целиком (общая полоса памяти — поэтому мерить, а не делить);
  5. частоты (*FREQUENCY, --freq-modes форм) на сетках не больше --freq-max-nodes узлов, при наибольшем T.
Расчёт, который по прогнозу (время меньшей сетки × (N/N_пред)^2) не успевает до конца срока, пропускается
(status skipped_time). Срок — как в speed_su2.py (--job-max-min от /proc/uptime или --budget-s).

  python3 speed_ccx.py run --case core130 --targets 60000,200000,600000 --threads 1,8,64 --pack 8 \\
          --nproc "$NPROC" --job-max-min 75 --out "$GTU_OUT"
Выход: <out>/speed_ccx_<case>.json; рабочие файлы ccx — в --work (не выгружаются). Код 0 — хотя бы одна статика
на каждой сетке прошла; 3 — нет; 2 — ошибка входа.
"""
from __future__ import annotations

import argparse
import datetime
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PLATFORM = os.path.dirname(HERE)
IMPELLER = os.path.join(PLATFORM, "geom", "impeller.py")
sys.path.insert(0, HERE)
from speed_su2 import machine_info, uptime_s, _git_head  # noqa: E402

E_PA, NU, RHO = 114.0e9, 0.34, 4430.0


# ------------------------------------------------------------------ сетка и вход ccx
def mesh_step(step: str, h: float, threads: int, inp: str) -> dict:
    """Сетка C3D10 по STEP с наибольшим шагом h; .inp пишет Gmsh (только физическая группа SOLID)."""
    import gmsh
    t0 = time.monotonic()
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.NumThreads", threads)
        gmsh.model.occ.importShapes(step)
        gmsh.model.occ.synchronize()
        vols = [v for _, v in gmsh.model.getEntities(3)]
        if not vols:
            raise RuntimeError("в STEP нет объёмов")
        gmsh.model.addPhysicalGroup(3, vols, name="SOLID")
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", h)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)  # HXT: параллельный
        gmsh.option.setNumber("Mesh.ElementOrder", 2)
        gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)
        gmsh.option.setNumber("Mesh.HighOrderOptimize", 0)
        gmsh.option.setNumber("Mesh.SaveAll", 0)
        gmsh.model.mesh.generate(3)
        etypes, etags, _ = gmsh.model.mesh.getElements(3)
        n_el = int(sum(len(t) for t in etags))
        n_nodes = len(gmsh.model.mesh.getNodes()[0])
        gmsh.write(inp)
    finally:
        gmsh.finalize()
    return dict(h=h, nodes=n_nodes, elements=n_el, element_types=[int(t) for t in etypes],
                mesh_s=round(time.monotonic() - t0, 2))


def read_inp_nodes(inp: str):
    xs = []
    with open(inp) as f:
        it = iter(f)
        for ln in it:
            if ln.strip().upper() == "*NODE":
                break
        for ln in it:
            if ln.startswith("*"):
                break
            p = ln.split(",")
            xs.append((int(p[0]), float(p[1]), float(p[2]), float(p[3])))
    return xs


def write_decks(inp: str, rpm: float, d2: float, freq_modes: int) -> dict:
    """Вход ccx: статика (CENTRIF) и частоты; закрепление у отверстия. Ось колеса — z (impeller.py)."""
    nodes = read_inp_nodes(inp)
    r = [math.hypot(x, y) for _, x, y, _ in nodes]
    rmin = min(r)
    tol = 0.012 * d2 / 0.435
    fix = [n[0] for n, rr in zip(nodes, r) if rr <= rmin + tol]
    om = rpm * 2 * math.pi / 60
    base = os.path.basename(inp)
    common = [f"*INCLUDE,INPUT={base}", "*NSET,NSET=NFIX"]
    common += [", ".join(map(str, fix[i:i + 12])) for i in range(0, len(fix), 12)]
    common += ["*MATERIAL,NAME=TI64", "*ELASTIC", f"{E_PA:.4g},{NU}", "*DENSITY", f"{RHO:g}",
               "*SOLID SECTION,ELSET=SOLID,MATERIAL=TI64"]
    d = os.path.dirname(inp)
    static = common + ["*STEP", "*STATIC", "*BOUNDARY", "NFIX,1,3", "*DLOAD",
                       f"SOLID,CENTRIF,{om * om:.8g},0.,0.,0.,0.,0.,1.", "*NODE FILE", "U", "*END STEP"]
    freq = common + ["*BOUNDARY", "NFIX,1,3", "*STEP", "*FREQUENCY", f"{freq_modes}", "*END STEP"]
    open(os.path.join(d, "static.inp"), "w").write("\n".join(static) + "\n")
    open(os.path.join(d, "freq.inp"), "w").write("\n".join(freq) + "\n")
    return dict(fixed_nodes=len(fix), r_min_m=rmin, fix_band_m=tol, omega_rad_s=om, rpm=rpm)


# ------------------------------------------------------------------ запуск и разбор ccx
def ccx_env(threads: int) -> dict:
    t = str(threads)
    return dict(os.environ, OMP_NUM_THREADS=t, CCX_NPROC_STIFFNESS=t, CCX_NPROC_EQUATION_SOLVER=t,
                CCX_NPROC_RESULTS=t)


def run_ccx(workdir: str, job: str, threads: int, limit_s: float, ccx: str = "ccx"):
    """Один ccx: время, процессорное время и пик памяти именно этого процесса (os.wait4)."""
    t0 = time.monotonic()
    log = open(os.path.join(workdir, f"{job}.log"), "w")
    p = subprocess.Popen([ccx, "-i", job], cwd=workdir, stdout=log, stderr=subprocess.STDOUT, env=ccx_env(threads))
    return dict(p=p, t0=t0, log=log, limit=limit_s, dir=workdir, job=job, threads=threads)


def wait_ccx(h: dict) -> dict:
    p, status = h["p"], "finished"
    deadline = h["t0"] + h["limit"]
    while True:
        pid, st, ru = os.wait4(p.pid, os.WNOHANG)
        if pid:
            break
        if time.monotonic() > deadline:
            p.kill()
            pid, st, ru = os.wait4(p.pid, 0)
            status = "timeout"
            break
        time.sleep(0.5)
    p.returncode = os.waitstatus_to_exitcode(st)
    h["log"].close()
    txt = open(os.path.join(h["dir"], f"{h['job']}.log"), errors="replace").read()
    errors = len(re.findall(r"\*ERROR", txt))
    if status == "finished" and (errors or "Job finished" not in txt):
        status = "failed"
    return dict(job=h["job"], threads=h["threads"], status=status, returncode=p.returncode,
                wall_s=round(time.monotonic() - h["t0"], 2), user_s=round(ru.ru_utime, 2), sys_s=round(ru.ru_stime, 2),
                maxrss_gb=round(ru.ru_maxrss / 1024 ** 2, 3), errors=errors,
                spooles_cpus=_last_int(r"Using up to (\d+) cpu\(s\) for spooles", txt),
                ccx_total_s=_last_float(r"Total CalculiX Time:\s*([0-9.Ee+-]+)", txt))


def _last_int(pat, txt):
    m = re.findall(pat, txt)
    return int(m[-1]) if m else None


def _last_float(pat, txt):
    m = re.findall(pat, txt)
    try:
        return float(m[-1]) if m else None
    except ValueError:
        return None


def frd_max_disp(path: str):
    """max|U| по блоку DISP файла .frd (формат: i3, i10, 6e12.5)."""
    best, in_disp, n = 0.0, False, 0
    try:
        with open(path, errors="replace") as f:
            for ln in f:
                if ln.startswith(" -4") and "DISP" in ln:
                    in_disp = True
                    continue
                if in_disp:
                    if ln.startswith(" -3"):
                        break
                    if ln.startswith(" -1"):
                        try:
                            u = [float(ln[13 + 12 * i:25 + 12 * i]) for i in range(3)]
                        except ValueError:
                            continue
                        n += 1
                        best = max(best, math.sqrt(u[0] ** 2 + u[1] ** 2 + u[2] ** 2))
    except OSError:
        return None, 0
    return (best if n else None), n


def dat_frequencies(path: str, nmax: int = 10):
    """Собственные частоты (Гц) из таблицы E I G E N V A L U E   O U T P U T файла .dat."""
    out = []
    try:
        txt = open(path, errors="replace").read()
    except OSError:
        return out
    m = re.search(r"E I G E N V A L U E   O U T P U T(.*?)(?:P A R T I C I P|\Z)", txt, re.S)
    if not m:
        return out
    for ln in m.group(1).splitlines():
        p = ln.split()
        if len(p) >= 4 and p[0].isdigit():
            try:
                out.append(float(p[3]))
            except ValueError:
                pass
    return out[:nmax]


def predict(prev: dict | None, nodes: int, exp: float = 2.0):
    if not prev or not prev.get("nodes") or prev.get("status") != "finished":
        return None
    return prev["wall_s"] * (nodes / prev["nodes"]) ** exp


# ------------------------------------------------------------------ главная
def cmd_run(a) -> int:
    t_start = time.monotonic()
    try:
        targets = [int(float(x)) for x in a.targets.split(",") if x.strip()]
        threads = sorted({min(int(x), a.nproc) for x in a.threads.split(",") if x.strip()}, reverse=True)
    except ValueError as e:
        print(f"[speed_ccx] {e}", file=sys.stderr)
        return 2
    ccx = shutil.which(a.ccx)
    if not ccx:
        print(f"[speed_ccx] нет {a.ccx} в PATH", file=sys.stderr)
        return 2
    end_rel = (a.job_max_min - 8 - a.reserve_min) * 60 - uptime_s() if a.job_max_min else a.budget_s
    deadline = t_start + end_rel
    os.makedirs(a.work, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    ver = subprocess.run([ccx, "-v"], capture_output=True, text=True).stdout.strip()  # ccx -v выходит с кодом 201
    doc = dict(method=("speed_ccx.py: колесо impeller.py (STEP) -> Gmsh C3D10 (HXT, SecondOrderLinear) -> ccx статика "
                       "CENTRIF и частоты; время и пик памяти процесса — os.wait4; упаковка — nproc/T одновременных копий"),
               date=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
               machine=machine_info(a.nproc), commit=_git_head(), ccx_version=ver,
               inputs=dict(case=a.case, targets=targets, threads=threads, pack=a.pack, freq_modes=a.freq_modes,
                           freq_max_nodes=a.freq_max_nodes, job_max_min=a.job_max_min, reserve_min=a.reserve_min,
                           material=dict(E_Pa=E_PA, nu=NU, rho=RHO, note="принято; замер скорости, не прочность")),
               step={}, meshes=[], skipped=[],
               notes=["решатель SPOOLES (пакет calculix-ccx Ubuntu); PARDISO/PaStiX дали бы другое время",
                      "закрепление и материал условные: max|U| служит только сверкой одинаковости расчёта при разных T"])
    out_json = os.path.join(a.out, f"speed_ccx_{a.case}.json")

    def save():
        doc["elapsed_s"] = round(time.monotonic() - t_start, 1)
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(doc, f, ensure_ascii=False, indent=1)

    # 1. STEP
    gdir = os.path.join(a.work, "geom")
    t0 = time.monotonic()
    with open(os.path.join(a.out, f"geom_{a.case}_step.log"), "w") as fh:
        rc = subprocess.call([sys.executable, IMPELLER, "--case", a.case, "--level", "coarse", "--no-mesh", "--out", gdir,
                              "--threads", str(min(a.nproc, 8))] + a.impeller_args.split(), stdout=fh, stderr=subprocess.STDOUT)
    try:
        ij = json.load(open(os.path.join(gdir, f"impeller_{a.case}_coarse.json"), encoding="utf-8"))
        step, rpm, d2 = ij["files"]["impeller_step"], float(ij["d1"]["rpm"]), float(ij["d1"]["D2"])
    except (OSError, KeyError, ValueError, json.JSONDecodeError) as e:
        doc["step"] = dict(rc=rc, error=f"{type(e).__name__}")
        doc["ok"] = False
        save()
        return 3
    doc["step"] = dict(rc=rc, step_s=round(time.monotonic() - t0, 1), rpm=rpm, D2_m=d2,
                       volume_m3=(ij.get("solid") or {}).get("volume_m3"))
    save()

    # 2. шаг сетки по пробной сетке
    h0 = d2 / 20
    cal = os.path.join(a.work, "cal")
    os.makedirs(cal, exist_ok=True)
    try:
        m0 = mesh_step(step, h0, min(a.nproc, 8), os.path.join(cal, "mesh.inp"))
    except Exception as e:  # noqa: BLE001 — сбой Gmsh записывается в итог
        doc["calibration"] = dict(error=f"{type(e).__name__}: {e}")
        doc["ok"] = False
        save()
        return 3
    doc["calibration"] = m0
    prev = {}  # потоки -> последний итог статики (для прогноза)
    ok_meshes = 0
    for tgt in targets:
        h = h0 * (m0["nodes"] / tgt) ** (1 / 3)
        wd = os.path.join(a.work, f"n{tgt}")
        shutil.rmtree(wd, ignore_errors=True)
        os.makedirs(wd)
        if time.monotonic() > deadline - 60:
            doc["skipped"].append(dict(target=tgt, what="mesh", reason="skipped_time"))
            continue
        try:
            mi = mesh_step(step, h, a.nproc, os.path.join(wd, "mesh.inp"))
        except Exception as e:  # noqa: BLE001
            doc["meshes"].append(dict(target=tgt, h=h, error=f"{type(e).__name__}: {e}"))
            save()
            continue
        mi.update(target=tgt, decks=write_decks(os.path.join(wd, "mesh.inp"), rpm, d2, a.freq_modes), static=[],
                  pack=None, freq=None)
        doc["meshes"].append(mi)
        save()
        print(f"[speed_ccx] сетка {tgt}: узлов {mi['nodes']}, элементов {mi['elements']}, {mi['mesh_s']} с", flush=True)
        # 3. статика при разных T
        for t in threads:
            remain = deadline - time.monotonic()
            pr = predict(prev.get(t), mi["nodes"])
            if remain < 30 or (pr is not None and pr > remain):
                mi["static"].append(dict(threads=t, status="skipped_time", predicted_s=pr, remain_s=round(remain)))
                continue
            sd = os.path.join(wd, f"static_t{t}")
            os.makedirs(sd)
            shutil.copy(os.path.join(wd, "mesh.inp"), sd)
            shutil.copy(os.path.join(wd, "static.inp"), sd)
            r = wait_ccx(run_ccx(sd, "static", t, min(a.run_limit_s, remain), ccx))
            r["max_disp_m"], r["disp_nodes"] = frd_max_disp(os.path.join(sd, "static.frd"))
            r["nodes"] = mi["nodes"]
            mi["static"].append(r)
            prev[t] = r
            print(f"[speed_ccx]   статика T={t}: {r['status']} {r['wall_s']} с, {r['maxrss_gb']} ГБ, "
                  f"max|U| {r['max_disp_m']}", flush=True)
            for fn in ("static.frd", "mesh.inp"):  # большие файлы — сразу удалить (диск)
                try:
                    os.remove(os.path.join(sd, fn))
                except OSError:
                    pass
            save()
        if any(s.get("status") == "finished" for s in mi["static"]):
            ok_meshes += 1
        # 4. упаковка
        if a.pack:
            t = min(a.pack, a.nproc)
            copies = max(1, a.nproc // t)
            one = next((s for s in mi["static"] if s.get("threads") == t and s.get("status") == "finished"), None)
            ref = one or next((s for s in mi["static"] if s.get("status") == "finished"), None)
            mem_total = (doc["machine"].get("mem_total_gb") or 0) * 0.8
            remain = deadline - time.monotonic()
            est = (one or ref or {}).get("wall_s")
            if ref is None:
                mi["pack"] = dict(threads=t, copies=copies, status="skipped_no_single")
            elif mem_total and ref["maxrss_gb"] * copies > mem_total:
                mi["pack"] = dict(threads=t, copies=copies, status="skipped_memory", need_gb=ref["maxrss_gb"] * copies)
            elif est is None or est * 3 > remain:
                mi["pack"] = dict(threads=t, copies=copies, status="skipped_time", remain_s=round(remain))
            else:
                hs, t0 = [], time.monotonic()
                for i in range(copies):
                    pd = os.path.join(wd, f"pack_{i}")
                    os.makedirs(pd)
                    shutil.copy(os.path.join(wd, "mesh.inp"), pd)
                    shutil.copy(os.path.join(wd, "static.inp"), pd)
                    hs.append(run_ccx(pd, "static", t, min(a.run_limit_s, remain), ccx))
                rs = [wait_ccx(hh) for hh in hs]
                wall = time.monotonic() - t0
                fin = [x for x in rs if x["status"] == "finished"]
                mi["pack"] = dict(threads=t, copies=copies, status="finished" if len(fin) == copies else "partial",
                                  batch_wall_s=round(wall, 2), finished=len(fin),
                                  cases_per_hour_machine=round(len(fin) * 3600 / wall, 3) if fin else None,
                                  wall_s_each=[x["wall_s"] for x in rs], maxrss_gb_each=[x["maxrss_gb"] for x in rs])
                for i in range(copies):
                    shutil.rmtree(os.path.join(wd, f"pack_{i}"), ignore_errors=True)
                print(f"[speed_ccx]   упаковка {copies}x{t}: {mi['pack']['status']} {mi['pack']['batch_wall_s']} с", flush=True)
            save()
        # 5. частоты
        if mi["nodes"] <= a.freq_max_nodes:
            t = threads[0]
            remain = deadline - time.monotonic()
            st = prev.get(t)
            if remain < 60 or (st and st.get("status") == "finished" and st["wall_s"] * 4 > remain):
                mi["freq"] = dict(threads=t, status="skipped_time", remain_s=round(remain))
            else:
                fd = os.path.join(wd, "freq")
                os.makedirs(fd)
                shutil.copy(os.path.join(wd, "mesh.inp"), fd)
                shutil.copy(os.path.join(wd, "freq.inp"), fd)
                r = wait_ccx(run_ccx(fd, "freq", t, min(a.run_limit_s, remain), ccx))
                r["freq_hz"] = dat_frequencies(os.path.join(fd, "freq.dat"), a.freq_modes)
                r["modes"] = a.freq_modes
                mi["freq"] = r
                print(f"[speed_ccx]   частоты T={t}: {r['status']} {r['wall_s']} с, {r['freq_hz'][:3]} Гц", flush=True)
                shutil.rmtree(fd, ignore_errors=True)
            save()
        try:
            os.remove(os.path.join(wd, "mesh.inp"))
        except OSError:
            pass
    doc["ok"] = ok_meshes == len(targets)
    save()
    print(f"[speed_ccx] итог: {out_json}; сеток со статикой {ok_meshes} из {len(targets)}", flush=True)
    return 0 if doc["ok"] else 3


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="замер на этой машине")
    r.add_argument("--case", default="core130")
    r.add_argument("--targets", default="60000,200000,600000", help="целевые числа узлов C3D10")
    r.add_argument("--threads", default="1,8,64", help="числа потоков ccx (обрезаются до --nproc)")
    r.add_argument("--pack", type=int, default=8, help="потоков на копию при упаковке (0 — без упаковки)")
    r.add_argument("--freq-modes", type=int, default=10)
    r.add_argument("--freq-max-nodes", type=int, default=250000)
    r.add_argument("--nproc", type=int, default=int(os.environ.get("NPROC") or os.cpu_count() or 1))
    r.add_argument("--impeller-args", default="")
    r.add_argument("--ccx", default="ccx")
    r.add_argument("--job-max-min", type=int, default=0)
    r.add_argument("--reserve-min", type=int, default=6)
    r.add_argument("--budget-s", type=int, default=3600)
    r.add_argument("--run-limit-s", type=int, default=1800, help="срок одного ccx, с")
    r.add_argument("--work", default="/tmp/gtu_speed_ccx")
    r.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    return cmd_run(a)


if __name__ == "__main__":
    sys.exit(main())
