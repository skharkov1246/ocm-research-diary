#!/usr/bin/env python3
"""Проверки второго пробного прогона раннера v2 (jobs.smoke2.json) — на машине AWS или локально.

  python  routes/gtu/platform/aws/smoke2.py py    # Python, cantera, numpy, scipy; точки цикла gtu_cycle.cycle в пуле NPROC
  python3 routes/gtu/platform/aws/smoke2.py ccx   # ccx -v; осесимметричный диск disk.py на грубой сетке

Пишет в $GTU_OUT (без него — ./smoke2_out): smoke_<режим>.json и az_check.json — встала ли машина
в зону, которую раннер выбрал (GTU_AZ), и была ли она самой дешёвой (GTU_AZ_CHEAPEST); зона машины —
GTU_AZ_MACHINE (снаружи контейнера, из IMDS) или IMDS напрямую.

Код выхода: 0 — всё прошло; 2 — сбой счёта (импорт, cantera, ccx, disk.py); 3 — счёт прошёл, но машина
не в самой дешёвой зоне (причина — в az_check.json; при сбое ёмкости раннер законно берёт следующую зону).
"""
import json
import os
import subprocess
import sys
import time
import urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", ".."))
OUT = os.path.abspath(os.environ.get("GTU_OUT") or "smoke2_out")

# Точка цикла для сверки: SGT-400-подобные КПД узлов (как калибровка в gtu_cycle.py), T41 = 1150 °C.
CYCLE_ARGS = dict(PR=16.8, T41_C=1150.0, eta_pc=0.895, eta_pt_gg=0.875, eta_pt_pt=0.895,
                  x_ngv=0.055, x_rot=0.055, x_leak=0.015, W1=39.4 / 1.006)
# Локальный расчёт той же точки 14.11.2026 (Python 3.13, cantera 3.2.0, numpy 2.5.3, scipy 1.18.1):
# другая версия cantera может дать малое расхождение — оно записывается, но счёт не валит.
CYCLE_REF = {"P_MW": 12.1664, "eff": 0.3623, "source": "локально 14.11.2026, cantera 3.2.0"}


def imds(path: str):
    """IMDSv2; None, если недоступен (из контейнера при hop limit 1, локально)."""
    direct = urllib.request.build_opener(urllib.request.ProxyHandler({}))   # IMDS — только напрямую
    try:
        req = urllib.request.Request("http://169.254.169.254/latest/api/token", method="PUT",
                                     headers={"X-aws-ec2-metadata-token-ttl-seconds": "60"})
        tok = direct.open(req, timeout=2).read().decode()
        req = urllib.request.Request(f"http://169.254.169.254/latest/meta-data/{path}",
                                     headers={"X-aws-ec2-metadata-token": tok})
        return direct.open(req, timeout=2).read().decode().strip()
    except Exception:  # noqa: BLE001
        return None


def nproc() -> int:
    try:
        return max(1, int(os.environ.get("NPROC", "")))
    except ValueError:
        return os.cpu_count() or 1


def az_check() -> dict:
    planned = os.environ.get("GTU_AZ", "")
    cheapest = os.environ.get("GTU_AZ_CHEAPEST", "")
    machine = os.environ.get("GTU_AZ_MACHINE", "") or imds("placement/availability-zone") or ""
    d = {"planned": planned, "cheapest": cheapest, "machine": machine,
         "machine_in_planned": bool(machine) and machine == planned,
         "planned_is_cheapest": bool(planned) and planned == cheapest}
    d["ok"] = d["machine_in_planned"] and d["planned_is_cheapest"]
    if not planned:
        d["why"] = "GTU_AZ пуст: раннер не знал цен по зонам (или запуск не раннером)"
    elif not d["planned_is_cheapest"]:
        d["why"] = "раннер взял не самую дешёвую зону: сбой ёмкости в ней или зона после прерывания — см. ledger"
    elif not d["machine_in_planned"]:
        d["why"] = "зона машины не совпала с запрошенной или не прочитана"
    return d


def _cycle_point(t41_c: float) -> dict:
    sys.path.insert(0, os.path.join(ROOT, "routes", "gtu"))
    import gtu_cycle  # noqa: E402  (меняет cwd на routes/gtu — пути здесь абсолютные)
    args = dict(CYCLE_ARGS, T41_C=t41_c)
    r = gtu_cycle.cycle(**args)
    return {k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()}


def check_py() -> tuple:
    out = {"python": sys.version.split()[0], "executable": sys.executable, "nproc": nproc()}
    import cantera
    import numpy
    import scipy
    out.update(cantera=cantera.__version__, numpy=numpy.__version__, scipy=scipy.__version__)
    try:
        import boto3
        out["boto3"] = boto3.__version__
    except ImportError:
        out["boto3"] = None
    t = time.time()
    point = _cycle_point(CYCLE_ARGS["T41_C"])
    out["cycle_point"] = point
    out["cycle_point_s"] = round(time.time() - t, 2)
    out["cycle_ref"] = CYCLE_REF
    out["cycle_ref_rel_diff"] = {k: round(abs(point[k] / CYCLE_REF[k] - 1), 6) for k in ("P_MW", "eff")}
    # Пул процессов на NPROC ядер (поле Д раннера): по точке на ядро, T41 от 1000 до 1250 °C.
    from concurrent.futures import ProcessPoolExecutor
    n = min(nproc(), 64)
    temps = [1000.0 + 250.0 * i / max(1, n - 1) for i in range(n)]
    t = time.time()
    with ProcessPoolExecutor(max_workers=n) as ex:
        pts = list(ex.map(_cycle_point, temps))
    out["pool"] = {"workers": n, "points": len(pts), "wall_s": round(time.time() - t, 2),
                   "P_MW_range": [round(min(p["P_MW"] for p in pts), 4), round(max(p["P_MW"] for p in pts), 4)]}
    return out, 0


def check_ccx() -> tuple:
    out = {"python": sys.version.split()[0], "nproc": nproc()}
    p = subprocess.run(["ccx", "-v"], capture_output=True, text=True, timeout=60)
    ver = next((ln.strip() for ln in (p.stdout + p.stderr).splitlines() if "Version" in ln), "")
    out["ccx_version"] = ver
    out["ccx_v_exit_code"] = p.returncode          # ccx -v выходит с кодом 201 — норма (Dockerfile платформы)
    if not ver:
        out["error"] = "ccx -v не напечатал версию"
        return out, 2
    sys.path.insert(0, os.path.join(ROOT, "routes", "gtu", "platform", "fea"))
    import disk  # noqa: E402
    case = dict(disk.example_cases()[0], mesh=dict(nr=30, nz=4))     # колесо компрессора схемы 130, грубо
    t = time.time()
    res = disk.run_case(case, workdir=os.path.join(OUT, "disk_run"))
    out["disk"] = {"case": case["name"], "mesh": case["mesh"], "wall_s": round(time.time() - t, 2),
                   "vm_max_node_MPa": round(res["vm_max_node_MPa"], 3),
                   "burst_margin": round(res["burst"]["margin"], 4)}
    return out, 0


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1 or argv[0] not in ("py", "ccx"):
        print(__doc__)
        return 2
    mode = argv[0]
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    try:
        res, rc = (check_py if mode == "py" else check_ccx)()
    except Exception as e:  # noqa: BLE001
        res, rc = {"error": f"{type(e).__name__}: {e}"}, 2
    res["wall_s"] = round(time.time() - t0, 2)
    az = az_check()
    with open(os.path.join(OUT, f"smoke_{mode}.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1)
    with open(os.path.join(OUT, "az_check.json"), "w", encoding="utf-8") as f:
        json.dump(az, f, ensure_ascii=False, indent=1)
    print(json.dumps({"mode": mode, "rc": rc, "az_ok": az["ok"], **{k: res.get(k) for k in
                      ("python", "cantera", "ccx_version", "wall_s", "error") if k in res}}, ensure_ascii=False))
    if rc == 0 and not az["ok"]:
        rc = 3
    return rc


if __name__ == "__main__":
    sys.exit(main())
