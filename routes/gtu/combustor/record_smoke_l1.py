"""Воронка камеры сгорания KV, уровень 1: запись местной обкатки шаблона кейса -> comb_l1_smoke_results.json.

Обкатка - кейсы comb_cfd.py, посчитанные apt-OpenFOAM 1912 на 4 общих ядрах этой среды (не v2512 и не c7a). Для
каждого кейса: итоги comb_post.analyse, время ячейки·итерации из журналов reactingFoam (между 10-й и последней
итерацией: ExecutionTime - процессорное время, нагрузку машины учитывает частично), диагностика сходимости - расход через
сечения по оси (внутренний поток phi) и диапазон локального шага псевдовремени (поле rDeltaT).
Числа обкатки - проверка шаблона и основа сметы, не результат по камере.
Запуск: python record_smoke_l1.py --inputs <l1_inputs.json> --case <папка>:<строка>:<нагрузка>:<метка> [...]
"""
import argparse, json, os, re, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import comb_post as CP          # noqa: E402
import make_jobs_l123 as MJ     # noqa: E402

OUT = os.path.join(HERE, "comb_l1_smoke_results.json")


def flux_profile(case, time=None, k=12):
    cs = CP.Case(case, time)
    m = cs.mesh
    phi = cs.field("phi")["internal"]
    nin = len(m.neigh)
    Sf, Cf = m.Sf[:nin], m.Cf[:nin]
    zf = np.abs(Sf[:, 2]) > 0.9 * np.linalg.norm(Sf, axis=1)
    zs = sorted(set(np.round(Cf[zf, 2], 7)))
    pick = zs[::max(1, len(zs) // k)] + [zs[-1]]
    prof = []
    for z in pick:
        sel = zf & (np.abs(Cf[:, 2] - z) < 1e-6)
        prof.append(dict(z=float(z), m_wedge=float((phi[sel] * np.sign(Sf[sel, 2])).sum())))
    rd = cs.field("rDeltaT")
    dt = None
    if rd is not None and not np.isscalar(rd["internal"]):
        v = 1.0 / rd["internal"]
        dt = dict(min_s=float(v.min()), max_s=float(v.max()), median_s=float(np.median(v)))
    return dict(time=cs.time, axial_mass_flux=prof, local_dt=dt)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--case", action="append", required=True, help="папка:строка входов:нагрузка:метка")
    a = ap.parse_args(argv)
    ents = json.load(open(a.inputs, encoding="utf-8"))
    rows = []
    for spec in a.case:
        path, row, li, tag = spec.split(":")
        e = ents["entries"][int(row)]
        r = CP.analyse(path, dict(T3=e["loads"][int(li)]["T3_K"], lhv=e["lhv"]))
        if "outlet" in r and "profile" in r["outlet"]:
            r["outlet"].pop("profile")
        desc = json.load(open(os.path.join(path, "case.json"), encoding="utf-8"))
        rows.append(dict(tag=tag, idx=e["idx"], scheme=e["concept"]["scheme"], chamber=e["concept"]["chamber"],
                         kind=desc["kind"], load=desc["load"], n_cells=desc["n_cells"], fuel_mode=desc.get("fuel_mode"),
                         analysis=r, diagnosis=flux_profile(path) if desc["kind"] == "cell" else None,
                         l0=dict(NOx=e["l0"]["NOx"], CO=e["l0"]["CO"], margin_lbo=e["l0"]["margin_lbo"])))
    timing = MJ.smoke_timing(os.path.commonpath([s.split(":")[0] for s in a.case]))
    res = dict(title="Камера KV-3, уровень 1: местная обкатка шаблона осесимметричного кейса (OpenFOAM 1912)",
               method="comb_cfd.py -> blockMesh, checkMesh, reactingFoam LTS (без химии, setFields, с химией) -> "
                      "comb_post.analyse; время - журналы reactingFoam; диагностика - поток phi через сечения и rDeltaT",
               software=f"OpenFOAM {os.environ.get('WM_PROJECT_VERSION', '1912')} (apt, Ubuntu 24.04), python "
                        f"{sys.version.split()[0]}, numpy {np.__version__}",
               date="2026-11-14", inputs=[dict(file=os.path.basename(a.inputs), fingerprint=ents["fingerprint"],
                                               from_="comb_l0_debug_results.json (отбор prepare_l1, 12 вариантов)")],
               status="обкатка шаблона (не результат по камере); 4 общих ядра под нагрузкой 17-34 (чужие расчёты)",
               cases=rows, timing=timing)
    json.dump(res, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=float)
    print(json.dumps(dict(timing={k: v for k, v in timing.items() if k != "rows"},
                          cases=[(x["tag"], x["analysis"].get("status"), x["analysis"].get("mass_imbalance")) for x in rows]),
                     ensure_ascii=False, default=float))


if __name__ == "__main__":
    main()
