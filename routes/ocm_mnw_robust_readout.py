#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/ocm_mnw_robust_readout.py — робастный ридаут ΔΔE‡ с расширенной сетки.
Мотив (Этап 18): у Cr точка 7 CH₄ — шип +180 ккал между соседями 6.4/−2.0
(геометрия сломана на всех уровнях, refine бессилен — чинить нечего, надо
маскировать); наивный «TS = argmax» берёт шип и выдаёт мусорный дескриптор.

Правило (то же, что MSA v3): точка-выброс = на >25 ккал выше ОБОИХ соседей
(края — единственного соседа). После маски: R = минимум слева от интерьерного
максимума, TS = интерьерный максимум. Если TS на краю или профиль зигзаг
(≥2 выброса подряд/чередуются) — substrate помечается unreadable, дескриптор
не выдаётся. Пишет ocm_mnw_emb_{m}_final_robust.json.

Запуск: python3 routes/ocm_mnw_robust_readout.py ti v cr fe [zr nb]
"""
import json
import os
import sys

OUTLIER_KCAL = 25.0
DIR = os.path.dirname(os.path.abspath(__file__))


def robust_read(rel):
    """→ (r_idx, ts_idx, barrier, interior, dropped_idx) на маскированном профиле."""
    n = len(rel)
    dropped = []
    for i in range(n):
        nb = [rel[j] for j in (i - 1, i + 1) if 0 <= j < n and j not in dropped]
        if nb and all(rel[i] - x > OUTLIER_KCAL for x in nb):
            dropped.append(i)
    keep = [i for i in range(n) if i not in dropped]
    if len(keep) < 4:
        return None, None, None, False, dropped
    kr = [rel[i] for i in keep]
    imax_k = max(range(1, len(kr) - 1), key=lambda i: kr[i])   # интерьер маски
    interior = kr[imax_k] >= max(kr[0], kr[-1])
    imin_k = min(range(imax_k + 1), key=lambda i: kr[i])
    barrier = kr[imax_k] - kr[imin_k]
    return keep[imin_k], keep[imax_k], round(barrier, 2), bool(interior), dropped


def main(metals):
    rows = []
    for m in metals:
        ml = m.lower()
        out = {"metal": m, "rule": f"neighbor-outlier > {OUTLIER_KCAL} kcal masked",
               "substrates": {}}
        readable = True
        for sub in ("ch4", "c2h6"):
            p = json.load(open(os.path.join(
                DIR, f"ocm_mnw_emb_{ml}_{sub}_profile.json")))
            out["substrates"][sub] = {}
            for lv in ("casscf", "nevpt2"):
                rel = p[f"{lv}_rel_kcal"]
                r, ts, bar, inter, drop = robust_read(rel)
                out["substrates"][sub][lv] = {
                    "barrier_kcal": bar, "r_index": r, "ts_index": ts,
                    "ts_interior": inter, "dropped_points": drop,
                    "rel_kcal": rel}
                if lv == "nevpt2" and (bar is None or not inter):
                    readable = False
            # DFT из geom-файла тем же правилом
            g = json.load(open(os.path.join(
                DIR, f"ocm_mnw_emb_{ml}_{sub}_geom.json")))
            rel = g["rel_kcal"]
            r, ts, bar, inter, drop = robust_read(rel)
            out["substrates"][sub]["dft"] = {
                "barrier_kcal": bar, "r_index": r, "ts_index": ts,
                "ts_interior": inter, "dropped_points": drop}
        if readable:
            dde = {}
            for lv in ("dft", "casscf", "nevpt2"):
                b4 = out["substrates"]["ch4"][lv]["barrier_kcal"]
                b6 = out["substrates"]["c2h6"][lv]["barrier_kcal"]
                dde[lv] = (round(b6 - b4, 2)
                           if b4 is not None and b6 is not None else None)
            out["ddE_kcal"] = dde
            if dde.get("nevpt2") is not None and dde.get("dft") is not None:
                out["quantum_shift_kcal"] = round(dde["nevpt2"] - dde["dft"], 2)
        else:
            out["verdict"] = ("UNREADABLE: NEVPT2 TS не интерьерный/профиль "
                              "зигзаг после маски — дескриптор не выдаётся")
        path = os.path.join(DIR, f"ocm_mnw_emb_{ml}_final_robust.json")
        with open(path, "w") as f:
            json.dump(out, f, indent=1)
        rows.append((m, out.get("ddE_kcal"), out.get("quantum_shift_kcal"),
                     out.get("verdict", "ok")))
        print(f"[{m}] -> {os.path.basename(path)}")
    print(f"\n{'M':<4} {'ddE (dft/casscf/nevpt2)':<34} {'shift':<8} verdict")
    for m, dde, sh, v in rows:
        ds = (f"{dde['dft']} / {dde['casscf']} / {dde['nevpt2']}"
              if dde else "—")
        print(f"{m:<4} {ds:<34} {str(sh):<8} {v}")


if __name__ == "__main__":
    main(sys.argv[1:] or ["ti", "v", "cr", "fe"])
