#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""routes/nh3_proton_kinetics.py — ДЕНЬ 4: кинетика протон-переноса на выживших
Li/Na-нитридах. Настоящая верификация «крадут оба» — не термодинамика (хрупкий CHE),
а БАРЬЕР: N-протонирование нитрида (нужный шаг выпуска NH₃) против H* на металле
(паразитный HER). Куда вулкан всё время указывал.

Метод (distinguished-coordinate, честный v1). Берём релаксированный мотив M₄N
низшего спина (из nh3_nitride_alloy_<compo>_results.json). Моделируем PCET как
согласованный перенос H-атома (нейтральный H — чистая парность заряда, прокси для
протон+электрон): подводим H к мишени и СКАНИРУЕМ замороженное расстояние d(H–мишень)
от 2.2 → 1.0 Å, релаксируя остальное (geomeTRIC freeze). Барьер = max(E) − E(2.2 Å).
Два канала: (N) мишень = азот нитрида → путь к NH₃; (M) мишень = вершина металла →
H* (зародыш HER). Селективность Δ‡ = барьер_M − барьер_N; >0 ⇒ N-протонирование
кинетически выигрывает ⇒ сплав отдаёт NH₃, а не жжёт водород.

Составы: Li2Na2, Li3Na (выжившие после гибрида) + Li (чистый эталон).

ОГРАНИЧЕНИЯ: нейтральный HAT-прокси — не полный электрохимический протон+электрон
с потенциалом/сольватацией; distinguished-coordinate ≠ истинное седло (верхняя
оценка барьера); каркас заморожен вне сканируемого атома; кластер, газовая фаза,
UKS-PBE/def2-SVP+DF. Меряет кинетическое ПРЕДПОЧТЕНИЕ N-vs-металл, не абсолютный ток.

Запуск:  COMPO=Li2Na2 python3 -u routes/nh3_proton_kinetics.py
Smoke:   --smoke (Li2Na2, канал N, одна точка d=1.4 — проверка пайплайна)
Выход:   routes/nh3_proton_kinetics_<compo>_results.json
"""
import os, sys, json, time, argparse, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
from pyscf import gto, dft, lib  # noqa: E402
from pyscf.geomopt.geometric_solver import optimize  # noqa: E402
from pyscf.data.elements import charge as Z  # noqa: E402

EV = 27.211386245988
DGRID = [2.2, 1.9, 1.7, 1.5, 1.3, 1.15, 1.0]     # d(H–мишень), Å
XC = os.environ.get("XC", "pbe")

def say(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
def outpath(c): return os.path.join(HERE, f"nh3_proton_kinetics_{c.lower()}_results.json")

def atomic_save(o, p):
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w") as f: json.dump(o, f, ensure_ascii=False, indent=1)
    os.replace(tmp, p)

def load_m4n(compo):
    # формат alloy: species.M4N.{xyz,spin}; формат ladder: M4N.lowest.spin -> M4N[spin].xyz
    for src in (os.path.join(HERE, f"nh3_nitride_alloy_{compo.lower()}_results.json"),
                os.path.join(HERE, f"nh3_nitride_{compo.lower()}_results.json")):
        if not os.path.exists(src): continue
        d = json.load(open(src))
        m4n = d.get("species", {}).get("M4N", {})
        if m4n.get("xyz"):
            return [[el, tuple(map(float, xyz))] for el, xyz in m4n["xyz"]], int(m4n["spin"])
        lo = d.get("M4N", {}).get("lowest", {})
        if lo.get("spin") and d.get("M4N", {}).get(lo["spin"], {}).get("xyz"):
            xyz = d["M4N"][lo["spin"]]["xyz"]
            return [[el, tuple(map(float, c)) ] for el, c in xyz], int(lo["spin"].split("=")[1])
    return None, None

def mkmf(mol):
    mf = dft.UKS(mol).density_fit()
    mf.xc = XC; mf.conv_tol = 1e-9; mf.max_cycle = 200; mf.level_shift = 0.2
    return mf

def scf(atoms, spin):
    mol = gto.M(atom=atoms, basis="def2-svp", spin=spin, charge=0, verbose=0)
    mf = mkmf(mol); e = mf.kernel()
    if not mf.converged:
        mfn = mf.newton(); mfn.max_cycle = 80
        e = mfn.kernel(mf.mo_coeff, mf.mo_occ); mf.converged = bool(mfn.converged); e = float(mfn.e_tot)
    return float(e), bool(mf.converged)

def constrained_relax(atoms, spin, i, j, d, maxsteps=30):
    """Заморозить расстояние атомов i,j = d (Å), релаксировать остальное."""
    a = [[el, tuple(xyz)] for el, xyz in atoms]
    # выставить мишень-H на нужную длину вдоль текущего направления
    p_i = np.array(a[i][1]); p_j = np.array(a[j][1])
    v = p_j - p_i; v = v / (np.linalg.norm(v) + 1e-9)
    a[j] = [a[j][0], tuple(p_i + v * d)]
    mol = gto.M(atom=a, basis="def2-svp", spin=spin, charge=0, verbose=0)
    mf = mkmf(mol)
    fd, cf = tempfile.mkstemp(dir=HERE, suffix=".txt")
    with os.fdopen(fd, "w") as f:
        f.write(f"$freeze\ndistance {i+1} {j+1}\n")     # geomeTRIC 1-индексация
    try:
        meq = optimize(mf, constraints=cf, maxsteps=maxsteps)
        e, conv = scf([[meq.atom_symbol(k),
                        tuple(float(x) for x in meq.atom_coord(k, unit="Angstrom"))]
                       for k in range(meq.natm)], spin)
        return e, conv
    except Exception as exc:
        say(f"      relax d={d} FAILED ({type(exc).__name__}: {str(exc)[:50]}) -> SP")
        return scf(a, spin)
    finally:
        try: os.remove(cf)
        except Exception: pass

def channel(atoms0, spin0, target_idx, label, store, save, grid):
    """Скан d(H–target); H = добавленный последний атом."""
    coords = np.array([x[1] for x in atoms0], float)
    t = coords[target_idx]; out = t - coords.mean(0); out = out / (np.linalg.norm(out) + 1e-9)
    base = atoms0 + [["H", tuple(t + out * 2.2)]]
    h_idx = len(base) - 1
    sp = (spin0 + 1) if (spin0 == 0) else (spin0 - 1)   # +H меняет чётность; берём разумный
    blk = store.setdefault(label, {"points": {}, "target_idx": target_idx})
    for d in grid:
        dk = f"{d:.2f}"
        if blk["points"].get(dk, {}).get("e") is not None: continue
        t0 = time.time()
        e, conv = constrained_relax(base, sp, target_idx, h_idx, d)
        blk["points"][dk] = {"e": round(e, 6), "converged": conv, "t_s": round(time.time()-t0, 1)}
        say(f"    [{label}] d={dk}Å: E={e:.5f} conv={conv} ({blk['points'][dk]['t_s']}s)")
        save()
    pts = sorted((float(k), v["e"]) for k, v in blk["points"].items() if v.get("converged"))
    if len(pts) >= 3:
        E = np.array([p[1] for p in pts]); e0 = E[0]
        blk["barrier_eV"] = round(float((E.max() - e0) * EV), 4)
        blk["d_at_max_A"] = pts[int(np.argmax(E))][0]
    return blk

def run(compo, smoke=False):
    atoms, spin = load_m4n(compo)
    if atoms is None:
        say(f"{compo}: нет M4N-геометрии (нужен nh3_nitride_alloy_{compo.lower()}_results.json)"); return 2
    n_idx = len(atoms) - 1                       # N — последний атом M4N
    OUT = outpath(compo)
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    res.setdefault("meta", {
        "compo": compo, "xc": XC,
        "what": "distinguished-coordinate H-atom-transfer barriers on M4N: channel N "
                "(nitride protonation -> NH3) vs channel M (H* on metal -> HER); "
                "selectivity = barrier_M - barrier_N; UKS-PBE/def2-SVP+DF, frozen d(H-target)",
        "limits": "neutral HAT proxy for PCET (no potential/solvent); distinguished "
                  "coordinate != true saddle (upper bound); frozen scaffold; cluster; "
                  "measures kinetic PREFERENCE N-vs-metal, not absolute rate",
    })
    sp = res.setdefault("channels", {})
    save = lambda: atomic_save(res, OUT)
    save()
    say(f"=== proton kinetics {compo} (xc={XC}) === threads={lib.num_threads()}  M4N spin 2S={spin}")
    grid = [1.4] if smoke else DGRID
    channel(atoms, spin, n_idx, "N", sp, save, grid)          # нитрид-N
    if not smoke:
        channel(atoms, spin, 0, "M", sp, save, DGRID)         # вершина металла (idx 0)
        bN = sp.get("N", {}).get("barrier_eV"); bM = sp.get("M", {}).get("barrier_eV")
        if bN is not None and bM is not None:
            res["selectivity_eV"] = round(bM - bN, 4)
            res["N_protonation_wins"] = bool(bM > bN)
        res["status"] = "ok" if res.get("selectivity_eV") is not None else "incomplete"
        save()
        say(f"  {compo}: барьер_N={bN} барьер_M={bM} Δ‡(M−N)={res.get('selectivity_eV')} "
            f"N-выигрывает={res.get('N_protonation_wins')}")
    else:
        res["status"] = "smoke"; save()
        p = sp.get("N", {}).get("points", {}).get("1.40", {})
        say(f"SMOKE {compo}: N d=1.4 E={p.get('e')} conv={p.get('converged')}")
    say(f"status={res['status']} -> {OUT}")
    return 0

def main(a):
    if a.smoke: return run("Li2Na2", smoke=True)
    c = os.environ.get("COMPO")
    for cc in ([c] if c else ["Li2Na2", "Li3Na", "Li"]):
        run(cc, smoke=False)
    return 0

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--smoke", action="store_true")
    sys.exit(main(ap.parse_args()))
