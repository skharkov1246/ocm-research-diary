#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""routes/cc_cogef_cas22.py — бирадикальная поправка к COGEF: CASSCF(2,2)
single-point'ы на сохранённых геометриях скана cc_cogef_results.json.

Зачем: в COGEF-скане ⟨S²⟩=0 на всей сетке — UKS-PBE схлопнулся в рестриктед-
подобное решение, а гомолиз C–C это σ/σ* бирадикал: PBE на растяжении завышает
энергию хвоста (или занижает — вопрос именно к этой пробе). CASSCF(2,2) в
σ(C–C)/σ*(C–C) окне — минимальная объективная мультиреференс-обработка разрыва.
Считаем E_cas(d) и NOON → n_u(d): где включается бирадикал и как меняется
ΔE‡(F) против чистого PBE.

Метод: для каждой точки скана (xyz сохранены в JSON) RHF/def2-SVP+DF →
AVAS не нужен: стартовые σ/σ* берём HOMO/LUMO RHF (для растянутой центральной
C–C это именно они; проверка — NOON: если 2.0/0.0 на равновесии и дробные на
хвосте, окно выбрано верно) → CASSCF(2,2), fix_spin ss=0. Без градиентов
(упавшая точка 2.80 Å не нужна — у неё нет xyz).

Ограничения: (2,2)-окно минимально (без гиперконъюгации/CH-орбиталей); геометрии —
PBE-релаксированные (вертикальная поправка); бутан-модель как в COGEF.

Запуск: OMP_NUM_THREADS=4 python3 -u routes/cc_cogef_cas22.py
Smoke:  OMP_NUM_THREADS=1 python3 routes/cc_cogef_cas22.py --smoke (~1 мин)
Выход:  routes/cc_cogef_cas22_results.json
"""
import os, sys, json, time, argparse, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
from pyscf import gto, scf, mcscf, fci, lib  # noqa: E402

SRC = os.path.join(HERE, "cc_cogef_results.json")
OUT = os.path.join(HERE, "cc_cogef_cas22_results.json")
EV = 27.211386245988
NN_PER_HA_BOHR = 82.38723498
A2B = 1.8897259886
FGRID_NN = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

def say(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

def atomic_save(obj):
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1)
    os.replace(tmp, OUT)

def cas22_at(atoms, verbose=0):
    mol = gto.M(atom=atoms, basis="def2-svp", spin=0, verbose=0)
    mf = scf.RHF(mol).density_fit()
    mf.conv_tol = 1e-9; mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        mfn = mf.newton(); mfn.max_cycle = 60
        mfn.kernel(mf.mo_coeff, mf.mo_occ)
        mf.mo_coeff, mf.mo_occ = mfn.mo_coeff, mfn.mo_occ
        mf.converged = bool(mfn.converged); mf.e_tot = float(mfn.e_tot)
    mc = mcscf.CASSCF(mf, 2, 2)
    mc.conv_tol = 1e-8; mc.max_cycle_macro = 60
    try: fci.addons.fix_spin_(mc.fcisolver, ss=0)
    except Exception: pass
    e = mc.kernel()[0]
    dm1 = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    noon = sorted(np.linalg.eigvalsh(dm1), reverse=True)
    n_u = float(sum(min(x, 2 - x) for x in noon))
    return float(mf.e_tot), float(e), bool(mc.converged), \
        [round(float(x), 4) for x in noon], round(n_u, 4)

def analyze(res):
    pts = sorted((float(k), v) for k, v in res["points"].items()
                 if v.get("cas_converged"))
    if len(pts) < 4: return
    D = np.array([p[0] for p in pts])
    E = np.array([p[1]["e_cas"] for p in pts]); E -= E.min()
    dEdd = np.gradient(E, D)
    rows = []
    for F in FGRID_NN:
        f_au = F / NN_PER_HA_BOHR * A2B
        Et = E - f_au * (D - D[0])
        i_ts = int(np.argmax(Et)); i_min = int(np.argmin(Et[:i_ts + 1])) if i_ts else 0
        rows.append({"F_nN": F, "barrier_eV": round(float((Et[i_ts] - Et[i_min]) * EV), 3)})
    res["cogef_cas22"] = {
        "F_max_nN": round(float(dEdd.max() * NN_PER_HA_BOHR / A2B), 2),
        "barrier_vs_force": rows,
        "onset_biradical_A": next((f"{p[0]:.2f}" for p in pts if p[1]["n_u"] > 0.5), "за сеткой"),
    }

def main(args):
    say(f"CASSCF(2,2)//COGEF-геометрии — threads={lib.num_threads()}")
    src = json.load(open(SRC))
    scan = {k: v for k, v in src["scan"].items() if v.get("converged") and v.get("xyz")}
    keys = sorted(scan, key=float)
    if args.smoke: keys = [keys[0], keys[-1]]
    res = json.load(open(OUT)) if os.path.exists(OUT) and not args.smoke else {}
    res.setdefault("meta", {
        "what": "CASSCF(2,2) sigma/sigma* single-points on the saved COGEF scan "
                "geometries (biradical-honest homolysis tail); RHF start "
                "HOMO/LUMO window; fix_spin ss=0",
        "source": "cc_cogef_results.json (UKS-PBE relaxed constrained scan)",
        "honesty": "(2,2) minimal window; vertical on PBE geometries; butane model",
    })
    res.setdefault("points", {})
    for k in keys:
        if k in res["points"] and res["points"][k].get("cas_converged"):
            say(f"  d={k} — есть (рестарт)"); continue
        t0 = time.time()
        try:
            e_hf, e_cas, conv, noon, n_u = cas22_at(scan[k]["xyz"])
            res["points"][k] = {"d_A": float(k), "e_hf": e_hf, "e_cas": e_cas,
                                "cas_converged": conv, "noon": noon, "n_u": n_u,
                                "t_s": round(time.time() - t0, 1)}
            say(f"  d={k} Å: E_cas={e_cas:.6f} conv={conv} NOON={noon} n_u={n_u}")
        except Exception as exc:
            res["points"][k] = {"d_A": float(k), "cas_converged": False,
                                "error": f"{type(exc).__name__}: {str(exc)[:80]}"}
            say(f"  d={k} Å FAILED: {exc}")
        if not args.smoke: atomic_save(res)
    if not args.smoke:
        analyze(res)
        done = sum(1 for v in res["points"].values() if v.get("cas_converged"))
        res["status"] = "ok" if done == len(keys) else f"partial({done}/{len(keys)})"
        atomic_save(res)
        say(f"status={res['status']} -> {OUT}")
    else:
        ok = all(v.get("cas_converged") for v in res["points"].values())
        eq, st = res["points"][keys[0]], res["points"][keys[-1]]
        say(f"SMOKE {'PASS' if ok else 'FAIL'}: n_u(eq)={eq.get('n_u')} -> n_u({keys[-1]})={st.get('n_u')}")
        return 0 if ok else 1
    return 0

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    sys.exit(main(a))
