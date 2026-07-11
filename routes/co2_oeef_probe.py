#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""routes/co2_oeef_probe.py — FLIP-эксперимент №1: ориентированное внешнее
электрическое поле (OEEF, Shaik-катализ) на барьер сочленения C–C (Cu₂+2*CO).

Мотивация из НАШИХ же чисел: полное 44q-пространство показало, что корреляция
узла слабая и насыщенная (n_u≈0.49, поправка к барьеру не растёт) → рычаг для
co2-сочленения НЕ квантовый, а ЭЛЕКТРОСТАТИЧЕСКИЙ. OEEF — «мало пробованный
индустриально» переворот (лит. якорь: Shaik et al., ориентированные поля как
катализатор/переключатель селективности; на электроде поле — родная ручка:
двойной слой ~0.1-1 В/Å). Вопрос пробы: НАСКОЛЬКО барьер сочленения управляется
полем вдоль оси C···C — dΔE‡/dF, знак и величина.

Протокол: геометрии реактанта (d=4.2 Å) и TS (d=2.6 Å) Cu₂ — 1:1 реконструкция
из ts_barrier.profile (build_ads из co2_to_fuels_ts.py, как в casscf-скриптах).
DF-RKS PBE/def2-SVP + однородное поле E вдоль x (ось C···C): H = H₀ + E·r
(finite-field через подмену get_hcore, common origin в центре масс). Поля:
0, ±0.002, ±0.005, ±0.010 а.е. (0.01 а.е. = 0.514 В/Å — масштаб двойного слоя).
Warm-chain по полю (dm предыдущей точки) для устойчивости ветки SCF.

Ограничения: кластер (не электрод), однородное поле (не двойной слой), PBE (мы
знаем его занижение барьера), замороженный металл, def2-SVP; проба измеряет
ЧУВСТВИТЕЛЬНОСТЬ dΔE‡/dF на кластерном уровне, не предсказывает ток/потенциал.

Запуск: OMP_NUM_THREADS=4 python3 -u routes/co2_oeef_probe.py
Smoke:  OMP_NUM_THREADS=1 python3 routes/co2_oeef_probe.py --smoke  (2×CO, ~1 мин)
Выход:  routes/co2_oeef_results.json (atomic, incremental, рестарт по точкам)
"""
import os, sys, json, time, argparse, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pyscf import gto, dft, lib  # noqa: E402

SITE = os.environ.get("OEEF_SITE", "Cu2")
AXIS = {"x":0,"y":1,"z":2}[os.environ.get("OEEF_AXIS","x")]
OUT = os.path.join(HERE, f"co2_oeef_{SITE.lower()}{'_'+os.environ['OEEF_AXIS'] if os.environ.get('OEEF_AXIS','x')!='x' else ''}_results.json" if (SITE!="Cu2" or os.environ.get("OEEF_AXIS","x")!="x") else "co2_oeef_results.json")
HARTREE_EV = 27.211386245988
AU_FIELD_V_PER_A = 51.42206747  # 1 а.е. поля = 51.422 В/Å
FIELDS = [0.0, 0.002, -0.002, 0.005, -0.005, 0.010, -0.010]  # а.е., вдоль x

def say(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

def atomic_save(obj, path=OUT):
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1)
    os.replace(tmp, path)

def scf_in_field(mol, efield_x, dm0=None, verbose=0):
    """DF-RKS PBE с однородным полем E вдоль x: H = H0 + E*x (finite field)."""
    mf = dft.RKS(mol).density_fit()
    mf.xc = "pbe"
    mf.verbose = verbose
    mf.conv_tol = 1e-9
    mf.max_cycle = 200
    if abs(efield_x) > 0:
        mol.set_common_orig(np.einsum("z,zx->x", mol.atom_charges(),
                                      mol.atom_coords()) / mol.atom_charges().sum())
        ao_dip = mol.intor_symmetric("int1e_r", comp=3)
        h0 = mf.get_hcore()
        E = np.zeros(3); E[AXIS] = efield_x
        mf.get_hcore = lambda *a: h0 + np.einsum("x,xij->ij", E, ao_dip)
    e = mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel()
    if not mf.converged:
        mfn = mf.newton(); mfn.max_cycle = 80
        e = mfn.kernel(mf.mo_coeff, mf.mo_occ)
        mf.mo_coeff, mf.mo_occ = mfn.mo_coeff, mfn.mo_occ
        mf.e_tot, mf.converged = float(mfn.e_tot), bool(mfn.converged)
    return mf, float(mf.e_tot), bool(mf.converged)

def run_geometry(tag, mol, res, args):
    blk = res.setdefault(tag, {})
    dm = None
    # warm-chain: 0 -> +F по возрастанию, отдельно 0 -> -F (сортируем по |F|, знаку)
    order = sorted(FIELDS, key=lambda f: (abs(f), -np.sign(f)))
    for f in order:
        key = f"{f:+.3f}"
        if key in blk and blk[key].get("converged"):
            say(f"  [{tag}] F={key} а.е. — уже есть (рестарт)")
            continue
        t0 = time.time()
        mf, e, conv = scf_in_field(mol, f, dm0=dm, verbose=args.verbose)
        dm = mf.make_rdm1()
        blk[key] = {"efield_au": f, "efield_V_per_A": round(f * AU_FIELD_V_PER_A, 3),
                    "e_pbe": round(e, 6), "converged": conv,
                    "t_s": round(time.time() - t0, 1)}
        say(f"  [{tag}] F={key} а.е. ({f*AU_FIELD_V_PER_A:+.2f} В/Å): "
            f"E={e:.6f} conv={conv} ({blk[key]['t_s']}s)")
        atomic_save(res)
    return blk

def summarize(res):
    r, t = res.get("reactant", {}), res.get("ts", {})
    rows = []
    for f in FIELDS:
        key = f"{f:+.3f}"
        if key in r and key in t and r[key].get("converged") and t[key].get("converged"):
            bar = (t[key]["e_pbe"] - r[key]["e_pbe"]) * HARTREE_EV
            rows.append({"efield_au": f, "efield_V_per_A": round(f * AU_FIELD_V_PER_A, 3),
                         "barrier_eV": round(bar, 4)})
    rows.sort(key=lambda x: x["efield_au"])
    res["barrier_vs_field"] = rows
    if len(rows) >= 3:
        F = np.array([x["efield_au"] for x in rows])
        B = np.array([x["barrier_eV"] for x in rows])
        slope = float(np.polyfit(F, B, 1)[0])          # эВ / а.е. поля
        res["summary"] = {
            "dBarrier_dField_eV_per_au": round(slope, 3),
            "dBarrier_dField_eV_per_VA": round(slope / AU_FIELD_V_PER_A, 4),
            "barrier_zero_field_eV": next((x["barrier_eV"] for x in rows
                                           if x["efield_au"] == 0.0), None),
            "best_field": min(rows, key=lambda x: x["barrier_eV"]),
            "note": ("dΔE‡/dF вдоль оси C···C; знак минус = поле в +x снижает "
                     "барьер. Кластер/PBE/однородное поле — чувствительность, "
                     "не предсказание потенциала."),
        }
    res["status"] = "ok" if len(rows) == len(FIELDS) else "incomplete"
    atomic_save(res)

def main_real(args):
    from co2_to_fuels_sites import METALS, _mol
    from co2_to_fuels_ts import build_ads
    say(f"OEEF probe — threads={lib.num_threads()}")
    with open(os.path.join(HERE, "co2_to_fuels_results.json")) as f:
        tb = json.load(f)["ts_barrier"][SITE]
    prof = {round(p["d"], 2): p for p in tb["profile"]}
    geoms = {"reactant": tb["d_reactant"], "ts": tb["d_ts"]}
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    res.setdefault("meta", {
        "what": f"OEEF probe ({SITE}): uniform E-field along the C...C axis vs coupling "
                "barrier, metal+2*CO cluster, DF-RKS PBE/def2-SVP, warm-chain in field",
        "fields_au": FIELDS, "axis": "x (C...C)",
        "geometry_source": "ts_barrier.profile of co2_to_fuels_results.json "
                           "(reactant d=4.2, TS d=2.6), build_ads reconstruction",
        "honesty": "cluster != electrode; uniform field != double layer; PBE "
                   "underestimates the barrier (known); frozen metal; def2-SVP. "
                   "This measures SENSITIVITY dE^/dF, not an operating potential.",
        "flip_motivation": "44q full-CAS showed weak/saturated correlation "
                           "(n_u~0.49) -> the lever for co2 coupling is "
                           "electrostatics, not quantum correlation",
    })
    atomic_save(res)
    for tag, d in geoms.items():
        p = prof[round(float(d), 2)]
        ads = build_ads(float(d), np.array(p["free"], dtype=float))
        mol = _mol(METALS[SITE], ads)
        say(f"[{tag}] d={d} Å, nao={mol.nao}")
        run_geometry(tag, mol, res, args)
    summarize(res)
    say(f"status={res['status']} -> {OUT}")

def main_smoke(args):
    say("SMOKE: 2xCO (без металла), поле ±0.005 а.е., реактант/TS-подобные d")
    res = {"meta": {"mode": "smoke"}}
    out = os.path.join(tempfile.mkdtemp(prefix="oeef_smoke_"), "smoke.json")
    global OUT, FIELDS
    OUT = out
    FIELDS = [0.0, 0.005, -0.005]
    def mk(d):
        a = (f"C {-d/2:.3f} 0 0; O {-d/2:.3f} 0 1.135; "
             f"C {d/2:.3f} 0 0; O {d/2:.3f} 0 1.135")
        return gto.M(atom=a, basis="def2-svp", spin=0, verbose=0)
    for tag, d in [("reactant", 4.2), ("ts", 2.6)]:
        run_geometry(tag, mk(d), res, args)
    summarize(res)
    ok = res["status"] == "ok" and all(
        abs(res[t][f"{f:+.3f}"]["e_pbe"]) > 100 for t in ("reactant", "ts")
        for f in FIELDS)
    say(f"SMOKE {'PASS' if ok else 'FAIL'}; поле сдвигает энергию: "
        f"{res['reactant']['+0.005']['e_pbe'] != res['reactant']['+0.000']['e_pbe']}")
    return 0 if ok else 1

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--verbose", type=int, default=0)
    a = ap.parse_args()
    sys.exit(main_smoke(a) if a.smoke else main_real(a))
