#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""routes/co2_cation_probe.py — ДОЖИМ ДОПУЩЕНИЯ №1 OEEF-результата:
«однородное поле ≈ реальный двойной слой». Вместо абстрактного поля ставим
НАСТОЯЩИЙ катион K⁺ около CuAl-сайта и меряем сдвиг барьера сочленения.

Физика моста: точечный заряд +1 на расстоянии r создаёт поле E≈14.4/r² В/Å
(6 Å → ~0.40 В/Å — ровно диапазон нашей OEEF-ручки). Литературный якорь:
катионные эффекты в CO2RR (K⁺/Cs⁺ ускоряют C2+) объясняют полем катиона у
*OCCO — наш прогноз слопа −0.384 эВ/(В/Å) должен воспроизвести это «настоящим»
ионом, если допущение однородного поля честно.

Протокол: геометрии реактанта (d=3.0) и TS (d=2.2) CuAl из ts_barrier.profile
(как в OEEF-пробе), + K⁺ в позициях вдоль оси C···C с обеих сторон на 6 и 8 Å
от середины связи C–C (z=1.5 Å, плоскость адсорбатов), заряд системы +1,
RKS-PBE/def2-SVP+DF. Сдвиг барьера vs голый сайт (0.9937 эВ) и vs предсказание
слопа по полю точечного заряда в середине C–C.
САМОТЕСТ: K⁺ на 20 Å (поле ~0.036 В/Å) обязан дать барьер ≈ голому (±0.02 эВ).

ЧЕСТНО: газофазный ион БЕЗ растворителя и противоиона — поле сверху-оценка
(реальный электролит экранирует); кластер, PBE, замороженный металл. Проба
отвечает на один вопрос: работает ли ручка от РЕАЛЬНОГО иона, а не только от
математического поля.

Запуск: OMP_NUM_THREADS=4 python3 -u routes/co2_cation_probe.py
Smoke:  --smoke (только самотест 20 Å на реактанте, ~10-20 мин)
Выход:  routes/co2_cation_probe_results.json (atomic, incremental)
"""
import os, sys, json, time, argparse, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pyscf import gto, dft, lib  # noqa: E402

OUT = os.path.join(HERE, "co2_cation_probe_results.json")
EV = 27.211386245988
BARE_BARRIER_EV = 0.9937          # наш OEEF-прогон, F=0, та же ветка что профиль
SLOPE_EV_PER_VA = -0.384          # наш чистый слоп (ось x, CuAl)
SITE = "CuAl"

def say(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

def atomic_save(obj):
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1)
    os.replace(tmp, OUT)

def geometries():
    from co2_to_fuels_sites import METALS, _mol  # noqa: F401
    from co2_to_fuels_ts import build_ads
    with open(os.path.join(HERE, "co2_to_fuels_results.json")) as f:
        tb = json.load(f)["ts_barrier"][SITE]
    prof = {round(p["d"], 2): p for p in tb["profile"]}
    out = {}
    for tag, d in (("reactant", tb["d_reactant"]), ("ts", tb["d_ts"])):
        p = prof[round(float(d), 2)]
        ads = build_ads(float(d), np.array(p["free"], dtype=float))
        atoms = [[el, tuple(map(float, xyz))] for el, xyz in METALS[SITE]]
        for el, xyz in zip(["C", "O", "C", "O"], ads):
            atoms.append([el, tuple(float(v) for v in xyz)])
        out[tag] = (atoms, float(d))
    return out

def scf_with_k(atoms, k_pos, verbose=0):
    full = atoms + [["K", tuple(k_pos)]] if k_pos is not None else atoms
    mol = gto.M(atom=full, basis="def2-svp", spin=0,
                charge=(1 if k_pos is not None else 0), verbose=0)
    mf = dft.RKS(mol).density_fit()
    mf.xc = "pbe"; mf.verbose = verbose
    mf.conv_tol = 1e-9; mf.max_cycle = 200
    e = mf.kernel()
    if not mf.converged:
        mfn = mf.newton(); mfn.max_cycle = 80
        e = mfn.kernel(mf.mo_coeff, mf.mo_occ)
        mf.converged = bool(mfn.converged); e = float(mfn.e_tot)
    return float(e), bool(mf.converged)

def field_at_cc_mid(k_pos, cc_mid):
    r = np.linalg.norm(np.array(k_pos) - np.array(cc_mid))
    ex = 14.39964 / r**2 * (cc_mid[0] - k_pos[0]) / r   # В/Å, x-компонента
    return round(float(ex), 3), round(float(r), 2)

def main(args):
    say(f"cation probe ({SITE}) — threads={lib.num_threads()}")
    g = geometries()
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    res.setdefault("meta", {
        "what": "explicit K+ near the CuAl coupling site vs the uniform-field "
                "OEEF prediction; system charge +1, RKS-PBE/def2-SVP+DF",
        "bare_barrier_eV": BARE_BARRIER_EV,
        "oeef_slope_eV_per_VA": SLOPE_EV_PER_VA,
        "honesty": "gas-phase cation, NO solvent/counter-ion screening -> "
                   "upper bound of the field; cluster, PBE, frozen metal; "
                   "answers only: does a REAL ion move the barrier as the "
                   "uniform-field slope predicts",
        "selftest": "K+ at 20 A must reproduce the bare barrier within ~0.02 eV",
    })
    res.setdefault("points", {})
    atomic_save(res)
    # середина C–C по реактанту (x=0 по построению build_ads, y=0), z≈средняя высота C
    zc = float(np.mean([g["reactant"][0][2][1][2], g["reactant"][0][4][1][2]]))
    cc_mid = (0.0, 0.0, zc)
    positions = {"selftest_20A": (-20.0, 0.0, zc)}
    if not args.smoke:
        for sgn, side in ((-1, "mx"), (+1, "px")):
            for R in (6.0, 8.0):
                positions[f"{side}_{R:.0f}A"] = (sgn * R, 0.0, zc)
    for name, kpos in positions.items():
        blk = res["points"].setdefault(name, {})
        if blk.get("barrier_eV") is not None:
            say(f"  {name}: есть (рестарт)"); continue
        exf, r = field_at_cc_mid(kpos, cc_mid)
        blk.update(k_pos=[round(v, 2) for v in kpos], r_from_ccmid_A=r,
                   field_x_at_ccmid_V_per_A=exf,
                   predicted_dbarrier_eV=round(SLOPE_EV_PER_VA * exf, 3))
        es = {}
        for tag in ("reactant", "ts"):
            t0 = time.time()
            e, conv = scf_with_k(g[tag][0], kpos)
            es[tag] = e
            blk[f"e_{tag}"] = round(e, 6); blk[f"conv_{tag}"] = conv
            say(f"  {name}/{tag}: E={e:.6f} conv={conv} ({time.time()-t0:.0f}s)")
            atomic_save(res)
        bar = (es["ts"] - es["reactant"]) * EV
        blk["barrier_eV"] = round(bar, 4)
        blk["dbarrier_vs_bare_eV"] = round(bar - BARE_BARRIER_EV, 4)
        say(f"  {name}: ΔE‡={bar:.4f} эВ (сдвиг {bar-BARE_BARRIER_EV:+.4f}; "
            f"поле {exf:+.3f} В/Å, прогноз слопа {blk['predicted_dbarrier_eV']:+.3f})")
        atomic_save(res)
    st = res["points"].get("selftest_20A", {})
    ok = st.get("barrier_eV") is not None and abs(st["dbarrier_vs_bare_eV"]) < 0.03
    res["selftest_pass"] = bool(ok)
    res["status"] = "smoke_ok" if args.smoke else (
        "ok" if all(p.get("barrier_eV") is not None for p in res["points"].values()) else "incomplete")
    atomic_save(res)
    say(f"selftest_pass={ok} status={res['status']} -> {OUT}")
    return 0 if (ok or not args.smoke) else 1

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    sys.exit(main(a))
