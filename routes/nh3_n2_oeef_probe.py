#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""routes/nh3_n2_oeef_probe.py — ШОРТКАТ-ЗОНД (ammonia): электрическая накачка
активации N₂. Кросс-рычаг «поле × N≡N» на минимальном нитрогеназном мотиве.

Идея: активация N₂ на Fe идёт через обратную донацию d(Fe)→π*(N₂); π* полярна
вдоль оси Fe–N–N → внешнее поле должно накачивать/сливать заселённость π*,
управляя степенью активации (удлинение R(N–N), перенос заряда на N₂). Если
dR_NN/dF заметен — поле = ручка активации N≡N при мягких условиях (электро-
статический шорткат вместо 200 атм/450°C; двойной слой на электроде даёт
0.1–1 В/Å бесплатно).

Модель v1 (строго говоря минимальная): [Fe(SH)₃(N₂)]⁻ — Fe(II) d6 на трёх тиолах
(мотив Fe–S окружения нитрогеназы), N₂ end-on аксиально. Геометрия идеализи-
рована по лит. мотивам (Fe–S 2.30 Å тригонально, Fe–N 1.80, N–N старт 1.12,
S–H 1.34; НЕ оптимизирована — фиксируем каркас, сканируем ТОЛЬКО R(N–N):
1.06→1.30 Å, 7 точек). UKS-PBE/def2-SVP+DF, спины 2S=4 и 2S=2 (низший берём).
Поля вдоль оси Fe–N–N (z): 0, ±0.005, ±0.010 а.е. (±0.26/±0.51 В/Å).
Дескрипторы: R_min(F) — положение минимума кривой E(R_NN) (параболическая
интерполяция трёх нижних точек), Малликен-заряд фрагмента N₂ в минимуме.
САМОТЕСТ: при F=0 R_min обязан быть 1.10–1.25 Å (активированный, но целый N₂).

Ограничения: каркас заморожен и идеализирован (не оптимум!) — сдвиги ПО ПОЛЮ при
фиксированном каркасе надёжны как чувствительность, абсолютный R_min — нет;
кластер без белка/растворителя; PBE; заряд −1 в вакууме. Зонд отвечает на один
вопрос: есть ли у поля первопорядковая ручка на активацию N₂.

Запуск: OMP_NUM_THREADS=4 python3 -u routes/nh3_n2_oeef_probe.py
Smoke:  --smoke (F=0, один спин, 3 точки R; ~10-20 мин)
Выход:  routes/nh3_n2_oeef_results.json (atomic, incremental)
"""
import os, json, time, argparse, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
from pyscf import gto, dft, lib  # noqa: E402

OUT = os.path.join(HERE, "nh3_n2_oeef_results.json")
EV = 27.211386245988
VA = 51.42206747
RGRID = [1.06, 1.09, 1.12, 1.16, 1.20, 1.25, 1.30]
FIELDS = [0.0, 0.005, -0.005, 0.010, -0.010]
SPINS = [4, 2]                      # 2S: квинтет-подобный d6 HS и триплет

def say(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

def atomic_save(obj):
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1)
    os.replace(tmp, OUT)

def build(r_nn):
    # Fe в нуле; 3×SH тригонально в плоскости xy чуть ниже; N2 вдоль +z (end-on)
    fs, fn = 2.30, 1.80
    tilt = np.deg2rad(100.0)        # S чуть под плоскостью (пирамидализация)
    atoms = [["Fe", (0.0, 0.0, 0.0)]]
    for k in range(3):
        a = 2 * np.pi * k / 3
        s = (fs * np.sin(tilt) * np.cos(a), fs * np.sin(tilt) * np.sin(a),
             -fs * np.cos(tilt))
        atoms.append(["S", s])
        sh = np.array(s); h = sh + np.array([0.9 * np.cos(a), 0.9 * np.sin(a), -1.0])
        h = sh + (h - sh) / np.linalg.norm(h - sh) * 1.34
        atoms.append(["H", tuple(h)])
    atoms.append(["N", (0.0, 0.0, fn)])
    atoms.append(["N", (0.0, 0.0, fn + r_nn)])
    return atoms

def scf_point(r_nn, f_au, spin, dm0=None):
    mol = gto.M(atom=build(r_nn), basis="def2-svp", charge=-1, spin=spin,
                verbose=0)
    mf = dft.UKS(mol).density_fit()
    mf.xc = "pbe"; mf.conv_tol = 1e-8; mf.max_cycle = 200; mf.level_shift = 0.2
    if abs(f_au) > 0:
        mol.set_common_orig((0.0, 0.0, 0.0))
        ao_dip = mol.intor_symmetric("int1e_r", comp=3)
        h0 = mf.get_hcore()
        E = np.array([0.0, 0.0, f_au])
        mf.get_hcore = lambda *a: h0 + np.einsum("x,xij->ij", E, ao_dip)
    e = mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel()
    if not mf.converged:
        mfn = mf.newton(); mfn.max_cycle = 80
        e = mfn.kernel(mf.mo_coeff, mf.mo_occ)
        mf.mo_coeff, mf.mo_occ = mfn.mo_coeff, mfn.mo_occ
        mf.converged = bool(mfn.converged); e = float(mfn.e_tot)
    # Малликен-заряд фрагмента N2 (последние 2 атома)
    try:
        pop = mf.mulliken_pop(verbose=0)[1]
        qn2 = float(pop[-1] + pop[-2])
    except Exception:
        qn2 = None
    return float(e), bool(mf.converged), qn2, mf.make_rdm1()

def parabolic_min(rs, es):
    i = int(np.argmin(es))
    if i == 0 or i == len(rs) - 1: return float(rs[i]), False
    x = np.array(rs[i-1:i+2]); y = np.array(es[i-1:i+2])
    c = np.polyfit(x, y, 2)
    return float(-c[1] / (2 * c[0])), True

def main(args):
    say(f"N2-activation OEEF probe — threads={lib.num_threads()}")
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    res.setdefault("meta", {
        "what": "[Fe(SH)3(N2)]- end-on; scan R(N-N) in axial field; descriptors "
                "R_min(F) (parabolic) and Mulliken q(N2) near min; UKS-PBE/"
                "def2-SVP+DF, spins 2S=4,2 (lowest taken)",
        "honesty": "idealized frozen scaffold (NOT optimized) — only FIELD-"
                   "induced shifts at fixed scaffold are honest; gas-phase "
                   "anion, no protein/solvent; PBE; v1 sensitivity probe",
        "selftest": "F=0: R_min must fall in 1.10-1.25 A",
    })
    res.setdefault("scan", {})
    fields = [0.0] if args.smoke else FIELDS
    rgrid = RGRID[1:6] if args.smoke else RGRID
    spins = [SPINS[0]] if args.smoke else SPINS
    for f_au in fields:
        fkey = f"{f_au:+.3f}"
        blk = res["scan"].setdefault(fkey, {})
        for spin in spins:
            skey = f"2S={spin}"
            sb = blk.setdefault(skey, {"points": {}})
            dm = None
            for r in rgrid:
                rkey = f"{r:.2f}"
                if rkey in sb["points"] and sb["points"][rkey].get("converged"):
                    continue
                t0 = time.time()
                try:
                    e, conv, qn2, dm = scf_point(r, f_au, spin, dm0=dm)
                    sb["points"][rkey] = {"e": e, "converged": conv,
                                          "q_N2": qn2,
                                          "t_s": round(time.time() - t0, 1)}
                    say(f"  F={fkey} {skey} R={rkey}: E={e:.5f} conv={conv} "
                        f"q(N2)={qn2 if qn2 is None else round(qn2,3)} "
                        f"({sb['points'][rkey]['t_s']}s)")
                except Exception as exc:
                    sb["points"][rkey] = {"converged": False,
                                          "error": f"{type(exc).__name__}: {str(exc)[:80]}"}
                    say(f"  F={fkey} {skey} R={rkey}: FAILED {exc}")
                atomic_save(res)
            ok = [(float(k), v["e"]) for k, v in sb["points"].items()
                  if v.get("converged")]
            if len(ok) >= 3:
                ok.sort()
                rmin, interior = parabolic_min([p[0] for p in ok],
                                               [p[1] for p in ok])
                sb["r_min_A"] = round(rmin, 4); sb["min_interior"] = interior
                i = int(np.argmin([p[1] for p in ok]))
                sb["e_min"] = ok[i][1]
                qk = f"{ok[i][0]:.2f}"
                sb["q_N2_at_min"] = sb["points"][qk].get("q_N2")
            atomic_save(res)
    # сводка: низший спин на каждом поле, R_min(F), q(F)
    rows = []
    for fkey, blk in sorted(res["scan"].items(), key=lambda kv: float(kv[0])):
        best = None
        for skey, sb in blk.items():
            if sb.get("e_min") is not None and (best is None or sb["e_min"] < best["e_min"]):
                best = {"spin": skey, **{k: sb[k] for k in
                        ("e_min", "r_min_A", "min_interior", "q_N2_at_min")}}
        if best:
            rows.append({"field_au": float(fkey),
                         "field_V_per_A": round(float(fkey) * VA, 3), **best})
    res["summary_rows"] = rows
    z = next((r for r in rows if r["field_au"] == 0.0), None)
    res["selftest_pass"] = bool(z and z["min_interior"] and 1.10 <= z["r_min_A"] <= 1.25)
    if len(rows) >= 3:
        F = np.array([r["field_V_per_A"] for r in rows])
        R = np.array([r["r_min_A"] for r in rows])
        res["dRmin_dF_A_per_VA"] = round(float(np.polyfit(F, R, 1)[0]), 5)
        q = [r["q_N2_at_min"] for r in rows if r["q_N2_at_min"] is not None]
        if len(q) == len(rows):
            res["dqN2_dF_e_per_VA"] = round(float(np.polyfit(F, np.array(q), 1)[0]), 4)
    res["status"] = "smoke_ok" if args.smoke else (
        "ok" if len(rows) == len(FIELDS) else "incomplete")
    atomic_save(res)
    say(f"selftest={res.get('selftest_pass')} status={res['status']} -> {OUT}")
    return 0

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    raise SystemExit(main(a))
