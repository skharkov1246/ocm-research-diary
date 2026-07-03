#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/h2_oer_quantum.py — Этап 3 водородного трека: мультиреференс-диагностика
и квантовая поправка ЛИМИТИРУЮЩЕЙ стадии OER на ОПТИМИЗИРОВАННЫХ геометриях.

Связывает Этап 1 (VQE↔CASCI фундамент) с Этапом 2 (ΔG-лестница): берёт из
h2_oer_ladder_results.json лимитирующую стадию и лучшие (низшие по спину)
оптимизированные геометрии её частиц, и на каждой считает:
  ROHF → AVAS(Ni 3d + O 2p адсорбата) → CASCI:
    • NOON → n_u («лишние» неспаренные e⁻) — мультиреференсна ли частица;
    • E_static = E_CASCI − E_ROHF — статическая корреляция в активном пространстве;
  поправка стадии: ΔΔE = E_static(продукт) − E_static(реагент)
  (газовые H₂O/H₂ закрытооболочечные — их статическая корреляция ≈ 0, что
   проверяется явно, а не постулируется).

Опционально (H2_VQE=1): UCCSD-VQE (PennyLane lightning.qubit) на активном
пространстве частицы с наибольшим n_u, если ≤ H2_VQE_MAXQ кубитов (деф. 16) —
PoC «квантовый компьютер решил бы этот кусок», сверка с CASCI.

ЧЕСТНО: ΔΔE — диагностическая поправка (CASCI-в-AVAS поверх ROHF, без
динамической корреляции; PBE-лестница уже содержит часть корреляции — прямое
сложение чисел из разных методов не делаем, поправку репортим ОТДЕЛЬНО).
Тот же подход, что Stage 4 cc-activation (там поправка вышла пренебрежимой).

Запуск: python3 -u routes/h2_oer_quantum.py [--source mono|di]
Результат: routes/h2_oer_quantum_results.json
"""
import os, sys, json, time, argparse
import numpy as np
from pyscf import gto, scf, mcscf
from pyscf.mcscf import avas

HARTREE_EV = 27.211386245988
KCAL = 627.5094740631
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "h2_oer_quantum_results.json")

STEP_SPECIES = {  # стадия → (реагент, продукт) среди кластерных частиц
    "dG1_bare_to_OH": ("bare", "oh"), "dG2_OH_to_O": ("oh", "o"),
    "dG3_O_to_OOH": ("o", "ooh"), "dG4_OOH_to_O2": ("ooh", "bare"),
}


def _mol_from_run(run, basis):
    """Восстановить молекулу из сохранённого рана (atoms | coords+species)."""
    if "atoms" in run:
        atoms = [(a[0], (a[1], a[2], a[3])) for a in run["atoms"]]
    else:
        sys.path.insert(0, HERE)
        from h2_oer_ladder import SPECIES
        syms = [a[0] for a in SPECIES[run["species"]]["atoms"]]
        atoms = [(s, tuple(c)) for s, c in zip(syms, run["coords"])]
    return gto.M(atom=atoms, basis=basis, charge=0, spin=run["spin"],
                 verbose=0, unit="Angstrom")


def _save_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1, ensure_ascii=False)
    os.replace(tmp, path)


def _rohf(mol):
    """ROHF: Newton → level-shift DIIS → warm-start Newton. Незасошедшийся
    SCF — ошибка (честно), а не тихий вход в AVAS/CASCI."""
    mf = scf.ROHF(mol); mf.verbose = 0; mf.conv_tol = 1e-8
    mfn = mf.newton()
    try:
        mfn.kernel()
    except Exception:
        mfn.converged = False
    if getattr(mfn, "converged", False):
        return mfn
    mf2 = scf.ROHF(mol); mf2.verbose = 0
    mf2.level_shift = 0.5; mf2.max_cycle = 400; mf2.conv_tol = 1e-7
    mf2.kernel()
    mfn = mf2.newton()
    try:
        mfn.kernel(mf2.make_rdm1())
    except Exception:
        mfn.converged = False
    if not getattr(mfn, "converged", False):
        raise RuntimeError("ROHF unconverged (Newton + level-shift + re-Newton)")
    return mfn


def casci_probe(run, basis):
    """ROHF→AVAS(Ni 3d + адсорбат-O 2p)→CASCI на сохранённой геометрии."""
    t0 = time.time()
    mol = _mol_from_run(run, basis)
    mf = _rohf(mol)
    labels = []
    has_ni = any(mol.atom_symbol(i) == "Ni" for i in range(mol.natm))
    if has_ni:
        labels.append("Ni 3d")
        # адсорбатные O: всё после первых 5 атомов ядра [Ni(OH)2] (моно) —
        # либо явные атомы после ядра (ди, длина ядра в run["core_len"])
        core_len = run.get("core_len", 5)
        labels += [f"{i} O 2p" for i in range(core_len, mol.natm)
                   if mol.atom_symbol(i) == "O"]
    else:
        # газовые контроли: связывающие+разрыхляющие, а не тривиально
        # заполненное O-2p (в нём корреляция нулевая по построению)
        labels = ["O 2p", "H 1s"]
    if mol.natm == 2 and set(mol.atom_symbol(i) for i in range(2)) == {"H"}:
        labels = ["H 1s"]
    # openshell_option=3: дефолт (2) на ROHF умеет терять SOMO-структуру —
    # CASCI выходил ВЫШЕ ROHF на +41 ккал/моль (нарушение вариационности)
    ncas, nelec, orbs = avas.avas(mf, labels, threshold=0.2,
                                  openshell_option=3, verbose=0)
    mc = mcscf.CASCI(mf, ncas, nelec); mc.verbose = 0
    mc.fcisolver.conv_tol = 1e-10
    mc.kernel(orbs)
    if mc.e_tot - mf.e_tot > 1e-3 / 627.5:   # вариационный сторож
        raise RuntimeError(f"CASCI above ROHF by {(mc.e_tot-mf.e_tot)*627.5:.2f} "
                           "kcal/mol — active space is inconsistent")
    dm = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
    noon = np.sort(np.linalg.eigvalsh(dm))[::-1]
    spin = mol.spin
    n_u_formal = spin
    n_u = float(np.sum(np.minimum(noon, 2 - noon)))
    return dict(e_rohf=float(mf.e_tot), e_casci=float(mc.e_tot),
                casci_converged=bool(mc.converged),
                e_static_ha=float(mc.e_tot - mf.e_tot),
                cas=[int(sum(nelec) if hasattr(nelec, "__len__") else nelec),
                     int(ncas)],
                qubits=2 * int(ncas),
                noon=[round(float(x), 3) for x in noon],
                n_unpaired=round(n_u, 3),
                n_unpaired_excess=round(n_u - n_u_formal, 3),
                spin=spin, seconds=round(time.time() - t0, 1))


def vqe_poc(run, basis, maxq, iters=400):
    """CASCI-интегралы → JW-кубитный гамильтониан → ворота корректности
    (низшее собственное значение H == E_CASCI) → UCCSD-VQE (lightning.qubit).
    Восстановленный рецепт Stage-1 NiO (тот, что гонялся на AWS)."""
    import pennylane as qml
    import openfermion as of
    from openfermion.chem.molecular_data import spinorb_from_spatial
    from pyscf import ao2mo
    from scipy.sparse.linalg import eigsh

    t0 = time.time()
    mol = _mol_from_run(run, basis)
    mf = _rohf(mol)
    core_len = run.get("core_len", 5)
    labels = ["Ni 3d"] + [f"{i} O 2p" for i in range(core_len, mol.natm)
                          if mol.atom_symbol(i) == "O"]
    ncas, nelec, orbs = avas.avas(mf, labels, threshold=0.2,
                                  openshell_option=3, verbose=0)
    nq = 2 * int(ncas)
    if nq > maxq:
        return {"skipped": f"{nq}q > H2_VQE_MAXQ={maxq}"}
    mc = mcscf.CASCI(mf, ncas, nelec); mc.verbose = 0; mc.kernel(orbs)
    e_casci = float(mc.e_tot)

    h1e, ecore = mc.get_h1eff()
    h2e = ao2mo.restore(1, mc.get_h2eff(), ncas)
    # PySCF chemist (pq|rs) → OpenFermion ⟨pq|rs⟩: transpose(0,2,3,1)
    one_s, two_s = spinorb_from_spatial(h1e, h2e.transpose(0, 2, 3, 1))
    ham = of.InteractionOperator(float(ecore), one_s, 0.5 * two_s)
    qubit_op = of.jordan_wigner(of.get_fermion_operator(ham))
    H = qml.qchem.import_operator(qubit_op, format="openfermion")

    # ворота корректности: диагонализация JW-гамильтониана В ПРАВИЛЬНОМ
    # (N, Sz)-СЕКТОРЕ vs CASCI. Полный фоковский минимум может лежать в чужом
    # секторе частиц/спина — сравнение с ним было бы ложной «ошибкой ворот».
    # На ≥18q sparse-сборка ест >12 ГБ (урок Stage 1: OOM на 2¹⁸) — там ворота
    # пропускаем: рецепт валидирован на меньших пространствах.
    na, nb = int(mc.nelecas[0]), int(mc.nelecas[1])
    ne_tot, sz = na + nb, 0.5 * (na - nb)
    gate_maxq = int(os.environ.get("H2_GATE_MAXQ", "16"))
    e_exact, gate_kcal = None, None
    if nq <= gate_maxq:
        from openfermion.linalg import jw_sz_restrict_operator
        sp_ham = of.get_sparse_operator(qubit_op, n_qubits=nq)
        sp_sector = jw_sz_restrict_operator(sp_ham, sz_value=sz,
                                            n_electrons=ne_tot, n_qubits=nq)
        if sp_sector.shape[0] <= 2:
            e_exact = float(np.linalg.eigvalsh(sp_sector.toarray())[0])
        else:
            e_exact = float(eigsh(sp_sector, k=1, which="SA",
                                  return_eigenvectors=False)[0])
        gate_kcal = abs(e_exact - e_casci) * KCAL

    # референс и возбуждения — в ТОМ ЖЕ секторе: openfermion/pennylane
    # интерливинг (чётные = α); α занимает na низших орбиталей, β — nb.
    # Это ROHF-детерминант при упорядоченных активных орбиталях (docc→socc).
    init_state = np.zeros(nq, dtype=int)
    init_state[[2 * i for i in range(na)]] = 1
    init_state[[2 * i + 1 for i in range(nb)]] = 1
    occ = [i for i in range(nq) if init_state[i]]
    virt = [i for i in range(nq) if not init_state[i]]
    szf = lambda i: 0.5 if i % 2 == 0 else -0.5
    singles = [[o, v] for o in occ for v in virt if szf(o) == szf(v)]
    doubles = [[o1, o2, v1, v2]
               for ii, o1 in enumerate(occ) for o2 in occ[ii + 1:]
               for jj, v1 in enumerate(virt) for v2 in virt[jj + 1:]
               if abs(szf(o1) + szf(o2) - szf(v1) - szf(v2)) < 1e-9]
    s_w, d_w = qml.qchem.excitations_to_wires(singles, doubles)
    dev = qml.device("lightning.qubit", wires=nq)

    @qml.qnode(dev, diff_method="adjoint")
    def cost(w):
        qml.UCCSD(w, wires=range(nq), s_wires=s_w, d_wires=d_w,
                  init_state=init_state)
        return qml.expval(H)

    w = qml.numpy.zeros(len(singles) + len(doubles), requires_grad=True)
    opt = qml.AdamOptimizer(0.05); e_prev = 1e9; hist = []
    for it in range(iters):
        w, e = opt.step_and_cost(cost, w)
        if it % 25 == 0:
            hist.append(round(float(e), 6))
        if it > 40 and abs(e - e_prev) < 1e-7:
            break
        e_prev = e
    e_vqe = float(cost(w))
    return dict(qubits=nq, cas=[ne_tot, int(ncas)], terms=len(qubit_op.terms),
                e_casci=e_casci, e_exact_jw=e_exact,
                gate_exact_vs_casci_kcal=(None if gate_kcal is None
                                          else round(gate_kcal, 4)),
                e_vqe=e_vqe, vqe_minus_casci_kcal=round((e_vqe - e_casci) * KCAL, 2),
                iters=it + 1, energy_trace=hist,
                seconds=round(time.time() - t0, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", choices=("mono", "di"), default="mono")
    ap.add_argument("--basis", default="def2-svp")
    ap.add_argument("--step", default=None, help="переопределить стадию")
    args = ap.parse_args()

    src = os.path.join(HERE, "h2_oer_ladder_results.json" if args.source == "mono"
                       else "h2_oer_dinuclear_results.json")
    with open(src) as f:
        data = json.load(f)
    ladder = data.get("ladder") or {}
    # лестница обязана быть того же базиса (и гео-опт для моно) — иначе
    # лимитирующая стадия могла прийти из другого расчёта
    if ladder and ladder.get("basis") != args.basis:
        print(f"[fatal] лестница в {src} посчитана в {ladder.get('basis')}, "
              f"а запрошен {args.basis}"); return
    if ladder and args.source == "mono" and not ladder.get("geom_opt", False):
        print("[fatal] лестница в", src, "не гео-оптимизирована"); return
    step = args.step or ladder.get("limiting")
    if not step:
        print("[fatal] нет лестницы в", src); return
    rname, pname = STEP_SPECIES[step]

    # лучшие раны частиц стадии + газовые контроли (только успешные,
    # только режим гео-опт — ':sp'-записи отсеиваются по ключу)
    best = {}
    for k, r in data["runs"].items():
        if not k.startswith(args.basis) or "e_ha" not in r or k.endswith(":sp"):
            continue
        sp = r.get("species") or k.split(":")[1]
        if args.source == "mono":
            if sp in ("h2o", "h2"):
                if r.get("opt") is not True:
                    continue
            elif r.get("opt") is not True:
                continue
        else:
            if sp == "bare":
                if not r.get("opt_ok"):
                    continue
            elif not r.get("relax_ok"):
                continue
        if sp not in best or r["e_ha"] < best[sp]["e_ha"]:
            best[sp] = r
    if args.source == "di":
        for r in best.values():
            r["core_len"] = 14

    out = {"meta": {"source": args.source, "basis": args.basis, "step": step,
                    "honest": "CASCI-в-AVAS диагностика; поправка репортится "
                              "отдельно от PBE-лестницы, не суммируется"},
           "species": {}}
    for sp in (rname, pname, "h2o", "h2"):
        if sp not in best:
            if sp in ("h2o", "h2") and args.source == "di":
                with open(os.path.join(HERE, "h2_oer_ladder_results.json")) as f:
                    mruns = json.load(f)["runs"]
                cand = [r for k, r in mruns.items()
                        if k.startswith(args.basis) and "e_ha" in r
                        and (r.get("species") or k.split(":")[1]) == sp]
                if cand:
                    best[sp] = min(cand, key=lambda r: r["e_ha"])
                    best[sp]["species"] = sp
                else:
                    print(f"[warn] нет {sp}"); continue
            else:
                print(f"[warn] нет {sp}"); continue
        print(f"[cas ] {sp} (2S={best[sp]['spin']}) …", flush=True)
        try:
            out["species"][sp] = casci_probe(best[sp], args.basis)
        except Exception as exc:
            out["species"][sp] = {"error": f"{type(exc).__name__}: {exc}"}
        _save_json(RESULTS, out)
        r = out["species"][sp]
        if "e_static_ha" in r:
            print(f"       CAS{tuple(r['cas'])}={r['qubits']}q  "
                  f"E_static={r['e_static_ha']*KCAL:.2f} kcal/mol  "
                  f"n_u={r['n_unpaired']} (сверх формального: "
                  f"{r['n_unpaired_excess']})", flush=True)
        else:
            print("       ", r.get("error"), flush=True)

    rr, pp = out["species"].get(rname, {}), out["species"].get(pname, {})
    if "e_static_ha" in rr and "e_static_ha" in pp:
        dd = pp["e_static_ha"] - rr["e_static_ha"]
        out["step_correction"] = {
            "step": step, "reactant": rname, "product": pname,
            "ddE_static_eV": round(dd * HARTREE_EV, 3),
            "ddE_static_kcal": round(dd * KCAL, 2),
            "note": "= E_static(prod) − E_static(react); диагностика того, "
                    "насколько многодетерминантность сдвигает лимитирующий шаг"}
        _save_json(RESULTS, out)
        print(json.dumps(out["step_correction"], indent=1, ensure_ascii=False))

    if os.environ.get("H2_VQE") == "1":
        maxq = int(os.environ.get("H2_VQE_MAXQ", "16"))
        target = max((s for s in (rname, pname) if s in out["species"]
                      and "n_unpaired" in out["species"][s]),
                     key=lambda s: out["species"][s]["n_unpaired"], default=None)
        if target:
            print(f"[vqe ] target={target}")
            try:
                out["vqe"] = vqe_poc(best[target], args.basis, maxq)
            except Exception as exc:
                out["vqe"] = {"error": f"{type(exc).__name__}: {exc}"}
            _save_json(RESULTS, out)


if __name__ == "__main__":
    main()
