#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""routes/msa_persulfate_echem.py — электрохимический дескриптор персульфатного
инициатора: ГЕЙТ поворота T2 (APPLIED_UNEXPECTED_TWISTS.md), трек E
WEEK_PLAN_300VCPU.md «электро-инициация как дескриптор (потенциал генерации
SO₄•⁻)».

Вопрос T2: замыкается ли серно-радикальный цикл SO₂ → электросинтез
персульфата → инициатор MSA-цепи ПО ПОТЕНЦИАЛАМ. Персульфат синтезируют
анодно при ~2.0–2.2 В (лит., BDD/Pt в H₂SO₄); мы считаем, при каком
потенциале сам S₂O₈²⁻ отдаёт рабочий радикал SO₄•⁻ (диссоциативный перенос
электрона), т.е. сколько «электрохимического запаса» несёт инициатор.

Дескрипторы (все — UKS-PBE0/def2-SVP + density fitting, сольватация
ddCOSMO(вода, eps 78.3553), геометрии pyberny; СКРИН-УРОВЕНЬ):
  1) EA_vert / EA_adiab:  S₂O₈²⁻ + e⁻ → [S₂O₈³•⁻]@geom / → SO₄²⁻ + SO₄•⁻
     U_gen(В vs SHE) = EA_adiab(эВ) − 4.44
     (абсолютная шкала SHE 4.44 эВ — ЛИТЕРАТУРНАЯ константа, Trasatti/IUPAC).
     Кросс-маршрут без 4.44: CHE-цикл SO₄•⁻ + H⁺ + e⁻ → HSO₄⁻,
     U_CHE = −ΔG/e (vs SHE при pH 0 ≡ RHE).
  2) BDE(O–O): дианион S₂O₈²⁻ → 2 SO₄•⁻ (электрохимически релевантный) и
     нейтральный H₂S₂O₈ → 2 •OSO₃H — прямой кросс-чек с якорем 10.9 ккал/моль
     из msa_initiator (поле anchor_check).
  3) ΔG HAT: SO₄•⁻ + CH₄ → HSO₄⁻ + CH₃• (термодинамика на DFT-уровне;
     кинетический якорь HAT-барьера 8.7 ккал/моль уже есть в msa_initiator).

ЧЕСТНОСТЬ (не прятать):
  * только электронные энергии — без ZPE/термики (±0.1–0.2 эВ);
  * ddCOSMO — линейный континуум; для ди-/трианионов это самое слабое звено,
    явных противоионов нет; def2-SVP без диффузных функций для анионов —
    скрин-уровень. Итоговая честная полоса потенциалов: ±0.3–0.5 В;
  * абсолютный потенциал SHE 4.44 эВ — литература, не наш расчёт;
  * анодные 2.0–2.2 В электросинтеза персульфата — литература (tier literature);
  * якорь «NEVPT2-скрина» 10.9/8.7 ккал/моль — это DFT-строка того скрина:
    NEVPT2-колонка для HSO4 там не сошлась (CAS conv=false, −20.77) — честно
    помечено и здесь, и там.

env: MSAE_STAGE_TIMEOUT (сек на частицу, default 3600),
     MSAE_SPECIES (CSV-подмножество частиц), MSAE_BASIS (default def2-svp).

Запуск: python3 routes/msa_persulfate_echem.py [run|sanity]   (default run)
Выход:  routes/msa_persulfate_echem_results.json (поточечные чекпойнты, resume)
"""
import json
import math
import os
import signal
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "msa_persulfate_echem_results.json")

HARTREE_EV = 27.211386245988
HARTREE_KCAL = 627.509474
ABS_SHE_EV = 4.44          # литература: Trasatti/IUPAC абсолютный потенциал SHE
LIT_ANODIC_V = (2.0, 2.2)  # литература: анодный электросинтез S2O8^2- (BDD/Pt)
ANCHOR_BDE_KCAL = 10.9     # msa_initiator_results.json (DFT-строка NEVPT2-скрина)
ANCHOR_HAT_KCAL = 8.7      # там же: барьер HAT HSO4* + CH4 (DFT-строка)
UNC_V = 0.5                # консервативная полуширина честной полосы (0.3-0.5)

BASIS = os.environ.get("MSAE_BASIS", "def2-svp")
XC = "pbe0"
STAGE_TIMEOUT = int(os.environ.get("MSAE_STAGE_TIMEOUT", "3600"))

try:
    from pyscf import gto, dft, lib, solvent
    from pyscf.geomopt.berny_solver import optimize as berny_optimize
    PYSCF_OK, PYSCF_ERR = True, None
except Exception as _e:                      # noqa: BLE001 — sanity без pyscf
    PYSCF_OK, PYSCF_ERR = False, f"{type(_e).__name__}: {_e}"


def say(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# ---------------------------------------------------------------- геометрии
def h2_atoms():
    return [("H", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, 0.74))]


def h2o_atoms():
    return [("O", (0.0, 0.0, 0.0)), ("H", (0.0, 0.76, 0.59)),
            ("H", (0.0, -0.76, 0.59))]


def ch4_atoms():
    d, th = 1.09, math.radians(109.47)
    a = [("C", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, d))]
    for k in range(3):
        phi = math.radians(120 * k)
        a.append(("H", (d * math.sin(th) * math.cos(phi),
                        d * math.sin(th) * math.sin(phi), d * math.cos(th))))
    return a


def ch3_atoms():
    a = [("C", (0.0, 0.0, 0.0))]
    for k in range(3):
        phi = math.radians(120 * k)
        a.append(("H", (1.08 * math.cos(phi), 1.08 * math.sin(phi), 0.0)))
    return a


def so4_atoms():
    """SO4 тетраэдр (для SO4^2- и SO4^•-), r(S-O)=1.49."""
    r = 1.49 / math.sqrt(3.0)
    dirs = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    return [("S", (0.0, 0.0, 0.0))] + [("O", (r * x, r * y, r * z))
                                       for x, y, z in dirs]


def hso4_atoms():
    """HSO4 (для HSO4^- и радикала •OSO3H) — релакс. геометрия радикала из
    msa_initiator_bde.json (округлена)."""
    return [("S", (0.0, -0.0669, 0.0559)), ("O", (0.0, -0.1082, 1.4867)),
            ("O", (1.1339, 0.6141, -0.6428)), ("O", (-1.1339, 0.6142, -0.6427)),
            ("O", (0.0, -1.5925, -0.4142)), ("H", (0.0, -1.6366, -1.3844))]


def s2o8_atoms():
    """S2O8^2-: скелет h2s2o8 из msa_initiator без двух H (O-O мост r~1.44)."""
    a = []
    for sgn in (1, -1):
        x0 = sgn * 2.25
        a += [("S", (x0, 0.0, 0.0)),
              ("O", (sgn * 0.72, 0.0, 0.35)),            # мостовой O
              ("O", (x0 + sgn * 0.6, 1.30, -0.55)),
              ("O", (x0 + sgn * 0.6, -1.30, -0.55)),
              ("O", (x0 + sgn * 1.35, 0.0, 1.05))]
    return a


def h2s2o8_atoms():
    """HO3S-O-O-SO3H — как в msa_initiator.py (якорный кросс-чек)."""
    a = []
    for sgn in (1, -1):
        x0 = sgn * 2.25
        a += [("S", (x0, 0.0, 0.0)),
              ("O", (sgn * 0.72, 0.0, 0.35)),
              ("O", (x0 + sgn * 0.6, 1.30, -0.55)),
              ("O", (x0 + sgn * 0.6, -1.30, -0.55)),
              ("O", (x0 + sgn * 1.35, 0.0, 1.05)),
              ("H", (x0 + sgn * 2.05, 0.6, 1.35))]
    return a


# частица: (геометрия, charge, spin=2S, оптимизировать?, комментарий)
SPECIES = {
    "h2":          (h2_atoms,     0, 0, True,  "референс CHE (1/2 H2 = H+ + e-)"),
    "h2o":         (h2o_atoms,    0, 0, True,  "референс среды (CHE, вода)"),
    "ch4":         (ch4_atoms,    0, 0, True,  "субстрат HAT"),
    "ch3_rad":     (ch3_atoms,    0, 1, True,  "дублет, продукт HAT"),
    "so4_2m":      (so4_atoms,   -2, 0, True,  "SO4^2-, синглет"),
    "so4_rad_m":   (so4_atoms,   -1, 1, True,  "SO4^•-, ДУБЛЕТ (честный спин)"),
    "hso4_m":      (hso4_atoms,  -1, 0, True,  "HSO4^-, синглет"),
    "s2o8_2m":     (s2o8_atoms,  -2, 0, True,  "персульфат-дианион, синглет"),
    "h2s2o8":     (h2s2o8_atoms, 0, 0, True,  "нейтральная кислота Маршалла"),
    "hso4_rad":    (hso4_atoms,   0, 1, True,  "•OSO3H, дублет (якорный BDE)"),
    # вертикальная EA: SP трианиона-дублета на геометрии s2o8_2m (без опт.)
    "s2o8_3m_vert": (None,       -3, 1, False, "S2O8^3•- @ geom(s2o8_2m), SP"),
}
ORDER = ["h2", "h2o", "ch4", "ch3_rad", "so4_2m", "so4_rad_m", "hso4_m",
         "s2o8_2m", "s2o8_3m_vert", "h2s2o8", "hso4_rad"]


def atomic_save(obj, path):
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1)
    os.replace(tmp, path)


def with_timeout(fn, seconds):
    def handler(signum, frame):
        raise TimeoutError(f"stage timeout {seconds}s")
    old = signal.signal(signal.SIGALRM, handler)
    signal.alarm(seconds)
    try:
        return fn()
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


def mkmol(atoms, charge, spin):
    return gto.M(atom=atoms, basis=BASIS, charge=charge, spin=spin,
                 verbose=0, max_memory=24000)


def mkmf(mol):
    mf = dft.UKS(mol).density_fit()
    mf.xc = XC
    mf.conv_tol = 1e-8
    mf.max_cycle = 200
    mf = solvent.ddCOSMO(mf)          # solvation=ddCOSMO (вода) — скрин-уровень
    mf.with_solvent.eps = 78.3553
    return mf


def scf_e(mol, dm0=None):
    mf = mkmf(mol)
    e = mf.kernel(dm0=dm0) if dm0 is not None else mf.kernel()
    if not mf.converged:
        mf.level_shift = 0.3
        mf.max_cycle = 300
        e = mf.kernel(mf.make_rdm1())
    return float(e), bool(mf.converged), mf


def relax_e(atoms, charge, spin, maxsteps=100):
    """pyberny-оптимизация в ddCOSMO; при провале — честный SP на старте."""
    mol = mkmol(atoms, charge, spin)
    mf = mkmf(mol)
    try:
        meq = berny_optimize(mf, maxsteps=maxsteps)
        e, conv, _ = scf_e(meq)
        geo = [[meq.atom_symbol(i),
                [round(float(x), 6) for x in meq.atom_coord(i, unit="Angstrom")]]
               for i in range(meq.natm)]
        return e, conv, True, geo
    except Exception as exc:                                 # noqa: BLE001
        say(f"    relax FAILED ({type(exc).__name__}: {str(exc)[:80]}) -> SP")
        e, conv, _ = scf_e(mol)
        geo = [[s, [round(float(x), 6) for x in xyz]] for s, xyz in atoms]
        return e, conv, False, geo


def run_species(res, save):
    sp = res["species"]
    only = [s.strip() for s in os.environ.get("MSAE_SPECIES", "").split(",")
            if s.strip()]
    for tag in ORDER:
        if only and tag not in only:
            continue
        if sp.get(tag, {}).get("e") is not None:
            say(f"  {tag}: есть E={sp[tag]['e']:.6f} (resume)")
            continue
        geom_fn, charge, spin, do_opt, note = SPECIES[tag]
        if tag == "s2o8_3m_vert":
            base = sp.get("s2o8_2m", {})
            if base.get("e") is None or not base.get("xyz"):
                say(f"  {tag}: нет s2o8_2m — пропуск (вертикальная EA)")
                continue
            atoms = [(el, tuple(xyz)) for el, xyz in base["xyz"]]
        else:
            atoms = geom_fn()
        t0 = time.time()
        say(f"  {tag}: q={charge} 2S={spin} opt={do_opt} ({note})")
        try:
            if do_opt:
                e, conv, relaxed, geo = with_timeout(
                    lambda: relax_e(atoms, charge, spin), STAGE_TIMEOUT)
            else:
                def _sp():
                    mol = mkmol(atoms, charge, spin)
                    e_, c_, _ = scf_e(mol)
                    return e_, c_, False, [[s, list(x)] for s, x in atoms]
                e, conv, relaxed, geo = with_timeout(_sp, STAGE_TIMEOUT)
            sp[tag] = {"e": round(e, 8), "converged": conv, "charge": charge,
                       "spin_2S": spin, "relaxed": relaxed, "note": note,
                       "xyz": geo, "wall_s": round(time.time() - t0, 1)}
            say(f"    -> E={e:.6f} conv={conv} relaxed={relaxed} "
                f"({sp[tag]['wall_s']}s)")
        except (TimeoutError, Exception) as exc:             # noqa: BLE001
            sp[tag] = {"e": None, "failed": f"{type(exc).__name__}: "
                                            f"{str(exc)[:120]}",
                       "charge": charge, "spin_2S": spin,
                       "wall_s": round(time.time() - t0, 1)}
            say(f"    -> FAILED: {sp[tag]['failed']}")
        save()   # поточечный чекпойнт после КАЖДОЙ частицы


def ee(sp, tag):
    v = sp.get(tag, {})
    return v.get("e") if v.get("e") is not None and v.get("converged") else None


def build_descriptors(res, save):
    sp = res["species"]
    d = {}
    e_s2o8, e_so42, e_so4r = ee(sp, "s2o8_2m"), ee(sp, "so4_2m"), ee(sp, "so4_rad_m")
    e_vert, e_hso4, e_h2 = ee(sp, "s2o8_3m_vert"), ee(sp, "hso4_m"), ee(sp, "h2")
    e_ch4, e_ch3 = ee(sp, "ch4"), ee(sp, "ch3_rad")
    e_h2s2o8, e_hso4r = ee(sp, "h2s2o8"), ee(sp, "hso4_rad")

    # (1) диссоциативный перенос электрона -> потенциал генерации SO4^•-
    if None not in (e_s2o8, e_so42, e_so4r):
        ea_ad = (e_s2o8 - e_so42 - e_so4r) * HARTREE_EV
        u = ea_ad - ABS_SHE_EV
        d["EA_adiabatic_dissoc_eV"] = round(ea_ad, 3)
        d["U_gen_SO4rad_V_vs_SHE"] = {
            "value": round(u, 3),
            "range_V": [round(u - UNC_V, 2), round(u + UNC_V, 2)],
            "uncertainty": "±0.3–0.5 В (скрин-уровень: ddCOSMO-анионы, "
                           "без ZPE, def2-SVP без диффузных)",
            "abs_SHE_scale_eV": ABS_SHE_EV,
            "abs_SHE_tier": "literature (Trasatti/IUPAC), не наш расчёт",
            "cycle": "S2O8^2- + e- -> SO4^2- + SO4^•- (адиабатический, "
                     "CHE-подобный: U = EA_adiab/e - 4.44)"}
    if None not in (e_s2o8, e_vert):
        d["EA_vertical_eV"] = round((e_s2o8 - e_vert) * HARTREE_EV, 3)
        d["EA_vertical_note"] = ("вертикальный захват e- (S2O8^3•- @ geom "
                                 "дианиона) — нижняя грань; разница с адиаб. = "
                                 "движущая сила диссоциации O-O")
    # кросс-маршрут CHE без абсолютной шкалы: SO4^•- + H+ + e- -> HSO4^-
    if None not in (e_hso4, e_so4r, e_h2):
        dg = (e_hso4 - e_so4r - 0.5 * e_h2) * HARTREE_EV
        d["U_CHE_SO4rad_HSO4_V_vs_SHE"] = {
            "value": round(-dg, 3),
            "cycle": "SO4^•- + H+ + e- -> HSO4^- (CHE, pH0: SHE==RHE)",
            "note": "окислительная сила SO4^•-; лит. ориентир 2.4–3.1 В — "
                    "кросс-чек маршрута через 4.44 эВ"}

    # (2) BDE O-O на том же уровне + якорный кросс-чек
    if None not in (e_so4r, e_s2o8):
        d["bde_OO_S2O8_2m_kcal"] = round(
            (2 * e_so4r - e_s2o8) * HARTREE_KCAL, 2)
    if None not in (e_hso4r, e_h2s2o8):
        bde_n = (2 * e_hso4r - e_h2s2o8) * HARTREE_KCAL
        d["bde_OO_H2S2O8_kcal"] = round(bde_n, 2)
        d["anchor_check"] = {
            "anchor_kcal": ANCHOR_BDE_KCAL,
            "anchor_source": "msa_initiator_results.json (скрин инициаторов)",
            "anchor_tier": "DFT-строка NEVPT2-скрина: NEVPT2-колонка для HSO4 "
                           "там НЕ сошлась (CAS conv=false, -20.77) — якорь "
                           "именно UKS-PBE0/газ; здесь тот же гомолиз в "
                           "ddCOSMO — сольватационный сдвиг ожидаем",
            "this_work_kcal": round(bde_n, 2),
            "delta_kcal": round(bde_n - ANCHOR_BDE_KCAL, 2),
            "ok": bool(abs(bde_n - ANCHOR_BDE_KCAL) <= 6.0)}

    # (3) термодинамика HAT SO4^•- + CH4 -> HSO4^- + CH3^•
    if None not in (e_hso4, e_ch3, e_so4r, e_ch4):
        d["dG_HAT_SO4rad_CH4_kcal"] = {
            "value": round((e_hso4 + e_ch3 - e_so4r - e_ch4) * HARTREE_KCAL, 2),
            "note": "электронная ΔE реакции (термодинамика, НЕ барьер); "
                    f"кинетический якорь барьера HAT {ANCHOR_HAT_KCAL} "
                    "ккал/моль (HSO4* + CH4, msa_initiator, DFT-строка "
                    "NEVPT2-скрина) уже есть — здесь только кросс-чек знака "
                    "движущей силы"}

    d["lit_anodic_persulfate_V"] = {
        "range": list(LIT_ANODIC_V), "tier": "literature",
        "note": "типовой анодный потенциал электросинтеза S2O8^2- "
                "(BDD/Pt, конц. сульфатная среда)"}

    # вердикт T2 по потенциалам
    ug = d.get("U_gen_SO4rad_V_vs_SHE", {})
    if ug:
        lo, hi = ug["range_V"]
        closes = bool(hi < LIT_ANODIC_V[0])
        hat = d.get("dG_HAT_SO4rad_CH4_kcal", {}).get("value")
        hat_ok = (hat is not None and hat < 5.0)
        d["verdict"] = {
            "cycle_closes_by_potentials": closes,
            "logic": "гейт T2: анод делает персульфат при 2.0–2.2 В; если "
                     "потенциал генерации SO4^•- (восстановительный запуск "
                     "инициатора) даже с честной полосой ±0.3–0.5 В лежит НИЖЕ "
                     "анодного окна — цикл не закорачивается: то, что сделали "
                     "на аноде, остаётся окислителем/инициатором в объёме",
            "U_gen_range_V_vs_SHE": [lo, hi],
            "anodic_window_V": list(LIT_ANODIC_V),
            "hat_thermo_ok": hat_ok,
            "text": ("ЦИКЛ СХОДИТСЯ по потенциалам" if closes and hat_ok else
                     "ЦИКЛ ПОД ВОПРОСОМ — см. числа; NO-VERDICT легален") +
                    " (скрин-уровень, полоса ±0.3–0.5 В; финал — NEVPT2 + "
                    "явная сольватация)"}
    else:
        d["verdict"] = {"cycle_closes_by_potentials": None,
                        "text": "NO VERDICT: не хватает сошедшихся частиц"}

    res["descriptors"] = d
    res["status"] = ("ok" if ug and "bde_OO_H2S2O8_kcal" in d else "incomplete")
    save()


def stage_sanity():
    print(f"pyscf: {'OK' if PYSCF_OK else 'MISSING (' + str(PYSCF_ERR) + ')'}")
    for tag in ORDER:
        geom_fn, charge, spin, do_opt, note = SPECIES[tag]
        nat = len(geom_fn()) if geom_fn else "geom(s2o8_2m)"
        line = f"[sanity] {tag}: q={charge} 2S={spin} atoms={nat} opt={do_opt}"
        if PYSCF_OK and geom_fn:
            try:
                mol = mkmol(geom_fn(), charge, spin)
                line += f" nao={mol.nao} nelec={mol.nelectron} OK"
            except Exception as e:                           # noqa: BLE001
                line += f" FAIL: {e}"
        print(line, flush=True)


def stage_run():
    if not PYSCF_OK:
        sys.exit(f"pyscf недоступен: {PYSCF_ERR}")
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    res.setdefault("purpose",
        "Гейт поворота T2 (SO2 -> электросинтез персульфата -> инициатор "
        "MSA-цепи): потенциал генерации SO4^•- из S2O8^2- vs анодное окно "
        "электросинтеза 2.0-2.2 В — сходится ли серно-радикальный цикл по "
        "потенциалам. Трек E WEEK_PLAN_300VCPU.")
    res.setdefault("honesty",
        "СКРИН-УРОВЕНЬ: UKS-PBE0/def2-SVP+DF, сольватация ddCOSMO (вода, "
        "линейный континуум) — для ди/трианионов слабое звено; без ZPE/термики; "
        "без явных противоионов; def2-SVP без диффузных для анионов; "
        "абсолютная шкала SHE 4.44 эВ и анодные 2.0-2.2 В — литература. "
        "Честная полоса потенциалов ±0.3-0.5 В. Спины честно: SO4^•- и "
        "S2O8^3•- — дублеты. Якорь 10.9/8.7 ккал/моль = DFT-строка "
        "NEVPT2-скрина msa_initiator (NEVPT2 для HSO4 там не сошёлся).")
    res.setdefault("level", {
        "method": "UKS-PBE0", "basis": BASIS, "density_fit": True,
        "solvation": "ddCOSMO", "solvent": "water(eps=78.3553)",
        "geometry": "pyberny", "tier": "screen"})
    res.setdefault("species", {})
    save = lambda: atomic_save(res, OUT)                     # noqa: E731
    save()
    say(f"=== msa_persulfate_echem === threads={lib.num_threads()} "
        f"basis={BASIS} stage_timeout={STAGE_TIMEOUT}s")
    run_species(res, save)
    build_descriptors(res, save)
    say(f"status={res['status']} -> {OUT}")
    dd = res.get("descriptors", {})
    say(f"U_gen={dd.get('U_gen_SO4rad_V_vs_SHE', {}).get('value')} В vs SHE  "
        f"BDE(S2O8^2-)={dd.get('bde_OO_S2O8_2m_kcal')}  "
        f"BDE(H2S2O8)={dd.get('bde_OO_H2S2O8_kcal')} (якорь {ANCHOR_BDE_KCAL})  "
        f"verdict={dd.get('verdict', {}).get('cycle_closes_by_potentials')}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    {"run": stage_run, "sanity": stage_sanity}[mode]()
