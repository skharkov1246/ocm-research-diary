#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/cler_ladder.py — гейт кейса N1-v2 (аноды Кольской ГМК).

Зачем. Главный электролизный передел Кольской ГМК (~120 кт Ni/год, ~500 ванн)
работает на ХЛОРНОМ анолите: анодная реакция — выделение Cl₂, а НЕ кислорода.
Норникель прямо сейчас (ванна №304, апрель–октябрь 2026) испытывает Ti-аноды
с покрытием из смешанных оксидов Ru-Pd-Ir, заменяя часть иридия собственным
палладием, с заявкой −5–10 % электропотребления. Наш кейс — не «продать идею»,
а дать расчётную поддержку ИХ живой программе и ответить на их же вопрос:
СКОЛЬКО Ir МОЖНО УБРАТЬ, не потеряв перенапряжение и селективность по Cl₂.
(Контекст: APPLIED_PRIOR_ART_SCAN.md, раздел 3 — переформулировка N1.)

Химия. Двухэлектронный путь Фольмера–Гейровского на рутиловом оксидном
центре M–O (computational chlorine electrode, аналог CHE):

    * + Cl⁻            → Cl* + e⁻        ΔG1
    Cl* + Cl⁻          → * + Cl₂ + e⁻    ΔG2
    ------------------------------------------------
    2Cl⁻ → Cl₂ + 2e⁻,  ΔG1 + ΔG2 = 2·1.36 эВ  (замыкание по U_eq)

Референс электрона/иона берём через Cl₂ (равновесие при U_eq = 1.36 В):
    G(Cl⁻) − G(e⁻)|_U=0 = ½G(Cl₂) + 1.36 эВ,
откуда с дескриптором связывания хлора

    ΔG_Cl* = G(Cl*) − G(*) − ½G(Cl₂)          ← классический вулкан ClER
    ΔG1 = ΔG_Cl* + 1.36 эВ
    ΔG2 = 1.36 эВ − ΔG_Cl*
    η_ClER = max(ΔG1, ΔG2)/e − 1.36 В = |ΔG_Cl*|

Замыкание ΔG1+ΔG2 = 2.72 эВ здесь ТОЖДЕСТВЕННО (не подгонка): оно встроено
в определение референса, а вся физика центра сидит в одном числе ΔG_Cl*.

Конкуренция с OER на ТОМ ЖЕ центре (паразитный O₂ съедает ток и разрушает
покрытие) — стандартная CHE-лестница, как в routes/oer_sulfate_ladder.py:
    ΔG_OH*  = G(OH*)  − G(*) − G(H₂O) + ½G(H₂)
    ΔG_O*   = G(O*)   − G(*) − G(H₂O) +  G(H₂)
    ΔG_OOH* = G(OOH*) − G(*) − 2G(H₂O) + 1.5G(H₂)
    шаги: ΔG_OH*, ΔG_O*−ΔG_OH*, ΔG_OOH*−ΔG_O*, 4.92−ΔG_OOH*;  η_OER = max/e − 1.23

ГЛАВНОЕ ЧИСЛО ДЖОБА — дескриптор селективности Cl₂ против O₂:

    dSel = ΔG_Cl* − ΔG_O*      [эВ]

Чем ОТРИЦАТЕЛЬНЕЕ dSel, тем прочнее центр держит хлор ОТНОСИТЕЛЬНО кислорода,
т.е. тем сильнее хлорный путь выигрывает у кислородного. Абсолютная шкала
кластерная и смещена — работаем РАНЖИРОМ центров, не абсолютом.

Центры (рутиловый мотив MO₂, минимальный узел; активный сайт — metals[0]):
    ru        [Ru]           чистый RuO₂ — промышленный референс DSA
    ir        [Ir]           IrO₂ — дорогой стабильный
    pd        [Pd]           «чистый PdO» — предел «максимум своего металла»
    ru_pd     [Ru, Pd]       Ru с соседним Pd — их состав
    ru_pd_ir  [Ru, Pd, Ir]   тройной — состав ванны №304

Геометрия строится генератором build_cluster(): цепочка рёберно-связанных
MO₆-узлов вдоль оси z (рутиловый мотив, M–M 3.05 Å), соседние металлы связаны
парой мостиков (μ-O + μ-OH), длины M–O 1.70–2.15 Å (оксо 1.70, гидроксо 1.95,
мостик 1.98, аква 2.15, M–Cl 2.28). Формальные заряды подобраны так, что ВСЕ
кластеры НЕЙТРАЛЬНЫ при M(IV):
    N=1: M(=O)(OH)₂(H₂O)         +4 −2 −1 −1 = 0
    N=2: M₀(=O)(OH)–(μO)(μOH)–M₁(=O)(H₂O)          +1 −3 +2 = 0
    N=3: M₀(=O)(OH)–(μO)(μOH)–M₁(OH)–(μO)(μOH)–M₂(=O)(H₂O)   +1 −3 +3 −3 +2 = 0
Адсорбат X (Cl*/O*/OH*/OOH*) садится на свободную апикальную позицию M₀ (−z).

Метод: UKS PBE(DF)/def2-SVP + def2-ECP для Ru(core 28)/Pd(core 28)/Ir(core 60)
(проверено: pyscf отдаёт def2-ECP по ecp="def2-svp"), pyberny-оптимизация
геометрии, скан спиновых состояний с автоматической проверкой ЧЁТНОСТИ по
реальному числу электронов (после ECP!), gradient-audit (gmax в JSON),
SCF-каскад Newton → level-shift DIIS → smearing → re-DIIS.

Честные оговорки (computed_here в meta):
  • электронные энергии и геометрии — реальный расчёт (computed_here: true);
  • ZPE−TS поправки — ЛИТЕРАТУРНЫЕ ЗАГЛУШКИ, НЕ считаны для этих кластеров
    (computed_here: false), источники в GCORR_SOURCE;
  • молекулярный кластер ≠ рутиловая поверхность (110): нет решётки, нет
    покрытия соседями, нет растворителя, двойного слоя и потенциала;
  • PBE переоценивает связь Cl₂ и занижает d-зазоры 4d/5d;
  ⇒ это ЗНАКОВЫЙ СКРИН РАНЖИРА центров, а НЕ абсолютные перенапряжения.

env-ручки:
  CLER_CENTERS        подмножество центров через запятую (default: все)
  CLER_BASIS          базис (default def2-svp)
  CLER_XC             функционал pyscf (default "pbe,pbe")
  CLER_STAGE_TIMEOUT  лимит секунд на (species,spin)-вариант (default 5400)
  CLER_JSON_SUFFIX    суффикс файла результатов ("_a"/"_b" при сплите)

Запуск:  python3 -u routes/cler_ladder.py            # полный
         python3 -u routes/cler_ladder.py --fast     # смок: single-point, без opt
Резюм:   повторный запуск досчитывает только недостающие (species,spin) —
         каждая сошедшаяся частица немедленно пишется в results-JSON; при
         таймауте стадии сохраняется ПРОМЕЖУТОЧНАЯ геометрия, и следующий
         запуск стартует с неё (важно для тяжёлого тройного центра).
Результат: routes/cler_ladder_results{CLER_JSON_SUFFIX}.json
"""
import os, sys, json, time, signal, argparse

import numpy as np
try:
    from pyscf import gto, dft
    from pyscf.scf import addons
except ImportError as _exc:                       # в сандбоксе pyscf может отсутствовать
    sys.exit("[cler] pyscf не установлен (%s). Скрипт предназначен для "
             "AWS-запуска через routes/cler_aws.py; локально: "
             "pip install numpy scipy pyscf pyberny" % _exc)

HARTREE_EV = 27.211386245988
HERE = os.path.dirname(os.path.abspath(__file__))
SUFFIX = os.environ.get("CLER_JSON_SUFFIX", "")
RESULTS = os.path.join(HERE, f"cler_ladder_results{SUFFIX}.json")
STAGE_TIMEOUT = int(os.environ.get("CLER_STAGE_TIMEOUT", "5400"))  # с на (species,spin)

U_EQ_CL = 1.36                 # В, равновесный потенциал Cl₂/Cl⁻ (станд. условия)
U_EQ_O = 1.23                  # В, равновесный потенциал OER
WATER_SPLIT_EV = 4.92          # 2H₂O → O₂ + 2H₂, эксперимент (computed_here: false)

# --- Литературные ZPE−TS поправки (эВ). computed_here: FALSE ------------------
# Конвенция та же, что в h2_oer_ladder/oer_sulfate_ladder: G ≈ E + (ZPE − TS),
# 298 K, 1 бар, без члена H(298)−H(0). Числа — ЗАГЛУШКИ из литературы, для
# ЭТИХ кластеров колебания НЕ считались.
GCORR_ADS_EV = {"bare": 0.0, "cl": +0.04, "oh": +0.35, "o": +0.05, "ooh": +0.41}
GCORR_REF_EV = {"h2o": -0.11, "h2": -0.14, "cl2": -0.65}
GCORR_SOURCE = {
    "oh/o/ooh/h2o/h2":
        "Nørskov et al. JPC B 108 (2004) 17886; Valdés et al. JPC C 112 (2008) "
        "9872 — стандартный набор CHE (ЗАГЛУШКА, не пересчитана здесь)",
    "cl":
        "Hansen/Man/Rossmeisl/Bligaard/Nørskov, PCCP 12 (2010) 283, "
        "'Electrochemical chlorine evolution at rutile oxide (110) surfaces' — "
        "ZPE−TS для Cl* на RuO2(110) ~+0.04 эВ (ЗАГЛУШКА)",
    "cl2":
        "NIST-JANAF: ν(Cl2)=559.7 см⁻¹ → ZPE=0.035 эВ; S°(Cl2,g)=223.1 "
        "Дж/(моль·К) → TS(298 K)=0.689 эВ; ZPE−TS=−0.654 эВ (ЗАГЛУШКА)",
}

# --- Геометрический генератор рутилового мотива -------------------------------
D_MM = 3.05        # M–M вдоль рёберно-связанной цепочки (RuO2 c ≈ 3.11 Å)
D_OXO = 1.70       # M=O терминальный оксо
D_OH = 1.95        # M–OH гидроксо
D_OH2 = 2.15       # M–OH2 аква
D_OBR = 1.98       # M–O мостиковый
D_MCL = 2.28       # M–Cl

# Активные центры: список металлов цепочки, активный сайт — metals[0].
CENTERS = {
    "ru":       dict(metals=["Ru"],
                     note="чистый RuO2 — промышленный референс DSA"),
    "ir":       dict(metals=["Ir"],
                     note="IrO2 — дорогой стабильный эталон"),
    "pd":       dict(metals=["Pd"],
                     note="'чистый PdO' в рутиловом мотиве — предельный случай "
                          "'максимум своего металла' (формальный Pd(IV))"),
    "ru_pd":    dict(metals=["Ru", "Pd"],
                     note="Ru с соседним Pd — состав Норникеля (Pd вместо части Ir)"),
    "ru_pd_ir": dict(metals=["Ru", "Pd", "Ir"],
                     note="тройной Ru-Pd-Ir — состав покрытия ванны №304"),
}
ADSORBATES = ("bare", "cl", "oh", "o", "ooh")
REFS = ("cl2", "h2o", "h2")


def _ads_atoms(ads):
    """Адсорбат на апикальной позиции M₀ (направление −z), слегка асимметрично."""
    if ads == "bare":
        return []
    if ads == "cl":
        return [("Cl", 0.05, 0.10, -D_MCL)]
    if ads == "o":
        return [("O", 0.04, 0.09, -1.68)]
    if ads == "oh":
        return [("O", 0.05, 0.10, -1.93), ("H", 0.60, 0.75, -2.48)]
    if ads == "ooh":
        return [("O", 0.05, 0.10, -1.96), ("O", 1.05, 0.45, -2.72),
                ("H", 0.85, 1.30, -3.10)]
    raise KeyError(ads)


def build_cluster(metals, ads):
    """Кластер рутилового мотива: цепочка рёберно-связанных узлов вдоль +z.

    metals[0] — активный сайт (адсорбат по −z). Соседние металлы связаны парой
    мостиков μ-O / μ-OH в плоскости xz (стороны чередуются, чтобы не копить
    искусственный диполь). Терминальные лиганды подобраны так, что кластер
    НЕЙТРАЛЕН при формальном M(IV) — см. docstring модуля.
    Возвращает список (symbol, x, y, z) в ангстремах.
    """
    n = len(metals)
    at = []
    for i, m in enumerate(metals):
        at.append((m, 0.0, 0.0, i * D_MM))

    # мостики между соседями
    dx = float(np.sqrt(max(D_OBR ** 2 - (D_MM / 2.0) ** 2, 0.25)))
    for i in range(n - 1):
        zm = (i + 0.5) * D_MM
        s = 1.0 if i % 2 == 0 else -1.0          # чередуем сторону μ-O
        at.append(("O", s * dx, 0.0, zm))                        # μ-O
        at.append(("O", -s * dx, 0.0, zm))                       # μ-OH
        at.append(("H", -s * (dx + 0.62), 0.70, zm + 0.12))

    # терминальные лиганды активного сайта M0
    at.append(("O", 0.0, D_OXO, -0.06))                          # M0=O
    at.append(("O", 0.0, -D_OH, 0.05))                           # M0-OH
    at.append(("H", 0.62, -(D_OH + 0.60), 0.32))
    if n == 1:
        # мостиков нет — свободные ±x достраиваем до нейтрального N=1 узла
        at.append(("O", D_OH, 0.05, 0.10))                       # второй OH
        at.append(("H", D_OH + 0.57, 0.68, 0.55))
        at.append(("O", -D_OH2, 0.05, -0.08))                    # аква
        at.append(("H", -(D_OH2 + 0.45), 0.55, 0.66))
        at.append(("H", -(D_OH2 + 0.47), -0.53, -0.52))

    # терминальные лиганды спектаторных металлов
    for i in range(1, n):
        z = i * D_MM
        if i == n - 1:                                # хвостовой: (=O)(H2O)
            at.append(("O", 0.06, 0.10, z + D_OXO))
            at.append(("O", 0.0, D_OH2, z + 0.05))
            at.append(("H", 0.55, D_OH2 + 0.50, z + 0.55))
            at.append(("H", -0.52, D_OH2 + 0.52, z - 0.45))
        else:                                         # средний: (OH)
            at.append(("O", 0.0, D_OH, z))
            at.append(("H", 0.60, D_OH + 0.58, z + 0.35))
    return at + _ads_atoms(ads)


REF_GEOM = {
    "cl2": [("Cl", 0.0, 0.0, 0.0), ("Cl", 1.99, 0.0, 0.0)],
    "h2o": [("O", 0.0, 0.0, 0.0), ("H", 0.76, 0.59, 0.0), ("H", -0.76, 0.59, 0.0)],
    "h2":  [("H", 0.0, 0.0, 0.0), ("H", 0.74, 0.0, 0.0)],
}
DIATOMIC = {"h2": ("H", "H", 0.62, 0.92), "cl2": ("Cl", "Cl", 1.80, 2.25)}


def species_atoms(species):
    if species in REF_GEOM:
        return REF_GEOM[species]
    center, ads = species.rsplit("_", 1)
    return build_cluster(CENTERS[center]["metals"], ads)


def _ecp_map(atoms, basis):
    """def2-* несёт ECP только для тяжёлых (Ru/Pd: core 28, Ir: core 60).
    Отдаём словарь только по металлам — чтобы pyscf не шумел по O/H/Cl."""
    if not basis.lower().startswith("def2"):
        return None
    heavy = {a[0] for a in atoms if gto.charge(a[0]) > 36}
    return {el: basis for el in heavy} or None


def _mol(species, basis, spin, coords=None):
    """coords (Å) — тёплый старт с сохранённой геометрии (резюм после таймаута)."""
    atoms = species_atoms(species)
    if coords and len(coords) == len(atoms):
        atoms = [(a[0], c[0], c[1], c[2]) for a, c in zip(atoms, coords)]
    return gto.M(atom=[(a[0], (a[1], a[2], a[3])) for a in atoms],
                 basis=basis, ecp=_ecp_map(atoms, basis), charge=0, spin=spin,
                 verbose=0, unit="Angstrom")


def spins_for(species, basis):
    """Спиновый скан с ЧЕСТНОЙ чётностью: число электронов считаем реально
    (после снятия ECP-остова), поэтому 2S = Nα−Nβ подбирается автоматически."""
    if species in REFS:
        return (0,)                               # Cl2/H2O/H2 — closed-shell
    mol = _mol(species, basis, None)              # spin=None → pyscf сам проставит
    odd = mol.nelectron % 2
    ads = species.rsplit("_", 1)[1]
    if odd:
        return (1, 3)
    return (0, 2, 4) if ads == "o" else (2, 0)


def _scf(mol, xc, dm0=None):
    """UKS(DF), Newton первичен; fallback: level-shift DIIS → smearing-DIIS →
    warm-start level-shift DIIS. Незасошедшийся SCF — ошибка, а не тихо
    возвращённая энергия (рецепт h2_oer_ladder/oer_sulfate_ladder)."""
    mf = dft.UKS(mol).density_fit(); mf.xc = xc; mf.verbose = 0
    mf.conv_tol = 1e-8
    mfn = mf.newton()
    try:
        mfn.kernel(dm0=dm0)
    except Exception:
        mfn.converged = False
    if getattr(mfn, "converged", False):
        return mfn
    mf2 = dft.UKS(mol).density_fit(); mf2.xc = xc; mf2.verbose = 0
    mf2.level_shift = 0.3; mf2.max_cycle = 400; mf2.conv_tol = 1e-8
    mf2.kernel(dm0=dm0)
    if mf2.converged:
        return mf2
    mfs = dft.UKS(mol).density_fit(); mfs.xc = xc; mfs.verbose = 0
    mfs = addons.smearing_(mfs, sigma=0.01, method="fermi")
    mfs.max_cycle = 400; mfs.conv_tol = 1e-6
    mfs.kernel(dm0=dm0)
    mf3 = dft.UKS(mol).density_fit(); mf3.xc = xc; mf3.verbose = 0
    mf3.level_shift = 0.3; mf3.max_cycle = 400; mf3.conv_tol = 1e-8
    mf3.kernel(dm0=mfs.make_rdm1())
    if not mf3.converged:
        raise RuntimeError("SCF unconverged after Newton/level-shift/smearing")
    return mf3


GMAX_TOL = 6e-4   # Ha/Bohr — gradient-audit: max-компонента градиента в финале


class StageTimeout(Exception):
    pass


def _alarm(_sig, _frm):
    raise StageTimeout(f"CLER_STAGE_TIMEOUT={STAGE_TIMEOUT}s exceeded")


def _opt_budget(natm):
    """Бюджет berny масштабируем по размеру кластера. Для тяжёлых центров берём
    КОРОТКИЕ прогоны с бОльшим числом рестартов: между прогонами геометрия
    чекпойнтится в results-JSON, поэтому STAGE_TIMEOUT режет максимум один
    прогон, а не всю стадию. Осознанный размен — числа уходят в JSON."""
    if natm <= 12:
        return 60, 2      # до 180 шагов, чекпойнт каждые 60
    if natm <= 16:
        return 25, 3      # до 100 шагов, чекпойнт каждые 25
    return 15, 5          # до 90 шагов, чекпойнт каждые 15


def _opt_diatomic(species, basis, xc):
    """berny не любит двухатомники — 1D-скан + доводка (рецепт h2_oer_ladder)."""
    a, b, lo, hi = DIATOMIC[species]
    ecp = _ecp_map([(a, 0, 0, 0)], basis)

    def e_at(d):
        m = gto.M(atom=[(a, (0, 0, 0)), (b, (float(d), 0, 0))], basis=basis,
                  ecp=ecp, spin=0, verbose=0, unit="Angstrom")
        return float(_scf(m, xc).e_tot)
    step = (hi - lo) / 15.0
    coarse = {round(float(d), 3): e_at(d) for d in np.arange(lo, hi, step)}
    d0 = min(coarse, key=coarse.get)
    fine = {round(float(d), 4): e_at(d)
            for d in np.arange(d0 - step, d0 + step * 1.01, step / 5.0)}
    d1 = min(fine, key=fine.get)
    return d1, fine[d1]


def optimize_species(species, spin, basis, xc, do_opt=True, start_coords=None,
                     progress=None):
    """Оптимизация геометрии (pyberny) для (species, spin). Сходимость
    проверяется ПРЯМО — нормой градиента в финальной точке (gradient-audit:
    gmax уходит в JSON).

    start_coords — тёплый старт с сохранённой геометрии; progress — mutable
    словарь, куда после КАЖДОГО berny-прогона кладётся текущая геометрия, чтобы
    таймаут стадии не обнулял работу (вызывающий пишет его в results-JSON)."""
    t0 = time.time()
    if do_opt and species in DIATOMIC:
        d1, e1 = _opt_diatomic(species, basis, xc)
        return dict(species=species, spin=spin, e_ha=e1, opt=True, gmax=0.0,
                    ss=0.0, coords=[[0, 0, 0], [d1, 0, 0]], geom_opt=True,
                    natm=2, maxsteps=None, restarts=None,
                    seconds=round(time.time() - t0, 1))
    mol = _mol(species, basis, spin, coords=start_coords)
    natm = mol.natm
    maxsteps, restarts = _opt_budget(natm)
    mf = _scf(mol, xc)
    e0 = float(mf.e_tot)
    coords = mol.atom_coords(unit="Angstrom").tolist()
    converged_opt, gmax = None, None
    if do_opt:
        from pyscf.geomopt.berny_solver import optimize as berny_opt
        try:
            cur_mf = mf
            for _attempt in range(restarts + 1):
                mol_eq = berny_opt(cur_mf, maxsteps=maxsteps, verbose=0)
                atoms = list(zip([mol.atom_symbol(i) for i in range(mol.natm)],
                                 mol_eq.atom_coords(unit="Angstrom")))
                if progress is not None:     # чекпойнт геометрии ДО дорогого SCF
                    progress["coords"] = mol_eq.atom_coords(
                        unit="Angstrom").tolist()
                    progress["berny_runs"] = _attempt + 1
                cur_mf = _scf(gto.M(atom=atoms, basis=basis,
                                    ecp=_ecp_map([(s, 0, 0, 0) for s, _c in atoms],
                                                 basis),
                                    charge=0, spin=spin, verbose=0,
                                    unit="Angstrom"), xc)
                gmax = float(np.abs(cur_mf.nuc_grad_method().kernel()).max())
                if progress is not None:
                    progress["gmax"] = gmax
                if gmax < GMAX_TOL:
                    break
            mf = cur_mf
            coords = mol_eq.atom_coords(unit="Angstrom").tolist()
            converged_opt = True if gmax < GMAX_TOL else \
                f"UNCONVERGED gmax={gmax:.5f} after {restarts+1} berny runs"
        except StageTimeout:
            raise                                # таймаут фиксируем выше, честно
        except Exception as exc:                 # неудачу оптимизации фиксируем
            converged_opt = f"FAILED: {type(exc).__name__}: {exc}"
    e_final = float(mf.e_tot)
    ss = float(mf.spin_square()[0])
    failed = isinstance(converged_opt, str) and converged_opt.startswith("FAILED")
    return dict(species=species, spin=spin,
                e_ha=e0 if failed else e_final,
                opt=converged_opt, gmax=gmax, ss=round(ss, 3), coords=coords,
                geom_opt=bool(do_opt), natm=natm, warm_start=bool(start_coords),
                maxsteps=maxsteps if do_opt else None,
                restarts=restarts if do_opt else None,
                seconds=round(time.time() - t0, 1))


def load_results():
    if os.path.exists(RESULTS):
        with open(RESULTS) as f:
            return json.load(f)
    return {"meta": {}, "runs": {}}


def save_results(data):
    tmp = RESULTS + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=1, ensure_ascii=False)
    os.replace(tmp, RESULTS)


def _gcorr(name):
    """ZPE−TS поправка (эВ): ru_pd_cl → cl, ir_bare → bare, cl2 → cl2."""
    if name in GCORR_REF_EV:
        return GCORR_REF_EV[name]
    return GCORR_ADS_EV[name.rsplit("_", 1)[1]]


def ladder(best, center):
    """ClER-лестница центра + конкурирующая OER-лестница на ТОМ ЖЕ центре."""
    need_cl = [f"{center}_bare", f"{center}_cl", "cl2"]
    if any(k not in best for k in need_cl):
        return None
    g = {}
    for k in [f"{center}_{a}" for a in ADSORBATES] + list(REFS):
        if k in best:
            g[k] = best[k]["e_ha"] * HARTREE_EV + _gcorr(k)
    if any(k not in g for k in need_cl):
        return None

    dg_cl = g[f"{center}_cl"] - g[f"{center}_bare"] - 0.5 * g["cl2"]
    dg1 = dg_cl + U_EQ_CL
    dg2 = U_EQ_CL - dg_cl
    out = {
        "center": center,
        "metals": CENTERS[center]["metals"],
        "note": CENTERS[center]["note"],
        "dG_Cl_star_eV": round(dg_cl, 3),
        "cler_steps_eV": {"dG1_bare_to_Cl": round(dg1, 3),
                          "dG2_Cl_to_Cl2": round(dg2, 3)},
        "cler_limiting": "dG1_bare_to_Cl" if dg1 >= dg2 else "dG2_Cl_to_Cl2",
        "eta_ClER_V": round(max(dg1, dg2) - U_EQ_CL, 3),
        "cler_sum_check_eV": round(dg1 + dg2, 3),      # обязан быть 2*1.36
        "best_spins": {k: best[k]["spin"] for k in
                       [f"{center}_{a}" for a in ADSORBATES] + list(REFS)
                       if k in best},
    }

    # --- конкурирующий OER на том же центре ---
    have_h = "h2o" in g and "h2" in g
    if have_h and f"{center}_oh" in g:
        out["dG_OH_star_eV"] = round(
            g[f"{center}_oh"] - g[f"{center}_bare"] - g["h2o"] + 0.5 * g["h2"], 3)
    if have_h and f"{center}_o" in g:
        dg_o = g[f"{center}_o"] - g[f"{center}_bare"] - g["h2o"] + g["h2"]
        out["dG_O_star_eV"] = round(dg_o, 3)
        # ГЛАВНОЕ ЧИСЛО: селективность Cl2 против O2 на одном и том же центре
        out["dSel_Cl_minus_O_eV"] = round(dg_cl - dg_o, 3)
        # смещение относительно скейлинга ClER↔OER (ΔG_Cl ≈ ½ΔG_O + const)
        out["scaling_offset_eV"] = round(dg_cl - 0.5 * dg_o, 3)
    if have_h and f"{center}_ooh" in g:
        out["dG_OOH_star_eV"] = round(
            g[f"{center}_ooh"] - g[f"{center}_bare"] - 2 * g["h2o"]
            + 1.5 * g["h2"], 3)
    if all(k in out for k in ("dG_OH_star_eV", "dG_O_star_eV", "dG_OOH_star_eV")):
        s1 = out["dG_OH_star_eV"]
        s2 = out["dG_O_star_eV"] - out["dG_OH_star_eV"]
        s3 = out["dG_OOH_star_eV"] - out["dG_O_star_eV"]
        s4 = WATER_SPLIT_EV - out["dG_OOH_star_eV"]
        steps = {"dG1_bare_to_OH": s1, "dG2_OH_to_O": s2,
                 "dG3_O_to_OOH": s3, "dG4_OOH_to_O2": s4}
        lim = max(steps, key=steps.get)
        out["oer_steps_eV"] = {k: round(v, 3) for k, v in steps.items()}
        out["oer_limiting"] = lim
        out["eta_OER_V"] = round(steps[lim] - U_EQ_O, 3)
        out["eta_OER_minus_ClER_V"] = round(out["eta_OER_V"]
                                            - out["eta_ClER_V"], 3)
    return out


def verdict(lads):
    """Ранжир по η_ClER и по dSel + прямой ответ на вопрос «сколько Ir убрать»."""
    v = {}
    by_eta = sorted((l for l in lads.values() if "eta_ClER_V" in l),
                    key=lambda l: l["eta_ClER_V"])
    v["rank_by_eta_ClER"] = [{"center": l["center"], "eta_ClER_V": l["eta_ClER_V"],
                              "dG_Cl_star_eV": l["dG_Cl_star_eV"]} for l in by_eta]
    by_sel = sorted((l for l in lads.values() if "dSel_Cl_minus_O_eV" in l),
                    key=lambda l: l["dSel_Cl_minus_O_eV"])
    v["rank_by_dSel_more_negative_is_more_Cl2_selective"] = [
        {"center": l["center"], "dSel_Cl_minus_O_eV": l["dSel_Cl_minus_O_eV"],
         "dG_Cl_star_eV": l["dG_Cl_star_eV"],
         "dG_O_star_eV": l["dG_O_star_eV"]} for l in by_sel]

    def _get(c, k):
        return lads.get(c, {}).get(k)

    ir_cut = {}
    for a, b, label in (("ru_pd_ir", "ru_pd", "убрать Ir из тройного (№304 → Ru-Pd)"),
                        ("ru_pd_ir", "ru", "убрать и Ir, и Pd (№304 → чистый RuO2)"),
                        ("ru_pd", "pd", "убрать Ru (Ru-Pd → чистый PdO, предел)"),
                        ("ir", "ru", "IrO2 → RuO2 (эталонная замена)")):
        ea, eb = _get(a, "eta_ClER_V"), _get(b, "eta_ClER_V")
        sa, sb = _get(a, "dSel_Cl_minus_O_eV"), _get(b, "dSel_Cl_minus_O_eV")
        if ea is None or eb is None:
            continue
        d_eta = eb - ea
        ir_cut[f"{b}_minus_{a}"] = {
            "what": label,
            "d_eta_ClER_V": round(d_eta, 3),
            "d_dSel_eV": (None if sa is None or sb is None else round(sb - sa, 3)),
            "verdict": ("выигрыш: упрощённый состав даже лучше по η"
                        if d_eta <= -0.05 else
                        "нейтрально: замена без потери (|Δη| < 0.05 В)"
                        if d_eta < 0.05 else
                        "потеря η — убранный металл нужен"),
        }
    v["ir_reduction"] = ir_cut
    if by_eta:
        v["best_by_eta"] = by_eta[0]["center"]
    if by_sel:
        v["best_by_dSel"] = by_sel[0]["center"]
    v["how_to_read"] = (
        "η_ClER = |ΔG_Cl*| (вулкан ClER, оптимум ΔG_Cl*=0). dSel = ΔG_Cl* − ΔG_O*: "
        "чем отрицательнее, тем сильнее хлорный путь выигрывает у паразитного O2 "
        "на ТОМ ЖЕ центре. Оба числа — КЛАСТЕРНЫЕ: сравнивать между центрами, "
        "не с экспериментальными вольтами.")
    return v


def _model_sizes(centers, basis):
    """Размеры кластеров (атомы/электроны после ECP/AO) — в JSON, чтобы читатель
    видел, НАСКОЛЬКО минимальна модель, и мог оценить бюджет счёта."""
    out = {}
    for c in centers:
        try:
            m = _mol(f"{c}_bare", basis, None)
            mc = _mol(f"{c}_cl", basis, None)
            out[c] = {"natm_bare": m.natm, "nelec_bare": m.nelectron,
                      "nao_bare": int(m.nao), "nao_Cl_star": int(mc.nao),
                      "spins_bare": list(spins_for(f"{c}_bare", basis))}
        except Exception as exc:                  # диагностика не должна ронять счёт
            out[c] = {"error": f"{type(exc).__name__}: {exc}"}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true",
                    help="смок: single-point, без оптимизации геометрии "
                         "(базис остаётся def2-SVP: sto-3g не имеет Ir/ECP)")
    ap.add_argument("--basis", default=None)
    ap.add_argument("--centers", default=None,
                    help="подмножество центров через запятую (приоритет над "
                         "CLER_CENTERS)")
    args = ap.parse_args()
    basis = args.basis or os.environ.get("CLER_BASIS") or "def2-svp"
    xc = os.environ.get("CLER_XC", "pbe,pbe")
    do_opt = not args.fast

    sel = args.centers or os.environ.get("CLER_CENTERS") or ""
    centers = [c.strip() for c in sel.split(",") if c.strip()] or list(CENTERS)
    unknown = [c for c in centers if c not in CENTERS]
    if unknown:
        sys.exit(f"[cler] неизвестные центры: {unknown}; "
                 f"допустимые: {sorted(CENTERS)}")
    # референсы Cl2/H2O/H2 нужны КАЖДОЙ ветке сплита — считаем всегда
    todo = [f"{c}_{a}" for c in centers for a in ADSORBATES] + list(REFS)

    data = load_results()
    data["meta"] = {
        "purpose": "N1-v2: расчётная поддержка ЖИВОЙ анодной программы Ru-Pd-Ir "
                   "Кольской ГМК (хлорный анолит, ванна №304, апр–окт 2026, "
                   "~120 кт Ni/год, ~500 ванн); цель — МИНИМУМ Ir при удержании "
                   "перенапряжения ClER и селективности Cl2 против O2",
        "case": "N1-v2 хлорный анод электроэкстракции Ni (не OER — анолит "
                "хлорный); контекст APPLIED_PRIOR_ART_SCAN.md разд. 3",
        "centers": {c: CENTERS[c] for c in centers},
        "method": f"UKS-{xc}(DF)/{basis} + def2-ECP (Ru/Pd core 28, Ir core 60)"
                  + (" + berny geom-opt" if do_opt else " single-point"),
        "descriptors": {
            "dG_Cl_star_eV": "G(Cl*) − G(*) − 1/2 G(Cl2) — дескриптор вулкана ClER",
            "dG1": "ΔG_Cl* + 1.36 эВ  (* + Cl- -> Cl* + e-)",
            "dG2": "1.36 эВ − ΔG_Cl*  (Cl* + Cl- -> Cl2 + e-)",
            "eta_ClER_V": "max(dG1,dG2)/e − 1.36 В = |ΔG_Cl*|",
            "closure": "dG1 + dG2 = 2·1.36 эВ тождественно (встроено в референс "
                       "Cl2, не подгонка)",
            "dG_O_star_eV": "G(O*) − G(*) − G(H2O) + G(H2)",
            "dG_OH_star_eV": "G(OH*) − G(*) − G(H2O) + 1/2 G(H2)",
            "dSel_Cl_minus_O_eV": "ΔG_Cl* − ΔG_O* — ГЛАВНОЕ ЧИСЛО: селективность "
                                  "Cl2 против паразитного O2 на одном центре; "
                                  "отрицательнее = хлор выигрывает",
            "eta_OER_V": "max(4 шага CHE)/e − 1.23 В — паразитный кислород",
        },
        "model_sizes": _model_sizes(centers, basis),
        "gcorr_eV": {**GCORR_REF_EV, **GCORR_ADS_EV},
        "gcorr_source": GCORR_SOURCE,
        "computed_here": {"electronic_energies": True, "geometries": True,
                          "gmax_gradient_audit": True,
                          "spin_scan": True,
                          "zpe_ts_gcorr": False,        # лит. заглушки, см. gcorr_source
                          "cl2_ref_energy": True,       # Cl2 считан явно (closed-shell)
                          "u_eq_1.36V": False,          # табличный станд. потенциал
                          "water_split_4.92eV": False}, # эксперимент
        "stage_timeout_s": STAGE_TIMEOUT,
        "honesty": "Молекулярный кластер рутилового мотива != рутиловая "
                   "поверхность (110): нет решётки, нет соседей по покрытию, нет "
                   "растворителя, двойного слоя и приложенного потенциала. "
                   "ZPE-TS поправки ЛИТЕРАТУРНЫЕ (computed_here: false), для этих "
                   "кластеров колебания не считались. PBE переоценивает связь Cl2 "
                   "и занижает d-зазоры 4d/5d; Pd(IV) в рутиловом мотиве — "
                   "формальная конструкция. Это ЗНАКОВЫЙ СКРИН РАНЖИРА центров, "
                   "а НЕ абсолютные перенапряжения: сравнивать центры между "
                   "собой, не с вольтами ванны №304.",
    }
    save_results(data)

    have_alarm = hasattr(signal, "SIGALRM")
    if have_alarm:
        signal.signal(signal.SIGALRM, _alarm)
    for sp in todo:
        for spin in spins_for(sp, basis):
            key = f"{basis}:{xc}:{sp}:s{spin}" + ("" if do_opt else ":sp")
            prev = data["runs"].get(key, {})
            # резюм: пропускаем только УСПЕШНЫЕ записи (ошибки пересчитываем)
            if "e_ha" in prev and prev.get("opt") in (True, None):
                print(f"[skip] {key} (готово)"); continue
            # тёплый старт: если прошлая попытка упала по таймауту/не сошлась,
            # но геометрия успела прочекпойнтиться — продолжаем с неё
            start = prev.get("coords") if isinstance(prev.get("coords"), list) \
                else None
            print(f"[run ] {key} (timeout {STAGE_TIMEOUT}s"
                  + (", warm start" if start else "") + ") …", flush=True)
            progress = {}
            if have_alarm:
                signal.alarm(STAGE_TIMEOUT)
            try:
                res = optimize_species(sp, spin, basis, xc, do_opt=do_opt,
                                       start_coords=start, progress=progress)
            except StageTimeout as exc:
                res = dict(species=sp, spin=spin, timeout=True, error=str(exc),
                           warm_start=bool(start), **progress,
                           resume_hint="геометрия прочекпойнчена — повторный "
                                       "запуск продолжит с неё, не с нуля")
            except Exception as exc:
                res = dict(species=sp, spin=spin,
                           error=f"{type(exc).__name__}: {exc}",
                           warm_start=bool(start), **progress)
            finally:
                if have_alarm:
                    signal.alarm(0)
            # поточечный чекпойнт: каждая частица пишется НЕМЕДЛЕННО
            data["runs"][key] = res
            save_results(data)
            print(f"       E={res.get('e_ha', float('nan')):.6f} Ha  "
                  f"opt={res.get('opt')}  gmax={res.get('gmax')}  "
                  f"⟨S²⟩={res.get('ss')}  natm={res.get('natm')}  "
                  f"{res.get('seconds', '?')}s", flush=True)

    # лучшие (низшие) по спину — только успешные И того же режима (opt vs :sp)
    best = {}
    for key, r in data["runs"].items():
        parts = key.split(":")
        if len(parts) < 4:
            continue
        b, kxc, sp, is_sp = parts[0], parts[1], parts[2], key.endswith(":sp")
        if b != basis or kxc != xc or "e_ha" not in r or is_sp == do_opt:
            continue
        if do_opt and r.get("opt") is not True:   # UNCONVERGED/FAILED — вон
            continue
        if sp not in best or r["e_ha"] < best[sp]["e_ha"]:
            best[sp] = r

    lads = {}
    for c in centers:
        lad = ladder(best, c)
        if lad is None:
            missing = [k for k in (f"{c}_bare", f"{c}_cl", "cl2") if k not in best]
            print(f"[warn] центр {c}: ClER-лестница не собрана, нет {missing}")
            continue
        lads[c] = lad
        print(json.dumps(lad, indent=1, ensure_ascii=False), flush=True)
    if lads:
        entry = {"basis": basis, "xc": xc, "geom_opt": do_opt,
                 "centers": lads, "verdict": verdict(lads)}
        # изоляция режимов: --fast НЕ трогает продакшен-слот "ladder"
        data.setdefault("ladders", {})[f"{basis}:{xc}"
                                       + ("" if do_opt else ":sp")] = entry
        if do_opt:
            data["ladder"] = entry
        save_results(data)
        print(json.dumps(entry["verdict"], indent=1, ensure_ascii=False),
              flush=True)


if __name__ == "__main__":
    main()
