#!/usr/bin/env python3
"""
RadShield screen — скрининг материалов радиационной защиты скафандра.

Что считает (всё — из первых принципов / арифметики состава, без подгонки):
  1. Водородная плотность: масс. доля H, моль H на грамм.
  2. Электронная плотность Z/A (главный фактор тормозной способности на грамм).
  3. Эффективный потенциал ионизации I_eff — брэгговская аддитивность
     (ICRU-37: элементные I, для H/C/N/O — рекомендованные "in compounds").
  4. Массовая тормозная способность протонов, формула Бете (реляционная,
     полный T_max, БЕЗ shell/density поправок — оговорка ниже).
  5. CSDA-пробег протонов (численное интегрирование 1/S от 0.5 МэВ).
  6. Энергия отсечки протонов для слоёв 0.5 / 1 / 2 г/см²
     (какой протон слой ещё останавливает).
  7. Захват тепловых нейтронов на ¹⁰B: атомы ¹⁰B на грамм, вероятность
     захвата на слое (после термализации! — см. оговорку).

Валидация (встроенная, роняет скрипт при расхождении):
  - вода @100 МэВ против NIST PSTAR ≈ 7.28–7.29 МэВ·см²/г
    (значение сверено WebSearch 2026-08-25, physics.nist.gov PSTAR);
  - I_eff(полиэтилен) против табличного ICRU-37 = 57.4 эВ.

Честные оговорки:
  - Без shell-поправки: ниже ~10 МэВ ошибка S растёт до нескольких %,
    ниже ~2 МэВ формула Бете ломается; вклад хвоста <2 МэВ в пробег
    при E>=10 МэВ — <1%, поэтому пробеги/отсечки честны в ~2-3%.
  - Тормозная способность — это SEP-протоны и захваченные пояса.
    ГКЛ (ГэВ/нуклон) слоем скафандра НЕ останавливаются — этот скрипт
    не считает дозовый транспорт (нужен HZETRN/Geant4, вне тулчейна).
  - Захват ¹⁰B применим к УЖЕ термализованным нейтронам (альбедо от тела,
    отражённые от корпуса); транспорт быстрых нейтронов не считается.

Запуск:  python3 calc/radshield_screen.py
Выход:   calc/radshield_screen_results.json
"""

import json
import math
import os

# ---------- константы (PDG) ----------
K = 0.307075          # МэВ·см²/моль (4π N_A r_e² m_e c²)
ME = 0.5109989        # МэВ, m_e c²
MP = 938.2721         # МэВ, m_p c²
NA = 6.02214076e23

# атомные массы (г/моль) и заряды
ELEMS = {
    "H": (1.008, 1), "Li": (6.94, 3), "B": (10.811, 5), "C": (12.011, 6),
    "N": (14.007, 7), "O": (15.999, 8), "Mg": (24.305, 12), "Al": (26.982, 13),
}

# I, эВ. ICRU-37: H,C,N,O — рекомендованные для конденсированных соединений;
# Li, B, Mg, Al — элементные табличные.
I_EV = {"H": 19.2, "Li": 40.0, "B": 76.0, "C": 81.0, "N": 82.0,
        "O": 106.0, "Mg": 156.0, "Al": 166.0}

B10_ABUNDANCE = 0.199        # доля ¹⁰B в природном боре (сверено 2026-08-25)
B10_SIGMA_TH = 3837e-24      # см², тепловой захват ¹⁰B (сверено 2026-08-25)

# ---------- материалы: стехиометрия {элемент: число атомов} ----------
MATERIALS = {
    "Al (корпус, базовая линия)":      {"Al": 1},
    "Вода H2O":                        {"H": 2, "O": 1},
    "Полиэтилен (CH2)n / СВМПЭ":       {"C": 1, "H": 2},
    "Кевлар (C14H10N2O2)n":            {"C": 14, "H": 10, "N": 2, "O": 2},
    "LiH":                             {"Li": 1, "H": 1},
    "MgH2":                            {"Mg": 1, "H": 2},
    "LiBH4":                           {"Li": 1, "B": 1, "H": 4},
    "Боразан NH3BH3 (AB)":             {"N": 1, "B": 1, "H": 6},
    "BN (нанотрубки BNNT)":            {"B": 1, "N": 1},
}

# композиты: весовые доли компонент из MATERIALS
COMPOSITES = {
    "ГБК-50: 50% AB + 50% СВМПЭ":
        [("Боразан NH3BH3 (AB)", 0.5), ("Полиэтилен (CH2)n / СВМПЭ", 0.5)],
    "ГБК-30: 30% AB + 70% СВМПЭ":
        [("Боразан NH3BH3 (AB)", 0.3), ("Полиэтилен (CH2)n / СВМПЭ", 0.7)],
}

E_REPORT_S = [10.0, 30.0, 100.0, 250.0]      # МэВ, точки отчёта S
E_REPORT_R = [10.0, 30.0, 100.0]             # МэВ, точки отчёта пробега
DEPTHS = [0.5, 1.0, 2.0]                     # г/см², слои для отсечки


def weight_fractions(stoich):
    m = {el: n * ELEMS[el][0] for el, n in stoich.items()}
    tot = sum(m.values())
    return {el: v / tot for el, v in m.items()}, tot


def mix_fractions(parts):
    w = {}
    for name, frac in parts:
        wf, _ = weight_fractions(MATERIALS[name])
        for el, v in wf.items():
            w[el] = w.get(el, 0.0) + frac * v
    return w


def z_over_a(wf):
    return sum(w * ELEMS[el][1] / ELEMS[el][0] for el, w in wf.items())


def i_eff(wf):
    num = sum(w * ELEMS[el][1] / ELEMS[el][0] * math.log(I_EV[el])
              for el, w in wf.items())
    return math.exp(num / z_over_a(wf))


def bethe_S(E, za, i_ev):
    """Массовая тормозная способность протона, МэВ·см²/г."""
    g = 1.0 + E / MP
    b2 = 1.0 - 1.0 / g**2
    tmax = 2 * ME * b2 * g**2 / (1 + 2 * g * ME / MP + (ME / MP) ** 2)
    i_mev = i_ev * 1e-6
    arg = 2 * ME * b2 * g**2 * tmax / i_mev**2
    return K * za / b2 * (0.5 * math.log(arg) - b2)


def build_range_table(za, i_ev, e_lo=0.5, e_hi=400.0, n=4000):
    """CSDA-пробег: трапеции по 1/S на лог-сетке. Хвост <0.5 МэВ отброшен."""
    es = [e_lo * (e_hi / e_lo) ** (i / (n - 1)) for i in range(n)]
    inv = [1.0 / bethe_S(e, za, i_ev) for e in es]
    r, rs = 0.0, [0.0]
    for i in range(1, n):
        r += 0.5 * (inv[i] + inv[i - 1]) * (es[i] - es[i - 1])
        rs.append(r)
    return es, rs


def interp(xs, ys, x):
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    t = (x - xs[lo]) / (xs[hi] - xs[lo])
    return ys[lo] + t * (ys[hi] - ys[lo])


def analyze(name, wf):
    za = z_over_a(wf)
    ie = i_eff(wf)
    h_w = wf.get("H", 0.0)
    b_w = wf.get("B", 0.0)
    es, rs = build_range_table(za, ie)
    b10_per_g = b_w / ELEMS["B"][0] * NA * B10_ABUNDANCE
    sigma_m = b10_per_g * B10_SIGMA_TH  # см²/г
    return {
        "material": name,
        "weight_fractions": {el: round(v, 5) for el, v in sorted(wf.items())},
        "H_wt_pct": round(100 * h_w, 2),
        "H_mol_per_g": round(h_w / ELEMS["H"][0], 5),
        "B_wt_pct": round(100 * b_w, 2),
        "Z_over_A": round(za, 5),
        "I_eff_eV": round(ie, 1),
        "S_p_MeV_cm2_g": {f"{e:g}MeV": round(bethe_S(e, za, ie), 3)
                          for e in E_REPORT_S},
        "csda_range_g_cm2": {f"{e:g}MeV": round(interp(es, rs, e), 4)
                             for e in E_REPORT_R},
        "proton_cutoff_MeV": {f"{d:g}g_cm2": round(interp(rs, es, d), 2)
                              for d in DEPTHS},
        "B10_atoms_per_g": f"{b10_per_g:.3e}",
        "thermal_n_capture_prob_1g_cm2":
            round(1.0 - math.exp(-sigma_m * 1.0), 4),
    }


def main():
    results = []
    for name, st in MATERIALS.items():
        wf, _ = weight_fractions(st)
        results.append(analyze(name, wf))
    for name, parts in COMPOSITES.items():
        results.append(analyze(name, mix_fractions(parts)))

    # ---- встроенная валидация ----
    wf_w, _ = weight_fractions(MATERIALS["Вода H2O"])
    s_water_100 = bethe_S(100.0, z_over_a(wf_w), i_eff(wf_w))
    wf_pe, _ = weight_fractions(MATERIALS["Полиэтилен (CH2)n / СВМПЭ"])
    i_pe = i_eff(wf_pe)
    validation = {
        "water_S_100MeV_computed": round(s_water_100, 3),
        "water_S_100MeV_PSTAR_anchor": "7.28-7.29 (NIST PSTAR, сверено WebSearch 2026-08-25)",
        "water_anchor_dev_pct": round(100 * (s_water_100 / 7.289 - 1), 2),
        "PE_I_eff_computed_eV": round(i_pe, 1),
        "PE_I_eff_ICRU37_eV": 57.4,
        "PE_I_dev_pct": round(100 * (i_pe / 57.4 - 1), 2),
        "note": ("Bethe без shell/density поправок; I по брэгговской "
                 "аддитивности. Дозовый транспорт (ГКЛ, нейтроны) НЕ "
                 "считается — нужен HZETRN/Geant4."),
    }
    assert abs(s_water_100 / 7.289 - 1) < 0.03, "PSTAR-якорь разошёлся >3%"
    assert abs(i_pe / 57.4 - 1) < 0.05, "I_eff(PE) разошёлся с ICRU-37 >5%"

    out = {
        "what": "RadShield screen: H-плотность, Z/A, Бете-тормозная способность "
                "протонов, CSDA-пробеги, отсечки по слоям, захват тепловых n на 10B",
        "date": "2026-08-25",
        "method": "состав-арифметика + Bethe (PDG, full Tmax, no shell/density) "
                  "+ Bragg additivity (ICRU-37) + sigma_th(10B)=3837 b",
        "validation": validation,
        "results": results,
    }
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "radshield_screen_results.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print(f"OK -> {path}")
    print(f"валидация: S(вода,100МэВ)={s_water_100:.3f} "
          f"(PSTAR 7.289, откл. {validation['water_anchor_dev_pct']}%); "
          f"I_eff(PE)={i_pe:.1f} эВ (ICRU 57.4, откл. {validation['PE_I_dev_pct']}%)")
    hdr = f"{'материал':38s} {'H,%':>6s} {'Z/A':>6s} {'S@100':>6s} {'Ec@1':>6s}"
    print(hdr)
    for r in results:
        print(f"{r['material']:38s} {r['H_wt_pct']:6.2f} {r['Z_over_A']:6.4f} "
              f"{r['S_p_MeV_cm2_g']['100MeV']:6.3f} "
              f"{r['proton_cutoff_MeV']['1g_cm2']:6.2f}")


if __name__ == "__main__":
    main()
