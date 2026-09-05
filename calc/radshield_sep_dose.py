#!/usr/bin/env python3
"""
RadShield этап 4-lite — доза кожи за слоем скафандра в солнечном протонном
событии (SEP), straight-ahead приближение. Отвечает на вопрос: «во сколько раз
дольше можно безопасно оставаться в EVA во время события в ГБК-костюме?»

Метод (первый порядок, честные границы применимости):
  - Спектр: worst-case событие августа-1972 (King 1974, AL-событие):
      J(>E) = J0·exp[(30−E)/E0], J0 = 7.9e9 см⁻², E0 = 26.5 МэВ
    (формула и параметры сверены WebSearch 2026-08-30, SPENVIS help/King).
    Плюс скан жёсткости E0 = 10 / 26.5 / 50 МэВ (мягкое/1972/жёсткое) —
    для отношений доз нормировка не важна.
  - Транспорт: CSDA + straight-ahead (протон идёт по прямой, замедляется по
    Бете из radshield_screen, валидированного по NIST PSTAR −0.03%).
    Пренебрегается: ядерные взаимодействия (для E<100 МэВ на 1–2 г/см² — на
    уровне ~1–2%), рассеяние, вторичные нейтроны, геометрия тела (плоский слой).
  - Доза: энерговыделение в первом 0.1 г/см² (1 мм) воды за слоем — «shallow
    skin dose». Кожа — лимитирующий орган EVA при SEP.
  - Ниже 1 МэВ остаточной энергии тормозная способность клампится на S(1 МэВ)
    (Бете там ломается); занижение хвоста Брэгг-пика одинаково для всех
    материалов — на отношения почти не влияет.

Что НЕ утверждается: абсолютные греи здесь — первый порядок, с BRYNTRN/Geant4
не сверены (первоисточник-таблицы недоступны из этого окружения); выводить в
дневник можно ТОЛЬКО отношения доз / коэффициенты времени.

Запуск:  python3 calc/radshield_sep_dose.py
Выход:   calc/radshield_sep_dose_results.json
"""

import json
import math
import os

import radshield_screen as rs

J0 = 7.9e9        # см⁻², King AL-1972
E30 = 30.0        # МэВ, опорная энергия формулы King
E0_1972 = 26.5    # МэВ
E0_SCAN = [10.0, 26.5, 50.0]

E_MAX = 400.0     # МэВ, верх интегрирования (вклад >400 МэВ в дозу кожи мал)
SKIN = 0.1        # г/см², слой воды для «shallow dose»
DEPTHS = [0.3, 1.0, 2.0]   # г/см² защиты
MEV_G_TO_GY = 1.602176634e-10  # 1 МэВ/г = 1.602e-10 Гр

COMPARE = [
    "Al (корпус, базовая линия)",
    "Полиэтилен (CH2)n / СВМПЭ",
    "ГБК-50: 50% AB + 50% СВМПЭ",
    "Боразан NH3BH3 (AB)",
]


def material_wf(name):
    if name in rs.MATERIALS:
        return rs.weight_fractions(rs.MATERIALS[name])[0]
    return rs.mix_fractions(rs.COMPOSITES[name])


class Medium:
    def __init__(self, name):
        self.name = name
        wf = material_wf(name)
        self.za = rs.z_over_a(wf)
        self.ie = rs.i_eff(wf)
        self.es, self.rt = rs.build_range_table(self.za, self.ie,
                                                e_lo=0.5, e_hi=E_MAX, n=6000)

    def S(self, e):
        return rs.bethe_S(max(e, 1.0), self.za, self.ie)  # кламп <1 МэВ

    def rng(self, e):
        return rs.interp(self.es, self.rt, e)

    def energy_at_range(self, r):
        if r <= 0:
            return 0.0
        return rs.interp(self.rt, self.es, r)


WATER = Medium("Вода H2O")


def skin_dep(e_out):
    """Энерговыделение (МэВ) протона с энергией e_out в 0.1 г/см² воды."""
    if e_out <= 0:
        return 0.0
    r = WATER.rng(e_out)
    if r <= SKIN:
        return e_out
    return e_out - WATER.energy_at_range(r - SKIN)


def dose_behind(med, depth, e0):
    """Доза (Гр) в слое кожи за depth г/см² материала med, спектр King(E0)."""
    ec = med.energy_at_range(depth) if depth > 0 else 0.5
    if ec >= E_MAX:
        return 0.0, ec
    n = 4000
    # лог-сетка по u = E - Ec: сгущение у отсечки, где Брэгг-пик
    u_lo, u_hi = 1e-3, E_MAX - ec
    us = [u_lo * (u_hi / u_lo) ** (i / (n - 1)) for i in range(n)]
    total = 0.0
    prev_u, prev_f = None, None
    for u in us:
        e = ec + u
        phi = (J0 / e0) * math.exp((E30 - e) / e0)          # p/(см²·МэВ)
        e_out = med.energy_at_range(med.rng(e) - depth)
        f = phi * skin_dep(e_out)
        if prev_u is not None:
            total += 0.5 * (f + prev_f) * (u - prev_u)
        prev_u, prev_f = u, f
    return total / SKIN * MEV_G_TO_GY, ec


def main():
    media = [Medium(name) for name in COMPARE]
    out_rows = []
    for depth in DEPTHS:
        doses = {}
        for m in media:
            d, ec = dose_behind(m, depth, E0_1972)
            doses[m.name] = (d, ec)
        d_al = doses[COMPARE[0]][0]
        for m in media:
            d, ec = doses[m.name]
            out_rows.append({
                "depth_g_cm2": depth,
                "material": m.name,
                "cutoff_MeV": round(ec, 2),
                "fluence_above_cutoff_cm2": f"{J0 * math.exp((E30 - ec) / E0_1972):.3e}",
                "skin_dose_Gy_1972_firstorder": round(d, 3),
                "dose_ratio_vs_Al": round(d / d_al, 3),
                "safe_time_factor_vs_Al": round(d_al / d, 2),
            })

    # скан жёсткости: ГБК-50 против Al
    hbc = media[2]
    al = media[0]
    scan = []
    for e0 in E0_SCAN:
        for depth in DEPTHS:
            d_h, _ = dose_behind(hbc, depth, e0)
            d_a, _ = dose_behind(al, depth, e0)
            scan.append({"E0_MeV": e0, "depth_g_cm2": depth,
                         "dose_ratio_HBC50_vs_Al": round(d_h / d_a, 3),
                         "safe_time_factor": round(d_a / d_h, 2)})

    # sanity-инварианты
    for depth in DEPTHS:
        row = {r["material"]: r for r in out_rows if r["depth_g_cm2"] == depth}
        assert row["ГБК-50: 50% AB + 50% СВМПЭ"]["dose_ratio_vs_Al"] < 1.0
        assert row["Боразан NH3BH3 (AB)"]["dose_ratio_vs_Al"] <= \
            row["Полиэтилен (CH2)n / СВМПЭ"]["dose_ratio_vs_Al"] + 1e-9

    out = {
        "what": "SEP-доза кожи за слоем (straight-ahead, King AL-1972) + "
                "коэффициент безопасного времени EVA во время события",
        "date": "2026-08-30",
        "spectrum": "J(>E)=7.9e9*exp((30-E)/26.5) cm^-2 (King 1974, AL-1972; "
                    "сверено 2026-08-30) + скан E0=10/26.5/50",
        "assumptions": "CSDA straight-ahead; без ядерных взаимодействий, "
                       "рассеяния, вторичных, геометрии тела; кожа=0.1 г/см² "
                       "воды; S клампится ниже 1 МэВ. Абсолютные Гр — первый "
                       "порядок, НЕ сверены с BRYNTRN; выводить только отношения.",
        "gcr_note": "Хроническая доза на поверхности Луны (ГКЛ+альбедо) "
                    "1369 мкЗв/день (Chang'E-4 LND, Sci.Adv. 2020, "
                    "10.1126/sciadv.aaz1334) слоем скафандра практически не "
                    "снижается — коэффициенты ниже относятся ТОЛЬКО к SEP.",
        "results_1972": out_rows,
        "hardness_scan_HBC50_vs_Al": scan,
    }
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "radshield_sep_dose_results.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print(f"OK -> {path}")
    print(f"{'глубина':>8s} {'материал':38s} {'Ec,МэВ':>7s} "
          f"{'D,Гр(1972)':>10s} {'D/D_Al':>7s} {'t_safe×':>8s}")
    for r in out_rows:
        print(f"{r['depth_g_cm2']:8.1f} {r['material']:38s} "
              f"{r['cutoff_MeV']:7.2f} {r['skin_dose_Gy_1972_firstorder']:10.3f} "
              f"{r['dose_ratio_vs_Al']:7.3f} {r['safe_time_factor_vs_Al']:8.2f}")
    print("\nскан жёсткости (ГБК-50 vs Al):")
    for s in scan:
        print(f"  E0={s['E0_MeV']:5.1f}  x={s['depth_g_cm2']:3.1f} г/см²  "
              f"D/D_Al={s['dose_ratio_HBC50_vs_Al']:.3f}  "
              f"t_safe×{s['safe_time_factor']:.2f}")


if __name__ == "__main__":
    main()
