#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/econ_ysz_flux.py — закрытие риска №3 реестра (METHANE_INVENTION_ECONOMICS):
несёт ли YSZ-стек поток O²⁻, нужный скиду I4, и что это стоит по площади/мощности.
Параметрика, НЕ CFD: проводимость 8YSZ из литературы (σ≈0.02 С/см при 800°C,
E_a≈0.9 эВ — учебные значения SOFC), j = σ·U/L (омическое приближение, тонкий
электролит, без поляризационных потерь — оптимистично, честно помечено).

Плюс: перевод флагманского ΔΔE‡(Cr)=+0.73 ккал в кинетический множитель
селективности при T OCM — ключевой «не обосраться»-риск дескриптора.
"""
import json
import math
import os

DIR = os.path.dirname(os.path.abspath(__file__))
KB_EV = 8.617e-5          # эВ/К
R_KCAL = 1.987e-3         # ккал/(моль·К)
F = 96485.0               # Кл/моль e-
YEAR_S = 3.156e7

# --- YSZ (8 мол.% Y2O3): σ(T) = σ0·exp(-Ea/kT); калибровка на 0.02 С/см @ 1073К
EA_EV = 0.90
SIGMA_1073 = 2.0          # С/м (=0.02 С/см)
SIGMA0 = SIGMA_1073 / math.exp(-EA_EV / (KB_EV * 1073))

def sigma(T):
    return SIGMA0 * math.exp(-EA_EV / (KB_EV * T))

# --- скид 1 MMscf/д: 7300 т CH4/год; OCM 2CH4 + O2 -> C2H4 + 2H2O
CH4_T_Y = 7300.0
Y = 0.35                                  # таргет-выход C2H4
C2H4_T_Y = CH4_T_Y * Y * (28.0 / (2 * 16.04))
MOL_C2H4_Y = C2H4_T_Y * 1e6 / 28.0        # моль/год
MOL_O2_Y = MOL_C2H4_Y                     # 1 O2 на 1 C2H4
MOL_O_Y = 2 * MOL_O2_Y                    # атомов O (=O2- через решётку)
I_FULL = MOL_O_Y * 2 * F / YEAR_S         # Ампер: полный кислород через YSZ

OUT = {"assumptions": {
    "sigma_1073K_S_per_m": SIGMA_1073, "Ea_eV": EA_EV,
    "electrolyte_um": 50, "U_V": 0.5,
    "note": "омическое приближение, без электродной поляризации (оптимистично); "
            "skid 1 MMscf/d, Y=0.35, 2CH4+O2->C2H4+2H2O"},
    "skid": {"CH4_t_y": CH4_T_Y, "C2H4_t_y": round(C2H4_T_Y, 0),
             "O_mol_y": round(MOL_O_Y, 0), "I_full_A": round(I_FULL, 0)},
    "scenarios": {}}

L = 50e-6                                 # м, тонкоплёночный электролит
U = 0.5                                   # В на ячейку
for T in (973, 1073, 1173):
    s = sigma(T)
    j = s * U / L                         # А/м^2 (омический предел)
    j_prac = min(j, 10000.0)              # практический потолок ~1 А/см^2 (SOFC)
    for share, tag in ((1.0, "full_pump"), (0.10, "trim_10pct"),
                       (0.01, "trim_1pct")):
        I_need = I_FULL * share
        area = I_need / j_prac            # м^2 стека
        power_kw = I_need * U / 1000.0    # электрическая мощность
        OUT["scenarios"][f"{T}K_{tag}"] = {
            "sigma_S_per_m": round(s, 3),
            "j_ohmic_A_per_m2": round(j, 0),
            "j_used_A_per_m2": round(j_prac, 0),
            "I_A": round(I_need, 0), "stack_area_m2": round(area, 1),
            "power_kW": round(power_kw, 1),
            "stack_cost_USD_at_3k_per_m2": round(area * 3000, 0)}

# --- дескриптор -> кинетика: множитель "щажения" C2H6 при T
DDE = 0.73                                # ккал/моль, Cr robust
OUT["descriptor_kinetics"] = {}
for T in (873, 973, 1073):
    k = math.exp(DDE / (R_KCAL * T))
    OUT["descriptor_kinetics"][f"{T}K"] = {
        "selectivity_multiplier": round(k, 2),
        "note": "exp(ddE/RT): во сколько раз медленнее центр рвёт C2H6 vs CH4"}
# сколько ddE нужно для 3х/10х при 1073К
for target in (3.0, 10.0):
    OUT["descriptor_kinetics"][f"ddE_needed_for_{target}x_1073K_kcal"] = round(
        math.log(target) * R_KCAL * 1073, 2)

with open(os.path.join(DIR, "econ_ysz_flux.json"), "w") as f:
    json.dump(OUT, f, indent=1)
print(json.dumps(OUT, indent=1, ensure_ascii=False))
