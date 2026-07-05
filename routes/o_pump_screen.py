#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/o_pump_screen.py — скрин №2 off-label карты: автомобильная лямбда-зонда
(YSZ, серийная в каждой машине) как электрохимический O²⁻-насос для OCM.

Идея. В OCM кислород co-feed'ом — перегорание неизбежно (газовый O₂ жжёт C₂).
Off-label: YSZ-элемент качает решёточный O²⁻ электрохимически, импульсами,
под потенциалом U. Эффективная энергия отдачи решёточного O носителем сдвигается:
    Ev_eff(U) = Ev − 2·F·U        (2 электрона на O²⁻; 1 эВ = 23.061 ккал/моль)
Т.е. НАПРЯЖЕНИЕ становится ручкой, двигающей носитель ПО вулкану Сабатье —
композиция больше не приговор. Наш вклад: Ev по 8 металлам уже посчитан
CASSCF/NEVPT2 (chemloop_vacancy_results.json) — там, где DFT врёт на десятки
ккал/моль, окно напряжений по DFT было бы просто неверным.

Селективное окно из нашего же волкана: рабочие chemical-looping носители
(Cr/Mn/Fe-band) легли в Ev_NEVPT2 ≈ [−10, +45] ккал/моль — берём его как
целевую полосу. Выход: для каждого металла диапазон U, вгоняющий Ev_eff в окно.
Чистая арифметика поверх готовых чисел — новых расчётов нет.
"""
import json
import os

DIR = os.path.dirname(os.path.abspath(__file__))
EV_KCAL = 23.061
WINDOW = (-10.0, 45.0)   # ккал/моль, селективная полоса из нашего NEVPT2-волкана


def main():
    src = json.load(open(os.path.join(DIR, "chemloop_vacancy_results.json")))
    rows = []
    for c in src["carriers"]:
        ev = c.get("Ev_nevpt2")
        ev_dft = c.get("Ev_dft")
        if ev is None:
            continue
        # Ev_eff = Ev - 2*23.061*U  ∈ [lo, hi]  =>  U ∈ [(Ev-hi)/46.1, (Ev-lo)/46.1]
        u_lo = (ev - WINDOW[1]) / (2 * EV_KCAL)
        u_hi = (ev - WINDOW[0]) / (2 * EV_KCAL)
        # то же по DFT — чтобы показать, насколько DFT-окно промахивается
        du_dft = ((ev_dft - ev) / (2 * EV_KCAL)) if ev_dft is not None else None
        rows.append({"metal": c["metal"], "Ev_nevpt2": ev, "Ev_dft": ev_dft,
                     "U_window_V": [round(u_lo, 2), round(u_hi, 2)],
                     "dft_window_error_V": round(du_dft, 2) if du_dft is not None
                     else None})
    out = {
        "screen": "off-label #6: automotive lambda-sensor (YSZ) as pulsed O2- pump for OCM",
        "relation": "Ev_eff(U) = Ev - 2*F*U; window = NEVPT2 selective band "
                    f"[{WINDOW[0]}, {WINDOW[1]}] kcal/mol from our chemloop volcano",
        "carriers": rows,
        "reading": {
            "knob": "voltage moves any carrier along the Sabatier volcano — "
                    "composition is no longer destiny; pulsed dosing decouples "
                    "CH4 activation from C2 overoxidation (no gas-phase O2)",
            "our_edge": "U-windows computed on NEVPT2 Ev; dft_window_error_V shows "
                        "how far a DFT-designed voltage window would miss (up to "
                        "volts for early metals) — DFT-blind design fails here",
            "hardware": "mass-produced YSZ lambda sensors / SOFC electrolyte plates; "
                        "wellhead skid compatible"},
    }
    path = os.path.join(DIR, "o_pump_screen_results.json")
    json.dump(out, open(path, "w"), indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
