#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/ocm_mnw_figure.py — фигура Этапа 15: ΔΔE‡ закрыт на едином CAS(9e,10o).
Левая панель — HAT-профили обоих субстратов на DFT (сплошные) и NEVPT2
(пунктир): каждая поверхность со своими R/TS. Правая — дескриптор ΔΔE‡ по
уровням теории против BEP-стены (0) и цели +4 ккал/моль (70% селективности).
Все числа читаются из routes/ocm_mnw_selectivity_results.json.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = os.path.dirname(os.path.abspath(__file__))
ASSETS = os.path.join(DIR, "..", "assets")
C = {"ch4": "#2a78d6", "c2h6": "#eda100"}   # категориальные слоты 1/3, валидировано
INK, MUT, GRID = "#0b0b0b", "#898781", "#e1e0d9"
NAME = {"ch4": "CH₄", "c2h6": "C₂H₆"}

with open(os.path.join(DIR, "ocm_mnw_selectivity_results.json")) as f:
    R = json.load(f)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.8, 5.6), width_ratios=[1.2, 1])
fig.patch.set_facecolor("#fcfcfb")

# --- панель 1: профили DFT (сплошные) и NEVPT2 (пунктир), свои R/TS
for s in ("ch4", "c2h6"):
    g = R["substrates"][s]
    xs = list(range(len(g["dft_rel_kcal"])))
    ax1.plot(xs, g["dft_rel_kcal"], "-o", color=C[s], lw=2, ms=5,
             label=f"{NAME[s]} DFT")
    ax1.plot(xs, g["nevpt2_rel_kcal"], "--s", color=C[s], lw=1.6, ms=5,
             alpha=0.75, label=f"{NAME[s]} NEVPT2")
    i = g["nevpt2_ts_index"]
    ax1.annotate(f"{g['nevpt2_barrier_kcal']:.1f}",
                 (i, g["nevpt2_rel_kcal"][i]), textcoords="offset points",
                 xytext=(0, 8), ha="center", fontsize=10, fontweight="bold",
                 color=C[s])
ax1.set_xlabel("точка HAT-пути (пин r(C–H)+r(O–H)) / pinned HAT point", fontsize=10)
ax1.set_ylabel("ΔE отн. точки 0, ккал/моль", fontsize=10)
ax1.set_title("HAT-профили: DFT (сплошные) и NEVPT2 (пунктир)\nкаждая поверхность — свои R и TS; подписи: барьер NEVPT2",
              fontsize=10.5, fontweight="bold", color=INK)
ax1.legend(frameon=False, fontsize=9, loc="upper left", ncol=2)
ax1.grid(axis="y", color=GRID, lw=0.7)
for sp in ("top", "right"):
    ax1.spines[sp].set_visible(False)

# --- панель 2: ΔΔE‡ по уровням теории
levels = ["dft", "casscf", "nevpt2"]
labels = ["DFT\n(PBE0)", "CASSCF\n(стат. корр.)", "NEVPT2\n(+ динам. корр.)"]
vals = [R["ddE_kcal"][l] for l in levels]
bars = ax2.bar(labels, vals, width=0.55, color=["#c3c2b7", "#86b6ef", "#2a78d6"],
               edgecolor="none")
for b, v in zip(bars, vals):
    ax2.annotate(f"{v:+.1f}", (b.get_x() + b.get_width() / 2, v),
                 textcoords="offset points", xytext=(0, 6 if v >= 0 else -16),
                 ha="center", fontsize=12, fontweight="bold", color=INK)
ax2.axhline(0, color=MUT, lw=1)
ax2.axhline(4, color="#0ca30c", lw=1.4, ls="--")
ax2.text(2.42, 4, " +4 = 70% C₂\n(карта потолка)", fontsize=9, color="#006300",
         va="center")
ax2.text(2.42, 0, " 0 = BEP-стена", fontsize=9, color=MUT, va="center")
ax2.set_xlim(-0.6, 3.4)
ax2.set_ylabel("ΔΔE‡ = барьер(C₂H₆) − барьер(CH₄), ккал/моль", fontsize=10)
ax2.set_title("Дескриптор селективности по уровням теории\n(единое CAS(9e,10o), оба субстрата досчитаны)",
              fontsize=10.5, fontweight="bold", color=INK)
ax2.grid(axis="y", color=GRID, lw=0.7)
for sp in ("top", "right"):
    ax2.spines[sp].set_visible(False)

shift = R["quantum_shift_kcal"]
fig.text(0.01, 0.005,
         f"Реальный прогон: MnO(секстет) — минимальная модель центра Mn–Na₂WO₄/SiO₂ (НЕ периодический катализатор). "
         f"UKS-PBE0/def2-SVP геометрия (полировка стабильностью + двунаправленные цепочки плотности); "
         f"ROHF→fixed-AVAS(Mn 3d, O 2p, H 1s)→CASSCF(9e,10o)→SC-NEVPT2 во всех 16 точках, орбитально-цепной ремонт выбросов. "
         f"Квантовый сдвиг дескриптора NEVPT2−DFT = {shift:+.1f} ккал/моль.",
         fontsize=7.3, style="italic", color=MUT, wrap=True)
plt.tight_layout(rect=(0, 0.055, 1, 1))
out = os.path.join(ASSETS, "ocm_mnw_dde.png")
plt.savefig(out, dpi=110, facecolor=fig.get_facecolor())
print("saved", out)
