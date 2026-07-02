#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/ocm_mnw_figure.py — фигура Этапа 15: ΔΔE‡ закрыт на едином CAS(9e,10o).
Левая панель — честные HAT-профили обоих субстратов (DFT-скан, те самые точки);
правая — дескриптор ΔΔE‡ по уровням теории (DFT / CASSCF / NEVPT2) против
нулевой линии BEP-стены и цели +4 ккал/моль (70% селективности).
Все числа читаются из routes/ocm_mnw_*_geom.json и *_selectivity_results.json.
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

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 5.4), width_ratios=[1.15, 1])
fig.patch.set_facecolor("#fcfcfb")

# --- панель 1: HAT-профили (rel к точке 0; R = минимум профиля, TS = максимум)
for s in ("ch4", "c2h6"):
    g = R["substrates"][s]
    xs = list(range(len(g["rel_kcal"])))
    ax1.plot(xs, g["rel_kcal"], "-o", color=C[s], lw=2, ms=6, label=NAME[s])
    i, j = g["ts_index"], g["r_index"]
    ax1.plot([j], [g["rel_kcal"][j]], "s", color=C[s], ms=9, mfc="none", mew=1.6)
    ax1.annotate(f"TS  ΔE‡={g['dft_barrier_kcal']:.1f}",
                 (i, g["rel_kcal"][i]), textcoords="offset points",
                 xytext=(0, 9), ha="center", fontsize=10, fontweight="bold",
                 color=C[s])
ax1.set_xlabel("точка HAT-пути (пин r(C–H)+r(O–H)) / pinned HAT point", fontsize=10)
ax1.set_ylabel("ΔE vs R, ккал/моль (UKS-PBE0)", fontsize=10)
ax1.set_title("Честные HAT-профили: подъём → TS → спад\n(каркас релаксирован, TS внутри сетки)",
              fontsize=11, fontweight="bold", color=INK)
ax1.legend(frameon=False, fontsize=10, loc="lower left")
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
              fontsize=11, fontweight="bold", color=INK)
ax2.grid(axis="y", color=GRID, lw=0.7)
for sp in ("top", "right"):
    ax2.spines[sp].set_visible(False)

shift = R["quantum_shift_kcal"]
fig.text(0.01, 0.01,
         f"Реальный прогон: MnO(секстет)-модель центра Mn–Na₂WO₄/SiO₂, UKS-PBE0/def2-SVP геометрия, "
         f"ROHF→AVAS(Mn 3d, O 2p, H 1s)→CASSCF(9e,10o)→NEVPT2, ЕДИНОЕ активное пространство во всех 4 структурах. "
         f"Квантовый сдвиг дескриптора NEVPT2−DFT = {shift:+.1f} ккал/моль. Минимальный кластер — не периодический катализатор.",
         fontsize=7.5, style="italic", color=MUT, wrap=True)
plt.tight_layout(rect=(0, 0.05, 1, 1))
out = os.path.join(ASSETS, "ocm_mnw_dde.png")
plt.savefig(out, dpi=110, facecolor=fig.get_facecolor())
print("saved", out)
