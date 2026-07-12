#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Фигура Этапа 21: (а) ΔΔE‡(F) по уровням — инверсия знака полевого отклика;
(б) ΔΔ_vol(M−Mn) — переворот ранжира летучести. Числа из results-JSON."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = os.path.dirname(os.path.abspath(__file__))
fld = json.load(open(os.path.join(DIR, "ocm_field_dde3_results.json")))
vol = json.load(open(os.path.join(DIR, "ocm_vol_results.json")))["desc"]["levels"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
F = [-0.257, 0.0, 0.257]
keys = ["-0.005", "+0.000", "+0.005"]
colors = {"dft": "#888888", "casscf": "#4477aa", "nevpt2": "#cc3311"}
for lv in ("dft", "casscf", "nevpt2"):
    y = [fld["ddE"][k][lv] for k in keys]
    ax1.plot(F, y, "o-", color=colors[lv], label=lv.upper(), lw=2)
ax1.axhline(0, color="k", lw=0.6)
ax1.set_xlabel("Поле вдоль Cr→O, В/Å")
ax1.set_ylabel("ΔΔE‡ (C₂H₆−CH₄), ккал/моль")
ax1.set_title("Этап 21-F: полевой отклик селективности\n"
              "DFT −4.75 vs NEVPT2 +0.37 ккал/(В/Å) — инверсия знака")
ax1.legend()
ax1.grid(alpha=0.3)

metals = ["Cr", "W"]
x = range(len(metals))
w = 0.25
for i, lv in enumerate(("dft", "casscf", "nevpt2")):
    dd = vol[lv]["ddE_vs_Mn_kcal"]
    ax2.bar([xx + (i - 1) * w for xx in x],
            [dd.get(f"{m}-Mn") for m in metals], w,
            color=colors[lv], label=lv.upper())
ax2.axhline(0, color="k", lw=0.6)
ax2.axhline(-20, color="r", ls="--", lw=1, label="kill-линия A2")
ax2.set_xticks(list(x))
ax2.set_xticklabels([f"{m}−Mn" for m in metals])
ax2.set_ylabel("ΔΔE_vol vs Mn, ккал/моль")
ax2.set_title("Этап 21-B: летучесть допанта\n"
              "DFT: красный флаг · NEVPT2: Cr>W>Mn по стойкости")
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3, axis="y")

fig.tight_layout()
out = os.path.join(DIR, "..", "assets", "ocm_stage21_field_vol.png")
fig.savefig(out, dpi=150)
print("saved", out)
