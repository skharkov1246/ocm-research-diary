#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/ocm_mnw_dopant_volcano_5metal.py — 5-metal anti-BEP volcano figure
(Stages 16–17b): Mn/Cr/Fe/Ti/V in embedded Na/W-кластере. Левая панель —
ΔΔE‡ по уровням DFT/CASSCF/NEVPT2 для пяти металлов; правая — квантовый
сдвиг (DFT→NEVPT2) стрелками, показывает metal-зависимую инверсию ранжирования.
Все числа из ocm_mnw_emb{,_cr,_fe,_ti,_v}_final.json.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = os.path.dirname(os.path.abspath(__file__))
ASSETS = os.path.join(DIR, "..", "assets")
INK, MUT, GRID = "#0b0b0b", "#898781", "#e1e0d9"
POS, NEG = "#2e8b57", "#c0392b"          # селективно (>0) / анти-селективно (<0)
CMET = {"Mn": "#898781", "Cr": "#7a3ea8", "Fe": "#b5651d", "Ti": "#2166ac", "V": "#d6604d"}

FILES = {"Mn": "ocm_mnw_emb_final.json", "Cr": "ocm_mnw_emb_cr_final.json",
         "Fe": "ocm_mnw_emb_fe_final.json", "Ti": "ocm_mnw_emb_ti_final.json",
         "V": "ocm_mnw_emb_v_final.json"}
D = {m: json.load(open(os.path.join(DIR, f))) for m, f in FILES.items()}
METALS = ["Mn", "Cr", "Fe", "Ti", "V"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.0, 5.6), width_ratios=[1.25, 1])
fig.patch.set_facecolor("#fcfcfb")

# --- панель 1: ΔΔE‡ по уровням для 5 металлов ---
lv = ["dft", "casscf", "nevpt2"]
lab = ["DFT", "CASSCF", "NEVPT2"]
w = 0.16
for j, m in enumerate(METALS):
    e = D[m]["ddE_kcal"]
    xs = [i + (j - 2) * w for i in range(len(lv))]
    ax1.bar(xs, [e[k] for k in lv], w, color=CMET[m],
            label=f"{m} (2S={D[m]['spin_2S']})")
    for x, k in zip(xs, lv):
        val = e[k]
        ax1.annotate(f"{val:+.1f}", (x, val), textcoords="offset points",
                     xytext=(0, 3 if val >= 0 else -11), ha="center",
                     fontsize=7, fontweight="bold", color=CMET[m])
ax1.axhline(0, color=INK, lw=1)
ax1.set_xticks(range(len(lv)))
ax1.set_xticklabels(lab, fontsize=10)
ax1.set_ylabel("ΔΔE‡ = барьер(C₂H₆) − барьер(CH₄), ккал/моль", fontsize=9.5)
ax1.set_title("5-металльный вулкан анти-BEP (embedded Na/W)\nΔΔE‡ по уровням теории; >0 = селективно к C₂",
              fontsize=10.5, fontweight="bold", color=INK)
ax1.legend(frameon=False, fontsize=9, loc="lower left", ncol=5)
ax1.grid(axis="y", color=GRID, lw=0.7)
for sp in ("top", "right"):
    ax1.spines[sp].set_visible(False)

# --- панель 2: квантовый сдвиг DFT→NEVPT2 стрелкой на металл ---
for j, m in enumerate(METALS):
    e = D[m]["ddE_kcal"]
    d0, d2 = e["dft"], e["nevpt2"]
    ax2.plot([j, j], [d0, d2], color=CMET[m], lw=2.2, zorder=1)
    ax2.annotate("", xy=(j, d2), xytext=(j, d0),
                 arrowprops=dict(arrowstyle="-|>", color=CMET[m], lw=2.2))
    ax2.scatter([j], [d0], s=42, color=MUT, zorder=3, edgecolor="white")
    up = d2 >= d0
    ax2.scatter([j], [d2], s=90, color=(POS if d2 >= 0 else NEG), zorder=3,
                edgecolor="white")
    ax2.annotate(f"{d2:+.1f}", (j, d2), textcoords="offset points",
                 xytext=(11, -3 if up else 3), ha="left",
                 fontsize=10, fontweight="bold", color=(POS if d2 >= 0 else NEG))
    ax2.annotate(f"DFT {d0:+.1f}", (j, d0), textcoords="offset points",
                 xytext=(11, 2), ha="left", fontsize=8, color=MUT)
    ax2.annotate(f"Δ {D[m]['quantum_shift_kcal']:+.1f}", (j, (d0 + d2) / 2),
                 textcoords="offset points", xytext=(-11, 0), ha="right",
                 fontsize=8.5, fontweight="bold", color=CMET[m])
ax2.axhline(0, color=INK, lw=1, ls="--", alpha=0.6)
ax2.text(4.35, 4.0, "селективно", fontsize=8, color=POS, ha="right")
ax2.text(4.35, -6.8, "анти-селективно", fontsize=8, color=NEG, ha="right")
ax2.set_xticks(range(len(METALS)))
ax2.set_xticklabels([f"{m}\n2S={D[m]['spin_2S']}" for m in METALS], fontsize=9.5)
ax2.set_ylabel("ΔΔE‡, ккал/моль (серое=DFT → цвет=NEVPT2)", fontsize=9)
ax2.set_title("Квантовый сдвиг металл-зависим: НЕВПТ2 ранжирование\nперевёрнуто vs DFT (Ti коррекция +8.6!)",
              fontsize=10, fontweight="bold", color=INK)
ax2.grid(axis="y", color=GRID, lw=0.7)
ax2.set_ylim(-8.5, 5)
for sp in ("top", "right"):
    ax2.spines[sp].set_visible(False)

fig.tight_layout()
out = os.path.join(ASSETS, "ocm_mnw_dopant_volcano_5metal.png")
fig.savefig(out, dpi=140, facecolor=fig.get_facecolor())
print("wrote", out)
