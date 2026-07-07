#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/ocm_mnw_dopant_figure.py — фигура Этапа 17: допант-панель анти-BEP на
EMBEDDED Na/W-кластере. Активный оксо-металл Mn(Э16) / Cr / Fe в ОДНОМ каркасе,
единый протокол ΔΔE‡. Левая панель — ΔΔE‡ по уровням DFT/CASSCF/NEVPT2 для трёх
металлов; правая — квантовый сдвиг (DFT→NEVPT2) как стрелка на каждый металл,
показывает металл-зависимость знака сдвига и ПЕРЕВОРОТ ранжирования (DFT «лучший»
Fe → NEVPT2 «лучший» Cr). Все числа — из ocm_mnw_emb{,_cr,_fe}_final.json.
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
CMET = {"Mn": "#898781", "Cr": "#7a3ea8", "Fe": "#b5651d"}

FILES = {"Mn": "ocm_mnw_emb_final.json", "Cr": "ocm_mnw_emb_cr_final.json",
         "Fe": "ocm_mnw_emb_fe_final.json"}
D = {m: json.load(open(os.path.join(DIR, f))) for m, f in FILES.items()}
METALS = ["Mn", "Cr", "Fe"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 5.4), width_ratios=[1.25, 1])
fig.patch.set_facecolor("#fcfcfb")

# --- панель 1: ΔΔE‡ по уровням для 3 металлов ---
lv = ["dft", "casscf", "nevpt2"]
lab = ["DFT", "CASSCF", "NEVPT2"]
w = 0.26
for j, m in enumerate(METALS):
    e = D[m]["ddE_kcal"]
    xs = [i + (j - 1) * w for i in range(len(lv))]
    ax1.bar(xs, [e[k] for k in lv], w, color=CMET[m],
            label=f"{m} (2S={D[m]['spin_2S']})")
    for x, k in zip(xs, lv):
        ax1.annotate(f"{e[k]:+.1f}", (x, e[k]), textcoords="offset points",
                     xytext=(0, 3 if e[k] >= 0 else -11), ha="center",
                     fontsize=8, fontweight="bold", color=CMET[m])
ax1.axhline(0, color=INK, lw=1)
ax1.set_xticks(range(len(lv)))
ax1.set_xticklabels(lab, fontsize=10)
ax1.set_ylabel("ΔΔE‡ = барьер(C₂H₆) − барьер(CH₄), ккал/моль", fontsize=9.5)
ax1.set_title("Допант-панель анти-BEP (embedded Na/W)\nΔΔE‡ по уровням теории; "
              ">0 = селективно к C₂", fontsize=10.5, fontweight="bold", color=INK)
ax1.legend(frameon=False, fontsize=9, loc="lower left", ncol=3)
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
    ax2.annotate(f"сдвиг {D[m]['quantum_shift_kcal']:+.1f}", (j, (d0 + d2) / 2),
                 textcoords="offset points", xytext=(-11, 0), ha="right",
                 fontsize=8.5, fontweight="bold", color=CMET[m])
ax2.axhline(0, color=INK, lw=1, ls="--", alpha=0.6)
ax2.text(2.35, 0.25, "селективно", fontsize=8, color=POS, ha="right")
ax2.text(2.35, -0.9, "анти-селективно", fontsize=8, color=NEG, ha="right")
ax2.set_xticks(range(len(METALS)))
ax2.set_xticklabels([f"{m}\n2S={D[m]['spin_2S']}" for m in METALS], fontsize=9.5)
ax2.set_ylabel("ΔΔE‡, ккал/моль (серое=DFT → цвет=NEVPT2)", fontsize=9)
ax2.set_title("Квантовый сдвиг металл-зависим\nCr/Mn ↑ к селективности, Fe ↓ — "
              "DFT-ранжирование перевёрнуто", fontsize=10, fontweight="bold",
              color=INK)
ax2.grid(axis="y", color=GRID, lw=0.7)
for sp in ("top", "right"):
    ax2.spines[sp].set_visible(False)

fig.tight_layout()
out = os.path.join(ASSETS, "ocm_mnw_dopant_dde.png")
fig.savefig(out, dpi=140, facecolor=fig.get_facecolor())
print("wrote", out)
