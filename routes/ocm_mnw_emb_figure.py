#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/ocm_mnw_emb_figure.py — фигура Этапа 16: ΔΔE‡ на EMBEDDED Na/W-кластере
Mn=O. Левая панель — HAT-профили CH4/C2H6 до TS (точки 0..2, интерьерный TS)
на DFT (сплошные) и NEVPT2 (пунктир); подписи — барьеры NEVPT2. Правая —
дескриптор ΔΔE‡ по уровням: bare (Этап 15) vs embedded (Этап 16), с квантовым
сдвигом. Все числа — из routes/ocm_mnw_emb_final.json (+ bare для сравнения).
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = os.path.dirname(os.path.abspath(__file__))
ASSETS = os.path.join(DIR, "..", "assets")
C = {"ch4": "#2a78d6", "c2h6": "#eda100"}
INK, MUT, GRID = "#0b0b0b", "#898781", "#e1e0d9"
NAME = {"ch4": "CH₄", "c2h6": "C₂H₆"}

R = json.load(open(os.path.join(DIR, "ocm_mnw_emb_final.json")))
sub = R["substrates"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.8, 5.4), width_ratios=[1.2, 1])
fig.patch.set_facecolor("#fcfcfb")

# --- панель 1: профили R->TS (точки 0..2, интерьерный TS=2) ---
for s in ("ch4", "c2h6"):
    g = sub[s]
    n = g["nevpt2_ts_index"] + 1              # R..TS (3 точки)
    xs = list(range(n))
    ax1.plot(xs, g["dft_rel_kcal"][:n], "-o", color=C[s], lw=2, ms=6,
             label=f"{NAME[s]} DFT")
    ax1.plot(xs, g["nevpt2_rel_kcal"][:n], "--s", color=C[s], lw=1.7, ms=6,
             alpha=0.8, label=f"{NAME[s]} NEVPT2")
    i = g["nevpt2_ts_index"]
    ax1.annotate(f"{g['nevpt2_barrier_kcal']:.1f}", (i, g["nevpt2_rel_kcal"][i]),
                 textcoords="offset points", xytext=(0, 9), ha="center",
                 fontsize=11, fontweight="bold", color=C[s])
ax1.set_xlabel("HAT-путь R→TS (пин r(C–H)+r(O–H)) — интерьерный TS", fontsize=10)
ax1.set_ylabel("ΔE отн. R, ккал/моль", fontsize=10)
ax1.set_title("Embedded Na/W: HAT-барьеры CH₄ vs C₂H₆\nDFT (сплошн.) / NEVPT2 (пунктир); подписи — барьер NEVPT2",
              fontsize=10.5, fontweight="bold", color=INK)
ax1.legend(frameon=False, fontsize=9, loc="upper left", ncol=2)
ax1.grid(axis="y", color=GRID, lw=0.7)
for sp in ("top", "right"):
    ax1.spines[sp].set_visible(False)

# --- панель 2: ΔΔE‡ по уровням, bare (Этап 15) vs embedded (Этап 16) ---
lv = ["dft", "casscf", "nevpt2"]
lab = ["DFT", "CASSCF", "NEVPT2"]
bare = {"dft": 4.7, "casscf": -9.5, "nevpt2": 11.6}     # Этап 15 (bare MnO)
emb = R["ddE_kcal"]
x = range(len(lv))
w = 0.38
ax2.bar([i - w/2 for i in x], [bare[k] for k in lv], w, color=MUT,
        label="bare MnO (Этап 15)")
ax2.bar([i + w/2 for i in x], [emb[k] for k in lv], w, color="#7a3ea8",
        label="embedded Na/W (Этап 16)")
ax2.axhline(0, color=INK, lw=1)
for i, k in enumerate(lv):
    ax2.annotate(f"{emb[k]:+.1f}", (i + w/2, emb[k]),
                 textcoords="offset points", xytext=(0, 4 if emb[k] >= 0 else -12),
                 ha="center", fontsize=9, fontweight="bold", color="#7a3ea8")
ax2.set_xticks(list(x)); ax2.set_xticklabels(lab, fontsize=10)
ax2.set_ylabel("ΔΔE‡ = барьер(C₂H₆) − барьер(CH₄), ккал/моль", fontsize=9.5)
ax2.set_title(f"Селективностный дескриптор ΔΔE‡\nкванто-сдвиг embedded +{R['quantum_shift_kcal']} "
              "(DFT→NEVPT2): анти-BEP направление выжило,\nно ослаблено (embedded ≈ 0)",
              fontsize=10, fontweight="bold", color=INK)
ax2.legend(frameon=False, fontsize=9, loc="lower right")
ax2.grid(axis="y", color=GRID, lw=0.7)
for sp in ("top", "right"):
    ax2.spines[sp].set_visible(False)

fig.tight_layout()
out = os.path.join(ASSETS, "ocm_mnw_emb_dde.png")
fig.savefig(out, dpi=140, facecolor=fig.get_facecolor())
print("wrote", out)
