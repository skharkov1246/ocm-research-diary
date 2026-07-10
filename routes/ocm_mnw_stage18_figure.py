#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/ocm_mnw_stage18_figure.py — фигура Этапа 18: стресс-тест допант-вулкана
расширенной сеткой (5→9 точек). Слева — ΔΔE‡(NEVPT2) «5 точек vs 9 точек
(robust)»: Cr/Fe воспроизведены точь-в-точь, Ti-«сенсация» схлопнулась к нулю,
V нечитаем. Справа — NEVPT2-профили CH₄, объясняющие почему: чистые перевалы
Cr/Fe (шип-выброс Cr замаскирован), монотонный подъём V (перевала нет),
шумная поверхность Ti. Числа: ocm_mnw_emb_*_final_robust.json + Этап 17/17b.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIR = os.path.dirname(os.path.abspath(__file__))
ASSETS = os.path.join(DIR, "..", "assets")
INK, MUT, GRID = "#0b0b0b", "#898781", "#e1e0d9"
POS, NEG = "#2e8b57", "#c0392b"
CMET = {"Cr": "#7a3ea8", "Fe": "#b5651d", "Ti": "#2166ac", "V": "#d6604d"}

R = {m: json.load(open(os.path.join(DIR, f"ocm_mnw_emb_{m.lower()}_final_robust.json")))
     for m in ("Ti", "V", "Cr", "Fe")}
OLD5 = {"Ti": 4.18, "V": 0.85, "Cr": 0.73, "Fe": -7.19}   # Этап 17/17b, 5 точек

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.0, 5.4), width_ratios=[1, 1.15])
fig.patch.set_facecolor("#fcfcfb")

# --- панель 1: 5-точечный vs 9-точечный robust ΔΔE‡ (NEVPT2) ---
metals = ["Cr", "Fe", "Ti", "V"]
w = 0.34
for j, m in enumerate(metals):
    old = OLD5[m]
    ax1.bar(j - w / 2, old, w, color=CMET[m], alpha=0.35,
            label="5 точек (Э17/17b)" if j == 0 else None)
    new = (R[m].get("ddE_kcal") or {}).get("nevpt2")
    if new is not None:
        ax1.bar(j + w / 2, new, w, color=CMET[m],
                label="9 точек robust (Э18)" if j == 0 else None)
        ax1.annotate(f"{new:+.2f}", (j + w / 2, new), textcoords="offset points",
                     xytext=(0, 4 if new >= 0 else -12), ha="center",
                     fontsize=9, fontweight="bold", color=CMET[m])
    else:
        ax1.text(j + w / 2, 1.5, "9 т.:\nнечитаем", ha="center", va="bottom",
                 fontsize=8.5, color=NEG, fontweight="bold")
    ax1.annotate(f"{old:+.2f}", (j - w / 2, old), textcoords="offset points",
                 xytext=(0, 4 if old >= 0 else -12), ha="center",
                 fontsize=8, color=CMET[m], alpha=0.75)
ax1.axhline(0, color=INK, lw=1)
ax1.set_xticks(range(len(metals)))
ax1.set_xticklabels(metals, fontsize=11)
ax1.set_ylabel("ΔΔE‡(NEVPT2), ккал/моль", fontsize=10)
ax1.set_title("Стресс-тест сеткой: что выжило\nCr/Fe воспроизведены, Ti-«сенсация» снята, V выбыл",
              fontsize=10.5, fontweight="bold", color=INK)
ax1.legend(frameon=False, fontsize=9, loc="lower left")
ax1.grid(axis="y", color=GRID, lw=0.7)
for sp in ("top", "right"):
    ax1.spines[sp].set_visible(False)

# --- панель 2: NEVPT2-профили CH₄ на 9 точках ---
for m in metals:
    rel = R[m]["substrates"]["ch4"]["nevpt2"]["rel_kcal"]
    drop = R[m]["substrates"]["ch4"]["nevpt2"]["dropped_points"]
    xs = list(range(len(rel)))
    xs_ok = [x for x in xs if x not in drop]
    ax2.plot(xs_ok, [rel[x] for x in xs_ok], "-o", color=CMET[m], lw=1.8,
             ms=4.5, label=m)
    for x in drop:
        y = min(rel[x], 70)          # шип за осью — крестим на кромке
        ax2.scatter([x], [y], marker="x", s=80, color=CMET[m], zorder=4)
        ax2.annotate(f"выброс +{rel[x]:.0f} → маска", (x, y),
                     textcoords="offset points", xytext=(-4, -12), ha="right",
                     fontsize=7.5, fontweight="bold", color=CMET[m])
ax2.axhline(0, color=INK, lw=0.8, ls="--", alpha=0.5)
ax2.set_xlabel("точка пути HAT (R → продукты)", fontsize=9.5)
ax2.set_ylabel("NEVPT2 rel, ккал/моль (CH₄)", fontsize=9.5)
ax2.set_title("Почему: Cr/Fe — чистые интерьерные перевалы;\nV растёт до края (перевала нет); Ti слаб (35 ккал)",
              fontsize=10, fontweight="bold", color=INK)
ax2.legend(frameon=False, fontsize=9, ncol=4, loc="upper left")
ax2.grid(axis="y", color=GRID, lw=0.7)
ax2.set_ylim(-40, 75)
for sp in ("top", "right"):
    ax2.spines[sp].set_visible(False)

fig.tight_layout()
out = os.path.join(ASSETS, "ocm_mnw_stage18_stress.png")
fig.savefig(out, dpi=140, facecolor=fig.get_facecolor())
print("wrote", out)
