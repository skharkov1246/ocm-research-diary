#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/ocm_mnw_7metal_map.py — финальная карта допант-вулкана (Этапы 15–19):
7 металлов (Mn/Cr/Fe 3d-окно, Ti/V ранние 3d, Zr/Nb 4d-гомологи) в одном
embedded Na/W-каркасе, robust NEVPT2-ридаут. Cr — единственный селективный.
V/Nb — неактивны (перевала нет); Ti/Zr — шумные ~0/минус. Числа:
final_robust.json (Ti/V/Cr/Fe/Zr/Nb) + Этап 16 FINAL (Mn, чистый 5-точечный).
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

# (metal, dde_nevpt2, статус, подпись)
ROWS = []
mn = json.load(open(os.path.join(DIR, "ocm_mnw_emb_final.json")))
ROWS.append(("Mn", mn["ddE_kcal"]["nevpt2"], "ok", "нейтрален"))
for m in ("Cr", "Fe", "Ti", "Zr"):
    r = json.load(open(os.path.join(DIR, f"ocm_mnw_emb_{m.lower()}_final_robust.json")))
    dde = r["ddE_kcal"]["nevpt2"]
    note = {"Cr": "селективен ✓ (2 сетки)", "Fe": "анти-селективен",
            "Ti": "≈0, шумно", "Zr": "минус, шумно"}[m]
    ROWS.append((m, dde, "ok", note))
for m, note in (("V", "неактивен: перевала нет (~59 ккал в гору)"),
                ("Nb", "неактивен: перевала нет (~67 ккал)")):
    ROWS.append((m, None, "inactive", note))
ORDER = ["Ti", "V", "Cr", "Mn", "Fe", "Zr", "Nb"]
ROWS.sort(key=lambda r: ORDER.index(r[0]))

fig, ax = plt.subplots(figsize=(10.5, 5.2))
fig.patch.set_facecolor("#fcfcfb")
for j, (m, dde, st, note) in enumerate(ROWS):
    if st == "ok":
        c = POS if dde > 0.3 else (MUT if abs(dde) <= 0.6 else NEG)
        ax.bar(j, dde, 0.55, color=c)
        ax.annotate(f"{dde:+.2f}", (j, dde), textcoords="offset points",
                    xytext=(0, 5 if dde >= 0 else -13), ha="center",
                    fontsize=10, fontweight="bold", color=c)
        ax.annotate(note, (j, min(dde, 0) - 1.9), ha="center", fontsize=7.6,
                    color=INK, alpha=0.8)
    else:
        ax.bar(j, 0, 0.55, color=MUT)
        ax.scatter([j], [0], marker="x", s=90, color=NEG, zorder=3)
        ax.annotate("неактивен\n(перевала нет)", (j, 0.6), ha="center",
                    fontsize=8, color=NEG, fontweight="bold")
ax.axhline(0, color=INK, lw=1)
ax.axhspan(0, 5.5, color=POS, alpha=0.05)
ax.text(6.42, 4.9, "селективно к C₂", fontsize=8.5, color=POS, ha="right")
ax.set_xticks(range(len(ROWS)))
ax.set_xticklabels([f"{m}\n{'3d' if m in ('Ti','V','Cr','Mn','Fe') else '4d'}"
                    for m, *_ in ROWS], fontsize=10.5)
ax.set_ylabel("ΔΔE‡(NEVPT2) robust, ккал/моль", fontsize=10)
ax.set_ylim(-10.5, 5.5)
ax.set_title("Финальная карта допант-вулкана (Этапы 15–19): 7 металлов, один каркас\n"
             "Cr — единственный подтверждённый селективный центр; группы 4/5 не работают в обоих рядах",
             fontsize=10.5, fontweight="bold", color=INK)
ax.grid(axis="y", color=GRID, lw=0.7)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.tight_layout()
out = os.path.join(ASSETS, "ocm_mnw_7metal_map.png")
fig.savefig(out, dpi=140, facecolor=fig.get_facecolor())
print("wrote", out)
