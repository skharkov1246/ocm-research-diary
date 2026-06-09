#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
routes/co2_to_fuels.py — фигуры вкладки co2-to-fuels (электровосстановление CO₂ → C2+).

Две схемы (по образцу routes/polyolefins.py):
  06_co2_to_fuels_molecules.png   — механизм C–C сочленения *CO–*CO → C2+ (схема, illustrative)
  06_co2_to_fuels_descriptor.png  — дескриптор ΔE‡(CASCI−DFT): профиль барьера + сито сайтов

Реальные числа (профиль CASSCF vs PBE, NOON в TS, ΔE‡ по сайтам) берутся из
routes/co2_to_fuels_results.json, который пишет routes/co2_to_fuels_calc.py (PySCF).
Всё, для чего нет расчёта, помечено «illustrative». Литературные якоря помечены «лит.».

Запуск:  python3 routes/co2_to_fuels.py
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.normpath(os.path.join(HERE, "..", "assets", "matterforge"))
RESULTS = os.path.join(HERE, "co2_to_fuels_results.json")

# ── палитра (в тон дневнику: фиолетовые заголовки, зелёный — цвет вкладки) ──
PURPLE = "#7c5cff"
GREEN  = "#48bb78"
INK    = "#1a202c"
GREY   = "#94a3b8"
RED    = "#e53e3e"
COPPER = "#d97a3a"
ALUM   = "#9aa3b2"
OXY    = "#e05a5a"
CARB   = "#2d3748"
plt.rcParams["font.family"] = "DejaVu Sans"


def _load_results():
    if os.path.exists(RESULTS):
        with open(RESULTS) as f:
            return json.load(f)
    return None


# ════════════════════════════════════════════════════════════════════════
#  ФИГУРА 1 — молекулы / механизм C–C сочленения
# ════════════════════════════════════════════════════════════════════════
def atom(ax, xy, color, r=0.16, label=None, lc="white", fs=12, ec=INK):
    ax.add_patch(Circle(xy, r, fc=color, ec=ec, lw=1.6, zorder=5))
    if label:
        ax.text(xy[0], xy[1], label, ha="center", va="center",
                color=lc, fontsize=fs, fontweight="bold", zorder=6)


def bond(ax, a, b, lw=3.2, color=CARB, ls="-", z=4):
    ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=lw, ls=ls,
            solid_capstyle="round", zorder=z)


def double_bond(ax, a, b, off=0.05, **kw):
    dx, dy = b[0]-a[0], b[1]-a[1]
    L = np.hypot(dx, dy); nx, ny = -dy/L*off, dx/L*off
    bond(ax, (a[0]+nx, a[1]+ny), (b[0]+nx, b[1]+ny), **kw)
    bond(ax, (a[0]-nx, a[1]-ny), (b[0]-nx, b[1]-ny), **kw)


def fig_molecules(res):
    fig, ax = plt.subplots(figsize=(14.2, 7.0), dpi=100)
    ax.set_xlim(0, 14.2); ax.set_ylim(0, 7.0); ax.axis("off")

    ax.text(7.1, 6.62, "CO₂ → C2+ (этилен/этанол): сочленение C–C через  *CO–*CO",
            ha="center", fontsize=18, fontweight="bold", color=INK)
    ax.text(7.1, 6.22, "узкое место всего маршрута — барьер димеризации *CO; делает (почти только) медь   (иллюстративно)",
            ha="center", fontsize=11.5, color=GREY)

    # ── ЗОНА A: CO₂ → *CO на Cu ─────────────────────────────────────────
    ax.text(2.05, 5.55, "1 · АКТИВАЦИЯ CO₂", ha="center", fontsize=12.5,
            fontweight="bold", color=PURPLE)
    # CO2 (O=C=O)
    cC=(2.05,4.7); oL=(1.35,4.7); oR=(2.75,4.7)
    double_bond(ax,cC,oL); double_bond(ax,cC,oR)
    atom(ax,oL,OXY,label="O"); atom(ax,oR,OXY,label="O"); atom(ax,cC,CARB,label="C")
    ax.text(2.05,4.18,"CO₂  (линейная, инертная)",ha="center",fontsize=10,color=INK)
    # arrow down
    ax.add_patch(FancyArrowPatch((2.05,3.95),(2.05,3.3),arrowstyle="-|>",
                 mutation_scale=18,color=PURPLE,lw=2.2))
    ax.text(2.55,3.62,"2e⁻, 2H⁺",fontsize=9,color=PURPLE)
    # *CO atop Cu
    cu=(2.05,2.55); c1=(2.05,3.05); o1=(2.05,3.5)
    double_bond(ax,c1,o1); bond(ax,cu,c1,lw=2.4,color=COPPER,ls=(0,(1,1)))
    atom(ax,o1,OXY,label="O"); atom(ax,c1,CARB,label="C")
    atom(ax,cu,COPPER,r=0.2,label="Cu")
    ax.text(2.05,2.0,"*CO  (хемосорбция на Cu)",ha="center",fontsize=10,color=INK)

    # ── ЗОНА B: C–C сочленение (узкое место) ────────────────────────────
    ax.text(7.0,5.95,"2 · C–C СОЧЛЕНЕНИЕ  ·  мультиреференсный TS",
            ha="center",fontsize=12.5,fontweight="bold",color=RED)
    # two *CO atop two Cu, C...C approaching
    cuA=(6.1,2.55); cuB=(7.9,2.55)
    cA=(6.35,3.15); cB=(7.65,3.15); oA=(6.0,3.62); oB=(8.0,3.62)
    bond(ax,cuA,cA,lw=2.4,color=COPPER,ls=(0,(1,1)))
    bond(ax,cuB,cB,lw=2.4,color=COPPER,ls=(0,(1,1)))
    double_bond(ax,cA,oA); double_bond(ax,cB,oB)
    # forming C...C bond (dashed red)
    bond(ax,cA,cB,lw=2.6,color=RED,ls=(0,(2,2)),z=4)
    ax.text(7.0,3.46,"C···C",ha="center",fontsize=10.5,color=RED,fontweight="bold")
    atom(ax,oA,OXY,label="O"); atom(ax,oB,OXY,label="O")
    atom(ax,cA,CARB,label="C"); atom(ax,cB,CARB,label="C")
    atom(ax,cuA,COPPER,r=0.2,label="Cu"); atom(ax,cuB,COPPER,r=0.2,label="Cu")
    ax.text(7.0,1.98,"*CO + *CO → *OCCO  (лимитирующая стадия)",ha="center",
            fontsize=10,color=INK)

    # diradicaloid TS callout — литературное число (поверхность Cu(100), 2025)
    noon_txt = "NOON(TS) ≈ 1.92 / 0.09 e  ·  дирадикалоид  (лит., поверхность)"
    box = FancyBboxPatch((5.0,4.55),4.0,0.95,boxstyle="round,pad=0.10,rounding_size=0.10",
                         fc="#fff5f5",ec=RED,lw=1.6,zorder=7)
    ax.add_patch(box)
    ax.text(7.0,5.18,"переходное состояние дирадикалоидно:",ha="center",
            fontsize=9.6,color=RED,fontweight="bold",zorder=8)
    ax.text(7.0,4.82,noon_txt,ha="center",fontsize=9.0,color=INK,zorder=8)
    ax.text(7.0,4.62,"DFT-PBE занижает барьер ~1 эВ (0.48 vs CASPT2 1.43 эВ, лит. 2025)",
            ha="center",fontsize=8.4,color=GREY,zorder=8,style="italic")

    # ── ЗОНА C: продукты C2+ ────────────────────────────────────────────
    ax.text(12.0,5.55,"3 · ПРОДУКТЫ C2+",ha="center",fontsize=12.5,
            fontweight="bold",color=GREEN)
    # ethylene CH2=CH2
    e1=(11.3,4.5); e2=(12.1,4.5)
    double_bond(ax,e1,e2); atom(ax,e1,CARB,label="C"); atom(ax,e2,CARB,label="C")
    ax.text(11.7,4.02,"этилен  CH₂=CH₂",ha="center",fontsize=10,color=INK)
    # ethanol CH3-CH2-OH
    a1=(11.2,2.9); a2=(11.9,2.9); a3=(12.6,3.2)
    bond(ax,a1,a2); bond(ax,a2,a3)
    atom(ax,a1,CARB,label="C"); atom(ax,a2,CARB,label="C"); atom(ax,a3,OXY,label="O")
    ax.text(11.9,2.42,"этанол  CH₃CH₂OH",ha="center",fontsize=10,color=INK)
    ax.text(12.0,1.85,"+ n-пропанол, ацетат …",ha="center",fontsize=9,color=GREY)

    # big arrows between zones
    ax.add_patch(FancyArrowPatch((3.5,3.0),(4.6,3.0),arrowstyle="-|>",
                 mutation_scale=22,color=GREY,lw=2.4))
    ax.add_patch(FancyArrowPatch((9.4,3.0),(10.5,3.0),arrowstyle="-|>",
                 mutation_scale=22,color=GREY,lw=2.4))

    # disclaimer
    ax.text(0.35,0.55,
            "Схема (механизм C–C сочленения, не расчёт). Лимитирующая стадия маршрута — "
            "димеризация *CO–*CO; именно её барьер DFT систематически занижает.",
            fontsize=10,color=GREY,style="italic")
    ax.text(13.9,0.2,"схема · illustrative",ha="right",fontsize=10,
            color=GREY,style="italic")
    fig.tight_layout()
    p=os.path.join(OUT,"06_co2_to_fuels_molecules.png")
    fig.savefig(p,dpi=100,facecolor="white"); plt.close(fig)
    print("wrote",p)


# ════════════════════════════════════════════════════════════════════════
#  ФИГУРА 2 — дескриптор ΔE‡(CASCI−DFT) и вопрос «от обратного»
# ════════════════════════════════════════════════════════════════════════
def fig_descriptor(res):
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(15.0, 7.6), dpi=100,
                                   gridspec_kw=dict(width_ratios=[1.05, 1.0]))
    fig.suptitle("Дескриптор сочленения C–C: расхождение CASCI − DFT(PBE)  —  и вопрос «от обратного»",
                 fontsize=16.5, fontweight="bold", color=INK, y=0.985)

    # ── Панель A: реальный профиль CASCI vs PBE + фронтир-NOON ──────────
    rows = (res or {}).get("occo_scan", {}).get("rows", [])
    if rows:
        R = np.array([r["R"] for r in rows])
        ep = (np.array([r["e_pbe"] for r in rows]) - rows[-1]["e_pbe"]) * 27.2114
        ec = (np.array([r["e_casci"] for r in rows]) - rows[-1]["e_casci"]) * 27.2114
        axA.plot(R, ep, "--", color="#3182ce", lw=2.4, label="DFT (PBE)")
        axA.plot(R, ec, "-", color=RED, lw=2.6, label="CASCI(12e,10o)")
        axA.fill_between(R, ep, ec, color=RED, alpha=0.12)
        axA.annotate("расхождение CASCI − DFT\n(расчёт, газофаза)",
                     xy=(1.9, 4.5), xytext=(2.45, 6.4), fontsize=9.5, color=RED,
                     ha="center", arrowprops=dict(arrowstyle="->", color=RED, lw=1.3))
    axA.set_xlabel("координата сочленения  C···C, Å", fontsize=11.5)
    axA.set_ylabel("энергия (отн. разведённых *CO), эВ", fontsize=11.5)
    axA.set_title("A · профиль на единой геометрии  (def2-SVP)", fontsize=12, color=INK)
    axA.grid(alpha=0.25); axA.legend(loc="lower right", fontsize=10)
    axA.invert_xaxis()

    # honest box on panel A (top-left, empty region)
    box = FancyBboxPatch((0.03, 0.60), 0.66, 0.34, transform=axA.transAxes,
                         boxstyle="round,pad=0.02,rounding_size=0.02",
                         fc="#fffaf0", ec="#dd6b20", lw=1.4, zorder=5)
    axA.add_patch(box)
    front = (res or {}).get("frontier", {})
    nu_lo = front.get("n_u_min", 0.04); nu_hi = front.get("n_u_max", 0.12)
    axA.text(0.05, 0.885,
             "Честно: жёсткий газофазный скан =\nотталкивательная стена, не TS (артефакт).",
             transform=axA.transAxes, fontsize=8.8, color="#dd6b20", fontweight="bold", zorder=6)
    axA.text(0.05, 0.655,
             f"Фронтир связи C–C почти closed-shell:\nn_u = {nu_lo:.2f}→{nu_hi:.2f}, strong=0. Нейтральная\n"
             "газофаза НЕ даёт дирадикалоид → нужен металл.",
             transform=axA.transAxes, fontsize=8.5, color=INK, zorder=6)

    # ── Панель B: Cu vs сплав CuAl — меняет ли квант порядок сайтов ──────
    axB.set_xlim(0, 10); axB.set_ylim(0, 10); axB.axis("off")
    axB.set_title("B · Cu vs сплав CuAl: меняет ли квант порядок сайтов?",
                  fontsize=12, color=INK)
    sites = (res or {}).get("sites", {})
    labels = {"Cu2": "чистый Cu₂", "CuAl": "сплав CuAl"}
    colors = {"Cu2": COPPER, "CuAl": "#6b7cff"}
    if all(k in sites for k in ("Cu2", "CuAl")):
        dp = {t: sites[t]["dEcouple_pbe_eV"] for t in ("Cu2", "CuAl")}
        dc = {t: sites[t]["dEcouple_casci_eV"] for t in ("Cu2", "CuAl")}
        allv = list(dp.values()) + list(dc.values())
        lo, hi = min(allv), max(allv); span = max(hi - lo, 0.4)
        def Y(v): return 2.3 + 4.8 * (v - lo) / span
        xL, xR = 3.2, 6.8
        flipped = sites.get("_summary", {}).get("order_flipped", False)
        verdict = ("КВАНТ МЕНЯЕТ порядок сайтов — DFT-волкан ломается"
                   if flipped else "порядок сайтов сохраняется")
        axB.text(5.0, 9.05, verdict, ha="center", fontsize=11,
                 color=("#c53030" if flipped else "#276749"), fontweight="bold")
        axB.text(xL, 7.9, "DFT (PBE)", ha="center", fontsize=10.5, color="#3182ce", fontweight="bold")
        axB.text(xR, 7.9, "CASCI(12e,10o)", ha="center", fontsize=10.5, color=RED, fontweight="bold")
        for t in ("Cu2", "CuAl"):
            axB.plot([xL, xR], [Y(dp[t]), Y(dc[t])], "-o", color=colors[t], lw=2.6, ms=11, zorder=4)
            axB.text(xL - 0.35, Y(dp[t]), f"{dp[t]:+.2f}", ha="right", va="center",
                     fontsize=9, color=colors[t])
            axB.text(xR + 0.35, Y(dc[t]), f"{labels[t]}  {dc[t]:+.2f} эВ", ha="left",
                     va="center", fontsize=9.5, color=colors[t], fontweight="bold")
        dcorr = {t: sites[t].get("quantum_corr_eV") for t in ("Cu2", "CuAl")}
        axB.text(5.0, 1.55,
                 f"мультиреференс-поправка ΔE_couple(CASCI−PBE):  Cu₂ {dcorr['Cu2']:+.2f} эВ,  "
                 f"CuAl {dcorr['CuAl']:+.2f} эВ",
                 ha="center", fontsize=9, color=INK)
        axB.text(5.0, 0.95,
                 "ось ↑ = ΔE_couple (E сочленённые − разведённые); ниже = легче сочленение.",
                 ha="center", fontsize=8, color=GREY)
    else:
        axB.text(5.0, 5.0, "ВОПРОС ОТ ОБРАТНОГО (расчёт идёт):\n"
                 "меняет ли ΔE‡(CASCI−DFT) порядок Cu vs CuAl?", ha="center",
                 fontsize=10, color="#276749", fontweight="bold")
    # compact literature anchor reference
    axB.text(5.0, 3.7, "лит. якорь (поверхность Cu(100)): барьер PBE 0.48 vs CASPT2 1.43 эВ; "
             "NOON(TS) 1.92/0.09 — дирадикалоид (2025)", ha="center", fontsize=8,
             color=GREY, style="italic")

    # footer
    fig.text(0.012, 0.012,
             "Панель B: ΔE_couple — РЕАКЦИОННАЯ энергия сочленения (релакс-эндпоинты, замороженный "
             "металл-кластер, def2-SVP), НЕ барьер ΔE‡ и НЕ фарадеевская эффективность (кластер ≠ "
             "электрод). Седло + constant-potential — следующий прогон на AWS. Контроль метрики: O₂ excess +0.09.",
             fontsize=8.3, color=GREY, style="italic")
    fig.text(0.988, 0.012, "схема · отчасти illustrative", ha="right", fontsize=10,
             color=GREY, style="italic")
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    p = os.path.join(OUT, "06_co2_to_fuels_descriptor.png")
    fig.savefig(p, dpi=100, facecolor="white"); plt.close(fig)
    print("wrote", p)


if __name__ == "__main__":
    res = _load_results()
    fig_molecules(res)
    fig_descriptor(res)
    print("results loaded:", bool(res))
