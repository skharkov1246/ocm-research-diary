#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MatterForge · вкладка h2o2-direct — генератор фигур (по образцу routes/polyolefins.py).

Две ИЛЛЮСТРАТИВНЫЕ схемы (не расчёт):
  07_h2o2_direct_molecules.png   — реагенты H2 + O2 (триплет) + развилка селективности на Pd-центре
  07_h2o2_direct_descriptor.png  — дескриптор-ландшафт: связать *OOH, НЕ разорвав O-O (2e- -> H2O2)

Каждая фигура несёт дисклеймер «схема · illustrative» — позиции отражают
литературные тренды (Pd1 single-atom, PdAu/PdSn), а не расчёт этого проекта.

Реальные числа (мультиреференсность O2 / разрыва O-O) живут в дневнике (updates[]),
а воспроизводятся функцией multireference_probe() ниже (PySCF, метрика NOON как в OCM/H2).

Запуск:  python3 routes/h2o2_direct.py            # генерит обе фигуры
         python3 routes/h2o2_direct.py --probe    # печатает реальный зонд (нужен pyscf)
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

OUT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "assets", "matterforge"))

INK   = "#1a202c"
GRAY  = "#718096"
PURP  = "#6b46c1"
ACC   = "#00b5d8"   # цвет вкладки h2o2-direct
GREEN = "#2f855a"
RED   = "#c53030"
OXY   = "#c53030"
HYD   = "#4a5568"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "savefig.dpi": 150,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})


# --------------------------------------------------------------------------- #
#  низкоуровневые рисовалки молекул (кружок = атом, линия = связь)
# --------------------------------------------------------------------------- #
def atom(ax, x, y, label, r=0.20, fc="white", ec=INK, fs=12, tc=INK, lw=2.0):
    ax.add_patch(Circle((x, y), r, facecolor=fc, edgecolor=ec, lw=lw, zorder=4))
    ax.text(x, y, label, ha="center", va="center", fontsize=fs, color=tc,
            fontweight="bold", zorder=5)


def bond(ax, p1, p2, double=False, color=INK, lw=2.2, off=0.05):
    (x1, y1), (x2, y2) = p1, p2
    if not double:
        ax.plot([x1, x2], [y1, y2], color=color, lw=lw, zorder=2,
                solid_capstyle="round")
        return
    dx, dy = x2 - x1, y2 - y1
    n = (dy, -dx)
    L = (n[0] ** 2 + n[1] ** 2) ** 0.5
    ox, oy = off * n[0] / L, off * n[1] / L
    for s in (+1, -1):
        ax.plot([x1 + s * ox, x2 + s * ox], [y1 + s * oy, y2 + s * oy],
                color=color, lw=lw, zorder=2, solid_capstyle="round")


def arrow(ax, p1, p2, color=PURP, lw=2.6, ms=18):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=ms,
                                 lw=lw, color=color, zorder=3,
                                 shrinkA=2, shrinkB=2))


def textbox(ax, x, y, s, fc="#fff5f5", ec=RED, tc=RED, fs=10.5, ha="center",
            va="center", w="bold"):
    ax.text(x, y, s, ha=ha, va=va, fontsize=fs, color=tc, fontweight=w,
            zorder=6, bbox=dict(boxstyle="round,pad=0.45", fc=fc, ec=ec, lw=1.6))


def illustrative_tag(fig):
    fig.text(0.985, 0.012, "схема · illustrative", ha="right", va="bottom",
             fontsize=11, style="italic", color=GRAY)


# --------------------------------------------------------------------------- #
#  ФИГУРА 1 — молекулы / механизм
# --------------------------------------------------------------------------- #
def fig_molecules(path):
    fig, ax = plt.subplots(figsize=(14.2, 7.1))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 9)
    ax.axis("off")

    ax.set_title(
        "Прямой синтез H₂O₂:  реагенты + развилка селективности на Pd-центре"
        "  (иллюстративно)",
        fontsize=15.5, fontweight="bold", color=INK, pad=14)

    # колоночные заголовки
    ax.text(2.6, 8.25, "РЕАГЕНТЫ", ha="center", fontsize=13.5,
            fontweight="bold", color=PURP)
    ax.text(7.5, 8.25, "АКТИВНЫЙ ЦЕНТР", ha="center", fontsize=13.5,
            fontweight="bold", color=PURP)
    ax.text(7.5, 7.78, "M = Pd₁ (single-atom) · PdAu · PdSn", ha="center",
            fontsize=11, color=GRAY)
    ax.text(12.8, 8.25, "РАЗВИЛКА СЕЛЕКТИВНОСТИ", ha="center", fontsize=13.5,
            fontweight="bold", color=PURP)

    # --- реагенты ---
    # H2
    ax.text(2.6, 6.95, "H₂", ha="center", fontsize=12, color=INK, fontweight="bold")
    atom(ax, 2.05, 6.4, "H", r=0.17, ec=HYD, tc=HYD, fs=10)
    atom(ax, 3.15, 6.4, "H", r=0.17, ec=HYD, tc=HYD, fs=10)
    bond(ax, (2.22, 6.4), (2.98, 6.4), color=HYD, lw=2.2)

    # O2 — триплет, два неспаренных e- на π*
    ax.text(2.6, 4.95, "O₂  (³Σg⁻ триплет)", ha="center", fontsize=12,
            color=INK, fontweight="bold")
    atom(ax, 2.05, 4.3, "O", r=0.22, ec=OXY, tc=OXY)
    atom(ax, 3.15, 4.3, "O", r=0.22, ec=OXY, tc=OXY)
    bond(ax, (2.27, 4.3), (2.93, 4.3), double=True, color=OXY, lw=2.2, off=0.07)
    # спиновые стрелки ↑↑
    for sx in (2.05, 3.15):
        ax.annotate("", xy=(sx, 4.92), xytext=(sx, 4.62),
                    arrowprops=dict(arrowstyle="-|>", color=OXY, lw=2.0))
    ax.text(2.6, 3.65, "2 неспаренных e⁻ на π*  (открытая оболочка)",
            ha="center", fontsize=9.6, color=GRAY, style="italic")

    arrow(ax, (3.7, 5.2), (5.3, 5.4), color=ACC)

    # --- активный центр: Pd + *OOH (O-O ЦЕЛА) ---
    atom(ax, 6.6, 5.4, "Pd", r=0.34, fc=PURP, ec=PURP, tc="white", fs=12, lw=2)
    # *OOH: Pd—O—O—H, связь O-O сохранена
    atom(ax, 7.55, 5.85, "O", r=0.21, ec=OXY, tc=OXY)
    atom(ax, 8.35, 5.4, "O", r=0.21, ec=OXY, tc=OXY)
    atom(ax, 9.05, 5.85, "H", r=0.16, ec=HYD, tc=HYD, fs=9.5)
    bond(ax, (6.92, 5.55), (7.36, 5.78), color=GRAY, lw=1.8)   # Pd···O
    bond(ax, (7.74, 5.78), (8.18, 5.49), color=OXY, lw=2.2)    # O–O (цела)
    bond(ax, (8.53, 5.5), (8.92, 5.74), color=HYD, lw=1.8)     # O–H
    ax.text(8.05, 6.5, "*OOH связан,  O–O ЦЕЛА", ha="center", fontsize=10.5,
            color=GREEN, fontweight="bold")
    textbox(ax, 7.5, 4.05,
            "КЛЮЧ: связать *OOH без\nдиссоциативной адсорбции O₂\n(не рвать O–O)",
            fc="#f0fff4", ec=GREEN, tc=GREEN, fs=10)

    # стрелка в развилку
    arrow(ax, (9.4, 5.4), (10.5, 5.4), color=INK)

    # --- развилка ---
    # верх: 2e- -> H2O2  (хорошо)
    ax.add_patch(FancyBboxPatch((10.7, 6.0), 4.9, 1.55,
                 boxstyle="round,pad=0.08", fc="#f0fff4", ec=GREEN, lw=2,
                 zorder=1))
    ax.text(11.0, 7.25, "2e⁻  →  H₂O₂", ha="left", fontsize=12.5,
            color=GREEN, fontweight="bold")
    # H-O-O-H (O-O цела)
    atom(ax, 11.6, 6.55, "H", r=0.14, ec=HYD, tc=HYD, fs=9)
    atom(ax, 12.25, 6.55, "O", r=0.19, ec=OXY, tc=OXY, fs=11)
    atom(ax, 12.95, 6.55, "O", r=0.19, ec=OXY, tc=OXY, fs=11)
    atom(ax, 13.6, 6.55, "H", r=0.14, ec=HYD, tc=HYD, fs=9)
    bond(ax, (11.73, 6.55), (12.07, 6.55), color=HYD, lw=1.7)
    bond(ax, (12.44, 6.55), (12.76, 6.55), color=OXY, lw=2.4)
    bond(ax, (13.14, 6.55), (13.47, 6.55), color=HYD, lw=1.7)
    ax.text(14.05, 6.55, "✓ O–O\nсохранена", ha="left", va="center",
            fontsize=9.3, color=GREEN, fontweight="bold")

    # низ: 4e- -> 2 H2O (плохо)
    ax.add_patch(FancyBboxPatch((10.7, 3.5), 4.9, 1.55,
                 boxstyle="round,pad=0.08", fc="#fff5f5", ec=RED, lw=2,
                 zorder=1))
    ax.text(11.0, 4.75, "4e⁻  →  2 H₂O", ha="left", fontsize=12.5,
            color=RED, fontweight="bold")
    # две молекулы воды (O-O разорвана)
    for cx in (12.1, 13.4):
        atom(ax, cx, 4.05, "O", r=0.19, ec=OXY, tc=OXY, fs=11)
        atom(ax, cx - 0.5, 3.85, "H", r=0.13, ec=HYD, tc=HYD, fs=8.5)
        atom(ax, cx + 0.5, 3.85, "H", r=0.13, ec=HYD, tc=HYD, fs=8.5)
        bond(ax, (cx - 0.16, 3.98), (cx - 0.37, 3.9), color=HYD, lw=1.6)
        bond(ax, (cx + 0.16, 3.98), (cx + 0.37, 3.9), color=HYD, lw=1.6)
    # «разрыв» между двумя O
    ax.plot([12.55, 12.95], [4.05, 4.05], color=RED, lw=1.4, ls=(0, (2, 2)))
    ax.text(12.75, 4.32, "✗ O–O\nразорвана", ha="center", va="bottom",
            fontsize=8.6, color=RED, fontweight="bold")

    # стрелки развилки
    arrow(ax, (10.55, 5.4), (10.85, 6.55), color=GREEN, lw=2.4, ms=15)
    arrow(ax, (10.55, 5.4), (10.85, 4.3), color=RED, lw=2.4, ms=15)

    textbox(ax, 12.9, 2.55,
            "+ паразитные пути: разложение и\nгидрирование уже готового H₂O₂ → H₂O",
            fc="#fffaf0", ec="#dd6b20", tc="#c05621", fs=9.6)

    fig.text(0.012, 0.028,
             "Схема (механизм, не расчёт).  Вся борьба — удержать 2e⁻-путь: "
             "связать *OOH и НЕ разорвать O–O.  O₂-триплет и супероксо/пероксо "
             "*OOH — открытая оболочка (реальный зонд мультиреференсности — в дневнике).",
             ha="left", va="bottom", fontsize=10, style="italic", color=GRAY)
    illustrative_tag(fig)

    fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.085)
    fig.savefig(path)
    plt.close(fig)
    print("wrote", path)


# --------------------------------------------------------------------------- #
#  ФИГУРА 2 — дескриптор-ландшафт
# --------------------------------------------------------------------------- #
def fig_descriptor(path):
    fig, ax = plt.subplots(figsize=(13.8, 8.4))

    ax.set_title(
        "Дескриптор селективности H₂O₂:  связать *OOH, НЕ разорвав O–O\n"
        "(иллюстративно, не расчёт)",
        fontsize=15, fontweight="bold", color=INK, pad=12)

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_xlabel("склонность к разрыву O–O  (диссоциативная адсорбция O₂)  →",
                  fontsize=12.5)
    ax.set_ylabel("прочность связывания *OOH  (2e⁻-дескриптор)  →", fontsize=12.5)
    ax.grid(True, ls=":", color="#cbd5e0", alpha=0.8)
    ax.set_axisbelow(True)

    # зона цели: низкий разрыв O-O + умеренно-сильное *OOH -> высокая H2O2-селективность
    ax.add_patch(Rectangle((0.4, 6.2), 3.7, 3.2, facecolor=ACC, alpha=0.12,
                           edgecolor=ACC, lw=1.4, ls="--", zorder=0))
    ax.text(2.25, 9.05, "зона высокой H₂O₂-селективности  (2e⁻)",
            ha="center", fontsize=11.5, color="#086f83", fontweight="bold")
    ax.text(0.5, 0.55, "←  держит O–O  (→ H₂O₂)", fontsize=10.5, color=GREEN,
            fontweight="bold")
    ax.text(9.5, 0.55, "рвёт O–O  (→ H₂O)  →", fontsize=10.5, color=RED,
            fontweight="bold", ha="right")

    # литературные точки (тренды, НЕ расчёт): (x=разрыв O-O, y=*OOH binding)
    pts = [
        # label, x, y, yerr, color, dx, dy, ha
        ("изолированный Pd₁\n(single-atom)", 1.5, 8.1, 0.5, "#2b6cb0", 0.35, 0.0, "left"),
        ("PdAu (разбавл. Pd)", 2.2, 7.3, 0.55, "#2f855a", 0.35, -0.1, "left"),
        ("PdSn", 1.9, 6.7, 0.5, "#38a169", 0.35, -0.05, "left"),
        ("чистый Pd\n(наночастица)", 5.6, 6.4, 0.6, "#805ad5", 0.35, 0.15, "left"),
        ("Pt (рвёт O–O)", 8.4, 5.2, 0.65, "#c05621", -0.4, 0.1, "right"),
    ]
    for lab, x, y, ye, c, dx, dy, ha in pts:
        ax.errorbar(x, y, yerr=ye, fmt="o", ms=15, color=c, ecolor=c,
                    elinewidth=1.6, capsize=4, mec="white", mew=1.4, zorder=5)
        ax.text(x + dx, y + dy, lab, fontsize=10.5, color=c, fontweight="bold",
                ha=ha, va="center", zorder=6)

    # ЦЕЛЬ дизайна: дешёвый 3d-центр с тем же дескриптором, что Pd1
    ax.scatter([2.6], [7.7], marker="*", s=620, color=ACC, edgecolor="#086f83",
               linewidths=1.8, zorder=7)
    ax.annotate("ЦЕЛЬ: дешёвый 3d-центр\n(Co/Ni/Fe–Nₓ single-atom)\n"
                "с тем же дескриптором, что Pd₁",
                xy=(2.6, 7.7), xytext=(4.6, 8.7), fontsize=10.5, color="#086f83",
                fontweight="bold", ha="left", va="center",
                arrowprops=dict(arrowstyle="-|>", color="#086f83", lw=1.8),
                bbox=dict(boxstyle="round,pad=0.4", fc="#e6fffa", ec=ACC, lw=1.5))

    # квант-edge (качественно — физика, не числа этого графика)
    textbox(ax, 6.95, 2.95,
            "квант-edge:  порядок связи O–O и спин на *OOH —\n"
            "мультиреференсны (O₂-триплет; супероксо ↔ пероксо).\n"
            "Одно-детерминантный DFT тут систематически ненадёжен →\n"
            "цель — многоконфигурационная (квантовая) точность.",
            fc="#f7fafc", ec=PURP, tc=PURP, fs=10, ha="center")

    # льюисова кислотность
    ax.text(1.35, 4.4, "льюисова кислотность\nсреды помогает 2e⁻-пути",
            fontsize=9.6, color=GRAY, style="italic", ha="left")
    ax.annotate("", xy=(1.7, 6.3), xytext=(1.9, 5.0),
                arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=1.3))

    # честный caveat (инженерия/экономика)
    textbox(ax, 6.6, 0.85,
            "Честно: дескриптор доказан для ORR-электрокатализа — перенос на "
            "термическую H₂+O₂ надо проверить.\n"
            "Главная стена — взрывное окно смеси H₂/O₂ (~4–96 %), и её квант НЕ двигает.",
            fc="#fffaf0", ec="#dd6b20", tc="#c05621", fs=9.8, ha="center")

    fig.text(0.012, 0.012,
             "Схема — позиции отражают известные литературные тренды "
             "(Pd₁ / PdAu / PdSn), НЕ расчёт этого проекта.  Та же дескриптор-логика, "
             "что в наших картах ΔΔG (OCM / полиолефины).",
             ha="left", va="bottom", fontsize=9.6, style="italic", color=GRAY)
    illustrative_tag(fig)

    fig.subplots_adjust(left=0.075, right=0.985, top=0.9, bottom=0.105)
    fig.savefig(path)
    plt.close(fig)
    print("wrote", path)


# --------------------------------------------------------------------------- #
#  РЕАЛЬНЫЙ зонд мультиреференсности (PySCF) — воспроизводит числа из дневника.
#  Метрика N_u = Σ nᵢ(2−nᵢ) по заселённостям натуральных орбиталей (базис-инвариантна,
#  как в OCM / H₂ / полиолефинах). Фигуры от него НЕ зависят (они иллюстративны).
# --------------------------------------------------------------------------- #
def multireference_probe():
    import numpy as np
    from pyscf import gto, scf, mcscf
    from pyscf.mcscf import avas

    def Nu(occ):
        occ = np.asarray(occ)
        return float(np.sum(occ * (2.0 - occ)))

    def cas_noon(mc):
        dm1 = mc.fcisolver.make_rdm1(mc.ci, mc.ncas, mc.nelecas)
        return np.linalg.eigvalsh(dm1)[::-1]

    print("== O2 основное состояние: триплет-бирадикал? ==")
    o2 = gto.M(atom="O 0 0 0; O 0 0 1.207", basis="def2-svp", spin=2,
               symmetry=False, verbose=0)
    mf = scf.ROHF(o2).run()
    mc = mcscf.CASSCF(mf, 6, 8); mc.verbose = 0; mc.kernel()
    occ = cas_noon(mc)
    print(f"  O2 триплет CAS(8e,6o)/def2-SVP  NOON={np.round(occ,3)}  N_u={Nu(occ):.2f}")

    print("== судьба связи O–O: цела (H2O2, 2e-) vs разрыв (H2O, 4e-) ==")

    def h2o2(Roo, roh=0.965, aooh=100.0, dih=113.7):
        a, f = np.deg2rad(aooh), np.deg2rad(dih)
        O1 = np.array([0, 0, 0.0]); O2 = np.array([0, 0, Roo])
        H1 = O1 + roh * np.array([np.sin(a), 0, np.cos(a)])
        H2 = O2 + roh * np.array([np.sin(np.pi - a) * np.cos(f),
                                  np.sin(np.pi - a) * np.sin(f),
                                  np.cos(np.pi - a)])
        at = [("O", O1), ("O", O2), ("H", H1), ("H", H2)]
        return "; ".join(f"{s} {p[0]:.4f} {p[1]:.4f} {p[2]:.4f}" for s, p in at)

    for Roo, tag in [(1.452, "равновесие (продукт H2O2)"),
                     (1.90, "растянута"),
                     (2.40, "к гомолизу O–O (цель 4e-)")]:
        m = gto.M(atom=h2o2(Roo), basis="def2-svp", spin=0, symmetry=False, verbose=0)
        mf = scf.RHF(m); mf.kernel()
        ncas, nelecas, mo = avas.avas(mf, ["O 2p"], canonicalize=True, verbose=0)
        mc = mcscf.CASSCF(mf, ncas, nelecas); mc.verbose = 0
        mc.max_cycle_macro = 100; mc.fix_spin_(ss=0); mc.kernel(mo)
        occ = cas_noon(mc)
        print(f"  R(O-O)={Roo:.2f}Å [{tag}] CAS({nelecas}e,{ncas}o) "
              f"{'conv' if mc.converged else 'NOTconv'}  N_u={Nu(occ):.2f}  "
              f"NOON={np.round(occ,3)}")

    print("== интермедиаты 2e--пути: супероксо O2^-• и гидропероксил HOO• (=*OOH) ==")
    om = gto.M(atom="O 0 0 0; O 0 0 1.28", basis="def2-svp", spin=1, charge=-1,
               symmetry=False, verbose=0)
    mf = scf.ROHF(om); mf.kernel()
    mc = mcscf.CASSCF(mf, 6, 9); mc.verbose = 0; mc.kernel()
    occ = cas_noon(mc)
    print(f"  супероксо O2^-• CAS(9e,6o) N_u={Nu(occ):.2f}  NOON={np.round(occ,3)}  "
          f"(газофазный анион в def2-SVP без диффузных — для тренда)")

    a = np.deg2rad(104.3); roh = 0.971; Roo = 1.331
    O2 = np.array([0, 0, Roo])
    H = O2 + roh * np.array([np.sin(np.pi - a), 0, np.cos(np.pi - a)])
    hoo = f"O 0 0 0; O 0 0 {Roo}; H {H[0]:.4f} {H[1]:.4f} {H[2]:.4f}"
    m = gto.M(atom=hoo, basis="def2-svp", spin=1, charge=0, symmetry=False, verbose=0)
    mf = scf.ROHF(m); mf.kernel()
    ncas, nelecas, mo = avas.avas(mf, ["O 2p"], canonicalize=True, verbose=0)
    mc = mcscf.CASSCF(mf, ncas, nelecas); mc.verbose = 0; mc.kernel(mo)
    occ = cas_noon(mc)
    print(f"  HOO• (*OOH) CAS({nelecas}e,{ncas}o) N_u={Nu(occ):.2f}  NOON={np.round(occ,3)}")


# --------------------------------------------------------------------------- #
def main():
    if "--probe" in sys.argv:
        multireference_probe()
        return
    os.makedirs(OUT, exist_ok=True)
    fig_molecules(os.path.join(OUT, "07_h2o2_direct_molecules.png"))
    fig_descriptor(os.path.join(OUT, "07_h2o2_direct_descriptor.png"))


if __name__ == "__main__":
    main()
