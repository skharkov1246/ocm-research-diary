#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cc-activation-recycling — schematic figures (REVERSE mirror of routes/polyolefins.py).

Generates two illustrative schematics (NOT data plots):
  08_cc_activation_recycling_molecules.png   -- the reverse mechanism: a late-TM pincer
       does oxidative addition (OA) into an inert C(sp3)-C(sp3) bond of polyethylene.
  08_cc_activation_recycling_descriptor.png  -- descriptor landscape: barrier E-dagger
       vs thermodynamic tilt dG(OA), target zone, literature anchors, inversion insight.

Every figure carries a "схема · illustrative" disclaimer. Literature anchors
(57 kJ/mol Rh-POCOP; 29 C Brookhart; 125-200 C tandem; Conk/Grubbs 87 %) are drawn
AS literature, never as our own number. If routes/cc_activation_recycling_results.json
exists, the real computed N_u(sigma->TS) is added as an honest annotation.
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, '..', 'assets', 'matterforge'))
ACCENT = '#9f7aea'      # tab color (purple)
ACCENT_D = '#6b46c1'
INK = '#2d3748'
MUTE = '#718096'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'svg.fonttype': 'none'})

RESULTS = os.path.join(HERE, 'cc_activation_recycling_results.json')
def load_results():
    if os.path.exists(RESULTS):
        try:
            return json.load(open(RESULTS))
        except Exception:
            return None
    return None

def disclaimer(fig, text):
    fig.text(0.012, 0.016, text, fontsize=9.6, style='italic', color=MUTE, ha='left',
             wrap=True)
    fig.text(0.988, 0.016, 'схема · illustrative', fontsize=9.6, style='italic',
             color=MUTE, ha='right')

# ---------------------------------------------------------------------------
def atom(ax, xy, label, r=0.16, fc='white', ec=INK, fs=12, tc=INK, lw=2.0):
    ax.add_patch(Circle(xy, r, fc=fc, ec=ec, lw=lw, zorder=3))
    ax.text(xy[0], xy[1], label, ha='center', va='center', fontsize=fs,
            color=tc, zorder=4, fontweight='bold')

def bond(ax, p, q, lw=2.4, color=INK, ls='-', z=2):
    ax.plot([p[0], q[0]], [p[1], q[1]], ls, lw=lw, color=color, zorder=z,
            solid_capstyle='round')

def dbond(ax, p, q, lw=2.4, color=INK, off=0.05):
    d = np.array(q) - np.array(p); n = np.array([-d[1], d[0]]); n = n/np.linalg.norm(n)*off
    bond(ax, np.array(p)+n, np.array(q)+n, lw=lw, color=color)
    bond(ax, np.array(p)-n, np.array(q)-n, lw=lw, color=color)

def arrow(ax, p, q, color=ACCENT, lw=3.0, mut=22):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle='-|>', mutation_scale=mut,
                 lw=lw, color=color, zorder=5))

# ===========================================================================
def fig_molecules():
    fig, ax = plt.subplots(figsize=(14.2, 6.3))
    ax.set_xlim(0, 14.2); ax.set_ylim(0, 6.3); ax.axis('off')

    fig.suptitle('C–C-активация: разрыв инертной цепи ПЭ late-TM пинцером — '
                 'reverse полиолефинов  (иллюстративно)',
                 fontsize=16, fontweight='bold', y=0.975, color='#1a202c')

    # column headers
    ax.text(2.45, 5.5, 'СУБСТРАТ (ПЭ)', fontsize=14, fontweight='bold', color=ACCENT_D, ha='center')
    ax.text(7.1, 5.5, 'late-TM пинцер   M = Rh / Ir', fontsize=14, fontweight='bold',
            color=ACCENT_D, ha='center')
    ax.text(11.75, 5.5, 'ФРАГМЕНТЫ', fontsize=14, fontweight='bold', color=ACCENT_D, ha='center')

    # --- LEFT: polyethylene chain, one inert C-C bond highlighted ---
    ax.text(2.45, 4.85, 'полиэтилен  …–CH₂–CH₂–…', fontsize=12, color=INK, ha='center')
    ys = 3.7
    xs = [0.7, 1.45, 2.2, 2.95, 3.7]
    pts = [(x, ys + (0.28 if i % 2 else -0.28)) for i, x in enumerate(xs)]
    for i in range(len(pts) - 1):
        col = '#e53e3e' if i == 2 else INK
        lw = 3.4 if i == 2 else 2.4
        bond(ax, pts[i], pts[i+1], lw=lw, color=col)
    for p in pts:
        atom(ax, p, 'C', r=0.15, fs=10)
    ax.annotate('инертная C(sp³)–C(sp³)\nBDE ≈ 85 ккал/моль · нет диполя/π',
                xy=pts[2], xytext=(2.2, 2.05), ha='center', fontsize=10.5, color='#c53030',
                arrowprops=dict(arrowstyle='->', color='#c53030', lw=1.6))
    ax.text(2.45, 1.25, 'нечего «ухватить» →\nсегодня: пиролиз 400–600 °C,\nрадикально, БЕЗ селективности',
            fontsize=10, color=MUTE, ha='center')

    # --- CENTER: pincer metal forming sigma-complex then OA ---
    M = (7.1, 3.55)
    # two P arms (pincer)
    p1 = (6.35, 4.35); p2 = (6.35, 2.75)
    bond(ax, M, p1, lw=2.4, color=ACCENT_D); bond(ax, M, p2, lw=2.4, color=ACCENT_D)
    atom(ax, p1, 'P', r=0.17, fc='#edf2f7', ec=ACCENT_D, tc=ACCENT_D, fs=11)
    atom(ax, p2, 'P', r=0.17, fc='#edf2f7', ec=ACCENT_D, tc=ACCENT_D, fs=11)
    ax.text(5.75, 3.55, 'пинцер\n(P–C–P)', fontsize=9.5, color=ACCENT_D, ha='center', va='center')
    atom(ax, M, 'M', r=0.30, fc=ACCENT, ec=ACCENT_D, tc='white', fs=14)
    # eta2 C-C approaching (sigma-complex), dotted to M
    c1 = (8.05, 3.95); c2 = (8.05, 3.15)
    bond(ax, c1, c2, lw=2.6, color=INK)
    atom(ax, c1, 'C', r=0.15, fs=10); atom(ax, c2, 'C', r=0.15, fs=10)
    ax.plot([M[0]+0.30, c1[0]-0.12], [M[1]+0.10, c1[1]-0.04], ':', lw=2.0, color=MUTE, zorder=2)
    ax.plot([M[0]+0.30, c2[0]-0.12], [M[1]-0.10, c2[1]+0.04], ':', lw=2.0, color=MUTE, zorder=2)
    ax.text(7.1, 2.15, 'σ-комплекс  M···(C–C)\nкинетика: τ(σ) — мишень дизайна',
            fontsize=10, color=INK, ha='center')
    ax.text(7.1, 1.35, 'M(I) d⁸ → OA → M(III) d⁶\n(разрыв σ + смена с.о. металла)',
            fontsize=9.5, color=MUTE, ha='center')

    # --- arrow OA ---
    arrow(ax, (8.9, 3.55), (10.2, 3.55), color=ACCENT, lw=3.2, mut=24)
    ax.text(9.55, 3.95, 'OA', fontsize=13, fontweight='bold', color=ACCENT_D, ha='center')
    ax.text(9.55, 3.2, 'окислит.\nприсоединение', fontsize=9, color=MUTE, ha='center')
    ax.text(9.55, 4.5, '(reverse внедрения)', fontsize=9, color=MUTE, ha='center', style='italic')

    # --- RIGHT: two M-C fragments ---
    M2 = (11.4, 3.55)
    pr1 = (10.75, 4.3); pr2 = (10.75, 2.8)
    bond(ax, M2, pr1, lw=2.2, color=ACCENT_D); bond(ax, M2, pr2, lw=2.2, color=ACCENT_D)
    atom(ax, pr1, 'P', r=0.15, fc='#edf2f7', ec=ACCENT_D, tc=ACCENT_D, fs=10)
    atom(ax, pr2, 'P', r=0.15, fc='#edf2f7', ec=ACCENT_D, tc=ACCENT_D, fs=10)
    atom(ax, M2, 'M', r=0.28, fc=ACCENT, ec=ACCENT_D, tc='white', fs=13)
    cf1 = (12.25, 4.2); cf2 = (12.25, 2.9)
    bond(ax, M2, cf1, lw=2.6, color=ACCENT_D); bond(ax, M2, cf2, lw=2.6, color=ACCENT_D)
    atom(ax, cf1, 'C', r=0.15, fs=10); atom(ax, cf2, 'C', r=0.15, fs=10)
    ax.text(11.95, 1.95, '2 × M–C\n(цепь разрезана)', fontsize=10.5, color=INK, ha='center')
    ax.text(11.95, 1.2, '→ короткие фрагменты /\nмономер (цель)', fontsize=10, color=MUTE, ha='center')

    disclaimer(fig, 'Схема (механизм, не расчёт). Reverse-рука polyolefins: разрыв C–C вместо '
                    'внедрения мономера. Риск: прямой OA эндотермичен → тандем.')
    fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.1)
    path = os.path.join(OUT, '08_cc_activation_recycling_molecules.png')
    fig.savefig(path, dpi=130); plt.close(fig)
    print('wrote', path)

# ===========================================================================
def fig_descriptor():
    res = load_results()
    fig, ax = plt.subplots(figsize=(14.6, 9.0))
    fig.suptitle('Descriptor landscape: барьер C–C OA × термодинамический наклон —\n'
                 'мишень дизайна не «горячее», а τ(σ)>τ(TS) при ΔG(OA)≈0  (иллюстративно, не расчёт)',
                 fontsize=15.5, fontweight='bold', y=0.985, color='#1a202c')

    ax.set_xlim(-34, 48); ax.set_ylim(-5, 45)
    ax.set_xlabel('термодинамика OA:  ΔG(OA), ккал/моль   (← экзо · термонейтрально · эндо →)',
                  fontsize=12.5)
    ax.set_ylabel('кинетический барьер  E‡(C–C OA), ккал/моль', fontsize=12.5)
    ax.grid(True, ls=':', color='#cbd5e0', alpha=0.6)
    ax.axvline(0, color=MUTE, lw=1.1, ls='--')
    ax.text(0, 30.5, 'ΔG≈0', fontsize=10, color=MUTE, ha='center')

    # --- target zone: low barrier (<=20) AND near-thermoneutral ---
    ax.add_patch(Rectangle((-13, 0), 26, 21, fc=ACCENT, ec='none', alpha=0.12, zorder=0))
    ax.text(0, 22.7, 'зона цели:  E‡≲20  ·  ΔG(OA)≈0  ·  τ(σ)>τ(TS-OA)',
            fontsize=11.5, fontweight='bold', color=ACCENT_D, ha='center')

    # --- literature anchors (drawn AS literature) ---
    def pt(x, y, color, ms=14):
        ax.plot([x], [y], 'o', ms=ms, color=color, mec='white', mew=1.6, zorder=6)
    def lab(x, y, text, color, ha='center'):
        ax.text(x, y, text, fontsize=9.8, fontweight='bold', color=color, ha=ha,
                va='center', zorder=6)
    # Rh-POCOP: intrinsic barrier MILD (57 kJ/mol=13.6 kcal) WHEN geometrically captured
    pt(4, 13.6, '#2b6cb0'); lab(4, 17.6, 'Rh-POCOP (захвачен)\nлит. ΔG‡≈57 кДж/моль (~13.6)', '#2b6cb0')
    # Brookhart Ir-PONOP: RE of CH4 at 29 C -> by microreversibility mild OA reachable
    pt(-9, 8.5, '#38a169')
    ax.annotate('Brookhart Ir-PONOP\nлит. RE CH₄ при 29 °C →\nмягкий OA (микрообрат.)',
                xy=(-9, 8.5), xytext=(-22, 15), fontsize=9.6, fontweight='bold',
                color='#2f855a', ha='center', zorder=6,
                arrowprops=dict(arrowstyle='-', color='#38a169', lw=1.0, alpha=0.6))
    # acyclic analog: ZERO OA (sigma-complex dissociates before OA)
    pt(40, 25, '#a0aec0'); lab(38, 20.8, 'ациклич. аналог\nлит. ≈0 OA\n(σ уходит раньше)', '#718096')
    # pyrolysis: brutal, radical, no selectivity
    pt(41, 13, '#e53e3e'); lab(33.5, 10.2, 'пиролиз 400–600 °C\nрадикально, БЕЗ селект.', '#c53030')

    # --- TANDEM box (top-left): the industrial bypass of C-C OA ---
    ax.add_patch(FancyBboxPatch((-33, 34.6), 29.5, 9.6, boxstyle='round,pad=0.3',
                 fc='#fffaf0', ec='#dd6b20', lw=1.6, zorder=4))
    ax.text(-18.2, 39.4, 'ТАНДЕМ (обход C–C OA):\nGoldman–Brookhart deH–метатезис 125–200 °C;\n'
            'Conk/Grubbs ПЭ→пропилен, выход 87 %\n→ через C–H / олефин, НЕ через C–C OA',
            fontsize=9.6, color='#9c4221', ha='center', va='center', zorder=5)

    # --- INVERSION box (top-right): the key design insight ---
    ax.add_patch(FancyBboxPatch((17, 32.5), 30.5, 11.7, boxstyle='round,pad=0.3',
                 fc='#fff5f5', ec='#e53e3e', lw=1.8, zorder=4))
    ax.text(32.2, 38.3, 'ИНВЕРСИЯ ДИЗАЙНА:\nмешает не высота TS, а то, что σ-комплекс\n'
            'диссоциирует РАНЬШЕ OA. Мишень — хелат/полость,\n'
            'удлиняющая τ(σ), + металл с ΔG(OA)≈0\n(микрообратимость = враг).',
            fontsize=9.8, color='#9b2c2c', ha='center', va='center', zorder=5, fontweight='bold')

    # --- QUANTUM-EDGE box (bottom-left): multireference + our measurement ---
    ax.add_patch(FancyBboxPatch((-33, -4.6), 29.5, 11.6, boxstyle='round,pad=0.3',
                 fc='#faf5ff', ec=ACCENT_D, lw=1.6, zorder=4))
    lines = ['КВАНТ-EDGE — наш расчёт (PySCF):']
    if res and 'headline' in res:
        h = res['headline']
        eth = res.get('ethane_curve', [])
        eth_max = max((row.get('nu_66') or 0) for row in eth) if eth else None
        if eth_max is not None:
            lines.append(f'• гомолиз C–C этана: N_u→{eth_max:.2f} (пиролиз = MR)')
        lines.append(f'• OA-TS на Rh-пинцере: N_u≈{h["N_u_TS"]} (слабо)')
    qc = (res or {}).get('quantum_correction')
    if qc and qc.get('dcorr_barrier_kcal') is not None:
        lines.append(f'• поправка к барьеру Δcorr‡≈{qc["dcorr_barrier_kcal"]:+.1f} ккал/моль → мала')
        lines.append('⇒ для барьера OA DFT достаточен; квант — в гомолизе')
    else:
        lines.append('(для RhCl₃ даже CCSD не рекомендован)')
    ax.text(-18.2, 1.2, '\n'.join(lines),
            fontsize=8.2, color=ACCENT_D, ha='center', va='center', zorder=5, linespacing=1.4)

    # --- ECONOMICS box (bottom-right) ---
    ax.add_patch(FancyBboxPatch((17, -4.6), 30.5, 10.2, boxstyle='round,pad=0.3',
                 fc='#fffaf0', ec='#dd6b20', lw=1.6, zorder=4))
    ax.text(32.2, 0.4, 'Экономика честно: катализатор — малая доля\nсебестоимости рециклинга. Рычаг — '
            'селективность /\nмягкая T / чистота мономера vs пиролиз, НЕ\nцена катализатора; '
            'feedstock доминирует.',
            fontsize=9.4, color='#9c4221', ha='center', va='center', zorder=5)

    disclaimer(fig, 'Позиции — литературные тренды (57 кДж/моль Rh-POCOP · 29 °C Brookhart · '
                    'тандем 125–200 °C · Conk 87 %), НЕ расчёт. Та же дескриптор-логика (барьер × ΔG).')
    fig.subplots_adjust(left=0.06, right=0.985, top=0.88, bottom=0.1)
    path = os.path.join(OUT, '08_cc_activation_recycling_descriptor.png')
    fig.savefig(path, dpi=130); plt.close(fig)
    print('wrote', path)

if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    fig_molecules()
    fig_descriptor()
    print('done.')
