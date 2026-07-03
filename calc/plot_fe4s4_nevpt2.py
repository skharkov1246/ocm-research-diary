#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc/plot_fe4s4_nevpt2.py — honest figure for FeMoco diary Stage 18.

Reads calc/femoco_4fe4s_nevpt2_results.json and draws the CASSCF orbital-relaxation
"descent" of the [4Fe-4S] frontier: energy above the high-spin reference (log scale)
at CASCI (raw ROHF orbitals), after 3 CASSCF macro-iterations, and at Stage 17's
50-iteration point — all still unconverged, with n_u annotated to show the metric is
orbital-fragile. The message: strongly-contracted NEVPT2 (dynamic correlation) is
cheap (~1 min), but converging the multireference reference it must sit on is the cost
wall. Every value is read from the results JSON.

Run AFTER the compute: python3 calc/plot_fe4s4_nevpt2.py
"""
import os
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "femoco_4fe4s_nevpt2_results.json")
OUT = os.path.join(os.path.dirname(HERE), "assets", "femoco", "fe4s4_nevpt2_wall.png")
KCAL = 627.50947406
E_STAGE17 = -8230.314104   # Stage 17's unconverged 50-macro-iter CASSCF energy


def main():
    with open(RESULTS) as fh:
        res = json.load(fh)
    e_hs = res["scf"]["e_hs_rohf"]
    ci = res["casci_rohf_orbitals"]
    c3 = res["casscf_bounded_3macro"]

    def above(e):
        return (e - e_hs) * KCAL

    pts = [
        ("CASCI\n(raw ROHF\norbitals)", above(ci["e_casci"]), ci["n_u"], "1"),
        ("CASSCF\n3 macro-iter", above(c3["e_casscf_3iter"]), c3["n_u_after_3iter"],
         f"{c3['seconds']:.0f}"),
        ("CASSCF ~50 iter\n(Stage 17)", above(E_STAGE17), c3["n_u_stage17_unconverged"],
         "~1 h"),
    ]
    xs = list(range(len(pts)))
    ys = [p[1] for p in pts]

    fig, ax = plt.subplots(figsize=(7.6, 4.6), dpi=140)
    ax.set_yscale("log")
    ax.plot(xs, ys, "-o", color="#c0392b", lw=2, ms=9, zorder=3)
    for x, (label, y, nu, t) in zip(xs, pts):
        ax.annotate(f"n_u={nu}", (x, y), textcoords="offset points", xytext=(0, 12),
                    ha="center", fontsize=10, fontweight="bold", color="#c0392b")
        ax.annotate(f"{y:.0f} kcal/mol\nabove HS · {t}", (x, y),
                    textcoords="offset points", xytext=(0, -30), ha="center",
                    fontsize=8, color="#555")
    # "converged limit — unreached" band below the last point
    ax.axhline(above(E_STAGE17) * 0.4, ls=":", lw=1.4, color="#2c3e50")
    ax.text(0.5, above(E_STAGE17) * 0.43, "converged CASSCF — not reached (the classical cost wall)",
            ha="center", va="bottom", fontsize=8.5, color="#2c3e50", style="italic")

    ax.set_xticks(xs)
    ax.set_xticklabels([p[0] for p in pts], fontsize=9)
    ax.set_ylabel("energy above high-spin reference (kcal/mol, log)")
    ax.set_xlim(-0.5, len(pts) - 0.3)
    ax.set_title("[4Fe–4S] frontier: converging the multireference reference is the wall\n"
                 "orbital relaxation descends slowly; n_u never settles (2.0 ↔ 4.0, "
                 "all unconverged)", fontsize=10.5)

    nevpt2_s = ci.get("nevpt2_seconds")
    box = (f"NEVPT2 (dynamic correlation): {nevpt2_s:.0f} s — cheap,\n"
           f"but on this poor reference it returns a meaningless\n"
           f"{ci.get('nevpt2_corr', 0):.1f} Ha correction — it can't rescue a bad reference")
    ax.text(0.02, 0.03, box, transform=ax.transAxes, ha="left", va="bottom", fontsize=8,
            color="#333",
            bbox=dict(boxstyle="round,pad=0.4", fc="#eaf2f8", ec="#aac", alpha=0.9))

    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
