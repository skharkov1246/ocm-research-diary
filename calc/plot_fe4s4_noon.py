#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc/plot_fe4s4_noon.py — honest figure for FeMoco diary Stage 17.

Reads calc/femoco_4fe4s_casscf_results.json (produced by femoco_4fe4s_casscf.py)
and draws the natural-orbital occupations (NOON) of the converged frontier CASSCF
of the [4Fe-4S] cubane. Deviation of the occupations from the closed-shell values
{2, 0} IS the multireference character that broken-symmetry DFT (Stage 15) cannot
represent; the Head-Gordon n_u = sum n_i (2 - n_i) is annotated.

Every value plotted is read from the results JSON — nothing is hand-entered.
Run AFTER the compute: python3 calc/plot_fe4s4_noon.py
"""
import os
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "femoco_4fe4s_casscf_results.json")
OUT = os.path.join(os.path.dirname(HERE), "assets", "femoco", "fe4s4_casscf_noon.png")


def main():
    with open(RESULTS) as fh:
        res = json.load(fh)
    fr = res.get("frontier_casscf", {})
    noon = fr.get("noon")
    if not noon:
        raise SystemExit(f"no frontier NOON in {RESULTS} (status={res.get('status')}, "
                         f"frontier={fr.get('error', fr)})")
    n_u = fr.get("n_u")
    ncas = fr.get("ncas")
    nelec = fr.get("nelec")
    qub = fr.get("qubits_jw")
    conv = fr.get("casscf_converged")
    scf_ok = fr.get("based_on_converged_scf")

    x = list(range(1, len(noon) + 1))
    # colour by distance from the nearest closed-shell integer (0 or 2)
    dev = [min(v, 2.0 - v) for v in noon]
    dmax = max(dev) or 1.0
    colours = [plt.cm.magma(0.25 + 0.55 * (d / dmax)) for d in dev]

    fig, ax = plt.subplots(figsize=(7.2, 4.3), dpi=140)
    ax.axhline(2.0, ls="--", lw=1, color="#888", zorder=0)
    ax.axhline(0.0, ls="--", lw=1, color="#888", zorder=0)
    ax.axhspan(0.15, 1.85, color="#ffd54f", alpha=0.12, zorder=0)  # "open" band
    bars = ax.bar(x, noon, color=colours, edgecolor="#222", linewidth=0.7, zorder=3)
    for xi, v in zip(x, noon):
        ax.text(xi, v + 0.05, f"{v:.2f}", ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x)
    ax.set_xlabel("frontier natural orbital (index)")
    ax.set_ylabel("occupation number")
    ax.set_ylim(0, 2.25)
    ax.set_xlim(0.4, len(noon) + 0.6)

    cas = f"CAS({nelec}e,{ncas}o) = {qub} qubits"
    title = (f"[4Fe–4S] frontier CASSCF NOON — measured multireference character\n"
             f"{cas}   n_u = {n_u}"
             + ("" if conv else "   (CASSCF not fully converged)")
             + ("" if scf_ok else "   (SCF not fully converged)"))
    ax.set_title(title, fontsize=10.5)
    ax.text(0.99, 0.02,
            "occupations away from {0, 2} = static (multireference) correlation\n"
            "— exactly what a single-determinant BS-DFT solution (Stage 15) cannot hold",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
            color="#444")

    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    print(f"wrote {OUT}  (n_u={n_u}, {cas}, casscf_conv={conv}, scf_conv={scf_ok})")


if __name__ == "__main__":
    main()
