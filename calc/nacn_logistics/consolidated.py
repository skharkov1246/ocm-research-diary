# -*- coding: utf-8 -*-
"""Consolidated: baseline vs optimised, with knocked-down tanks. Final numbers."""
import sys, json, math
sys.path.insert(0,"/home/user/ocm-research-diary/calc/nacn_logistics")
from packages import CONT, RATE
P=[tuple(x) for x in json.load(open("/home/user/ocm-research-diary/calc/nacn_logistics/pkg.json"))]
D={p[0]:p for p in P}
m=lambda t:D[t][2]; v=lambda t:D[t][3]*D[t][4]*D[t][5]; c=lambda t:D[t][6]*1000.0
BOX=6300.; OOGC=23100.; USE_V=CONT["40HC"][3]*0.72; USE_M=24.1; LEV=0.0245; KD=0.18
TANKS=["T1","T2","T3","T4","T5","T6"]
MODS=["M1","M2","M3","M4","M5","M6"]
LOCAL_FINAL=["T4","C1","C2","Z2"]            # fully local
VOIDFILLED_T, VOIDFILLED_V = 25.5, 72.9

def ship_vol(t): return v(t)*KD if t in TANKS else v(t)
ship=[t for t in D if t not in LOCAL_FINAL and t!="G2"]
loose=[t for t in ship if t not in MODS]
lm=sum(m(t) for t in loose)-VOIDFILLED_T + 8.0     # +C2 liner rolls
lv=sum(ship_vol(t) for t in loose)-VOIDFILLED_V + 22.0
n=math.ceil(max(lv/USE_V, lm/USE_M))
print("=== OPTIMISED SHIPMENT (scenario F) ===")
print(f"  6 ISO-frame process modules, each <= 12.192 x 2.438 x 2.896 m, gross 20.9-25.8 t")
print(f"     113.0 t equipment + {VOIDFILLED_T} t void-fill = {113+VOIDFILLED_T:.1f} t inside the frames")
print(f"  loose crated cargo {lm:.1f} t / {lv:.1f} m3 -> {n} x 40'HC "
      f"(volume-driven {lv/USE_V:.2f}, mass-driven {lm/USE_M:.2f})")
print(f"  Pt/Rh gauze 3.2 kg -> air freight")
print(f"  TOTAL UNITS = {6+n}   |  OOG pieces = 0  |  abnormal-load convoys = 0")

fob=sum(c(t) for t in ship)+c("G2")+c("C2")*0.6
fr=6*(BOX+800)+n*BOX+2804
ins=(fob+fr)*RATE["insurance_pct"]
lev_ex=(fob+fr+ins)*LEV; lev_nex=(fob+fr+ins)*0.0745
logF_ex=fr+ins+lev_ex; logF_nex=fr+ins+lev_nex
loc=sum(c(t)*1.30 for t in ["T4"])+c("C1")+c("C2")+c("Z2")
erect_kd=sum(m(t) for t in TANKS if t!="T4")*900
print(f"\n  China FOB shipped      ${fob:>10,.0f}")
print(f"  ocean+port+inland      ${fr:>10,.0f}")
print(f"  marine insurance 1.1%  ${ins:>10,.0f}")
print(f"  levies exempt / not    ${lev_ex:>10,.0f} / ${lev_nex:,.0f}")
print(f"  LOGISTICS TOTAL        ${logF_ex:>10,.0f} (exempt) / ${logF_nex:,.0f} (non-exempt)")
print(f"  Ghana field erection of KD tanks  ${erect_kd:>8,.0f}")

A_log_ex, A_log_nex, A_fob = 435444., 659077., 4146800.
print(f"\n=== VS BASELINE A (9 OOG + 11 box = 20 units) ===")
print(f"  logistics  ${A_log_ex:,.0f} -> ${logF_ex:,.0f}   saving ${A_log_ex-logF_ex:,.0f} "
      f"({100*(1-logF_ex/A_log_ex):.0f}%)   [exempt]")
print(f"  logistics  ${A_log_nex:,.0f} -> ${logF_nex:,.0f}   saving ${A_log_nex-logF_nex:,.0f} "
      f"({100*(1-logF_nex/A_log_nex):.0f}%)   [non-exempt]")
print(f"  units      20 -> {6+n}  ({100*(1-(6+n)/20):.0f}% fewer transport places)")
print(f"  ocean volume {sum(v(t) for t in D if t!='C1'):.0f} m3 -> {sum(ship_vol(t) for t in loose)+6*76.3:.0f} m3")
print(f"  as % of $4.2M CAPEX: logistics {100*A_log_ex/4.2e6:.1f}% -> {100*logF_ex/4.2e6:.1f}% (exempt); "
      f"{100*A_log_nex/4.2e6:.1f}% -> {100*logF_nex/4.2e6:.1f}% (non-exempt)")

print("\n=== LOCAL CONTENT (Ghana value-added), scenario F ===")
tot_val=sum(c(t) for t in D)
gh_material = c("C1")+c("C2")+c("Z2")+c("T4")*1.30
gh_labour   = erect_kd + 100.0*1000   # module hook-up, piping tie-in, E&I termination, commissioning support
# install labour estimate: 6 modules x ~450 mh + tanks + civils supervision
print(f"  Ghana-supplied material+civil   ${gh_material:>10,.0f}")
print(f"  Ghana labour (erection/tie-in)  ${gh_labour:>10,.0f}")
print(f"  Ghana total value-added         ${gh_material+gh_labour:>10,.0f}  "
      f"= {100*(gh_material+gh_labour)/(tot_val+gh_labour+logF_ex):.0f}% of delivered+installed cost")
print(f"  by shipped VOLUME avoided: {100*(1-(sum(ship_vol(t) for t in loose)+6*76.3)/sum(v(t) for t in D if t!='C1')):.0f}% of the "
      f"original ocean volume never leaves China (KD nesting + local T4/C1/C2)")

print("\n=== SENSITIVITY: what if Ghana fabricators quote better? ===")
for k in (1.10,1.20,1.30,1.44):
    extra=[t for t in ["S1","S2","S3","S4","S5","PP3","PP4","E4","T1","T3"]
           if c(t)*k < c(t)+max(ship_vol(t)/USE_V,m(t)/USE_M)*BOX*(1+LEV)]
    saved=sum(max(ship_vol(t)/USE_V,m(t)/USE_M)*BOX*(1+LEV)-c(t)*(k-1) for t in extra)
    print(f"   k={k:4.2f}: additionally localise {extra if extra else '(nothing)'} -> ${saved:,.0f}")
