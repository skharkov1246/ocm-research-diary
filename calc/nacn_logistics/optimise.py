# -*- coding: utf-8 -*-
"""Per-item localisation break-even + lever decomposition (waterfall)."""
import sys, json, math
sys.path.insert(0,"/home/user/ocm-research-diary/calc/nacn_logistics")
from packages import CONT, RATE
P=[tuple(x) for x in json.load(open("/home/user/ocm-research-diary/calc/nacn_logistics/pkg.json"))]
D={p[0]:p for p in P}
m=lambda t:D[t][2]; v=lambda t:D[t][3]*D[t][4]*D[t][5]; c=lambda t:D[t][6]*1000.0
K = {"S1":1.44,"S2":1.44,"S3":1.44,"S4":1.44,"S5":1.44,"T1":1.30,"T2":1.30,"T3":1.25,
     "T4":1.30,"T5":1.35,"T6":1.30,"PP3":1.44,"PP4":1.44,"E3":1.10,"E4":1.48,"E5":1.15}
BOX_VOL=CONT["40HC"][3]; BOX_PAY=CONT["40HC"][4]; STOW=0.72; ROAD=28.0
BOX_COST=RATE["ocean_40HC"]+RATE["dest_port_40"]+RATE["inland_40_tarkwa"]   # 6300
LEV_EX, LEV_NEX = 0.0245, 0.0745

print("PER-ITEM LOCALISATION BREAK-EVEN  (marginal slot cost method)")
print(f"  1 x 40'HC delivered Tarkwa = ${BOX_COST:,.0f}; usable {BOX_VOL*STOW:.1f} m3 or {min(BOX_PAY,ROAD-3.9):.1f} t")
print(f"{'tag':4} {'m,t':>6} {'v,m3':>7} {'slots':>6} {'freight$':>9} {'CN dlvd$':>9} {'GH dlvd$':>9} {'k':>5} {'verdict':>10}")
print("-"*80)
loc_yes=[]; loc_no=[]
for t in sorted(K, key=lambda x:-(v(x))):
    slots = max(v(t)/(BOX_VOL*STOW), m(t)/min(BOX_PAY,ROAD-3.9))
    fr    = slots*BOX_COST
    cn_d  = c(t) + fr + (c(t)+fr)*LEV_EX
    gh_d  = c(t)*K[t]
    good  = gh_d < cn_d
    (loc_yes if good else loc_no).append(t)
    print(f"{t:4} {m(t):6.1f} {v(t):7.1f} {slots:6.2f} {fr:9,.0f} {cn_d:9,.0f} {gh_d:9,.0f} {K[t]:5.2f} "
          f"{'LOCALISE' if good else 'import':>10}  {'+' if good else '-'}${abs(cn_d-gh_d):,.0f}")
print("-"*80)
print("  LOCALISE:", " ".join(loc_yes))
print("  keep imported:", " ".join(loc_no))
sav = sum(max(v(t)/(BOX_VOL*STOW), m(t)/min(BOX_PAY,ROAD-3.9))*BOX_COST*(1+LEV_EX) - (c(t)*K[t]-c(t))
          for t in loc_yes)
print(f"  net saving from optimal localisation only (exempt case): ${sav:,.0f}")
sav_n = sum(max(v(t)/(BOX_VOL*STOW), m(t)/min(BOX_PAY,ROAD-3.9))*BOX_COST*(1+LEV_NEX)
            + c(t)*LEV_NEX - (c(t)*K[t]-c(t)) for t in loc_yes)
print(f"  net saving, non-exempt (localised goods dodge 7.45% levies too): ${sav_n:,.0f}")

# ---- share of scope localised ----
tot_m = sum(m(t) for t in D); tot_v=sum(v(t) for t in D); tot_c=sum(c(t) for t in D)
ship_m= sum(m(t) for t in D if t!="C1"); 
for lbl, S in [("optimal (break-even) set", loc_yes),
               ("all Ghana-capable set", list(K)+["C1","C2","Z2"])]:
    lm=sum(m(t) for t in S); lv=sum(v(t) for t in S); lc=sum(c(t) for t in S)
    print(f"\n  LOCAL SHARE, {lbl}:")
    print(f"    by mass   {lm:7.1f} / {tot_m:7.1f} t   = {100*lm/tot_m:5.1f}%   "
          f"(excl. civil concrete: {100*sum(m(x) for x in S if x!='C1')/(tot_m-m('C1')):.1f}%)")
    print(f"    by volume {lv:7.1f} / {tot_v:7.1f} m3  = {100*lv/tot_v:5.1f}%")
    print(f"    by value  {lc/1000:7.0f} / {tot_c/1000:7.0f} k$  = {100*lc/tot_c:5.1f}%")
