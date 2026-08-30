# -*- coding: utf-8 -*-
"""Lever waterfall + break-even k + sensitivity."""
import sys, json, math
sys.path.insert(0,"/home/user/ocm-research-diary/calc/nacn_logistics")
from packages import CONT, RATE
P=[tuple(x) for x in json.load(open("/home/user/ocm-research-diary/calc/nacn_logistics/pkg.json"))]
D={p[0]:p for p in P}
m=lambda t:D[t][2]; v=lambda t:D[t][3]*D[t][4]*D[t][5]; c=lambda t:D[t][6]*1000.0
K={"S1":1.44,"S2":1.44,"S3":1.44,"S4":1.44,"S5":1.44,"T1":1.30,"T2":1.30,"T3":1.25,
   "T4":1.30,"T5":1.35,"T6":1.30,"PP3":1.44,"PP4":1.44,"E3":1.10,"E4":1.48,"E5":1.15}
BOX_VOL=CONT["40HC"][3]; BOX_PAY=CONT["40HC"][4]; STOW=0.72; ROAD=28.0
BOX=RATE["ocean_40HC"]+RATE["dest_port_40"]+RATE["inland_40_tarkwa"]
OOG=RATE["ocean_40FR_OOG"]+RATE["dest_port_OOG"]+RATE["inland_OOG_tarkwa"]
LEV=0.0245

print("A. BREAK-EVEN LOCAL PRICE RATIO k*  (localise if Ghana quote <= k* x China FOB)")
print(f"{'tag':4} {'item':44} {'k*':>6} {'k assumed':>10} {'margin':>8}")
for t in sorted(K, key=lambda x:-(max(v(x)/(BOX_VOL*STOW), m(x)/24.1)*BOX/c(x))):
    fr=max(v(t)/(BOX_VOL*STOW), m(t)/min(BOX_PAY,ROAD-3.9))*BOX
    kstar=1+fr*(1+LEV)/c(t)
    print(f"{t:4} {D[t][1][:44]:44} {kstar:6.2f} {K[t]:10.2f} {kstar-K[t]:+8.2f}")

print("\nB. LEVER WATERFALL  (delta on total logistics + procurement, exempt case)")
lev=[]
# L1: design every module & vessel into ISO gauge -> 9 OOG pieces become in-gauge slots
d1 = 9*OOG - (6*BOX + 3*BOX)
lev.append(("L1  design-to-gauge: 9 OOG pieces -> 9 in-gauge slots", d1))
# L2: modules built ON ISO corner-casting frames (frame = container)
#     saves box hire/tare/stuffing/de-stuffing/extra lift; tare 3.9t x6 freed
d2 = 6*(1250) # box hire+stuffing+destuffing+one extra crane lift, SRC=ASSUMPTION $1250/module
lev.append(("L2  ISO-frame modules (frame is the container)", d2))
# L3: void-fill 25.5 t / 72.9 m3 into module free space -> boxes avoided
d3 = 3*BOX*(1+LEV)
lev.append(("L3  void-fill 25.5 t / 72.9 m3 inside modules  (3 boxes avoided)", d3))
# L4: safety-driven scope deletion (cross-node): 60 t NH3 -> urea
nh3_fob=310000.; nh3_ship=2*OOG+1*BOX
d4 = nh3_ship + nh3_fob*LEV - 0  # silo already in inventory as T3
lev.append(("L4  NH3 bullets (2x25 t, D3.0 m = OOG) deleted by urea route", d4))
# L5: product buffer 7 d -> 1 d : 2x75 m3 tanks (D4.5 OOG) -> 2x25 m3 in-gauge
d5 = 2*OOG - 2*(0.6*BOX)
lev.append(("L5  product buffer 7 d -> 1 d: 2 OOG tanks deleted", d5))
# L6: optimal localisation (T1,T3,T4,T5 only)
d6 = 58338.
lev.append(("L6  localise ONLY the bulky low-value vessels (T1,T3,T4,T5)", d6))
# L7: blanket localisation of the rest (negative!)
d7 = -sum(c(t)*(K[t]-1) - max(v(t)/(BOX_VOL*STOW), m(t)/24.1)*BOX*(1+LEV)
          for t in K if t not in ("T1","T3","T4","T5"))
lev.append(("L7  blanket-localise everything else (DO NOT: net negative)", d7))
tot=0
for n,d in lev:
    if not n.startswith("L7"): tot+=d
    print(f"   {n:62} {d:+10,.0f}")
print(f"   {'NET SAVING (L1..L6)':62} {tot:+10,.0f}")

print("\nC. SENSITIVITY")
for oc in (2500,4000,6000,8000):
    B=oc+RATE["dest_port_40"]+RATE["inland_40_tarkwa"]
    O=oc*3+RATE["dest_port_OOG"]+RATE["inland_OOG_tarkwa"]
    s = 9*O-9*B + 6*1250 + 3*B*(1+LEV) + (2*O+B) + nh3_fob*LEV + (2*O-1.2*B)
    print(f"   ocean 40'HC ${oc:>5,} -> in-gauge/void/scope levers save ${s:>9,.0f}")
for kk in (1.10,1.25,1.44,1.75):
    s=sum(max(v(t)/(BOX_VOL*STOW), m(t)/24.1)*BOX*(1+LEV)-c(t)*(kk-1) for t in K)
    print(f"   uniform local ratio k={kk:4.2f} -> blanket localisation nets ${s:>9,.0f}")
