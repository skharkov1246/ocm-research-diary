# -*- coding: utf-8 -*-
"""Package inventory + localisation decision + freight model."""
import sys, json, math
sys.path.insert(0,"/home/user/ocm-research-diary/calc/nacn_logistics")
from packages import CONT, RATE, GH, CN, DIAG_40HC

# origin: CN = must ship from China ; GH = fabricable/procurable in Ghana ; MIX
# wetted: 1 = pressure boundary or wetted part in HCN/CN- /NH3 service (safety-critical)
# fields: (tag, name, mass_t, L, W, H, fob_kUSD, origin, wetted, note)
P = [
 # ---------- A. process modules (frame-mounted, China shop) ----------
 ("M1","Reaction module: Pt/Rh gauze burner D1.05m, mixer, flame arrestors, WHB, steam drum, start-up burner",
   16.5, 11.8, 2.30, 2.60, 470, "CN", 1, "Pt-gauze burner + HCN-side WHB: pressure boundary, no field joints"),
 ("M2","Quench + NH3 recovery module: quench cooler, NH3 absorber D0.8x10m, MP/AP circuit, pumps",
   19.0, 11.8, 2.30, 2.62, 385, "CN", 1, "column lies flat 10.0 m < 12.03 m internal"),
 ("M3","HCN absorption + NaCN reactor module: HCN absorber D0.8x9m, NaOH dosing, product cooler, pH loop",
   18.5, 11.8, 2.30, 2.62, 405, "CN", 1, "pH>11.5 interlock; all wetted CS/SS certified 3.1"),
 ("M4","Utility module: process air blower 1600 Nm3/h, IA package, chiller, PHEs, CW pumps",
   17.0, 11.8, 2.30, 2.60, 300, "CN", 0, "rotating equipment + PHEs; densest module"),
 ("M5","E-house / control module: MCC, VSDs, PLC+SIS, UPS, HVAC, gas-detection panel",
   21.0, 11.8, 2.30, 2.60, 335, "CN", 0, "certified SIL loops; walk-in, blast-rated wall to process side"),
 ("M6","Abatement module: thermal oxidiser, caustic scrubber, emergency scrubber, ID fan",
   21.0, 11.8, 2.30, 2.62, 265, "CN", 1, "HCN-bearing off-gas; monolithic scrubber shell"),
 # ---------- B. tankage / storage ----------
 ("T1","NaOH 50% storage 120 m3 (14 d), CS, heat-traced, D6.0x4.5m", 14.0, 6.0,6.0,4.6,  62,"GH",0,
   "field-erected CS tank; NaOH service, standard local fab scope"),
 ("T2","NaCN 30% product tank 2x25 m3 (1 d buffer), CS, D3.0x3.8m",  7.2, 3.2,3.2,4.0,  46,"GH",1,
   "CS is standard for alkaline NaCN; local fab ONLY with 3rd-party ITP + full RT + no Zn/Cu/Al contact"),
 ("T3","Urea silo 280 m3 (30 d) + hydrolyser feed hopper, D7.0x8.0m",16.5, 7.2,7.2,8.2,  74,"GH",0,
   "bulk solid, non-hazardous: pure structural scope -> Ghana"),
 ("T4","Firewater + deluge tank 250 m3, CS bolted",                  11.0, 8.0,8.0,5.0,  38,"GH",0,
   "commodity; local"),
 ("T5","Spent-scrubber / cyanide-destruction holding 60 m3, HDPE-lined CS", 5.5,4.5,4.5,4.2,29,"GH",1,
   "wetted CN-; local fab with imported liner + certified welding procedure"),
 ("T6","Sumps, catch pits, drum/IBC decant station (steel)",          4.0, 3.0,2.2,2.2,  18,"GH",1,""),
 # ---------- C. structural steel ----------
 ("S1","Pipe rack + module interconnect steel",                      22.0, 12.0,2.3,2.3, 35,"GH",0,
   "plain sections, galvanised: Tema/Takoradi fabricators"),
 ("S2","Platforms, stairs, ladders, handrail, grating",              18.0, 6.0,2.3,2.3,  29,"GH",0,""),
 ("S3","Equipment support frames, tank saddles, silo legs",          12.5, 6.0,2.3,2.3,  20,"GH",0,""),
 ("S4","Weather canopy / roof over process modules + cladding",       9.5, 6.0,2.3,2.3,  15,"GH",0,
   "Ghana coastal-humid: local design better than shipped sheet"),
 ("S5","Stack 30 m x DN900, 3 sections + guys",                       7.5,10.5,1.1,1.1, 17,"GH",0,
   "rolled plate cylinder, post-abatement duty: local fab OK"),
 # ---------- D. piping bulk ----------
 ("PP1","HCN / NH3 / CN- process spools: SS316L + CS, shop-prefab, 100% RT, PMI",
   9.5, 11.8,2.3,1.2, 245,"CN",1,"safety-critical wetted: prefab spools ONLY, field joints minimised"),
 ("PP2","Valves, safety valves, rupture discs, flame arrestors (process)",
   6.0, 4.0,2.0,1.5, 285,"CN",1,"certified; densest crate in the shipment (~500 kg/m3)"),
 ("PP3","Low-class piping: cooling water, service water, plant air, drains, firewater ring",
   26.0, 12.0,2.3,2.3, 65,"GH",0,"non-hazardous utilities -> Ghana supply+install"),
 ("PP4","Pipe supports, hangers, U-bolts, shoes",                     6.5, 6.0,2.3,2.3, 14,"GH",0,""),
 # ---------- E. electrical & instrumentation ----------
 ("E1","Field instruments, HCN/NH3 gas detectors, analysers, SIS field devices",
   3.0, 3.0,2.0,1.6, 265,"CN",0,"certified Ex/SIL; air-freightable if schedule slips"),
 ("E2","Instrument + fire-resistant + fibre cable, glands, JBs",       6.5, 4.0,2.0,2.0, 78,"CN",0,""),
 ("E3","LV power cable, earthing, lighting, small power",             14.0, 4.0,2.0,2.0, 96,"GH",0,
   "Tema cable makers (Tropical Cable, Nexans Kabelmetal) - IEC 60502, local content"),
 ("E4","Cable ladder/tray, supports, transformers pad, DB boards",     11.0, 6.0,2.3,2.0, 42,"GH",0,""),
 ("E5","Diesel gen-set standby 500 kVA + ATS",                         6.5, 4.0,2.0,2.3, 68,"GH",0,
   "widely stocked in Ghana with local service network - buy local for uptime, not for freight"),
 # ---------- F. civil ----------
 ("C1","Foundations, plinths, slabs, roads (concrete+rebar)",        520.0, 0,0,0,  185,"GH",0,
   "never shipped; 100% local by definition - listed to size the local share honestly"),
 ("C2","Bunding: RC walls + HDPE geomembrane + sump, cyanide-tight",  95.0, 0,0,0,   72,"GH",1,
   "ICMI Code Std 4: liner+sump imported roll goods, civil works local"),
 # ---------- G. safety / misc ----------
 ("G1","Deluge/foam skid, safety showers, escape sets, antidote stations, HCN refuge",
   5.0, 5.0,2.3,2.3, 88,"CN",0,"certified life-safety: import"),
 ("G2","Pt/Rh gauze pack + spare pack (~3.2 kg PGM)",                0.02,0.6,0.6,0.4,155,"CN",0,
   "AIR FREIGHT, high-value, zero sea-logistics impact; lease/buy-back candidate"),
 ("Z1","Commissioning spares + 2-yr ZIP (seals, bearings, gaskets, packing, refractory, elements)",
   12.0, 0,0,0, 128,"CN",0,"VOID-FILL candidate: dense crates, no dedicated container"),
 ("Z2","Erection consumables: welding, bolts, paint, insulation, scaffolding",
   9.0, 0,0,0,  46,"GH",0,"local; scaffolding hired not bought"),
]

def vol(p):
    _,_,m,L,W,H,*_ = p
    return L*W*H

hdr=f"{'tag':4} {'mass_t':>7} {'vol_m3':>7} {'FOB_k$':>7} {'org':>4} {'wet':>3}  name"
print(hdr); print("-"*110)
tot_m=tot_v=tot_c=0
for p in P:
    tag,name,m,L,W,H,c,org,wet,note = p
    v=L*W*H
    tot_m+=m; tot_v+=v; tot_c+=c
    print(f"{tag:4} {m:7.1f} {v:7.1f} {c:7.0f} {org:>4} {wet:>3}  {name[:60]}")
print("-"*110)
print(f"TOTAL {tot_m:7.1f} {tot_v:7.1f} {tot_c:7.0f}")
json.dump([list(p) for p in P], open("/home/user/ocm-research-diary/calc/nacn_logistics/pkg.json","w"),
          ensure_ascii=False, indent=1)
