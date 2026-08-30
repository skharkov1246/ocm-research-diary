# -*- coding: utf-8 -*-
"""
NaCN 2500 t/y Andrussow modular skid, Tarkwa (Ghana).
NODE: PLUG-AND-PLAY ERECTION & COMMISSIONING.
Part 1 - field joint census and welding / NDT / testing man-hours.

Modules M1..M6 (113.0 t total, heaviest 21.0 t, envelope 11.8 x 2.3 x 2.6 m)
from calc/nacn_logistics/pkg.json (this repo).
Process flows from calc/nacn_logistics/balance.py:
  burner feed 1596 Nm3/h, gas to HCN absorber 1946 m3/h @60C,
  product 1041.7 kg/h of 30% NaCN, NaOH 50% 510 kg/h, NH3 159.7 kg/h.

CODE BASIS (primary, ASME B31.3-2022 Process Piping):
  HCN = "Category M Fluid Service" (B31.3 300.2 definition: a fluid such that a
  single exposure to a very small quantity, caused by leakage, can produce
  serious irreversible harm). Chapter VIII then applies:
   - M302.2.4/M308: threaded and unlisted joints restricted; welded construction
     preferred, flanged joints minimised;
   - M341.4: 100% examination - all girth/miter groove welds radiographed (or UT),
     100% visual, no "random" sampling allowed;
   - M345.8: SENSITIVE LEAK TEST required in ADDITION to the hydrostatic /
     pneumatic leak test; sensitivity not less than 1e-3 atm*mL/s.
  These three lines are what makes an HCN tie-in ~2.5x slower than a utility one.

All labour rates flagged SRC=ASSUMPTION - vendor/contractor quotation required.
No primary productivity source could be retrieved (web budget exhausted).
"""
import math

# ---------------------------------------------------------------
# 1. LINE SIZING for the inter-module tie-ins (from balance.py flows)
# ---------------------------------------------------------------
def d_from_flow(q_m3s, v):
    A = q_m3s / v
    return math.sqrt(4*A/math.pi)*1000  # mm

print("="*72)
print("1. SIZING OF INTER-MODULE PROCESS LINES  (D = sqrt(4*Q/(pi*v)))")
print("="*72)
# reactor off-gas after WHB, ~250 C
q1 = 1596/3600 * (273+250)/273
print(f"M1->M2 reactor off-gas: Q={1596} Nm3/h @250C = {q1:.3f} m3/s, v=15 m/s "
      f"-> D={d_from_flow(q1,15):.0f} mm -> DN250")
q2 = 1946/3600
print(f"M2->M3 gas to HCN abs : Q=1946 m3/h @60C = {q2:.3f} m3/s, v=15 m/s "
      f"-> D={d_from_flow(q2,15):.0f} mm -> DN200")
q3 = 1946/3600*0.88   # HCN + part of water removed
print(f"M3->M6 tail gas (H2)  : Q={q3*3600:.0f} m3/h, v=18 m/s "
      f"-> D={d_from_flow(q3,18):.0f} mm -> DN200")
q4 = 1155.4/3600*(1.013/2.5)   # blower discharge ~1.5 barg
print(f"M4->M1 process air    : Q=1155 Nm3/h @1.5 barg = {q4:.3f} m3/s, v=20 m/s "
      f"-> D={d_from_flow(q4,20):.0f} mm -> DN150")
q5 = 230.5/3600*(1.013/4.0)
print(f"BL->M1 natural gas    : Q=230 Nm3/h @3 barg = {q5:.4f} m3/s, v=20 m/s "
      f"-> D={d_from_flow(q5,20):.0f} mm -> DN50")
q6 = 1041.7/1157/3600  # rho 30% NaCN ~1157 kg/m3
print(f"M3->T2 product 30%    : Q={q6*3600:.2f} m3/h, v=1.5 m/s "
      f"-> D={d_from_flow(q6,1.5):.0f} mm -> DN25 (DN40 chosen for fouling margin)")

# ---------------------------------------------------------------
# 2. FIELD JOINT CENSUS
#    (tag, service, DN, n_butt_welds, n_flanged, catM?)
# ---------------------------------------------------------------
JOINTS = [
 # --- Category M (HCN / CN- / NH3-bearing) --------------------------------
 ("M1-M2 reactor off-gas",        "HCN gas 250C", 250, 2, 2, 1),
 ("M2-M3 gas to HCN absorber",    "HCN gas 60C",  200, 2, 2, 1),
 ("M3-M6 tail gas to oxidiser",   "HCN traces",   200, 2, 1, 1),
 ("M2-M1 NH3 recycle (AP/MP)",    "NH3 aq",        50, 2, 1, 1),
 ("M2 NH3 feed ex hydrolyser",    "NH3 gas",       25, 2, 1, 1),
 ("M1 PSV/rupture-disc header",   "HCN relief",   100, 3, 2, 1),
 ("M2 relief header to M6",       "HCN relief",   100, 3, 2, 1),
 ("M3 relief header to M6",       "HCN relief",   100, 3, 2, 1),
 ("M3-M6 emergency scrubber vent","HCN emerg",    150, 2, 1, 1),
 ("M3-T2 product line",           "NaCN 30%",      40, 4, 4, 1),
 ("T2 loading/return loop",       "NaCN 30%",      40, 4, 4, 1),
 ("M3 analyser/sample tie-ins",   "NaCN/HCN",      15, 6, 6, 1),
 ("M6 scrubber blowdown to T5",   "CN- liquor",    50, 3, 2, 1),
 ("M1/M2/M3 closed drains to T6", "CN- liquor",    50, 6, 3, 1),
 ("M3 NaOH dosing ex T1",         "NaOH 50%",      25, 3, 2, 0),
 ("M6 caustic make-up ex T1",     "NaOH 50%",      25, 3, 2, 0),
 # --- utility / non-Cat-M --------------------------------------------------
 ("M4-all cooling water S+R",     "CW",           150, 8, 8, 0),
 ("M4-all instrument air",        "IA",            50, 4, 6, 0),
 ("N2 header to M1/M2/M3/M6",     "N2",            50, 5, 6, 0),
 ("BL-M1 natural gas",            "NG",            50, 3, 3, 0),
 ("M4-M1 process air",            "air",          150, 2, 2, 0),
 ("M1-M2 BFW / steam / condensate","steam",         50, 6, 4, 0),
 ("Firewater ring + deluge",      "FW",           150, 10, 8, 0),
 ("Service/potable water, drains", "SW",            50, 8, 6, 0),
 ("M6-S5 stack duct (thin-wall, 0.5 mh/DI - see note)","flue",900,3,1,0),
]

print()
print("="*72)
print("2. FIELD JOINT CENSUS (inter-module + battery-limit tie-ins)")
print("="*72)
print(f"{'tie-in':34}{'DN':>5}{'weld':>6}{'flg':>5}{'DI':>7}  catM")
tot = {"di_m":0.0,"di_u":0.0,"w_m":0,"w_u":0,"f_m":0,"f_u":0}
for tag, svc, dn, nw, nf, m in JOINTS:
    di_in = dn/25.4
    di = nw*di_in
    print(f"{tag:34}{dn:5}{nw:6}{nf:5}{di:7.1f}  {'YES' if m else '-'}")
    if m: tot["di_m"]+=di; tot["w_m"]+=nw; tot["f_m"]+=nf
    else: tot["di_u"]+=di; tot["w_u"]+=nw; tot["f_u"]+=nf
print("-"*72)
print(f"Category M : {tot['w_m']:3} butt welds, {tot['di_m']:6.1f} dia-inch, {tot['f_m']:3} flanged joints")
print(f"Utility    : {tot['w_u']:3} butt welds, {tot['di_u']:6.1f} dia-inch, {tot['f_u']:3} flanged joints")
print(f"TOTAL      : {tot['w_m']+tot['w_u']:3} butt welds, "
      f"{tot['di_m']+tot['di_u']:6.1f} dia-inch, {tot['f_m']+tot['f_u']:3} flanged joints")

# ---------------------------------------------------------------
# 3. MAN-HOURS  (SRC=ASSUMPTION rates, band given)
# ---------------------------------------------------------------
R = {
 "mh_per_di_catM": 2.9,  # SS316L GTAW root+fill, Ar purge dams, 100% RT support,
                         # fit-up, alignment, grinding.  band 2.3-3.6
 "mh_per_di_util": 1.3,  # CS SMAW/FCAW, 5-10% RT.       band 1.0-1.8
 "mh_flange_catM": 3.0,  # PCC-1 controlled bolt-up, new gasket, cross-pattern
                         # torque in 4 passes + record.  band 2.2-4.0
 "mh_flange_util": 1.1,
 "mh_rt_per_weld": 1.6,  # RT shot + processing + interpretation + exclusion-zone
                         # set-up (2-man NDT crew).      band 1.2-2.2
 "rt_frac_util":   0.10, # B31.3 341.4.1 normal fluid service: 5% random; 10% used
 "mh_helium_per_joint": 0.35,  # sniffing at 1e-3 atm*mL/s per B31.3 M345.8
}
mh_weld_m = tot["di_m"]*R["mh_per_di_catM"]
mh_weld_u = tot["di_u"]*R["mh_per_di_util"]
mh_flg_m  = tot["f_m"]*R["mh_flange_catM"]
mh_flg_u  = tot["f_u"]*R["mh_flange_util"]
mh_rt_m   = tot["w_m"]*R["mh_rt_per_weld"]
mh_rt_u   = tot["w_u"]*R["rt_frac_util"]*R["mh_rt_per_weld"]
n_catM_joints = tot["w_m"]+tot["f_m"]
mh_he     = n_catM_joints*R["mh_helium_per_joint"]
# repair allowance: reject rate on field SS welds
REJ = 0.06   # 6% RT reject, band 3-12%
mh_rep = (mh_weld_m+mh_rt_m)*REJ*2.2   # cut out, re-prep, re-weld, re-shoot

print()
print("="*72)
print("3. MECHANICAL COMPLETION MAN-HOURS (piping)")
print("="*72)
for lbl,v,expl in [
  ("Cat-M welding", mh_weld_m, f"{tot['di_m']:.1f} DI x {R['mh_per_di_catM']} mh/DI"),
  ("Utility welding", mh_weld_u, f"{tot['di_u']:.1f} DI x {R['mh_per_di_util']} mh/DI"),
  ("Cat-M flange bolt-up", mh_flg_m, f"{tot['f_m']} x {R['mh_flange_catM']} mh"),
  ("Utility flange bolt-up", mh_flg_u, f"{tot['f_u']} x {R['mh_flange_util']} mh"),
  ("RT 100% Cat-M (M341.4)", mh_rt_m, f"{tot['w_m']} welds x {R['mh_rt_per_weld']} mh"),
  ("RT 10% utility", mh_rt_u, f"{tot['w_u']}x{R['rt_frac_util']} x {R['mh_rt_per_weld']} mh"),
  ("Helium sensitive leak test (M345.8)", mh_he, f"{n_catM_joints} joints x {R['mh_helium_per_joint']} mh"),
  ("Weld repair allowance", mh_rep, f"{REJ*100:.0f}% reject x 2.2 rework factor"),
]:
    print(f"  {lbl:38}{v:8.0f} mh   [{expl}]")
mh_pipe = mh_weld_m+mh_weld_u+mh_flg_m+mh_flg_u+mh_rt_m+mh_rt_u+mh_he+mh_rep
print(f"  {'SUBTOTAL interconnect piping':38}{mh_pipe:8.0f} mh")
import json
json.dump({"joints":tot,"mh_pipe":mh_pipe,"n_catM_joints":n_catM_joints,
           "mh_rt_m":mh_rt_m,"mh_he":mh_he},
          open("/home/user/ocm-research-diary/calc/nacn_erection/joints_out.json","w"),indent=1)
