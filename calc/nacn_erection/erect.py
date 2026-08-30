# -*- coding: utf-8 -*-
"""
Part 2 - civil readiness, crane/setting, E&I, testing, commissioning durations.
Follows joints.py. All rates SRC=ASSUMPTION unless marked.
"""
import math
MH_PIPE = 1002.0   # from joints.py (duct rate corrected)

print("="*74); print("4. CIVIL READINESS (CLIENT SCOPE, BEFORE MODULES ARRIVE)"); print("="*74)
# C1 = 520 t concrete+rebar, C2 = 95 t bunding (pkg.json)
rho_rc, rebar_frac = 2.45, 0.09   # t/m3 reinforced concrete; rebar mass fraction
V_c1 = 520*(1-rebar_frac)/rho_rc
V_c2 = 95*(1-rebar_frac)/rho_rc
print(f"C1 foundations/slabs/roads: 520 t -> V = 520*(1-{rebar_frac})/{rho_rc} = {V_c1:.0f} m3 concrete, "
      f"{520*rebar_frac:.0f} t rebar")
print(f"C2 cyanide-tight bunding  :  95 t -> V = {V_c2:.0f} m3 + HDPE liner")
V_tot = V_c1+V_c2
# module plinths only (the part that gates module setting)
V_plinth = 6*(11.8*2.3*0.9*0.35)  # 6 modules, strip plinths ~35% of footprint x0.9 m deep
print(f"of which MODULE PLINTHS (the schedule-critical part): 6 x 11.8x2.3 strip "
      f"plinths, 0.9 m deep, 35% coverage = {V_plinth:.0f} m3")
# pour rate & cure
pour_rate = 35.0    # m3/day, one mixer truck fleet + pump, small site
d_pour = V_tot/pour_rate
print(f"pour duration = {V_tot:.0f} / {pour_rate} m3/d = {d_pour:.1f} working days (not continuous - "
      f"formwork/rebar cycles dominate)")
print("civil cycle (survey+excavate+blind+rebar+form+pour+strip) at 3 pours/week -> "
      "typ 6-8 weeks for this scope in Ghana with one crew")
# cure: ACI 318-19 (primary) - f'c at 28 d; 70% of f'c typically at 7 d for Type I OPC
print("ACI 318-19: design strength f'c at 28 d. Anchor-bolt/plinth loading normally")
print("  released at >=70% f'c ~ 7 d (Type I OPC) or 3-4 d with 52.5R / early-strength mix.")
print("  => MODULE SETTING GATE = plinth pour + 7 d cure + survey of anchor bolts.")
d_cure = 7

print(); print("="*74); print("5. MODULE SETTING - CRANE SIZING AND LIFT CYCLE"); print("="*74)
m_max, m_spread, m_sling = 21.0, 2.5, 0.6
DAF = 1.25   # dynamic amplification, planning value
W = (m_max+m_spread+m_sling)*DAF
print(f"heaviest module M5 (e-house) = {m_max} t; spreader beam {m_spread} t; rigging {m_sling} t")
print(f"gross load = ({m_max}+{m_spread}+{m_sling}) x DAF {DAF} = {W:.1f} t")
R_lift, util = 14.0, 0.75
print(f"set radius {R_lift} m (crane hardstand outside bund), utilisation limit {util:.0%}")
print(f"required chart capacity = {W:.1f}/{util} = {W/util:.1f} t @ {R_lift} m")
print("  -> 130-160 t class all-terrain (e.g. LTM 1130-5.1 / LTM 1160-5.2) on main boom.")
print("  -> a 100 t class is marginal at 14 m: DO NOT plan on it. SRC=ASSUMPTION (no")
print("     load chart retrieved this session); confirm against actual chart + ground")
print("     bearing pressure (outrigger pads on 40 t/m2 min, Tarkwa laterite needs mats).")
lift_cycle = 0.75+0.5+1.25   # rig / lift+travel+land / shim+survey+derig, hours
n_lifts = 6 + 6 + 3 + 4      # modules, tanks/silo, stack sections, misc steel
print(f"lift cycle = 0.75 rig + 0.50 lift&land + 1.25 shim/survey/derig = {lift_cycle:.1f} h")
print(f"lifts = 6 modules + 6 tanks/silo + 3 stack + 4 misc = {n_lifts}")
h_crane = n_lifts*lift_cycle
print(f"crane hours = {n_lifts} x {lift_cycle:.1f} = {h_crane:.0f} h -> at 9 h/day = {h_crane/9:.1f} crane-days")
print("  + 1 d mobilise/assemble + 1 d demobilise -> book the crane for 5 days.")
# setting tolerance / grout
print("setting tolerance: +-6 mm elevation, +-10 mm plan (else the 44 Cat-M spool")
print("  fit-ups go out of tolerance and become field re-cuts).")
print("grout: non-shrink cementitious, 24-48 h before load transfer, 7 d before hydrotest")
print("  of module-supported piping -> grout cure is ON the critical path, 3 d practical.")

print(); print("="*74); print("6. E&I FIELD SCOPE - THE REAL MODULAR SAVING"); print("="*74)
# I/O estimate for this plant
n_io = 420; n_loops = 185; n_sif = 26
print(f"estimated I/O {n_io}, control loops {n_loops}, safety instrumented functions {n_sif}")
print("STICK-BUILT: every field device home-runs to the marshalling cabinet.")
L_home = 55.0
cable_stick = n_io*L_home
print(f"  {n_io} home runs x {L_home} m avg = {cable_stick:,.0f} m of instrument cable")
print("MODULAR: devices are shop-wired to module JBs; field scope = trunk multicores only.")
n_jb = 6*4 + 4      # 4 JBs per module + 4 field JBs (tanks, bund, gas detection)
n_pwr = 22          # LV feeders from module M5 to motors outside modules + BL
L_trunk = 38.0
cable_mod = (n_jb+n_pwr)*L_trunk
print(f"  ({n_jb} JB trunks + {n_pwr} LV feeders) x {L_trunk} m = {cable_mod:,.0f} m")
print(f"  cable pulled in field: {cable_mod/cable_stick:.0%} of stick-built "
      f"({cable_stick-cable_mod:,.0f} m avoided)")
mh_pull_per_100m = 6.5   # tray-mounted multicore, 2-man
mh_term_per_core = 0.20
cores = (n_jb*24 + n_pwr*4)
mh_ei = cable_mod/100*mh_pull_per_100m + cores*2*mh_term_per_core
print(f"  pulling {cable_mod:.0f} m x {mh_pull_per_100m} mh/100 m = "
      f"{cable_mod/100*mh_pull_per_100m:.0f} mh")
print(f"  terminations {cores} cores x2 ends x {mh_term_per_core} mh = {cores*2*mh_term_per_core:.0f} mh")
mh_loop = n_loops*0.8 + n_sif*6.0
print(f"loop checking: {n_loops} loops x 0.8 mh (shop-FAT'd, field = trunk continuity + "
      f"end-to-end sample) + {n_sif} SIF x 6.0 mh full end-to-end proof test")
print(f"  IEC 61511-1 cl.15 requires FIELD validation of every SIF after installation")
print(f"  (shop FAT does not discharge it) -> {mh_loop:.0f} mh")
mh_ei_tot = mh_ei+mh_loop+180  # +180 earthing, lighting, small power, tray
print(f"E&I FIELD SUBTOTAL = {mh_ei:.0f} + {mh_loop:.0f} + 180 (earthing/lighting/tray) "
      f"= {mh_ei_tot:.0f} mh")

print(); print("="*74); print("7. OTHER FIELD MECHANICAL SCOPE"); print("="*74)
OTH = [("set + align 6 modules, shim, grout, anchor-bolt torque", 6*26),
       ("tanks T1-T6 + silo T3 field erect / set (GH-fabricated)", 380),
       ("stack S5 30 m assemble+erect+guy", 110),
       ("pipe rack S1, platforms S2, frames S3, canopy S4 (86.5 t steel)", 86.5*7.5),
       ("pipe supports PP4, insulation, painting, touch-up", 240),
       ("deluge/foam skid G1, safety showers, refuge, antidote stations", 150),
       ("Pt/Rh gauze pack G2 install + burner internals set (2-man, clean-room disc.)", 24)]
mh_oth = sum(v for _,v in OTH)
for k,v in OTH: print(f"  {k:62}{v:7.0f} mh")
print(f"  {'SUBTOTAL':62}{mh_oth:7.0f} mh")

MH_CONSTR = MH_PIPE + mh_ei_tot + mh_oth
print(); print(f"TOTAL FIELD CONSTRUCTION MAN-HOURS (mech.completion) = "
      f"{MH_PIPE:.0f} + {mh_ei_tot:.0f} + {mh_oth:.0f} = {MH_CONSTR:.0f} mh")

import json; json.dump({"MH_CONSTR":MH_CONSTR,"MH_PIPE":MH_PIPE,"mh_ei_tot":mh_ei_tot,
  "mh_oth":mh_oth,"h_crane":h_crane,"cable_stick":cable_stick,"cable_mod":cable_mod,
  "n_io":n_io,"n_loops":n_loops,"n_sif":n_sif,"V_tot":V_tot,"d_cure":d_cure},
  open("/home/user/ocm-research-diary/calc/nacn_erection/erect_out.json","w"),indent=1)
