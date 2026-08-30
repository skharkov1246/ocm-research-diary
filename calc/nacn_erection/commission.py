# -*- coding: utf-8 -*-
"""Part 3 - pressure testing, purging, dry-out, hot commissioning: CALENDAR items."""
import math, json
E = json.load(open("/home/user/ocm-research-diary/calc/nacn_erection/erect_out.json"))

print("="*74); print("8. PRESSURE / LEAK TESTING (ASME B31.3 para 345)"); print("="*74)
P_des = {"HCN gas":2.0,"NH3":16.0,"CW":6.0,"IA":10.0,"NG":6.0,"FW":12.0,"product":6.0}
print("345.4.2 hydrostatic: Pt = 1.5*P_des*(St/S); hold >= 10 min, then examine at P_des")
print("345.5.5 pneumatic  : Pt = 1.1*P_des; step to 0.5*Pt, then 1/10 steps, hold 10 min")
print("M345.8 SENSITIVE LEAK TEST (Cat M only): sensitivity >= 1e-3 atm*mL/s,")
print("   BPVC Sec.V Art.10 gas&bubble or He mass-spec. IN ADDITION to 345.4/345.5.")
for k,v in P_des.items():
    print(f"   {k:9} P_des {v:5.1f} barg -> hydro {1.5*v:5.1f} barg / pneum {1.1*v:5.1f} barg")
n_pack = 14
h_pack = 9.0   # blind, fill/pressurise, hold, inspect, drain, dry, reinstate
print(f"test packages: {n_pack} (7 Cat-M + 7 utility); {h_pack} h each incl. reinstatement")
print(f"  -> {n_pack*h_pack:.0f} test-hours; 2 packs/day in parallel -> {n_pack/2:.0f} working days")
print("  Cat-M packs additionally get He mass-spec sniffing of all 77 joints (27 mh)")
print("  and 100% RT film review sign-off BEFORE pressurisation (hold point).")

print(); print("="*74); print("9. DRYING AND INERTING THE HCN CIRCUIT"); print("="*74)
V = {"reactor/burner D1.05x3.0":math.pi/4*1.05**2*3.0,
     "WHB shell + quench":2.0,
     "NH3 absorber D0.8x10":math.pi/4*0.8**2*10,
     "HCN absorber D0.8x9":math.pi/4*0.8**2*9,
     "thermal oxidiser":4.0,"scrubbers 2x":6.0,
     "DN250 x20 m":math.pi/4*0.25**2*20,"DN200 x40 m":math.pi/4*0.20**2*40,
     "small-bore + relief hdr":1.2}
for k,v in V.items(): print(f"   {k:26}{v:6.2f} m3")
Vt = sum(V.values()); Vd = Vt*1.25
print(f"   HCN-wetted volume = {Vt:.1f} m3; +25% dead legs/jackets = {Vd:.1f} m3")
n_id = math.log(20.9/0.5)
print(f"dilution purge (perfect mixing): n = ln(C0/C) = ln(20.9/0.5) = {n_id:.2f} volume changes")
n_real = 7
print(f"real geometry (dead legs, gauze pack, packed columns) -> use n = {n_real}")
Q_n2 = Vd*n_real
print(f"N2 required = {Vd:.1f} x {n_real} = {Q_n2:.0f} Nm3 = {Q_n2*1.25:.0f} kg")
print(f"  = {Q_n2/(2*22.4*1000/28):.2f} standard 50 L/200 bar packs -> ONE liquid-N2 dewar")
print("  (180 kg LN2 = 145 Nm3) covers a purge; keep 2 on site + a 500 L bulk vessel for")
print("  the continuous seal/analyser purge (est. 4-8 Nm3/h -> 100-200 Nm3/day).")
# drying
m_res = 20.0   # kg residual water after drain of hydrotested Cat-M circuit
w_out = 0.005  # kg H2O per kg N2 carried at 60 C sweep, not saturated (conservative)
m_n2 = m_res/w_out; V_n2 = m_n2/1.25
q_dry = 100.0
print(f"dry-out after hydrotest: {m_res} kg residual water; sweep gas picks up "
      f"{w_out} kg H2O/kg N2")
print(f"  N2 = {m_res}/{w_out} = {m_n2:.0f} kg = {V_n2:.0f} Nm3; at {q_dry:.0f} Nm3/h "
      f"-> {V_n2/q_dry:.0f} h = {V_n2/q_dry/24:.1f} days")
print("  WHY IT MATTERS: residual water + HCN + alkali -> HCN polymerisation (black")
print("  polymer) and hydrolysis to formamide/formate. Wet start-up fouls the gauze")
print("  and the absorber packing on day one. Dry-out is NOT optional and NOT fast.")
print("  Mitigation: pneumatic-test the Cat-M gas circuit instead of hydro where B31.3")
print("  345.5 permits -> removes 32 h from the critical path (but 345.5 needs a written")
print("  hazard review; Cat M pneumatic testing also needs owner's approval).")

print(); print("="*74); print("10. REFRACTORY DRY-OUT (thermal oxidiser M6 + WHB)"); print("="*74)
steps = [("ambient->150 C @25 C/h",(150-30)/25),("hold 150 C",12),
         ("150->350 C @25 C/h",(350-150)/25),("hold 350 C",12),
         ("350->550 C @30 C/h",(550-350)/30),("hold 550 C",8),
         ("controlled cool to standby",6)]
tot=0
for k,v in steps: print(f"   {k:26}{v:6.1f} h"); tot+=v
print(f"   TOTAL refractory dry-out = {tot:.0f} h = {tot/24:.1f} days, CONTINUOUS,")
print("   operator-attended, needs fuel gas + ID fan + stack + DCS live.")
print("   SRC=ASSUMPTION: generic castable curve (~1 h hold per 25 mm at each plateau);")
print("   the vendor's actual curve governs. NOT FOUND: a published curve for this unit.")

print(); print("="*74); print("11. HOT COMMISSIONING SEQUENCE (Andrussow specifics)"); print("="*74)
seq = [
 ("N2 leak/tightness at operating P, gas detection live",           8),
 ("air blower + oxidiser on fuel gas only, stack in service",      16),
 ("caustic circulation M3/M6, pH loop tuned to >11.5, deluge test",24),
 ("Pt/Rh gauze pack installed (last item - $155k, theft/damage)",   6),
 ("start-up burner light, gauze preheat to ~600 C on CH4/air",     10),
 ("NH3 cut-in, ratio ramp to Andrussow mix, LIGHT-OFF ~1100 C",     8),
 ("gauze conditioning ('cauliflower' restructuring), yield climbs", 48),
 ("ramp 40%->70%->100% load, 3 plateaus x 12 h",                   36),
 ("stabilise, tune, snag close-out before guarantee run",          24),
]
t=0
for k,v in seq: print(f"   {k:58}{v:4.0f} h"); t+=v
print(f"   subtotal hot commissioning = {t:.0f} h = {t/24:.1f} days elapsed (24/7 shift work)")
print("   NOTE: gauze conditioning (48 h) is why a 72 h guarantee run cannot start on")
print("   the day of light-off - HCN yield on a fresh Pt/Rh gauze is below steady state.")

print(); print("="*74); print("12. GUARANTEE / PERFORMANCE TEST"); print("="*74)
rate_nacn = 312.5   # kg/h 100% NaCN (balance.py)
for h in (72,168,720):
    print(f"   {h:4} h run -> {rate_nacn*h/1000:6.1f} t NaCN(100%) = "
          f"{rate_nacn*h/0.30/1000:6.1f} t of 30% solution")
print("   72 h @ nameplate is the standard mechanical/performance guarantee.")
print("   For a SAFETY-BENCHMARK plant add a 30-day AVAILABILITY test (720 h): it is the")
print("   only test that exposes trip-rate, scrubber saturation and unmanned-night")
print("   behaviour. Cost of that: 225 t NaCN(100%) - which the mine consumes anyway.")
print("   KEY PLUG-AND-PLAY ADVANTAGE: there is NO off-spec disposal problem. Commissioning")
print("   liquor at any strength goes into the leach circuit provided pH>11.5 and free-CN")
print("   is assayed. The customer IS the sink. (Contrast: a merchant plant must dispose.)")

json.dump({"Vd":Vd,"Q_n2":Q_n2,"V_n2_dry":V_n2,"h_dry":V_n2/q_dry,"h_refr":tot,
 "h_hot":t,"n_pack":n_pack},
 open("/home/user/ocm-research-diary/calc/nacn_erection/comm_out.json","w"),indent=1)
