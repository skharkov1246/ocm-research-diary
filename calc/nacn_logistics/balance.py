# -*- coding: utf-8 -*-
"""Mass balance -> equipment size -> shipping envelope. Andrussow, 2500 t/y NaCN."""
MW = {"NaCN":49.007,"HCN":27.026,"NH3":17.031,"CH4":16.043,"NaOH":39.997,"O2":31.999,"urea":60.06}
CAP_T_Y = 2500.0          # t/y NaCN, 100% basis
HOURS   = 8000.0          # h/y  (91.3% availability)
Y_NH3   = 0.68            # Andrussow selectivity on NH3 (industry 0.60-0.72)
Y_CH4   = 0.62            # on CH4 (0.55-0.65)
O2_CH4  = 1.05            # mol O2 per mol CH4 fed

kmol_NaCN_h = CAP_T_Y*1000/MW["NaCN"]/HOURS
hcn_kg_h    = kmol_NaCN_h*MW["HCN"]
nh3_kmol_h  = kmol_NaCN_h/Y_NH3
ch4_kmol_h  = kmol_NaCN_h/Y_CH4
o2_kmol_h   = ch4_kmol_h*O2_CH4
air_kmol_h  = o2_kmol_h/0.2095
tot_kmol_h  = ch4_kmol_h+nh3_kmol_h+air_kmol_h
nm3_h       = tot_kmol_h*22.414
naoh_kg_h   = kmol_NaCN_h*MW["NaOH"]
sol30_kg_h  = CAP_T_Y*1000/HOURS/0.30

print(f"NaCN            {kmol_NaCN_h:8.3f} kmol/h = {CAP_T_Y*1000/HOURS:7.1f} kg/h")
print(f"HCN made        {hcn_kg_h:8.1f} kg/h")
print(f"NH3 feed        {nh3_kmol_h:8.3f} kmol/h = {nh3_kmol_h*MW['NH3']:7.1f} kg/h = {nh3_kmol_h*MW['NH3']*HOURS/1000:6.0f} t/y")
print(f"CH4 feed        {ch4_kmol_h:8.3f} kmol/h = {ch4_kmol_h*22.414:7.1f} Nm3/h = {ch4_kmol_h*22.414*HOURS/1e6:5.2f} MNm3/y")
print(f"Air feed        {air_kmol_h:8.3f} kmol/h = {air_kmol_h*22.414:7.1f} Nm3/h")
print(f"TOTAL burner feed {nm3_h:7.0f} Nm3/h  ({nm3_h/3600:.3f} Nm3/s)")
print(f"NaOH 100%       {naoh_kg_h:8.1f} kg/h -> 50% sol {naoh_kg_h/0.5:7.1f} kg/h = {naoh_kg_h/0.5*HOURS/1000:5.0f} t/y")
print(f"Product 30% sol {sol30_kg_h:8.1f} kg/h = {sol30_kg_h*HOURS/1000:6.0f} t/y")

# --- gauze burner sizing: superficial velocity at gauze ---
# reaction at ~1100 C, ~1.3 bara -> volumetric expansion
T=1373.0; P=1.3
m3_s_hot = tot_kmol_h*22.414*(T/273.15)*(1.013/P)/3600
for v in (1.5, 2.0, 2.5):
    A = m3_s_hot/v
    d = (4*A/3.14159)**0.5
    print(f"  gauze: v_hot={v} m/s -> A={A:.4f} m2, D={d*1000:.0f} mm")

# --- HCN absorber sizing (packed, NaOH circulation) ---
# design gas velocity 1.2 m/s at ~60 C exit of WHB/quench
T2=333.0
g_m3_s = tot_kmol_h*22.414*(T2/273.15)/3600
for v in (1.0,1.2,1.5):
    A=g_m3_s/v; d=(4*A/3.14159)**0.5
    print(f"  absorber: v={v} m/s -> D={d*1000:.0f} mm")
print(f"  gas to absorber {g_m3_s*3600:.0f} m3/h @60C")

# --- storage volumes ---
rho30 = 1160.0  # kg/m3, NaCN 30% solution ~1.16 (alkaline)
for days in (1,3,7):
    v = sol30_kg_h*24*days/rho30
    print(f"  product buffer {days} d -> {sol30_kg_h*24*days/1000:6.1f} t sol -> {v:6.1f} m3")
# NaOH 50% storage, 14 d
print(f"  NaOH 50% 14 d -> {naoh_kg_h/0.5*24*14/1000:5.1f} t -> {naoh_kg_h/0.5*24*14/1525:5.1f} m3 (rho 1525)")
# urea route: urea -> NH3 (CO(NH2)2 + H2O -> 2NH3 + CO2), 30 d silo
urea_kg_h = nh3_kmol_h/2*MW["urea"]
print(f"  urea route: {urea_kg_h:.1f} kg/h = {urea_kg_h*HOURS/1000:.0f} t/y; 30 d silo = {urea_kg_h*24*30/1000:.1f} t -> {urea_kg_h*24*30/750:.1f} m3 (bulk 750 kg/m3)")
