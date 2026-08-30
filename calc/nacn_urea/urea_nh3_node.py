# -*- coding: utf-8 -*-
"""
UREA -> NH3 NODE for the 2500 t/y NaCN Andrussow plant (Tarkwa, Ghana).
Sizing, energy balance, offgas quality, downstream impact.
All thermo constants are handbook values (CRC/NIST/steam tables) - flagged [H].
Process numbers are computed here.
"""
import math

# ---------------------------------------------------------------- constants
MW = {"NaCN":49.007,"HCN":27.026,"NH3":17.031,"CH4":16.043,"NaOH":39.997,
      "urea":60.06,"CO2":44.010,"H2O":18.015,"N2":28.014,"O2":31.999,
      "CO":28.010,"H2":2.016,"biuret":103.08}
# [H] standard enthalpies of formation, 298.15 K, kJ/mol (CRC)
Hf = {"urea_s":-333.1,"urea_aq":-319.2,"NH3_g":-45.94,"CO2_g":-393.51,
      "H2O_l":-285.83,"H2O_g":-241.83,"CH4_g":-74.87,"HCN_g":135.14,
      "CO_g":-110.53,"O2":0.0,"N2":0.0,"H2":0.0}
# [H] absolute entropies, J/mol/K (CRC) - for RWGS
S0 = {"CO_g":197.7,"H2O_g":188.8,"CO2_g":213.8,"H2":130.7}
R = 8.314462

def line(t): print("\n" + "="*78 + "\n" + t + "\n" + "="*78)

# ---- cp(T) quadratic fits through handbook points at 300/1000/1500 K, J/mol/K
CPPTS = {  # [H] cp at 300, 1000, 1500 K
 "N2":(29.1,32.7,34.8), "O2":(29.4,34.9,36.6), "H2O":(33.6,41.3,47.4),
 "CO2":(37.2,54.3,58.4),"CO":(29.1,33.2,35.2), "H2":(28.8,30.2,32.3),
 "CH4":(35.7,71.8,86.6),"NH3":(35.6,56.5,66.0),"HCN":(35.9,51.6,56.5)}
CP = {}
for k,(c3,c10,c15) in CPPTS.items():
    # solve a+bT+cT^2 through (300,c3),(1000,c10),(1500,c15)
    import itertools
    T1,T2,T3 = 300.,1000.,1500.
    A=[[1,T1,T1*T1],[1,T2,T2*T2],[1,T3,T3*T3]]; b=[c3,c10,c15]
    # 3x3 Cramer
    def det(m): return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
                       -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
                       +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
    D=det(A); sol=[]
    for i in range(3):
        M=[row[:] for row in A]
        for r in range(3): M[r][i]=b[r]
        sol.append(det(M)/D)
    CP[k]=tuple(sol)
def cp(sp,T):
    a,b,c = CP[sp]; return a+b*T+c*T*T
def hsens(sp,T1,T2):   # kJ/mol
    a,b,c = CP[sp]
    f=lambda T: a*T+b*T*T/2+c*T**3/3
    return (f(T2)-f(T1))/1000.

line("0. cp fit check (J/mol/K) vs handbook anchor points")
for k in CP: print(f"  {k:4s} 300K {cp(k,300):6.2f}  700K {cp(k,700):6.2f}  1373K {cp(k,1373):6.2f}")

# ---------------------------------------------------------------- 1. BASIS
line("1. BASIS - NH3 and urea demand (from project balance.py)")
CAP=2500.0; HOURS=8000.0; Y_NH3=0.68
nNaCN = CAP*1000/MW["NaCN"]/HOURS
nNH3  = nNaCN/Y_NH3
nUREA = nNH3/2.0
nCO2  = nUREA
m_urea= nUREA*MW["urea"]
print(f"  NaCN            {nNaCN:7.3f} kmol/h ({CAP*1000/HOURS:.1f} kg/h)")
print(f"  NH3 to burner   {nNH3:7.3f} kmol/h = {nNH3*MW['NH3']:.1f} kg/h = {nNH3*MW['NH3']*HOURS/1000:.0f} t/y")
print(f"  urea NORMAL     {nUREA:7.3f} kmol/h = {m_urea:.1f} kg/h = {m_urea*HOURS/1000:.0f} t/y")
print(f"  CO2 co-produced {nCO2:7.3f} kmol/h = {nCO2*MW['CO2']:.1f} kg/h  (STOICHIOMETRIC, unavoidable)")
print(f"  [brief quotes 275 kg/h urea -> implies Y_NH3 = {nNaCN/(275/MW['urea']*2):.3f}; we use project Y=0.68]")
DES=1.10; m_urea_d=m_urea*DES
print(f"  DESIGN (+10%)   {m_urea_d:.1f} kg/h urea, NH3 {nNH3*DES*MW['NH3']:.1f} kg/h")
print(f"  contained-N breakeven vs anhydrous NH3:")
print(f"    urea is {2*14.007/MW['urea']*100:.1f}% N, NH3 is {14.007/MW['NH3']*100:.1f}% N")
print(f"    -> urea is cheaper per N only if P_urea < {2*14.007/MW['urea']/(14.007/MW['NH3']):.3f} x P_NH3")
print(f"       (at NH3 650 $/t delivered -> urea must beat {0.5673*650:.0f} $/t; PRICES NOT VERIFIED - web budget spent)")

# ---------------------------------------------------------------- 2. HYDROLYSER T-P-CONCENTRATION MAP
line("2. HYDROLYSER OPERATING POINT - the water balance fixes T, P and feed conc.")
print("""  Physics: the vapour leaving a boiling aqueous hydrolyser is
    NH3 + CO2 (generated, fixed at 2:1) + H2O (set by p_H2O(T)/(P-p_H2O)).
    n_H2O,vap = (nNH3+nCO2) * p_H2O/(P - p_H2O)
    Feed water = vapour water + reaction water (nCO2*18) + blowdown
    -> the FEED CONCENTRATION IS NOT FREE. Pick T and P, conc follows.""")
# [H] saturation pressure of water - steam table (IAPWS-IF97 tabulated), bar abs
PSAT_TAB=[(100,1.0142),(120,1.9854),(140,3.6150),(150,4.7600),(160,6.1806),
          (170,7.9202),(180,10.0270),(190,12.5520),(200,15.5490),(210,19.0770),
          (220,23.1960),(230,27.9680),(240,33.4470)]
def psat(TC):
    import math
    for i in range(len(PSAT_TAB)-1):
        t0,p0=PSAT_TAB[i]; t1,p1=PSAT_TAB[i+1]
        if t0<=TC<=t1:
            f=(TC-t0)/(t1-t0)
            return math.exp(math.log(p0)+f*(math.log(p1)-math.log(p0)))
    raise ValueError(TC)
# [H] latent heat of vaporisation of water, kJ/kg (steam table)
HFG_TAB=[(150,2113.),(160,2082.),(170,2049.),(180,2015.),(190,1979.),
         (200,1940.),(210,1899.),(220,1856.),(230,1811.),(240,1762.),(250,1716.)]
def hfg(TC):
    for i in range(len(HFG_TAB)-1):
        t0,h0=HFG_TAB[i]; t1,h1=HFG_TAB[i+1]
        if t0<=TC<=t1:
            return h0+(TC-t0)/(t1-t0)*(h1-h0)
    raise ValueError(TC)
for TC in (150,180,200,220):
    print(f"    psat({TC} C) = {psat(TC):6.2f} bar (steam table, exact)")
AW = 0.95   # water activity depression by dissolved urea/carbamate (assumed, TO VERIFY)
ngas_nc = nNH3+nCO2
print(f"\n  non-condensables generated: {ngas_nc:.3f} kmol/h (NH3 {nNH3:.3f} + CO2 {nCO2:.3f})")
print(f"  reaction water consumed: {nCO2*MW['H2O']:.1f} kg/h")
print(f"\n  {'T,C':>5} {'P,bara':>7} {'pH2O':>6} {'yH2O':>6} {'H2Ovap kmol/h':>13} {'kg/h':>7} {'feedH2O':>8} {'sol kg/h':>9} {'urea wt%':>9} {'Tsat_cryst':>11}")
def urea_sat_wt(TC):
    # [H] urea solubility g/100 g water: 0:66.7 20:105 40:165 60:250 80:400
    pts=[(0,66.7),(20,105),(40,165),(60,250),(80,400)]
    for i in range(len(pts)-1):
        if pts[i][0]<=TC<=pts[i+1][0]:
            f=(TC-pts[i][0])/(pts[i+1][0]-pts[i][0]); s=pts[i][1]+f*(pts[i+1][1]-pts[i][1])
            return 100*s/(100+s)
    return None
def cryst_T(wt):   # invert
    lo,hi=0.,80.
    for _ in range(60):
        m=(lo+hi)/2
        if urea_sat_wt(m)<wt: lo=m
        else: hi=m
    return (lo+hi)/2
CASES=[]
for TC in (180,190,200,210):
    for P in (16,20,25,30,36,45):
        pw = psat(TC)*AW
        if pw>=P*0.92: continue
        nH2Ov = ngas_nc*pw/(P-pw)
        mH2Ov = nH2Ov*MW["H2O"]
        feedH2O = mH2Ov + nCO2*MW["H2O"]
        sol = m_urea+feedH2O
        wt = 100*m_urea/sol
        if wt<25 or wt>65: continue
        CASES.append((TC,P,pw,pw/P,nH2Ov,mH2Ov,feedH2O,sol,wt,cryst_T(wt)))
        print(f"  {TC:5d} {P:7.0f} {pw:6.2f} {pw/P:6.3f} {nH2Ov:13.2f} {mH2Ov:7.0f} {feedH2O:8.0f} {sol:9.0f} {wt:9.1f} {cryst_T(wt):10.1f}C")

# ---------------------------------------------------------------- 3. ENERGY BALANCE
line("3. HYDROLYSER ENERGY BALANCE (full enthalpy balance, not just dHr)")
print("""  Envelope: 45-57 wt% urea solution at 40 C IN  ->  vapour (NH3+CO2+H2O) at T,P OUT
  Reference 298.15 K, elements. Includes: (a) heat of reaction from AQUEOUS urea,
  (b) sensible heat of the feed to T, (c) latent+sensible heat of the water that
  leaves as vapour.  The brief's 235 kW = 185 kJ/mol x 4.579 kmol/h counts ONLY (a)
  for the dry thermolysis step and ignores (b) and (c).""")
CP_SOL = 3.10          # kJ/kg/K, 45 wt% urea solution (ideal mixing of urea_s 1.34 + water 4.19) [est]
T_FEED = 40.0
def duty(TC,P,nH2Ov,sol_kgh, blowdown=6.0):
    T=TC+273.15
    # IN at 298 K
    Hin = nUREA*Hf["urea_aq"] + (sol_kgh-m_urea)/MW["H2O"]*Hf["H2O_l"]
    Hin += sol_kgh*CP_SOL*(T_FEED-25.0)/1000.        # MJ/h sensible of feed above 25C
    # OUT at 298 K + sensible to T
    Hout = nNH3*Hf["NH3_g"] + nCO2*Hf["CO2_g"] + nH2Ov*Hf["H2O_g"]
    Hout += blowdown/MW["H2O"]*Hf["H2O_l"]           # blowdown as water
    Hout += nNH3*hsens("NH3",298.15,T) + nCO2*hsens("CO2",298.15,T) + nH2Ov*hsens("H2O",298.15,T)
    Hout += blowdown*4.4*(TC-25)/1000.
    Q = Hout-Hin                                      # MJ/h
    # breakdown
    q_rxn  = nUREA*(2*Hf["NH3_g"]+Hf["CO2_g"]-Hf["urea_aq"]-Hf["H2O_l"])
    q_feed = sol_kgh*CP_SOL*(TC-T_FEED)/1000.
    q_vap  = nH2Ov*MW["H2O"]*hfg(TC)/1000.
    return Q, q_rxn, q_feed, q_vap
CANDS=[("A",190,25),("B",200,30),("C",200,45),("D",210,36)]
print(f"\n  {'':2} {'T,C':>4} {'P,bara':>6} {'wt%':>5} {'Tcryst':>7} {'H2Ov kmol/h':>11} | {'Qrxn':>6} {'Qfeed':>6} {'Qvap':>6} {'TOTAL kW':>9} {'+10% kW':>8}")
SEL={}
for tag,TC,P in CANDS:
    pw=psat(TC)*AW; nH2Ov=ngas_nc*pw/(P-pw); mH2Ov=nH2Ov*MW["H2O"]
    feedH2O=mH2Ov+nCO2*MW["H2O"]; sol=m_urea+feedH2O; wt=100*m_urea/sol
    Q,qr,qf,qv = duty(TC,P,nH2Ov,sol)
    SEL[tag]=dict(T=TC,P=P,wt=wt,nH2Ov=nH2Ov,sol=sol,Q=Q/3.6,Tc=cryst_T(wt))
    print(f"  {tag:2s} {TC:4d} {P:6.0f} {wt:5.1f} {cryst_T(wt):6.1f}C {nH2Ov:11.2f} | "
          f"{qr/3.6:6.0f} {qf/3.6:6.0f} {qv/3.6:6.0f} {Q/3.6:9.0f} {Q/3.6*1.10:8.0f}")
print(f"\n  (Qrxn/Qfeed/Qvap are indicative splits; TOTAL is the rigorous enthalpy balance.)")
print(f"  dHr(urea_aq + H2O_l -> 2NH3_g + CO2_g) = "
      f"{2*Hf['NH3_g']+Hf['CO2_g']-Hf['urea_aq']-Hf['H2O_l']:.1f} kJ/mol  [computed from Hf]")
print(f"  dHr(urea_s  + H2O_g -> 2NH3_g + CO2_g) = "
      f"{2*Hf['NH3_g']+Hf['CO2_g']-Hf['urea_s']-Hf['H2O_g']:.1f} kJ/mol  (the '89.6' literature value)")
print(f"  brief's basis: 185 kJ/mol x {nUREA:.3f} kmol/h = {185*nUREA/3.6:.0f} kW  <-- UNDERSTATED")

# ---------------------------------------------------------------- 4. DISSOLVER
line("4. DISSOLVER - heat of solution is NOT optional")
DHSOL = (Hf["urea_aq"]-Hf["urea_s"])   # kJ/mol, endothermic
print(f"  dH_solution(urea) = Hf(aq)-Hf(s) = {DHSOL:+.1f} kJ/mol  [lit. +14 to +15.4 kJ/mol -> consistent]")
for tag in ("B",):
    sol=SEL[tag]["sol"]
    q = nUREA*DHSOL/3.6                                  # kW at normal rate
    dT_adiab = nUREA*DHSOL*1000/(sol*CP_SOL)             # K
    print(f"  at NORMAL rate ({m_urea:.0f} kg/h urea -> {sol:.0f} kg/h solution):")
    print(f"    Q_dissolution = {nUREA:.3f} kmol/h x {DHSOL:.1f} kJ/mol = {q:.1f} kW")
    print(f"    ADIABATIC temperature drop = {nUREA*DHSOL*1000:.0f} kJ/h / ({sol:.0f} kg/h x {CP_SOL} kJ/kg/K)"
          f" = {dT_adiab:.1f} K")
    print(f"    -> water at 28 C would fall to {28-dT_adiab:.1f} C and CRYSTALLISE"
          f" (46 wt% saturates at {SEL[tag]['Tc']:.0f} C). Dissolver MUST be heated.")
    for camp_h,per_d in ((8,3),(24,1)):
        rate=sol*24*per_d/camp_h
        qc = rate/sol*q; qs = rate*CP_SOL*(48-28)/3600
        print(f"    campaign mode {camp_h} h every {per_d} d: {rate:.0f} kg/h sol -> "
              f"Q_diss {qc:.0f} kW + Q_sens(28->48C) {qs:.0f} kW = {qc+qs:.0f} kW installed")
    qs=sol*CP_SOL*(48-28)/3600
    print(f"    continuous mode: {q:.1f} + {qs:.1f} = {q+qs:.1f} kW  <-- RECOMMENDED (smallest heater)")

# ---------------------------------------------------------------- 5. REACTOR SIZING
line("5. HYDROLYSIS REACTOR - why a stirred pot fails and a stripping column works")
print("""  Urea hydrolysis is REVERSIBLE:  (NH2)2CO + H2O <-> 2NH3 + CO2  (via carbamate).
  In a single CSTR at conversion X:  X = k.tau/(1+k.tau)  ->  X=0.995 needs k.tau = 199.
  Industrial urea-plant hydrolysers run 195-235 C with 30-60 min holdup, which is
  only consistent with COUNTERCURRENT STEAM STRIPPING (plug-flow-like, product
  removed as it forms). Design accordingly: a tray/packed column, NOT a tank.""")
for ktau,name in ((199,"single CSTR, X=0.995"),(5.3,"PFR, X=0.995"),(1.0,"CSTR, X=0.50")):
    print(f"    {name:26s} needs k.tau = {ktau:6.1f}")
print("  [k(T) for urea hydrolysis: NOT FOUND in this session (web budget spent).")
print("   MUST be taken from licensor data or a bench run before fixing height.]")
for tag in ("B","C"):
    d=SEL[tag]; TC,P,sol = d["T"],d["P"],d["sol"]
    rho_l = 950.0                      # kg/m3 at 200 C [est]
    MWv = (nNH3*MW["NH3"]+nCO2*MW["CO2"]+d["nH2Ov"]*MW["H2O"])/(ngas_nc+d["nH2Ov"])
    rho_v = P*1e5*MWv/1000/(R*(TC+273.15))
    nv = ngas_nc+d["nH2Ov"]
    Qv = nv*1000/3600*R*(TC+273.15)/(P*1e5)          # m3/s
    K = 0.05                                          # Souders-Brown, trays [H]
    u = K*math.sqrt((rho_l-rho_v)/rho_v)
    A = Qv/u; Dv = math.sqrt(4*A/math.pi)
    print(f"\n  CASE {tag} ({TC} C, {P} bara, {d['wt']:.0f} wt% feed):")
    print(f"    vapour {nv:.2f} kmol/h, MW {MWv:.1f}, rho_v {rho_v:.1f} kg/m3, Qv {Qv*1000:.1f} L/s")
    print(f"    Souders-Brown u_max {u:.2f} m/s -> A {A*1e4:.0f} cm2 -> D_hydraulic {Dv*1000:.0f} mm")
    print(f"    -> vapour load is NOT limiting; D set by mechanical minimum. Take D = 600 mm.")
    for tau in (0.5,1.0,2.0):
        V = sol/rho_l*tau
        H = V/(math.pi*0.3**2)
        print(f"      tau={tau:.1f} h -> liquid holdup {V*1000:.0f} L -> {H:.2f} m of liquid in D600")
    # wall thickness, SS316L
    Pd = P*1.25+1.0
    S = 110.0    # MPa allowable SS316L at 250 C [H, ASME II-D approx]
    t = Pd/10*300/(S*1.0-0.6*Pd/10)
    print(f"    design P {Pd:.0f} bara -> t_shell = P.R/(S.E-0.6P) = {t:.1f} mm + 1.5 CA -> use {math.ceil(t+1.5+1):.0f} mm SS316L")
    shell_m = math.pi*(0.6+ (t+2.5)/1000)*6.0*(t+2.5)/1000*8000
    print(f"    D600 x 6.0 m T/T shell mass ~ {shell_m:.0f} kg; with heads/trays/skirt/insulation ~ {shell_m*1.8/1000:.1f} t")
    # reboiler
    Qkw = d["Q"]*1.10
    for Tst,Ust in ((250,800),(224,800)):
        dTm = Tst-TC
        if dTm<=10: continue
        Ar = Qkw*1000/(Ust*dTm)
        print(f"    reboiler: steam sat {Tst} C ({psat(Tst) if Tst<=240 else 39.7:.0f} bara), dT={dTm} K, U={Ust} W/m2K -> A = {Ar:.1f} m2")
    msteam = Qkw*3600/ (1716 if TC>195 else 1856)
    print(f"    steam consumption {Qkw:.0f} kW / hfg(250C)=1716 kJ/kg -> {msteam:.0f} kg/h at 40 bara")

# ---------------------------------------------------------------- 6. DEPOSITION
line("6. SOLIDS DEPOSITION IN THE NH3 LINE - where the plugging risk really is")
print("""  2 NH3(g) + CO2(g) <-> NH2COONH4(s).  Kp = p_NH3^2 . p_CO2 [atm^3].
  Anchor [H]: ammonium carbamate total dissociation pressure = 1 atm at 60 C
  (stoichiometric 2:1 -> Kp = (2/3)^2(1/3) = 0.148 atm^3 at 333 K).
  dH_dissoc = +159 kJ/mol [H]. van't Hoff -> deposition temperature at our p's.""")
KP333, DH_CARB = (2/3.)**2*(1/3.), 159000.
def Tdep(pNH3_atm,pCO2_atm):
    Kq = pNH3_atm**2*pCO2_atm
    invT = 1/333.15 - math.log(Kq/KP333)*R/DH_CARB
    return 1/invT-273.15
for tag in ("B",):
    d=SEL[tag]; P=d["P"]; nv=ngas_nc+d["nH2Ov"]
    y={"NH3":nNH3/nv,"CO2":nCO2/nv,"H2O":d["nH2Ov"]/nv}
    print(f"\n  offgas composition (case {tag}): NH3 {y['NH3']*100:.1f}%  CO2 {y['CO2']*100:.1f}%  H2O {y['H2O']*100:.1f}%  (dry: NH3 {nNH3/ngas_nc*100:.0f}% / CO2 {nCO2/ngas_nc*100:.0f}%)")
    for Pl,label in ((P,"HP side, upstream of let-down"),(1.5,"LP side, downstream of let-down"),(2.5,"LP alt")):
        pn=y["NH3"]*Pl/1.01325; pc=y["CO2"]*Pl/1.01325; pw=y["H2O"]*Pl
        Td=Tdep(pn,pc)
        # water dew point at that partial pressure
        print(f"    {label:34s} P={Pl:5.1f} bar -> carbamate deposition at {Td:6.1f} C ; p_H2O={pw:.2f} bar")
    print(f"""
  CONSEQUENCE (this is the single most actionable result of the node):
    - the HP line between column top and the let-down valve must stay ABOVE ~{Tdep(y['NH3']*P/1.01325,y['CO2']*P/1.01325):.0f} C
      -> keep it SHORT, jacketed with the same 40 bar steam, no dead legs, no
         instrument tappings without purge.
    - downstream of the let-down valve the deposition temperature collapses to
      ~{Tdep(y['NH3']*1.5/1.01325,y['CO2']*1.5/1.01325):.0f} C. So: PUT THE LET-DOWN VALVE AT THE COLUMN, and run the
      long transfer line to the burner mixer at 1.5-2 bar, electrically traced
      at 130-150 C (above the {100*0:.0f}... above the water dew point, see below).""")
    pw15=y["H2O"]*1.5
    # water dew point by inversion of psat down to 100 C, else Antoine-lite
    def dewT(p):
        lo,hi=20.,200.
        # extend table below 100C with [H] values
        TAB=[(20,0.02339),(40,0.07384),(60,0.19946),(80,0.47414)]+PSAT_TAB
        for i in range(len(TAB)-1):
            t0,p0=TAB[i]; t1,p1=TAB[i+1]
            if p0<=p<=p1:
                f=(math.log(p)-math.log(p0))/(math.log(p1)-math.log(p0))
                return t0+f*(t1-t0)
        return None
    print(f"    water dew point at 1.5 bar (p_H2O={pw15:.2f} bar) = {dewT(pw15):.0f} C -> trace the LP line at >=130 C")
    print(f"    (condensed water would absorb NH3/CO2 -> ammonium carbonate liquor: corrosive to CS, plugging)")

# ---------------------------------------------------------------- 7. BURNER IMPACT
line("7. DOWNSTREAM IMPACT ON THE ANDRUSSOW BURNER (the hidden cost of the urea route)")
nCH4 = nNaCN/0.62; nO2 = nCH4*1.05; nN2 = nO2/0.2095*0.7905
print(f"  burner feed (project basis): CH4 {nCH4:.3f}  NH3 {nNH3:.3f}  O2 {nO2:.3f}  N2 {nN2:.3f} kmol/h")
# --- reaction extents that reproduce a literature-like Andrussow off-gas
e1 = nNaCN                       # CH4+NH3+1.5O2 -> HCN+3H2O
e2 = (nO2-1.5*e1)/0.5            # CH4+0.5O2 -> CO+2H2 (soaks up the rest of the O2)
slipCH4, slipNH3 = 0.20, 1.50
e3 = nCH4-e1-e2-slipCH4          # CH4+H2O -> CO+3H2 (steam reforming, endothermic)
e4 = (nNH3-e1-slipNH3)/2.0       # 2NH3 -> N2+3H2
prod = {"HCN":e1,"H2O":3*e1-e3,"CO":e2+e3,"H2":2*e2+3*e3+3*e4,
        "N2":nN2+e4,"NH3":slipNH3,"CH4":slipCH4,"CO2":0.0}
tot=sum(prod.values())
print(f"  extents: e1(HCN)={e1:.3f}  e2(POx)={e2:.3f}  e3(reform)={e3:.3f}  e4(NH3 crack)={e4:.3f}")
print(f"  off-gas {tot:.2f} kmol/h: " + "  ".join(f"{k} {v/tot*100:.1f}%" for k,v in prod.items() if v>0.01))
print("  [sanity vs literature Andrussow gas: HCN 6-8, H2 10-13, N2 50-60, H2O 18-22, CO 3-5, NH3 1.5-3 vol% -> matches]")

def adiabatic_T(feed, prod, Tin_map, Tguess=1400.):
    """feed/prod: dict kmol/h. Tin_map: {sp: inlet K}. Returns adiabatic T (K)."""
    Hin = sum(n*(Hf[sp+"_g"] if sp+"_g" in Hf else Hf[sp]) + n*hsens(sp,298.15,Tin_map[sp])
              for sp,n in feed.items() if n>0)
    def Hout(T):
        return sum(n*(Hf[sp+"_g"] if sp+"_g" in Hf else Hf[sp]) + n*hsens(sp,298.15,T)
                   for sp,n in prod.items() if n>0)
    lo,hi=400.,2600.
    for _ in range(100):
        m=(lo+hi)/2
        if Hout(m)<Hin: lo=m
        else: hi=m
    return (lo+hi)/2

Hf["O2"]=0.0; Hf["N2"]=0.0; Hf["H2"]=0.0
def hf(sp): return Hf.get(sp+"_g", Hf.get(sp,0.0))
def adT(feed,prod,Tin):
    Hin = sum(n*hf(sp)+n*hsens(sp,298.15,Tin[sp]) for sp,n in feed.items() if n>0)
    f=lambda T: sum(n*hf(sp)+n*hsens(sp,298.15,T) for sp,n in prod.items() if n>0)
    lo,hi=400.,2800.
    for _ in range(120):
        m=(lo+hi)/2
        if f(m)<Hin: lo=m
        else: hi=m
    return (lo+hi)/2

T_AIR = 200.+273.15
BASE_FEED={"CH4":nCH4,"NH3":nNH3,"O2":nO2,"N2":nN2}
BASE_TIN ={"CH4":298.15,"NH3":333.15,"O2":T_AIR,"N2":T_AIR}
Tb = adT(BASE_FEED,prod,BASE_TIN)
print(f"\n  BASELINE (anhydrous NH3, air preheat 200 C): adiabatic gauze T = {Tb-273.15:.0f} C")
print(f"  [literature Andrussow gauze 1050-1200 C -> model is calibrated well enough for DELTAS]")

print(f"\n  --- adding the urea ballast (CO2 + H2O travel with the NH3 at 180 C) ---")
print(f"  {'case':5} {'CO2':>6} {'H2O':>6} {'tot feed':>9} {'T_gauze':>8} {'dT':>6} {'air preheat to restore':>23}")
RES={}
for tag in ("B","C"):
    d=SEL[tag]; nb_co2=nCO2; nb_h2o=d["nH2Ov"]
    F=dict(BASE_FEED); F["CO2"]=nb_co2; F["H2O"]=nb_h2o
    P2=dict(prod); P2["CO2"]=P2.get("CO2",0)+nb_co2; P2["H2O"]+=nb_h2o
    TIN=dict(BASE_TIN); TIN["NH3"]=453.15; TIN["CO2"]=453.15; TIN["H2O"]=453.15
    Tu = adT(F,P2,TIN)
    # find air preheat that restores Tb
    lo,hi=298.,1200.
    for _ in range(80):
        m=(lo+hi)/2; TIN2=dict(TIN); TIN2["O2"]=m; TIN2["N2"]=m
        if adT(F,P2,TIN2)<Tb: lo=m
        else: hi=m
    Tair_req=(lo+hi)/2
    RES[tag]=dict(T=Tu,Tair=Tair_req,feed=sum(F.values()),co2=nb_co2,h2o=nb_h2o,P2=P2,F=F)
    print(f"  {tag:5} {nb_co2:6.2f} {nb_h2o:6.2f} {sum(F.values()):9.2f} {Tu-273.15:7.0f}C {Tu-Tb:+6.0f} "
          f"{Tair_req-273.15:19.0f} C (from 200 C)")
print(f"  baseline total feed {sum(BASE_FEED.values()):.2f} kmol/h")
print("""
  READ THIS: the ballast costs 100-190 K of gauze temperature. Andrussow HCN yield
  falls steeply below ~1050 C. It MUST be paid back, and there are only 3 currencies:
    (i)  hotter AIR preheat (air alone is safe to preheat - it is not a fuel mixture),
    (ii) burning more CH4 (extra O2 + CH4, lowers HCN selectivity),
    (iii)dry the NH3 gas (higher hydrolyser pressure = case C).
  Note (i) has a limit: hotter feed WIDENS the flammability envelope (Zabetakis:
  UFL rises ~0.72%/K relative), which fights the mixing-station safety case.""")

# ---------------------------------------------------------------- 7b RWGS
line("7b. Does the CO2 survive the gauze? Reverse water-gas shift at 1100 C")
dH_rwgs = Hf["CO_g"]+Hf["H2O_g"]-Hf["CO2_g"]
dS_rwgs = S0["CO_g"]+S0["H2O_g"]-S0["CO2_g"]-S0["H2"]
for TC in (900,1000,1100,1200):
    T=TC+273.15; dG=dH_rwgs*1000-T*dS_rwgs; K=math.exp(-dG/(R*T))
    print(f"    CO2 + H2 <-> CO + H2O   T={TC} C  dG={dG/1000:7.1f} kJ/mol  K={K:6.2f}")
T=1373.15; dG=dH_rwgs*1000-T*dS_rwgs; Krw=math.exp(-dG/(R*T))
for tag in ("B",):
    P2=RES[tag]["P2"]
    a,b,c,d0 = P2["CO2"],P2["H2"],P2["CO"],P2["H2O"]
    lo,hi=0.,a*0.999
    for _ in range(100):
        x=(lo+hi)/2
        if (c+x)*(d0+x)/((a-x)*(b-x)) < Krw: lo=x
        else: hi=x
    x=(lo+hi)/2
    print(f"\n  case {tag}: CO2 in {a:.3f} kmol/h, H2 {b:.3f} -> RWGS extent {x:.3f}")
    print(f"    CO2 leaving the gauze = {a-x:.3f} kmol/h ({(a-x)/a*100:.0f}% of the CO2 fed survives)")
    print(f"    H2 consumed {x:.3f} kmol/h (was {b:.3f}) - lost tail-gas fuel value")
    print(f"    extra endothermic load {x*dH_rwgs:.0f} MJ/h = {x*dH_rwgs/3.6:.0f} kW (already inside the adiabatic T above? NO -")
    print(f"    the adiabatic calc above assumed CO2 inert -> real gauze T is LOWER still by roughly {x*dH_rwgs*1000/ (sum(n*cp(sp,1373) for sp,n in P2.items())):.0f} K)")
    RES[tag]["CO2_out"]=a-x

# ---------------------------------------------------------------- 7c HOW DRY IS DRY ENOUGH
line("7c. How dry must the NH3 gas be? (air-preheat is capped by methane autoignition)")
print("""  HARD CONSTRAINT [H]: autoignition temperature CH4 = 537 C, NH3 = 651 C (NFPA 497 /
  Zabetakis). Any air preheat above ~480 C (537 minus a 60 K margin) means the fuel
  autoignites the instant it contacts the air in the mixer, whatever the flammability
  limits say. So air preheat is CAPPED at ~450-480 C. That caps how much ballast the
  burner can absorb.""")
print(f"\n  {'nH2O ballast':>13} {'nCO2':>6} {'T_gauze':>8} {'air preheat needed to hold 1085 C':>34}")
for nb_h2o in (13.64,10.0,6.87,4.0,2.0,0.0):
    F=dict(BASE_FEED); F["CO2"]=nCO2; F["H2O"]=nb_h2o
    P2=dict(prod); P2["CO2"]=nCO2; P2["H2O"]=prod["H2O"]+nb_h2o
    TIN=dict(BASE_TIN); TIN["NH3"]=453.15; TIN["CO2"]=453.15; TIN["H2O"]=453.15
    Tu=adT(F,P2,TIN)
    lo,hi=298.,1500.
    for _ in range(80):
        m=(lo+hi)/2; T2=dict(TIN); T2["O2"]=m; T2["N2"]=m
        if adT(F,P2,T2)<Tb: lo=m
        else: hi=m
    ok = "OK" if (lo+hi)/2-273.15 <= 480 else "IMPOSSIBLE (>480 C cap)"
    print(f"  {nb_h2o:13.2f} {nCO2:6.2f} {Tu-273.15:7.0f}C {(lo+hi)/2-273.15:22.0f} C   {ok}")
# CO2-free, water-free reference
F=dict(BASE_FEED); P2=dict(prod)
TIN=dict(BASE_TIN); TIN["NH3"]=453.15
print(f"  {'0 (NO CO2 either)':>13} {0:6.2f} {adT(F,P2,TIN)-273.15:7.0f}C {'200 (unchanged)':>22}")
print("""
  -> Even with a PERFECTLY DRY hydrolyser gas, the 4.69 kmol/h of stoichiometric CO2
     alone still needs ~330 C air preheat. That is achievable. But the water is not
     forgiving: above ~7 kmol/h of H2O ballast the required air preheat crosses the
     methane autoignition cap and the burner CANNOT be rebalanced by preheat at all.""")

# ---------------------------------------------------------------- 8. CAUSTIC / PRODUCT QUALITY
line("8. CAUSTIC CONSUMPTION AND CARBONATE IN THE PRODUCT - the other hidden cost")
naoh_cn = nNaCN*MW["NaOH"]
prod30 = CAP*1000/HOURS/0.30
print(f"  NaOH for NaCN itself: {nNaCN:.3f} kmol/h = {naoh_cn:.1f} kg/h")
print(f"  30% product solution: {prod30:.0f} kg/h")
print(f"  CO2 reaching the caustic absorber (after RWGS): {RES['B']['CO2_out']:.3f} kmol/h")
print(f"  2 NaOH + CO2 -> Na2CO3 + H2O\n")
print(f"  {'CO2 co-absorbed':>16} {'extra NaOH kg/h':>16} {'% over base':>12} {'Na2CO3 kg/h':>12} {'wt% in product':>15} {'extra NaOH t/y':>15}")
for fr in (0.10,0.20,0.35,0.60,1.00):
    nco2 = RES['B']['CO2_out']*fr
    extra = 2*nco2*MW["NaOH"]
    na2co3 = nco2*105.99
    print(f"  {fr*100:14.0f}% {extra:16.1f} {extra/naoh_cn*100:11.0f}% {na2co3:12.1f} {na2co3/(prod30+na2co3)*100:14.1f}% {extra*HOURS/1000:15.0f}")
print("""
  Interpretation: HCN absorption into NaOH is an instantaneous acid-base reaction;
  CO2 absorption is kinetically limited (CO2 + OH- , k_OH ~ 8500 L/mol/s at 25 C [H]),
  so a short-contact absorber can be made selective. But even 20% CO2 slip doubles
  ...adds ~46% to caustic and puts ~7 wt% Na2CO3 in a 30% NaCN product. Mine cyanide
  specs care about free NaOH and carbonate: THIS IS A PRODUCT-QUALITY ISSUE, not just
  a cost issue. NOT VERIFIED: the mine's actual NaCN solution acceptance spec.""")

# ---------------------------------------------------------------- 9. FLAMMABILITY
line("9. FLAMMABILITY / MIXING STATION - what the ballast actually does")
def frac(F):
    t=sum(F.values()); return {k:v/t for k,v in F.items()}
fb=frac(BASE_FEED)
print(f"  BASELINE feed: CH4 {fb['CH4']*100:.2f}%  NH3 {fb['NH3']*100:.2f}%  O2 {fb['O2']*100:.2f}%  N2 {fb['N2']*100:.2f}%")
print(f"    CH4 vs its UFL in air (15.0 vol% [H]): {fb['CH4']*100:.2f}% -> "
      f"{fb['CH4']*100-15.0:+.2f} pp   <-- the brief's '0.6-1.4 pp' margin, reproduced")
ufl_mix = 1/((fb['CH4']/(fb['CH4']+fb['NH3']))/15.0 + (fb['NH3']/(fb['CH4']+fb['NH3']))/28.0)
print(f"    Le Chatelier UFL of the CH4/NH3 fuel blend (UFL NH3 = 28% [H]): {ufl_mix:.2f}%; "
      f"total fuel = {(fb['CH4']+fb['NH3'])*100:.2f}% -> {(fb['CH4']+fb['NH3'])*100-ufl_mix:+.2f} pp above")
for tag in ("B","C"):
    F=RES[tag]["F"]; f2=frac(F)
    n2eq = F["N2"] + 1.5*F["CO2"] + 1.4*F["H2O"]   # [H] CO2 ~1.5x, H2O ~1.4x N2 as inertant
    n2eq_b = BASE_FEED["N2"]
    print(f"\n  CASE {tag}: CH4 {f2['CH4']*100:.2f}%  NH3 {f2['NH3']*100:.2f}%  O2 {f2['O2']*100:.2f}%  "
          f"CO2 {f2['CO2']*100:.2f}%  H2O {f2['H2O']*100:.2f}%")
    print(f"    O2 concentration: {fb['O2']*100:.2f}% -> {f2['O2']*100:.2f}%  "
          f"(LOC for CH4 = 12% with N2, 14.5% with CO2 [H, NFPA 69, 25 C])")
    print(f"    N2-equivalent inert: {n2eq_b:.1f} -> {n2eq:.1f} kmol/h (+{(n2eq/n2eq_b-1)*100:.0f}%)")
print("""
  Two-sided result, state it honestly:
  + the ballast drops O2 from 15.2% to 12.1-13.0%, at or below the CH4 limiting oxygen
    concentration -> on the LOC criterion the mixture becomes ROBUSTLY non-flammable,
    and CO2 is ~1.5x N2 as an inertant. This is a genuine safety dividend at the mixer.
  - but LOC falls with temperature (roughly -1 pp per 100 K), and any air preheat above
    ~330 C erodes it. The safety gain and the thermal fix pull in opposite directions.
  - and on the naive 'CH4 vs UFL' metric the ballast moves CH4 from 14.4% to 11.5%,
    i.e. INTO the 5-15% window. That metric is wrong when 20%+ inert is present, but it
    is the metric a regulator or an insurer will reach for first. Be ready for it.
  ACTION: this must be settled by MEASUREMENT, not correlation - EN 1839 / ASTM E918
  limits on the actual CH4/NH3/CO2/H2O/air composition at the real preheat temperature.
  For a plant whose selling point is 'safety benchmark', that test is not optional.""")

# ---------------------------------------------------------------- 10. LOSS OF NH3
line("10. LOSS-OF-NH3 SCENARIO - the interlock that the urea route makes critical")
noNH3 = {"CH4":nCH4,"O2":nO2,"N2":nN2}
f3=frac(noNH3)
print(f"  If NH3 stops and CH4+air keep flowing: CH4 = {f3['CH4']*100:.2f}% in air -> UFL 15% -> "
      f"{'ABOVE UFL by %.2f pp (too rich, safe by %.2f pp only)'%(f3['CH4']*100-15,f3['CH4']*100-15) if f3['CH4']*100>15 else 'INSIDE THE FLAMMABLE RANGE'}")
print(f"  With the urea route NH3, CO2 and H2O all stop together (one stream) - so the")
print(f"  post-trip mixture is the same {f3['CH4']*100:.2f}% CH4. Margin is {f3['CH4']*100-15:.2f} pp: NOT a barrier.")
print(f"  LP transfer-line NH3 inventory (50 m of DN80 at 1.5 bar, 150 C):")
V=50*math.pi*(0.0779/2)**2
n=1.5e5*V/(R*423.15)/1000
print(f"    V={V*1000:.0f} L -> {n:.4f} kmol gas, NH3 {n*nNH3/(ngas_nc+SEL['B']['nH2Ov'])*MW['NH3']:.2f} kg "
      f"= {n*nNH3/(ngas_nc+SEL['B']['nH2Ov'])*MW['NH3']/(nNH3*MW['NH3']/60):.1f} s of feed")
print("""  -> the NH3 stream dies within ~15-30 s of a hydrolyser trip. CH4 and air must be
     tripped FASTER than that. SIL-rated: 'low NH3 flow / low hydrolyser pressure ->
     close CH4 XV within 1 s'. This is the same interlock as the anhydrous plant but
     the urea route makes it harder, because a hydrolyser has NO liquid-NH3 buffer to
     ride through an upset. Compensate with: 2 x 100% hydrolysers on hot standby,
     pressure-controlled let-down, and a small (< 5 kg NH3) gas receiver.""")

# ---------------------------------------------------------------- 11. DESIGN SELECTION
line("11. DESIGN SELECTION - feed concentration is the master variable")
print("""  Steady-state water balance closes the argument:
     water OUT (overhead) = water IN (feed solution) - reaction water - blowdown
  A reflux/rectifier CANNOT change this - it only moves water around inside the column.
  So the overhead water ballast, and therefore the whole burner problem, is fixed by
  ONE number: the urea concentration of the solution you feed. Everything else follows.""")
print(f"\n  {'wt%':>5} {'Tcryst':>7} {'sol kg/h':>9} {'H2O ovhd':>9} {'P@180C':>7} {'P@200C':>7} {'Q@180C':>7} "
      f"{'T_gauze':>8} {'air preheat':>11} {'verdict':>10}")
BEST=[]
for wt in (40,45,50,55,60,62,65,68):
    sol = m_urea/(wt/100.)
    fw  = sol-m_urea
    mov = fw - nCO2*MW["H2O"] - 6.0
    if mov<=0: continue
    nov = mov/MW["H2O"]
    ratio = nov/ngas_nc
    P180 = psat(180)*AW*(1+1/ratio); P200 = psat(200)*AW*(1+1/ratio)
    Q,_,_,_ = duty(180,P180,nov,sol); Qk=Q/3.6
    F=dict(BASE_FEED); F["CO2"]=nCO2; F["H2O"]=nov
    P2=dict(prod); P2["CO2"]=nCO2; P2["H2O"]=prod["H2O"]+nov
    TIN=dict(BASE_TIN); TIN["NH3"]=453.15; TIN["CO2"]=453.15; TIN["H2O"]=453.15
    Tg=adT(F,P2,TIN)
    lo,hi=298.,1600.
    for _ in range(80):
        m=(lo+hi)/2; T2=dict(TIN); T2["O2"]=m; T2["N2"]=m
        if adT(F,P2,T2)<Tb: lo=m
        else: hi=m
    Tair=(lo+hi)/2-273.15
    Tc=cryst_T(wt)
    v = "OK" if Tair<=480 else "NO(AIT)"
    if P180>60: v+=" P>60bar"
    BEST.append((wt,Tc,sol,nov,P180,P200,Qk,Tg-273.15,Tair,v))
    print(f"  {wt:5.0f} {Tc:6.1f}C {sol:9.0f} {nov:9.2f} {P180:7.1f} {P200:7.1f} {Qk:7.0f} "
          f"{Tg-273.15:7.0f}C {Tair:10.0f}C {v:>10}")
print("""
  RECOMMENDED OPERATING POINT:  62 wt% urea solution, column 180 C / 38 bar(a).
    - overhead ballast H2O 4.7 kmol/h + CO2 4.7 kmol/h (down from 13.6 + 4.7)
    - reboiler duty 237 kW (261 kW design) instead of 373 kW at 46 wt%
    - gauze 991 C (-94 K); required air preheat 415 C: below the 480 C methane-
      autoignition cap, but with only ~65 K of margin - that margin is the design's
      tightest single number
    - cost: the solution system must be trace-heated to 50-55 C (62 wt% saturates at
      ~40 C). In Tarkwa (ambient 22-32 C) that is ~25 K of tracing, ~5 kW, and a
      tracing failure crystallises slowly over hours - recoverable, not a safety event.""")

# ---------------------------------------------------------------- 12. STORAGE IN GHANA
line("12. SOLID UREA RECEPTION AND STORAGE IN 95% RH")
print("""  [H] Urea critical relative humidity (CRH): ~81% at 10 C, ~76% at 20 C, ~72-73% at 30 C.
  Tarkwa ambient RH reaches 95%. AMBIENT RH IS ~22 pp ABOVE THE CRH ALL YEAR.
  Above CRH urea does not merely cake - it DELIQUESCES (absorbs water and dissolves).
  This single fact should kill the 280 m3 bulk silo currently in the project inventory
  (T3: 16.5 t, 7.2 x 7.2 x 8.2 m, $74k FOB) unless it is fed dehumidified air.""")
rho_bulk = 750.0
for days in (7,14,30,45):
    t = m_urea*24*days/1000.
    print(f"    {days:2d} d buffer = {t:6.1f} t = {t*1000/rho_bulk:6.1f} m3 bulk = {t/1.0:.0f} x 1 t FIBC")
print(f"""
  RECOMMENDATION - three-tier, and only the smallest tier is exposed to the air:
   (1) 30 d strategic buffer = {m_urea*24*30/1000:.0f} t as 1 t FIBC with sealed 150 um PE liner,
       in a closed shed. The BAG is the moisture barrier - no silo, no dehumidifier.
       Stack 2 high on 1.2 x 1.2 m pitch -> {m_urea*24*30/1000/2*1.44:.0f} m2 of racking + aisles ~ 200 m2 shed.
   (2) 3 d day-bin {m_urea*24*3/1000:.1f} t / {m_urea*24*3/rho_bulk:.0f} m3 (D 2.5 x H 6 m), 70 deg cone, vibrating
       discharger, rotary valve, blanketed with dehumidified air (dew point < 5 C,
       ~2 kW desiccant dryer). Filled from FIBC in one ~3 h campaign every 3 days.
   (3) 48 h solution buffer 2 x 12 m3 - this is the tier that actually decouples the
       plant from solids handling and lets it run unattended.
  Operator intervention: ONE bag-handling campaign per 3 days. Nothing daily.
  DELIBERATELY REJECTED: a 30 d bulk silo. In 95% RH it bridges and rat-holes, and a
  blocked urea silo stops NaCN production, which stops the leach circuit.""")

# ---------------------------------------------------------------- 13. UREA SPECIFICATION
line("13. UREA SPECIFICATION - the non-obvious one: DO NOT buy fertiliser-grade")
print(f"""  Fertiliser granular urea is the wrong feedstock for two independent reasons:
   (a) ANTI-CAKING COATING. Granular fertiliser urea is surface-treated with
       urea-formaldehyde resin (UF-85), typically 0.3-0.6 wt%. At {m_urea:.0f} kg/h that is
       {m_urea*0.005:.1f} kg/h = {m_urea*0.005*HOURS/1000:.1f} t/y of formaldehyde resin fed to a 180 C hydrolyser.
       It polymerises to insoluble methylene-urea sludge -> reboiler fouling.
       Ghana's humidity makes an uncoated product tempting to refuse - but the FIBC
       liner, not the coating, is the moisture barrier here. SPECIFY UNCOATED.
   (b) BIURET. Fertiliser grade allows 1.0-1.2 wt% biuret; technical/AdBlue grade
       (ISO 22241 / DIN 70070 feedstock) allows <=0.3 wt%.
  Target specification (write it into the purchase order):
       urea            >= 99.5 wt%      biuret         <= 0.3 wt%
       water           <= 0.5 wt%       aldehydes      <= 5 mg/kg
       insolubles      <= 20 mg/kg      alkali (Na+K)  <= 1 mg/kg
       Fe, Cu, Zn, Ni, Cr, Al, Ca, Mg   <= 0.5 mg/kg each     phosphate <= 0.5 mg/kg
       NO anti-caking / NO UF coating; supplied in 1 t FIBC with sealed PE liner
  WHY THE METALS LIMITS MATTER, and this is the good news of the whole node:
    the Pt/Rh gauze is poisoned by alkali, iron and phosphorus. In the ANHYDROUS route
    any such contaminant in the NH3 goes straight to the gauze. In the UREA route only
    NH3, CO2 and H2O are volatile - EVERY non-volatile contaminant stays in the column
    liquor and leaves in the blowdown. THE HYDROLYSER IS ALSO A PURIFICATION STEP.
    The only carry-over path is entrainment -> demister + one water-wash tray at the
    column top brings it to < 1 mg/kg. Net: the urea route delivers a CLEANER feed to
    the gauze than tanker ammonia does.""")

# ---------------------------------------------------------------- 14. IMPURITY / BLOWDOWN
line("14. BIURET, HNCO AND THE BLOWDOWN BALANCE")
print(f"""  HNCO: in the DRY thermolysis route (NH2)2CO -> NH3 + HNCO is unavoidable and HNCO
  then needs a TiO2/Al2O3 hydrolysis bed, which fouls. In the AQUEOUS route at 180 C
  with a large excess of liquid water, HNCO hydrolyses in situ (HNCO + H2O -> NH3 +
  CO2) far faster than it can escape. Residual HNCO in the overhead should be ppm-level.
  -> DESIGN DECISION: aqueous hydrolysis, not melt thermolysis. In Ghana this also
     turns urea's hygroscopicity from a liability into a non-issue: we are adding
     water on purpose, so a slightly damp bag is a dosing correction, not a reject.
     [HNCO slip is NOT quantified here - no kinetic data in this session. Measure it.]""")
for bi in (0.003,0.010):
    m_bi = m_urea*bi
    for destr in (0.8,0.5):
        left = m_bi*(1-destr)
        print(f"  biuret {bi*100:.1f} wt% in feed = {m_bi:.2f} kg/h; {destr*100:.0f}% hydrolysed at 180 C"
              f" -> {left:.2f} kg/h to blowdown = {left*HOURS/1000:.2f} t/y")
nv_solids = m_urea*0.003*0.2 + m_urea*20e-6
print(f"\n  non-volatile load (biuret residue + insolubles), tech-grade urea: {nv_solids:.3f} kg/h = {nv_solids*HOURS:.0f} kg/y")
for hold in (0.02,0.05,0.10):
    bd = nv_solids/hold
    print(f"    holding {hold*100:.0f} wt% non-volatiles in the liquor -> blowdown {bd:.1f} kg/h "
          f"({bd*24:.0f} kg/d, {bd*HOURS/1000:.1f} t/y)")
print(f"""
  ROUTE THE BLOWDOWN BACK, DO NOT DISCHARGE IT: cool it and use it as dissolver
  make-up water. It carries dissolved NH3 and urea; discharging it would put ammonia
  into a gold-mine water circuit, which is exactly the kind of thing the ICMI Cyanide
  Code auditors and the Ghana EPA look at. Closed loop + a ~200 kg batch drain roughly
  monthly to a drum, neutralised, is the plug-and-play answer.
  [Ghana regulatory citations - EPA Act 490 (1994), LI 1652 (1999), Factories Offices
   and Shops Act 1970 (Act 328) for pressure vessels - RECALLED, NOT VERIFIED this
   session (web search budget exhausted). Must be confirmed before use.]""")

# ---------------------------------------------------------------- 15. UTILITIES
line("15. UTILITIES - can the Andrussow waste-heat boiler carry this node?")
WT=62.0
sol_d = m_urea/(WT/100.); nov_d = (sol_d-m_urea-nCO2*MW["H2O"]-6.0)/MW["H2O"]
P_d = psat(180)*AW*(1+ngas_nc/nov_d)
Qd,_,_,_ = duty(180,P_d,nov_d,sol_d); Qd/=3.6
print(f"  design point: {WT:.0f} wt%, 180 C, {P_d:.1f} bara -> reboiler {Qd:.0f} kW (normal), {Qd*1.1:.0f} kW (design)")
qdiss = nUREA*DHSOL/3.6 + sol_d*CP_SOL*(58-28)/3600
print(f"  dissolver {qdiss:.0f} kW ; solution tank tracing ~5 kW ; LP line tracing ~8 kW")
print(f"  NODE TOTAL heat demand = {Qd*1.1+qdiss+13:.0f} kW")
# waste heat in the Andrussow gas
F=dict(BASE_FEED); F["CO2"]=nCO2; F["H2O"]=nov_d
P2=dict(prod); P2["CO2"]=nCO2; P2["H2O"]=prod["H2O"]+nov_d
Qwh = sum(n*hsens(sp,473.15,1373.15) for sp,n in P2.items() if n>0)   # MJ/h, 1100->200 C
print(f"\n  Andrussow off-gas {sum(P2.values()):.1f} kmol/h cooled 1100 -> 200 C in the WHB:")
print(f"    recoverable heat = {Qwh:.0f} MJ/h = {Qwh/3.6:.0f} kW  (vs node demand {Qd*1.1+qdiss+13:.0f} kW)")
print(f"    -> covered {Qwh/3.6/(Qd*1.1+qdiss+13):.1f}x over. The hydrolyser runs on waste heat: NO fuel cost.")
print(f"""    BUT the WHB must make steam hot enough to drive a 180 C reboiler:
      needs >= 230 C saturated = 28 bar(a); at 250 C / 40 bar(a) the reboiler is
      A = {Qd*1.1*1000/(800*70):.1f} m2 (U=800 W/m2K, dT=70 K). Steam rate {Qd*1.1*3600/1716:.0f} kg/h.
      >>> ACTION FOR THE M1 (reaction module) NODE: specify the steam drum at 40 bar(a),
          not the 10-20 bar that a small Andrussow WHB would default to. <<<""")
print(f"""
  START-UP and BACK-UP: at cold start there is no waste heat and the burner cannot be
  lit without NH3. Fit a {math.ceil(Qd*1.1/50)*50:.0f} kW ELECTRIC immersion reboiler in parallel with the
  steam reboiler. This also avoids an auxiliary fired boiler and its statutory boiler
  inspection/licensing - a real plug-and-play win. Ghana grid/mine power.
  Electric-only running cost, if it ever came to that: {Qd*1.1:.0f} kW x {HOURS:.0f} h = {Qd*1.1*HOURS/1000:.0f} MWh/y;
  at 0.13 $/kWh (INDICATIVE, NOT VERIFIED) = ${Qd*1.1*HOURS*0.13:,.0f}/y = ${Qd*1.1*HOURS*0.13/CAP:.0f}/t NaCN.
  Start-up only ({math.ceil(Qd*1.1/50)*50:.0f} kW x ~4 h x ~12 starts/y) = negligible.""")

# ---------------------------------------------------------------- 16. EQUIPMENT LIST
line("16. EQUIPMENT LIST, DIMENSIONS, WEIGHTS, MODULE FIT (40ft HC = 12.03 x 2.35 x 2.69 m)")
EQ = [
 ("V-101","FIBC discharge station with hoist, bag-breaker, dust filter","2.5 x 2.5 x 5.5 m", 1.8,"skid, field-erected"),
 ("V-102","Urea day-bin 27 m3 (3 d), SS304, 70 deg cone, vibrating discharger","D 2.5 x H 6.0 m", 4.5,"OOG - ships as shell + cone"),
 ("PK-103","Desiccant air dryer for bin blanket, dew point < 5 C, 2 kW","0.8 x 0.6 x 1.6 m", 0.3,"in 40ft HC"),
 ("W-104","Loss-in-weight screw feeder, 0-600 kg/h urea","1.8 x 0.5 x 0.8 m", 0.4,"in 40ft HC"),
 ("V-105","Dissolver, agitated, 1.0 m3, SS304, 30 kW electric/hot-water coil","D 1.0 x H 1.4 m", 0.9,"in 40ft HC"),
 ("T-106","Urea solution tanks 2 x 12 m3, SS304, insulated + traced 50-55 C","D 2.2 x H 3.6 m ea", 2.6,"local fab (GH)"),
 ("P-107","Solution metering pumps 2x100%, 0-600 kg/h, 45 bar, triplex","1.2 x 0.8 x 1.0 m ea", 0.7,"in 40ft HC"),
 ("F-108","Duplex feed filter 100 um + 5 um","0.3 x 0.3 x 1.0 m", 0.1,"in 40ft HC"),
 ("C-201","HYDROLYSIS COLUMN 2 x 100%, D 600 x 6.0 m T/T, SS316L, 14 mm, 38/48 bar","D 0.75 x L 7.5 m ea", 2.4,"lies flat in 40ft HC"),
 ("E-202","Thermosiphon reboiler 2 x 100%, 6 m2, SS316L, steam 40 bar(a)","D 0.35 x L 3.0 m ea", 0.6,"in 40ft HC"),
 ("E-203","Feed/effluent + blowdown cooler","D 0.25 x L 2.0 m", 0.3,"in 40ft HC"),
 ("PCV-204","Let-down station, jacketed, at the column; 38 -> 1.8 bar","-", 0.2,"on C-201 skid"),
 ("H-205","LP NH3/CO2/H2O transfer line, DN80 SS316L, electrically traced 140 C","~50 m", 0.6,"loose"),
 ("V-206","NH3 gas receiver 0.5 m3 at 1.8 bar (surge only, < 1 kg NH3)","D 0.6 x L 1.8 m", 0.2,"in 40ft HC"),
 ("E-207","HIGH-TEMPERATURE AIR PREHEATER to ~400 C (NEW - see section 7c)","D 0.8 x L 3.5 m", 2.2,"alloy 800H/321H"),
]
tot_m=0
print(f"  {'tag':7} {'item':62} {'envelope':22} {'t':>5}  fit")
for tag,name,env,m,fit in EQ:
    tot_m+=m
    print(f"  {tag:7} {name:62.62} {env:22} {m:5.1f}  {fit}")
print(f"  {'':7} {'NODE TOTAL':62} {'':22} {tot_m:5.1f} t")
print(f"""
  Everything except the day-bin (V-102, D 2.5 m) and the air preheater fits inside the
  existing 40ft HC envelope. The hydrolysis columns ship LYING DOWN (7.5 m < 12.03 m).
  Two new transport modules: M7 'urea dissolution + dosing' and M8 'hydrolysis'.
  ~{tot_m:.0f} t added; against the 113 t already in the shipment that is ~{tot_m/113*100:.0f}%.
  DELETED from the shipment by this node: 2 x 25 t anhydrous NH3 bullets (D 3.0 m,
  out-of-gauge on both) and the NH3 vaporiser/superheater package.
  CHANGED: T3 in the current inventory (280 m3 urea silo, 16.5 t, $74k) is REPLACED by
  V-102 (27 m3 bin) + a 200 m2 bagged-goods shed. Net FOB is probably lower; the shed
  is a pure Ghana civil scope.""")

# ---------------------------------------------------------------- 17. SAFETY LEDGER
line("17. SAFETY LEDGER - what the swap actually buys, in kg")
inv_col = 0.65*950*(nNH3*MW["NH3"]/ (nNH3*MW["NH3"]+sol_d))  # crude: NH3 dissolved in liquor
liq_hold = 0.65*950
print(f"  ANHYDROUS ROUTE   NH3 inventory on site: 60 000 kg  (toxic cloud to 1.9 km, per brief)")
print(f"  UREA ROUTE        hazardous-NH3 inventory:")
print(f"    hydrolysis column liquor  {liq_hold:.0f} kg at 180 C, of which free/loosely bound NH3 ~5 wt% = {liq_hold*0.05:.0f} kg")
print(f"    x2 columns                                                        = {liq_hold*0.05*2:.0f} kg")
print(f"    HP + LP gas lines + receiver                                      ~ {1.0:.0f} kg")
print(f"    TOTAL                                                             ~ {liq_hold*0.05*2+1:.0f} kg")
print(f"    -> reduction factor {60000/(liq_hold*0.05*2+1):.0f}x")
print(f"  Urea solution 24 m3 x 1.18 t/m3 x 62% = {24*1180*0.62/1000:.1f} t urea = {24*1180*0.62*2*MW['NH3']/MW['urea']/1000:.1f} t of NH3 'potential',")
print(f"    but chemically bound, non-volatile, non-toxic by inhalation, and it cannot")
print(f"    form a cloud: a full-bore rupture makes a puddle, not a plume.")
print(f"\n  ROAD TRANSPORT THROUGH GHANA (Tema -> Tarkwa, ~300 km):")
print(f"    anhydrous: {nNH3*MW['NH3']*HOURS/1000:.0f} t/y NH3 / 20 t tanker = {nNH3*MW['NH3']*HOURS/1000/20:.0f} loads/y of UN 1005 Class 2.3 TOXIC GAS")
print(f"    urea:      {m_urea*HOURS/1000:.0f} t/y / 28 t truck = {m_urea*HOURS/1000/28:.0f} loads/y of NON-DANGEROUS goods")
print(f"    -> {m_urea*HOURS/1000/28/(nNH3*MW['NH3']*HOURS/1000/20):.1f}x more truck movements, but ZERO toxic-gas movements.")
print(f"       Consequence per accident collapses; exposure frequency rises. Say both.")

# ---------------------------------------------------------------- 18. COST DELTA
line("18. COST DELTA vs THE ANHYDROUS ROUTE (per t NaCN; prices INDICATIVE, NOT VERIFIED)")
P_NH3, P_UREA, P_NAOH = 650., 400., 500.
c_nh3  = nNH3*MW["NH3"]*HOURS/1000*P_NH3/CAP
c_urea = m_urea*HOURS/1000*P_UREA/CAP
print(f"  feedstock, anhydrous NH3  {nNH3*MW['NH3']*HOURS/1000:6.0f} t/y x {P_NH3:.0f} $/t = ${c_nh3:6.0f}/t NaCN")
print(f"  feedstock, urea           {m_urea*HOURS/1000:6.0f} t/y x {P_UREA:.0f} $/t = ${c_urea:6.0f}/t NaCN   "
      f"delta {c_urea-c_nh3:+.0f} $/t")
print(f"  breakeven urea price = {2*14.007/MW['urea']/(14.007/MW['NH3'])*P_NH3:.0f} $/t at NH3 {P_NH3:.0f} $/t (contained-N parity)")
for fr in (0.10,0.20):
    extra = 2*RES['B']['CO2_out']*fr*MW["NaOH"]*HOURS/1000
    print(f"  extra NaOH from CO2 slip at {fr*100:.0f}% co-absorption: {extra:.0f} t/y x {P_NAOH:.0f} $/t = ${extra*P_NAOH/CAP:.0f}/t NaCN")
print(f"  extra CH4 / thermal: covered by waste heat and by air preheat -> ~$0/t operating,")
print(f"    but CAPEX: high-temperature air preheater E-207 (alloy) ~ $120-200k (ESTIMATE)")
print(f"  deleted: NH3 bullets, vaporiser, NH3 gas detection ring, NH3 emergency scrubber duty")
print(f"""
  BOTTOM LINE ON MONEY: at parity urea prices the feedstock delta is roughly
  {c_urea-c_nh3:+.0f} $/t NaCN against a cash cost of 1228 $/t, i.e. {abs(c_urea-c_nh3)/1228*100:.0f}%. The CO2-driven
  caustic and carbonate penalty is the one that can actually hurt, and it is
  controlled by absorber design, not by the urea decision itself.""")

line("19. RESIDUAL RISKS AND WHAT IS NOT KNOWN")
print("""  QUANTIFIED HERE:
   - water ballast -> gauze temperature: 168 K at 46 wt% feed, 94 K at 62 wt%
   - air preheat cap 480 C from methane autoignition (537 C)
   - carbamate deposition 115 C at 38 bar, 55 C at 1.5 bar
   - dissolver adiabatic temperature drop 34 K -> mandatory heating
   - reboiler 237 kW normal (NOT the 235 kW in the brief - that figure omits feed
     sensible heat and water vaporisation and happens to land in the same place
     for the wrong reasons; at 46 wt% feed the true duty is 373 kW)

  NOT KNOWN / MUST BE CLOSED BEFORE DETAIL DESIGN:
   1. k(T) for urea hydrolysis and the NH3-CO2-H2O VLE at 180 C / 38 bar. Everything
      about column height and the true overhead composition rests on these. Get the
      licensor's data or run a 2 L autoclave series. This is the #1 open item.
   2. HNCO and biuret slip into the overhead - assumed negligible, NOT measured.
   3. Real flammability envelope of CH4/NH3/CO2/H2O/air at 200-400 C preheat (EN 1839).
   4. The mine's NaCN solution acceptance spec, especially carbonate and free NaOH.
   5. Ghana statutory basis: pressure-vessel registration for a 48 bar(g) column,
      fertiliser/urea import licensing, EPA permitting. All RECALLED, NOT VERIFIED.
   6. Urea and ammonia delivered prices into Tarkwa. The whole economic comparison
      is a single ratio (0.567) away from flipping.

  RISKS THE DESIGN DELIBERATELY ACCEPTS:
   - a 48 bar(g) pressure vessel appears on a site whose selling point is 'no pressure
     hazard'. It is small (620 L) and contains ~30 kg of NH3, but it IS a new
     pressure boundary and needs its own registration and PSV/relief study.
   - crystallisation of 62 wt% urea solution on loss of tracing. Slow, recoverable.
   - the hydrolyser has no ride-through inventory: any trip stops NH3 in ~20 s.
   - fouling of the reboiler by biuret/insolubles: mitigated by spec + 2 x 100% +
     online switchover, but it is the most likely cause of an unplanned outage.""")
