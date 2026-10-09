# NaCN Nigeria transfer model. Bases: routes/nacn_s1/audit_cost_out.txt, final_econ_out.txt,
# own_feedstock_out.txt, calc/nacn_hub/scale_out.txt, calc/nacn_series/series_verify_out.txt
def npv(rate, flows): return sum(f/(1+rate)**t for t,f in enumerate(flows))
def irr(flows):
    lo,hi=-0.99,1.0
    if npv(lo+1e-9,flows)*npv(hi,flows)>0: return None
    for _ in range(200):
        m=(lo+hi)/2
        if npv(lo,flows)*npv(m,flows)<=0: hi=m
        else: lo=m
    return m
def proj(capex, ebitda, life=15, build=(0.4,0.6), ramp=0.6):
    f=[-capex*build[0], -capex*build[1]]
    for y in range(life):
        f.append(ebitda*(ramp if y==0 else 1.0))
    return f

def option_b(Q, case):  # Q kt/y solid NaCN in Delta, export to region
    lo = case=='opt'
    s=Q/2.5
    # CAPEX, $M
    base = 10.22*s**0.62 if lo else 11.64*s**0.72*(1.18 if Q>2.5 else 1.0)
    solid = base*(0.35 if lo else 0.60)
    loc = 1.51*s**0.7*(0.8 if lo else 1.5)
    gascond = (1.0 if lo else 2.5)*s**0.6
    avgMW = 0.19*s + 0.05*Q
    power = avgMW*2.5*(0.9 if lo else 1.4)
    sec = 0.3 if lo else 0.8
    tie = 0.5 if lo else 2.0
    capex = base+solid+loc+gascond+power+sec+tie
    # variable $/t
    urea = (0.929*400) if lo else (1.071*520)
    alk = ((0.857+0.094)*680) if lo else ((0.90+0.188)*830)
    gas = (37*0.25) if lo else (42*1.0)
    pwr = (700*0.015) if lo else (850*0.025)
    water = 10 if lo else 30
    cat = 20 if lo else 40
    ins = 50 if lo else 100
    icmi = 40 if lo else 78
    eff = 22 if lo else 38
    pack = 50 if lo else 110
    var = urea+alk+gas+pwr+water+cat+ins+icmi+eff+pack
    # fixed $M/y
    nloc = (18 if lo else 24)*s**0.4
    staff = nloc*(0.012 if lo else 0.025) + (1*0.30 if lo else 2*0.40)
    secop = (0.4 if lo else 1.0)*s**0.3
    toir = capex*(0.035 if lo else 0.05)
    fixed = staff+secop+toir
    # price/netback
    price = 2240*0.97+130 if lo else 2036*0.92
    freight = 240 if lo else 500
    netback = price-freight
    util = 0.95 if lo else 0.80
    sold = Q*1000*util
    opex = sold*var/1e6 + fixed
    community = 0.0 if lo else 0.03*opex
    ebitda = sold*netback/1e6 - opex - community
    cash_per_t = (opex+community)*1e6/sold
    return dict(Q=Q,case=case,capex=capex,base=base,solid=solid,var=var,fixed=fixed,staff=staff,secop=secop,toir=toir,
                cash_per_t=cash_per_t,price=price,netback=netback,rev=sold*netback/1e6,gross_rev=sold*price/1e6,
                opex=opex+community,ebitda=ebitda,sold=sold)

def option_a(case):  # 400 t/y Ca(CN)2/NaCN solution by pipe to Segilola CIL, gas by CNG
    lo = case=='opt'
    Q = 0.30 if not lo else 0.40   # kt/y, Segilola need 250-400
    capex = 4.0 if lo else 7.0     # scaled 3.2-3.3 by n=0.62-0.72, + non-scaling safety/permits
    urea = 0.929*(480+40) if lo else 1.071*(480+60)
    alk = 0.572*1000*0.27*1.3 if lo else (0.90+0.188)*(770+80)  # lime route (opt) vs NaOH (pess)
    gas = 27.5*8.0 if lo else 27.5*11.0
    pwr = 640*0.10 if lo else 640*0.16
    other = 20+35+30+22 if lo else 40+55+48+38
    var = urea+alk+gas+pwr+other
    staff = 14*0.012+0.25 if lo else 14*0.025+0.40
    toir = capex*(0.035 if lo else 0.05)
    fixed = staff+toir
    price = 2160+300 if lo else 1675   # gate (import CIF NG + 125-217 delivery) + avoided after-gate
    sold = Q*1000
    ebitda = sold*price/1e6 - sold*var/1e6 - fixed
    return dict(Q=Q,capex=capex,var=var,fixed=fixed,staff=staff,cash_per_t=(sold*var/1e6+fixed)*1e6/sold,
                price=price,rev=sold*price/1e6,ebitda=ebitda)

print("=== (a) Segilola, раствор по трубе, газ CNG ===")
for c in ('opt','pess'):
    r=option_a(c); f=proj(r['capex'],r['ebitda'],life=4)
    print(c, {k:round(v,3) for k,v in r.items()}, 'NPV15(4y life)=%.2f'%npv(0.15,f))
print()
print("=== (b) твёрдый NaCN в Дельте на экспорт ===")
for Q in (2.5,5,10):
    for c in ('opt','pess'):
        r=option_b(Q,c); f=proj(r['capex'],r['ebitda'])
        pb = r['capex']/r['ebitda'] if r['ebitda']>0 else float('inf')
        i=irr(f)
        print(f"Q={Q:>4} {c:4} capex={r['capex']:.2f} (base {r['base']:.2f}, solid {r['solid']:.2f}) var={r['var']:.0f} $/t "
              f"fixed={r['fixed']:.2f} (staff {r['staff']:.2f}, sec {r['secop']:.2f}, toir {r['toir']:.2f}) cash={r['cash_per_t']:.0f} $/t "
              f"price={r['price']:.0f} netback={r['netback']:.0f} rev={r['rev']:.2f} grossrev={r['gross_rev']:.2f} EBITDA={r['ebitda']:.2f} "
              f"payback={pb:.1f} NPV15={npv(0.15,f):.2f} NPV20={npv(0.20,f):.2f} IRR={'n/a' if i is None else round(i*100,1)}")
print()
# breakeven netback for 10kt opt to get payback 5y
for Q in (2.5,5,10):
  for c in ('opt','pess'):
    r=option_b(Q,c)
    need = (r['capex']/5 + r['opex'])*1e6/r['sold']
    need0 = r['opex']*1e6/r['sold']
    print(f"Q={Q} {c}: netback для EBITDA=0: {need0:.0f}; для окупаемости 5 лет: {need:.0f}; цена CIF-экв. (+фрахт {240 if c=='opt' else 500}): {need+(240 if c=='opt' else 500):.0f}")
print()
# gas use
for Q in (2.5,5,10):
    mmbtu = Q*1000*np_ if (np_:=37) else 0
    for g in (37,42):
        mmscfd = Q*1000*g/1.05/1000/365
        print(f"Q={Q} kt: газ {g} MMBtu/t -> {mmscfd:.2f} MMscf/d; доля 1/5/15 MMscf/d: {mmscfd/1*100:.0f}% / {mmscfd/5*100:.0f}% / {mmscfd/15*100:.1f}%")
print()
# benchmarks per 1 MMscf/d
ph=265.1
for hp in (27.7,50):
    gross=ph*hp*365/1e6; opex=4.84-3.5
    print(f"майнинг hashprice {hp}: выручка {gross:.2f}, EBITDA ~{gross-opex:.2f} $M/y, capex 9.02, простая окуп. {9.02/max(gross-opex,1e-9):.1f} г (ASIC ~2 г)")
gwh=4.64*8760*0.85/1000
for p in (0.12,0.15,0.20):
    e=gwh*(p-0.02-0.006)
    print(f"captive power {p}$/kWh: {gwh:.1f} GWh/y, EBITDA {e:.2f} $M/y, capex 4.6-7.0 -> окуп {4.6/e:.1f}-{7.0/e:.1f} г")
# NH3 own lever at 10kt
nh3 = 10000*0.95*0.511/1000, 10000*0.95*0.589/1000
print("NH3 нужда 10 кт NaCN, т/сут:", [round(x*1000/350,1) for x in nh3])
