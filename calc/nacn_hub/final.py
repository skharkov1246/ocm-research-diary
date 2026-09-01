import math
Q0,CAPEX0=2500.0,10.22e6
ITEMS=[("мочевина",506,584,1.0),("едкий натр",512,537,1.0),("NaOH на CO2",56,112,1.0),
 ("газ",104,302,1.0),("электро",57,122,1.0),("вода",9,29,1.0),("Pt/Rh",20,40,1.0),
 ("стоки",22,38,1.0),("персонал",166,300,0.30),("ТОиР",143,204,None),
 ("страхование",35,55,None),("ICMI/лаб",30,48,0.25)]
CRF=lambda r,n: r/(1-(1+r)**-n); K=CRF(0.12,15); AF=1/K
def cx(Q,n,pen=0.0): return CAPEX0*(Q/Q0)**n*(1+(pen if Q>6000 else 0))
def cash(Q,n,idx):
    s=0
    for _,l,h,e in ITEMS: s+=(l if idx==0 else h)*(Q/Q0)**((n if e is None else e)-1)
    return s
# узел налива автоцистерн — то, что было СНЯТО ПО БЕЗОПАСНОСТИ, и возвращается вместе с плечом
LOAD_LO,LOAD_HI=1.46e6,3.10e6
def trip(km,lo):
    d=2*km; f=d*0.35*(1.05 if lo else 1.25); c=20 if lo else 35
    e=(50 if lo else 90)+d*(0.15 if lo else 0.30); t=(60 if lo else 110)+d*(0.35 if lo else 0.60)
    return (f+c+e+t)*(1.18 if lo else 1.30)/6.78
# конфигурации: (имя, Q, [(кт, км)] внеплощадочных объёмов)
CFG=[("A база: только Тарква, доля рынка",2500,[]),
     ("B Тарква одна, плечо 0 км",        6000,[]),
     ("C +Идуаприем 10 км",               8000,[(2.0,10)]),
     ("D +Даманг 30 км, Васса 35 км",    10000,[(2.0,10),(1.8,30),(0.5,35)]),
     ("E +Эдикан 85 км",                 12000,[(2.0,10),(1.8,30),(0.5,35),(1.9,85)])]
print(f"=== ИНТЕГРИРОВАННАЯ МОДЕЛЬ (аннуитет 12 проц. / 15 лет, CRF={K:.5f}) ===")
hdr=f"{'конфигурация':34s}{'Q':>6s}{'капитал млн$':>15s}{'достав.$/т':>12s}{'ПОЛНАЯ+дост.$/т':>18s}"
print(hdr); print("-"*len(hdr))
R={}
for nm,Q,offs in CFG:
    cl=cx(Q,0.62); ch=cx(Q,0.72,0.18)
    if offs: cl+=LOAD_LO; ch+=LOAD_HI
    dl=sum(kt*1000*trip(km,True)  for kt,km in offs)/Q
    dh=sum(kt*1000*trip(km,False) for kt,km in offs)/Q
    fl=cash(Q,0.62,0)+cl*K/Q+dl; fh=cash(Q,0.72,1)+ch*K/Q+dh
    R[nm]=(Q,cl,ch,fl,fh,dl,dh,cash(Q,0.62,0),cash(Q,0.72,1))
    print(f"{nm:34s}{Q:6d}{cl/1e6:8.1f}-{ch/1e6:6.1f}{dl:7.0f}-{dh:4.0f}{fl:12.0f}-{fh:5.0f}")
base=R["A база: только Тарква, доля рынка"]
print("\n=== ЭФФЕКТ ОТНОСИТЕЛЬНО БАЗЫ 2500 т/год ===")
for nm,(Q,cl,ch,fl,fh,dl,dh,_,_) in R.items():
    if Q==2500: continue
    print(f"  {nm:34s} эффект {fl-base[3]:+6.0f} ... {fh-base[4]:+6.0f} $/т | капитал {(cl-CAPEX0)/1e6:+5.1f}...{(ch-CAPEX0)/1e6:+5.1f} млн")

print("\n=== ЦЕНА КОНЦЕНТРАЦИИ: сдвиг реализуемой цены (переговорный излишек) ===")
PAR,BRW=2207.0,2596.0; S=BRW-PAR
for lab,slo,shi,q in (("2500 т (13% нужд одного рудника)",0.40,0.60,2500),
                      ("8-12 кт (42-63% рынка Ганы, 1-2 покупателя)",0.15,0.30,10000)):
    print(f"  {lab:44s} цена {PAR+slo*S:.0f}-{PAR+shi*S:.0f} $/т")
print(f"  => потеря реализуемой цены при масштабе: {PAR+0.40*S-(PAR+0.15*S):.0f}...{PAR+0.60*S-(PAR+0.30*S):.0f} $/т [Д]")

print("\n=== НЕТТО: эффект масштаба МИНУС цена концентрации ===")
for nm in ("B Тарква одна, плечо 0 км","C +Идуаприем 10 км","D +Даманг 30 км, Васса 35 км","E +Эдикан 85 км"):
    Q,cl,ch,fl,fh,*_=R[nm]
    print(f"  {nm:34s} {fl-base[3]+97:+6.0f} ... {fh-base[4]+117:+6.0f} $/т нетто")

print("\n=== ОКУПАЕМОСТЬ (простая, EBITDA = цена - денежная себестоимость) ===")
for price,plab in ((2207,"паритет с импортом"),(2757,"greenfield-середина")):
    print(f"  цена {price} $/т ({plab}):")
    for nm,(Q,cl,ch,fl,fh,dl,dh,ml,mh) in R.items():
        for cap,cc,dd,tag in ((cl,ml,dl,"низ"),(ch,mh,dh,"верх")):
            eb=(price-cc-dd)*Q
            pb=cap/eb if eb>0 else float('inf')
            print(f"     {nm:34s} {tag:4s} EBITDA {eb/1e6:6.2f} млн  окуп. {('%.1f г'%pb) if pb<60 else 'НЕ ОКУПАЕТСЯ'}")
