import math
# NIST Shomate (J/mol/K), t=T/1000
CH4=[(298,1300,[-0.703029,108.4773,-42.52157,5.862788,0.678565,-76.84376,158.7163,-74.8731])]
H2=[(298,1000,[33.066178,-11.363417,11.432816,-2.772874,-0.158558,-9.980797,172.707974,0.0]),
    (1000,2500,[18.563083,12.257357,-2.859786,0.268238,1.977990,-1.147438,156.288133,0.0])]
C=[(298,6000,[21.17510,-0.812428,0.448537,-0.043256,-0.013103,710.3470,183.8734,716.6700])] # gas! not graphite
def HS(sh,T):
    for lo,hi,c in sh:
        if lo-1<=T<=hi+1:
            A,B,Cc,D,E,F,G,Hh=c;t=T/1000
            H=A*t+B*t*t/2+Cc*t**3/3+D*t**4/4-E/t+F-Hh   # kJ/mol (H-H298)
            S=A*math.log(t)+B*t+Cc*t*t/2+D*t**3/3-E/(2*t*t)+G
            return H,S
    raise ValueError(T)
# graphite: Cp ~ approx (Butland & Maddison) use integrated simple fit: Cp=16.86+0.00477T-8.54e5/T^2 (J/molK, 298-2500)
def Hgr(T):
    f=lambda T:16.86*T+0.00477*T*T/2+8.54e5/T
    return (f(T)-f(298.15))/1000
def Sgr(T):
    g=lambda T:16.86*math.log(T)+0.00477*T+8.54e5/(2*T*T)
    return 5.74+(g(T)-g(298.15))
dHf_CH4=-74.87
S298={'CH4':186.25,'H2':130.68}
for T in (1323,1373):
    hc,sc=HS(CH4,min(T,1300)); # extrapolate CH4 above 1300 linearly by Cp
    if T>1300:
        h1,_=HS(CH4,1299);h2,_=HS(CH4,1300);cp=(h2-h1)*1000
        hc=h2+cp*(T-1300)/1000; sc=sc+cp*math.log(T/1300)
    hh,sh=HS(H2,T)
    dH=74.87+Hgr(T)+2*hh-hc
    dS=Sgr(T)+2*sh-sc
    dG=dH-T*dS/1000
    K=math.exp(-dG*1000/(8.314*T))
    print(f"T={T} dH={dH:.1f} kJ dS={dS:.1f} dG={dG:.1f} Kp={K:.0f} bar; CH4 sensible 298->T {hc:.1f} kJ")
    for P in (1,5,7.5,10,20):
        a=K/P; X=math.sqrt(a/(4+a)); print(f"   P={P} bar Xeq={X:.3f}")
    for X in (0.75,0.87,0.95):
        # heat absorbed, cold feed 298K, products leave at T, no recovery
        Q=X*dH+hc  # heat CH4 to T then react fraction X at T
        # with feed/effluent recovery: assume feed preheated to 873K by hot product gas
        h873,_=HS(CH4,873); Qr=Q-h873
        print(f"   X={X} Q_abs no-recov={Q:.0f} kJ/molCH4; feed preheat 600C: {Qr:.0f}; fuel at eta 0.5/0.7/0.85: {Q/0.5:.0f}/{Q/0.7:.0f}/{Q/0.85:.0f}; w/preheat eta0.7: {Qr/0.7:.0f}")
# Direct genset MW per MMscf/d
mol_s=7300e6/16.04/(365*86400); print('mol/s',mol_s,'LHV MW',mol_s*802.3/1000,'x0.38',mol_s*802.3/1000*0.38)
# 1 MMscf/d of CH4 at 60F,14.696 psia
lbmol=1e6/379.48; t_d=lbmol*16.04*0.45359/1000; print('t CH4/d',t_d,'t/y',t_d*365)
print('MW LHV true 1MMscfd',t_d*365*1e6/16.04*802.3/1e3/31.536e6*1e-3*1e3, )
