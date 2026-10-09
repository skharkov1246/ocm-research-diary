# Peng-Robinson isothermal flash, pure python. Checks C3 recovery claim for refrigeration LTS.
import math
R=8.314
# Tc K, Pc bar, omega
P={'N2':(126.2,33.9,0.039),'CO2':(304.1,73.8,0.225),'C1':(190.6,46.0,0.011),'C2':(305.3,48.7,0.099),
   'C3':(369.8,42.5,0.152),'iC4':(408.1,36.5,0.181),'nC4':(425.1,38.0,0.200),'iC5':(460.4,33.8,0.227),
   'nC5':(469.7,33.7,0.251),'C6':(507.6,30.2,0.301),'C7':(540.2,27.4,0.349)}
def kij(i,j):
    s={i,j}
    if 'CO2' in s and len(s)==2 and s!={'CO2','N2'}: return 0.12
    if 'N2' in s and len(s)==2 and s!={'CO2','N2'}: return 0.03 if 'C1' in s else 0.08
    return 0.0
def ab(T,c):
    Tc,Pc,w=P[c]; Pc*=1e5
    m=0.37464+1.54226*w-0.26992*w*w
    a=0.45724*R**2*Tc**2/Pc*(1+m*(1-math.sqrt(T/Tc)))**2
    b=0.07780*R*Tc/Pc
    return a,b
def cubic_roots(c2,c1,c0):
    # z^3+c2 z^2+c1 z+c0
    q=(3*c1-c2*c2)/9; r=(9*c2*c1-27*c0-2*c2**3)/54
    D=q**3+r*r
    if D>=0:
        s=math.copysign(abs(r+math.sqrt(D))**(1/3),r+math.sqrt(D)); t=math.copysign(abs(r-math.sqrt(D))**(1/3),r-math.sqrt(D))
        return [s+t-c2/3]
    th=math.acos(r/math.sqrt(-q**3))
    return [2*math.sqrt(-q)*math.cos((th+2*k*math.pi)/3)-c2/3 for k in range(3)]
def lnphi(x,T,p,phase):
    comps=list(x); A={};B={}
    for c in comps: A[c],B[c]=ab(T,c)
    am=sum(x[i]*x[j]*math.sqrt(A[i]*A[j])*(1-kij(i,j)) for i in comps for j in comps)
    bm=sum(x[i]*B[i] for i in comps)
    Am=am*p/(R*T)**2; Bm=bm*p/(R*T)
    roots=[z for z in cubic_roots(-(1-Bm),Am-3*Bm*Bm-2*Bm,-(Am*Bm-Bm*Bm-Bm**3)) if z>Bm]
    z=max(roots) if phase=='v' else min(roots)
    out={}
    for i in comps:
        s=sum(x[j]*math.sqrt(A[i]*A[j])*(1-kij(i,j)) for j in comps)
        out[i]=B[i]/bm*(z-1)-math.log(z-Bm)-Am/(2*math.sqrt(2)*Bm)*(2*s/am-B[i]/bm)*math.log((z+(1+math.sqrt(2))*Bm)/(z+(1-math.sqrt(2))*Bm))
    return out
def flash(zf,T,p):
    K={c:P[c][1]*1e5/p*math.exp(5.373*(1+P[c][2])*(1-P[c][0]/T)) for c in zf}
    for it in range(500):
        lo,hi=0.0,1.0
        f=lambda V: sum(zf[c]*(K[c]-1)/(1+V*(K[c]-1)) for c in zf)
        if f(0)<0: return 0.0,{c:zf[c] for c in zf},None
        if f(1)>0: return 1.0,None,zf
        for _ in range(100):
            m=(lo+hi)/2
            if f(m)>0: lo=m
            else: hi=m
        V=m
        x={c:zf[c]/(1+V*(K[c]-1)) for c in zf}; y={c:K[c]*x[c] for c in zf}
        sx=sum(x.values()); sy=sum(y.values()); x={c:x[c]/sx for c in x}; y={c:y[c]/sy for c in y}
        pl=lnphi(x,T,p,'l'); pv=lnphi(y,T,p,'v')
        Kn={c:math.exp(pl[c]-pv[c]) for c in zf}
        err=sum(abs(Kn[c]/K[c]-1) for c in zf); K=Kn
        if err<1e-9: break
    return V,x,y
rich=dict(N2=0.10,CO2=0.46,C1=82.72,C2=6.73,C3=5.68,iC4=1.56,nC4=1.60,iC5=0.40,nC5=0.26,C6=0.16,C7=0.20)
mid=dict(C1=88.0,C2=4.6,C3=2.3,iC4=0.6,nC4=0.6,iC5=0.35,nC5=0.35)
for name,zz in [('rich',rich),('mid',mid)]:
    s=sum(zz.values()); z={c:v/s for c,v in zz.items()}
    for Tc_ in (-20,-30,-35,-40):
        for pb in (27.5,40):
            V,x,y=flash(z,273.15+Tc_,pb*1e5)
            L=1-V
            rec={c:round(L*x[c]/z[c],2) for c in ('C1','C2','C3','iC4','nC4','nC5')} if x else None
            print(name,Tc_,'C',pb,'bar','L frac',round(L,4),'liq recov',rec)
