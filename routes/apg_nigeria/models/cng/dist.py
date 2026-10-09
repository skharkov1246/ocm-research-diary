import csv,math,sys
f=sys.argv[1]
rows=[r for r in csv.DictReader(open(f)) if r['Year']=='2025' and float(r['BCM'] or 0)>0]
# CNG demand hubs (industrial / transport), coordinates approx [Д]
hubs={'PortHarcourt':(4.815,7.035),'Onne':(4.71,7.15),'Aba':(5.12,7.37),'Owerri':(5.48,7.03),'Uyo':(5.03,7.93),
'Warri':(5.52,5.75),'Sapele':(5.89,5.68),'Benin':(6.34,5.63),'Asaba/Onitsha':(6.16,6.77),'Nnewi':(6.02,6.92),'Enugu':(6.45,7.51),'Lagos':(6.6,3.35)}
def hav(a,b,c,e):
    R=6371;p=math.radians;dl=p(c-a);dn=p(e-b)
    h=math.sin(dl/2)**2+math.cos(p(a))*math.cos(p(c))*math.sin(dn/2)**2
    return 2*R*math.asin(math.sqrt(h))
ROAD=1.35  # коэффициент извилистости дорог [Д 1.25-1.5]
res=[]
for r in rows:
    la,lo=float(r['lat']),float(r['lon']);q=float(r['MMSCFD'])
    d={h:hav(la,lo,*c)*ROAD for h,c in hubs.items()}
    res.append((r,q,d))
on=[x for x in res if x[0]['Location']=='ONSHORE']
print('onshore',len(on),round(sum(x[1] for x in on),1))
for thr in [1,5]:
  for R in [50,100,150,200]:
    for hubset,name in [(['PortHarcourt','Onne','Aba','Owerri','Uyo','Warri','Sapele','Benin','Asaba/Onitsha','Nnewi','Enugu'],'любой хаб Юга'),(['Asaba/Onitsha','Nnewi','Enugu','Benin'],'Онича/Энугу/Бенин')]:
      y=[x for x in on if x[1]>=thr and min(x[2][h] for h in hubset)<=R]
      print(f'>= {thr} MMscfd, дорога <= {R} км до [{name}]: n={len(y)}, vol={sum(x[1] for x in y):.1f}')
# median road km from >=1 onshore flares to nearest hub and to Onitsha, Lagos
import statistics as st
x=[t for t in on if t[1]>=1]
print('n>=1',len(x),'median km to nearest south hub',round(st.median(min(v for k,v in t[2].items() if k!='Lagos') for t in x)))
print('median km to Onitsha',round(st.median(t[2]['Asaba/Onitsha'] for t in x)),'to Lagos',round(st.median(t[2]['Lagos'] for t in x)))
# by size class
for lo_,hi_ in [(1,5),(5,10),(10,100)]:
    z=[t for t in on if lo_<=t[1]<hi_]
    print(f'class {lo_}-{hi_}: n={len(z)} vol={sum(t[1] for t in z):.0f}')
