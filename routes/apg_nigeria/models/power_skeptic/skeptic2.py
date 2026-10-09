import io,contextlib,sys
sys.path.insert(0,'.')
with contextlib.redirect_stdout(io.StringIO()):
    import skeptic as K
C=K.C; b=C.S['base']
# 1) buyer BATNA: own gas engine on trucked LNG/CNG
FX=1330
for usd in (4.90,6.4,9.8):
    for capcrf,hrs in ((1500*0.171,5256),(1500*0.199,5256)):
        fuel=usd*8.98; om=(15,25); cap=capcrf/hrs*1000
        lo=(fuel+om[0]+cap)*FX/1000; hi=(fuel+om[1]+cap+5)*FX/1000
        print(f'LNG/CNG ${usd}/MMBtu CRF {capcrf/1500:.3f}: N{lo:.0f}-{hi:.0f}/kWh')
print('Greenville N230/scm ->',230/0.0353/FX,'$/MMBtu; CNG truck N450 ->',450/0.0353/FX)
# 2) price scan with combined soft corrections (as prior skeptic) + generator-net price
soft=dict(lag=2,extra_opex=0.3,top=True,extra_capex=1.0,wc_days=60)
for pr in (200,210,220,250,280,320):
    p=dict(b);p['price_ngn']=pr
    print('net price',pr,'clean',K.run(p),' soft',K.run(p,**soft))
# tripartite via PHED: net price = BATNA ceiling - DUoS/margin, DisCo bad debt 15%
for ceil,duos in ((250,30),(250,50),(280,30)):
    p=dict(b);p['price_ngn']=ceil-duos;p['baddebt']=0.15
    print('tripartite ceil',ceil,'duos',duos,K.run(p),K.run(p,**soft))
# FX: 10%/y dep, 50% indexation on top of soft
print('soft+fx10 idx50',K.run(b,fx_dep=0.10,index=0.5,**soft))
# gas price claim
for g in (0.25,2.68):
    p=dict(b);p['gas']=g;print('gas',g,K.run(p))
# break-even price under soft for NPV15/20=0 ~ from prior: 360 / 430
# capital efficiency vs mining
print('EBITDA/capex gen base',4.14/14.61,'soft',3.73/15.61,'mining $50',3.5/9.02,'$27.7',1.94/9.02)
