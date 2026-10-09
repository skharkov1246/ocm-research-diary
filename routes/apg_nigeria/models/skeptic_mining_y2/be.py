import sys,io,contextlib,json
sys.path.insert(0,'/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/compute')
sys.path.insert(0,'.')
with contextlib.redirect_stdout(io.StringIO()):
    import mining_ai as M
    import hyb_y2 as H
def bisect(f):
    lo,hi=1,1000
    for _ in range(60):
        m=(lo+hi)/2
        if f(m)<0: lo=m
        else: hi=m
    return round(m,1)
res={}
# гибрид: безубыточный hashprice (NPV15=0), база
for prof in ['flat','two_shift','one_shift_plus_base']:
    for decl in [0.0,0.05]:
        for hm in [10,0]:
            res[f'hyb_BE|{prof}|decl{decl}|halv_m{hm}']=bisect(lambda hp:H.run(hp,prof,decl=decl,halv_month=hm)['npv15'])
# гибрид оптимизм-параметры на двухсменном профиле
opt=dict(usd_th=9,landed=1.15,gas=0.25,vom=12,erosion=0.0,halv=0.70,pool=0.025,tor=0.035)
for prof in ['flat','two_shift']:
    for hp in [39.87,50]:
        res[f'hyb_opt|{prof}|hp{hp}']=H.run(hp,prof,S=1.4,**opt)['npv15']
# автономный майнинг, 4 г, старт ПОСЛЕ халвинга (задержка 12-36 мес. из контекста)
for hm in [10,0]:
    M.HALVING_MONTH=hm
    for sc in ['base','opt']:
        r={hp:M.mining(hp,'A_13.5JTH',sc) for hp in [27.7,39.87,50.0]}
        res[f'solo4y_{sc}_halv_m{hm}']={hp:(v['rows'][0]['ebitda'],v['rows'][1]['ebitda'],v['npv15']) for hp,v in r.items()}
        res[f'solo4y_{sc}_BE15_halv_m{hm}']=bisect(lambda hp:M.mining(hp,'A_13.5JTH',sc)['npv15'])
M.HALVING_MONTH=10
# год 2, база, $39.87: DSCR при 50% долга на 5 лет под 14%
r=M.mining(39.87,'A_13.5JTH','base'); debt=0.5*r['capex']; ann=debt*0.14/(1-1.14**-5)
res['dscr_base_39.87']=[round(x['ebitda']/ann,2) for x in r['rows']]; res['debt_service']=round(ann,2)
# роль 2: ASIC на всей станции, пока нет клиента (BOP утоплен), пост-халвинг, маржинально
for hp in [27.7,39.87,50]:
    th=4.14e6/13.5; capex=th*13*1.2/1e6
    mwh=4.14*8760*0.88; rev=th/1000*hp*0.6*365*0.88/1e6
    eb=rev-mwh*22.5/1e6-0.04*rev-0.04*capex
    res[f'role2_marg_hp{hp}']=dict(capex=round(capex,2),ebitda_y=round(eb,2),simple_pb=round(capex/eb,1) if eb>0 else None)
json.dump(res,open('be_out.json','w'),indent=1)
for k,v in res.items(): print(k,v)
