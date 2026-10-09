import sys,io,contextlib
sys.path.insert(0,'/tmp/claude-0/-home-user-ocm-research-diary/50ead64b-47e2-5bdd-b34b-09f6667b6ca9/scratchpad/multi')
with contextlib.redirect_stdout(io.StringIO()):
    import multi_model as M
def s(Q,cfg,**ov):
    r=M.site(Q,'base',cfg,ov=ov); return {k:r[k] for k in ('capex','npv15','npv20','irr','payback')}|{'e1':r['y1']['ebitda'],'MW':r['MW_gross']}
P=('P',)
print('ref 5 P',s(5,P)); print('ref 15 P',s(15,P))
for Q in (5,15):
  for pr in (320,280,250,230,210):
    print(Q,'tariff',pr,s(Q,P,p_dir=pr),' +delay1',s(Q,P,p_dir=pr,delay=1))
# gas 1/4 of declared (NGFCP round-1 lesson)
for Q in (5,15):
  for av in (0.6,0.4,0.25):
    print(Q,'avail',av,s(Q,P,avail=av))
# combined realistic: tariff 250, delay 1, cm 2.6
for Q in (5,15):
  print(Q,'COMB 250+delay1',s(Q,P,p_dir=250,delay=1))
  print(Q,'COMB 250+delay1+FX(1330->1600 unindexed)',s(Q,P,p_dir=250*1330/1600,delay=1))
  print(Q,'COMB 280+delay1+LF0.6',s(Q,P,p_dir=280,delay=1,LF=0.6))
# expected value power-only 5: 50/50 base/pess
b=M.site(5,'base',P); p=M.site(5,'pess',P)
print('EV 5 P npv15',0.5*b['npv15']+0.5*p['npv15'],'npv20',0.5*b['npv20']+0.5*p['npv20'])
# take-or-pay on whole available flare (GSA volume = site), power only
def top_extra(Q):
    r=M.site(Q,'base',P); G=Q*0.85; ex=(G-r['g_pow'])*1000*365*(0.6*1.05+0.15)/1e6
    return round(ex,2)
print('ToP on unused flare gas $M/y: 5',top_extra(5),'15',top_extra(15))
# share of flare actually extinguished by power-only
for Q in (5,15):
    r=M.site(Q,'base',P); print(Q,'flare utilised share',round(r['g_pow']/(Q*0.85),2), 'with LF/av burn', round(r['g_pow']*0.72*0.92/(Q*0.85),2))
# full line scenarios EV
for Q in (5,15):
    b=M.site(Q,'base',('P','L','C')); p=M.site(Q,'pess',('P','L','C'))
    print(Q,'full line EV npv15',0.5*b['npv15']+0.5*p['npv15'])
print('screen 5 full',M.site(5,'screen',('P','L','C'))['y1']['ebitda'], [M.site(5,'base',('P','L','C'),ov=dict(cm=c))['y1']['ebitda'] for c in (2.4,2.8)], [M.site(5,'base',('P','L','C'),ov=dict(cm=c))['capex'] for c in (2.4,2.8)])
