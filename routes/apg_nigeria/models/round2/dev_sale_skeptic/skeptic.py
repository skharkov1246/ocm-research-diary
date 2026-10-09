# Скептик к dev_sale: (1) стратег в рамках обязательных поправок аудита NaCN; (2) найровый PPA без индексации;
# (3) ценность «опциона» продажи на R2 против своего FID; (4) отзыв разрешений 12.2026 до конца ворот 1.
import os, sys, io, contextlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'dev_sale'))
with contextlib.redirect_stdout(io.StringIO()):
    import dev_sale as D
sn = '5_UghelliEast_25km'
p0 = dict(D.SITES[sn], gas=D.gas_from_allin(1.0))
def BE(extra, tl=None, key='npv15', **kw):
    tl = tl or D.TL; kw = dict(D.kw_e(True), **kw)
    return round(D.be(dict(p0, **extra), D.S['dopgas'], tl, key, **kw))
print('=== 1. Безубыточный тариф стратега при поправках аудита (EDTI, стройка 1 г) BE15/BE20')
for cm in (2.2, 2.4, 2.6, 2.8):
    for sm in (1.4, 1.8):
        ex = dict(cm=cm, staff_mult=sm, dev=0.3)
        r = D.val(D.SITES[sn], 405, True, 1, ex)
        print(f'  cm {cm} staff {sm}: BE15 {BE(ex)} BE20 {BE(ex, key="npv20")}  | T405: capex {r["capex"]} NPV15 {r["npv15"]} NPV20 {r["npv20"]} '
              f'премия 4% {0.04*r["capex"]:.2f} / alpha .35 {0.35*r["npv15"]:.2f}')
print('  наш (cm2.6, staff1.8, dev0.8): BE15', BE({}), 'BE20', BE({}, key='npv20'))

print('\n=== 2. PPA в найре: индексация idx (1 = долларовый, 0 = фикс найра), девальвация 8/12/16%')
for dep in (0.08, 0.12, 0.16):
    for idx in (1.0, 0.5, 0.0):
        tl = dict(D.TL, dep=dep, idx=idx)
        for nm, ex in (('стратег x2.2', D.BUY), ('стратег x2.6 экспат', dict(cm=2.6, staff_mult=1.8, dev=0.3))):
            r = D.run(dict(p0, p_dir=405, **ex), D.S['dopgas'], tl, **D.kw_e(True))
            P = max(0, min(0.35*r['npv15'], 0.04*r['capex']))
            print(f'  dep {dep} idx {idx} {nm}: NPV15 {r["npv15"]:6.2f} NPV20 {r["npv20"]:6.2f} премия база {P:.2f}  BE15 {BE(ex, tl)}')

print('\n=== 3. Ценность встроенного опциона «продать на R2 вместо своего FID» (база)')
sc = D.SC['base']
bv = D.buyer[f'{sn}|405|E|1y']
P = D.premium('R2', sc, bv)
pR2 = sc['p0']*sc['p1']*sc['p2']
for pfc_own in (0.65, 0.40, 0.30):
    # если свой FID не состоялся (1-pfc_own), пытаемся продать; корреляция кредитора: стратегу FC даётся с pfc*k
    for corr in (1.0, 0.5):
        sale_net = sc['ps_R2']*((P*sc['upfront']-max(0.03, P*sc['sale_cost']))/1.15**0.3 + sc['pfc']*corr*P*(1-sc['upfront'])/1.15**1.0)
        opt = pR2*(1-pfc_own)*sale_net/1.15**1.0
        print(f'  pfc_own {pfc_own} корреляция {corr}: ценность опциона на t0 = ${opt*1e3:.1f} тыс. (P(R2)={pR2:.3f}, чистая продажа|R2 {sale_net:.3f})')
print('  стоимость опциона: юрвопросы +$3-8 тыс. + зондаж $0-5 тыс. [из оценки]')

print('\n=== 4. Своя стройка при pfc_own ниже, чем у стратега (нет баланса), база ₦405')
o = D.buyer[f'{sn}|405|E|1y']
def gc(r, k=0.5):
    probs = [1.0, sc['p0'], sc['p0']*sc['p1']]
    return sum(pp*(D.stage_cost(g, 0)*(1-k)+D.stage_cost(g, 1)*k)/(1+r)**tt for g, pp, tt in zip(('g0','g1','g2'), probs, D.TIMES['R2']['c']))
for pf in (0.65, 0.4, 0.3, 0.2):
    print(f'  pfc_own {pf}: EV15 своя {pR2*pf*o["npv15_own"]/1.15 - gc(0.15):+.3f}  EV20 {pR2*pf*o["npv20_own"]/1.2 - gc(0.20):+.3f};  продажа R2 EV15 {D.dev_ev("R2", sc, bv, 0.15, 0.5)[0]:+.3f}; стоп после G0: {-D.stage_cost("g0",0)*0.5-D.stage_cost("g0",1)*0.5:+.3f}')

print('\n=== 5. Отзыв разрешений ~12.12.2026 (t=0.18 г): G1 (опцион с держателем) заканчивается t~0.5')
for p_keep in (1.0, 0.6, 0.4):
    s2 = dict(sc, p1=sc['p1']*p_keep)
    print(f'  P(держатель сохраняет разрешение) {p_keep}: R2 EV15 {D.dev_ev("R2", s2, bv, 0.15, 0.5)[0]:+.3f}, R1 EV15 {D.dev_ev("R1", s2, bv, 0.15, 0.5)[0]:+.3f}')

print('\n=== 6. Ughelli East реально 6.5 км до Ughelli, модель 25 км: стратег NPV15 при трассе 8 км')
for km in (8, 25):
    r = D.val(dict(D.SITES[sn], line_km=km, loss=0.03*km/25), 405, True, 1, D.BUY)
    print(f'  трасса {km} км: capex {r["capex"]} NPV15 {r["npv15"]} премия 4% {0.04*r["capex"]:.2f}')

print('\n=== 7. Аукцион без потолка и EV при стратеге в рамках аудита (cm 2.6, экспат), ₦405')
ra = D.val(D.SITES[sn], 405, True, 1, dict(cm=2.6, staff_mult=1.8, dev=0.3))
bva = dict(npv15_buy=ra['npv15'], capex_buy=ra['capex'])
for scn, s in D.SC.items():
    sa = dict(s, cap_R1=1.0, cap_R2=1.0)
    print(f'  {scn}: R2 EV15 с потолком {D.dev_ev("R2", s, bva, 0.15, s["cost_k"])[0]:+.3f}; аукцион P={s["alpha_R2"]*ra["npv15"]:.2f} EV15 {D.dev_ev("R2", sa, bva, 0.15, s["cost_k"])[0]:+.3f}; '
          f'R1 аукцион EV15 {D.dev_ev("R1", sa, bva, 0.15, s["cost_k"])[0]:+.3f}')
