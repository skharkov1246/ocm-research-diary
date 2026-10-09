# Предельный благоприятный стек: all3 + эскроу + капитал по Seplat x2.6 + опекс по верхней границе аудита (скрин x2.4)
# вместо снизу-вверх (x4.7 к скрину). Показывает, держится ли вердикт даже при всех поправках в нашу пользу. [О]
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gas_to_pipe'))
import g2p
TR = 0.80/1.05
lean = dict(staff=0.0, ovh=0.0, sec=0.0, comm=0.0, tor=0.0, ins=0.0)   # bottom-up обнулён -> работает max(.., screen*mult)
out = {}
for site, Q, km in (('Ughelli East', 6.7, 3.7), ('Ughelli East qmin', 6.3, 3.7), ('generic5_4km', 5, 4), ('generic15_4km', 15, 4)):
    for coll_case, cov in (('coll0.85', {}), ('escrow', dict(coll=1.0, fxloss=0.0, dso=45))):
        for b in (None, 0.36, 0.22):
            ov = dict(price=2.18+TR, up=1.0, gas_in=0.25, opx_screen_mult=2.4, **lean, **cov)
            if b:
                pk = g2p.capex(Q, km, dict(g2p.SC['base'], **ov))['pkg']; s = b*Q/pk
                ov.update(usd_hp=1250*s, dehy=0.35*s, meter=0.40*s)
            r = g2p.run(Q, km, 'base', ov=ov, contract=True)
            out[f'{site}|{coll_case}|capex_{b or "bottomup"}'] = dict(cap=r['capex']['total'], fixed=r['fixed_opex'], e2=r['y2']['ebitda'], npv15=r['npv15'], npv20=r['npv20'], irr=r['irr'])
            # то же до контракта (спрос /2)
            r = g2p.run(Q, km, 'base', ov=ov, contract=False)
            out[f'{site}|{coll_case}|capex_{b or "bottomup"}|pre'] = dict(npv15=r['npv15'], irr=r['irr'])
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sk2_out.json'), 'w'), indent=1, ensure_ascii=False)
for k, v in out.items(): print(k, v)
