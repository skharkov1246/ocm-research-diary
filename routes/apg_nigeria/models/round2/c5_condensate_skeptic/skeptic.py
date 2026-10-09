# Раунд 2, скептик (оптика «второй год и деньги») к оценке c5_condensate.
# Проверки [О]: (1) крупнейшие наземные факелы VIIRS 2025 против порога третьей стороны 19.6 MMscf/д (оптимизм);
# (2) «оператор сам» с налогом оператора (HT 30% + CIT 30% по PIA для PML [В, не проверено по каждой OML];
#     PPT 85% старого режима JV [В]) и с обязательными поправками аудита NaCN (капитал x2.4-2.8, опекс >= скрин x1.7-2.4);
# (3) оптимизм третьей стороны «наполовину»: по одному снимаем оптимистичные допущения (Brent, hp, падение, капитал).
# Запуск: python3 skeptic.py -> skeptic_out.json
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'c5_condensate'))
import econ as E
out = {}
# (1) площадки выше порога
fl = {'Egbema West (NEPL) 2025': 35.8, 'Egbema West min23-25': 1.4, 'Odidi (NEPL) 2025': 20.6, 'Odidi min': 5.0,
      'Obiafu-Obrikom 2025': 18.0, 'Agu (Sterling) 2025': 17.4, 'max single flare 27.6': 27.6}
out['big_flares'] = {}
for k, q in fl.items():
    for gas in ('mid', 'rich'):
        r = E.run(q, gas, 'opt'); rb = E.run(q, gas, 'base')
        out['big_flares'][f'{k}|{gas}'] = dict(opt_npv15=r['npv15'], opt_npv20=r['npv20'], opt_irr=r['irr'], opt_capex=r['capex'],
                                              opt_y2_ebitda=r['y2']['ebitda'], base_npv20=rb['npv20'])
# (1b) сохранение факела: если через 2 года факел вернулся к минимуму 2023-25 (Egbema West 1.4) — opt с падением 50%/г
r = E.run(35.8, 'mid', 'opt', ov=dict(decl=0.5)); out['egbema_collapse_decl50'] = {k: r[k] for k in ('npv15', 'npv20', 'irr')}
# (2) оператор сам: налог и поправки аудита
def op_run(Q, comp, key, sc, tax, audit):
    E.TAX = tax
    if audit:  # оператор = третья сторона по капиталу/опексу, но без охраны и общины (они уже есть), персонал половина
        ov = dict(hp=1.0, comp=comp, jt_key=key, sec=0.0, comm=0.0, staff=E.SC[sc]['staff']/2, ovh=E.SC[sc]['ovh']/2,
                  tariff=0.0, deal=1.0)
        r = E.run(Q, 'mid', sc, ov=ov)
    else:
        r = E.run(Q, 'mid', sc, ov=dict(hp=1.0, comp=comp, jt_key=key), operator=True)
    E.TAX = 0.34
    return r
out['operator'] = {}
for Q in (5, 15, 20.6, 35.8):
    for comp, key in (('SPE_AG6', '30->3'), ('NJTD_A', '50->3')):
        for sc in ('base', 'opt'):
            for tname, tax in (('34', 0.34), ('HT30+CIT30=51', 0.51), ('PPT85', 0.85)):
                for audit in (False, True):
                    r = op_run(Q, comp, key, sc, tax, audit)
                    out['operator'][f'Q{Q}|{comp}_{key}|{sc}|tax{tname}|audit{int(audit)}'] = dict(
                        capex=r['capex'], fixed=r['fixed'], y2E=r['y2']['ebitda'], npv15=r['npv15'], npv20=r['npv20'], irr=r['irr'])
# порог Q для оператора, NPV20=0, оптимизм, с налогом 51% и аудитом
out['op_threshold'] = {}
for comp, key in (('SPE_AG6', '30->3'), ('NJTD_A', '50->3')):
    for tname, tax in (('34', 0.34), ('51', 0.51), ('85', 0.85)):
        for audit in (False, True):
            out['op_threshold'][f'{comp}_{key}|tax{tname}|audit{int(audit)}'] = dict(
                opt=E.solve(lambda x: op_run(x, comp, key, 'opt', tax, audit)['npv20'], 0.3, 400),
                base=E.solve(lambda x: op_run(x, comp, key, 'base', tax, audit)['npv20'], 0.3, 400))
# (3) оптимизм третьей стороны на 15 и 20 MMscf/д: снимаем по одному допущению
out['opt_strip'] = {}
for Q in (15, 20.6, 27.6):
    for lab, ov in (('opt', {}), ('brent70', dict(brent=70)), ('hp0.6', dict(hp=0.6)), ('decl8', dict(decl=0.08)),
                    ('cm2.8', dict(cm=2.8)), ('roy17', dict(roy=0.17)), ('loss8_tariff3', dict(loss=0.08, tariff=3.0)),
                    ('lf0.78', dict(lf=0.78)), ('meoh_hi_900', dict(meoh_hi=True, meoh_usd_t=900)),
                    ('opx2.0_staff0.4_sec0.25', dict(opx_mult=2.0, staff=0.4, sec=0.25)), ('dso120_coll0.85', dict(dso=120, coll=0.85)),
                    ('build2', dict(build=2)), ('precontract/2', dict()),):
        r = E.run(Q, 'mid', 'opt', ov=ov, contract=(lab != 'precontract/2')) if lab != 'precontract/2' else E.run(Q, 'mid', 'opt', ov=dict(deal=0.5), contract=False)
        out['opt_strip'][f'Q{Q}|{lab}'] = dict(y2E=r['y2']['ebitda'], npv15=r['npv15'], npv20=r['npv20'], irr=r['irr'])
json.dump(out, open(os.path.join(HERE, 'skeptic_out.json'), 'w'), ensure_ascii=False, indent=1)
for s, d in out.items():
    print('==', s)
    if isinstance(d, dict) and all(isinstance(v, dict) for v in d.values()):
        for k, v in d.items(): print(k, v)
    else: print(d)
# (4) BOOT-доля: оператор платит нам долю нетбэка (остальное — его маржа за титул/роялти/риск) [Д 0.5-0.8]
out['boot_share'] = {}
for Q in (20.6, 27.6, 35.8):
    for sh in (1.0, 0.8, 0.65, 0.5):
        p = E.SC['opt']; net = (p['brent']-p['diff'])*(1-p['roy'])*(1-p['loss'])-p['tariff']
        # эквивалентная цена Brent, дающая net*sh
        b = (net*sh + p['tariff'])/((1-p['roy'])*(1-p['loss'])) + p['diff']
        r = E.run(Q, 'mid', 'opt', ov=dict(brent=b))
        out['boot_share'][f'Q{Q}|share{sh}'] = dict(net_usd_bbl=round(net*sh, 1), npv15=r['npv15'], npv20=r['npv20'], irr=r['irr'])
json.dump(out, open(os.path.join(HERE, 'skeptic_out.json'), 'w'), ensure_ascii=False, indent=1)
for k, v in out['boot_share'].items(): print(k, v)
