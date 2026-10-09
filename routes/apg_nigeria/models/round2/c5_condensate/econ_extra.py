# Раунд 2. Дополнения к econ.py: потолок выручки при 100% извлечении C5+, оператор сам с глубоким JT,
# пороги «оператор сам». Запуск: python3 econ_extra.py -> econ_extra_out.json
import json, os
import econ as E
HERE = os.path.dirname(os.path.abspath(__file__))
out = dict(ceiling={}, operator_deep={}, op_thresholds={}, us_best={})
# 1) Потолок: весь C5+ пробы, загрузка 0.78, без вычетов и без капитала; валовая выручка при Brent $70 [О]
for comp in ('SPE_AG1', 'SPE_AG2', 'SPE_AG4', 'SPE_AG6', 'NJTD_A', 'lean_ND'):
    pot = E.PROC['potential'][comp]['C5p_bbl_per_MMscf']
    for Q in (1, 5, 15):
        bbl_d = Q*0.78*pot
        out['ceiling'][f'{comp}_Q{Q}'] = dict(C5p_bbl_per_MMscf=pot, bbl_d=round(bbl_d, 1),
            gross_usd_m_70=round(bbl_d*365*70/1e6, 3), net_usd_m_base=round(bbl_d*365*48.16/1e6, 3))
# 2) Оператор делает сам, газ ВД 100%, глубокий JT с метанолом (верхняя граница выхода)
for Q in (1, 5, 15):
    for gas, comp, key in (('mid', 'SPE_AG2', '30->3'), ('mid', 'SPE_AG6', '30->3'), ('rich', 'SPE_AG4', '30->3'), ('rich', 'NJTD_A', '50->3')):
        for sc in ('base', 'opt'):
            r = E.run(Q, gas, sc, ov=dict(hp=1.0, comp=comp, jt_key=key), operator=True)
            out['operator_deep'][f'Q{Q}_{comp}_{key}_{sc}'] = {k: r[k] for k in ('yield_bbl_per_MMscf', 'meoh_kg_per_MMscf', 'net_usd_bbl', 'capex', 'fixed', 'y1', 'y2', 'npv15', 'npv20', 'irr', 'payback')}
for comp, key in (('SPE_AG2', 'noinh:30'), ('SPE_AG2', '30->3'), ('SPE_AG6', '30->3'), ('NJTD_A', '50->3')):
    out['op_thresholds'][f'{comp}_{key}'] = dict(
        Q_npv20_base=E.solve(lambda x: E.run(x, 'mid', 'base', ov=dict(hp=1.0, comp=comp, jt_key=key), operator=True)['npv20'], 0.3, 400),
        Q_npv20_opt=E.solve(lambda x: E.run(x, 'mid', 'opt', ov=dict(hp=1.0, comp=comp, jt_key=key), operator=True)['npv20'], 0.3, 400),
        Q_npv20_us_base=E.solve(lambda x: E.run(x, 'mid', 'base', ov=dict(hp=1.0, comp=comp, jt_key=key))['npv20'], 0.3, 400))
# 3) Мы (третья сторона), лучшее сочетание на базовых ценах: газ ВД 100%, глубокий JT, богатая проба
for Q in (5, 15):
    for comp, key in (('SPE_AG6', '30->3'), ('NJTD_A', '50->3')):
        r = E.run(Q, 'mid', 'base', ov=dict(hp=1.0, comp=comp, jt_key=key))
        out['us_best'][f'Q{Q}_{comp}_{key}_base'] = {k: r[k] for k in ('yield_bbl_per_MMscf', 'capex', 'fixed', 'y1', 'y2', 'npv15', 'npv20', 'irr')}
json.dump(out, open(os.path.join(HERE, 'econ_extra_out.json'), 'w'), ensure_ascii=False, indent=1)
for s, d in out.items():
    print('==', s)
    for k, v in d.items(): print(k, v)
