# Скептик (нигерийская оптика) к «Конденсат C5+ с JT-сепарации». Метки: [О] расчёт, [Д] допущение, [В] вторичный.
# Повторно использует econ.py оценщика без изменений; меняет только входы.
import json, os, sys, importlib.util
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', 'c5_condensate')
sys.path.insert(0, SRC); os.chdir(SRC)
import econ as E
out = {}
# 0) Воспроизведение ключевых чисел оценщика
for Q in (1, 5, 15):
    r = E.run(Q, 'mid', 'base'); o = E.run(Q, 'mid', 'opt')
    out[f'repro_Q{Q}'] = dict(capex=r['capex'], y2_bbl_d=r['y2']['bbl_d'], ebitda_y2=r['y2']['ebitda'], npv15=r['npv15'], npv15_opt=o['npv15'], irr_opt=o['irr'])
# 1) Рациональная база: метанол дешевле упущенного выхода (ингибитор vs без ингибитора), AG2, 30->3
y_no = E.NOINH['SPE_AG2']['30']['yield_stab_bbl_per_MMscf']; y_me = E.jt_yield('SPE_AG2', '30->3')
mk = E.meoh_kg('SPE_AG2', '30->3', hi=True)
out['meoh_vs_yield_per_MMscf'] = dict(extra_bbl=round(y_me-y_no, 2), extra_rev_usd=round((y_me-y_no)*48.16, 1),
                                     meoh_cost_usd_hi=round(mk*0.9, 1))   # $900/т, верх расхода
# 2) Лучшая мыслимая третья сторона: O&M и охрана оператора, газ ВД 100%, глубокий JT AG6/NJTD, капитал x2.4 (низ аудита)
for comp, key in (('SPE_AG6', '30->3'), ('NJTD_A', '50->3')):
    for Q in (15, 27.6):
        r = E.run(Q, 'mid', 'base', ov=dict(hp=1.0, comp=comp, jt_key=key, staff=0.10, sec=0.05, comm=0.05, ovh=0.05, cm=2.4, dev=0.3, opx_mult=1.7))
        out[f'lean3p_Q{Q}_{comp}_{key}'] = {k: r[k] for k in ('capex', 'fixed', 'y2', 'npv15', 'npv20', 'irr')}
# 3) Нигерийская реальность на единственных площадках с NPV>=0 в оптимизме (Obiafu-Obrikom 18.0, Oredo 15.5):
#    потери Ebocha-Brass/TNP/TFP 2021-22 (десятки %) [В], дебиторка NNPC-JV 180 сут [Д], доля нетбэка нам 60% [Д],
#    простои трубы (форс-мажоры) готовность 0.85 [Д]
for name, Q in (('Obiafu-Obrikom', 18.0), ('Oredo', 15.5)):
    base_opt = E.run(Q, 'mid', 'opt')
    res = dict(opt=base_opt['npv15'])
    for lab, ov in (('loss20', dict(loss=0.20)), ('loss30', dict(loss=0.30)), ('dso180', dict(dso=180)), ('up0.85', dict(up=0.85)),
                    ('fee60pct', None), ('all', None)):
        if lab == 'fee60pct':
            p = E.SC['opt']; net = (p['brent']-p['diff'])*(1-p['roy'])*(1-p['loss'])-p['tariff']
            ov = dict(tariff=p['tariff'] + 0.4*net)
        if lab == 'all':
            p = E.SC['opt']; net = (p['brent']-p['diff'])*(1-p['roy'])*(1-0.20)-p['tariff']
            ov = dict(loss=0.20, dso=180, up=0.85, tariff=p['tariff'] + 0.4*net)
        res[lab] = E.run(Q, 'mid', 'opt', ov=ov)['npv15']
    out[f'ng_reality_opt_{name}'] = res
# 4) Оператор сам: налог оператора не 34%, а HT 30% + CIT 30% (HT вычитается) = 51% (PML на суше) [В PIA s.260/267; Д]
for tax in (0.34, 0.51):
    E.TAX = tax
    out[f'op_thr_tax{tax}'] = {f'{c}_{k}': dict(base=E.solve(lambda x: E.run(x, 'mid', 'base', ov=dict(hp=1.0, comp=c, jt_key=k), operator=True)['npv20'], 0.3, 400),
                                                opt=E.solve(lambda x: E.run(x, 'mid', 'opt', ov=dict(hp=1.0, comp=c, jt_key=k), operator=True)['npv20'], 0.3, 400))
                               for c, k in (('SPE_AG6', '30->3'), ('NJTD_A', '50->3'))}
E.TAX = 0.34
# 5) Потолок с C4: если весь C4+ закачать в нефть (RVP пула позволяет при малых объёмах) [О]
P = E.PROC['potential']
out['c4_note'] = {c: dict(C5p=P[c]['C5p_bbl_per_MMscf'], C3p=P[c]['C3p_bbl_per_MMscf']) for c in ('SPE_AG2', 'SPE_AG4', 'SPE_AG6', 'NJTD_A')}
json.dump(out, open(os.path.join(HERE, 'skeptic_ng_out.json'), 'w'), ensure_ascii=False, indent=1)
for k, v in out.items(): print(k, v)
