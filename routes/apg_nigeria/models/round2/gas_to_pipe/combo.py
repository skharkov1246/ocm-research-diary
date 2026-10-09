# Надстройка «газ в трубу» к станции «энергия якорю» на том же факеле (идея критика №5):
# станция берёт 11/(4.40*0.96) = 2.60 MMscf/д [О, исправленная мощность 4.40 МВт на MMscf/д],
# остаток доступного газа (загрузка 0.78 по VIIRS) уходит компрессией в трубу до ТЭС.
# Общие с энергоблоком: площадка, охрана, девелопмент, вход/сепарация. Приростной опекс — только персонал компрессии,
# ТОиР и страховка на приростной капитал. Метки: [О] расчёт, [Д] допущение.
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g2p
HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
for site, Q, km in (('Ughelli East', 6.7, 3.7), ('Oredo (qmin 8.5)', 8.5, 16.4), ('Oredo (2025 15.5)', 15.5, 16.4), ('generic 5', 5.0, 4.0)):
    g_pow = 11/(4.40*0.96)
    for sc in ('opt', 'base', 'pess'):
        lf = g2p.SC[sc]['lf']
        rem_nameplate = max(0.0, Q - g_pow/lf)          # остаток в «номинальных» MMscf/д (до загрузки)
        if rem_nameplate <= 0.05:
            out[f'{site}|{sc}'] = 'остатка нет'; continue
        ov = dict(dev=0.1, staff=0.12, sec=0.0, comm=0.03, ovh=0.05, opx_screen_mult=1.0)
        r = g2p.run(rem_nameplate, km, sc, ov=ov, contract=True)
        out[f'{site}|{sc}'] = dict(rem_nameplate=round(rem_nameplate, 2), capex=r['capex']['total'], fixed=r['fixed_opex'],
                                   ebitda_y2=r['y2']['ebitda'], npv15=r['npv15'], npv20=r['npv20'], irr=r['irr'])
json.dump(out, open(os.path.join(HERE, 'combo_out.json'), 'w'), indent=1, ensure_ascii=False)
for k, v in out.items(): print(k, v)
