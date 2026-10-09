# Скептик (нигерийская оптика) к варианту «газ генераторам покупателя». Импортирует модель оценщика
# round2/gas_to_genset/g2g.py без изменений и пересчитывает спорные места. Метки: [О] расчёт, [Д] допущение, [В] вторичный.
import sys, os, json, copy
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'gas_to_genset'))
import g2g
from g2g import run, solve, D_ENG, D_DUAL, SC, npv, irr, HHV
out = {}

# 1. Поправка аудита «спрос /1.3-2.7 до контракта» обязательна: покупателя нет ни в одной конфигурации.
#    Оценщик в заголовке брал строки GSA (контракт подписан). Пересчёт тех же кейсов в режиме pre.
CASES = {
 'UghelliE->Ughelli 8km eng': (6.3, D_ENG, 1, 6.5*1.3),
 'UghelliE->Warri 27km eng':  (6.3, D_ENG, 1, 20.8*1.3),
 'Oredo->Sapele 28km eng':    (8.5, D_ENG, 1, 21.8*1.3),
 'Utorogu->Ughelli 17km eng': (3.1, min(D_ENG, 3.1*0.78*0.92*0.945), 1, 13.2*1.3),
 'Q5 26km eng':               (5, D_ENG, 1, 26),
 'Q5 26km dual':              (5, D_DUAL, 1, 26),
 'UghelliE->Ughelli 8km dual':(6.3, D_DUAL, 1, 6.5*1.3),
 'Q15 26km cluster_big':      (15, 6.0, 8, 26),
}
r1 = {}
for k, (Q, dm, nb, km) in CASES.items():
    row = {}
    for sc in ('opt', 'base', 'pess'):
        g = run(Q, dm, nb, km, sc, contract=True); p = run(Q, dm, nb, km, sc, contract=False)
        row[sc] = dict(cap=g['design']['total'], npv15_GSA=g['npv15'], npv15_pre=p['npv15'], npv20_pre=p['npv20'],
                       irr_GSA=g['irr'], irr_pre=p['irr'], e2_pre=p['y2']['ebitda'])
    row['p20_base_GSA'] = solve(lambda x: run(Q, dm, nb, km, 'base', ov=dict(price=x))['npv20'], 0.5, 120)
    row['p20_base_pre'] = solve(lambda x: run(Q, dm, nb, km, 'base', ov=dict(price=x), contract=False)['npv20'], 0.5, 200)
    row['p15_base_pre'] = solve(lambda x: run(Q, dm, nb, km, 'base', ov=dict(price=x), contract=False)['npv15'], 0.5, 200)
    row['p20_opt_pre'] = solve(lambda x: run(Q, dm, nb, km, 'opt', ov=dict(price=x), contract=False)['npv20'], 0.5, 200)
    r1[k] = row
out['1_precontract'] = r1

# 2. Честный к оценщику пересчёт: перенос убытков (в g2g налог считается без переноса убытков). Ставим верхнюю оценку
#    эффекта: убираем налог вовсе (tax=0) -- это потолок любого налогового улучшения, включая EDTI.
g2g.TAX = 0.0
r2 = {}
for k in ('UghelliE->Ughelli 8km eng', 'Q5 26km eng', 'Q15 26km cluster_big'):
    Q, dm, nb, km = CASES[k]
    r2[k] = dict(npv15_base_GSA_notax=run(Q, dm, nb, km, 'base')['npv15'],
                 npv15_base_pre_notax=run(Q, dm, nb, km, 'base', contract=False)['npv15'],
                 p20_base_GSA_notax=solve(lambda x: run(Q, dm, nb, km, 'base', ov=dict(price=x))['npv20'], 0.5, 120))
g2g.TAX = 0.34
out['2_tax_ceiling'] = r2

# 3. Нигерийские трения, которых в базе нет: блокада/врезка трубы на 3 мес. во 2-й год (как у судьи для ЛЭП),
#    и сочетание pre + блокада + сбор 85% (неплатежи промпокупателя в найре) для лучшей площадки.
def run_block(Q, dm, nb, km, sc, months=3, year=2, **kw):
    # блокада: в году year продаж нет months/12, фиксированный опекс остаётся. Реализуем через 2 прогона и разность.
    base = run(Q, dm, nb, km, sc, **kw)
    r = dict(base)
    # потеря маржи года: (выручка - переменный газ) * months/12, дисконт на (b + year)
    y = base['y2'] if year == 2 else base['y5']
    b = SC[sc]['build']
    lost = (y['rev'] - y['gas']) * months/12 * (1-0.34)
    r['npv15'] = round(base['npv15'] - lost/(1.15)**(b+year-1), 2)
    r['npv20'] = round(base['npv20'] - lost/(1.20)**(b+year-1), 2)
    return r
Q, dm, nb, km = CASES['UghelliE->Ughelli 8km eng']
out['3_frictions_UghelliE_8km'] = dict(
  block3m_GSA=run_block(Q, dm, nb, km, 'base')['npv15'],
  coll85_GSA=run(Q, dm, nb, km, 'base', ov=dict(coll=0.85))['npv15'],
  pre_coll85_block=run_block(Q, dm, nb, km, 'base', ov=dict(coll=0.85), contract=False)['npv15'],
  opt_pre_coll85_block=run_block(Q, dm, nb, km, 'opt', ov=dict(coll=0.85), contract=False)['npv15'],
  p20_GSA_coll85_block=solve(lambda x: run_block(Q, dm, nb, km, 'base', ov=dict(price=x, coll=0.85))['npv20'], 0.5, 120),
)

# 4. Потолок цены, который покупатель реально видит, против порогов. Покупатель на дизеле с dual-fuel комплектом
#    сравнивает наш газ не с дизелем, а с привозным CNG/СПГ в тот же комплект: $8-10.75 CNG, $4.7-6 СПГ [В фаз 1-2].
#    Трубный газ: регулируемый ориентир commercial $2.68/MMBtu с 01.04.2026 [В: Leadership, NMDPRA] + доставка НЕ НАЙДЕНА.
out['4_ceiling'] = dict(cng=[8.0, 10.75], lng_greenville=[4.7, 6.0], commercial_DBP_2026=2.68,
                        note='порог NPV20=0 любой одиночной конфигурации в базе > потолка CNG')

# 5. Метановое число жирного попутного газа (Kubesh 1992: MON по составу, MN = 1.624*MON - 119.1) [О]
COMP = {}
src = open(os.path.join(HERE, '..', 'c5_condensate', 'process.py')).read()
import re
for name in ('SPE_AG2', 'SPE_AG4', 'SPE_AG5', 'SPE_AG6', 'NJTD_A', 'lean_ND'):
    m = re.search(r"'%s':\s*dict\(([^)]*)\)" % name, src); d = {}
    for kv in m.group(1).split(','):
        if '=' in kv: kk, vv = kv.split('='); d[kk.strip()] = float(vv)
    COMP[name] = d
def mn(d, drop_c5=False):
    t = {k: v for k, v in d.items()}
    if drop_c5:
        for k in ('iC5', 'nC5', 'C6', 'C7', 'C8', 'C9', 'C10'): t[k] = 0.0
        t['iC4'] = t.get('iC4', 0)*0.85; t['nC4'] = t.get('nC4', 0)*0.7
    s = sum(t.values()); x = {k: v/s for k, v in t.items()}
    c4p = sum(x.get(k, 0) for k in ('iC4', 'nC4', 'iC5', 'nC5', 'C6', 'C7', 'C8', 'C9', 'C10'))
    mon = 137.78*x.get('C1', 0) + 29.948*x.get('C2', 0) - 18.193*x.get('C3', 0) - 167.062*c4p + 181.233*x.get('CO2', 0) + 26.994*x.get('N2', 0)
    return round(1.624*mon - 119.1, 1)
out['5_methane_number'] = {k: dict(raw=mn(v), after_dewpoint=mn(v, True)) for k, v in COMP.items()}

json.dump(out, open(os.path.join(HERE, 'sk_out.json'), 'w'), indent=1, ensure_ascii=False)
for k, v in out.items():
    print('==', k)
    if isinstance(v, dict):
        for kk, vv in v.items(): print(' ', kk, vv)
