# Скептик к g2p.py: (1) проверка пересчётом; (2) «стальной человек» — снимаем возможные двойные счёты
# в пользу проекта; (3) обратные поправки (падение дебита 2025 -> старт 2028), которых в g2p нет.
# Метки: [О] расчёт, [Д] допущение. Запуск: python3 sk_g2p.py
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gas_to_pipe'))
import g2p
HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
def show(tag, r): 
    out[tag] = dict(cap=r['capex']['total'], fixed=r['fixed_opex'], E1=r['y1']['ebitda'], E2=r['y2']['ebitda'],
                    npv15=r['npv15'], npv20=r['npv20'], irr=r['irr'], Gin2=r['y2']['Gin'])
    print(f"{tag:70s} cap {r['capex']['total']:5.1f} F {r['fixed_opex']:.2f} E2 {r['y2']['ebitda']:6.2f} Gin2 {r['y2']['Gin']:.2f} N15 {r['npv15']:6.1f} IRR {r['irr']}")
SITES = {'UghelliEast': (6.7, 3.7), 'Afam': (5.6, 0.9), 'Oben': (8.7, 2.0), 'Q5_4km': (5, 4), 'Q15_4km': (15, 4)}
# Стальной человек (всё в пользу проекта, не нарушая обязательных поправок NaCN):
#  up=1.0 — VIIRS уже среднегодовой, простой апстрима в нём учтён [Д];
#  ToиР 3.5% (нижняя граница аудита) и страховка 1.0%; транспортная составляющая +$0.5/MMBtu к DBP за прямой отвод [Д, тариф НЕ НАЙДЕН];
#  собираемость 1.0 и DSO 45 без валютных потерь (эскроу), газ на входе $0.25 (пол NGFCP), капитал x2.4.
steel = dict(up=1.0, tor=0.035, ins=0.010, transport=-0.5, coll=1.0, dso=45, fxloss=0.0, gas_in=0.25, cm=2.4)
for s, (Q, km) in SITES.items():
    show(f'{s}|base GSA (как в g2p)', g2p.run(Q, km, 'base', contract=True))
    show(f'{s}|base GSA, up=1.0 (снят двойной счёт простоя)', g2p.run(Q, km, 'base', ov=dict(up=1.0), contract=True))
    show(f'{s}|base GSA, стальной человек', g2p.run(Q, km, 'base', ov=steel, contract=True))
    # обратная поправка: старт эксплуатации 2028 -> дебит 2025 x 0.92^2 (FID после ворот 0 в 12.2026, стройка 1 г) [Д]
    show(f'{s}|base GSA, дебит к старту 2028 (x0.92^2)', g2p.run(Q*0.92**2, km, 'base', ov=dict(), contract=True))
# Порог цены в «стальном человеке» для Ughelli East
pr = g2p.solve(lambda x: g2p.run(6.7, 3.7, 'base', ov=dict(steel, price=x), contract=True)['npv20'], 0.5, 40)
out['steel_price_threshold_UghelliEast_usd_mmbtu'] = pr; print('порог цены (стальной человек) Ughelli East:', pr)
pr = g2p.solve(lambda x: g2p.run(6.7, 3.7, 'base', ov=dict(steel, staff=x/1.3*0+0.0, sec=0, comm=0, ovh=0, opx_screen_mult=1.7, price=2.18+0*x), contract=True)['npv20'], 0, 1)
# Гипотетика: нулевой персонал/охрана/община (только ТОиР+страховка), стальной человек:
r = g2p.run(6.7, 3.7, 'base', ov=dict(steel, staff=0, sec=0, comm=0, ovh=0, opx_screen_mult=0), contract=True)
show('UghelliEast|стальной + без персонала/охраны/общины', r)
# Единичная экономика на номинальный MMscf/д (база, год 2): маржа vs платёж на капитал
m = 0.78*0.92*0.95*0.92*365*1000*(0.935*2.18*1.05*0.85*(1-0.12*120/365) - 0.78)/1e6
c_unit = (g2p.capex(20, 0, g2p.SC['base'])['total'] - g2p.capex(10, 0, g2p.SC['base'])['total'])/10
crf = 0.15/(1-1.15**-10)
out['unit'] = dict(margin_per_nominal=round(m, 3), marginal_capex_per_mmscfd=round(c_unit, 2), capital_charge=round(c_unit*(crf+0.0425+0.015), 3))
print(out['unit'])
json.dump(out, open(os.path.join(HERE, 'sk_g2p_out.json'), 'w'), indent=1, ensure_ascii=False)
# Компрессоры под фактический поток (Q x 0.78), а не под номинал VIIRS x1.5 [Д], поверх стального человека:
for s, (Q, km) in SITES.items():
    show(f'{s}|стальной + компрессия под Q*0.78', g2p.run(Q*0.78, km, 'base', ov=dict(steel, lf=1.0), contract=True))
json.dump(out, open(os.path.join(HERE, 'sk_g2p_out.json'), 'w'), indent=1, ensure_ascii=False)
