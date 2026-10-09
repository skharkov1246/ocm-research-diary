# Скептик к девелоперской модели (dev_sale.py). Запуск: python3 skeptic_dev_sale.py > skeptic_dev_sale_out.txt
# Все результаты [О]. Проверяем: (1) стратег с капиталом x2.2 нарушает обязательную поправку аудита x2.4-2.8;
# (2) «LOI на ₦45-55 ниже» убивает запасной путь своей стройки; (3) ценность продажи как встроенного опциона;
# (4) спрос /1.3-2.7 при стратеге по аудиту.
import os, sys, io, contextlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'dev_sale'))
with contextlib.redirect_stdout(io.StringIO()):
    import dev_sale as D
S = D.S; TL = D.TL
SN = '5_UghelliEast_25km'
def v(T, ex, site=SN, d=1, **over):
    return D.val(dict(D.SITES[site], **over), T, True, d, ex)
def be(ex, key, site=SN):
    p = dict(D.SITES[site], gas=D.gas_from_allin(1.0), **ex)
    return round(D.be(p, S['dopgas'], TL, key, **D.kw_e(True)))

print('=== 1. Стратег по поправке аудита (капитал x2.4-2.8; персонал 1.4 без экспата / 1.8 с экспатом), ₦405, EDTI, 1 г')
for site in (SN, '15_Oredo_30km_Q8.5'):
    for cm in (2.2, 2.4, 2.6, 2.8):
        for sm in (1.4, 1.8):
            ex = dict(cm=cm, staff_mult=sm, dev=0.3)
            b = v(405, ex, site)
            P2 = min(0.35*b['npv15'], 0.04*b['capex']) if b['npv15'] > 0 else 0
            print(f'  {site} cm {cm} staff {sm}: capex {b["capex"]}, NPV15 {b["npv15"]}, NPV20 {b["npv20"]}, IRR {b["irr"]}; '
                  f'BE15/BE20 {be(ex,"npv15",site)}/{be(ex,"npv20",site)}; BE20 с премией 1.0 {be(dict(ex, dev=1.3),"npv20",site)}; премия R2 база {P2:.2f}')
print('  наш (cm 2.6, staff 1.8, dev 0.8): BE15/BE20', be(D.OWN, 'npv15'), '/', be(D.OWN, 'npv20'))

print('\n=== 2. Дев. EV (база/пессимизм/оптимизм) при стратеге по аудиту, ₦405')
for cm, sm in ((2.2, 1.4), (2.4, 1.4), (2.6, 1.8)):
    b = v(405, dict(cm=cm, staff_mult=sm, dev=0.3))
    bv = dict(npv15_buy=b['npv15'], capex_buy=b['capex'])
    row = []
    for scn in ('base', 'pess', 'opt'):
        sc = D.SC[scn]
        e15 = D.dev_ev('R2', sc, bv, 0.15, sc['cost_k'])[0]
        ea = D.dev_ev('R2', dict(sc, cap_R2=1.0), bv, 0.15, sc['cost_k'])[0]
        row.append(f'{scn}: EV15 {e15:+.3f}, аукцион {ea:+.3f} (P {sc["alpha_R2"]*b["npv15"]:.2f})')
    print(f'  cm {cm} staff {sm}:', '; '.join(row))

print('\n=== 3. «LOI на ₦45-55 ниже»: что будет с нашим запасным путём своей стройки')
for T in (300, 345, 350, 377, 405):
    o = v(T, D.OWN); b1 = v(T, D.BUY); b2 = v(T, dict(cm=2.4, staff_mult=1.4, dev=0.3)); b3 = v(T, dict(cm=2.6, staff_mult=1.8, dev=0.3))
    print(f'  ₦{T}: наша станция NPV15 {o["npv15"]:+.2f} NPV20 {o["npv20"]:+.2f}; стратег x2.2 NPV20 {b1["npv20"]:+.2f}; '
          f'стратег x2.4 NPV20 {b2["npv20"]:+.2f}; стратег как мы NPV20 {b3["npv20"]:+.2f}')

print('\n=== 4. Ценность продажи как встроенного опциона в воротах судьи (база, ₦405, R2 достигнут)')
# В состоянии R2: решаем «строим сами» (FC с вероятностью pfc_own) или «продаём» (стратег покупает ps, закрывает pfc).
# Рационально: строим сами, если FC достижимо; продажа ценна только в ветке, где наш FC не состоялся.
# Опцион = P(R2) * (1-pfc_own) * ps * pfc_buy|our_fail * чистая премия PV.  Корреляция: если мы не закрылись из-за
# слабости проекта (якорь/PPA), стратег тоже закрывается реже -> множитель rho [Д 0.3-1.0].
sc = D.SC['base']
o = v(405, D.OWN); b = v(405, D.BUY)
P = min(sc['alpha_R2']*b['npv15'], sc['cap_R2']*b['capex'])
net_sale = (P*sc['upfront']*(1-sc['sale_cost'])/1.15**0.3 + P*(1-sc['upfront'])/1.15**1.0)*0.85   # налог/расходы ~15% [Д]
pR2 = sc['p0']*sc['p1']*sc['p2']
for pfc_own in (0.65, 0.35, 0.20):
    for rho in (1.0, 0.5, 0.3):
        build = pfc_own*o['npv15']
        opt_val = (1-pfc_own)*sc['ps_R2']*sc['pfc']*rho*net_sale
        # только продажа (как в оценке): ps*pfc*... vs своя стройка
        print(f'  pfc_own {pfc_own} rho {rho}: в R2 своя стройка {build:.2f}, продажа-как-запаска добавляет {opt_val:.3f} '
              f'(в EV на t=0: {pR2*opt_val/1.15:.4f}) ; продажа вместо стройки {sc["ps_R2"]*net_sale:.2f}')
print('  премия R2', round(P, 2), 'чистыми PV ~', round(net_sale, 2), 'P(R2)=', pR2)
# если LOI взяли на ₦345 (под стратега), своя стройка отрицательна по NPV20 -> запаски нет, EV пути = только продажа
o345 = v(345, D.OWN); b345 = v(345, D.BUY)
P345 = min(sc['alpha_R2']*b345['npv15'], sc['cap_R2']*b345['capex'])
print(f'  LOI ₦345: наша NPV15 {o345["npv15"]:+.2f} NPV20 {o345["npv20"]:+.2f}; стратег x2.2 NPV15 {b345["npv15"]:+.2f}; премия {P345:.2f}; '
      f'EV15 продажи {D.dev_ev("R2", sc, dict(npv15_buy=b345["npv15"], capex_buy=b345["capex"]), 0.15, 0.5)[0]:+.3f}')
b345a = v(345, dict(cm=2.4, staff_mult=1.4, dev=0.3)); b345c = v(345, dict(cm=2.6, staff_mult=1.8, dev=0.3))
print(f'  LOI ₦345, стратег x2.4: NPV15 {b345a["npv15"]:+.2f} NPV20 {b345a["npv20"]:+.2f}; стратег как мы: NPV15 {b345c["npv15"]:+.2f} NPV20 {b345c["npv20"]:+.2f}')

print('\n=== 5. Спрос /1.3 и /2.7 при стратеге по аудиту (x2.4, 1.4), ₦405')
for div in (1.0, 1.3, 2.7):
    b = v(405, dict(cm=2.4, staff_mult=1.4, dev=0.3), anchor=11.0/div)
    print(f'  LOI/{div}: {11/div:.1f} МВт, capex {b["capex"]}, NPV15 {b["npv15"]}, NPV20 {b["npv20"]}')

print('\n=== 6. Сколько P(готово и продано) нужно при премии в потолке 4% и при аукционе, стратег x2.4')
b = v(405, dict(cm=2.4, staff_mult=1.4, dev=0.3)); bv = dict(npv15_buy=b['npv15'], capex_buy=b['capex'])
def solve(f, lo, hi):
    if f(lo)*f(hi) > 0: return None
    for _ in range(60):
        m = (lo+hi)/2
        if f(lo)*f(m) <= 0: hi = m
        else: lo = m
    return round((lo+hi)/2, 3)
for lab, s2 in (('потолок 4%', dict(sc)), ('аукцион 0.5xNPV15', dict(sc, alpha_R2=0.5, cap_R2=1.0))):
    f = lambda x: D.dev_ev('R2', dict(s2, p0=x, p1=1, p2=1, ps_R2=1), bv, 0.15, 0.5)[0]
    print(f'  {lab}: нужна P(готово и продано) >= {solve(f, 0.0005, 1.0)}; модель даёт {sc["p0"]*sc["p1"]*sc["p2"]*sc["ps_R2"]:.3f}')
