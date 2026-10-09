# Законный для ТЭС стек: DBP $2.18 (PIA s.167(5)) вместо $2.68, и вариант DBP + транспортный тариф $0.80 (NERC 2014 [В], статус 2026 НЕ НАЙДЕН)
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'power_plus_pipe'))
import bundle as B
OUT = {}
for nm in ('Ughelli East 2025 (6.7) -> Transcorp Ughelli', 'Oredo 2025 (15.5) -> Sapele PS', 'типовой 15 MMscf/д, отвод 8 км'):
    Q, lk, pk = B.SITES[nm]
    for sc in ('opt', 'base'):
        p0 = B.model(Q, lk, pk, sc, 405, True, 'power')
        for cn, price in (('DBP 2.18', 2.18), ('DBP+транспорт 2.98', 2.98)):
            for cl, pov in (('эскроу 0.98/45', dict(coll=0.98, dso=45)), ('сбор 0.50/300', dict(coll=0.50, dso=300))):
                b = B.model(Q, lk, pk, sc, 405, True, 'bundle', price=price, gas_allin=0.40/1.05, pov=pov)
                k = f'{sc}|{nm[:22]}|{cn}|{cl}'
                OUT[k] = dict(inc15=round(b['npv15']-p0['npv15'], 2), inc20=round(b['npv20']-p0['npv20'], 2))
                print(k, OUT[k])
json.dump(OUT, open(os.path.join(HERE, 'legal_stack_out.json'), 'w'), indent=1, ensure_ascii=False)
