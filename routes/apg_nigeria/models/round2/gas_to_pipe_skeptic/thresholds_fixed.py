# Пороги с методическими поправками скептика (fix=ALL) [О]. Запуск: python3 thresholds_fixed.py
import json, os
from g2p_skeptic import run, solve, ALL
out = {}
for nm, Q, km in (('Ughelli East', 6.7, 3.7), ('Afam Umuosi', 5.6, 0.9), ('Q5_4km', 5, 4), ('Q15_8km', 15, 8), ('Oredo_qmin', 8.5, 16.4)):
    out[nm] = dict(
        fee_from_operator=(lambda g: None if g is None else round(-g, 2))(solve(lambda x: run(Q, km, ov=dict(gas_in=x), fix=ALL)['npv20'], -15, 3)),
        Q_threshold_base=solve(lambda x: run(x, km, fix=ALL)['npv20'], 0.2, 300),
        Q_threshold_floor_escrow=solve(lambda x: run(x, km, ov=dict(gas_in=0.25, coll=1.0, dso=45, fxloss=0.0), fix=ALL)['npv20'], 0.2, 300))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'thresholds_fixed_out.json'), 'w'), indent=1, ensure_ascii=False)
for k, v in out.items(): print(k, v)
