# Доп. прогон к bundle.py: «лучший реалистичный стек» условий для трубной надстройки и разбор приращения.
# Запуск: python3 stack.py -> stack_out.json. Все числа [О] на bundle.model; входы см. bundle.py.
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import bundle as B
out = {}
SITES = {k: B.SITES[k] for k in ('Ughelli East 2025 (6.7) -> Transcorp Ughelli', 'Oredo min23-25 (8.5) -> Sapele PS',
                                 'Oredo 2025 (15.5) -> Sapele PS', 'типовой 5 MMscf/д, отвод 4 км', 'типовой 15 MMscf/д, отвод 8 км')}
STACK = {
  'база': {},
  'стек A: $2.68 + эскроу ТЭС + газ по полу NGFCP': dict(price=2.68, gas_allin=0.40/1.05, pov=dict(coll=0.98, dso=45)),
  'стек A + плата оператора $1.0': dict(price=2.68, gas_allin=0.40/1.05, pov=dict(coll=0.98, dso=45), fee_op=1.0),
  'стек A + плата оператора $1.75': dict(price=2.68, gas_allin=0.40/1.05, pov=dict(coll=0.98, dso=45), fee_op=1.75),
  'стек A + ТОиР трубопровода 1.5% (вне правила аудита)': dict(price=2.68, gas_allin=0.40/1.05, pov=dict(coll=0.98, dso=45, tor=0.025)),
}
for nm, (Q, lk, pk) in SITES.items():
    for e in (True, False):
        T = 405 if e else 460          # порог LOI раунда 2 для DoP-газ + эскроу (summary_r2_out.txt)
        p0 = B.model(Q, lk, pk, 'base', T, e, 'power')
        for sn, kw in STACK.items():
            r = B.model(Q, lk, pk, 'base', T, e, 'bundle', **kw)
            out[f'{nm}|{"EDTI" if e else "noEDTI"}|T{T}|{sn}'] = dict(
                capex=r['capex'], capex_pipe=r['capex_pipe'], Qd=r['Qd_pipe'], hp=r['hp'], pipe_km=r['pipe_len'], fixed_pipe=r['fixed_pipe'],
                y2_pipe=r['y2_pipe'], npv15=r['npv15'], npv20=r['npv20'], inc15=round(r['npv15']-p0['npv15'], 2),
                inc20=round(r['npv20']-p0['npv20'], 2), power_npv15=p0['npv15'], power_npv20=p0['npv20'])
        # пороговая плата оператора поверх стека A (приращение NPV20 = 0)
        kwA = STACK['стек A: $2.68 + эскроу ТЭС + газ по полу NGFCP']
        f = B.solve(lambda x: B.model(Q, lk, pk, 'base', T, e, 'bundle', fee_op=x, **kwA)['npv20'] - p0['npv20'], 0, 15)
        out[f'{nm}|{"EDTI" if e else "noEDTI"}|T{T}|порог платы оператора поверх стека A, $/Mscf'] = f
# объёмы газа (база, 5-й год): сколько сжигает станция, сколько уходит в трубу, сколько остаётся в факеле
for nm, (Q, lk, pk) in SITES.items():
    for m in ('power', 'bundle'):
        r = B.model(Q, lk, pk, 'base', 405, True, m)
        out[f'объёмы|{nm}|{m}'] = {f'y{rw["t"]}': dict(avail=rw['avail'], burned=rw['burned'], slack=rw['slack'], pipe=rw['pipe_take'])
                                  for rw in r['rows'] if rw['t'] in (1, 2, 5, 10)}
json.dump(out, open(os.path.join(HERE, 'stack_out.json'), 'w'), indent=1, ensure_ascii=False)
for k, v in out.items(): print(k, v)
