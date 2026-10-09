# Регрессия: deal_r2 при параметрах раунда 1 (OLD, TL_OLD) должен повторить judge_check.py / deal2.py.
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from deal_r2 import *
S = STRUCTS
def r(pr, s, **kw):
    x = run(dict(OLD, p_dir=pr), S[s], TL_OLD, debt=(0.6, 0.14, 7), r1_round=True, **kw)
    return {k: x[k] for k in ('capex', 'y1', 'y2', 'npv15', 'npv20', 'irr', 'min_dscr')}
exp = {  # judge_check.py, раунд 1
 ('dopmargin', 377, True): (12.9, 4.8), ('dopmargin', 377, False): (7.5, -0.1), ('esc', 377, True): (8.7, 1.4),
 ('esc', 377, False): (3.6, -3.2), ('none', 377, True): (6.6, -0.4), ('dopmargin', 350, True): (9.3, 1.9)}
ok = True
for (s, pr, e), (a, b) in exp.items():
    x = r(pr, s, edti=e, depyrs=5 if e else 10)
    good = abs(x['npv15']-a) < 0.051 and abs(x['npv20']-b) < 0.051
    ok &= good; print('OK ' if good else 'BAD', s, pr, 'EDTI' if e else 'noEDTI', x, 'ожидалось', a, b)
for s, e, exp_be in [('dopmargin', True, (290, 332)), ('dopmargin', False, (316, 377)), ('esc', True, (317, 362)),
                     ('esc', False, (345, 411)), ('none', True, (330, 380)), ('none', False, (351, 420))]:
    kw = dict(edti=e, depyrs=5 if e else 10, r1_round=True)
    a = round(be(OLD, S[s], TL_OLD, 'npv15', **kw)); b = round(be(OLD, S[s], TL_OLD, 'npv20', **kw))
    good = abs(a-exp_be[0]) <= 1 and abs(b-exp_be[1]) <= 1   # раунд 1 округлял NPV до 0.1 внутри бисекции -> ±1 ₦
    if s == 'none' and not e:   # раунд 1 (deal3.py) считал эту ячейку с неплатежами 2% вместо 5%
        a2 = round(be(OLD, dict(S[s], bad=0.02), TL_OLD, 'npv15', **kw)); b2 = round(be(OLD, dict(S[s], bad=0.02), TL_OLD, 'npv20', **kw))
        good = (a2, b2) == exp_be
        print('   ячейка раунда 1 воспроизводится только при bad=0.02:', (a2, b2), '; при bad=0.05 верно', (a, b))
    ok &= good
    print('OK ' if good else 'BAD', 'BE', s, 'EDTI' if e else 'noEDTI', (a, b), 'ожидалось', exp_be)
x = run(dict(OLD, p_dir=377), S['dopmargin'], dict(TL_OLD, delay=2), edti=False, depyrs=10, r1_round=True)
print('delay2 377 noEDTI', x['npv15'], x['npv20'], '(судья: 2.4 / -5.2)')
print('REGRESSION', 'PASS' if ok else 'FAIL')
