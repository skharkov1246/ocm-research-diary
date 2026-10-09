# IRR девелопера в состоянии «всё прошло, продали, закрылись» (вероятности = 1), квартальная сетка.
import os, sys, io, contextlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):
    import dev_sale as D
def irr_t(flows):
    f = lambda r: sum(c/(1+r)**t for t, c in flows)
    lo, hi = -0.99, 20.0
    if f(lo)*f(hi) > 0: return None
    for _ in range(200):
        m = (lo+hi)/2
        if f(lo)*f(m) <= 0: hi = m
        else: lo = m
    return round(m*100, 1)
bv = D.buyer['5_UghelliEast_25km|405|E|1y']
for lvl in ('R1', 'R2'):
    gs = ['g0', 'g1'] if lvl == 'R1' else ['g0', 'g1', 'g2']
    for scn, sc in D.SC.items():
        k = sc['cost_k']; P = D.premium(lvl, sc, bv); t = D.TIMES[lvl]
        fl = [(tt, -(D.stage_cost(g, 0)*(1-k)+D.stage_cost(g, 1)*k)) for g, tt in zip(gs, t['c'])]
        sale_c = max(0.03, P*sc['sale_cost'])
        fl += [(t['sale'], P*sc['upfront']-sale_c), (t['fc'], P*(1-sc['upfront']))]
        basis = -sum(c for _, c in fl[:len(gs)]); gain = P - sale_c - basis
        fl_tax = fl + [(t['fc'], -max(0, gain)*0.34)]
        print(lvl, scn, 'премия', round(P, 2), 'затраты', round(basis, 3), 'IRR до налога', irr_t(fl), '% после налога', irr_t(fl_tax), '%')
