# Доп. пороги к dev_sale.py: премия при вероятностях = 1, и общая вероятность успеха при «аукционной» премии.
import os, sys, io, contextlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):
    import dev_sale as D
def solve(f, lo, hi):
    if f(lo)*f(hi) > 0: return None
    for _ in range(60):
        m = (lo+hi)/2
        if f(lo)*f(m) <= 0: hi = m
        else: lo = m
    return round((lo+hi)/2, 3)
bv = D.buyer['5_UghelliEast_25km|405|E|1y']
for lvl in ('R1', 'R2'):
    for scn in ('opt', 'base', 'pess'):
        sc = D.SC[scn]
        one = dict(sc, p0=1, p1=1, p2=1, ps_R1=1, ps_R2=1, pfc=1)
        for r in (0.15, 0.20):
            fP = lambda P: D.dev_ev(lvl, dict(one, **{f'alpha_{lvl}': 1.0, f'cap_{lvl}': P/bv['capex_buy']}), dict(bv, npv15_buy=1e9), r, sc['cost_k'])[0]
            Pmin = solve(fP, 0.0, 50.0)
            # вероятность успеха (общая p0*p1*p2*ps, pfc из сценария) при премии = 0.5 x NPV15 стратега (аукцион, без потолка)
            Pa = 0.5*bv['npv15_buy']
            def fp(x):
                s2 = dict(sc, p0=x, p1=1, p2=1, ps_R1=1, ps_R2=1, **{f'alpha_{lvl}': 0.5, f'cap_{lvl}': 1.0})
                return D.dev_ev(lvl, s2, bv, r, sc['cost_k'])[0]
            px = solve(fp, 0.0005, 1.0)
            print(f'{lvl} {scn} r={r}: премия для EV=0 при вероятностях 1 = ${Pmin} млн ({Pmin/bv["capex_buy"]*100 if Pmin else 0:.1f}% капзатрат стратега); '
                  f'при аукционной премии ${Pa:.2f} млн нужна P(готово и продано) >= {px}')
print('базовая P(готово и продано): R1', round(0.2*0.55*0.30, 3), 'R2', round(0.2*0.55*0.5*0.5, 3))
