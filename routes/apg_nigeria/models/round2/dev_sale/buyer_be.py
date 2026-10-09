# Безубыточный тариф стратега против нашего (₦/кВт·ч), DoP-газ + эскроу, газ $1.0 всё вкл.
import os, sys, io, contextlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):
    import dev_sale as D
for sn in ('5_UghelliEast_25km', '15_Oredo_30km_Q8.5'):
    for e in (True, False):
        for d in (1, 2):
            kw = D.kw_e(e)
            if d == 2: kw['capex_split'] = [0.5, 0.5]
            tl = dict(D.TL, delay=d)
            r = {}
            for who, ex in (('наш', D.OWN), ('стратег', D.BUY)):
                p = dict(D.SITES[sn], gas=D.gas_from_allin(1.0), **ex)
                r[who] = (round(D.be(p, D.S['dopgas'], tl, 'npv15', **kw)), round(D.be(p, D.S['dopgas'], tl, 'npv20', **kw)))
            # стратег + премия девелопера $1.0 млн в капитале
            p = dict(D.SITES[sn], gas=D.gas_from_allin(1.0), **dict(D.BUY, dev=0.3+1.0))
            r['стратег+премия1.0'] = (round(D.be(p, D.S['dopgas'], tl, 'npv15', **kw)), round(D.be(p, D.S['dopgas'], tl, 'npv20', **kw)))
            print(sn, 'EDTI' if e else 'noEDTI', f'{d}г', 'BE15/BE20:', r)
