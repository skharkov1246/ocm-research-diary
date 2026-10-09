# Поправка аудита «спрос делить на 1.3-2.7 до контракта»: LOI на 11 МВт -> 8.5 / 4.1 МВт по договору.
import os, sys, io, contextlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):
    import dev_sale as D
for sn in ('5_UghelliEast_25km', '15_Oredo_30km_Q8.5'):
    for div in (1.0, 1.3, 2.7):
        mw = 11.0/div
        for T in (405, 440):
            b = D.val(dict(D.SITES[sn], anchor=mw), T, True, 1, D.BUY)
            o = D.val(dict(D.SITES[sn], anchor=mw), T, True, 1, D.OWN)
            legal = 'да' if mw >= 8.33 else 'НЕТ (ниже 8.33 МВт)'
            print(f'{sn} LOI/{div}: {mw:.1f} МВт, законно: {legal}; T={T}: NPV15 стратега {b["npv15"]}, капекс стратега {b["capex"]}; '
                  f'наш NPV15 {o["npv15"]} NPV20 {o["npv20"]} капекс {o["capex"]}')
