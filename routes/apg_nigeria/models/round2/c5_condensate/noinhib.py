# Раунд 2. Вариант «без реагентов»: холодный сепаратор ВД держим на 3 C выше температуры гидратообразования
# (Motiee 1991 [В]), ингибитор не нужен. Выход стабильного конденсата при (Th+3, P1) [О, PR-флэш, process.py].
# Запуск: <venv thermo>/bin/python noinhib.py -> noinhib_out.json
import json, os, process as P
HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
for n, d in P.COMPS.items():
    r = {}
    for P1 in (15, 30, 50):
        Th, g = P.hydrate_T(d, P1); Tx = Th + 3.0
        y = P.cold_sep_yield(d, Tx, P1)
        r[f'{P1}'] = dict(hydrate_T=Th, Tx=round(Tx, 1), yield_stab_bbl_per_MMscf=y['stab'], gas_gravity=g,
                          recovery_vs_C5p=round(y['stab']/P.potential(d)['C5p_bbl_per_MMscf'], 2))
    out[n] = r; print(n, r, flush=True)
json.dump(out, open(os.path.join(HERE, 'noinhib_out.json'), 'w'), ensure_ascii=False, indent=1)
