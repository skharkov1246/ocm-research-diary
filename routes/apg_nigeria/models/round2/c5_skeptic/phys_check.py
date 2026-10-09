# Скептик раунда 2, «Конденсат C5+ (JT)». Проверка физики [О]:
# 1) повтор выхода без ингибитора (Th+3 C) на 30 бар; 2) верхняя граница продукта, если весь C4 жидкости
#    оставить в продукте (спайкинг бутана в нефть под давлением насоса) — стабилизатор оценки держит только 15/30%;
# 3) адиабатический JT без рекуперации. Запуск: <venv thermo>/bin/python phys_check.py
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'c5_condensate'))
import process as P
out = {}
for n in ('SPE_AG2', 'SPE_AG4', 'SPE_AG6', 'NJTD_A'):
    d = P.COMPS[n]; r = {}
    for P1 in (30, 50):
        Th, g = P.hydrate_T(d, P1); Tx = Th + 3.0
        base = P.cold_sep_yield(d, Tx, P1)['stab']
        old = dict(P.SPLIT)
        P.SPLIT.update(iC4=1.0, nC4=1.0)              # весь C4 жидкости в продукт [Д верх]
        allc4 = P.cold_sep_yield(d, Tx, P1)['stab']
        P.SPLIT.clear(); P.SPLIT.update(old)
        r[f'{P1}bar'] = dict(Th=Th, Tx=round(Tx, 1), y_stab=base, y_allC4=allc4)
    r['jt_adiab_30to3_from35'] = P.jt_T(d, 35, 30, 3)
    out[n] = r; print(n, r, flush=True)
json.dump(out, open(os.path.join(HERE, 'phys_check_out.json'), 'w'), ensure_ascii=False, indent=1)
