# Скептик связки «энергия + труба» через нигерийскую оптику: платёжная дисциплина ТЭС, охрана трубы, отбор ТЭС.
# Использует bundle.model без изменений; меняются только параметры через pov/shr/price/fee_op.
# Источники поправок:
#  - Transcorp Power, I кв. 2025 [В: BusinessDay 26.04.2025]: выручка N105.4 млрд, дебиторка N352.7 млрд (+N54.3 млрд за квартал),
#    долг поставщикам газа N128.7 млрд. Отсюда [О]: собрано (105.4-54.3)/105.4 = 48% квартального счёта; DSO = 352.7/(105.4*4)*365 = 305 дн.
#  - APGC, март 2026 [В]: gencos должны поставщикам газа ~N3.3 трлн; поставки газа <43% потребности; поставщики грозят отключением.
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'power_plus_pipe'))
import bundle as B
coll_tc = (105.4-54.3)/105.4; dso_tc = 352.7/(105.4*4)*365
SITES = {k: B.SITES[k] for k in ('Ughelli East 2025 (6.7) -> Transcorp Ughelli', 'Oredo min23-25 (8.5) -> Sapele PS',
                                 'Oredo 2025 (15.5) -> Sapele PS', 'типовой 15 MMscf/д, отвод 8 км')}
NG = dict(coll=0.50, dso=300)                                  # [О из В: Transcorp I кв. 2025]
SEC = dict(sec=0.30)                                           # охрана трассы [Д 0.2-0.5; в связке взято 0.10]
CASES = {
  'база оценщика': {},
  'сбор 0.50, DSO 300 (Transcorp I кв.25)': dict(pov=NG),
  '+ охрана 0.30 $млн/г': dict(pov=NG, shr=SEC),
  '+ ТЭС берёт газ 0.60 времени (газ <43% потребности, отключения)': dict(pov=NG, shr=SEC, bo=0.60),
  'стек A: $2.68 + эскроу 0.98/45 + пол NGFCP': dict(price=2.68, gas_allin=0.40/1.05, pov=dict(coll=0.98, dso=45)),
  'стек A без эскроу (сбор 0.50/300)': dict(price=2.68, gas_allin=0.40/1.05, pov=NG),
}
def run(nm, sc, kw, T=405, e=True):
    Q, lk, pk = B.SITES[nm]
    kw = dict(kw); bo = kw.pop('bo', None)
    saved = dict(B.SCN[sc])
    if bo is not None: B.SCN[sc]['bo'] = bo
    try:
        p0 = B.model(Q, lk, pk, sc, T, e, 'power')
        b = B.model(Q, lk, pk, sc, T, e, 'bundle', **kw)
        fee = B.solve(lambda x: B.model(Q, lk, pk, sc, T, e, 'bundle', fee_op=x, **kw)['npv20'] - p0['npv20'], 0.0, 30)
    finally:
        B.SCN[sc] = saved
    return dict(npv15_pow=p0['npv15'], npv15_bun=b['npv15'], inc15=round(b['npv15']-p0['npv15'], 2),
                inc20=round(b['npv20']-p0['npv20'], 2), y2_pipe=b['y2_pipe'], fee_be_npv20=fee)
OUT = {'collection_transcorp': round(coll_tc, 3), 'dso_transcorp': round(dso_tc)}
print(f'Transcorp I кв.2025: собрано {coll_tc:.2f} счёта, DSO {dso_tc:.0f} дн. [О из В]')
for sc in ('base', 'opt'):
    for nm in SITES:
        for cn, kw in CASES.items():
            r = run(nm, sc, kw); OUT[f'{sc}|{nm}|{cn}'] = r
            print(f'{sc:4s} {nm[:30]:30s} {cn:58s} NPV15 связка {r["npv15_bun"]:6.1f} (энергия {r["npv15_pow"]:5.1f}) '
                  f'прир.15 {r["inc15"]:6.1f} прир.20 {r["inc20"]:6.1f} E2 трубы {r["y2_pipe"]:5.2f} плата-порог {r["fee_be_npv20"]}')
# сравнение «в 9-14 раз дешевле»: только компрессорный пакет на MMscf/д (без трубы, осушки, x2.6) против Seplat $0.22-0.36
import g2p as G
for sc in ('opt', 'base'):
    pp = G.SC[sc]; hpq = G.hp_per_mmscfd(pp['p_suc'], pp['p_dis'])
    comp = hpq*1.5*pp['usd_hp']/1e6
    OUT[f'compr_pkg_per_mmscfd|{sc}'] = round(comp, 3)
    print(f'{sc}: компрессорный пакет 3x50% = {comp:.2f} $млн на MMscf/д; с x{pp["cm"]}: {comp*pp["cm"]:.2f}')
json.dump(OUT, open(os.path.join(HERE, 'skeptic_ng_out.json'), 'w'), indent=1, ensure_ascii=False)
