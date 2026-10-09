# Скептик BOOT, оптика «нигерийская реальность»: плательщик, неплатежи, доступ, права, безопасность.
# Берёт boot.py как есть (поправки аудита NaCN там уже стоят) и добавляет нигерийские слои [Д]:
#  N1 кредит NEPL/JV: сбор 0.80, дебиторка 240 дн (история cash call NNPC ~$8.5 млрд к 2016 [В]);
#  N2 риск разрыва контракта (смена владельца, отзыв лицензии, NGFCP-конфликт, блокада): hazard h/год, актив
#     застревает (передача за $0 и так), денежный поток лет эксплуатации x (1-h)^(t-1);
#  N3 охрана/общины: sec 0.8 млн/г вместо 0.35 (эскорт, JTF, GMoU) и готовность 0.85 вместо 0.92 (блокады);
#  N4 отбор газа оператором (вход завода/экспорт в домашний рынок с неплатежами GenCo): загрузка x0.9;
#  N5 тариф в найре, если долларовый внутренний договор с NEPL оспорят (CBN Act, аудит NPDC 2023 [В]): 8-16%/г.
# Запуск: python3 sk_ng.py -> sk_ng_out.json
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'compression_boot'))
import boot as B

def npv(r, cf): return sum(c/(1+r)**t for t, c in enumerate(cf))

def run_ng(Q, fee, sc='base', ov=None, h=0.0, ngn=None, top=None, conf='boost', km=None):
    r = B.run(Q, conf, sc, fee, ov=ov, ngn=ngn, top=top, km=km)
    b = B.SC[sc]['build'] if not (ov and 'build' in ov) else ov['build']
    cf = list(r['cf'])
    for k in range(b, len(cf)):
        cf[k] *= (1-h)**(k-b)
    return dict(npv15=round(npv(.15, cf), 2), npv20=round(npv(.20, cf), 2), e2=r['y2']['ebitda'], capex=r['capex']['total'])

def thr(Q, **kw):
    return B.solve(lambda x: run_ng(Q, x, **kw)['npv20'], 0.01, 60)

N1 = dict(coll=0.80, dso=240)
N3 = dict(sec=0.80, up=0.85)
N4 = dict(lf=0.78*0.9)
SEPLAT = dict(cm=0.6, staff=0.25, ovh=0.10, dev=0.2)          # «подрядчик уровня оператора» из оценки
ENERFLEX = dict(cm=1.6, staff=0.25, ovh=0.10)

def merge(*ds):
    o = {}
    for d in ds: o.update(d)
    return o

out = {}
def pr(k, v): out[k] = v; print(f'{k:60s} {v}')

for Q in (5, 15):
    for lab, ov, h, ngn in (
        ('model_base', None, 0.0, None),
        ('N1_credit', N1, 0.0, None),
        ('N2_hazard5', None, 0.05, None), ('N2_hazard10', None, 0.10, None),
        ('N3_security_up', N3, 0.0, None),
        ('N4_offtake', N4, 0.0, None),
        ('NG_real(N1+N2_7%+N3+N4)', merge(N1, N3, N4), 0.07, None),
        ('NG_real+naira12%', merge(N1, N3, N4), 0.07, 0.12),
        ('seplat_like_model', SEPLAT, 0.0, None),
        ('seplat_like+NG_real', merge(SEPLAT, N1, N3, N4), 0.07, None),
        ('seplat_like+NG_real+naira12%', merge(SEPLAT, N1, N3, N4), 0.07, 0.12),
        ('enerflex_like+NG_real', merge(ENERFLEX, N1, N3, N4), 0.07, None),
    ):
        pr(f'Q{Q}|{lab}|fee1.75', run_ng(Q, 1.75, ov=ov, h=h, ngn=ngn))
        pr(f'Q{Q}|{lab}|fee_npv20=0', thr(Q, ov=ov, h=h, ngn=ngn))

# Именованные площадки NEPL, газ по минимуму 2023-25, NG_real
for name in ('Oredo (NEPL, IGHF)', 'Odidi (NEPL)', 'Ughelli East (NEPL)', 'Obiafu-Obrikom (Oando)'):
    q, qmin, conf, km = B.NAMED[name]
    ov = merge(N1, N3, dict(lf=0.78*0.9*qmin/q))
    pr(f'NAMED {name}|qmin|NG_real|fee_npv20=0', B.solve(lambda x: run_ng(q, x, ov=ov, h=0.07, conf=conf, km=km)['npv20'], 0.01, 80))

# Мотив оператора: годовой штраф по 9 площадкам и доля NEPL (VIIRS 2025, [О])
vol = {k: v[0] for k, v in B.NAMED.items()}
tot = sum(vol.values()); nepl = sum(v for k, v in vol.items() if 'NEPL' in k)
pr('sites_total_MMscfd_2025', round(tot, 1)); pr('NEPL_share', round(nepl/tot, 3))
for p in (1.6, 1.9, 3.5):
    pr(f'penalty_all9_$M_per_y_at_{p}', round(tot*365*p/1000, 1))
# Средний штраф NUPRC 2025: N521.87 млрд / 203.97 млрд scf [В]; курс 2025 [Д] 1 450-1 600
for fx in (1330, 1450, 1530, 1600):
    pr(f'avg_penalty_$per_Mscf_fx{fx}', round(521.87e9/203.97e6/fx, 2))
# Если всё сожжённое облагалось бы $3.50, доля фактически оплаченного объёма:
for fx in (1450, 1530):
    pr(f'implied_billed_share_at_3.50_fx{fx}', round(521.87e9/fx/3.5/203.97e6, 2))
json.dump(out, open(os.path.join(HERE, 'sk_ng_out.json'), 'w'), indent=1, ensure_ascii=False)
