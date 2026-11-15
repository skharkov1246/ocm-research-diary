"""Фикстуры теста эскизов: геометрия уровня 1 для классов формы, которых может не оказаться в отладочном топ-50.

Каждая фикстура - настоящая концепция выборки уровня 0 (funnel_space.concept по idx, n, kv, seed прогона; из первых
--per-class концепций класса - наименее нарушающая правила), её расчётная точка
уровня 1 без оптимизации (funnel_l1.evaluate_point) и компоновка funnel_l1.layout - та же функция, что пишет поле geometry у
каждого из 50. Ничего не придумывается: размеры - из кода уровня 1; правила допустимости точка может не проходить (это не
вариант отбора, а образец формы для теста генератора эскизов). В рейтинг и в 50 фикстуры не входят.

Запуск (python с cantera и scipy - тот же, что для funnel_l1.py), воспроизводимо, без случайности сверх зерна выборки:
    python routes/gtu/funnel/sketch/make_fixtures.py            # -> sketch/fixtures_l1_layout.json
"""
import argparse, json, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
FUNNEL = os.path.dirname(HERE)
if FUNNEL not in sys.path:
    sys.path.insert(0, FUNNEL)
import funnel_space as FS   # noqa: E402
import funnel_l1 as L1      # noqa: E402
import merge_l1 as M1       # noqa: E402

KV, N = "KV-3", 2000          # выборка отладочного прогона AWS уровня 0 (вход funnel_l1_top50.json)
CLASSES = [   # (ключ, описание, условие на дискретные оси концепции)
    ("recup_integral", "печатный рекуператор в корпусе (кольцевой)", lambda c: c["recup"] == "печатный" and c["recup_integral"] == 1),
    ("recup_side", "печатный рекуператор блоком модулей сбоку", lambda c: c["recup"] == "печатный" and c["recup_integral"] == 0),
    ("cluster4", "кластер из 4 ядер", lambda c: c["cores"] == 4),
    ("cluster2", "кластер из 2 ядер", lambda c: c["cores"] == 2),
    ("comb_ext", "выносная камера сгорания", lambda c: c["comb"] == "выносная"),
    ("comb_wave", "волновой ротор", lambda c: c["comb"] == "ВР"),
    ("comb_rdc", "камера с ростом давления (RDC)", lambda c: c["comb"] == "RDC"),
    ("radial_1shaft_oil", "радиальная турбина, один вал, масляные подшипники, камера не выносная",
     lambda c: c["turb"] == "Р" and c["rotors"] == 1 and c["bear"] == "масло" and c["comb"] != "выносная"),
    ("axial_centrifugal", "осевые ступени + центробежная (ОЦ)", lambda c: c["comp"] == "ОЦ"),
    ("diagonal", "диагональное колесо", lambda c: c["comp"] == "Д"),
    ("double_inlet_1shaft", "двусторонний вход, один вал (выходной вал со стороны входа)",
     lambda c: c["inlet"] == "двух" and c["rotors"] == 1),
]


def entry_of(c, ev, key, text):
    lay = L1._clean(ev["layout"])
    mods = [dict(name=m["name"], mass_kg=round(float(m["mass_kg"]), 1)) for m in ev["models"]["sv"]["modules"]]
    env = lay["envelope"]
    r = dict(c)
    return dict(fixture=True, fixture_class=key, fixture_text=text, div_rank=None, idx=int(c["idx"]), kv=KV, arch=c["arch"],
                archetype=FS.ARCH[c["arch"]]["name"], scheme=M1.scheme(r),
                concept={k: c[k] for k in ("idx", "n", "seed", "kv", "arch", "rotors", "comp", "n_ax", "inlet", "turb", "nt", "npt",
                                           "cool", "recup", "comb", "bear", "cores", "rotor_join", "m_comb", "hot_integral",
                                           "recup_integral", "conn", "cartridge", "eps_rec")},
                size_m=dict(L=round(env["L"], 3), W=round(env["W"], 3), H=round(env["H"], 3)),
                mass_kg_estimate=round(float(lay["mass_kg_estimate"]), 0), modules=mods, geometry=lay,
                rules_failed=sorted({b["code"] for b in ev["bad"]}))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=os.path.join(HERE, "fixtures_l1_layout.json"))
    ap.add_argument("--scan", type=int, default=N, help="сколько idx выборки просматривать по порядку")
    ap.add_argument("--per-class", type=int, default=25, help="из скольких первых концепций класса брать наименее нарушающую правила")
    a = ap.parse_args(argv)
    t0 = time.time()
    L1._prep(N, FS.SEED, KV)
    tgt = FS.target(KV)
    found, best, tried = {}, {}, {}
    for idx in range(a.scan):
        if all(tried.get(k, 0) >= a.per_class for k, _, _ in CLASSES):
            break
        c = FS.concept(idx, N, KV, FS.SEED)
        keys = [(k, t) for k, t, cond in CLASSES if tried.get(k, 0) < a.per_class and cond(c)]
        if not keys:
            continue
        ev = L1.evaluate_point(c, tgt)
        for key, text in keys:
            tried[key] = tried.get(key, 0) + 1
            if ev["ok"] and (key not in best or (ev["viol"], idx) < best[key][0]):
                best[key] = ((ev["viol"], idx), entry_of(c, ev, key, text))
    for key, (rank, e) in best.items():
        e["viol"] = round(float(rank[0]), 5)
        found[key] = e
    out = dict(title="Фикстуры теста эскизов: компоновка уровня 1 для классов формы",
               method="funnel_space.concept(idx, n, kv, seed) -> funnel_l1.evaluate_point (без оптимизации) -> funnel_l1.layout; "
                      f"из первых {a.per_class} по idx концепций класса - с наименьшим превышением правил (viol; при равенстве - "
                      "меньший idx), у которой расчётная точка строится",
               software=f"python {sys.version.split()[0]}", date="2026-11-14", kv=KV, n=N, seed=FS.SEED, code_sha=L1.code_sha(),
               note="образцы формы для теста генератора эскизов; не варианты отбора, правила уровня 1 могут не проходить (rules_failed)",
               classes=[dict(key=k, text=t, found=k in found, concepts_tried=tried.get(k, 0)) for k, t, _ in CLASSES],
               entries=[found[k] for k, _, _ in CLASSES if k in found], t_s=round(time.time() - t0, 1))
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"{a.out}: {len(found)} из {len(CLASSES)} классов, {out['t_s']} с")
    for cl in out["classes"]:
        print(f"  {cl['key']:22s} {'есть' if cl['found'] else 'НЕТ '} (попыток {cl['concepts_tried']})")


if __name__ == "__main__":
    main()
