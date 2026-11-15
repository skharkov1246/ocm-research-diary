"""Блок «thermal» для routes/gtu/organic/feed_to_funnels.json: какие оси, значения и модели добавить в воронку
двигателя (routes/gtu/funnel) и воронку камеры (routes/gtu/combustor). Все числа — из results-JSON этой папки.
Файл feed_to_funnels.json общий для узлов organic: блок пишется под своим ключом, чужие ключи не трогаются.
python3 feed_thermal.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
FEED = HERE.parent / "feed_to_funnels.json"
REL = "routes/gtu/organic/thermal/"


def load(name):
    p = HERE / name
    return json.loads(p.read_text()) if p.exists() else None


def vascular_block(v):
    if v is None:
        return None
    ngv = v["cases"]["ngv"]
    rows = ngv["rows"]
    table = {}
    for r in rows:
        if r["feasible"]:
            table.setdefault(r["T41_C"], {})[r["net"]] = round(r["x_ngv_pct"], 3)
    best = {T: min(d.items(), key=lambda kv: kv[1]) for T, d in table.items()}
    lin = v["cases"]["liner"]["rows"]
    liner = {r["net"]: dict(coolant_kg_s_m2=(round(r["coolant_per_m2_kg_s"], 4) if r.get("coolant_per_m2_kg_s") else None),
                            T_max_full_dp_C=round(r["at_full_dp"]["T_max"] - 273.15, 1),
                            spread_full_dp_K=round(r["at_full_dp"]["T_spread"], 1))
             for r in lin}
    return dict(
        x_ngv_pct_by_T41=table,
        best_network_by_T41={str(T): dict(net=b[0], x_ngv_pct=b[1]) for T, b in best.items()},
        T_lim_surface_C=v["assumptions"]["ngv"]["T_lim_C"], N_vanes_estimate=ngv["N_vanes_estimate"],
        liner_panel=liner, liner_T_lim_C=v["assumptions"]["liner"]["T_lim_C"])


def tpms_block(t):
    if t is None:
        return None
    rc = t["recuperator"]["cases"]
    out = {}
    for key, c in rc.items():
        out[key] = {}
        for name, d in c["designs"].items():
            out[key][name] = dict(mass_core_kg=round(d["mass_core"]), mass_with_headers_kg=round(d["mass_total"]),
                                  V_core_m3=round(d["V"], 3), A_front_m2=round(d["A_fr"], 2), L_m=round(d["L"], 3),
                                  modules_BLT_S1500=d["modules"], feasible=d["feasible"])
        out[key]["cycles_plate_lpbf"] = c["duty"]["cycles_lpbf"]
        if c["duty"].get("cycles_foil"):
            out[key]["cycles_foil"] = c["duty"]["cycles_foil"]     # фольга (organic/cycles) — для сравнения массы
    tr = t["transport"]["mac_summary"]
    trans = {topo: [dict(level=r["level"], eps=round(r["eps"], 3), fRe=round(r["fRe"], 2), Nu_H1=round(r["Nu_H1_stokes"], 3))
                    for r in rows] for topo, rows in tr.items()}
    return dict(designs=out, transport_stokes=trans, crosscheck=t["transport"].get("crosscheck_mac_vs_lbm"),
                inertia=t["transport"].get("inertia_used"))


def topo_block(tp):
    if tp is None:
        return None
    out = {}
    for case, c in tp["cases"].items():
        out[case] = dict(n_out=c["setup"]["n_out"], nonuniformity_max_min_at_equal_back_pressure=c["nonuniformity"],
                         balanced_dissipation_ratio_vs_conventional=dict(
                             bp_plus_orifices=round(c["balanced_ratio_bp_vs_conventional"], 3),
                             warm_start_uniform=round(c["balanced_ratio_warm_vs_conventional"], 3),
                             cold_start_uniform=round(c["balanced_ratio_uniform_vs_conventional"], 3)),
                         note="Φ_balanced — рассеяние при равных расходах + идеальные дроссели до общего давления выходов; "
                              "2D, Стокс, безразмерно — показатель, не число для KV-3")
    return out


def main():
    v = load("vascular_results.json")
    t = load("tpms_results.json")
    tp = load("topo_results.json")
    vb, tb, pb = vascular_block(v), tpms_block(t), topo_block(tp)
    block = dict(
        source=[REL + "vascular_results.json", REL + "tpms_results.json", REL + "topo_results.json"],
        date="2026-11-14",
        status="оценка уровня облика: 1D-сеть + 2D-стенка; переносы TPMS — предел Стокса (воксели 32³) с поправкой "
               "LBM; оптимизация топологии — 2D, малое Re; всё — к проверке трёхмерным расчётом (jobs_thermal.json)",
        engine_funnel=dict(
            target="routes/gtu/funnel/funnel_space.py (CAT, CONT), parts_model.py, make_model.py",
            new_axes=[
                dict(axis="cool", kind="CAT", add_value="P (печатный сопловой аппарат с каналами у стенки)",
                     note="x_ngv по T41 — таблица vascular.x_ngv_pct_by_T41 (предел поверхности — оценка); "
                          "лучшая сеть: при умеренной нагрузке — одно-два ветвления, при высокой — прямые каналы. "
                          "Изготовимость (organic/print_review): каналы 0,6 мм × 48 мм (L/D 80) под обшивкой 0,6 мм — ниже "
                          "надёжного выноса порошка (≥ 1 мм, по выдержке поиска), нужна проливка каждого канала; таблица "
                          "посчитана для гладких каналов — шероховатость печати съедает запас перепада "
                          "(print_review_results.json → roughness_ngv): при 1050 °C прямым каналам нужно 181 кПа из 470 "
                          "гладкими и 300–431 кПа при ε 30–60 мкм (держат); при 1100 °C — 285 из 456 кПа гладкими, а при "
                          "ε 30–60 мкм предел 950 °C держится только на границе или не держится. Строки 1100 и 1150 °C "
                          "годны лишь для доведённых каналов (ε ≤ 15 мкм: абразивная проливка, оценка). В точке KV-3 "
                          "(T41 1000 °C, класс U) охлаждение "
                          "не нужно: значение P окупается только вместе с подъёмом T41",
                     model="vascular.min_flow_for_limit(ngv_case(T41), сеть) → x_ngv = 2·N·ṁ/G"),
                dict(axis="ngv_network", kind="CAT", values=["прямые", "8 деревьев K=1", "4 дерева K=2"],
                     note="глубокие деревья (K ≥ 3) и растущее дерево на равном объёме предел не держат — в перебор не "
                          "включать"),
                dict(axis="recup_core", kind="CAT",
                     values=["пластины LPBF (cycles)", "фольга (cycles)", "гироид LPBF, переменная ячейка",
                             "Schwarz-D LPBF, переменная ячейка"],
                     note="масса и объём сердцевины — tpms.designs; в ламинарном режиме (предел Стокса) гироид "
                          "не легче плоских каналов при той же стенке; выигрыш TPMS — только с инерцией (Re ≥ 100), "
                          "это задача AWS thermal-tpms-lbm-*. Изготовимость (organic/print_review): печатные сердцевины "
                          "тяжелее фольги (cycles_plate_lpbf и cycles_foil в tpms.designs), стенка 0,35 / 0,25 мм ниже "
                          "наименьшего элемента 0,38 мм (по выдержке поиска), плотность тонкой стенки и вынос порошка из "
                          "каналов Dh 1,5 мм длиной 0,1–0,26 м не подтверждены — базовое значение оси: фольга"),
                dict(axis="manifold", kind="CAT", values=["гребёнка/спицы из труб", "печатный (оптимизация топологии)"],
                     note="печатный коллектор — 1 деталь вместо N труб + коллектора + 2N соединений; равномерность "
                          "раздачи и потери — topo (2D, малое Re: показатель, не число для KV-3). Разброс topo — при "
                          "одинаковом противодавлении выходов без дросселей; где выходы дросселируются (форсунки, сопла), "
                          "выигрыш — только в числе деталей и соединений"),
            ],
            models=dict(vascular=vb, tpms=tb, manifold=pb)),
        combustor_funnel=dict(
            target="routes/gtu/combustor (targets.json, comb_space.py)",
            new_axes=[
                dict(axis="liner_cooling", kind="CAT", values=["плёночное (как есть)", "печатные каналы у стенки"],
                     note="при перепаде 3 % p3 предел 900 °C не держит ни одна сеть (vascular.liner_panel: coolant_kg_s_m2 "
                          "= null, лучшие — прямые каналы, T_max при полном перепаде в T_max_full_dp_C) — значение «печатные "
                          "каналы» допустимо только вместе с плёнкой или покрытием (оценка; излучение пламени — ε = 0,20)"),
                dict(axis="fuel_manifold", kind="CAT", values=["трубы", "печатный с равномерной раздачей"],
                     note="разброс расходов по форсункам брать из structure/manifold_results.json (с перепадом форсунки "
                          "4 % p3: кольцо 1,8 %, дерево 0 %; от допусков печати p95 2,3 % и 0,32 %; блок structure этого "
                          "файла). topo.nonuniformity (61 % у гребёнки) — при одинаковом противодавлении без форсунок: к "
                          "топливному коллектору KV-3 не применять — перепад форсунки (~55 кПа) в ~30 раз больше перепада "
                          "каналов (~2 кПа) и сам выравнивает раздачу"),
            ]),
        aws=dict(jobs=REL + "jobs_thermal.json",
                 note="5 задач, проверено командой plan без AWS; прогноз пакета — поле estimate (верхняя граница); "
                      "трёхмерные CFD/CHT — блок planned_3d (кода сборки кейсов нет)"),
    )
    (HERE / "feed_thermal.json").write_text(json.dumps(block, ensure_ascii=False, indent=1))
    feed = json.loads(FEED.read_text()) if FEED.exists() else {}
    feed["thermal"] = block
    FEED.write_text(json.dumps(feed, ensure_ascii=False, indent=1))
    print("записано: thermal в", FEED.name, "; ключи файла:", list(feed))


if __name__ == "__main__":
    main()
