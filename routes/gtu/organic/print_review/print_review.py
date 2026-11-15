#!/usr/bin/env python3
"""Скептик изготовимости печатью и пользы для владельца: генеративные формы organic/structure и organic/thermal.

Что считается (всё — из файлов узлов, внешние правила печати — с источником и статусом):
  1. STL узлов: доля поверхности с нависанием круче 45° к плите (лучшая из 6 ориентаций по осям), число отдельных
     оболочек и замкнутых внутренних полостей (порошок из них не выходит без отверстий).
  2. Наименьшие стенки, отверстия, каналы и отношение длины к диаметру — из results-JSON и констант узлов — против
     правил печати и выноса порошка.
  3. Чувствительность «сосудов» сопловой лопатки (thermal/vascular.py) к шероховатости печатного канала: тот же
     решатель, трение Churchill с ε/D, теплоотдача Gnielinski × (f/f_гл)^n (Norris), предел 950 °C.
  4. Время печати сердцевин рекуператора (thermal/tpms_results.json) по скорости камер funnel/make_model.PRINTERS.
  5. Собираемость оболочки «одной деталью» (structure/shell_results.json → parts): ротор из 2 колёс и 3 ступеней
     турбины — какой статор запирает какое колесо (frame3), сколько частей корпуса нужно на самом деле; поправка
     счёта деталей по funnel/parts_model.py.
  6. Цельная горячая часть (жаровая труба + сопловой аппарат + корпус): несовместность температурных деформаций.
  7. Посадка в камеры ГИП и проход лучей КТ через стенки оболочки.
  8. Польза для владельца по каждой форме (детали, масса, КПД) и решение: оставить / с условиями / отсеять.
python3 print_review.py                # ~4–20 мин (п. 3 — 42 решения vascular), пишет print_review_results.json рядом
python3 print_review.py --reuse-rough  # без п. 3: строки шероховатости берутся из прошлого прогона
"""
import json
import math
import struct
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ORG = HERE.parent
GTU = ORG.parent
ROOT = GTU.parent.parent
sys.path.insert(0, str(ORG / "thermal"))
sys.path.insert(0, str(GTU / "funnel"))

# ------------------------------------------------------------------ внешние правила (источник и статус)
RULES = {
    "overhang_deg": dict(value=45.0, text="самонесущее нависание — не круче 45° к плите (у 316L и бронз — 55°)",
                         src="https://www.protolabs.com/resources/design-tips/how-to-3d-print-metal-parts/",
                         status="по выдержке поиска"),
    "hole_horizontal_max_mm": dict(value=8.0, text="горизонтальное круглое отверстие до 8 мм строится без поддержек; "
                                   "больше — каплевидный или ромбовидный профиль", src="https://www.makerverse.com/resources/"
                                   "the-laser-powder-bed-fusion-l-pbf-design-guide/; Protolabs DMLS design guide",
                                   status="по выдержке поиска"),
    "min_feature_in718_mm": dict(value=0.381, text="наименьший элемент IN718 при обычном разрешении 0,015 in",
                                 src="https://www.protolabs.com/media/bn1dvyuy/pl_3dp_dmls_design_guide_en.pdf",
                                 status="по выдержке поиска"),
    "depowder_ideal_min_mm": dict(value=1.0, text="канал от 1 мм — удаление порошка надёжно; меньше — зависит от формы и "
                                  "порошка (провайдер очищал 0,4 мм)", src="https://3dadept.com/designing-for-depowdering-"
                                  "optimizing-3d-printed-parts-from-the-start/", status="по выдержке поиска"),
    "depowder_long_mm": dict(value=[3.0, 5.0, 50.0], text="канал ≥ 3 мм, при длине больше 50 мм — ≥ 5 мм (правило "
                             "поставщиков, один первоисточник на нескольких сайтах)",
                             src="https://www.assemblymag.com/articles/97510-design-for-additive-manufacturing-part-3",
                             status="по выдержке поиска"),
    "escape_holes": dict(value=2, text="не меньше двух отверстий выноса порошка на полость", src="3dadept.com (та же статья)",
                         status="по выдержке поиска"),
    "roughness_channels": dict(text="печатные каналы DMLS: трение растёт сильнее теплоотдачи, особенно в мини-каналах малого Dh; "
                               "средняя шероховатость Ra плохо предсказывает потери — нужны Rq и асимметрия профиля",
                               src="Stimpson, Snyder, Thole, Mongillo, J. Turbomach. 2016 (doi:10.1115/1.4032167); "
                                   "Wildgoose et al., J. Turbomach. 2024", status="по выдержке поиска (аннотации)"),
    "fatigue_ti64_lpbf": dict(text="усталость LPBF Ti-6Al-4V: поверхность после печати снижает предел выносливости ~45 %; "
                              "обработка поднимает его примерно втрое, ГИП — ещё на 10–20 %; база 10⁷ — 200–500 МПа "
                              "в зависимости от дефектов", src="Kahlin et al. (обзор); Ednie et al. JMRT 2024; "
                              "ASM IMAT 2026 abstract 65732", status="по выдержке поиска"),
    "hip_largest": dict(value_m=[2.05, 1.60, 0.88],
                        text="наибольшая камера ГИП в мире — Ø2,05 × 4,2 м (Metal Technology Co., Япония); Ø1,6 × 3,5 м "
                             "(MTC Powder Solutions, Швеция); серийная QIH 200 URC — Ø0,88 × 2,15 м",
                        src="https://quintustechnologies.com/products/material-densification/hot-isostatic-presses/large-hip-systems/; "
                            "https://www.prweb.com/releases/addition-of-super-sized-quintus-hot-isostatic-press-marks-major-"
                            "expansion-for-mtc-powder-solutions-842093421.html", status="по выдержке поиска"),
    "ct_penetration": dict(text="при 450 кВ и ниже просвечивание крупных плотных деталей затруднено (нужен ускоритель "
                           "9 МэВ); никель при 200 кВ — порядка 2,1 см", src="https://4nsi.com/scanning-service/high-energy-x-ray/; "
                           "https://npl.co.uk/products-services/dimensional/mpp-xct", status="по выдержке поиска"),
    "printer_boxes": dict(text="BLT-S1500: 1500×1500×1200 мм (DirectIndustry: Z 1200) или ×1500 (TCT, реселлер) — "
                          "источники расходятся, принята меньшая; Farsoon FS1521M-U: 1530×1530×1650 мм (Metal AM, 2025); "
                          "EP-M2050: 2050×2050×1100 мм (funnel/make_model.PRINTERS)",
                          src="https://www.tctmagazine.com/blt-guide/; https://www.metal-am.com/farsoon-fs1521m-u-now-"
                              "supports-up-to-32-lasers/", status="по выдержке поиска"),
}

# ------------------------------------------------------------------ 1. STL
STLS = [  # путь, вид (solid — металл, fluid — канал как тело), что это
    ("structure/stl/comp2_profile_eq.stl", "solid", "колесо 2, гладкая задняя поверхность (проверено ccx)"),
    ("structure/stl/comp2_A_mass_eq.stl", "solid", "колесо 2, воксельный проект A с полостями (в ccx не проходит)"),
    ("structure/stl/comp2_B_mass_eq_gyroid_sector40.stl", "solid", "колесо 2, сектор 40° с гироидом (не проходит)"),
    ("structure/stl/turb1_profile_eq.stl", "solid", "диск турбины 1, гладкий профиль eq (проверено)"),
    ("structure/stl/turb1_profile_lim.stl", "solid", "диск турбины 1, «предел» (нужна доводка)"),
    ("structure/stl/turb1_profile_baseline_frame3.stl", "solid", "диск турбины 1, исходный литой профиль"),
    ("structure/stl/shell_frame_kv3.stl", "solid", "рама опоры подшипника (оптимизация топологии, 2D, выдавлено)"),
    ("thermal/tpms_gyroid_optimum_block.stl", "solid", "лист гироида рекуператора (упрощённая копия)"),
    ("thermal/tpms_gyroid_graded_illustration.stl", "solid", "иллюстрация переменной ячейки (не оптимум)"),
    ("thermal/topo_comb6.stl", "fluid", "коллектор 1→6 (2D, выдавлено, каналы телом)"),
    ("thermal/topo_fan8.stl", "fluid", "коллектор 1→8 (2D, выдавлено, каналы телом)"),
    ("thermal/vascular_ngv_straight.stl", "fluid", "каналы лопатки СА: прямые"),
    ("thermal/vascular_ngv_tree_r8.stl", "fluid", "каналы лопатки СА: 8 деревьев"),
    ("thermal/vascular_ngv_tree_r4.stl", "fluid", "каналы лопатки СА: 4 дерева"),
    ("thermal/vascular_ngv_cco.stl", "fluid", "каналы лопатки СА: растущее дерево"),
]


def read_stl(path):
    b = Path(path).read_bytes()
    n = struct.unpack("<I", b[80:84])[0]
    if 84 + 50 * n == len(b):
        rec = np.frombuffer(b[84:84 + 50 * n], dtype=np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")]))
        return rec["v"].astype(float)
    tri, cur = [], []                                          # ASCII
    for line in b.decode(errors="ignore").splitlines():
        s = line.strip().split()
        if s and s[0] == "vertex":
            cur.append([float(x) for x in s[1:4]])
            if len(cur) == 3:
                tri.append(cur); cur = []
    return np.array(tri, float)


def stl_audit(path, kind):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    V = read_stl(path)
    e1, e2 = V[:, 1] - V[:, 0], V[:, 2] - V[:, 0]
    cr = np.cross(e1, e2)
    A2 = np.linalg.norm(cr, axis=1)
    ok = A2 > 0
    V, cr, A2 = V[ok], cr[ok], A2[ok]
    N = cr / A2[:, None]
    area = A2 / 2
    vol_total = float(np.einsum("ij,ij->i", V[:, 0], np.cross(V[:, 1], V[:, 2])).sum() / 6.0)
    inverted = vol_total < 0                               # нормали внутрь (обход по часовой) — дефект файла
    if inverted:
        N = -N
        V = V[:, ::-1]
    if kind == "fluid":
        N = -N                                             # стенка металла смотрит внутрь канала
    lo, hi = V.reshape(-1, 3).min(0), V.reshape(-1, 3).max(0)
    span = hi - lo
    tol = 1e-4 * float(span.max())
    best = None
    rows = {}
    for ax in range(3):
        for sgn in (1, -1):
            nz = sgn * N[:, ax]
            zmin = (sgn * V[:, :, ax]).min()
            on_plate = np.all(np.abs(sgn * V[:, :, ax] - zmin) < tol, axis=1)
            down = (nz < -math.cos(math.radians(RULES["overhang_deg"]["value"]))) & ~on_plate
            frac = float(area[down].sum() / area.sum())
            key = f"{'+' if sgn > 0 else '-'}{'xyz'[ax]}"
            rows[key] = round(frac, 4)
            if best is None or frac < best[1]:
                best = (key, frac)
    # оболочки и полости: склейка вершин, связность по общим вершинам, знак объёма оболочки
    q = np.round(V.reshape(-1, 3) / max(tol, 1e-9)).astype(np.int64)
    _, inv = np.unique(q, axis=0, return_inverse=True)
    inv = inv.reshape(-1, 3)
    nt = len(inv)
    r = np.repeat(np.arange(nt), 3)
    g = coo_matrix((np.ones(3 * nt), (r, inv.ravel())), shape=(nt, inv.max() + 1)).tocsr()
    adj = g @ g.T
    ncomp, lab = connected_components(adj, directed=False)
    vol = np.zeros(ncomp)
    np.add.at(vol, lab, np.einsum("ij,ij->i", V[:, 0], np.cross(V[:, 1], V[:, 2])) / 6.0)
    big = np.abs(vol) > 1e-6 * np.abs(vol).max()
    voids = int(np.sum(big & (vol < 0)))              # оболочка с нормалями внутрь — полость внутри тела
    bodies = int(np.sum(big & (vol > 0)))             # у fluid — отдельные куски каналов (трубки сегментов)
    return dict(normals_inverted=bool(inverted), triangles=int(nt), bbox=[round(float(x), 3) for x in span], area=round(float(area.sum()), 2),
                overhang_frac_by_build_dir=rows, best_build_dir=best[0], overhang_frac_best=round(best[1], 4),
                shells=int(ncomp), bodies=bodies, enclosed_voids=voids)


# ------------------------------------------------------------------ 2. наименьшие элементы форм
def features():
    v = json.loads((ORG / "thermal" / "vascular_results.json").read_text())
    a = v["assumptions"]["ngv"]
    al = v["assumptions"]["liner"]
    import vascular as VS
    Lx = a["side_len_over_c"] * a["chord"]
    pt = v["inputs"][0]["used"]
    H = pt["h_blade"]
    nets = {}
    s = VS.net_straight(Lx, H, a["N_straight"], a["d_straight"])
    Vt = a["N_straight"] * math.pi * a["d_straight"] ** 2 / 4 * Lx
    nets["straight"] = s
    for name, nr, ex in (("tree_r8", 8, 3.0), ("tree_r4", 4, 3.0), ("tree_r1", 1, 3.0)):
        K = int(round(math.log2(a["N_straight"] / nr)))
        t = VS.net_delta(Lx, H, K, 1e-3, ex, name, n_root=nr)
        t.scale_d(math.sqrt(Vt / t.volume()))
        nets[name] = t
    ngv = {}
    for k, n in nets.items():
        d = [s_[2] for s_ in n._ch()]
        ngv[k] = dict(d_min_mm=round(1e3 * min(d), 3), d_max_mm=round(1e3 * max(d), 3),
                      path_len_mm=round(1e3 * Lx, 1), L_over_dmin=round(Lx / min(d), 1))
    lat = json.loads((ORG / "structure" / "lattice_results.json").read_text())
    sheet = lat.get("sheet_thickness_mm")
    t = json.loads((ORG / "thermal" / "tpms_results.json").read_text())
    pr = t["recuperator"]["print_rules"]
    rc = t["recuperator"]["cases"]["R90-ic-RH-ibc-none-none"]["designs"]
    man = json.loads((ORG / "structure" / "manifold_results.json").read_text())
    sh = json.loads((ORG / "structure" / "shell_results.json").read_text())["geometry"]
    rmin = RULES
    out = dict(
        ngv_vascular=dict(nets=ngv, skin_mm=1e3 * a["t_skin"], wall_mm=1e3 * a["t_w"], T_surface_lim_C=a["T_lim_C"],
                          verdict="прямые каналы 0,6 мм × 48 мм (L/D 80) и обшивка 0,6 мм: ниже надёжного выноса порошка "
                                  f"({rmin['depowder_ideal_min_mm']['value']} мм) и далеко от правила «≥ 5 мм при длине > 50 мм»; "
                                  "печатаемо только мелким порошком с проверкой проливкой каждого канала (оценка); "
                                  "деревья K ≥ 1 сужают окончания (см. d_min) — ещё хуже"),
        liner_vascular=dict(d_straight_mm=1e3 * al["d_straight"], skin_mm=1e3 * al["t_skin"]),
        recuperator=dict(wall_mm=[1e3 * pr["t_w"], 1e3 * pr["t_w_thin"]], Dh_min_mm=1e3 * pr["Dh_min"],
                         flow_len_m={k: round(d["L"], 3) for k, d in rc.items()},
                         verdict=f"стенка 0,35 мм ниже наименьшего элемента IN718 {rmin['min_feature_in718_mm']['value']} мм "
                                 "по выдержке, 0,25 мм — тем более; плотность тонкой стенки (сквозные поры → перетечка воздуха "
                                 "8 бар в газ 1 бар) не проверена; Dh 1,5 мм при длине пути 0,1–0,26 м — в 30–70 раз длиннее "
                                 "правила выноса порошка (оценка по правилу поставщиков)"),
        rotor_gyroid=dict(sheet_mm=sheet, verdict="лист 0,39–0,52 мм внутри обшивки 2 мм ротора: замкнутая полость, вынос "
                                                  "порошка только через отверстия в самом нагруженном месте — форма снята узлом"),
        fuel_manifold_tree=dict(D_mm=[round(1e3 * x, 1) for x in man["tree"]["D_m"]], wall_band_mm=1e3 * sh["t_band"],
                                eps_rough_um=1e6 * man["assumptions"]["EPS_ROUGH_m"],
                                verdict="ствол Ø42 мм и ветви Ø32 мм не помещаются в стенку 40 мм рамы оболочки — нужен прилив "
                                        "(≥ 48 мм, оценка); горизонтальные каналы > 8 мм — каплевидный профиль; принятые "
                                        "15 мкм — «после обработки», а внутренний канал не обрабатывается (узел сам пишет 20–50 мкм)"),
        shell_channels=dict(radial_mm=[round(1e3 * c[1], 1) for c in sh["radial_channels"]],
                            axial_mm=[round(1e3 * c[1], 1) for c in sh["axial_channels"]],
                            verdict="радиальные каналы 12–22 мм в стойках лежат горизонтально при печати осью вверх — "
                                    "каплевидный профиль; осевые 16–24 мм вертикальны — без поддержек; масляные каналы "
                                    "печати — промывка до класса чистоты масла и контроль частиц обязательны (оценка)"),
    )
    return out


# ------------------------------------------------------------------ 3. шероховатость «сосудов»
def roughness_sensitivity(ks_list=(0.0, 15e-6, 30e-6, 60e-6), T41s=(1050.0, 1100.0), net_names=("straight", "tree_r8", "tree_r4"),
                          augs=("norris", "none")):
    import vascular as VS
    f0 = VS.f_churchill
    f0v = VS.f_churchill_v
    m0 = VS._m_of_dp
    nu0 = VS.nu_pipe
    hy0 = VS.hydraulics
    st = dict(ks=0.0, D=None, aug="norris")

    def fv(Re, rr=0.0):
        Re = np.asarray(Re, float)
        if st["ks"] > 0 and st["D"] is not None and np.shape(st["D"]) == np.shape(Re):
            return f0v(Re, st["ks"] / st["D"])
        return f0v(Re, rr)

    def mdp(dp, rho, mu, D, A, L, Kloc=None):
        st["D"] = np.asarray(D, float)
        return m0(dp, rho, mu, D, A, L, Kloc)

    def hyd(net, *a, **k):
        r = hy0(net, *a, **k)
        st["D"] = r["D"]
        r["f"] = fv(r["Re"])
        r["_d"] = r["D"]
        return r

    def nu(Re, Pr, f):
        """Gnielinski гладкой трубы × (f/f_гл)^n, n = 0,68·Pr^0,215, отношение f не больше 4 (Norris 1971 — по памяти, оценка)"""
        base = nu0(Re, Pr, f)
        if Re <= 2300 or st["ks"] <= 0 or st["aug"] == "none":
            return base
        fs = f0(Re)
        ratio = min(max(f / fs, 1.0), 4.0)
        return base * ratio ** (0.68 * Pr ** 0.215)

    VS.f_churchill_v, VS._m_of_dp, VS.nu_pipe, VS.hydraulics = fv, mdp, nu, hyd
    try:
        pt = VS.kv3_point()
        a = VS.A["ngv"]
        panel0, p_in, p_out, N_v = VS.ngv_case(pt, pt["T41"])
        Lx, H = panel0.Lx, panel0.H
        nets = {"straight": VS.net_straight(Lx, H, a["N_straight"], a["d_straight"])}
        Vt = a["N_straight"] * math.pi * a["d_straight"] ** 2 / 4 * Lx
        for name, nr, ex in (("tree_r8", 8, 3.0), ("tree_r4", 4, 3.0)):
            K = int(round(math.log2(a["N_straight"] / nr)))
            t = VS.net_delta(Lx, H, K, 1e-3, ex, name, n_root=nr)
            t.scale_d(math.sqrt(Vt / t.volume()))
            nets[name] = t
        rows = []
        for T41 in T41s:
            panel, p_in, p_out, _ = VS.ngv_case(pt, T41 + 273.15)
            for name in net_names:
                net = nets[name]
                for ks, aug in [(0.0, "smooth")] + [(k, g) for k in ks_list if k > 0 for g in augs]:
                    st["ks"], st["aug"] = ks, aug
                    t0 = time.time()
                    r = VS.min_flow_for_limit(panel, net, p_out, pt["T3"], a["T_lim_C"] + 273.15, p_in)
                    row = dict(T41_C=T41, net=name, ks_um=round(ks * 1e6, 1), nu_augmentation=aug, feasible=r["feasible"],
                               dp_available_kPa=round((p_in - p_out(Lx)) / 1e3, 1))
                    if r["feasible"]:
                        hydr = r["_hyd"]
                        fr = hydr["f"] / f0v(hydr["Re"])
                        row.update(x_ngv_pct=round(100 * 2 * N_v * r["m_c"] / pt["G"], 4), dp_kPa=round(r["dp"] / 1e3, 1),
                                   Re_max=round(r["Re_max"]), f_ratio_max=round(float(fr.max()), 3),
                                   T_spread_K=round(r["T_spread"], 1), Mach_max=round(r["Mach_max"], 3),
                                   at_backflow_floor=bool(r.get("at_backflow_floor", False)))
                    else:
                        row.update(T_max_full_dp_C=round(r["at_full_dp"]["T_max"] - 273.15, 1))
                    row["seconds"] = round(time.time() - t0, 1)
                    rows.append(row)
                    print("rough", row, flush=True)
        st["ks"] = 0.0
        return dict(rows=rows, N_vanes=N_v, G_kg_s=pt["G"],
                    method="thermal/vascular.py без изменений, подменены: трение Churchill с ε/D по диаметру сегмента; "
                           "теплоотдача — две границы: norris — Nu_Gnielinski,гл·min(f/f_гл, 4)^(0,68·Pr^0,215) (Norris 1971, "
                           "по памяти — оценка, оптимистично для печатных каналов), none — без приращения (осторожно); "
                           "наименьший расход при пределе поверхности 950 °C (min_flow_for_limit)",
                    note="ε — эквивалентная песочная шероховатость; 15 мкм ~ обработанный/полированный канал, 30–60 мкм — "
                         "печать без обработки малых каналов (оценка; у Stimpson 2016 приращение трения в мини-каналах "
                         "больше, чем по ε/D — по выдержке)")
    finally:
        VS.f_churchill_v, VS._m_of_dp, VS.nu_pipe, VS.hydraulics = f0v, m0, nu0, hy0


# ------------------------------------------------------------------ 4. время печати рекуператора
def recup_print_time():
    import make_model as MK
    t = json.loads((ORG / "thermal" / "tpms_results.json").read_text())
    rho = t["recuperator"]["print_rules"]["rho_m"]
    P = {p["name"]: p for p in MK.PRINTERS}
    out = {}
    for case, c in t["recuperator"]["cases"].items():
        out[case] = {}
        for name, d in c["designs"].items():
            vol_cm3 = d["mass_core"] / rho * 1e6
            out[case][name] = dict(metal_m3=round(vol_cm3 / 1e6, 3),
                                   machine_h_BLT_S1500=round(vol_cm3 / P["BLT-S1500"]["rate"]),
                                   machine_h_EP_M2050=round(vol_cm3 / P["Eplus3D EP-M2050"]["rate"]))
        foil = c["duty"].get("cycles_foil")
        if foil:
            out[case]["cycles_foil_mass_kg"] = foil["mass_kg"]
    return dict(rows=out, rates_cm3_h={k: P[k]["rate"] for k in ("BLT-S1500", "Eplus3D EP-M2050")},
                rate_basis={k: P[k]["rate_basis"] for k in ("BLT-S1500", "Eplus3D EP-M2050")},
                note="машинное время = объём металла сердцевины / объёмная скорость камеры (без коллекторов, поддержек, "
                     "охлаждения камеры и ГИП); у тонких стенок скорость ниже паспортной (контур) — числа снизу (оценка)")


# ------------------------------------------------------------------ 5. собираемость оболочки
def shell_assembly():
    import parts_model as PM
    f3 = json.loads((GTU / "frame3" / "frame3_point_results.json").read_text())
    sh = json.loads((ORG / "structure" / "shell_results.json").read_text())
    st1, st2 = f3["stage1"], f3["stage2"]
    tu = f3["turbine"]["options"]["3"]
    Dm = tu["Dm"]
    S = tu["stages"]
    geo = sh["geometry"]
    D_bear = 2 * geo["r_bear_o"]
    rotor = [("колесо 1 (двусторонний вход)", st1["D2_m"]), ("колесо 2", st2["D2_m"])]
    for i, s in enumerate(S):
        rotor.append((f"колесо турбины {i + 1}", Dm + s["h_exit_m"]))
    # статорные элементы между колёсами: отверстие (наименьший диаметр, через который проходит вал/ступица)
    stator = [("вход + передний корпус подшипника", D_bear),
              ("обратный канал к глазу колеса 2", st2["D1_tip_m"]),
              ("внутренний корпус камеры + СА 1", Dm - S[0]["h_rotor_in_m"]),
              ("СА 2", Dm - 0.5 * (S[0]["h_exit_m"] + S[1]["h_rotor_in_m"])),
              ("СА 3", Dm - 0.5 * (S[1]["h_exit_m"] + S[2]["h_rotor_in_m"])),
              ("задний корпус подшипника + выхлоп", D_bear)]
    rows = []
    for k, (nm, b) in enumerate(stator):
        left = rotor[k - 1][1] if k > 0 else None
        right = rotor[k][1] if k < len(rotor) else None
        blocks = [x for x in (left, right) if x is not None and x > b]
        rows.append(dict(stator=nm, bore_m=round(b, 3), D_left_m=None if left is None else round(left, 3),
                         D_right_m=None if right is None else round(right, 3), traps=len(blocks) > 0))
    need = sum(1 for r in rows if r["traps"])            # каждый запирающий статор — отдельная часть
    claimed = 3                                           # оболочка + горячая часть + выхлоп (shell_one_piece)
    ref = PM.count(PM.REF_BOM)
    one = sh["parts"]["shell_one_piece"]
    extra_parts = need - claimed
    per_bolt = math.ceil(math.pi * PM.REF_BOM["case"]["D"] / PM.BOLT_PITCH)
    corr = dict(parts=one["parts"] + extra_parts,
                fasteners_clamp=one["fasteners"] + 2, fasteners_bolts=one["fasteners"] + per_bolt)
    return dict(rotor=[dict(element=n, D_m=round(d, 3)) for n, d in rotor], stator=rows,
                stator_pieces_needed=need, stator_pieces_in_shell_count=claimed,
                parts=dict(reference=ref["parts"], shell_claimed=one["parts"], shell_assemblable=corr["parts"],
                           removed_claimed=ref["parts"] - one["parts"], removed_assemblable=ref["parts"] - corr["parts"]),
                fasteners=dict(reference=ref["fasteners"], shell_claimed=one["fasteners"],
                               shell_assemblable=[corr["fasteners_clamp"], corr["fasteners_bolts"]],
                               removed_assemblable=[ref["fasteners"] - corr["fasteners_bolts"],
                                                    ref["fasteners"] - corr["fasteners_clamp"]]),
                assumption="ротор собирается в корпусе по месту (колесо 1 — спереди, колесо 2 и турбина — сзади, стяжка "
                           "через обратный канал); добавочные части: передняя крышка с корпусом подшипника, кольца СА 2 и СА 3 "
                           "(цельные печатные кольца, вставляются между колёсами); новый наружный стык — один (передняя "
                           "крышка): хомут (2) или болты по правилу parts_model (π·D/0,1 м); D отверстий — frame3 и "
                           "shell_results.geometry (оценка уровня облика)",
                consequence="картридж ротора целиком не вынимается; балансировка ротора в сборе — до разборки и сборка по "
                            "меткам; замена колеса турбины 1 — через выхлоп, колёса 3, 2 и кольца СА 3, 2")


# ------------------------------------------------------------------ 6. цельная горячая часть
def hot_integral_mismatch(T_liner=900.0, T_case=None, r_liner=0.45):
    m = json.loads((GTU / "platform" / "fea" / "materials.json").read_text())["in713c"]
    f3 = json.loads((GTU / "frame3" / "frame3_point_results.json").read_text())
    if T_case is None:
        T_case = f3["cycle"]["frame3"]["T3_C"]
    a_tab = np.array(m["alpha_mean_1e6"], float)
    E_tab = np.array(m["E_GPa"], float)
    al = lambda T: float(np.interp(T, a_tab[:, 0], a_tab[:, 1])) * 1e-6
    eps = al(T_liner) * (T_liner - 20) - al(T_case) * (T_case - 20)
    E = float(np.interp(0.5 * (T_liner + T_case), E_tab[:, 0], E_tab[:, 1])) * 1e3
    return dict(T_liner_C=T_liner, T_case_C=round(T_case, 1), eps_mismatch_pct=round(100 * eps, 3),
                dr_free_mm_at_r=dict(r_m=r_liner, dr_mm=round(1e3 * eps * r_liner, 2)),
                sigma_elastic_biaxial_full_constraint_MPa=round(E * eps / (1 - m["nu"])),
                material_proxy="IN-713C (materials.json, α и E по выдержкам поиска) вместо Hastelloy X / Haynes 282 — "
                               "их свойств в репозитории нет; порядок величины тот же (оценка)",
                verdict="жёстко сросшиеся жаровая труба и корпус должны разойтись на ~3 мм по радиусу: упругое напряжение "
                        "при полном стеснении в разы выше предела текучести любых жаропрочных сплавов при 900 °C — пластическое "
                        "деформирование каждый пуск, малоцикловая усталость; нужна податливая связь (печатные пружинные "
                        "подвески/сильфонный участок) — в счёте деталей и расчётах её нет")


# ------------------------------------------------------------------ 7. ГИП и КТ
def hip_ct():
    sh = json.loads((ORG / "structure" / "shell_results.json").read_text())
    rows = []
    hips = RULES["hip_largest"]["value_m"]
    for r in sh["printer_fit"]["rows"]:
        for k, f in r["fit"].items():
            D = f["D_m"]
            rows.append(dict(frame=r["frame"], kD=k, D_m=D, L_m=f["L_m"],
                             hip_vessels_fit=[h for h in hips if h >= D + 0.05]))
    R, t = sh["geometry"]["R_o"], sh["geometry"]["t_band"]
    chord = 2 * math.sqrt(R ** 2 - (R - t) ** 2)
    return dict(hip=rows, hip_margin_m=0.05,
                ct=dict(shell_R_m=round(R, 3), wall_m=t, tangential_chord_through_wall_m=round(chord, 3),
                        verdict="луч по касательной проходит ~0,46 м металла стенки оболочки (Ti-6Al-4V или никелевый "
                                "сплав): трубкой 450 кВ не просвечивается, нужен ускоритель 9 МэВ с грубым вокселем; каналы "
                                "0,6–1 мм лопаток СА в кольце Ø0,54 м КТ всего кольца не разрешит — контроль проливкой "
                                "каждого канала, свидетели-купоны в той же печати, КТ купонов (оценка)"))


def benefit(res):
    """польза для владельца по каждой форме: детали, масса, КПД — из results-JSON узлов и этого файла"""
    st = json.loads((ORG / "structure" / "structure_results.json").read_text())
    v = json.loads((ORG / "thermal" / "vascular_results.json").read_text())
    tp = json.loads((ORG / "thermal" / "tpms_results.json").read_text())
    man = json.loads((ORG / "structure" / "manifold_results.json").read_text())
    shf = json.loads((ORG / "structure" / "shell_results.json").read_text())
    rc = tp["recuperator"]["cases"]["R90-ic-RH-ibc-none-none"]
    xs = {(r["T41_C"], r["net"]): r.get("x_ngv_pct") for r in v["cases"]["ngv"]["rows"]}
    sa = res["shell_assembly"]
    voids = {k: d["enclosed_voids"] for k, d in res["stl"].items()}
    K, C, F = "оставить", "оставить с условиями", "отсеять: красота без пользы"
    return [
        dict(form="колесо 2 — гладкая задняя поверхность (structure, profile_eq)", parts=0, mass="−14,0 %", J="−22,8 %",
             efficiency="косвенно (пуск)", make="поковка + точение, печать не нужна", decision=K,
             basis="structure_results.json → compressor_wheel_stage2.profile_eq"),
        dict(form="диск турбины 1 — гладкий профиль eq", parts=0, mass="−5,3 %", efficiency=0, make="литьё + точение",
             decision=K, basis="structure_results.json → turbine_disk_stage1.profile_eq"),
        dict(form="диск турбины 1 — «до предела»", parts=0, mass="−21 % (верхняя оценка)", efficiency=0,
             make="литьё + точение", decision="оставить как верхнюю оценку: ресурс 18 929 ч < 25 000 ч на сетке 160×24"),
        dict(form="воксельные проекты A колеса и диска (полости, тороиды)", parts=0, mass="−14…−36 % только в вокселях",
             efficiency=0, make=f"замкнутых полостей в STL колеса: {voids.get('structure/stl/comp2_A_mass_eq.stl')}; "
                                "вынос порошка — отверстиями в нагруженном полотне; в CalculiX загрузка 1,28–1,45",
             decision=F),
        dict(form="ротор с гироидом (вариант B)", parts=0, mass="−0,7 % диск; колесо не проходит", efficiency=0,
             make="замкнутая полость под обшивкой 2 мм, лист 0,39–0,52 мм", decision=F),
        dict(form="корпус-оболочка «одной деталью» + горячая часть", parts=f"−{sa['parts']['removed_assemblable']} "
             f"(заявлено −{sa['parts']['removed_claimed']})", fasteners=f"−{sa['fasteners']['removed_assemblable'][0]}…"
             f"−{sa['fasteners']['removed_assemblable'][1]}", mass="не считалась", efficiency="не считался",
             make="6 частей статора вместо 3; податливая связь жаровой трубы; ГИП только в камерах Ø1,6 / 2,05 м; КТ — "
                  "ускоритель", decision=C),
        dict(form="рама опоры подшипника (оптимизация топологии)", parts=0,
             stiffness=f"×{shf['stiffness_ratio_worst']:.2f} при той же массе", efficiency="потери стоек не считались",
             make="призма, печать осью вверх без нависаний (STL: 0 %)", decision="оставить как второстепенное (жёсткость)"),
        dict(form="топливный коллектор-дерево по Мюррею в стенке", parts=f"−{man['ring']['parts']} деталей, "
             f"−{man['ring']['joints'] - man['tree']['joints']} соединений",
             uniformity=f"{man['ring']['spread_pct']:.1f} % → {man['tree']['spread_pct']:.0f} %", efficiency=0,
             make="ствол Ø42 мм не помещается в стенку 40 мм — прилив; каплевидные каналы", decision=C),
        dict(form="сопловой аппарат: прямые каналы / 8 деревьев у стенки (класс P)", parts="кольцо вместо сегментов",
             coolant=f"{xs.get((1050.0, 'straight'))} % против 3 % класса N при 1050 °C", efficiency="только при T41 > 1000 °C",
             make="каналы 0,6 × 48 мм под обшивкой 0,6 мм; при 1100 °C нужна доводка каналов (roughness_ngv)", decision=C),
        dict(form="сопловой аппарат: 4 дерева (K = 2)", coolant=f"{xs.get((1000.0, 'tree_r4'))} против "
             f"{xs.get((1000.0, 'straight'))} % только при 1000 °C, где KV-3 не охлаждается (класс U); при 1050 °C хуже прямых",
             efficiency=0, decision=F),
        dict(form="сопловой аппарат: деревья K ≥ 3 и растущее дерево (CCO)", efficiency="предел не держат", decision=F),
        dict(form="жаровая труба: каналы у стенки без плёнки", efficiency="предел 900 °C не держит ни одна сеть",
             decision="отсеять как самостоятельное решение; только с плёнкой или покрытием"),
        dict(form="рекуператор: гироид / Schwarz-D", mass=f"ε 0,90: {rc['designs']['gyroid_graded']['mass_core'] / 1e3:.2f} / "
             f"{rc['designs']['diamond_graded']['mass_core'] / 1e3:.2f} т против пластин Dh 1,5 мм "
             f"{rc['designs']['plate_NuH1_Dh15']['mass_core'] / 1e3:.2f} т и фольги {rc['duty']['cycles_foil']['mass_kg'] / 1e3:.2f} т",
             make="стенка 0,35 / 0,25 мм ниже 0,38 мм; вынос порошка из Dh 1,5 мм длиной 0,1–0,26 м не подтверждён",
             decision=F + " (пока нет расчёта с инерцией, Re ≥ 100 — задача AWS)"),
        dict(form="переменная ячейка TPMS (иллюстрация)", mass="оптимизатор выбрал постоянную ячейку", decision=F),
        dict(form="коллектор оптимизацией топологии 2D (comb6, fan8)", parts="1 деталь вместо N труб",
             uniformity="61 % → 3,9 % только без дросселей на выходах", efficiency="для топлива — нет (форсунки выравнивают)",
             make="2D-выдавливание, нормали STL внутрь", decision="отсеять для топлива (там дерево Мюррея); оставить как метод "
                                                               "для раздачи без дросселей (воздух охлаждения) — после 3D"),
    ]


def main():
    t0 = time.time()
    reuse = "--reuse-rough" in sys.argv
    res = dict(title="Скептик изготовимости печатью и пользы: генеративные формы organic/structure и organic/thermal",
               method=__doc__.strip().split("\n\n")[1] if "\n\n" in __doc__ else __doc__.strip(),
               software=dict(python=sys.version.split()[0], numpy=np.__version__),
               date="2026-11-14",
               run_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               status="оценка уровня облика; внешние правила печати, ГИП, КТ, усталость — по выдержкам поиска (RULES)",
               inputs=[dict(file="routes/gtu/organic/structure/stl/*.stl, routes/gtu/organic/thermal/*.stl", use="STL"),
                       dict(file="routes/gtu/organic/thermal/vascular_results.json", field="assumptions.ngv, liner; inputs"),
                       dict(file="routes/gtu/organic/thermal/vascular.py", field="решатель сети и стенки (без правок)"),
                       dict(file="routes/gtu/organic/thermal/tpms_results.json", field="recuperator.print_rules, cases.*.designs"),
                       dict(file="routes/gtu/organic/structure/lattice_results.json", field="sheet_thickness_mm"),
                       dict(file="routes/gtu/organic/structure/manifold_results.json", field="tree.D_m, assumptions"),
                       dict(file="routes/gtu/organic/structure/shell_results.json", field="geometry, parts, printer_fit"),
                       dict(file="routes/gtu/frame3/frame3_point_results.json", field="stage1/2.D2_m, D1_tip_m; turbine.options.3"),
                       dict(file="routes/gtu/funnel/parts_model.py", field="count, REF_BOM, BOLT_PITCH"),
                       dict(file="routes/gtu/funnel/make_model.py", field="PRINTERS (rate)"),
                       dict(file="routes/gtu/platform/fea/materials.json", field="in713c.alpha_mean_1e6, E_GPa, nu")],
               rules=RULES)
    res["stl"] = {}
    for rel, kind, what in STLS:
        p = ORG / rel
        if p.exists():
            d = stl_audit(p, kind)
            d.update(kind=kind, what=what)
            res["stl"][rel] = d
            print("stl", rel, d["overhang_frac_best"], d["best_build_dir"], d["enclosed_voids"], flush=True)
    res["features"] = features()
    res["recup_print_time"] = recup_print_time()
    res["shell_assembly"] = shell_assembly()
    res["hot_integral"] = hot_integral_mismatch()
    res["hip_ct"] = hip_ct()
    prev = HERE / "print_review_results.json"
    if reuse and prev.exists() and "roughness_ngv" in json.loads(prev.read_text()):
        res["roughness_ngv"] = json.loads(prev.read_text())["roughness_ngv"]
        res["roughness_ngv"]["reused_from_previous_run"] = True
    else:
        res["roughness_ngv"] = roughness_sensitivity()
    res["benefit"] = benefit(res)
    res["runtime_s"] = round(time.time() - t0, 1)
    (HERE / "print_review_results.json").write_text(json.dumps(res, ensure_ascii=False, indent=1))
    print("записано", HERE / "print_review_results.json", res["runtime_s"], "с")


if __name__ == "__main__":
    main()
