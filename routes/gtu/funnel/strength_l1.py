"""Воронка двигателя-заменителя, уровень 1: быстрая прочность ротора. Всё - оценка.

  * Колёса компрессора и радиальной турбины (заодно с лопатками) - подобие расчётам CalculiX ../platform/fea/disk_results.json
    (тот же условный профиль, масштаб по D): напряжения ∝ ρ·U², запас по разрушающей частоте (критерий средних окружных
    напряжений, k = 1 - как disk.py, верхняя оценка) ∝ (U_пр/U)·√(σв(T)/σв_пр(T_пр)·ρ_пр/ρ). Колесо сборного ротора - с
    центральным отверстием под стяжку (пример comp_T60_1C_ti), цельное (заодно, сварка) - без отверстия (comp_T60_1C_ti_solid).
  * Диск осевой турбины (сборный ротор: диск + лопатки) - аналитика диска постоянной толщины с нагрузкой на ободе от лопаток:
    σθ,ср = ρω²R²/3 + p_обода (равновесие полудиска); блиск - так же с лопатками в массе обода.
  * Лопатка осевой турбины: напряжение в корне σ = K_T·ρ·AN²·π/1800, K_T подобран по ступени frame3 (469 МПа при AN² 6,028e7,
    IN-713C - frame3_fea_results.json, как правило AN² уровня 0).
  * Ползучесть - Ларсон - Миллер (C = 20) по блоку lmp IN-713C в materials.json (единственный сплав с данными); K418 -
    близкий аналог IN713C (make_model); прочие сплавы - та же кривая со сдвигом температуры dT_cap (make_model; оценка).
    Металл: радиальное колесо - обод при T41 - U²/(2cp) (как пример turb130), осевая лопатка - относительная температура
    торможения на входе ротора; охлаждаемое колесо - минус DT_COOL. Требование RUPTURE_H - оценка.
  * Керамика (Si3N4): вероятность разрушения по Вейбуллу P = 1 - exp(-(V_eff/V0)·(σ/σ0)^m): средняя прочность 594 МПа при
    1200 °C - Kyocera SN282A (по выдержке поиска, global.kyocera.com/prdct/fc/technologies/021.html); модуль Вейбулла m в
    источнике не найден - WEIBULL_M оценка; σ0 = σ_ср/Γ(1 + 1/m); V0 и доля объёма колеса под пиком - оценка.
  * CalculiX для лучших: run_disk_fea() строит вход ../platform/fea/disk.py (условный профиль примера в масштабе D, материалы
    materials.json) и считает осесимметричный диск (ccx) - Мизес, σθ,ср, запас по разрушающей частоте, ресурс (IN-713C).
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
FEA = os.path.join(GTU, "platform", "fea")
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import make_model as MK          # noqa: E402

K_BURST = 1.0                 # как disk.py: верхняя оценка (настоящий k < 1 - разгонные испытания)
OVERSPEED = 1.20              # запас по разрушающей частоте проверяется при 120 % частоты (задача владельца)
RUPTURE_H = 25000.0           # требуемый ресурс по длительной прочности, ч (оценка: межремонтный интервал промышленных ГТУ
#                               порядка 25-30 тыс. ч - не проверено)
DT_COOL = 175.0               # снижение температуры металла охлаждаемого колеса / лопатки, K (оценка: класс K допускает
#                               T41 1200 против 1000-1050 °C неохлаждаемого - make_model.t41_limit)
WEIBULL_M = 15.0              # модуль Вейбулла Si3N4 - оценка (в выдержке источника нет; принято 15 - середина 10-20)
SN_MEAN_MPA = 594.0           # средняя прочность SN282A при 1200 °C (Kyocera, по выдержке поиска)
WEIBULL_V0_MM3 = 1.0          # эффективный объём образца, мм³ (оценка: изгиб образца B)
WEIBULL_VFRAC = 0.02          # доля объёма колеса под пиковым напряжением (оценка)
PF_MAX = 1e-4                 # допустимая вероятность разрушения керамического колеса (оценка)
_D = {}


def _data():
    if not _D:
        mats = json.load(open(os.path.join(FEA, "materials.json"), encoding="utf-8"))
        cases = {c["case"]["name"]: c for c in json.load(open(os.path.join(FEA, "disk_results.json"), encoding="utf-8"))["cases"]}
        f3 = json.load(open(os.path.join(GTU, "frame3", "frame3_fea_results.json"), encoding="utf-8"))
        _D.update(mats=mats, cases=cases, f3=f3, lmp=_lmp_fit(mats["in713c"]["lmp"]))
    return _D


def _interp(tab, x):
    xs, ys = [t[0] for t in tab], [t[1] for t in tab]
    if x <= xs[0]:
        return ys[0], x >= xs[0] - 1e-9
    if x >= xs[-1]:
        return ys[-1], False
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            return ys[i] + (ys[i + 1] - ys[i]) * (x - xs[i]) / (xs[i + 1] - xs[i]), True
    return ys[-1], False


def _lmp_fit(lmp):
    """lg σ = a + b·P + c·P², P = T[K]·(C + lg t) - МНК по точкам materials.json (как disk.py)"""
    C = lmp["C"]
    P = [(T + 273.15) * (C + math.log10(t)) for T, t, s, _ in lmp["data"]]
    y = [math.log10(s) for T, t, s, _ in lmp["data"]]
    n = len(P)
    S = [[sum(p ** (i + j) for p in P) for j in range(3)] for i in range(3)]
    r = [sum(yy * p ** i for p, yy in zip(P, y)) for i in range(3)]
    # решение 3x3 (Крамер)
    def det(M):
        return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    D = det(S)
    coef = []
    for k in range(3):
        M = [row[:] for row in S]
        for i in range(3):
            M[i][k] = r[i]
        coef.append(det(M) / D)
    return dict(C=C, a=coef[0], b=coef[1], c=coef[2], P_range=(min(P), max(P)), n=n)


def rupture_hours(sigma_MPa, T_C, dT_cap=0.0):
    """ресурс по длительной прочности (ч) кривой IN-713C при σ и температуре металла T - dT_cap (сдвиг сплава);
    возврат: (часы или None, экстраполяция за данные?)"""
    if sigma_MPa <= 0:
        return math.inf, False
    L = _data()["lmp"]
    y = math.log10(sigma_MPa)
    a, b, c = L["a"] - y, L["b"], L["c"]
    disc = b * b - 4 * c * a
    if disc < 0:
        return None, True
    roots = [(-b + s * math.sqrt(disc)) / (2 * c) for s in (1, -1)] if abs(c) > 1e-30 else [-a / b]
    # убывающая ветвь lg σ(P): производная b + 2cP < 0
    P = next((p for p in roots if b + 2 * c * p < 0), max(roots))
    T = T_C - dT_cap + 273.15
    t = 10 ** (P / T - L["C"])
    extrap = not (L["P_range"][0] * 0.97 <= P <= L["P_range"][1] * 1.03)
    return t, extrap


def creep_material(mat):
    """кривая IN-713C и сдвиг температуры для сплава колеса: (сдвиг K, основание)"""
    if mat == "K418":
        return 0.0, "K418 - близкий аналог IN713C (make_model): кривая IN-713C"
    if mat == "SC":
        return 75.0, "монокристалл: кривая IN-713C со сдвигом +75 K (оценка; в make_model dT_cap нет)"
    d = MK.MAT.get(mat, {}).get("dT_cap")
    if d is None:
        return 0.0, f"{mat}: данных нет - кривая IN-713C без сдвига (оценка)"
    return float(d), f"{mat}: кривая IN-713C со сдвигом dT_cap {d:+.0f} K (make_model, оценка)"


def _uts(mat_key, T):
    return _interp(_data()["mats"][mat_key]["uts_MPa"], T)


def _ys(mat_key, T):
    return _interp(_data()["mats"][mat_key]["ys_MPa"], T)


def wheel(kind, mat, rho, U, T_mean_C, solid=True):
    """колесо заодно с лопатками по подобию примеру disk_results.json. kind: comp | rturb. mat - ключ materials.json или
    имя сплава (никелевые - свойства IN-713C, оценка). Возврат: vm, σθ,ср, запас разрушения (100 % и 120 %), отношение к σт"""
    D = _data()
    if kind == "comp":
        ex = D["cases"]["comp_T60_1C_ti_solid" if solid else "comp_T60_1C_ti"]
        ex_mat, ex_rho = "ti64", 4430.0
    else:
        ex = D["cases"]["turb130_in713c"]
        ex_mat, ex_rho = "in713c", 7916.0
    r = ex["result"]
    mk = "ti64" if mat in ("Ti64", "ti64") else ("in713c" if mat not in ("Al", "Si3N4") else None)
    k = (rho / ex_rho) * (U / r["U_tip"]) ** 2
    out = dict(example=ex["case"]["name"], scale_stress=k, vm_MPa=r["vm_max_node_MPa"] * k,
               hoop_mean_MPa=r["section"]["hoop_mean_ip_MPa"] * k, U=U, rho=rho, T_mean_C=T_mean_C, material_data=mk)
    if mk is None:
        out.update(burst_margin=None, burst_margin_120=None, yield_ratio=None,
                   note=f"{mat}: нет свойств в materials.json - запас не считается")
        return out
    uts_ex, _ = _uts(ex_mat, r["section"]["T_mean_C"])
    uts, inr = _uts(mk, T_mean_C)
    m100 = r["burst"]["margin"] * math.sqrt(K_BURST) * (r["U_tip"] / U) * math.sqrt((uts / uts_ex) * (ex_rho / rho))
    out.update(burst_margin=m100, burst_margin_120=m100 / OVERSPEED, uts_MPa=uts, uts_in_data=inr)
    if "ys_MPa" in D["mats"][mk] and len(D["mats"][mk]["ys_MPa"]) > 1:
        ys, yin = _ys(mk, T_mean_C)
        out.update(yield_ratio=out["vm_MPa"] / ys, ys_MPa=ys, ys_in_data=yin)
    else:
        out.update(yield_ratio=None)
    return out


def sigma_allow(T_C, hours=None, dT_cap=0.0):
    """допускаемое напряжение по длительной прочности (МПа) кривой IN-713C за hours при T - dT_cap; (σ, экстраполяция?)"""
    L = _data()["lmp"]
    hours = hours or RUPTURE_H
    P = (T_C - dT_cap + 273.15) * (L["C"] + math.log10(hours))
    extrap = not (L["P_range"][0] * 0.97 <= P <= L["P_range"][1] * 1.03)
    return 10 ** (L["a"] + L["b"] * P + L["c"] * P * P), extrap


DISC_RATIO_MAX = 6.0          # толщина ступицы / толщина обода равнопрочного диска (оценка: технологичность поковки)


def disc_axial(rho, n, Dm, H, m_blades, T_rim_C, mat="IN718", dT_cap=0.0):
    """диск осевой ступени, спроектированный равнопрочным (Стодола, как frame3): расчётное напряжение
    σ_d = min(k·σв(T)/1,2², длительная прочность за RUPTURE_H) - тогда запас по разрушению при 120 % >= 1 по построению;
    проверяется, можно ли такой диск сделать: толщина обода под нагрузку лопаток t_об = F/(2πRσ_d) и отношение толщин
    ступицы и обода exp(ρω²R²/(2σ_d)) (<= DISC_RATIO_MAX). Свойства - IN-713C (IN718 в materials.json нет; оценка)"""
    om = 2 * math.pi * n / 60
    R = max(0.05, 0.5 * (Dm - H) - 0.03 * Dm)               # живой обод под замками (оценка; frame3 - 0,156 при Dm 0,537)
    F = m_blades * om ** 2 * (0.5 * Dm)
    uts, inr = _uts("in713c", T_rim_C)
    s_b = K_BURST * uts / OVERSPEED ** 2
    s_c, ext = sigma_allow(T_rim_C, dT_cap=dT_cap)
    s_d = min(s_b, s_c)
    t_rim = F / (2 * math.pi * R * s_d * 1e6)
    ratio = math.exp(rho * om ** 2 * R ** 2 / (2 * s_d * 1e6))
    return dict(R=R, sigma_d_MPa=s_d, sigma_burst_MPa=s_b, sigma_creep_MPa=s_c, creep_extrapolated=ext, governed_by=(
        "разрушение при 120 %" if s_b <= s_c else "длительная прочность"), t_rim=t_rim, t_hub=t_rim * ratio, ratio=ratio,
        burst_margin_120=math.sqrt(K_BURST * uts / s_d) / OVERSPEED, uts_MPa=uts, uts_in_data=inr, U_rim=om * R, material=mat)


def k_taper():
    """K_T корня лопатки по ступени frame3 (turb3): σ = K_T·ρ·AN²·π/1800"""
    f3 = _data()["f3"]
    p = f3["parts"]["turb3"]["geometry"]
    an2 = float(p["AN2_point_m2rpm2"])
    sig = float(p["sigma_root_MPa"])
    rho = float(_data()["mats"]["in713c"]["density"])
    return sig * 1e6 / (rho * an2 * math.pi / 1800), dict(AN2=an2, sigma_root_MPa=sig, rho=rho)


def blade_root(rho, AN2):
    k, _ = k_taper()
    return k * rho * AN2 * math.pi / 1800 / 1e6


def weibull_pf(sigma_MPa, D_wheel, m=None, mean_MPa=None):
    """вероятность разрушения керамического колеса: V_eff = WEIBULL_VFRAC·V колеса (V ≈ 0,0825·D³, как масса уровня 0)"""
    m = m or WEIBULL_M
    mean = mean_MPa or SN_MEAN_MPA
    s0 = mean / math.gamma(1 + 1 / m)
    V = 0.0825 * D_wheel ** 3 * 1e9 * WEIBULL_VFRAC
    pf = 1 - math.exp(-(V / WEIBULL_V0_MM3) * (max(sigma_MPa, 0.0) / s0) ** m)
    return dict(pf=pf, sigma_MPa=sigma_MPa, sigma0_MPa=s0, m=m, V_eff_mm3=V,
                pf_m10=1 - math.exp(-(V / WEIBULL_V0_MM3) * (max(sigma_MPa, 0.0) / (mean / math.gamma(1.1))) ** 10),
                pf_m20=1 - math.exp(-(V / WEIBULL_V0_MM3) * (max(sigma_MPa, 0.0) / (mean / math.gamma(1.05))) ** 20))


def facts():
    k, src = k_taper()
    L = _data()["lmp"]
    return dict(K_BURST=K_BURST, OVERSPEED=OVERSPEED, RUPTURE_H=RUPTURE_H, DT_COOL=DT_COOL, WEIBULL_M=WEIBULL_M,
                SN_MEAN_MPA=SN_MEAN_MPA, WEIBULL_V0_MM3=WEIBULL_V0_MM3, WEIBULL_VFRAC=WEIBULL_VFRAC, PF_MAX=PF_MAX,
                DISC_RATIO_MAX=DISC_RATIO_MAX, k_taper=round(k, 4), k_taper_source=src,
                lmp_fit={x: (float(f"{v:.8g}") if isinstance(v, float) else v) for x, v in L.items()},
                sources=dict(
                    disk="../platform/fea/disk_results.json (CalculiX, условные профили; масштаб по ρU²)",
                    lmp="../platform/fea/materials.json, in713c.lmp (по выдержкам поиска - см. статусы файла)",
                    sn282a="Kyocera SN282A: 594 МПа при 1200 °C - по выдержке поиска (global.kyocera.com/prdct/fc/technologies/021.html); "
                           "модуль Вейбулла - не найден, оценка",
                    frame3="../frame3/frame3_fea_results.json, parts.turb3 - корень лопатки"))


# =========================================================================== CalculiX: disk.py для лучших
def disk_case(name, kind, mat_key, rpm, D, T_in_C, T_out_C, solid=True):
    """вход disk.py: условный профиль примера (comp_T60_1C_ti / turb130_in713c) в масштабе D - как disk.example_cases"""
    if FEA not in sys.path:
        sys.path.insert(0, FEA)
    import disk as DK                                           # noqa: E402
    ex = {c["name"]: c for c in DK.example_cases()}
    base = ex["comp_T60_1C_ti_solid" if (kind == "comp" and solid) else ("comp_T60_1C_ti" if kind == "comp" else "turb130_in713c")]
    R0 = base["profile"]["r"][-1]
    k = (D / 2) / R0
    prof = {key: [x * k for x in base["profile"][key]] for key in ("r", "t", "hb")}
    rin = prof["r"][0]
    return dict(name=name, material=mat_key, rpm=float(rpm), profile=prof, profile_interp="pchip", blade_fill=base["blade_fill"],
                T_C=DK._quad_T(rin, D / 2, T_in_C, T_out_C), burst_k=1.0, mesh=dict(nr=40, nz=6),
                note=f"уровень 1: условный профиль {base['name']} в масштабе {k:.3f} (оценка)")


def run_disk_fea(case, workdir):
    if FEA not in sys.path:
        sys.path.insert(0, FEA)
    import disk as DK                                           # noqa: E402
    r = DK.run_case(case, workdir=workdir, keep=False)
    life = r.get("creep") or {}
    return dict(vm_max_MPa=r["vm_max_node_MPa"], hoop_mean_MPa=r["section"]["hoop_mean_ip_MPa"], burst_margin=r["burst"]["margin"],
                burst_margin_120=r["burst"]["margin"] / OVERSPEED, yield_ratio=r.get("yield_ratio_max"),
                life_local_h=life.get("local_min_h"), life_s1_h=life.get("local_s1_min_h"), equilibrium_ok=r["quality"]["passed"],
                seconds=r["run"]["seconds"], flags=r.get("flags"))


def run_cases(path, out_dir, workers=1):
    """задача AWS «CalculiX для лучших 300»: cases.json (merge_l1 --fea-cases-out) -> результаты disk.py по каждому случаю"""
    import tempfile
    from concurrent.futures import ThreadPoolExecutor
    doc = json.load(open(path, encoding="utf-8"))
    tmp = tempfile.mkdtemp(prefix="l1_fea_")

    def one(item):
        try:
            return dict(key=item["key"], result=run_disk_fea(item["case"], os.path.join(tmp, item["case"]["name"])))
        except Exception as e:                       # noqa: BLE001 - отказ ccx пишется, задача не рушится
            return dict(key=item["key"], error=f"{type(e).__name__}: {e}")
    with ThreadPoolExecutor(max(1, workers)) as ex:
        got = list(ex.map(one, doc["cases"]))
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "funnel_l1_fea_results.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(dict(method="disk.py (CalculiX ccx, осесимметричный диск, условный профиль в масштабе D)", date="2026-11-14",
                       input=os.path.basename(path), n=len(got), results=got), fh, ensure_ascii=False, indent=1)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="CalculiX для лучших вариантов уровня 1 (задача AWS)")
    ap.add_argument("--run-cases", required=True)
    ap.add_argument("--out", default=os.environ.get("GTU_OUT", "."))
    ap.add_argument("--workers", type=int, default=int(os.environ.get("NPROC", "1") or 1))
    a = ap.parse_args()
    print(run_cases(a.run_cases, a.out, a.workers))
