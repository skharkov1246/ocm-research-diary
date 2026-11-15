"""Массовая 3D-прочность: запись варианта уровня 1 (поля concept, geometry, meanline, strength) -> роторные детали для
трёхмерного расчёта, их нагрузки и источники возбуждения.

Деталь (dict): id (comp1, turb1, pt1 ...), kind (wheel | axial), role (comp | turb | pt), geom (спецификация для geom3d),
rpm (расчётная частота детали), material (сплав), материал диска вставных лопаток, поле температур (T_spec), источники
возбуждения (eo_sources) и рабочий диапазон частоты вращения (op_range), неподдерживаемое - в pending с причиной.

Правила (все - оценка, причины - в тексте):
  * колесо компрессора - центробежная ступень: impeller.py по D2, b2, D1t, D1h, β2л, Z (+ полулопатки), двусторонний вход;
    осевые и диагональные ступени - pending (нет генератора, у осевой ступени в записи нет одномерной модели);
  * радиальная турбина - тот же генератор: D2 <- D4, b2 <- b4, вход колеса <- эксдьюсер (D5t, D5h), β2л = 0 (радиальные
    лопатки на входе), углы лопатки эксдьюсера - по треугольнику выхода (cm5/U4), осевая длина Lz одномерной модели;
  * осевая турбина и СТ - ступени meanline: Dm, H3 (высота на выходе ротора - как AN² уровня 1), хорда и осевая хорда ротора,
    углы β2, β3 и α2 на среднем радиусе -> свободный вихрь по высоте; диск - равнопрочный профиль strength_l1.disc_axial
    (t_об, σ_d) с ободом под лопаткой этой ступени; блиск - по правилу funnel_space._blisk (литьё до Dtip 0,50 м, печать - по
    камере принтера), иначе - вставные лопатки и диск IN718 (как уровень 1);
  * цельное колесо (ротор заодно, сварка) - технологическое отверстие (geom3d.A core_bore_r1h), сборный ротор - отверстие
    под стяжку;
  * поле температур (уровень 1, по радиусу): колесо компрессора - T_вх у отверстия -> T_вых на r2 по квадрату радиуса (как
    disk.py); радиальная турбина - обод T_rim_C (уровень 1), центр на 240 K холоднее (как merge_l1.fea); осевая ступень -
    перо и полка при температуре металла лопатки (относительная температура торможения на входе ротора минус охлаждение -
    как strength_l1), обод диска на 150 K холоднее (strength_l1.disc_axial), ступица ещё на 200 K (frame3_fea, funnel_l2);
  * источники возбуждения: число лопаток СА × 1 и × 2 (ступень турбины и СТ), число горелок (1-я ступень турбины газогенератора),
    поворотный ВНА × 1 и × 2 (1-я ступень компрессора), низкие порядки 1..LOW_EO; число лопаток диффузора в записи нет -
    для колёс компрессора считается, какие числа лопаток диффузора 9..40 дают резонанс (выбор - за конструктором);
  * рабочий диапазон: одновальная машина (генератор) - 0,95..1,05; газогенератор двухвальной - 0,80..1,02; СТ генераторного
    исполнения - 0,95..1,05 (диапазон приводного исполнения 0,70..1,05 - справочно).
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
FUNNEL = os.path.join(GTU, "funnel")

LOW_EO = 4                         # «низкие порядки»: 1..4 (дисбаланс, овальность, неравномерность входа; оценка)
OP_RANGE = dict(single=(0.95, 1.05), gg=(0.80, 1.02), pt=(0.95, 1.05), pt_drive=(0.70, 1.05))
OP_BASIS = ("одновальная (генератор): 0,95..1,05 n; газогенератор двухвальной: 0,80..1,02 n (частичная нагрузка - уровень 1 "
            "считает 50-100 %); СТ генераторного исполнения KV-3: 0,95..1,05; приводное исполнение 0,70..1,05 - справочно (оценка)")
DIFF_Z_RANGE = (9, 40)             # числа лопаток диффузора, проверяемые на резонанс (оценка)
T_HUB_DROP_RAD = 240.0             # радиальная турбина: центр холоднее обода (как merge_l1.fea, funnel_l2)
T_DISC_DROP = 150.0                # обод диска осевой ступени холоднее лопатки (strength_l1.disc_axial: Tb - 150)
T_HUB_DROP_AX = 200.0              # ступица диска осевой ступени холоднее обода (frame3_fea, funnel_l2.T_HUB_DROP)
T_HUB_RAD_T = dict(t_hub=0.014, t_tip=0.005)   # толщина лопатки радиальной турбины у втулки / на периферии, доля D4 (оценка)
DISC_RATIO_MAX = 6.0               # как strength_l1.DISC_RATIO_MAX: профиль диска ограничен по толщине ступицы

_FS = None


def _blisk(m_tw, mat, Dtip, Dm):
    """правило funnel_space._blisk (импорт; запасной путь - то же правило по make_model)"""
    global _FS
    try:
        if _FS is None:
            if FUNNEL not in sys.path:
                sys.path.insert(0, FUNNEL)
            import funnel_space as FS     # noqa: E402
            _FS = FS
        return bool(_FS._blisk(m_tw, mat, Dtip, Dm)), "funnel_space._blisk"
    except Exception as e:                 # noqa: BLE001
        if mat == "SC":
            return False, f"запасное правило ({type(e).__name__})"
        if m_tw == "печать":
            return True, f"запасное правило ({type(e).__name__})"
        return bool(m_tw == "литьё" and Dtip <= 0.50), f"запасное правило ({type(e).__name__})"


def _stage_inlets(ml):
    """полные T и p на входе ступеней компрессора (как funnel_l2._stage_inlets)"""
    T, p, out = 288.15, 101325.0, []
    for st in ml:
        out.append((T, p))
        pr, ep = float(st["PRs"]), float(st.get("eta_p") or 0.85)
        T = T * pr ** (0.4 / (1.4 * ep))
        p = p * pr
    return out


def T_quad(r_in, r_out, T_in, T_out):
    return dict(type="quad", r_in=r_in, r_out=r_out, T_in=T_in, T_out=T_out)


def comp_parts(v):
    out, pending = [], []
    c, geo, ml, stg = v["concept"], v["geometry"]["compressor"], v["meanline"]["compressor"], v["strength"]["comp"]
    cores = int(c.get("cores") or 1)
    G = float(v["G"]) / cores
    inl = _stage_inlets(ml)
    solid = c.get("rotor_join") in ("заодно", "сварка")
    vigv = (v["geometry"].get("inlet") or {}).get("vigv")
    k_ml = 0
    for i, g in enumerate(geo):
        sid = f"comp{i + 1}"
        if g.get("type") != "центробежная":
            pending.append(dict(id=sid, what=f"колесо компрессора ({g.get('type')})",
                                why="генератор 3D - только центробежное колесо; у осевой ступени в записи нет одномерной модели"))
            continue
        st = ml[k_ml]
        w = next((x for x in stg if x.get("stage") == i + 1), stg[k_ml] if k_ml < len(stg) else {})
        k_ml += 1
        mat = w.get("mat") or c.get("mat_cw")
        alloy = c.get("mat_tw") if mat == "Ni" else mat
        T01, p01 = inl[min(i, len(inl) - 1)]
        geom = dict(name=sid, m=round(G, 5), PR=float(st["PRs"]), rpm=float(v["n_rpm"]), D2=float(st["D2"]),
                    Z=int(round(float(st.get("Zfull") or g.get("blades")))), splitters=bool(g.get("splitters")),
                    b2=float(st["b2"]), beta2b=float(st["beta2b_deg"]), D1s=float(st["D1t"]), nu=float(st["D1h"]) / float(st["D1t"]),
                    gap=float(st.get("cl") or 3e-4), double_sided=bool(st.get("dbl")), T01=T01, p01=p01,
                    bore="core" if solid else "tie")
        T_in = float(w.get("T_in_C") if w.get("T_in_C") is not None else T01 - 273.15)
        T_out = float(w.get("T_out_C") if w.get("T_out_C") is not None else T_in + 120.0)
        eo = [dict(eo=k, kind="low", source=f"низкий порядок {k}") for k in range(1, LOW_EO + 1)]
        if i == 0 and vigv:
            zi = int(vigv.get("vanes_per_inlet") or 0)
            if zi:
                eo += [dict(eo=zi, kind="igv", source=f"ВНА {zi} лопаток × 1"), dict(eo=2 * zi, kind="igv2", source=f"ВНА {zi} лопаток × 2")]
        out.append(dict(id=sid, kind="wheel", role="comp", stage=i + 1, geom=geom, rpm=float(v["n_rpm"]), alloy=alloy,
                        T_spec=T_quad(None, geom["D2"] / 2, T_in, T_out), eo_sources=eo, free_stator=list(DIFF_Z_RANGE),
                        op=_op(v, "comp"), solid=solid))
    return out, pending


def _op(v, role):
    rotors = int(v["concept"].get("rotors") or 1)
    if role == "pt":
        return dict(range=OP_RANGE["pt"], info_range=OP_RANGE["pt_drive"], basis=OP_BASIS)
    if rotors == 1:
        return dict(range=OP_RANGE["single"], info_range=None, basis=OP_BASIS)
    return dict(range=OP_RANGE["gg"], info_range=None, basis=OP_BASIS)


def radial_turbine_parts(v):
    out, pending = [], []
    c, tm, t = v["concept"], v["meanline"]["turbine"], v["strength"].get("turb") or {}
    cores = int(c.get("cores") or 1)
    G = float(v["G"]) / cores
    solid = c.get("rotor_join") in ("заодно", "сварка")
    zn, zt = int(c.get("zn") or 0), int(c.get("zt") or 0)
    ninj = int((v["geometry"].get("combustor") or {}).get("injectors") or 0)
    T_rim = float(t.get("T_rim_C")) if t.get("T_rim_C") is not None else None
    for i, st in enumerate(tm.get("stages") or []):
        sid = f"turb{i + 1}"
        D4, D5t, D5h = float(st["D4"]), float(st["D5t"]), float(st["D5h"])
        cm_u = float(st.get("cm5_U4") or 0.5)
        geom = dict(name=sid, m=round(G, 5), PR=2.0, rpm=float(v["n_rpm"]), D2=D4, Z=zt, splitters=False, b2=float(st["b4"]),
                    beta2b=0.0, D1s=D5t, nu=D5h / D5t, Lz=float(st["Lz"]),
                    beta1s=math.degrees(math.atan((D5t / D4) / cm_u)), beta1h=math.degrees(math.atan((D5h / D4) / cm_u)),
                    T01=float(st.get("T05") or 1000.0), p01=float(st.get("p05") or 2e5), double_sided=False,
                    bore="core" if solid else "tie", **T_HUB_RAD_T)
        Tr = T_rim if (i == 0 and T_rim is not None) else (float(st.get("T_rel_in", 1100.0)) - 273.15
                                                           - (175.0 if T_rim is not None and i > 0 and _cooled(v) else 0.0))
        eo = [dict(eo=k, kind="low", source=f"низкий порядок {k}") for k in range(1, LOW_EO + 1)]
        if zn:
            eo += [dict(eo=zn, kind="ngv", source=f"СА {zn} лопаток × 1"), dict(eo=2 * zn, kind="ngv2", source=f"СА {zn} лопаток × 2")]
        if i == 0 and ninj:
            eo.append(dict(eo=ninj, kind="inj", source=f"горелки камеры: {ninj}"))
        out.append(dict(id=sid, kind="wheel", role="turb", stage=i + 1, geom=geom, rpm=float(v["n_rpm"]), alloy=c.get("mat_tw"),
                        T_spec=T_quad(None, D4 / 2, Tr - T_HUB_DROP_RAD, Tr), eo_sources=eo, free_stator=None, op=_op(v, "turb"),
                        solid=solid))
    return out, pending


def _cooled(v):
    t = (v["strength"].get("turb") or {})
    s = (v["meanline"].get("turbine") or {}).get("stages") or []
    if t.get("T_blade_C") is not None and s and s[0].get("T_rel_in"):
        return float(s[0]["T_rel_in"]) - 273.15 - float(t["T_blade_C"]) > 50.0
    if t.get("T_rim_C") is not None:
        return False
    return False


def _free_vortex_sections(st, n_rpm, r_h, r_t):
    """углы пера на корне, середине и периферии по свободному вихрю: cθ2·r = const, на выходе без закрутки (ax_turbine)"""
    Dm = float(st["Dm"])
    rm = Dm / 2
    Um = math.pi * Dm * n_rpm / 60.0
    cx = float(st["phi"]) * Um
    ct2m = cx * math.tan(math.radians(float(st["alpha2_deg"])))
    secs = []
    bx = float(st["bx_R"])
    for r, tc in ((r_h, None), (rm, None), (r_t, None)):
        U = Um * r / rm
        b_in = math.degrees(math.atan((ct2m * rm / r - U) / cx))
        b_out = math.degrees(math.atan(U / cx))
        f = (r - r_h) / max(r_t - r_h, 1e-9)
        tmc = 0.20 + (0.08 - 0.20) * f
        secs.append(dict(r=r, beta_in=b_in, beta_out=b_out, c_ax=bx, t_max_c=tmc))
    return secs


def _disc_profile(disc, rpm, rho, R):
    """равнопрочный диск (strength_l1.disc_axial): t(r) = t_об·exp(ρω²(R² - r²)/(2σ_d)), отверстие 0,2 R, толщина ступицы
    не больше DISC_RATIO_MAX·t_об"""
    om = 2 * math.pi * rpm / 60
    s_d = float(disc["sigma_d_MPa"]) * 1e6
    t_rim = float(disc["t_rim"])
    rs = [R * (0.20 + 0.80 * k / 12) for k in range(13)]
    ts = [min(DISC_RATIO_MAX * t_rim, t_rim * math.exp(rho * om ** 2 * (R ** 2 - r ** 2) / (2 * s_d))) for r in rs]
    return rs, ts


def axial_parts(v, role):
    """ступени осевой турбины газогенератора (role turb) или СТ (role pt)"""
    out, pending = [], []
    c = v["concept"]
    key = "turbine" if role == "turb" else "power_turbine"
    tm = v["meanline"].get(key)
    sk = "turb" if role == "turb" else "pt"
    t = v["strength"].get(sk) or {}
    if not tm or not t:
        return out, pending
    zn, zt = int(c.get("zn") or 0), int(c.get("zt") or 0)
    ninj = int((v["geometry"].get("combustor") or {}).get("injectors") or 0)
    if role == "turb":
        rpm = float(v["n_rpm"])
        m_tw, alloy = c.get("m_tw"), c.get("mat_tw")
        cooled = _cooled(v)
    else:
        rpm = float((v["geometry"].get("power_turbine") or {}).get("n_rpm") or v["n_rpm"])
        m_tw = c.get("m_tw") if c.get("m_tw") != "керамика" else "литьё"
        alloy = c.get("mat_tw") if m_tw == "печать" else "K418"
        cooled = False
    disc = t.get("disc") or {}
    for i, st in enumerate(tm.get("stages") or []):
        sid = f"{sk}{i + 1}"
        Dm, H3 = float(st["Dm"]), float(st["H3"])
        r_h, r_t = 0.5 * (Dm - H3), 0.5 * (Dm + H3)
        blisk, basis = _blisk(m_tw, alloy, r_t * 2, Dm)
        R = max(0.05, 0.5 * (Dm - H3) - 0.03 * Dm)
        if R >= r_h:
            R = 0.92 * r_h
        disc_alloy = alloy if blisk else "IN718"
        rho_disc = _rho(disc_alloy)
        if not disc or not disc.get("sigma_d_MPa"):
            pending.append(dict(id=sid, what="диск ступени", why="в записи уровня 1 нет диска (strength.disc)"))
            continue
        Tb = float(st["T_rel_in"]) - 273.15 - (175.0 if cooled else 0.0)
        if i == 0 and t.get("T_blade_C") is not None:
            Tb = float(t["T_blade_C"])
        Tdisc = Tb - T_DISC_DROP
        geom = dict(name=sid, Z=zt, blisk=blisk, blisk_basis=basis, r_hub=r_h, r_tip=r_t, R_live=R, chord=float(st["chord_R"]),
                    bx=float(st["bx_R"]), rpm=rpm, sigma_d_MPa=float(disc["sigma_d_MPa"]), rho_disc=rho_disc, rho_blade=_rho(alloy),
                    t_rim_L1=float(disc.get("t_rim") or 0.0), sections=_free_vortex_sections(st, rpm, r_h, r_t))
        eo = [dict(eo=k, kind="low", source=f"низкий порядок {k}") for k in range(1, LOW_EO + 1)]
        if zn:
            eo += [dict(eo=zn, kind="ngv", source=f"СА {zn} лопаток × 1"), dict(eo=2 * zn, kind="ngv2", source=f"СА {zn} лопаток × 2")]
        if role == "turb" and i == 0 and ninj:
            eo.append(dict(eo=ninj, kind="inj", source=f"горелки камеры: {ninj}"))
        out.append(dict(id=sid, kind="axial", role=role, stage=i + 1, geom=geom, rpm=rpm, alloy=alloy, disc_alloy=disc_alloy,
                        T_spec=dict(type="axial", r_bore=0.20 * R, R=R, r_root=r_h, T_hub=Tdisc - T_HUB_DROP_AX, T_rim=Tdisc, T_blade=Tb),
                        eo_sources=eo, free_stator=None, op=_op(v, role), solid=False))
    return out, pending


def _rho(alloy):
    try:
        if FUNNEL not in sys.path:
            sys.path.insert(0, FUNNEL)
        import make_model as MK           # noqa: E402
        return float(MK.MAT.get(alloy, {}).get("rho") or 8190.0)
    except Exception:                      # noqa: BLE001
        return 8190.0


def variant_parts(v):
    """все роторные детали варианта; pending - что не считается и почему"""
    parts, pending = [], []
    a, b = comp_parts(v)
    parts += a; pending += b
    c = v["concept"]
    if c.get("turb") == "Р":
        a, b = radial_turbine_parts(v)
    else:
        a, b = axial_parts(v, "turb")
    parts += a; pending += b
    if v["meanline"].get("power_turbine"):
        a, b = axial_parts(v, "pt")
        parts += a; pending += b
    if c.get("comb") == "ВР":
        pending.append(dict(id="wave_rotor", what="барабан волнового ротора", why="нет генератора геометрии"))
    return parts, pending


def compact(v):
    """запись уровня 1 -> поля, нужные 3D-расчёту (вход задач AWS - файл случаев)"""
    keep_c = ("idx", "kv", "arch", "rotors", "comp", "turb", "nt", "npt", "cool", "cores", "rotor_join", "m_cw", "m_tw", "mat_cw", "mat_tw",
              "zc", "zt", "zn", "comb", "inlet", "recup", "bear")
    g = v["geometry"]
    return dict(j=v["j"], idx=v["idx"], kv=v.get("kv"), arch=v.get("arch") or v["concept"].get("arch"), niche=v.get("niche"),
                score=v.get("score"), feasible=v.get("feasible"), n_rpm=v["n_rpm"], G=v["G"],
                concept={k: v["concept"].get(k) for k in keep_c},
                geometry=dict(compressor=g.get("compressor"), turbine=g.get("turbine"), power_turbine=g.get("power_turbine"),
                              combustor=g.get("combustor"), inlet=dict(vigv=(g.get("inlet") or {}).get("vigv"))),
                meanline=v["meanline"], strength=v["strength"])
