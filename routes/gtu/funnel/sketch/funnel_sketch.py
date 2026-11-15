"""Эскизы внешнего вида двигателя-заменителя по геометрии уровня 1 воронки (оценка, не конструкторская документация).

Решение владельца 14.11.2026: к каждому из 50 вариантов уровня 1 - эскизы внешнего вида с трёх ракурсов.

Вход - funnel_l1_top50.json: поле geometry каждого варианта (funnel_l1.layout), схема (concept, scheme), модули замены
(modules), масса (mass_kg_estimate). Из funnel_l1_results.json - только подпись расчётной точки (target); полей варианта
там нет, поэтому вариант без geometry не рисуется (причина - в sketch_dims.json). Ничего не досочиняется: узел, размера
которого в данных нет, рисуется условно (штрихом) и попадает в список условных узлов варианта.

Модель - упрощённая трёхмерная модель внешнего вида из тел вращения и параллелепипедов по размерам узлов, м:
  ось x - по валу от входного фланца (x = 0) к выхлопу, y - поперёк, z - вверх от рамы (z = 0);
  ось двигателя - на высоте 0,5·H (правило funnel_l1.layout: высота оси пакета - нет данных, запрос заказчику);
  габарит с обвязкой (штрихпунктир с двумя точками) СОБИРАЕТСЯ из размеров узлов по правилу компоновки уровня 1 (корпус
  компрессора + 0,5 м, 1,3·Dtip силовой турбины + 0,5 м, рекуператор в корпусе по объёму, блок рекуператора сбоку, выносная
  камера, ядра кластера) - это не копия поля envelope; тест сверяет его и габарит всей модели с JSON в пределах 1 %.
  Узел, выходящий за габарит, и пересечение узлов разных ядер не подрезаются, а записываются в issues варианта.
Виды (один масштаб для всего набора, параллельная проекция):
  1 - вид сбоку: ось вала горизонтально, ортогональная проекция (смотрим с -y, вход слева);
  2 - вид с торца со стороны выходного вала, ортогональная проекция;
  3 - аксонометрия 3/4 (параллельная проекция, взгляд спереди-слева-сверху: азимут 55°, возвышение 30°).
Рендер - pyvista (VTK; без экрана - OSMesa), надписи, размеры, выноски, легенда, штамп - matplotlib (Liberation Serif -
метрический аналог Times New Roman). Стиль - ведомость на белой бумаге: монохром, светло-серая заливка, чёрные контуры.

Запуск (интерпретатор с pyvista и matplotlib):
  python funnel_sketch.py --input ../funnel_l1_top50.json --out "$GTU_OUT/sketch" --dims sketch_dims.json
  python funnel_sketch.py --no-render ...            # только модель и sketch_dims.json (numpy, без pyvista)
"""
import argparse, datetime, hashlib, json, math, os, sys, textwrap

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FUNNEL = os.path.dirname(HERE)
DATE = "2026-11-14"
CAPTION = "Эскиз по размерам уровня 1 (оценка), не конструкторская документация"

# =========================================================================== лист (пиксели)
SHEET_W, SHEET_H = 2800, 1750
AREA = (40, 40, 2040, 1470)                 # поле вида x0, y0, x1, y1
PANEL = (2060, 40, 2760, 1470)              # легенда
TITLE = (40, 1490, 2760, 1710)              # штамп
TITLE_SPLIT = 2040                          # граница граф штампа
SS = 2                                      # суперсэмплинг рендера
SEG = 96                                    # сегментов по окружности тела вращения
LW_CONTOUR, LW_EDGE, LW_DASH = 2.4, 1.5, 1.5        # толщины линий, пикс. листа
DASH_PX, GAP_PX = 13.0, 8.0                 # штрих условного узла, пикс.
MOD_SHADES = ["#c4c4c4", "#a0a0a0", "#dcdcdc", "#8e8e8e", "#b4b4b4", "#d0d0d0", "#a8a8a8", "#989898"]
AMBIENT, DIFFUSE = 0.72, 0.34                # тень не темнее #5a: чёрные рёбра на заливке не ниже 3:1 (WCAG, нетекстовые)
NEUTRAL = "#f2f2f2"                         # узлы вне модулей замены: вход, выхлоп, переходники, вал
FONT = "Liberation Serif"
MARGIN_SIDE = dict(x=135, top=125, bottom=200)      # поля под размеры и выноски (подбор масштаба)
MARGIN_AXO = dict(x=150, top=125, bottom=230)
AXO_AZ, AXO_EL = 55.0, 30.0
VIEWS = ("side", "end", "axo")
VIEW_TITLE = {"side": "Вид сбоку: ось вала горизонтально, ортогональная проекция",
              "end": "Вид с торца со стороны выходного вала, ортогональная проекция",
              "axo": "Аксонометрия 3/4 (параллельная проекция)"}
VIEW_FILE = {"side": "1_side", "end": "2_end", "axo": "3_axo"}
TOL_ENV = 0.01                              # допуск сверки габарита, доля


# =========================================================================== числа по-русски
def fnum(x, nd=3):
    s = f"{float(x):.{nd}f}"
    if s.startswith("-") and float(s) == 0.0:
        s = s[1:]
    return s.replace(".", ",")


def fint(x):
    n = int(round(float(x)))
    return f"{n:,}".replace(",", " ") if abs(n) >= 10000 else str(n)


# =========================================================================== примитивы
def _basis(u):
    u = np.asarray(u, float)
    u = u / np.linalg.norm(u)
    a = np.array([0.0, 0.0, 1.0]) if abs(u[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(a, u)
    e1 /= np.linalg.norm(e1)
    return u, e1, np.cross(u, e1)


def rev_mesh(profile, p0, u, seg=SEG):
    """тело вращения: профиль [(s, r)] вдоль оси u от точки p0, торцы закрыты; треугольники наружу"""
    u, e1, e2 = _basis(u)
    p0 = np.asarray(p0, float)
    th = np.linspace(0.0, 2 * math.pi, seg, endpoint=False)
    ring = np.outer(np.cos(th), e1) + np.outer(np.sin(th), e2)
    V = np.vstack([p0 + s * u + r * ring for s, r in profile])
    i = np.arange(seg)
    j = (i + 1) % seg
    F = []
    for k in range(len(profile) - 1):
        if profile[k][1] <= 1e-9 and profile[k + 1][1] <= 1e-9:
            continue
        a, b, c, d = k * seg + i, k * seg + j, (k + 1) * seg + j, (k + 1) * seg + i
        F += [np.stack([a, b, c], 1), np.stack([a, c, d], 1)]
    caps = []
    for k, sign in ((0, -1), (len(profile) - 1, 1)):
        s, r = profile[k]
        if r > 1e-9:
            caps.append((len(V) + len(caps), k, sign, p0 + s * u))
    if caps:
        V = np.vstack([V] + [cp[3][None, :] for cp in caps])
    for ci, k, sign, _ in caps:
        a, b = k * seg + i, k * seg + j
        cc = np.full(seg, ci)
        F.append(np.stack([cc, b, a], 1) if sign < 0 else np.stack([cc, a, b], 1))
    return V, np.vstack(F).astype(np.int64)


def box_mesh(lo, hi):
    (x0, y0, z0), (x1, y1, z1) = lo, hi
    V = np.array([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                  [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]], float)
    F = np.array([[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4],
                  [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]], np.int64)
    return V, F


def box_edges(lo, hi):
    V, _ = box_mesh(lo, hi)
    E = [(0, 1), (3, 2), (0, 3), (1, 2), (4, 5), (7, 6), (4, 7), (5, 6), (0, 4), (1, 5), (2, 6), (3, 7)]   # все рёбра по +оси
    return [V[[a, b]] for a, b in E]


def rev_edges(profile, p0, u, v, seg=SEG):
    """контур условного тела вращения на виде v: окружности в точках профиля и две образующие силуэта"""
    u, e1, e2 = _basis(u)
    v_ = np.asarray(v, float)
    r_ = np.cross(v_, [0.0, 0.0, 1.0])
    if np.linalg.norm(r_) > 1e-6 and np.dot(e1, r_) < -1e-9:   # одна фаза штриха у зеркальных узлов (газоходы ядер)
        e1 = -e1
    if np.dot(np.cross(e1, e2), v_) > 0:
        e2 = -e2
    p0 = np.asarray(p0, float)
    th = np.linspace(0.0, 2 * math.pi, seg + 1)
    ring = np.outer(np.cos(th), e1) + np.outer(np.sin(th), e2)
    out = [p0 + s * u + r * ring for s, r in profile if r > 1e-6]
    d = np.cross(u, np.asarray(v, float))
    if np.linalg.norm(d) > 1e-6:
        d /= np.linalg.norm(d)
        for sg in (1.0, -1.0):
            out.append(np.array([p0 + s * u + sg * r * d for s, r in profile]))
    return out


def dash_segments(polylines, dash, gap):
    """штрихи по длине ломаных (в метрах): массив (N, 2, 3)"""
    segs = []
    for P in polylines:
        P = np.asarray(P, float)
        if len(P) < 2:
            continue
        L = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
        if L[-1] <= 1e-9:
            continue
        t = 0.0
        while t < L[-1]:
            t1 = min(t + dash, L[-1])
            ts = np.r_[t, L[(L > t) & (L < t1)], t1]
            pts = np.array([np.array([np.interp(tt, L, P[:, k]) for k in range(3)]) for tt in ts])
            segs += [pts[m:m + 2] for m in range(len(pts) - 1)]
            t = t1 + gap
    return np.array(segs) if segs else np.zeros((0, 2, 3))


# =========================================================================== модель варианта
class Model:
    """узлы (parts), оси тел (bodies) для поиска радиуса по x, точки подключения, условные узлы, габарит"""

    def __init__(self, entry):
        self.entry = entry
        self.parts, self.bodies, self.conns, self.hidden, self.cond_notes, self.issues = [], [], [], [], [], []
        self.modules = []
        self.enclosed = {}          # номер модуля -> чем закрыт снаружи (на виде сбоку - скрытым контуром)
        self.scale_note = None      # масштаб набора больше наибольшего для варианта (--allow-overflow)

    # ---- узлы
    def _add(self, name, kind, prim, module=None, inst=0, body=None, fit=False, note=None):
        if prim[0] == "rev":
            V, F = rev_mesh(prim[1], prim[2], prim[3])
        else:
            V, F = box_mesh(prim[1], prim[2])
        p = dict(name=name, kind=kind, prim=prim, module=module, inst=inst, body=body, fit=fit, V=V, F=F,
                 lo=V.min(0), hi=V.max(0))
        self.parts.append(p)
        if body is not None and prim[0] == "rev" and abs(prim[3][0]) > 0.99:
            x0 = prim[2][0]
            self.bodies[body]["segs"].append((x0 + prim[1][0][0], x0 + prim[1][-1][0], prim[1]))
        if kind == "cond" and note:
            self.cond_notes.append(note)
        return p

    def rev_x(self, name, kind, profile, body, module=None, inst=0, fit=False, note=None):
        b = self.bodies[body]
        return self._add(name, kind, ("rev", [(s, r) for s, r in profile], (0.0, b["y"], b["z"]), (1.0, 0.0, 0.0)),
                         module, inst, body, fit, note)

    def rev(self, name, kind, profile, p0, u, module=None, inst=0, fit=False, note=None):
        return self._add(name, kind, ("rev", profile, tuple(p0), tuple(u)), module, inst, None, fit, note)

    def box(self, name, kind, lo, hi, module=None, inst=0, fit=False, note=None):
        return self._add(name, kind, ("box", tuple(lo), tuple(hi)), module, inst, None, fit, note)

    def new_body(self, name, y, z, inst=0):
        self.bodies.append(dict(name=name, y=y, z=z, inst=inst, segs=[]))
        return len(self.bodies) - 1

    def r_at(self, body, x):
        """наибольший радиус тел вращения оси body в сечении x (или None)"""
        best = None
        for x0, x1, prof in self.bodies[body]["segs"]:
            if x0 - 1e-9 <= x <= x1 + 1e-9:
                xs = [x0 - prof[0][0] + s for s, _ in prof]
                rs = [r for _, r in prof]
                r = float(np.interp(x, xs, rs)) if xs[-1] > xs[0] else max(rs)
                best = r if best is None else max(best, r)
        return best

    def body_at(self, x, prefer=None):
        """ось, тела которой покрывают x: сначала prefer, потом остальные по порядку"""
        order = ([prefer] if prefer is not None else []) + [i for i in range(len(self.bodies)) if i != prefer]
        for i in order:
            if self.r_at(i, x) is not None:
                return i
        return None


def _shade(no):
    return MOD_SHADES[(no - 1) % len(MOD_SHADES)] if no else NEUTRAL


def build_model(e):
    """трёхмерная модель внешнего вида варианта по его geometry (размеры уровня 1, м)"""
    g = e.get("geometry")
    if not g:
        raise ValueError(f"у варианта idx {e.get('idx')} нет поля geometry - эскиз не строится")
    M = Model(e)
    con = e.get("concept") or {}
    mods = e.get("modules") or []
    names = [m["name"] for m in mods]
    cores = int((g.get("cores") or {}).get("n", 1) or 1)
    cluster = cores > 1
    comb = g["combustor"]
    ctype = comb["type"]
    pt, turb, exh, shaft, rec = g.get("power_turbine"), g["turbine"], g["exhaust"], g["output_shaft"], g.get("recuperator")
    stages, inlet = g["compressor"], g["inlet"]
    hot_out = not str(shaft["output_end"]).startswith("холодный")

    def no(name):
        return names.index(name) + 1 if name in names else None

    rotor_no = no("картридж ротора") or no("ротор")
    hot_no, comb_no, core_no, pt_no, rec_no = no("горячая часть"), no("камера"), no("ядро"), no("силовая турбина"), no("модуль рекуператора")
    comb_in_hot = comb_no is None and not cluster and ctype in ("кольцевая", "RDC")
    if cluster:
        rotor_no = hot_no = comb_no = core_no
    elif comb_in_hot:
        comb_no = hot_no
    for k, m in enumerate(mods):
        M.modules.append(dict(no=k + 1, name=m["name"], mass_kg=m.get("mass_kg"), shade=_shade(k + 1), drawn=False))

    # ---- размеры узлов -> габарит по правилу компоновки уровня 1 (funnel_l1.layout), не из поля envelope
    cc = [s for s in stages if s["type"] != "осевая"]
    k_case = cc[0]["D_case"] / cc[0]["D2"] if cc else 1.0
    Dc = max([s.get("D2", 0.0) for s in cc] + [s.get("D", 0.0) for s in stages if s["type"] == "осевая"])
    D_case_g = k_case * Dc                                   # funnel_space.derive: D_case = kD_case · max D ступеней
    D_pt_case = 1.3 * pt["Dtip"] if pt else None             # корпус СТ - 1,3·Dtip (funnel_space: D_env и стыки модуля)
    L_core = exh["x0"] + exh["L"]
    rec_side = rec is not None and not str(rec.get("shape", "")).startswith("кольцевой")
    if cluster:
        D_cc = D_case_g + 0.5
        kq = math.ceil(math.sqrt(cores))
        rows = math.ceil(cores / kq)
        W, H = D_cc, D_cc
        if ctype == "выносная":
            H += 1.2 * comb["D"]
        W, H = max(W, kq * D_cc), max(H, rows * D_cc)
        D_core = D_cc
    else:
        D_core = D_case_g + 0.5
        if pt:
            D_core = max(D_core, 1.3 * pt["Dtip"] + 0.5)
        if rec is not None and not rec_side:
            D_core = math.sqrt(D_core ** 2 + 4 * rec["V_m3"] / (math.pi * max(D_core, 0.5)))
        W = H = D_core
        if rec_side:
            nrow = max(1, int(L_core // max(rec["module_side"], 1e-6)))
            rows_al = min(nrow, rec["modules"])
            ncol = math.ceil(rec["modules"] / nrow)
            W += ncol * rec["module_side"] + 0.3
        if ctype == "выносная":
            H += 1.2 * comb["D"]
    z_ax = 0.5 * H                                           # правило funnel_l1.layout: ось на 0,5·H
    Lc = shaft["coupling_L"]
    x_lo = -Lc if not hot_out else 0.0
    L_env = L_core + Lc
    y_lo = -0.5 * W if cluster else -0.5 * D_core       # блок рекуператора сбоку - на стороне +y (за двигателем на виде сбоку)
    env = dict(lo=np.array([x_lo, y_lo, 0.0]), hi=np.array([x_lo + L_env, y_lo + W, H]))
    M.env = env
    M.z_ax, M.hot_out, M.cluster, M.cores = z_ax, hot_out, cluster, cores
    M.rule = dict(D_case=D_case_g, D_core=D_core, D_pt_case=D_pt_case, L_core=L_core, L=L_env, W=W, H=H)
    # рекуператор в корпусе: длина и D - из компоновки, положение по длине не задано - по центру (условно)
    drum = None
    if rec is not None and not rec_side:
        xd0 = 0.5 * (L_core - rec["L"])
        drum = dict(x0=xd0, x1=xd0 + rec["L"], r=0.5 * rec["D_outer"])

    def under_drum(x, r):
        """точка оси двигателя (x, радиус r) закрыта барабаном рекуператора"""
        return drum is not None and not cluster and drum["x0"] - 1e-9 <= x <= drum["x1"] + 1e-9 and r < drum["r"] - 1e-6

    # ---- оси тел
    if cluster:
        kq = math.ceil(math.sqrt(cores))
        rows = math.ceil(cores / kq)
        core_bodies = []
        for i in range(cores):
            col, row = i % kq, i // kq
            yc = (col - (kq - 1) / 2) * D_cc
            zc = z_ax + ((rows - 1) / 2 - row) * D_cc
            core_bodies.append(M.new_body(f"ядро {i + 1}", yc, zc, inst=i))
        main = M.new_body("силовая часть", 0.0, z_ax)
    else:
        main = M.new_body("двигатель", 0.0, z_ax)
        core_bodies = [main]
    M.main = main

    # ---- газогенератор: вход, компрессор, камера, турбина - на каждой оси ядра
    for ib, b in enumerate(core_bodies):
        by, bz = M.bodies[b]["y"], M.bodies[b]["z"]
        inst = ib
        if str(inlet["type"]).startswith("двусторонний"):
            Li, Wi, Hi = inlet["L"], inlet["W"], inlet["H"]
            a1, a2 = 0.46 * Li, 0.54 * Li
            for k, (xa, xb) in enumerate(((0.0, a1), (a2, Li))):
                M.box(f"входная камера {k + 1} (двусторонний вход)", "solid", (xa, by - Wi / 2, bz - Hi / 2), (xb, by + Wi / 2, bz + Hi / 2),
                      inst=inst)
                A = inlet["flange_area_m2"] / 2 / cores
                ly = min(Wi, A / max(xb - xa, 1e-6))
                M.box(f"фланец входа B{k + 1}", "cond", (xa, by - ly / 2, bz + Hi / 2), (xb, by + ly / 2, bz + Hi / 2 + 0.02), inst=inst,
                      note="фланцы входа: площадь - из расчёта (25 м/с), форма и место на камере - условно")
                xm, zt = 0.5 * (xa + xb), bz + Hi / 2 + 0.02
                M.conns.append(dict(letter=f"B{k + 1}", kind=f"вход воздуха {k + 1} (двусторонний вход)", x=xm, inst=inst, body=b,
                                    pts=[np.array([xm, by, zt]), np.array([xm, by - ly / 2, zt]), np.array([xm, by - Wi / 2, zt - 0.02]),
                                         np.array([xm, by + Wi / 2, zt - 0.02])]))
            M.cond_notes.append("двусторонний вход: две радиальные камеры, деление длины входа пополам - условно")
            # колесо 1-й ступени с двумя входами - между камерами («радиальная камера с двух сторон колеса»); в компоновке уровня 1
            # ступень 1 стоит за входом (x0 = L входа), поэтому колесо между камерами - условно: Ø по ступени 1, ширина - зазор камер
            s1 = stages[0]
            Dw = s1.get("D2") or s1.get("D")
            M.rev_x("колесо 1-й ступени между камерами (двусторонний вход)", "cond", [(a1, Dw / 2), (a2, Dw / 2)], b, module=rotor_no,
                    inst=inst, note=f"двусторонний вход: колесо 1-й ступени с двумя входами между камерами - условно (Ø {fnum(Dw)} м - "
                                    f"по ступени 1, ширина - зазор между камерами); ступень 1 по данным - за входом, x0 = {fnum(s1['x0'])} м")
            if ib == 0:
                M.issues.append(dict(code="double_inlet_order", node="двусторонний вход",
                                     text=(f"двусторонний вход - «радиальная камера с двух сторон колеса», а компоновка уровня 1 ставит обе "
                                           f"камеры перед ступенью 1 (вход 0…{fnum(Li)} м, ступень 1 с x0 = {fnum(s1['x0'])} м): колесо между "
                                           f"камерами на эскизе - условно штрихом, длина двигателя с таким колесом - нет данных")))
            # тело вращения на месте входа - для радиуса r(x) (не рисуется)
            M.bodies[b]["segs"].append((0.0, Li, [(0.0, Hi / 2), (Li, Hi / 2)]))
        else:
            Din, Li = inlet["D"], inlet["L"]
            M.rev_x("входной коллектор", "solid", [(0.0, Din / 2), (Li, Din / 2)], b, inst=inst)
            A1 = inlet["flange_area_m2"] / cores
            Deq = inlet["flange_D_equiv"] / math.sqrt(cores)
            if hot_out:                              # торец входа свободен - фланец на торце коллектора
                big = (f"; Ø фланца {fnum(Deq)} м больше коллектора (D = 1,6·D1t = {fnum(Din)} м) - переход с воздуховода на "
                       f"коллектор не задан") if Deq > Din + 1e-6 else ""
                M.rev_x("фланец входа B", "cond", [(0.0, Deq / 2)], b, inst=inst,
                        note="фланец входа: площадь - из расчёта (25 м/с), круглая форма на торце - условно" + big)
                M.conns.append(dict(letter="B", kind="вход воздуха (осевой, с коллектором)", x=0.0, inst=inst, body=b,
                                    pts=[np.array([0.0, by, bz + 0.25 * Din])] + [np.array([0.0, by - 0.5 * Din * math.cos(t), bz + 0.5 * Din * math.sin(t)])
                                                                                 for t in np.linspace(0.2, 2.9, 6)]))
            else:                                    # привод с холодного конца: на торце входа - вал и муфта, воздух - к коллектору сверху
                ly = min(Din, A1 / max(Li, 1e-6))
                zt = bz + Din / 2
                fits = Li * ly >= 0.99 * A1
                M.box("фланец входа B (сверху коллектора)", "cond", (0.0, by - ly / 2, zt), (Li, by + ly / 2, zt + 0.02), inst=inst, fit=True,
                      note="привод с холодного конца: торец входа занят выходным валом и муфтой - подвод воздуха к коллектору показан "
                           "сверху (радиально), место и форма - условно; площадь фланца - из расчёта (25 м/с)"
                           + ("" if fits else f", на коллекторе L × D = {fnum(Li)} × {fnum(Din)} м помещается {fnum(Li * ly, 2)} м² из "
                                              f"{fnum(A1, 2)} м² - переход с воздуховода не задан"))
                xm = 0.5 * Li
                M.conns.append(dict(letter="B", kind="вход воздуха (к коллектору сверху: привод с холодного конца)", x=xm, inst=inst, body=b,
                                    pts=[np.array([xm, by, zt + 0.02]), np.array([xm, by - ly / 2, zt + 0.02]), np.array([xm, by - Din / 2, bz])]))
        for s in stages:
            x0, L = s["x0"], s["L"]
            if s["type"] == "осевая":
                prof = [(x0, s["D"] / 2), (x0 + L, s["D"] / 2)]
            else:
                Lax = max(L - 0.3 * s["D2"], 0.3 * L)
                prof = [(x0, s["D1t"] / 2), (x0 + 0.55 * Lax, s["D2"] / 2), (x0 + Lax, s["D_case"] / 2), (x0 + L, s["D_case"] / 2)]
            M.rev_x(f"ступень компрессора {s['stage']} ({s['type']})", "solid", prof, b, module=rotor_no, inst=inst)
        x0c, Lcb, Dcb = comb["x0"], comb["L"], comb["D"]
        if ctype in ("кольцевая", "RDC"):
            M.rev_x("камера сгорания" + (" RDC" if ctype == "RDC" else ""), "solid", [(x0c, Dcb / 2), (x0c + Lcb, Dcb / 2)], b,
                    module=comb_no, inst=inst)
            fuel_pts = [np.array([x0c + 0.5 * Lcb, by, bz + Dcb / 2])]
        elif ctype == "ВР":
            M.rev_x("камера сгорания (волновой ротор)", "solid", [(x0c, Dcb / 2), (x0c + Lcb, Dcb / 2)], b, module=comb_no, inst=inst)
            x_d1 = turb["x0"]
            M.rev_x("барабан волнового ротора", "cond", [(x0c + Lcb, Dcb / 2), (x_d1, Dcb / 2)], b, module=comb_no, inst=inst, fit=True,
                    note="барабан волнового ротора: длина - из компоновки (0,6·Dc), диаметр не задан - условно D камеры")
            fuel_pts = [np.array([x0c + 0.5 * Lcb, by, bz + Dcb / 2])]
        else:                                            # выносная: ось вертикальная над корпусом (funnel_l1.layout)
            x_t = turb["x0"]
            Dh = turb["D_hot_case"]
            M.rev_x("переходный корпус под выносную камеру", "cond", [(x0c, Dh / 2), (x_t, Dh / 2)], b, inst=inst, fit=True,
                    note="переходный корпус под выносной камерой: длина - из компоновки (0,25·Dc), диаметр не задан - условно D горячей части")
            xc = 0.5 * (x0c + x_t)
            zc = bz + float(comb.get("z_center", 0.5 * D_case_g + 0.6 * Dcb))
            M.rev("выносная камера сгорания", "solid", [(0.0, Dcb / 2), (Lcb, Dcb / 2)], (xc, by, zc - Lcb / 2), (0.0, 0.0, 1.0),
                  module=comb_no, inst=inst)
            r_loc = Dh / 2
            if zc - Lcb / 2 > bz + r_loc:
                M.rev("патрубки выносной камеры", "cond", [(0.0, 0.25 * Dcb), (zc - Lcb / 2 - bz - r_loc, 0.25 * Dcb)],
                      (xc, by, bz + r_loc), (0.0, 0.0, 1.0), inst=inst, fit=True,
                      note="патрубки выносной камеры: размеров нет - условно")
            fuel_pts = [np.array([xc, by, zc + Lcb / 2])] + [np.array([xc + 0.49 * Dcb * math.sin(t), by - 0.49 * Dcb * math.cos(t), zc + Lcb / 2 - 0.01])
                                                             for t in np.linspace(0, 2 * math.pi, 12, endpoint=False)]
            M.ext_comb = dict(xc=xc, zc=zc, L=Lcb, D=Dcb)
        nf = sum(cn["n"] for cn in g["connections"] if cn["kind"] == "топливо")
        fuel_kind = f"топливо ({nf} форс." + (", " + str(comb.get("placement")) if ctype == "выносная" else "") + ")"
        if ctype != "выносная" and under_drum(fuel_pts[0][0], Dcb / 2):
            fuel_pts = [np.array([fuel_pts[0][0], by, bz + drum["r"]])]
            fuel_kind += "; камера под барабаном рекуператора, буква на поверхности барабана"
            M.cond_notes.append("топливо и выхлоп под барабаном рекуператора: буквы D и C - на наружной поверхности над точкой "
                                "подключения, вывод через корпус рекуператора - условно")
        M.conns.append(dict(letter="D", kind=fuel_kind, x=x0c, inst=inst, body=b, pts=fuel_pts))
        tx0, tL, Dh = turb["x0"], turb["L"], turb["D_hot_case"]
        M.rev_x(f"турбина газогенератора ({turb['type']}, {turb['stages']} ст.)", "solid", [(tx0, Dh / 2), (tx0 + tL, Dh / 2)], b,
                module=hot_no, inst=inst)

    # ---- силовая часть, выхлоп, вал
    t_end = turb["x0"] + turb["L"]
    if pt:
        Dp = D_pt_case
        if cluster:
            for ib, b in enumerate(core_bodies):
                p0 = np.array([t_end, M.bodies[b]["y"], M.bodies[b]["z"]])
                p1 = np.array([pt["x0"], 0.0, z_ax])
                d = p1 - p0
                M.rev(f"газоход ядра {ib + 1} к силовой турбине", "cond",
                      [(0.0, turb["D_hot_case"] / 2), (float(np.linalg.norm(d)), turb["D_hot_case"] / 2)],
                      p0, d / np.linalg.norm(d), inst=ib, fit=True, note="газоходы от ядер к силовой турбине: размеров нет - условно")
        else:
            M.rev_x("переходный канал к силовой турбине", "solid", [(t_end, turb["D_hot_case"] / 2), (pt["x0"], Dp / 2)], main)
        M.rev_x(f"силовая турбина ({pt['stages']} ст.)", "solid", [(pt["x0"], Dp / 2), (pt["x0"] + pt["L"], Dp / 2)], main, module=pt_no)
        M.cond_notes.append("корпус силовой турбины - 1,3·Dtip (правило funnel_space: габарит и стыки модуля), не расчёт корпуса")
    De = exh["D"]
    M.rev_x("выхлопной коллектор (радиальный)", "solid", [(exh["x0"], De / 2), (exh["x0"] + exh["L"], De / 2)], main)
    A = exh["flange_area_m2"]
    xe0, xe1 = exh["x0"], exh["x0"] + exh["L"]
    top = z_ax + De / 2
    c_kind, c_note = "выхлоп (радиальный коллектор, вверх)", ""
    if drum is not None and drum["x1"] > xe0 + 1e-9 and drum["r"] > De / 2:
        ov = min(drum["x1"], xe1) - xe0                  # барабан по центру закрывает начало выхлопного коллектора
        if xe1 - drum["x1"] >= max(0.2 * exh["L"], 0.1):  # фланец - на открытой части коллектора
            xe0 = drum["x1"]
            c_note = (f"; барабан рекуператора закрывает коллектор на {fnum(ov)} м из {fnum(exh['L'])} м - фланец показан на открытой "
                      f"части коллектора")
            c_kind = "выхлоп (радиальный коллектор, вверх; фланец на открытой части коллектора за барабаном рекуператора)"
        else:                                            # коллектор закрыт почти целиком - фланец на поверхности барабана
            top = z_ax + drum["r"]
            c_note = "; коллектор закрыт барабаном рекуператора - фланец показан на поверхности барабана"
            c_kind = "выхлоп (радиальный коллектор под барабаном рекуператора; фланец вверх на поверхности барабана)"
    xm_e = 0.5 * (xe0 + xe1)
    ly = min(De, A / max(xe1 - xe0, 1e-6))
    stub = max(0.0, min(0.1 * De, H - top))
    M.box("фланец выхлопа C", "cond", (xe0, -ly / 2, top - 0.05 * De), (xe1, ly / 2, top + stub), fit=True,
          note="фланец выхлопа: площадь - из расчёта (40 м/с), направление вверх и форма - условно (данных пакета нет)" + c_note)
    M.conns.append(dict(letter="C", kind=c_kind, x=xm_e, inst=0, body=main,
                        pts=[np.array([xm_e, 0.0, top + stub]), np.array([xm_e, -ly / 2, top + stub]),
                             np.array([xm_e, 0.0, top]), np.array([xm_e, -0.5 * De * math.sin(1.2), z_ax + 0.5 * De * math.cos(1.2)])]))
    d_sh, D_fl = shaft["d_shaft"], shaft["coupling_flange_D"]
    t_fl = min(0.2 * D_fl, 0.3 * Lc)
    if hot_out:
        xs = L_core
        prof = [(xs, d_sh / 2), (xs + Lc - t_fl, d_sh / 2), (xs + Lc - t_fl, D_fl / 2), (xs + Lc, D_fl / 2)]
        x_end = xs + Lc
    else:
        xs = 0.0
        prof = [(xs - Lc, D_fl / 2), (xs - Lc + t_fl, D_fl / 2), (xs - Lc + t_fl, d_sh / 2), (xs, d_sh / 2)]
        x_end = xs - Lc
    M.rev_x("выходной вал и фланец муфты", "solid", prof, main)
    M.cond_notes.append("толщина фланца муфты - условно (задан только диаметр фланца и длина муфты)")
    M.shaft_end = np.array([x_end, 0.0, z_ax])
    M.conns.append(dict(letter="A", kind=f"выходной вал: {fint(shaft['n_rpm'])} об/мин, {fnum(shaft['torque_kNm'], 2)} кН·м, вал Ø{fnum(d_sh)} м, "
                                         f"фланец муфты Ø{fnum(D_fl)} м", x=x_end, inst=0, body=main,
                        pts=[M.shaft_end + np.array([0.0, -0.3 * D_fl * math.cos(t), 0.3 * D_fl * math.sin(t)]) for t in (1.2, 0.6, 2.0)] + [M.shaft_end]))

    # ---- рекуператор
    if rec is not None:
        if rec_side:
            side = rec["module_side"]
            nrow = max(1, int(L_core // max(side, 1e-6)))
            rows_al = min(nrow, rec["modules"])
            xs0 = 0.5 * (L_core - rows_al * side)
            y0 = 0.5 * D_core + 0.3                    # сторона +y: на виде сбоку и в аксонометрии блок за двигателем, а не перед ним
            for m_ in range(rec["modules"]):
                r_, q_ = m_ % rows_al, m_ // rows_al
                M.box(f"модуль рекуператора {m_ + 1}", "solid", (xs0 + r_ * side, y0 + q_ * side, 0.0),
                      (xs0 + (r_ + 1) * side, y0 + (q_ + 1) * side, side), module=rec_no, inst=0)
            M.box("воздуховоды и газоходы рекуператора", "cond", (xs0, 0.45 * D_case_g, 0.35 * side), (xs0 + rows_al * side, y0, 0.65 * side),
                  fit=True, note="рекуператор сбоку: сторона (+y, за двигателем на виде сбоку), положение блока по длине и высоте и "
                                 "воздуховоды - условно (в данных - ряды, сторона модуля, зазор 0,3 м)")
        else:
            Lr, Dr = rec["L"], rec["D_outer"]
            xs0 = drum["x0"]
            M.rev_x(f"рекуператор в корпусе ({rec['modules']} мод.)", "solid", [(xs0, Dr / 2), (xs0 + Lr, Dr / 2)], main, module=rec_no)
            M.cond_notes.append("рекуператор в корпусе: D и длина - из компоновки, положение по длине (по центру) и разбивка на модули - условно")
            # что закрыто барабаном: модули целиком внутри (на виде сбоку - скрытым контуром), вход и выхлоп - на сколько
            if not cluster:
                for mm in M.modules:
                    ps = [p for p in M.parts if p["module"] == mm["no"] and p["kind"] == "solid"]
                    if mm["no"] == rec_no or not ps:
                        continue
                    if all(p["body"] == main and p["lo"][0] >= drum["x0"] - 1e-6 and p["hi"][0] <= drum["x1"] + 1e-6
                           and max(r for _, r in p["prim"][1]) < drum["r"] - 1e-6 for p in ps):
                        M.enclosed[mm["no"]] = "барабан рекуператора"
                free = exh["x0"] - inlet["L"]
                cov_in = max(0.0, inlet["L"] - drum["x0"])
                cov_ex = max(0.0, min(drum["x1"], exh["x0"] + exh["L"]) - exh["x0"])
                if cov_in > 1e-6 or cov_ex > 1e-6:
                    enc = ", ".join(str(n) for n in sorted(M.enclosed))
                    M.issues.append(dict(code="recup_overlap", node=f"рекуператор в корпусе ({rec['modules']} мод.)",
                                         by_m=round(float(Lr - free), 3),
                                         text=(f"рекуператор в корпусе: длина 0,8·L_core = {fnum(Lr)} м больше участка между входом и выхлопом "
                                               f"({fnum(free)} м), положение по длине в компоновке не задано; по центру барабан закрывает "
                                               + ", ".join(t for t in (f"конец входа на {fnum(cov_in)} м" if cov_in > 1e-6 else "",
                                                                       f"выхлопной коллектор на {fnum(cov_ex)} м из {fnum(exh['L'])} м"
                                                                       if cov_ex > 1e-6 else "") if t)
                                               + ((f"; модуль {enc} - внутри барабана" if len(M.enclosed) == 1 else f"; модули {enc} - внутри барабана")
                                                  + " (на виде сбоку - скрытым контуром)" if enc else ""))))

    # ---- опоры (положение - из компоновки, форма - условно)
    for s in g.get("supports") or []:
        xs_ = s["x"]
        bodies = [b for b in range(len(M.bodies)) if M.r_at(b, xs_) is not None]
        zmin = min((M.bodies[b]["z"] for b in bodies), default=None)
        for b in bodies:
            if abs(M.bodies[b]["z"] - zmin) > 1e-6:
                continue
            r = M.r_at(b, xs_)
            by, bz = M.bodies[b]["y"], M.bodies[b]["z"]
            if bz - 0.92 * r > 0.02:
                M.box(s["name"], "cond", (xs_ - 0.04, by - 0.55 * r, 0.0), (xs_ + 0.04, by + 0.55 * r, bz - 0.92 * r), fit=True,
                      note="опоры: положение по оси - из компоновки уровня 1 (доли длины), форма и размер - условно")

    # ---- точки подключения из geometry.connections (положение по оси - из компоновки, по окружности - условно)
    # у каждой точки своя буква: масло у передней и задней опоры - E1 и E2, два воздушных отбора - G1 и G2, электрика - всегда F1…
    cnt = {L_: sum(1 for cn in g["connections"] if cn["kind"].startswith(p_)) for L_, p_ in (("E", "масло"), ("G", "воздух"))}
    lettered = dict(E=0, F=0, G=0)
    for cn in g["connections"]:
        kind = cn["kind"]
        if kind == "топливо":
            continue
        n_txt = f", {cn['n']} шт." if cn.get("n", 1) > 1 else ""
        if kind.startswith("масло"):
            lettered["E"] += 1
            letter, pref = "E" + (str(lettered["E"]) if cnt["E"] > 1 else ""), "side"
            n_txt = f", {cn['n']} линии" if cn.get("n", 1) > 1 else ""       # n - линии у одной опоры (подвод и слив), не число опор
        elif kind.startswith("электрика"):
            lettered["F"] += 1
            letter, pref = f"F{lettered['F']}", "top"
        elif kind.startswith("воздух"):
            lettered["G"] += 1
            letter, pref = "G" + (str(lettered["G"]) if cnt["G"] > 1 else ""), "top"
        else:
            letter, pref = "?", "top"
        b = M.body_at(cn["x"], prefer=core_bodies[0])
        if b is None:
            continue
        M.conns.append(dict(letter=letter, kind=f"{kind}{n_txt} ({cn['where']})", x=cn["x"], inst=0, body=b, ring=(cn["x"], b, pref)))
    M.cond_notes.append("точки подключения: положение по оси - из компоновки, по окружности - условно")

    # ---- модули без размеров: блок датчиков, привод ВНА - условная коробка у точки подключения
    for mm in M.modules:
        nm = mm["name"]
        x = None
        if nm == "блок датчиков и жгут":
            x = next((cn["x"] for cn in g["connections"] if cn["kind"].startswith("электрика: жгут")), None)
        elif nm == "привод ВНА":
            x = next((cn["x"] for cn in g["connections"] if "ВНА" in cn["kind"]), None)
            if x is None:
                x = 0.5 * inlet["L"]
        if x is None:
            continue
        b = M.body_at(x, prefer=core_bodies[0])
        if b is None:
            continue
        r = M.r_at(b, x)
        by, bz = M.bodies[b]["y"], M.bodies[b]["z"]
        a = 0.11
        yc_ = max(by - r - a - 0.01, M.env["lo"][1] + a)     # снаружи корпуса со стороны вида сбоку, но в габарите
        M.box(nm, "cond", (x - a, yc_ - a, bz - a), (x + a, yc_ + a, bz + a), module=mm["no"], fit=True,
              note="блок датчиков, привод ВНА: размер и место не заданы - условная коробка у точки подключения")

    # ---- скрытые узлы: корпуса магнитных подшипников (на виде сбоку - пунктиром с точками)
    for mb in g.get("magnetic_bearings") or []:
        for b in core_bodies:
            M.hidden.append(dict(name=mb["name"], x=mb["x"], L=mb["L"], D=mb["D"], y=M.bodies[b]["y"], z=M.bodies[b]["z"]))
    for mm in M.modules:
        mm["drawn"] = any(p["module"] == mm["no"] for p in M.parts)
    M.cond_notes = list(dict.fromkeys(M.cond_notes))
    check_model(M)
    return M


# =========================================================================== габарит и проверки
def nodes_bbox(M, include_fit=True):
    ps = [p for p in M.parts if include_fit or not p["fit"]]
    return np.min([p["lo"] for p in ps], 0), np.max([p["hi"] for p in ps], 0)


def model_bbox(M):
    lo, hi = nodes_bbox(M)
    return np.minimum(lo, M.env["lo"]), np.maximum(hi, M.env["hi"])


def check_model(M):
    """узлы за габаритом (кроме условных, вписанных в габарит) и пересечение узлов разных ядер"""
    lo_e, hi_e = M.env["lo"], M.env["hi"]
    ext = hi_e - lo_e
    for p in M.parts:
        if p["fit"]:
            if np.any(p["lo"] < lo_e - 1e-6) or np.any(p["hi"] > hi_e + 1e-6):
                M.issues.append(dict(code="fit_outside", node=p["name"], text=f"условный узел {p['name']} вне габарита (ошибка эскиза)"))
            continue
        over_hi = p["hi"] - hi_e
        over_lo = lo_e - p["lo"]
        for k, ax in enumerate("xyz"):
            for over, side in ((over_hi[k], "+"), (over_lo[k], "-")):
                if over > max(0.005 * ext[k], 0.003):
                    M.issues.append(dict(code="protrudes", node=p["name"], axis=ax + side, by_m=round(float(over), 3),
                                         text=f"{p['name']} выходит за габарит уровня 1 по {ax}{side} на {fnum(over)} м"))
    if hasattr(M, "ext_comb"):
        ec = M.ext_comb
        H_add = 1.2 * ec["D"]
        top = ec["zc"] + ec["L"] / 2
        M.issues.append(dict(code="ext_comb_rule", node="выносная камера сгорания", by_m=round(float(top - hi_e[2]), 3),
                             text=(f"выносная камера: ось вертикальная, L = {fnum(ec['L'])} м = {fnum(ec['L'] / ec['D'], 1)}·D, а компоновка уровня 1 "
                                   f"добавила к H только 1,2·D = {fnum(H_add)} м и центр камеры поставила на 0,5·D корпуса + 0,6·D над осью: "
                                   f"верх камеры {fnum(top)} м при H = {fnum(hi_e[2])} м, низ камеры заходит в корпус - несогласованность "
                                   f"funnel_l1.layout (эскиз рисует по данным, не подрезая)")))
    if M.cluster:
        solids = [p for p in M.parts if p["kind"] == "solid" and p["body"] is None and p["name"].startswith("выносная")]
        for p in solids:
            for q in M.parts:
                if q["kind"] != "solid" or q["inst"] == p["inst"] or q["body"] is None:
                    continue
                if np.all(p["lo"] < q["hi"] - 1e-6) and np.all(q["lo"] < p["hi"] - 1e-6):
                    M.issues.append(dict(code="intersect", node=p["name"], other=q["name"],
                                         text=f"{p['name']} ядра {p['inst'] + 1} пересекает узел «{q['name']}» ядра {q['inst'] + 1}"))
                    break


def dims_of(M):
    """габарит JSON, габарит модели и узлов, отклонения - для sketch_dims.json и теста"""
    e = M.entry
    envj = e["geometry"]["envelope"]
    lo, hi = model_bbox(M)
    nlo, nhi = nodes_bbox(M, include_fit=False)
    model = dict(L=float(hi[0] - lo[0]), W=float(hi[1] - lo[1]), H=float(hi[2] - lo[2]))
    rule = dict(L=M.rule["L"], W=M.rule["W"], H=M.rule["H"])
    nodes = dict(L=float(nhi[0] - nlo[0]), W=float(nhi[1] - nlo[1]), H=float(nhi[2] - nlo[2]))
    dev = {k: (model[k] / envj[k] - 1.0) for k in "LWH"}
    dev_rule = {k: (rule[k] / envj[k] - 1.0) for k in "LWH"}
    return dict(json={k: round(float(envj[k]), 4) for k in "LWH"}, model={k: round(v, 4) for k, v in model.items()},
                rule={k: round(v, 4) for k, v in rule.items()}, nodes_solid={k: round(v, 4) for k, v in nodes.items()},
                dev_pct={k: round(100 * v, 3) for k, v in dev.items()}, dev_rule_pct={k: round(100 * v, 3) for k, v in dev_rule.items()},
                ok=all(abs(v) <= TOL_ENV for v in dev.values()))


# =========================================================================== виды и масштаб
def view_basis(M, view):
    if view == "side":
        v = np.array([0.0, 1.0, 0.0])
    elif view == "end":
        v = np.array([-1.0, 0.0, 0.0]) if M.hot_out else np.array([1.0, 0.0, 0.0])
    else:
        az, el = math.radians(AXO_AZ), math.radians(AXO_EL)
        v = np.array([math.cos(el) * math.cos(az), math.cos(el) * math.sin(az), -math.sin(el)])
    up = np.array([0.0, 0.0, 1.0])
    right = np.cross(v, up)
    right /= np.linalg.norm(right)
    up2 = np.cross(right, v)
    return v, right, up2


def env_corners(M):
    lo, hi = M.env["lo"], M.env["hi"]
    return np.array([[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1]) for z in (lo[2], hi[2])])


def bbox_corners(M):
    """углы габарита всей модели (габарит с обвязкой и выступающие узлы)"""
    lo, hi = model_bbox(M)
    return np.array([[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1]) for z in (lo[2], hi[2])])


def view_extent(M, view):
    """размах модели на виде, м (ширина, высота)"""
    v, right, up2 = view_basis(M, view)
    C = bbox_corners(M)
    a, b = C @ right, C @ up2
    return float(a.max() - a.min()), float(b.max() - b.min())


def max_scale(M):
    """наибольший масштаб, пикс./м, при котором все три вида варианта входят в поле (с полями под размеры и выноски)"""
    Wd, Hd = AREA[2] - AREA[0], AREA[3] - AREA[1]
    s = 1e9
    for view in VIEWS:
        w, h = view_extent(M, view)
        m = MARGIN_AXO if view == "axo" else MARGIN_SIDE
        s = min(s, (Wd - 2 * m["x"]) / max(w, 1e-6), (Hd - m["top"] - m["bottom"]) / max(h, 1e-6))
    return s


def common_scale(models):
    """один масштаб для всего набора: наибольший целый десяток пикс./м, при котором любой вид любого варианта входит в поле"""
    return max(10, int(min(max_scale(M) for M in models) // 10) * 10)


class Proj:
    """параллельная проекция вида на лист: мир -> пиксели листа"""

    def __init__(self, M, view, s):
        self.v, self.right, self.up = view_basis(M, view)
        self.s = s
        Wd, Hd = AREA[2] - AREA[0], AREA[3] - AREA[1]
        self.cx, self.cy = AREA[0] + Wd / 2, AREA[1] + Hd / 2
        C = bbox_corners(M)
        c3 = 0.5 * (C.min(0) + C.max(0))
        a, b = (C - c3) @ self.right * s, -(C - c3) @ self.up * s
        mx, my = 0.5 * (a.max() + a.min()), 0.5 * (b.max() + b.min())
        m = MARGIN_AXO if view == "axo" else MARGIN_SIDE
        dy = 0.5 * (m["top"] - m["bottom"])
        self.f = c3 - self.right * ((0.0 - mx) / s) + self.up * ((dy - my) / s)

    def __call__(self, p):
        p = np.atleast_2d(np.asarray(p, float)) - self.f
        return np.stack([self.cx + p @ self.right * self.s, self.cy - p @ self.up * self.s], 1)

    def dist(self, p, cam):
        return (np.atleast_2d(p) - cam) @ self.v


# =========================================================================== рендер вида
def _pv_mesh(pv, V, F):
    return pv.PolyData(V, np.hstack([np.full((len(F), 1), 3, np.int64), F]).ravel())


def _pv_lines(pv, segs):
    if len(segs) == 0:
        return None
    P = segs.reshape(-1, 3)
    n = len(segs)
    L = np.hstack([np.full((n, 1), 2, np.int64), np.arange(2 * n).reshape(n, 2)]).ravel()
    return pv.PolyData(P, lines=L)


def render_image(M, view, s):
    """pyvista: заливка, контуры, рёбра, штрихи условных узлов; возврат: картинка поля вида, глубина, проекция, камера"""
    import pyvista as pv
    P = Proj(M, view, s)
    Wd, Hd = AREA[2] - AREA[0], AREA[3] - AREA[1]
    pl = pv.Plotter(off_screen=True, window_size=(Wd * SS, Hd * SS), lighting="none")
    pl.set_background("white")
    D = 60.0
    cam = P.f - P.v * D
    pl.add_light(pv.Light(position=tuple(P.f - P.v * 10 + P.up * 9 - P.right * 6), focal_point=tuple(P.f), color="white",
                          intensity=0.55, light_type="scene light"))
    pl.add_light(pv.Light(position=tuple(cam), focal_point=tuple(P.f), color="white", intensity=0.45, light_type="scene light"))
    lift = -P.v * (1.5 / s)
    for p in M.parts:
        if p["kind"] != "solid":
            continue
        mesh = _pv_mesh(pv, p["V"], p["F"])
        smooth = p["prim"][0] == "rev"
        pl.add_mesh(mesh, color=_shade(p["module"]), ambient=AMBIENT, diffuse=DIFFUSE, specular=0.0, smooth_shading=smooth,
                    split_sharp_edges=smooth, feature_angle=35.0, silhouette=dict(color="black", line_width=LW_CONTOUR * SS))
        fe = mesh.extract_feature_edges(feature_angle=35.0, boundary_edges=False, non_manifold_edges=False, feature_edges=True,
                                        manifold_edges=False)
        if fe.n_points:
            fe.points = fe.points + lift
            pl.add_mesh(fe, color="black", line_width=LW_EDGE * SS)
    segs = []
    for p in M.parts:
        if p["kind"] != "cond":
            continue
        pr = p["prim"]
        pls = box_edges(pr[1], pr[2]) if pr[0] == "box" else rev_edges(pr[1], pr[2], pr[3], P.v)
        sg = dash_segments(pls, DASH_PX / s, GAP_PX / s)
        if len(sg):
            segs.append(sg + lift * 1.5)
    dashed = pl.add_mesh(_pv_lines(pv, np.vstack(segs)), color="black", line_width=LW_DASH * SS) if segs else None
    pl.enable_parallel_projection()
    pl.camera.position = tuple(cam)
    pl.camera.focal_point = tuple(P.f)
    pl.camera.up = tuple(P.up)
    pl.camera.parallel_scale = (Hd / 2) / s
    pl.camera.clipping_range = (1.0, 2 * D)
    if dashed is not None:                      # глубина - по телам, без штрихов условных узлов
        dashed.SetVisibility(False)
    pl.screenshot(return_img=True)
    depth = pl.get_image_depth(fill_value=np.nan, reset_camera_clipping_range=False)
    if dashed is not None:
        dashed.SetVisibility(True)
        pl.render()
    img = pl.screenshot(return_img=True)
    pl.close()
    from PIL import Image
    im = Image.fromarray(img).resize((Wd, Hd), Image.LANCZOS)
    return np.asarray(im), -np.asarray(depth), P, cam


# =========================================================================== надписи, размеры, выноски
def _enclosed_parts(M, no):
    return [p for p in M.parts if p["module"] == no and p["kind"] == "solid" and p["prim"][0] == "rev" and p["body"] is not None]


def hidden_outlines(M):
    """скрытые контуры модулей внутри барабана рекуператора на виде сбоку: ломаные (мир, м) в плоскости оси"""
    out = []
    for no in sorted(M.enclosed):
        for p in _enclosed_parts(M, no):
            prof, (_, by, bz) = p["prim"][1], p["prim"][2]
            top = [np.array([s, by, bz + r]) for s, r in prof]
            bot = [np.array([s, by, bz - r]) for s, r in prof]
            out += [np.array(top), np.array(bot), np.array([bot[0], top[0]]), np.array([bot[-1], top[-1]])]
    return out


def hidden_anchor(M, no):
    """точка выноски модуля внутри барабана: середина верхней образующей его наибольшего узла"""
    p = max(_enclosed_parts(M, no), key=lambda q: max(r for _, r in q["prim"][1]))
    prof, (_, by, bz) = p["prim"][1], p["prim"][2]
    xm = 0.5 * (prof[0][0] + prof[-1][0])
    r = float(np.interp(xm, [s for s, _ in prof], [r for _, r in prof]))
    return np.array([xm, by, bz + r])


def _visible(P, cam, depth, p, tol=0.04):
    X, Y = P(p)[0]
    i, j = int((Y - AREA[1]) * SS), int((X - AREA[0]) * SS)
    h, w = depth.shape
    if not (0 <= i < h and 0 <= j < w):
        return False
    win = depth[max(0, i - 3):i + 4, max(0, j - 3):j + 4]
    win = win[np.isfinite(win)]
    if win.size == 0:
        return True
    return float(P.dist(p, cam)[0]) <= float(win.max()) + tol


def anchors(M, P, cam, depth, view):
    """точки выносок: модули (по экземпляру) и подключения - видимые на этом виде"""
    labels = []
    groups = {}
    for p in M.parts:
        if p["module"]:
            groups.setdefault((p["module"], p["inst"] if p["name"].find("рекуператора") < 0 else 0), []).append(p)
    for (no, inst), ps in sorted(groups.items()):
        if no in M.enclosed:                     # модуль целиком внутри барабана: выноска - к скрытому контуру на виде сбоку
            if view == "side":
                w = hidden_anchor(M, no)
                labels.append(dict(text=str(no), kind="mod", anchor=P(w)[0], world=w, hidden=True))
            continue
        cand = []
        for p in ps:
            V, F = p["V"], p["F"]
            C = V[F].mean(1)
            if p["kind"] == "cond":                  # центр, центры граней, вершины
                lo_, hi_ = p["lo"], p["hi"]
                c_ = 0.5 * (lo_ + hi_)
                fc = [np.where(np.arange(3) == k, side, c_) for k in range(3) for side in (lo_, hi_)]
                C = np.vstack([c_] + fc + list(V))
            step = max(1, len(C) // 300)
            cand.append(C[::step])
        C = np.vstack(cand)
        vis = np.array([_visible(P, cam, depth, c) for c in C])
        if not vis.any():
            continue
        Cv = C[vis]
        S = P(Cv)
        med = np.median(S, 0)
        k = int(np.argmin(np.linalg.norm(S - med, axis=1)))
        labels.append(dict(text=str(no), kind="mod", anchor=S[k], world=Cv[k]))
    shown = set()
    for cn in M.conns:
        if cn["letter"] in shown:          # одна буква на вид (у ядер кластера - у первого видимого)
            continue
        if "ring" in cn:
            x, b, pref = cn["ring"]
            r = M.r_at(b, x) or 0.1
            by, bz = M.bodies[b]["y"], M.bodies[b]["z"]
            ang0 = math.pi / 2 if pref == "top" else math.pi
            angs = sorted(np.linspace(0, 2 * math.pi, 24, endpoint=False), key=lambda a: abs((a - ang0 + math.pi) % (2 * math.pi) - math.pi))
            pts = [np.array([x, by + r * math.cos(a), bz + r * math.sin(a)]) for a in angs]
        else:
            pts = cn["pts"]
        for p in pts:
            if _visible(P, cam, depth, p):
                a = P(p)[0]
                near = next((lb for lb in labels if lb["kind"] == "conn" and np.linalg.norm(lb["anchor"] - a) < 8.0), None)
                if near is not None:               # точки совпали на виде (камеры двустороннего входа с торца) - одна выноска
                    near["text"] += ", " + cn["letter"]
                else:
                    labels.append(dict(text=cn["letter"], kind="conn", anchor=a, world=p))
                shown.add(cn["letter"])
                break
    return labels


def _label_w(lb):
    return 36 if lb["kind"] == "mod" else 16 + 13 * len(lb["text"])


def place_rails(labels, ex0, ey0, ex1, ey1, bottom_y=None):
    """выноски на две «полки» - над и под видом; без наложений, порядок по x сохраняется"""
    out = []
    xmin, xmax = AREA[0] + 30, AREA[2] - 30
    mid = 0.5 * (ey0 + ey1)
    for side in ("top", "bottom"):
        L = sorted([lb for lb in labels if (lb["anchor"][1] <= mid) == (side == "top")], key=lambda lb: lb["anchor"][0])
        if not L:
            continue
        need = sum(_label_w(lb) + 12 for lb in L)
        levels = 2 if need > (xmax - xmin) * 0.55 else 1
        for lev in range(levels):
            sub = L[lev::levels]
            xs = [lb["anchor"][0] for lb in sub]
            for i in range(1, len(sub)):
                xs[i] = max(xs[i], xs[i - 1] + (_label_w(sub[i - 1]) + _label_w(sub[i])) / 2 + 12)
            if xs and xs[-1] > xmax:
                xs[-1] = xmax
                for i in range(len(sub) - 2, -1, -1):
                    xs[i] = min(xs[i], xs[i + 1] - (_label_w(sub[i]) + _label_w(sub[i + 1])) / 2 - 12)
            xs = [max(xmin, x) for x in xs]
            if side == "top":
                y = ey0 - 46 - 44 * lev
            else:
                y = (bottom_y if bottom_y is not None else ey1) + 46 + 44 * lev
            for lb, x in zip(sub, xs):
                out.append(dict(lb, pos=(x, y)))
    return out


def _arrow(ax, tip, back, size=13, half=4.2):
    d = np.asarray(tip) - np.asarray(back)
    d = d / (np.linalg.norm(d) + 1e-9)
    n = np.array([-d[1], d[0]])
    tri = np.array([tip, tip - d * size + n * half, tip - d * size - n * half])
    ax.fill(tri[:, 0], tri[:, 1], color="black", lw=0, zorder=6)


def dim(ax, p1, p2, n, d, text, fs=15.5):
    """размер: выносные линии от p1, p2 вдоль n на d пикс., размерная линия со стрелками, текст снаружи"""
    p1, p2, n = np.asarray(p1, float), np.asarray(p2, float), np.asarray(n, float)
    n = n / np.linalg.norm(n)
    a, b = p1 + n * d, p2 + n * d
    lw = 0.8 * 0.72
    for p, q in ((p1, a), (p2, b)):
        e0 = p + n * 6
        ax.plot([e0[0], (q + n * 9)[0]], [e0[1], (q + n * 9)[1]], color="black", lw=lw, zorder=5)
    ax.plot([a[0], b[0]], [a[1], b[1]], color="black", lw=lw, zorder=5)
    _arrow(ax, a, b)
    _arrow(ax, b, a)
    mid = 0.5 * (a + b) + n * 15
    ang = math.degrees(math.atan2(-(b - a)[1], (b - a)[0]))
    if ang > 90:
        ang -= 180
    if ang < -90:
        ang += 180
    ax.text(mid[0], mid[1], text, fontsize=fs, family=FONT, ha="center", va="center", rotation=ang, rotation_mode="anchor", zorder=8.5,
            bbox=dict(facecolor="white", edgecolor="none", pad=1.0))


def _phantom(ax, a, b, lw=0.9):
    ax.plot([a[0], b[0]], [a[1], b[1]], color="black", lw=lw * 0.72, dashes=(14, 3, 1.5, 3, 1.5, 3), zorder=4)


def _centerline(ax, a, b):
    ax.plot([a[0], b[0]], [a[1], b[1]], color="black", lw=0.7 * 0.72, dashes=(18, 3, 2, 3), zorder=4)


def _wrap_chars(width_px, fs_pt):
    """символов в строке для Liberation Serif заданного кегля (средняя ширина знака 0,46 кегля)"""
    return max(20, int(width_px / (0.46 * fs_pt * 100 / 72)))


def _neutral_name(name):
    """«входная камера 2 (двусторонний вход)» -> «входная камера»"""
    return name.split(" (")[0].rstrip("0123456789 ")


def axo_axis_scale(M, s):
    """пикс./м вдоль осей модели X, Y, Z на аксонометрии (масштаб s - в плоскости проекции)"""
    _, right, up = view_basis(M, "axo")
    return [float(np.hypot(e @ right, e @ up) * s) for e in np.eye(3)]


def scale_text(M, view, s, short=False):
    if view != "axo":
        return f"1 м = {s} пикс."
    kx, ky, kz = (int(round(k)) for k in axo_axis_scale(M, s))
    if short:
        return f"1 м = {s} пикс. в плоскости проекции; по осям X, Y, Z — {kx}, {ky}, {kz} пикс./м"
    return (f"1 м = {s} пикс. в плоскости проекции (перпендикулярно взгляду); вдоль осей модели: X — {kx}, Y — {ky}, "
            f"Z — {kz} пикс./м")


def legend_lines(M, view, s, kv_label):
    e, g = M.entry, M.entry["geometry"]
    L = [("h", "Модули замены (номера-выноски, оттенки серого)")]
    for mm in M.modules:
        cnt = ""
        if mm["name"] == "ядро":
            cnt = f", {M.cores} шт."
        if mm["name"] == "модуль рекуператора" and g.get("recuperator"):
            cnt = f", {g['recuperator']['modules']} шт."
        per = "по " if cnt else ""
        mass = f", {per}{fnum(mm['mass_kg'], 1)} кг" if mm.get("mass_kg") is not None else ""
        if not mm["drawn"]:
            tail = " — место и размер не заданы, на эскизе не показан"
        elif all(p["kind"] == "cond" for p in M.parts if p["module"] == mm["no"]):
            tail = " — условная коробка у точки подключения"
        elif mm["no"] in M.enclosed:
            tail = f" — внутри: закрыт снаружи ({M.enclosed[mm['no']]}), на виде сбоку — скрытым контуром (пунктир)"
        else:
            tail = ""
        L.append(("m", f"{mm['no']} — {mm['name']}{cnt}{mass}{tail}", mm))
    neutral = list(dict.fromkeys(_neutral_name(p["name"]) for p in M.parts if p["kind"] == "solid" and not p["module"]))
    if neutral:
        L.append(("m", "без номера — узлы вне модулей замены: " + ", ".join(neutral), dict(no=0, shade=NEUTRAL, drawn=True)))
    L.append(("h", "Точки подключения (буквы)"))
    seen = set()
    dbl = str(g["inlet"]["type"]).startswith("двусторонний")
    for cn in M.conns:
        if cn["letter"] in seen:
            continue
        seen.add(cn["letter"])
        extra = ""
        if cn["letter"].startswith("B"):
            A_in = g["inlet"]["flange_area_m2"]
            n_fl = (2 if dbl else 1) * M.cores
            if n_fl > 1:
                extra = (f"; фланец {fnum(A_in / n_fl, 2)} м² каждый ({n_fl} шт.), всего {fnum(A_in, 2)} м² на двигатель (25 м/с)"
                         + ("; по входу на каждом ядре" if M.cluster else ""))
            else:
                extra = f"; фланец {fnum(A_in, 2)} м² (25 м/с)"
        if cn["letter"] == "C":
            ex = g["exhaust"]
            extra = f"; фланец {fnum(ex['flange_area_m2'], 2)} м² (40 м/с), {fint(ex['T_C'])} °C"
        if cn["letter"] == "D" and M.cluster:
            extra = "; на каждом ядре"
        L.append(("t", f"{cn['letter']} — {cn['kind']}{extra}".replace(" - ", " — ")))
    L.append(("h", "Ступени компрессора, м (на эскизе — точки профиля корпуса)"))
    for st in g["compressor"]:
        if st["type"] == "осевая":
            L.append(("t", f"ст. {st['stage']}, осевая: D {fnum(st['D'])}, x {fnum(st['x0'])}…{fnum(st['x0'] + st['L'])}"))
        else:
            L.append(("t", f"ст. {st['stage']}, {st['type']}: D1t {fnum(st['D1t'])} → D2 {fnum(st['D2'])}, корпус {fnum(st['D_case'])}, "
                           f"x {fnum(st['x0'])}…{fnum(st['x0'] + st['L'])}"))
    L.append(("h", "Обозначения"))
    hid = [t for t, ok in (("корпуса магнитных подшипников", M.hidden), ("модули внутри барабана рекуператора", M.enclosed)) if ok]
    L.append(("t", "сплошная — контур узла по расчётным размерам уровня 1; штрих — условно (размера нет в данных); "
                   "штрихпунктир с двумя точками — габарит с обвязкой (правило компоновки уровня 1); штрихпунктир — ось"
                   + (f"; пунктир — скрытые {' и '.join(hid)} (на виде сбоку)" if hid else "")))
    if M.scale_note:
        L.append(("h", "Масштаб"))
        L.append(("t", "— " + M.scale_note))
    L.append(("h", "Условно в этом варианте"))
    for t in M.cond_notes:
        L.append(("t", "— " + t.replace(" - ", " — ")))
    if M.issues:
        L.append(("h", "Замечания к компоновке уровня 1"))
        seen = set()
        for it in M.issues:
            if it["text"] in seen:
                continue
            seen.add(it["text"])
            L.append(("t", "— " + it["text"].replace(" - ", " — ")))
    return L


def draw_panel(ax, M, view, s, kv_label, P):
    x0, y0, x1, y1 = PANEL
    ax.add_patch(_rect(x0, y0, x1, y1))
    lines = legend_lines(M, view, s, kv_label)
    if view == "axo":
        tri = np.array([[v @ P.right, -(v @ P.up)] for v in np.eye(3)]) * s
        bottom_reserved = (max(0.0, -tri[:, 1].min()) + 146)
    else:
        tri = None
        bottom_reserved = 150
    fs_b, fs_h = 13.0, 14.5
    lh_b, lh_h = 1.30, 2.0
    for _ in range(10):
        rows = []
        for ln in lines:
            if ln[0] == "h":
                rows.append(("h", ln[1], None, 0))
                continue
            ind = 34 if ln[0] == "m" else 0
            wrapped = textwrap.wrap(ln[1], _wrap_chars(x1 - x0 - 40 - ind - 18, fs_b)) or [""]
            for k, w in enumerate(wrapped):
                rows.append((ln[0] if k == 0 else "c", w, ln[2] if ln[0] == "m" else None, ind if k == 0 else ind + 18))
        hpx = sum((fs_h * 1.389 * lh_h if r[0] == "h" else fs_b * 1.389 * lh_b) for r in rows)
        if hpx <= (y1 - y0) - bottom_reserved - 30:
            break
        fs_b *= 0.94
        fs_h *= 0.95
    y = y0 + 10
    for kind, text, mm, ind in rows:
        if kind == "h":
            y += fs_h * 1.389 * 0.8
            ax.text(x0 + 18, y, text, fontsize=fs_h, family=FONT, va="top", ha="left")
            y += fs_h * 1.389 * 1.2
            continue
        if kind == "m":
            sw = fs_b * 1.389 * 0.8
            cond_only = bool(mm.get("no")) and mm["drawn"] and all(p["kind"] == "cond" for p in M.parts if p["module"] == mm["no"])
            sq = _rect(x0 + 20, y + 3, x0 + 20 + sw, y + 3 + sw, fc=mm["shade"] if mm["drawn"] and not cond_only else "white", lw=0.6)
            if cond_only:
                sq.set_linestyle((0, (2.5, 1.5)))
            ax.add_patch(sq)
        ax.text(x0 + 18 + ind, y, text, fontsize=fs_b, family=FONT, va="top", ha="left")
        y += fs_b * 1.389 * lh_b
    yb = y1 - bottom_reserved + 20
    if tri is not None:
        kx, ky, kz = (int(round(k)) for k in axo_axis_scale(M, s))
        ax.text(x0 + 18, yb - 26, f"Масштаб единый для всех листов набора: 1 м = {s} пикс. в плоскости проекции;", fontsize=13,
                family=FONT, va="top")
        ax.text(x0 + 18, yb, f"вдоль осей модели X — {kx}, Y — {ky}, Z — {kz} пикс./м (аксонометрия)", fontsize=13, family=FONT, va="top")
        ax.text(x0 + 18, yb + 26, "Отрезки 1 м вдоль осей модели:", fontsize=12.5, family=FONT, va="top")
    else:
        ax.text(x0 + 18, yb, f"Масштаб единый для всех листов набора: 1 м = {s} пикс.", fontsize=13, family=FONT, va="top")
    if tri is not None:
        o = np.array([x0 + 110 + max(0.0, -tri[:, 0].min()), y1 - 40.0])      # место под подпись «Y: 1 м» слева
        for d, nm in zip(tri, ("X", "Y", "Z")):
            e_ = o + d
            ax.plot([o[0], e_[0]], [o[1], e_[1]], color="black", lw=1.6 * 0.72, zorder=6)
            m_ = o + 0.5 * d
            ax.plot([m_[0]], [m_[1]], marker="o", ms=2.5, color="black")
            ax.text(e_[0] + (8 if d[0] >= 0 else -8), e_[1] - 6, f"{nm}: 1 м", fontsize=12.5, family=FONT, va="center",
                    ha="left" if d[0] >= 0 else "right")
    else:
        ax.text(x0 + 18, yb + 26, "Масштабная линейка:", fontsize=12.5, family=FONT, va="top")
        bx, by = x0 + 30, yb + 60
        ax.add_patch(_rect(bx, by, bx + s / 2, by + 9, fc="black", lw=0.8))
        ax.add_patch(_rect(bx + s / 2, by, bx + s, by + 9, fc="white", lw=0.8))
        for t, lab in ((0, "0"), (0.5, "0,5"), (1.0, "1 м")):
            ax.text(bx + s * t, by + 18, lab, fontsize=12.5, family=FONT, ha="center", va="top")


def _rect(x0, y0, x1, y1, fc="none", lw=0.8):
    from matplotlib.patches import Rectangle
    return Rectangle((x0, y0), x1 - x0, y1 - y0, fill=fc != "none", facecolor=fc if fc != "none" else "none", edgecolor="black",
                     lw=lw * 0.72, zorder=3)


def draw_title(ax, M, view, s, sheet_no, kv_label, n_total, src):
    e = M.entry
    x0, y0, x1, y1 = TITLE
    xs = TITLE_SPLIT
    ax.add_patch(_rect(x0, y0, x1, y1))
    ax.plot([xs, xs], [y0, y1], color="black", lw=0.8 * 0.72)
    ax.text(x0 + 18, y0 + 16, CAPTION, fontsize=19, family=FONT, va="top")
    rank = e.get("div_rank")
    head = (f"Вариант № {rank} из {n_total}" if rank else f"Фикстура теста «{e.get('fixture_text', e.get('fixture_class', ''))}»")
    ax.text(x0 + 18, y0 + 52, f"{head} · {kv_label} · idx {e.get('idx')} · {e.get('archetype', e.get('arch', ''))}",
            fontsize=14.5, family=FONT, va="top")
    wl = _wrap_chars(xs - x0 - 40, 13)
    sch = textwrap.wrap("Схема: " + str(e.get("scheme", "")), wl)
    for k, w in enumerate(sch[:2]):
        ax.text(x0 + 18, y0 + 84 + 23 * k, w, fontsize=13, family=FONT, va="top")
    d = dims_of(M)
    eta = e.get("eta_th_shaft_pct")
    sc = e.get("score")
    tail = (f"Габарит L × W × H = {fnum(d['json']['L'])} × {fnum(d['json']['W'])} × {fnum(d['json']['H'])} м (оценка); "
            f"масса {fint(e.get('mass_kg_estimate') or e['geometry'].get('mass_kg_estimate'))} кг (оценка)"
            + (f"; тепловой КПД на валу {fnum(eta, 1)} % (ISO, оценка)" if eta is not None else "")
            + (f"; взвешенная оценка {fnum(sc, 2)}" if sc is not None else "")
            + f"; высота оси выходного вала h = {fnum(M.z_ax)} м — принято 0,5·H (высота оси пакета — нет данных)")
    for k, w in enumerate(textwrap.wrap(tail, wl)[:3]):
        ax.text(x0 + 18, y0 + 136 + 23 * k, w, fontsize=13, family=FONT, va="top")
    ax.text(xs + 18, y0 + 16, f"Лист {sheet_no} из 3", fontsize=17, family=FONT, va="top")
    wr = _wrap_chars(x1 - xs - 40, 13.5)
    for k, w in enumerate(textwrap.wrap(VIEW_TITLE[view], wr)):
        ax.text(xs + 18, y0 + 50 + 23 * k, w, fontsize=13.5, family=FONT, va="top")
    sc_lines = textwrap.wrap(f"Масштаб (один для всего набора): {scale_text(M, view, s, short=True)}", _wrap_chars(x1 - xs - 40, 12.5))
    for k, w in enumerate(sc_lines):
        ax.text(xs + 18, y0 + 104 + 21 * k, w, fontsize=12.5, family=FONT, va="top")
    ys = y0 + 104 + 21 * len(sc_lines) + 8
    for k, w in enumerate(textwrap.wrap(f"Источник: {src}, поле geometry (funnel_l1.layout), {DATE}; модель и рендер — funnel_sketch.py",
                                        _wrap_chars(x1 - xs - 40, 12))):
        ax.text(xs + 18, ys + 21 * k, w, fontsize=12, family=FONT, va="top")


def alt_text(M, view):
    """текстовая замена листа (для сайта и читалок экрана)"""
    e, d = M.entry, dims_of(M)
    who = f"вариант № {e['div_rank']}" if e.get("div_rank") else f"фикстура «{e.get('fixture_class')}»"
    mods = "; ".join(f"{m['no']} — {m['name']}" for m in M.modules if m["drawn"])
    return (f"Эскиз, {VIEW_TITLE[view].lower()}: {who}, {e.get('archetype', e.get('arch', ''))}. Габарит с обвязкой "
            f"{fnum(d['json']['L'])} × {fnum(d['json']['W'])} × {fnum(d['json']['H'])} м, масса "
            f"{fint(e.get('mass_kg_estimate') or e['geometry'].get('mass_kg_estimate'))} кг (оценка уровня 1). Модули замены: {mods}. "
            f"Условно (штрихом): {', '.join(sorted({p['name'] for p in M.parts if p['kind'] == 'cond'}))}."
            + (f" Замечания: {'; '.join(dict.fromkeys(it['text'] for it in M.issues))}." if M.issues else ""))


def render_sheet(job):
    """один лист: вид модели + размеры, выноски, легенда, штамп -> PNG"""
    M, view, s, path, kv_label, n_total, src = job
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    img, depth, P, cam = render_image(M, view, s)
    fig = plt.figure(figsize=(SHEET_W / 100, SHEET_H / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, SHEET_W)
    ax.set_ylim(SHEET_H, 0)
    ax.axis("off")
    ax.imshow(img, extent=(AREA[0], AREA[2], AREA[3], AREA[1]), interpolation="none", zorder=0)
    ax.add_patch(_rect(20, 20, SHEET_W - 20, SHEET_H - 20, lw=1.4))
    lo, hi = M.env["lo"], M.env["hi"]
    C = env_corners(M)
    S = P(C)
    ex0, ey0 = S.min(0)
    ex1, ey1 = S.max(0)
    Sm = P(bbox_corners(M))                     # с выступающими узлами - для полок выносок
    rx0, ry0 = Sm.min(0)
    rx1, ry1 = Sm.max(0)
    # габарит с обвязкой
    E = [(a, b) for a in range(8) for b in range(a + 1, 8) if np.sum(np.abs(C[a] - C[b]) > 1e-9) == 1]
    far = int(np.argmax(P.dist(C, cam)))
    for a, b in E:
        if view == "axo" and far in (a, b):
            continue
        _phantom(ax, S[a], S[b])
    # оси
    pad = 0.12
    for bd in M.bodies:
        segs = [s_ for s_ in bd["segs"]]
        if not segs:
            continue
        xa = min(s_[0] for s_ in segs) - pad
        xb = max(s_[1] for s_ in segs) + pad
        if bd is M.bodies[M.main]:
            xa, xb = min(xa, lo[0] - pad), max(xb, hi[0] + pad)
        if view == "end":
            c = P([xa, bd["y"], bd["z"]])[0]
            _centerline(ax, c - [22, 0], c + [22, 0])
            _centerline(ax, c - [0, 22], c + [0, 22])
        else:
            a_, b_ = P([[xa, bd["y"], bd["z"]], [xb, bd["y"], bd["z"]]])
            _centerline(ax, a_, b_)
    # скрытые корпуса магнитных подшипников
    if view == "side":
        for hb in M.hidden:
            q = P([[hb["x"] - hb["L"] / 2, hb["y"], hb["z"] - hb["D"] / 2], [hb["x"] + hb["L"] / 2, hb["y"], hb["z"] + hb["D"] / 2]])
            xa, ya = q.min(0)
            xb, yb = q.max(0)
            ax.plot([xa, xb, xb, xa, xa], [ya, ya, yb, yb, ya], color="black", lw=0.9 * 0.72, dashes=(2, 2.5), zorder=4)
        for pl_ in hidden_outlines(M):           # модули внутри барабана рекуператора
            q = P(pl_)
            ax.plot(q[:, 0], q[:, 1], color="black", lw=0.9 * 0.72, dashes=(2, 2.5), zorder=4)
    # размеры
    out_right = (P(M.shaft_end)[0][0] > 0.5 * (ex0 + ex1))
    if view == "side":
        labels = anchors(M, P, cam, depth, view)
        placed = place_rails(labels, rx0, ry0, rx1, ry1)
        nb = max([lb["pos"][1] for lb in placed if lb["pos"][1] > ey1] + [ry1])
        dim(ax, (ex0, ey1), (ex1, ey1), (0, 1), nb - ey1 + 52, f"L = {fnum(hi[0] - lo[0])} м")
        hs = P([lo[0], 0, 0])[0][0] if out_right else P([hi[0], 0, 0])[0][0]
        dim(ax, (hs, ey1), (hs, ey0), (-1, 0) if out_right else (1, 0), 62, f"H = {fnum(hi[2] - lo[2])} м")
        sp = P(M.shaft_end)[0]
        dim(ax, (sp[0], ey1), sp, (1, 0) if out_right else (-1, 0), 62, f"h = {fnum(M.z_ax)} м")
    elif view == "end":
        labels = anchors(M, P, cam, depth, view)
        placed = place_rails(labels, rx0, ry0, rx1, ry1)
        nb = max([lb["pos"][1] for lb in placed if lb["pos"][1] > ey1] + [ry1])
        dim(ax, (ex0, ey1), (ex1, ey1), (0, 1), nb - ey1 + 52, f"W = {fnum(hi[1] - lo[1])} м")
        dim(ax, (ex0, ey1), (ex0, ey0), (-1, 0), 62, f"H = {fnum(hi[2] - lo[2])} м")
        sp = P(M.shaft_end)[0]
        dim(ax, (sp[0], ey1), sp, (1, 0), (ex1 - sp[0]) + 62, f"h = {fnum(M.z_ax)} м")
    else:
        c2 = S.mean(0)

        def edge_dim(i, j, text):
            a_, b_ = S[i], S[j]
            t = b_ - a_
            n = np.array([-t[1], t[0]])
            n /= np.linalg.norm(n)
            if np.dot(n, 0.5 * (a_ + b_) - c2) < 0:
                n = -n
            dim(ax, a_, b_, n, 44, text)

        def ci(x, y, z):
            return int(np.argmin(np.linalg.norm(C - np.array([x, y, z]), axis=1)))
        edge_dim(ci(lo[0], lo[1], 0), ci(hi[0], lo[1], 0), f"L = {fnum(hi[0] - lo[0])} м")
        edge_dim(ci(lo[0], lo[1], 0), ci(lo[0], hi[1], 0), f"W = {fnum(hi[1] - lo[1])} м")
        vert = [(ci(x, y, 0), ci(x, y, hi[2])) for x in (lo[0], hi[0]) for y in (lo[1], hi[1])]
        a_, b_ = min(vert, key=lambda ij: S[ij[0]][0])
        dim(ax, S[a_], S[b_], (-1, 0), 50, f"H = {fnum(hi[2] - lo[2])} м")
        labels = anchors(M, P, cam, depth, view)
        placed = place_rails(labels, rx0, ry0, rx1, ry1, bottom_y=max(ey1, ry1) + 60)
    for lb in placed:
        x, y = lb["pos"]
        a = lb["anchor"]
        w = _label_w(lb)
        r = 17
        d = np.array(a) - np.array([x, y])
        dn = d / (np.linalg.norm(d) + 1e-9)
        if lb["kind"] == "mod":
            start = np.array([x, y]) + dn * r
        else:
            hw, hh = w / 2, 15
            t = min(hw / (abs(dn[0]) + 1e-9), hh / (abs(dn[1]) + 1e-9))
            start = np.array([x, y]) + dn * t
        ax.plot([start[0], a[0]], [start[1], a[1]], color="black", lw=0.75 * 0.72, zorder=8)
        ax.plot([a[0]], [a[1]], marker="o", ms=3.6, color="black", zorder=9)
        if lb["kind"] == "mod":
            from matplotlib.patches import Circle
            ax.add_patch(Circle((x, y), r, facecolor="white", edgecolor="black", lw=0.9 * 0.72, zorder=9))
        else:
            ax.add_patch(_rect(x - w / 2, y - 15, x + w / 2, y + 15, fc="white", lw=0.9))
        ax.text(x, y + 1, lb["text"], fontsize=15.5, family=FONT, ha="center", va="center", zorder=10)
    draw_panel(ax, M, view, s, kv_label, P)
    draw_title(ax, M, view, s, VIEWS.index(view) + 1, kv_label, n_total, src)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=100, metadata={"Title": f"{CAPTION}: {VIEW_TITLE[view]}", "Description": alt_text(M, view),
                                         "Software": "funnel_sketch.py"})
    plt.close(fig)
    return path


# =========================================================================== набор
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def code_sha():
    return sha256(os.path.abspath(__file__))[:16]


def kv_label_of(top, results):
    kv = top.get("kv") or ""
    t = (results or {}).get("target") or {}
    if t.get("model"):
        return f"{kv} ({t['model']}" + (f", {t['variant']}" if t.get("variant") else "") + ")"
    return kv


def tag_of(e):
    return f"{int(e['div_rank']):02d}" if e.get("div_rank") else f"fx_{e.get('fixture_class', e.get('idx'))}"


def sheet_paths(out, e):
    t = tag_of(e)
    return {v: os.path.join(out, t, f"{t}_{VIEW_FILE[v]}.png") for v in VIEWS}


def _r4(x):
    return [round(float(v), 4) for v in np.atleast_1d(x)]


def _node_rec(p):
    """размеры узла так, как он нарисован: тело вращения по оси x - профиль [x, D]; иначе - ось и размеры или углы коробки"""
    pr = p["prim"]
    rec = dict(name=p["name"], module=p["module"])
    if pr[0] == "rev" and p["body"] is not None:
        rec["profile_x_D"] = [_r4([s, 2 * r]) for s, r in pr[1]]
    elif pr[0] == "rev":
        rec.update(p0=_r4(pr[2]), axis=_r4(pr[3]), L=round(float(pr[1][-1][0] - pr[1][0][0]), 4), D=round(float(2 * max(r for _, r in pr[1])), 4))
    else:
        rec.update(lo=_r4(pr[1]), hi=_r4(pr[2]))
    return rec


def wheels_drawn(M):
    """колёса компрессора по нарисованному профилю (ядро 1): D1t, D2, корпус - для сверки эскиза с geometry.compressor"""
    out = []
    for st in M.entry["geometry"]["compressor"]:
        p = next(q for q in M.parts if q["inst"] == 0 and q["name"].startswith(f"ступень компрессора {st['stage']} "))
        prof = p["prim"][1]
        rec = dict(stage=st["stage"], type=st["type"], x0=round(float(prof[0][0]), 4), x1=round(float(prof[-1][0]), 4))
        if st["type"] == "осевая":
            rec["D"] = round(float(2 * prof[0][1]), 4)
        else:
            rec.update(D1t=round(float(2 * prof[0][1]), 4), D2=round(float(2 * prof[1][1]), 4), D_case=round(float(2 * prof[2][1]), 4))
        out.append(rec)
    return out


def dims_record(M, out_rel=None):
    e = M.entry
    d = dims_of(M)
    rec = dict(div_rank=e.get("div_rank"), idx=e.get("idx"), arch=e.get("arch"), archetype=e.get("archetype"),
               niche=e.get("niche"), scheme=e.get("scheme"),
               envelope_json_m=d["json"], envelope_model_m=d["model"], envelope_rule_m=d["rule"], nodes_solid_bbox_m=d["nodes_solid"],
               dev_model_vs_json_pct=d["dev_pct"], dev_rule_vs_json_pct=d["dev_rule_pct"], envelope_ok_1pct=d["ok"],
               mass_kg_estimate=e.get("mass_kg_estimate") or e["geometry"].get("mass_kg_estimate"),
               mass_basis="funnel_l1 (оценка: массы узлов × L0.MASS_K)",
               axis_height_m=round(M.z_ax, 4), axis_height_basis="0,5·H (funnel_l1.layout: высота оси пакета - нет данных, запрос заказчику)",
               output_end=e["geometry"]["output_shaft"]["output_end"],
               cores=M.cores, recuperator=(e["geometry"].get("recuperator") or {}).get("shape"),
               inlet=e["geometry"]["inlet"]["type"], combustor=e["geometry"]["combustor"]["type"],
               solid_nodes=[p["name"] for p in M.parts if p["kind"] == "solid"],
               conditional_nodes=sorted({p["name"] for p in M.parts if p["kind"] == "cond"}),
               conditional_notes=M.cond_notes,
               hidden_nodes=sorted({h["name"] for h in M.hidden}),
               modules=[dict(no=m["no"], name=m["name"], mass_kg=m["mass_kg"], drawn=m["drawn"],
                             enclosed_by=M.enclosed.get(m["no"])) for m in M.modules],
               compressor_wheels_drawn_m=wheels_drawn(M),
               drawn_nodes_m=[_node_rec(p) for p in M.parts if p["kind"] == "solid" and p["inst"] == 0],
               connections=[dict(letter=c["letter"], kind=c["kind"], x_m=round(float(c["x"]), 4)) for c in M.conns],
               issues=[{k: v for k, v in it.items()} for it in M.issues])
    rec["alt_ru"] = {v: alt_text(M, v) for v in VIEWS}
    if out_rel is not None:
        rec["sheets"] = out_rel
    return rec


def load_entries(path):
    top = json.load(open(path, encoding="utf-8"))
    if "selected" in top:
        return top, top["selected"]
    return top, top["entries"]


class ScaleOverflow(ValueError):
    """заданный --px-per-m больше, чем помещается у части вариантов"""


def run(input_path, out, dims_path=None, results_path=None, render=True, px_per_m=None, workers=1, only=None, allow_overflow=False):
    top, entries = load_entries(input_path)
    results = json.load(open(results_path, encoding="utf-8")) if results_path and os.path.exists(results_path) else None
    kv_label = kv_label_of(top, results)
    models, skipped = [], []
    for e in entries:
        if only and e.get("div_rank") not in only and tag_of(e) not in only:
            continue
        try:
            models.append(build_model(e))
        except (ValueError, KeyError, TypeError) as ex:
            skipped.append(dict(div_rank=e.get("div_rank"), idx=e.get("idx"), reason=f"{type(ex).__name__}: {ex}"))
    s = int(px_per_m) if px_per_m else common_scale(models)
    # заданный масштаб проверяется: вид, не входящий в поле, режется и наезжает на легенду - без --allow-overflow это ошибка
    overflow = [dict(div_rank=M.entry.get("div_rank"), idx=M.entry.get("idx"), fixture_class=M.entry.get("fixture_class"),
                     own_max_px_per_m=round(max_scale(M), 1)) for M in models if max_scale(M) < s]
    if overflow and not allow_overflow:
        raise ScaleOverflow(f"масштаб 1 м = {s} пикс. не помещает виды у {len(overflow)} вариантов (наибольший общий - "
                            f"{common_scale(models)} пикс./м): {overflow}; уменьшите --px-per-m или добавьте --allow-overflow")
    for M in models:
        M.scale_note = None
        if max_scale(M) < s:
            M.scale_note = (f"масштаб набора {s} пикс./м больше наибольшего для этого варианта ({fnum(max_scale(M), 1)} пикс./м): "
                            f"вид может выходить за поле и наезжать на легенду (--allow-overflow)")
    n_total = len(entries)
    src = os.path.basename(input_path)
    jobs, recs = [], []
    for M in models:
        paths = sheet_paths(out, M.entry)
        recs.append(dims_record(M, {v: os.path.relpath(p, out) for v, p in paths.items()}))
        for v in VIEWS:
            jobs.append((M, v, s, paths[v], kv_label, n_total, src))
    done = []
    if render:
        if workers > 1 and len(jobs) > 1:
            import multiprocessing as mp
            with mp.get_context("spawn").Pool(workers) as pool:
                done = pool.map(render_sheet, jobs, chunksize=1)
        else:
            done = [render_sheet(j) for j in jobs]
    res = dict(title=f"Эскизы вариантов уровня 1: габариты, масса, условные узлы - {top.get('kv', '')}",
               method=("упрощённая трёхмерная модель внешнего вида по полю geometry (funnel_l1.layout): тела вращения и параллелепипеды по "
                       "размерам узлов; узел без размера - условно (штрих), список - conditional_nodes; габарит модели собран из размеров "
                       "узлов по правилу компоновки уровня 1 и сверен с полем envelope (допуск 1 %); три вида в одном масштабе: вид сбоку, "
                       "вид с торца со стороны выходного вала (ортогональные), аксонометрия 3/4 (параллельная); рендер pyvista (OSMesa), "
                       "надписи matplotlib"),
               software=_software(render), date=DATE, code_sha=code_sha(),
               inputs=[dict(file=os.path.relpath(input_path, FUNNEL) if input_path.startswith(FUNNEL) else os.path.basename(input_path),
                            sha256=sha256(input_path))] + ([dict(file=os.path.basename(results_path), sha256=sha256(results_path),
                                                                  used="подпись расчётной точки (target)")] if results else []),
               kv=top.get("kv"), kv_label=kv_label, n_entries=n_total, n_drawn=len(models), skipped=skipped,
               scale_px_per_m=s, scale_forced=bool(px_per_m), scale_overflow=overflow,
               axo_axis_px_per_m=(dict(zip("XYZ", (round(k, 1) for k in axo_axis_scale(models[0], s)))) if models else None),
               sheet_px=[SHEET_W, SHEET_H],
               scale_limited_by=(None if not models else (lambda M: dict(div_rank=M.entry.get("div_rank"), idx=M.entry.get("idx"),
                                                                       arch=M.entry.get("arch"), own_max_px_per_m=round(max_scale(M), 1)))(
                   min(models, key=max_scale))), views=[dict(key=v, title=VIEW_TITLE[v], file_suffix=VIEW_FILE[v]) for v in VIEWS],
               axo=dict(azimuth_deg=AXO_AZ, elevation_deg=AXO_EL, projection="параллельная"),
               tolerance_envelope=TOL_ENV, caption=CAPTION,
               note=("все размеры - оценка уровня 1; полноразмерные PNG - вне репозитория (папка --out), в репозитории - три примера листа "
                     "(examples/); высота оси выходного вала - 0,5·H, данных пакета нет"),
               variants=recs, rendered=len(done))
    if dims_path:
        with open(dims_path, "w", encoding="utf-8") as f:
            json.dump(res, f, ensure_ascii=False, indent=1)
    return res, models


def _software(render):
    parts = [f"python {sys.version.split()[0]}", f"numpy {np.__version__}"]
    if render:
        try:
            import pyvista, vtk, matplotlib
            parts += [f"pyvista {pyvista.__version__}", f"vtk {vtk.vtkVersion.GetVTKVersion()}", f"matplotlib {matplotlib.__version__}"]
        except ImportError:
            pass
    return ", ".join(parts)


def shrink_png(src, dst, limit=200_000):
    """пример листа для репозитория: оттенки серого, палитра, при необходимости - уменьшение до limit байт"""
    from PIL import Image
    im = Image.open(src).convert("L")
    scale = 1.0
    while True:
        w, h = int(im.width * scale), int(im.height * scale)
        q = im.resize((w, h), Image.LANCZOS) if scale < 1.0 else im
        q = q.quantize(colors=24, dither=Image.Dither.NONE)
        q.save(dst, optimize=True)
        if os.path.getsize(dst) <= limit or scale < 0.3:
            return dict(file=os.path.basename(dst), px=[w, h], bytes=os.path.getsize(dst))
        scale *= 0.9


def main(argv=None):
    ap = argparse.ArgumentParser(description="эскизы внешнего вида вариантов уровня 1 воронки (три вида, один масштаб)")
    ap.add_argument("--input", default=os.path.join(FUNNEL, "funnel_l1_top50.json"))
    ap.add_argument("--results", default=os.path.join(FUNNEL, "funnel_l1_results.json"))
    ap.add_argument("--out", default=os.path.join(os.environ.get("GTU_OUT", os.path.join(HERE, "_out")), "sketch"))
    ap.add_argument("--dims", default=os.path.join(HERE, "sketch_dims.json"), help="результат: габариты и масса, условные узлы ('' - не писать)")
    ap.add_argument("--px-per-m", type=float, default=None, help="масштаб (по умолчанию - общий наибольший для набора); "
                                                                 "больше наибольшего у какого-либо варианта - ошибка")
    ap.add_argument("--allow-overflow", action="store_true", help="рисовать и при заданном масштабе, не помещающем виды (пометка на листе)")
    ap.add_argument("--workers", type=int, default=int(os.environ.get("NPROC", "0") or 0) or min(4, os.cpu_count() or 1))
    ap.add_argument("--only", default="", help="номера вариантов через запятую (div_rank) или метки фикстур")
    ap.add_argument("--no-render", action="store_true", help="только модель и sketch_dims.json")
    ap.add_argument("--examples", default="", help="папка примеров для репозитория (листы <=200 КБ)")
    ap.add_argument("--example-sheets", default="1:side,1:axo,2:axo", help="какие листы в примеры: номер:вид")
    a = ap.parse_args(argv)
    only = None
    if a.only:
        only = [int(t) if t.isdigit() else t for t in a.only.split(",")]
    try:
        res, models = run(a.input, a.out, a.dims or None, a.results, render=not a.no_render, px_per_m=a.px_per_m, workers=a.workers,
                          only=only, allow_overflow=a.allow_overflow)
    except ScaleOverflow as ex:
        print(f"ОШИБКА: {ex}", file=sys.stderr)
        sys.exit(2)
    print(f"вариантов {res['n_entries']}, нарисовано {res['n_drawn']}, листов {res['rendered']}, масштаб 1 м = {res['scale_px_per_m']} пикс.")
    for ov in res["scale_overflow"]:
        print(f"  ВНИМАНИЕ: масштаб не помещает виды варианта {ov} (--allow-overflow)", file=sys.stderr)
    for r in res["variants"]:
        flag = "ок" if r["envelope_ok_1pct"] else "РАСХОЖДЕНИЕ"
        print(f"  {r['div_rank'] or r['idx']}: {r['arch']}: L×W×H JSON {r['envelope_json_m']} модель {r['envelope_model_m']} "
              f"отклонение % {r['dev_model_vs_json_pct']} {flag}; масса {r['mass_kg_estimate']} кг; условных узлов "
              f"{len(r['conditional_nodes'])}; замечаний {len(r['issues'])}")
    for sk in res["skipped"]:
        print(f"  пропущен {sk}")
    if a.examples and not a.no_render:
        os.makedirs(a.examples, exist_ok=True)
        ex = []
        for item in a.example_sheets.split(","):
            rk, v = item.split(":")
            e = next((M.entry for M in models if str(M.entry.get("div_rank")) == rk), None)
            if e is None:
                continue
            src = sheet_paths(a.out, e)[v]
            dst = os.path.join(a.examples, f"sketch_{tag_of(e)}_{VIEW_FILE[v]}.png")
            ex.append(shrink_png(src, dst))
        print("примеры:", ex)


if __name__ == "__main__":
    main()
