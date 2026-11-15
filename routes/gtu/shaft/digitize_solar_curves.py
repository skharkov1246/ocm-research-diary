"""Оцифровка графика «Available Power» из открытых листов данных Solar Taurus 60 (GS и MD).

Зачем. Паспорт даёт одну точку ISO; модель вала калибруется по ней. Для проверки вне калибровки
нужны другие точки — лист данных изготовителя содержит векторный график полной мощности по
температуре наружного воздуха (−30…+45 °C). Это открытый документ (тот же, что в fleet/ как источник).

Как: страница 2 PDF -> SVG (pdftocairo из poppler-utils); рамка графика (прямоугольник с Z) даёт оси:
левый край — −30 °C, правый — +45 °C, низ и верх — подписи шкалы мощности (заданы ниже, сверены с
картинкой страницы); кривая мощности — первый подпуть толстого пути (Безье) — вычисляется точно.
Скрипт ничего не скачивает и не выполняет из PDF.

Повторить (из корня репозитория):
    d=$(mktemp -d)
    curl -sS -o "$d"/gs.pdf https://s7d2.scene7.com/is/content/Caterpillar/CM20150703-52095-49890
    curl -sS -o "$d"/md.pdf https://s7d2.scene7.com/is/content/Caterpillar/C10550209
    python3 -I routes/gtu/shaft/digitize_solar_curves.py "$d"/gs.pdf "$d"/md.pdf
Результат — data/solar_available_power_taurus60.json рядом (каталог data/ узла).
"""
import hashlib, json, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SHEETS = {   # шкалы — подписи на графике (страница 2), сверены с картинкой страницы
    "GS": dict(doc="DS60PG/0320, Taurus 60 Gas Turbine Generator Set",
               url="https://s7d2.scene7.com/is/content/Caterpillar/CM20150703-52095-49890",
               chart="DS60GS-002M", y_axis="OUTPUT POWER, kWe", y_bottom=4000.0, y_top=8000.0, unit="кВт эл. (клеммы)",
               basis="Nominal rating – per ISO; at 15°C, sea level; no inlet/exhaust losses; relative humidity 60%; "
                     "natural gas LHV = 35 MJ/Nm3; no accessory losses"),
    "MD": dict(doc="DS60MD/713, Taurus 60 Gas Turbine Mechanical-Drive Package",
               url="https://s7d2.scene7.com/is/content/Caterpillar/C10550209",
               chart="DS60MD-002M", y_axis="OUTPUT POWER, kW (hp)", y_bottom=4500.0, y_top=6500.0, unit="кВт (вал)",
               basis="Nominal Rating – ISO; at 15°C, sea level; no inlet/exhaust losses; relative humidity 60%; "
                     "natural gas LHV = 35 MJ/nm3; optimum power turbine speed; AC-driven accessories"),
}
X_LEFT_C, X_RIGHT_C = -30.0, 45.0
T_EVAL = [-30.0, -22.5, -15.0, -7.5, 0.0, 7.5, 15.0, 22.5, 30.0, 37.5, 45.0]


def _num(s):
    return [float(v) for v in re.findall(r"-?\d+(?:\.\d+)?(?:e-?\d+)?", s)]


def _mat(p):
    t = re.search(r'transform="matrix\(([^)]+)\)"', p)
    return _num(t.group(1)) if t else [1, 0, 0, 1, 0, 0]


def _ap(m, x, y):
    return m[0] * x + m[2] * y + m[4], m[1] * x + m[3] * y + m[5]


def parse_svg(svg):
    frame, curve = None, None
    for p in re.findall(r"<path[^>]*>", svg):
        d = re.search(r' d="([^"]+)"', p)
        sw = re.search(r'stroke-width="([\d.]+)"', p)
        if not d or not sw:
            continue
        d, sw, m = d.group(1), float(sw.group(1)), _mat(p)
        if "Z" in d and 0.9 <= sw <= 1.1 and frame is None:
            pts = [_ap(m, *xy) for xy in zip(_num(d)[0::2], _num(d)[1::2])]
            xs, ys = [q[0] for q in pts], [q[1] for q in pts]
            if max(xs) - min(xs) > 200 and max(ys) - min(ys) > 150:
                frame = (min(xs), min(ys), max(xs), max(ys))
        if "C" in d and sw >= 2.5 and curve is None:
            sub = d.split("M")[1]                       # первый подпуть — кривая мощности
            toks = re.findall(r"[MC]|-?\d+(?:\.\d+)?", "M" + sub)
            nums = [float(t) for t in toks if t not in ("M", "C")]
            pts = [_ap(m, nums[i], nums[i + 1]) for i in range(0, len(nums), 2)]
            curve = pts                                 # P0, (c1, c2, P1), (c1, c2, P2) ...
    return frame, curve


def bezier_samples(pts, k=400):
    out = []
    for s in range(0, len(pts) - 1, 3):
        p0, c1, c2, p1 = pts[s], pts[s + 1], pts[s + 2], pts[s + 3]
        for i in range(k + 1):
            t = i / k
            a, b, c, dd = (1 - t) ** 3, 3 * t * (1 - t) ** 2, 3 * t * t * (1 - t), t ** 3
            out.append((a * p0[0] + b * c1[0] + c * c2[0] + dd * p1[0], a * p0[1] + b * c1[1] + c * c2[1] + dd * p1[1]))
    return out


def digitize(pdf, key):
    sh = SHEETS[key]
    with tempfile.TemporaryDirectory() as td:
        base = os.path.join(td, "p2.svg")
        subprocess.run(["pdftocairo", "-svg", "-f", "2", "-l", "2", pdf, base], check=True)
        frame, curve = parse_svg(open(base).read())
    x0, y0, x1, y1 = frame                               # y растёт вниз
    toT = lambda x: X_LEFT_C + (x - x0) / (x1 - x0) * (X_RIGHT_C - X_LEFT_C)
    toP = lambda y: sh["y_bottom"] + (y1 - y) / (y1 - y0) * (sh["y_top"] - sh["y_bottom"])
    s = [(toT(x), toP(y)) for x, y in bezier_samples(curve)]
    s.sort()
    Tmin, Tmax = s[0][0], s[-1][0]
    pts = []
    for T in T_EVAL:
        if Tmax < T <= Tmax + 0.5:      # кривая кончается чуть раньше края рамки: продолжение по последнему 1 °C
            (ta, pa) = next(q for q in s if q[0] >= Tmax - 1.0)
            pb = s[-1][1]
            pts.append(dict(T_amb_C=T, P_kW=round(pb + (pb - pa) / (Tmax - ta) * (T - Tmax), 1), extrapolated_C=round(T - Tmax, 2)))
            continue
        if T < Tmin - 1e-6 or T > Tmax + 1e-6:
            continue
        for (ta, pa), (tb, pb) in zip(s, s[1:]):
            if ta <= T <= tb:
                pts.append(dict(T_amb_C=T, P_kW=round(pa + (pb - pa) * (T - ta) / (tb - ta) if tb > ta else pa, 1)))
                break
    ppt = (sh["y_top"] - sh["y_bottom"]) / (y1 - y0)
    return dict(sheet=sh, pdf_sha256=hashlib.sha256(open(pdf, "rb").read()).hexdigest(),
                frame_pt=[round(v, 3) for v in frame], curve_T_range_C=[round(Tmin, 2), round(Tmax, 2)],
                kW_per_pt=round(ppt, 3), points=pts)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    out = dict(method="оцифровка векторного графика «Available Power»: рамка графика -> оси (−30…+45 °C; шкала "
                      "мощности по подписям), кривая — точное вычисление Безье первого подпути толстого пути "
                      "(pdftocairo -svg). Точность ограничена самим графиком: это сглаженная кривая чертёжника; "
                      "толщина линии 3 pt ≈ ±1,5 pt по вертикали — см. kW_per_pt. Значения — ПРЕДВАРИТЕЛЬНО (оцифровка)",
               software="poppler pdftocairo, python3", date=__import__("datetime").date.today().isoformat(),
               GS=digitize(sys.argv[1], "GS"), MD=digitize(sys.argv[2], "MD"))
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    dst = os.path.join(HERE, "data", "solar_available_power_taurus60.json")
    json.dump(out, open(dst, "w"), ensure_ascii=False, indent=1)
    for k in ("GS", "MD"):
        print(k, out[k]["frame_pt"], out[k]["curve_T_range_C"], [(p["T_amb_C"], p["P_kW"]) for p in out[k]["points"]])
