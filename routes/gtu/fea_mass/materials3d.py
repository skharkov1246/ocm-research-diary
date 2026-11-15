"""Массовая 3D-прочность: свойства материалов роторных деталей - те же источники и правила, что у уровня 1.

  * Ti-6Al-4V и литой IN-713C - ../platform/fea/materials.json (статус каждого числа - там; сняты с выдержек поиска).
  * Никелевые сплавы колёс и лопаток турбины (K418, CM247LC, Haynes 282, IN718, K438, монокристалл) - как уровень 1
    (funnel/strength_l1.py): упругость, КТР и σв - таблицы IN-713C; длительная прочность - кривая Ларсона - Миллера
    IN-713C со сдвигом температуры dT_cap сплава (funnel/make_model.py, strength_l1.creep_material); плотность - make_model.MAT.
    Это оценка: настоящие кривые сплавов в репозитории не собраны.
  * σ0,2(T) никелевых сплавов: в materials.json только 20 °C (689 МПа); при нагреве - σ0,2(20)·σв(T)/σв(20) (оценка).
  * Алюминиевый сплав колеса и нитрид кремния - в materials.json их нет; ниже - типовые справочные значения по памяти,
    НЕ ПРОВЕРЕНО (оценка), у каждого числа - пометка.
  * Малоцикловая усталость - констант Мэнсона - Коффина в репозитории нет (навык calculix-rotor, шаг 7). Оценка - методом
    модифицированных универсальных наклонов (Muralidharan, Manson, 1988; формула общеизвестна, первоисточник в сессии не
    открыт): Δε = 0,623·(σв/E)^0,832·N^-0,09 + 0,0196·εf^0,155·(σв/E)^-0,53·N^-0,56, εf = ln(1/(1-ψ)); относительное
    сужение ψ - принятое (не проверено). Поэтому число циклов - оценка низкой достоверности, в вердикте - «условно».
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.dirname(HERE)
FUNNEL = os.path.join(GTU, "funnel")
FEA = os.path.join(GTU, "platform", "fea")

_MJ = None
_ST = None


def _mats_json():
    global _MJ
    if _MJ is None:
        _MJ = json.load(open(os.path.join(FEA, "materials.json"), encoding="utf-8"))
    return _MJ


def strength_l1():
    """funnel/strength_l1.py (Ларсон - Миллер IN-713C, сдвиги сплавов) - импорт без правки кода воронки"""
    global _ST
    if _ST is None:
        if FUNNEL not in sys.path:
            sys.path.insert(0, FUNNEL)
        import strength_l1 as ST          # noqa: E402
        _ST = ST
    return _ST


def interp(tab, x):
    xs, ys = [t[0] for t in tab], [t[1] for t in tab]
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            return ys[i] + (ys[i + 1] - ys[i]) * (x - xs[i]) / (xs[i + 1] - xs[i])
    return ys[-1]


# относительное сужение ψ для оценки МЦУ - НЕ ПРОВЕРЕНО (типовые значения по памяти справочника, оценка)
PSI = dict(ti64=0.25, in713c=0.08, al=0.05, si3n4=None)
PSI_BASIS = ("ψ: Ti-6Al-4V отожжённый 0,25; литой IN-713C 0,08; литой алюминиевый сплав 0,05 - типовые справочные значения по "
             "памяти, первоисточник не открыт (оценка, не проверено)")

# алюминиевый сплав колеса (литьё ZL114A / печать AlSi10Mg; деформируемый 2618 - в impeller.py плотность 2760): оценка
AL = dict(name="алюминиевый сплав колеса (оценка свойств)", density=2800.0, nu=0.33,
          E_GPa=[[20, 72.0], [200, 66.0]], alpha_mean_1e6=[[20, 21.5], [200, 23.0]],
          uts_MPa=[[20, 400.0], [150, 360.0], [200, 300.0]], ys_MPa=[[20, 340.0], [150, 310.0], [200, 250.0]],
          status="НЕ ПРОВЕРЕНО: типовые значения деформируемого 2618-T61 по памяти справочника (оценка); плотность 2800 - "
                 "make_model.MAT['Al'] (gtu_explore.finalize). Литой сплав слабее")
# нитрид кремния: упругость и КТР - оценка; прочность - Вейбулл уровня 1 (strength_l1: SN282A 594 МПа при 1200 °C)
SI3N4 = dict(name="нитрид кремния (оценка свойств)", density=3200.0, nu=0.27, E_GPa=[[20, 300.0], [1200, 280.0]],
             alpha_mean_1e6=[[20, 3.2], [1200, 3.6]],
             status="НЕ ПРОВЕРЕНО: E и КТР - типовые значения спечённого Si3N4 по памяти (оценка); плотность 3200 - make_model.MAT")


class Material:
    """свойства для ccx и критериев. key - ключ данных (ti64 | in713c | al | si3n4); alloy - сплав записи уровня 1"""

    def __init__(self, key, alloy=None, rho=None):
        self.key, self.alloy = key, alloy or key
        self.T_ref = 20.0
        self.flags = []
        if key in ("ti64", "in713c"):
            m = _mats_json()[key]
            self.src = f"routes/gtu/platform/fea/materials.json[{key}]"
        elif key == "al":
            m = AL
            self.src = "fea_mass/materials3d.py AL (оценка, не проверено)"
        elif key == "si3n4":
            m = dict(SI3N4, uts_MPa=None, ys_MPa=None)
            self.src = "fea_mass/materials3d.py SI3N4 (оценка, не проверено)"
        else:
            raise ValueError(f"нет данных материала {key!r}")
        self.name = m["name"]
        self.rho = float(rho if rho is not None else m["density"])
        self.nu = float(m["nu"])
        self.E_tab = [(float(T), float(E)) for T, E in m["E_GPa"]]
        self.alpha_tab = [(float(T), float(a)) for T, a in m["alpha_mean_1e6"]]
        self.uts_tab = m.get("uts_MPa")
        ys = m.get("ys_MPa")
        self.ys_tab = ys if ys and len(ys) > 1 else None
        self.ys20 = ys[0][1] if ys else None
        self.psi = PSI.get(key)
        # длительная прочность: только никелевые сплавы (кривая IN-713C со сдвигом)
        self.creep = key == "in713c"
        self.dT_cap, self.creep_basis = (0.0, None)
        if self.creep:
            self.dT_cap, self.creep_basis = strength_l1().creep_material(self.alloy)
        if alloy and key == "in713c" and alloy not in ("in713c", "K418", "Ni"):
            self.flags.append(f"{alloy}: упругость, КТР и σв - таблицы IN-713C (оценка, как уровень 1)")

    def E(self, T):
        return interp(self.E_tab, T) * 1e9

    def uts(self, T):
        """σв(T), МПа (за краем таблицы - крайнее значение, как disk.py)"""
        return interp(self.uts_tab, T) if self.uts_tab else None

    def ys(self, T):
        """σ0,2(T), МПа: таблица; у никелевых - σ0,2(20)·σв(T)/σв(20) (оценка)"""
        if self.ys_tab:
            return interp(self.ys_tab, T)
        if self.ys20 and self.uts_tab:
            return self.ys20 * self.uts(T) / self.uts_tab[0][1]
        return None

    def ys_basis(self):
        if self.ys_tab:
            return f"таблица {self.src}"
        if self.ys20:
            return "σ0,2(20 °C)·σв(T)/σв(20 °C) (оценка: при нагреве данных σ0,2 нет)"
        return "нет данных"

    def rupture_h(self, sigma_MPa, T_C):
        """ресурс по длительной прочности, ч (Ларсон - Миллер, strength_l1.rupture_hours); (часы | None, экстраполяция)"""
        if not self.creep:
            return None, None
        if sigma_MPa <= 0:
            return math.inf, False
        t, ext = strength_l1().rupture_hours(sigma_MPa, T_C, self.dT_cap)
        if t is None:
            return 0.0, True            # напряжение выше кривой - ресурса нет
        return t, ext

    def lcf_cycles(self, d_eps, T_C):
        """оценка числа циклов до зарождения трещины по размаху деформации (модифицированные универсальные наклоны)"""
        if self.psi is None or not self.uts_tab or d_eps <= 0:
            return None
        su_E = self.uts(T_C) * 1e6 / self.E(T_C)
        ef = math.log(1.0 / (1.0 - self.psi))

        def de(N):
            return 0.623 * su_E ** 0.832 * N ** -0.09 + 0.0196 * ef ** 0.155 * su_E ** -0.53 * N ** -0.56
        lo, hi = 1.0, 1e12
        if de(hi) >= d_eps:
            return hi
        if de(lo) <= d_eps:
            return lo
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            if de(mid) > d_eps:
                lo = mid
            else:
                hi = mid
            if hi / lo < 1.0001:
                break
        return math.sqrt(lo * hi)

    def describe(self):
        return dict(key=self.key, alloy=self.alloy, name=self.name, rho=self.rho, nu=self.nu, source=self.src,
                    ys_basis=self.ys_basis(), creep=("Ларсон - Миллер IN-713C, " + self.creep_basis) if self.creep else "нет данных",
                    lcf=(f"универсальные наклоны, ψ = {self.psi} (оценка, не проверено)" if self.psi else "нет данных"),
                    flags=list(self.flags))


def for_part(alloy):
    """материал детали по сплаву записи уровня 1: Ti64 | Al | Si3N4 | никелевые (K418, CM247LC, Haynes282, IN718, K438, SC, Ni)"""
    if alloy in ("Ti64", "ti64"):
        return Material("ti64", "Ti64")
    if alloy in ("Al", "al"):
        return Material("al", "Al")
    if alloy in ("Si3N4",):
        return Material("si3n4", "Si3N4")
    rho = None
    try:
        MK = strength_l1().MK
        rho = MK.MAT.get(alloy, {}).get("rho")
    except Exception:            # noqa: BLE001 - плотность из materials.json
        rho = None
    return Material("in713c", alloy, rho=rho)
