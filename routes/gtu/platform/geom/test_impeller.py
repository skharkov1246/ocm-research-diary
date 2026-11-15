"""Тесты impeller.py.

Запуск:  python -m pytest -q test_impeller.py   (или  python test_impeller.py)
Строит оба случая на уровне coarse через CLI во временную папку и проверяет записанные файлы
независимо от Gmsh: файл SU2 разбирается заново, объёмы и качество ячеек пересчитываются numpy.
"""
import math
import os
import shutil
import tempfile
import unittest

import numpy as np

import impeller as I

KEEP = os.environ.get("IMPELLER_TEST_KEEP")  # не удалять временные файлы


def read_su2(path):
    """минимальный разбор SU2 (ASCII): элементы, узлы, маркеры"""
    with open(path) as fh:
        L = fh.read().split("\n")
    out, i = {"markers": {}}, 0
    while i < len(L):
        s = L[i].strip()
        if s.startswith("NDIME="):
            out["ndime"] = int(s.split("=")[1])
        elif s.startswith("NELEM="):
            n = int(s.split("=")[1])
            out["elem"] = np.array([l.split() for l in L[i + 1:i + 1 + n]], dtype=np.int64)
            i += n
        elif s.startswith("NPOIN="):
            n = int(s.split("=")[1].split()[0])
            out["X"] = np.array([l.split()[:3] for l in L[i + 1:i + 1 + n]], dtype=float)
            i += n
        elif s.startswith("NMARK="):
            out["nmark"] = int(s.split("=")[1])
        elif s.startswith("MARKER_TAG="):
            tag = s.split("=", 1)[1].strip()
            n = int(L[i + 1].split("=")[1])
            out["markers"][tag] = np.array([l.split() for l in L[i + 2:i + 2 + n]], dtype=np.int64)
            i += n + 1
        i += 1
    return out


def tet_volumes(X, T):
    a, b, c, d = (X[T[:, k]] for k in range(4))
    return np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0


def tet_mean_ratio(X, T, V):
    """мера формы 12 (3V)^(2/3) / сумма квадратов рёбер: 1 у правильного тетраэдра, 0 у вырожденного"""
    P = [X[T[:, k]] for k in range(4)]
    l2 = sum(((P[i] - P[j]) ** 2).sum(1) for i in range(4) for j in range(i + 1, 4))
    return 12.0 * np.cbrt(3.0 * np.abs(V)) ** 2 / l2 * np.sign(V)


def face_keys(F, n):
    F = np.sort(F, axis=1)
    return (F[:, 0] * n + F[:, 1]) * n + F[:, 2]


class _CaseTest:
    case = None
    args = ()        # дополнительные ключи CLI
    over = {}        # те же переопределения для Impeller
    inputs = None    # исходные данные Impeller для case custom (у именованных случаев - case_inputs)

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix=f"impeller_{cls.case}_")
        cls.res = I.main(["--case", cls.case, "--level", "coarse", "--out", cls.tmp, *cls.args])
        cls.su2 = read_su2(cls.res["files"]["passage_su2"])
        cls.imp = I.Impeller(**(cls.inputs or I.case_inputs(cls.case)), **cls.over)

    @classmethod
    def tearDownClass(cls):
        if not KEEP:
            shutil.rmtree(cls.tmp, ignore_errors=True)

    # -------------------------------------------------------------- исходные данные
    def test_1d_inputs(self):
        d = self.res["d1"]
        lo, hi = self.D2_range
        self.assertTrue(lo <= d["D2"] <= hi, d["D2"])
        self.assertAlmostEqual(d["PR"], self.PR, places=6)
        if self.rpm_range:
            self.assertTrue(self.rpm_range[0] <= d["rpm"] <= self.rpm_range[1], d["rpm"])
        self.assertLessEqual(d["U2"], 620.0 + 1.0)
        self.assertAlmostEqual(d["gap"], max(0.3e-3, 0.0007 * d["D2"]), places=12)
        self.assertTrue(0.0 < d["beta2b"] < 60.0)
        if self.res["params"]["beta2b"] is not None:  # загиб задан: подбора под psi нет
            self.assertEqual(d["beta2b"], self.res["params"]["beta2b"])
            return
        self.assertAlmostEqual(d["psi_euler"], d["psi"], places=6)
        # независимо от кода: Визнер, sigma = 1 - Cслип/U2 = 1 - sqrt(cos b2b)/Z^0,7 (Z - лопатки на выходе),
        # при входе без закрутки Ctheta2/U2 = sigma - phi2 tg b2b = psi
        b = math.radians(d["beta2b"])
        Zex = self.res["params"]["Z"] * (2 if self.res["params"]["splitters"] else 1)
        sigma = 1.0 - math.sqrt(math.cos(b)) / Zex ** 0.7
        self.assertAlmostEqual(sigma - d["phi2"] * math.tan(b), d["psi"], places=6)

    def test_hub_volume_pappus(self):
        """объём втулки (OCC) против теоремы Паппа-Гульдина по меридиональному контуру, без OCC"""
        if "solid" not in self.res:
            self.skipTest("твёрдое тело не строилось")
        d, s, imp = self.res["d1"], self.res["solid"], self.imp
        hub = I._bez(imp.Ph, np.linspace(0.0, 1.0, 200001))
        poly = np.vstack([[s["z_nose"], s["r_bore"]], [s["z_nose"], d["r1h"]], hub,
                          [d["Lz"] + s["t_back"], d["r2"]], [d["Lz"] + s["t_back"], s["r_bore"]],
                          [s["z_nose"], s["r_bore"]]])
        z, r = poly[:, 0], poly[:, 1]
        V = abs(math.pi * np.sum(np.diff(z) * (r[1:] ** 2 + r[1:] * r[:-1] + r[:-1] ** 2) / 3.0))
        self.assertLess(abs(V / s["volume_parts_m3"]["hub"] - 1.0), 1e-6)

    # -------------------------------------------------------------- STEP
    def test_step_solid(self):
        import gmsh
        if "solid" not in self.res:
            self.skipTest("твёрдое тело не строилось")
        path = self.res["files"]["impeller_step"]
        self.assertTrue(os.path.getsize(path) > 0)
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.occ.importShapes(path)
            gmsh.model.occ.synchronize()
            vols = gmsh.model.getEntities(3)
            self.assertEqual(len(vols), 1, "колесо должно быть одним телом")
            # объём прочитанного STEP (OCC) против объёма по частям при построении:
            # втулка + Z (лопатки - их часть внутри втулки). У двустороннего колеса getMass OCC сплавленного тела
            # точнее 1e-4 не бывает: на ступени 1 KV-3 сплавленная сторона давала +0,49 %, её зеркальная копия
            # (та же геометрия) +0,36 %, обе стороны вместе -0,01 % - ошибка интегрирования по граням втулки с
            # вырезами (README); там допуск как у объёма сетки ниже - 2e-3, он и есть независимая проверка
            V = gmsh.model.occ.getMass(3, vols[0][1])
            tol_occ = 2e-3 if self.res["d1"].get("double_sided") else 1e-4
            self.assertLess(abs(V / self.res["solid"]["volume_m3"] - 1.0), tol_occ)
            gmsh.option.setNumber("Mesh.MeshSizeMax", 0.01)
            gmsh.option.setNumber("Mesh.MeshSizeMin", 0.0004)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 16)
            gmsh.option.setNumber("Mesh.Algorithm3D", 10)
            gmsh.model.mesh.generate(3)
            tags, X, _ = gmsh.model.mesh.getNodes()
            X = X.reshape(-1, 3)
            r = np.hypot(X[:, 0], X[:, 1])
            r2 = self.res["d1"]["r2"]
            self.assertLess(r.max(), r2 * (1 + 1e-6))
            self.assertGreater(r.max(), r2 * (1 - 1e-4))
            # третий способ: объём тетраэдрической сетки прочитанного STEP (без интегрирования OCC)
            idx = np.zeros(int(tags.max()) + 1, dtype=np.int64)
            idx[tags.astype(np.int64)] = np.arange(len(tags))
            _, _, en = gmsh.model.mesh.getElements(3)
            Vt = np.abs(tet_volumes(X, idx[np.concatenate(en).astype(np.int64)].reshape(-1, 4))).sum()
            self.assertLess(abs(Vt / self.res["solid"]["volume_m3"] - 1.0), 2e-3)
        finally:
            gmsh.finalize()

    # -------------------------------------------------------------- сетка SU2
    def test_markers_present(self):
        m = self.su2
        self.assertEqual(m["ndime"], 3)
        self.assertEqual(m["nmark"], len(I.MARKERS))
        self.assertEqual(set(m["markers"]), set(I.MARKERS))
        for k, e in m["markers"].items():
            self.assertGreater(len(e), 0, k)
            self.assertTrue(np.all(e[:, 0] == 5), f"{k}: ожидались треугольники")

    def test_cells_nondegenerate(self):
        E, X = self.su2["elem"], self.su2["X"]
        self.assertTrue(np.all(E[:, 0] == 10), "ожидались только тетраэдры")
        T = E[:, 1:5]
        V = tet_volumes(X, T)
        q = tet_mean_ratio(X, T, V)
        print(f"\n[{self.case}] ячеек {len(T)}, узлов {len(X)}; объём ячейки min {V.min():.3e} м3, "
              f"мера формы min {q.min():.3f}; minSICN (Gmsh) {self.res['mesh']['quality']['minSICN']:.3f}")
        self.assertTrue(np.all(V > 0), "есть ячейки с неположительным объёмом")
        self.assertGreater(q.min(), 0.0)
        self.assertGreater(self.res["mesh"]["quality"]["minSICN"], 0.0)
        # сумма объёмов ячеек = объём области OCC (сетка заполняет область; линейные тетраэдры
        # на кривых стенках дают отклонение порядка 1e-3)
        Vf = self.res["mesh"]["volume_fluid_occ_m3"]
        print(f"[{self.case}] объём сетки / объём области OCC - 1 = {V.sum() / Vf - 1.0:.2e}")
        self.assertLess(abs(V.sum() / Vf - 1.0), 5e-3)
        # канал до вычитания лопаток = 1/Z тела вращения (проверка выбора объёма и сектора)
        self.assertLess(abs(self.res["mesh"]["passage_volume_vs_body_over_Z"]), 1e-3)

    def test_boundary_closed_and_tagged(self):
        E, X = self.su2["elem"], self.su2["X"]
        n = len(X)
        T = E[:, 1:5]
        F = np.vstack([T[:, [0, 1, 2]], T[:, [0, 1, 3]], T[:, [0, 2, 3]], T[:, [1, 2, 3]]])
        k, cnt = np.unique(face_keys(F, n), return_counts=True)
        self.assertTrue(np.all(cnt <= 2), "грань принадлежит более чем двум ячейкам")
        bnd = set(k[cnt == 1].tolist())
        tagged = []
        for e in self.su2["markers"].values():
            tagged.extend(face_keys(e[:, 1:4], n).tolist())
        self.assertEqual(len(tagged), len(set(tagged)), "треугольник в двух маркерах")
        self.assertEqual(set(tagged), bnd, "граница области и маркеры не совпадают")

    def test_periodic_rotation(self):
        Z = self.res["params"]["Z"]
        X = self.su2["X"]
        p1, p2 = self.su2["markers"]["periodic_1"][:, 1:4], self.su2["markers"]["periodic_2"][:, 1:4]
        self.assertEqual(len(p1), len(p2))
        n1, n2 = np.unique(p1), np.unique(p2)
        self.assertEqual(len(n1), len(n2))
        a = 2 * math.pi / Z
        R = np.array([[math.cos(a), -math.sin(a), 0], [math.sin(a), math.cos(a), 0], [0, 0, 1]])
        Y = X[n1] @ R.T
        from scipy.spatial import cKDTree
        dist, j = cKDTree(X[n2]).query(Y)
        D2 = self.res["d1"]["D2"]
        snap = self.res["mesh"]["periodic_snap"]
        print(f"\n[{self.case}] periodic: {len(n1)} узлов, наибольшее расхождение после поворота на 360/{Z}: "
              f"{dist.max():.2e} м (поправка узлов при записи: {snap['n_nodes']} узлов, до {snap['max_shift_m']:.2e} м)")
        # поправка - порядка допуска OCC, а не перестановка граней
        self.assertLess(snap["max_shift_m"], 1e-6 * D2)
        self.assertLess(dist.max(), 1e-10 * D2)
        mp = np.full(len(X), -1, dtype=np.int64)
        mp[n1] = n2[j]
        self.assertEqual(len(set(mp[n1].tolist())), len(n1), "отображение узлов не взаимно однозначно")
        n = len(X)
        self.assertEqual(set(face_keys(mp[p1], n).tolist()), set(face_keys(p2, n).tolist()),
                         "треугольники periodic_1 после поворота не совпали с periodic_2")

    def test_marker_geometry(self):
        X, imp = self.su2["X"], self.imp
        d = self.res["d1"]
        tol = 1e-5 * d["D2"]

        def nodes(k):
            P = X[np.unique(self.su2["markers"][k][:, 1:4])]
            return P, np.hypot(P[:, 0], P[:, 1])

        P, r = nodes("inlet")
        self.assertLess(np.abs(P[:, 2] - imp.z_in).max(), tol)
        P, r = nodes("outlet")
        self.assertLess(np.abs(r - imp.r_out).max(), tol)
        P, r = nodes("hub")
        self.assertLess(I._seg_dist(np.stack([P[:, 2], r], 1), imp.hub_poly).max(), tol)
        P, r = nodes("shroud")
        self.assertLess(I._seg_dist(np.stack([P[:, 2], r], 1), imp.sh_poly).max(), tol)
        P, r = nodes("blade")
        self.assertLess(r.max(), d["r2"] + tol)
        # ни одна грань входного участка не попала в blade (классификация по остатку): лопатка начинается
        # у z = 0, наклон торца входной кромки - не больше половины наибольшей толщины
        self.assertGreater(P[:, 2].min(), -0.5 * self.res["params"]["t_hub"] * d["D2"])
        self.assertGreater(r.min(), d["r1h"] - tol)
        ds = I._seg_dist(np.stack([P[:, 2], r], 1), imp.sh_poly)
        print(f"\n[{self.case}] радиальный зазор по узлам лопатки: min {ds.min() * 1e3:.4f} мм "
              f"(задан {d['gap'] * 1e3:.4f} мм)")
        self.assertGreater(ds.min(), d["gap"] - tol)
        self.assertLess(ds.min(), d["gap"] + tol)


class TestCore130(_CaseTest, unittest.TestCase):
    """схема 130: колесо 0,435 м, 22 200 об/мин, PR 4,25"""
    case = "core130"
    D2_range = (0.434, 0.436)
    rpm_range = (22150.0, 22250.0)
    PR = 4.25


class TestTaurus60(_CaseTest, unittest.TestCase):
    """пакет Taurus 60: одно колесо, PR 6,75, D2 около 0,6 м"""
    case = "taurus60"
    D2_range = (0.57, 0.63)
    rpm_range = None
    PR = 6.75


class TestFewBlades(_CaseTest, unittest.TestCase):
    """граничный случай: Z = 3 (шаг 120 град). До исправления развёртки угла сектора область строилась
    не между периодическими поверхностями и падала с «нет граней для маркеров ['blade']»"""
    case = "core130"
    args = ("--z", "3", "--no-step")
    over = dict(Z=3)
    D2_range = (0.434, 0.436)
    rpm_range = (22150.0, 22250.0)
    PR = 4.25


# SYNTHETIC: придуманное двустороннее колесо с горячим входом (как у второй ступени) и заданными размерами -
# проверка ключей --double-sided, --t01/--p01, --d1s/--nu/--b2/--lz/--beta1h/--beta1s/--beta2b, --r-bore 0.
# Значения не из расчёта трека, подобраны только для дешёвой сетки (b2/D2 около 0,045).
SYNTH_DBL = dict(m=12.0, PR=3.0, rpm=30000.0, D2=0.36, T01=400.0, p01=1.5e5)        # SYNTHETIC
SYNTH_DBL_OVER = dict(Z=8, double_sided=True, D1s=0.20, nu=0.36, b2=0.016, Lz=0.095,  # SYNTHETIC
                      beta1h=30.0, beta1s=58.0, beta2b=45.0, r_bore=0.0)


class TestDoubleSided(_CaseTest, unittest.TestCase):
    """двусторонний вход: канал и сетка - одной стороны на половину расхода, STEP - обе стороны зеркально;
    полные параметры на входе и размеры заданы извне (SYNTHETIC)"""
    case = "custom"
    s, o = SYNTH_DBL, SYNTH_DBL_OVER
    args = ("--m", str(s["m"]), "--pr", str(s["PR"]), "--rpm", str(s["rpm"]), "--d2", str(s["D2"]),
            "--t01", str(s["T01"]), "--p01", str(s["p01"]), "--double-sided", "--z", str(o["Z"]),
            "--d1s", str(o["D1s"]), "--nu", str(o["nu"]), "--b2", str(o["b2"]), "--lz", str(o["Lz"]),
            "--beta1h", str(o["beta1h"]), "--beta1s", str(o["beta1s"]), "--beta2b", str(o["beta2b"]),
            "--r-bore", str(o["r_bore"]))
    inputs = dict(SYNTH_DBL)
    over = dict(SYNTH_DBL_OVER)
    D2_range = (0.36, 0.36)
    rpm_range = (30000.0, 30000.0)
    PR = 3.0
    del s, o

    def test_given_inputs(self):
        """заданные величины дошли до 1D-данных без изменений; сторона - на половину расхода"""
        d, s, o = self.res["d1"], SYNTH_DBL, SYNTH_DBL_OVER
        self.assertTrue(d["double_sided"])
        self.assertEqual(d["m"], s["m"])
        self.assertEqual(d["m_side"], 0.5 * s["m"])
        self.assertAlmostEqual(self.res["su2"]["mass_flow_passage_kg_s"], 0.5 * s["m"] / o["Z"], places=12)
        self.assertEqual((d["T01"], d["p01"]), (s["T01"], s["p01"]))
        self.assertEqual(self.res["su2"]["inlet_total"], dict(p01_Pa=s["p01"], T01_K=s["T01"]))
        self.assertAlmostEqual(2 * d["r1s"], o["D1s"], places=12)
        self.assertAlmostEqual(d["r1h"], o["nu"] * o["D1s"] / 2, places=12)
        for k in ("b2", "Lz", "beta1h", "beta1s", "beta2b"):
            self.assertAlmostEqual(d[k], o[k], places=12, msg=k)
        self.assertAlmostEqual(d["Ls"], o["Lz"] - o["b2"], places=12)

    def test_mesh_sizes_from_file(self):
        """размеры прочитаны из самого файла SU2, не из JSON: вход r1h..r1s, высота выхода b2, D2 по лопатке,
        шаг периодичности 360/Z"""
        X, o, d = self.su2["X"], SYNTH_DBL_OVER, self.res["d1"]
        tol = 1e-5 * d["D2"]

        def nodes(k):
            return X[np.unique(self.su2["markers"][k][:, 1:4])]
        P = nodes("inlet")
        r = np.hypot(P[:, 0], P[:, 1])
        self.assertLess(abs(r.min() - o["nu"] * o["D1s"] / 2), tol)
        self.assertLess(abs(r.max() - o["D1s"] / 2), tol)
        P = nodes("outlet")
        self.assertLess(abs((P[:, 2].max() - P[:, 2].min()) - o["b2"]), tol)
        P = nodes("blade")
        self.assertLess(abs(2 * np.hypot(P[:, 0], P[:, 1]).max() - SYNTH_DBL["D2"]), 2 * tol)
        # угол между соответственными узлами periodic_1 и periodic_2 - по их среднему углу
        a1 = np.unwrap(np.arctan2(nodes("periodic_1")[:, 1], nodes("periodic_1")[:, 0])).mean()
        a2 = np.unwrap(np.arctan2(nodes("periodic_2")[:, 1], nodes("periodic_2")[:, 0])).mean()
        self.assertAlmostEqual(math.degrees(a2 - a1) % 360.0, 360.0 / o["Z"], places=6)

    def test_mirror_solid(self):
        """STEP двустороннего колеса: две стороны, зеркальные относительно средней плоскости общего диска"""
        import gmsh
        s = self.res["solid"]
        self.assertEqual(s["sides"], 2)
        self.assertAlmostEqual(s["z_mid_plane"], SYNTH_DBL_OVER["Lz"] + s["t_back"], places=12)
        self.assertEqual(s["r_bore"], 0.0)
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.occ.importShapes(self.res["files"]["impeller_step"])
            gmsh.model.occ.synchronize()
            gmsh.option.setNumber("Mesh.MeshSizeMax", 0.01)
            gmsh.model.mesh.generate(1)  # узлы на рёбрах: z-размах тела без грубой оценки рамки OCC
            X = gmsh.model.mesh.getNodes()[1].reshape(-1, 3)
        finally:
            gmsh.finalize()
        zc = s["z_mid_plane"]
        self.assertLess(abs((X[:, 2].min() + X[:, 2].max()) - 2 * zc), 1e-6)
        self.assertLess(abs(X[:, 2].min() - s["z_nose"]), 1e-6)
        # отражение узлов рёбер относительно средней плоскости ложится на узлы рёбер (с точностью до шага)
        from scipy.spatial import cKDTree
        Y = X.copy(); Y[:, 2] = 2 * zc - Y[:, 2]
        dist, _ = cKDTree(X).query(Y)
        self.assertLess(np.median(dist), 1e-9)


class TestInputs(unittest.TestCase):
    """явный отказ на недопустимых исходных данных (без Gmsh)"""

    def test_rejects(self):
        for kw in (dict(m=-1.0, PR=3.0), dict(m=5.0, PR=0.9), dict(m=5.0, PR=3.0, Z=2),
                   dict(m=5.0, PR=3.0, gap=-1e-4), dict(m=5.0, PR=3.0, m_split=1.2),
                   dict(m=5.0, PR=3.0, beta2b=80.0), dict(m=5.0, PR=3.0, D2=-0.4),
                   dict(m=5.0, PR=3.0, T01=-10.0), dict(m=5.0, PR=3.0, p01=0.0),
                   dict(m=5.0, PR=3.0, D2=0.4, D1s=0.45), dict(m=5.0, PR=3.0, b2=0.02, Lz=0.015),
                   dict(m=5.0, PR=3.0, beta1s=95.0),
                   dict(m=50.0, PR=3.0, rpm=30000.0, D2=0.36, D1s=0.05)):  # вход запирается
            with self.assertRaises(ValueError, msg=str(kw)):
                I.Impeller(**kw)

    def test_double_sided_equals_half_flow(self):
        """двустороннее колесо с расходом 2m - сторона та же, что одностороннее с m (D2 и частота заданы)"""
        a = I.Impeller(m=10.0, PR=4.0, rpm=15000.0, D2=0.6, double_sided=True)   # SYNTHETIC
        b = I.Impeller(m=5.0, PR=4.0, rpm=15000.0, D2=0.6)                       # SYNTHETIC
        for k in ("b2", "r1s", "r1h", "beta1h", "beta1s", "beta2b", "Lz", "phi", "C1", "M1s_rel", "psi"):
            self.assertAlmostEqual(a.d1[k], b.d1[k], places=12, msg=k)
        self.assertEqual(a.d1["m_side"], 5.0)

    def test_defaults_unchanged(self):
        """без новых ключей вход стандартный (gtu_explore.T0, P0) и одна сторона на весь расход"""
        E = I.load_explore()
        imp = I.Impeller(m=5.0, PR=3.0)                                          # SYNTHETIC
        self.assertEqual((imp.d1["T01"], imp.d1["p01"]), (E.T0, E.P0))
        self.assertEqual(imp.d1["m_side"], 5.0)
        self.assertFalse(imp.d1["double_sided"])
        self.assertEqual(imp.d1["given"], [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
