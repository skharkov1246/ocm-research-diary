"""Проверки frame3_geom.py и frame3_geom_results.json (без построения геометрии, секунды).

Запуск:  python3 -m pytest -q test_frame3_geom.py
Сетки и STEP строит сам frame3_geom.py (минуты); здесь - перенос чисел точки в impeller и итоговый JSON.
"""
import hashlib
import json
import os
import unittest

import frame3_geom as F

HERE = os.path.dirname(os.path.abspath(__file__))


class TestKwargs(unittest.TestCase):
    """числа точки доходят до одномерных данных impeller без изменений (без Gmsh)"""

    @classmethod
    def setUpClass(cls):
        cls.point = json.load(open(F.POINT, encoding="utf-8"))

    def test_stage_mapping(self):
        for key in F.STAGES:
            s = self.point[key]
            imp = F.I.Impeller(**F.stage_kwargs(s))
            d = imp.d1
            self.assertEqual(d["double_sided"], s["inlet"] == "двусторонний", key)
            self.assertAlmostEqual(d["m_side"], s["m_per_side_kg_s"], places=9, msg=key)
            # окружная скорость: impeller считает по D2 и частоте, точка - своим расчётом
            self.assertAlmostEqual(d["U2"], s["U2_m_s"], delta=1e-3, msg=key)
            self.assertAlmostEqual(2 * d["r1h"], s["D1_hub_m"], places=9, msg=key)
            self.assertAlmostEqual(2 * d["r1s"], s["D1_tip_m"], places=9, msg=key)
            for k, f in (("b2", "b2_m"), ("Lz", "axial_length_m"), ("gap", "tip_clearance_m"),
                         ("beta2b", "beta2_blade_from_radial_deg"), ("beta1h", "beta1_hub_deg"), ("beta1s", "beta1_tip_deg")):
                self.assertAlmostEqual(d[k], s[f], places=9, msg=f"{key}.{k}")
            self.assertAlmostEqual(d["T01"], s["T01_C"] + 273.15, places=9)
            self.assertAlmostEqual(d["p01"], s["p01_kPa"] * 1e3, places=6)
            self.assertEqual(imp.p["Z"], s["blades_full"])
            # свой треугольник скоростей impeller (без загромождения входа) - в пределах 2 град от точки
            self.assertLess(abs(d["beta1s_flow"] - s["beta1_tip_deg"]), 2.0)
            self.assertEqual(imp.warnings, [], key)


@unittest.skipUnless(os.path.exists(F.RESULTS), "frame3_geom_results.json ещё не записан")
class TestResults(unittest.TestCase):
    """итоговый JSON: обязательные поля, проверки пройдены, вход не менялся после прогона"""

    @classmethod
    def setUpClass(cls):
        cls.r = json.load(open(F.RESULTS, encoding="utf-8"))

    def test_fields(self):
        for k in ("method", "software", "date", "inputs", "stages", "checks", "not_done"):
            self.assertIn(k, self.r)
        self.assertEqual(self.r["date"], "2026-11-14")
        for inp in self.r["inputs"]:
            self.assertIn("file", inp)
            self.assertIn("field", inp)

    def test_checks_ok(self):
        for key in F.STAGES:
            c = self.r["checks"][key]
            self.assertIn("coarse", c, key)
            for level, v in c.items():
                self.assertTrue(v["ok"], f"{key} {level}: {v}")
            st = self.r["stages"][key]
            for level, L in st["levels"].items():
                self.assertGreater(L["minSICN"], 0.0)
                self.assertEqual(set(L["check_su2"]["markers"]), set(F.I.MARKERS))
                self.assertAlmostEqual(L["check_su2"]["periodic"]["angle_deg"], 360.0 / st["params_from_point"]["Z"])

    def test_input_unchanged(self):
        """точка не менялась после прогона (иначе геометрию надо перестроить)"""
        inp = next(i for i in self.r["inputs"] if i["file"].endswith("frame3_point_results.json"))
        self.assertEqual(inp["sha256"], hashlib.sha256(open(F.POINT, "rb").read()).hexdigest(),
                         "frame3_point_results.json изменился после прогона frame3_geom.py - перестроить")

    def test_big_files_where_present(self):
        """большие файлы (вне репозитория или в out/ из .gitignore): если они на месте - SHA-256 совпадает"""
        n = 0
        for key in F.STAGES:
            st = self.r["stages"][key]
            infos = [st["solid_file"]] + [f for L in st["levels"].values() for f in L["files"].values()]
            infos += [f["stored_as"] for f in infos if "stored_as" in f]
            for f in infos:
                if os.path.exists(f["path"]) and f["bytes"] < 2e8:
                    self.assertEqual(F.sha256(f["path"]), f["sha256"], f["path"])
                    n += 1
        if n == 0:
            self.skipTest("большие файлы не найдены (scratchpad очищен)")


if __name__ == "__main__":
    unittest.main(verbosity=2)
