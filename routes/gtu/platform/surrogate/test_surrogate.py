"""Тесты surrogate.py на заменителе дорогого расчёта.

Заменитель — центробежная ступень gtu_explore.comp_stage (титановое колесо), обращённая по
окружной скорости: входы U2 [м/с], m [кг/с], T1 [K] -> КПД eta, степень сжатия PR, диаметр
колеса D2 [м]. Степень сжатия находится brentq так, чтобы comp_stage дала заданную U2.
К eta добавлен шум СКО 0,003 (0,3 п. п.), к PR — 0,5 % (имитация шума сходимости и сетки),
D2 — без шума (геометрия).

Критерии задачи:
  * 200 точек плана -> ошибка на отложенной выборке не более 2 % по КПД;
  * оптимум, найденный по модели, проверяется прямым расчётом: отличие не более 1 п. п.

Запуск:  python -m pytest -q -s test_surrogate.py   (или  python test_surrogate.py)
SURROGATE_RECORD=1 — записать итоги в test_surrogate_results.json рядом.
"""
import json
import os
import sys
import tempfile
import time
import unittest

import numpy as np
from scipy.optimize import brentq, differential_evolution

HERE = os.path.dirname(os.path.abspath(__file__))
GTU = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, GTU)
_cwd = os.getcwd()
import gtu_explore as G  # noqa: E402  (модуль при импорте переходит в свою папку)
os.chdir(_cwd)
import surrogate as S  # noqa: E402

JOBS = int(os.environ.get("SURROGATE_JOBS", "4"))
RECORD = os.environ.get("SURROGATE_RECORD") == "1"
OUT = ["eta", "PR", "D2"]
SPACE = S.Space({"U2": [350.0, 650.0], "m": [1.5, 25.0], "T1": [243.15, 318.15]})
NOISE = dict(eta=0.003, PR_rel=0.005)
T_ISO = 288.15
MLP = dict(members=8, hidden=(32, 32), alpha="auto")
RESULTS = {}

# случаи оптимизации: максимум КПД при ограничениях, T1 = ISO
CASES = {
    "A_core130": dict(cons=["PR>=4.25"], note="степень сжатия ядра схемы 130"),
    "B_core130_D040": dict(cons=["PR>=4.25", "D2<=0.40"], note="то же, колесо не больше 0,40 м"),
    "C_taurus60": dict(cons=["PR>=6.5", "D2<=0.60", "U2<=620"], note="одно колесо пакета Taurus 60, титан 620 м/с"),
}


def stage(U2, m, T1):
    """прямой расчёт: (eta, PR, D2) ступени с окружной скоростью U2"""
    f = lambda pr: G.comp_stage(T1, G.P0, pr, m, "ti")["U2"] - U2  # noqa: E731
    pr = brentq(f, 1.05, 15.0, xtol=1e-12)
    r = G.comp_stage(T1, G.P0, pr, m, "ti")
    return r["eta"], pr, r["D2"]


def truth(X):
    return np.array([stage(*x) for x in np.atleast_2d(X)])


def noisy(Yt, seed):
    rng = np.random.default_rng(seed)
    Y = Yt.copy()
    Y[:, 0] += rng.normal(0.0, NOISE["eta"], len(Y))
    Y[:, 1] *= 1.0 + rng.normal(0.0, NOISE["PR_rel"], len(Y))
    return Y


def true_optimum(cons, T1=T_ISO, lo=None, hi=None):
    """истинный оптимум прямым расчётом: differential_evolution по (U2, m) без модели"""
    cons = [S.parse_constraint(c) for c in cons]
    lo = SPACE.lo if lo is None else lo
    hi = SPACE.hi if hi is None else hi

    def f(z):
        eta, pr, d2 = stage(z[0], z[1], T1)
        v = dict(eta=eta, PR=pr, D2=d2, U2=z[0], m=z[1])
        pen = sum(max(0.0, (c - v[n]) if op == ">=" else (v[n] - c)) for n, op, c in cons)
        return -eta + 10.0 * pen

    r = differential_evolution(f, [(lo[0], hi[0]), (lo[1], hi[1])], tol=1e-10, maxiter=300, popsize=25,
                               polish=True, **S._rng_kw(differential_evolution, 3))
    eta, pr, d2 = stage(r.x[0], r.x[1], T1)
    return dict(U2=float(r.x[0]), m=float(r.x[1]), eta=float(eta), PR=float(pr), D2=float(d2))


def check_cons(v, cons, rel=0.01):
    """ограничения выполнены прямым расчётом с запасом rel (ошибка модели по PR и D2 — доли процента)"""
    for n, op, c in (S.parse_constraint(c) for c in cons):
        if op == ">=" and v[n] < c * (1 - rel):
            return False
        if op == "<=" and v[n] > c * (1 + rel):
            return False
    return True


class SurrogateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        t0 = time.time()
        cls.X, info = S.design(SPACE, 200, "lhs", seed=1)
        cls.Yt = truth(cls.X)
        cls.Y = noisy(cls.Yt, seed=7)
        cls.Xh, _ = S.design(SPACE, 256, "sobol", seed=2)      # отложенная выборка — другой план
        cls.Yht = truth(cls.Xh)
        cls.Yh = noisy(cls.Yht, seed=8)
        t1 = time.time()
        cls.sur = S.Surrogate(SPACE, OUT, kinds=("mlp", "gp"), mlp=MLP, log=("D2",), seed=0, n_jobs=JOBS)
        cls.sur.fit(cls.X, cls.Y)
        cls.t_fit = time.time() - t1
        RESULTS["setup"] = dict(plan=info, holdout_n=len(cls.Xh), noise=NOISE, space=SPACE.bounds(),
                                truth_range={o: [float(cls.Yt[:, j].min()), float(cls.Yt[:, j].max())]
                                             for j, o in enumerate(OUT)},
                                fit_seconds=round(cls.t_fit, 1), data_seconds=round(t1 - t0, 1),
                                model=cls.sur.describe())
        print(f"\nобучение MLP+GP на 200 точках: {cls.t_fit:.1f} с; alpha MLP: {cls.sur.models['mlp'].alpha_}; "
              f"шум по GP: { {k: round(v, 5) for k, v in cls.sur.noise_sd().items()} }")

    @classmethod
    def tearDownClass(cls):
        if RECORD:
            RESULTS["versions"] = cls.sur.describe()["versions"]
            RESULTS["date"] = time.strftime("%Y-%m-%d")
            with open(os.path.join(HERE, "test_surrogate_results.json"), "w", encoding="utf-8") as fh:
                json.dump(RESULTS, fh, ensure_ascii=False, indent=1)

    # ------------------------------------------------------------------ план
    def test_01_design(self):
        Xs, info = S.design(SPACE, 64, "sobol", seed=0)
        self.assertEqual(Xs.shape, (64, 3))
        self.assertTrue(((Xs >= SPACE.lo) & (Xs <= SPACE.hi)).all())
        Xr, inf_r = S.design(SPACE, 64, "random", seed=0)
        self.assertLess(info["discrepancy_cd"], inf_r["discrepancy_cd"])     # Соболь равномернее случайного
        Xl, _ = S.design(SPACE, 50, "lhs", seed=0)
        U = SPACE.to_unit(Xl)
        for j in range(3):  # латинский гиперкуб: ровно одна точка в каждой из 50 полос по каждой оси
            self.assertEqual(sorted(np.floor(U[:, j] * 50).astype(int).tolist()), list(range(50)))
        g = S.block_groups(Xl, SPACE, n_blocks=5, mode="order")
        self.assertEqual(g.tolist(), sorted(g.tolist()))

    # ------------------------------------------------------------------ отложенная выборка
    def test_02_holdout(self):
        rep_true = {k: S.metrics(self.Yht, self.sur.predict(self.Xh, kind=k)[0], names=OUT) for k in self.sur.kinds}
        rep_noisy = S.validate(self.sur, self.Xh, self.Yh)
        RESULTS["holdout_vs_truth"] = rep_true
        RESULTS["holdout_vs_noisy"] = rep_noisy
        for k in self.sur.kinds:
            for name, rep in (("без шума", rep_true[k]), ("с шумом", rep_noisy[k])):
                print(S._fmt_metrics(f"отложенная выборка 256 точек, {k}, эталон {name}:", rep))
        for k in self.sur.kinds:
            e = rep_true[k]["eta"]
            # критерий задачи: ошибка по КПД не более 2 % (берётся максимум по 256 точкам, не среднее)
            self.assertLessEqual(e["max_rel_pct"], 2.0, k)
            # модель ближе к истине, чем один зашумлённый расчёт (СКО шума 0,003)
            self.assertLess(e["rmse"], NOISE["eta"], k)
            self.assertLessEqual(rep_true[k]["PR"]["max_rel_pct"], 2.0, k)
        # калибровка: с шумом в СКО прогноза доля |ошибка| <= 2 СКО не ниже 0,8
        self.assertGreaterEqual(rep_noisy["gp"]["eta"]["cover_2sd"], 0.8)
        # GP восстанавливает уровень шума по КПД (заложено 0,003)
        self.assertLess(abs(self.sur.noise_sd()["eta"] - NOISE["eta"]) / NOISE["eta"], 0.35)

    # ------------------------------------------------------------------ перекрёстная по блокам
    def test_03_block_cv(self):
        t0 = time.time()
        cv, g = S.block_cv(self.X, self.Y, SPACE, OUT, n_folds=5, groups="kmeans", n_blocks=20,
                           kinds=("mlp", "gp"), mlp=MLP, log=("D2",), seed=0, n_jobs=JOBS)
        dt = time.time() - t0
        RESULTS["block_cv"] = dict(folds=5, blocks=20, groups="kmeans", seconds=round(dt, 1),
                                   metrics={k: v["metrics"] for k, v in cv.items()})
        for k in cv:
            print(S._fmt_metrics(f"перекрёстная по блокам (k-means 20 блоков, 5 частей, эталон с шумом), {k}:",
                                 cv[k]["metrics"]))
        print(f"  {dt:.1f} с")
        self.assertEqual(len(np.unique(g)), 20)
        for k in cv:
            self.assertLessEqual(cv[k]["metrics"]["eta"]["max_rel_pct"], 2.0, k)

    # ------------------------------------------------------------------ оптимум по модели
    def _run_case(self, name, kind):
        c = CASES[name]
        r = S.optimize(self.sur, "eta", True, c["cons"], fixed={"T1": T_ISO}, kind=kind, tol=0.005, seed=0)
        x = r["x"]
        eta, pr, d2 = stage(x["U2"], x["m"], x["T1"])
        direct = dict(eta=float(eta), PR=float(pr), D2=float(d2), U2=x["U2"], m=x["m"])
        opt = true_optimum(c["cons"])
        model_eta = r["pred"][kind]["eta"]["mean"]
        out = dict(case=c["note"], constraints=c["cons"], kind=kind, x=x, model=r["pred"][kind], direct=direct,
                   true_opt=opt, gap_model_vs_direct_pp=100 * abs(model_eta - eta),
                   gap_direct_vs_true_opt_pp=100 * (opt["eta"] - eta), flags=r["flags"], at_edge=r["at_edge"],
                   notes=r["notes"], nn_ratio=r["nn_ratio"], de=r["de"])
        RESULTS.setdefault("optimize", {})[f"{name}/{kind}"] = out
        print(f"\n{name} [{kind}] {c['cons']}: U2={x['U2']:.1f} m={x['m']:.2f} | модель eta={model_eta:.4f} "
              f"±{r['pred'][kind]['eta']['sd']:.4f}; прямой eta={eta:.4f} PR={pr:.3f} D2={d2:.3f} | "
              f"истинный оптимум eta={opt['eta']:.4f} (U2={opt['U2']:.1f} m={opt['m']:.2f}) | "
              f"флаги {r['flags']} край {r['at_edge']} | DE {r['de']['nfev']} выз., {r['de']['seconds']} с")
        # критерий задачи: прогноз модели в найденной точке и прямой расчёт — не более 1 п. п.
        self.assertLessEqual(out["gap_model_vs_direct_pp"], 1.0)
        # найденная точка хуже истинного оптимума не более чем на 1 п. п.
        self.assertLessEqual(out["gap_direct_vs_true_opt_pp"], 1.0)
        self.assertTrue(check_cons(direct, c["cons"]), direct)
        self.assertFalse(r["uncertain"], r["notes"])
        return out

    def test_04_optimize_mlp(self):
        for name in CASES:
            self._run_case(name, "mlp")

    def test_05_optimize_gp(self):
        for name in CASES:
            self._run_case(name, "gp")

    # ------------------------------------------------------------------ только изученная область
    def test_06_studied_region(self):
        """данные только при U2 >= 520 м/с; истинный оптимум случая A (U2 около 505) снаружи"""
        sel = self.X[:, 0] >= 520.0
        sub = S.Surrogate(SPACE, OUT, kinds=("mlp", "gp"), mlp=MLP, log=("D2",), seed=0, n_jobs=JOBS)
        sub.fit(self.X[sel], self.Y[sel])
        r = S.optimize(sub, "eta", True, CASES["A_core130"]["cons"], fixed={"T1": T_ISO}, kind="mlp", tol=0.005)
        x = r["x"]
        eta, pr, d2 = stage(x["U2"], x["m"], x["T1"])
        opt = true_optimum(CASES["A_core130"]["cons"])
        # на краю и вне данных неопределённость больше, чем в середине изученной области
        _, sd_out = sub.predict([[400.0, 12.0, T_ISO]], kind="gp")
        _, sd_in = sub.predict([[585.0, 12.0, T_ISO]], kind="gp")
        RESULTS["studied_region"] = dict(n=int(sel.sum()), x=x, direct=dict(eta=eta, PR=pr, D2=d2), true_opt=opt,
                                         at_edge=r["at_edge"], flags=r["flags"], notes=r["notes"],
                                         gp_sd_eta_outside=float(sd_out[0, 0]), gp_sd_eta_inside=float(sd_in[0, 0]))
        print(f"\nизученная область U2>={sub.box_lo[0]:.1f}: найдено U2={x['U2']:.1f}, прямой eta={eta:.4f}, "
              f"истинный оптимум eta={opt['eta']:.4f} при U2={opt['U2']:.1f}; край {r['at_edge']}; "
              f"СКО GP вне/внутри {sd_out[0, 0]:.4f}/{sd_in[0, 0]:.4f}")
        self.assertGreaterEqual(x["U2"], sub.box_lo[0])
        self.assertIn("U2:низ", r["at_edge"])
        self.assertGreater(sd_out[0, 0], 2.0 * sd_in[0, 0])
        # заведомо малый допуск -> флаг неопределённости
        r2 = S.optimize(sub, "eta", True, CASES["A_core130"]["cons"], fixed={"T1": T_ISO}, kind="mlp", tol=1e-6)
        self.assertTrue(r2["uncertain"])
        self.assertIn("sigma", r2["flags"])
        with self.assertRaises(ValueError):  # закреплять вход вне изученной области нельзя
            S.optimize(sub, "eta", True, [], fixed={"U2": 400.0})

    # ------------------------------------------------------------------ следующая порция по EI
    def test_07_next_batch(self):
        cons = CASES["A_core130"]["cons"]
        Xn, rows = S.propose(self.sur, 8, "eta", True, cons, fixed={"T1": T_ISO}, kind="gp", seed=0)
        self.assertEqual(Xn.shape, (8, 3))
        self.assertTrue(np.allclose(Xn[:, 2], T_ISO))
        self.assertTrue(((Xn >= SPACE.lo - 1e-9) & (Xn <= SPACE.hi + 1e-9)).all())
        self.assertEqual(len({tuple(np.round(x, 6)) for x in Xn}), 8)
        self.assertGreater(rows[0]["acq"], 0.0)
        self.assertEqual(rows[0]["mode"], "ei")
        # порция не вырождается в повторы одной точки
        Ub = self.sur._ub(Xn)
        dmin = min(np.linalg.norm(Ub[i] - Ub[j]) for i in range(8) for j in range(i))
        self.assertGreaterEqual(dmin, 0.5 * self.sur.nn_med * 0.999)
        Yn_t = truth(Xn)
        # порция ложится у границы ограничения PR = 4,25, где и лежит оптимум
        Xr, _ = S.design(SPACE, 8, "random", seed=11)
        Xr[:, 2] = T_ISO
        d_ei = float(np.median(np.abs(Yn_t[:, 1] - 4.25)))
        d_rand = float(np.median(np.abs(truth(Xr)[:, 1] - 4.25)))
        # дорасчёт: добавить порцию (с шумом) и снова найти оптимум
        X2 = np.vstack([self.X, Xn])
        Y2 = np.vstack([self.Y, noisy(Yn_t, seed=9)])
        sur2 = S.Surrogate(SPACE, OUT, kinds=("mlp", "gp"), mlp=MLP, log=("D2",), seed=0, n_jobs=JOBS).fit(X2, Y2)
        before = RESULTS.get("optimize", {}).get("A_core130/gp") or self._run_case("A_core130", "gp")
        r = S.optimize(sur2, "eta", True, cons, fixed={"T1": T_ISO}, kind="gp", tol=0.005)
        eta2 = stage(r["x"]["U2"], r["x"]["m"], T_ISO)[0]
        RESULTS["next_batch"] = dict(k=8, rows=rows, min_dist_unit=float(dmin), nn_med_unit=self.sur.nn_med, true=Yn_t.tolist(), median_abs_PR_minus_target_ei=d_ei,
                                     median_abs_PR_minus_target_random=d_rand,
                                     after=dict(x=r["x"], model_eta=r["pred"]["gp"]["eta"]["mean"], direct_eta=eta2),
                                     before_direct_eta=before["direct"]["eta"])
        for rw, yt in zip(rows, Yn_t):
            print(f"  {rw['mode']:4s} U2={rw['x']['U2']:6.1f} m={rw['x']['m']:5.2f} acq={rw['acq']:.3g} "
                  f"pof={rw['pof']:.2f} | прямой eta={yt[0]:.4f} PR={yt[1]:.3f}")
        print(f"порция EI: |PR-4,25| медиана {d_ei:.3f} (случайные точки {d_rand:.3f}); "
              f"после дорасчёта прямой eta {eta2:.4f} (было {before['direct']['eta']:.4f}), "
              f"модель {r['pred']['gp']['eta']['mean']:.4f}")
        self.assertLess(d_ei, d_rand)
        self.assertLessEqual(100 * abs(r["pred"]["gp"]["eta"]["mean"] - eta2), 1.0)
        self.assertGreaterEqual(eta2, before["direct"]["eta"] - 0.001)

    # ------------------------------------------------------------------ командная строка
    def test_08_cli(self):
        with tempfile.TemporaryDirectory() as td:
            prob = os.path.join(td, "problem.json")
            with open(prob, "w") as fh:
                json.dump(dict(inputs=SPACE.bounds(), outputs=OUT), fh)
            plan = os.path.join(td, "plan.csv")
            self.assertEqual(S.main(["plan", "--problem", prob, "--n", "64", "--method", "sobol", "--out", plan]), 0)
            Xp, _, _, _ = S.read_table(plan, SPACE.names, [])
            self.assertEqual(Xp.shape, (64, 3))
            runs = os.path.join(td, "runs.csv")
            Yp = noisy(truth(Xp), seed=5)
            S.write_table(runs, Xp, SPACE.names, {o: Yp[:, j] for j, o in enumerate(OUT)})
            with open(runs, "a") as fh:  # неудачный прогон: пустой показатель
                fh.write("500,10,288.15,,4.1,0.43\n")
            model = os.path.join(td, "model.joblib")
            rep = os.path.join(td, "report.json")
            self.assertEqual(S.main(["fit", "--problem", prob, "--data", runs, "--kinds", "gp", "--folds", "3",
                                     "--n-blocks", "6", "--log", "D2", "--jobs", "1", "--out", model,
                                     "--report", rep]), 0)
            with open(rep) as fh:
                R = json.load(fh)
            self.assertEqual(R["dropped"], 1)
            self.assertEqual(R["n"], 64)
            opt = os.path.join(td, "opt.json")
            self.assertEqual(S.main(["optimize", "--model", model, "--objective", "eta", "--kind", "gp",
                                     "--constraint", "PR>=4.25", "--fixed", "T1=288.15", "--out", opt]), 0)
            with open(opt) as fh:
                O = json.load(fh)
            self.assertAlmostEqual(O["x"]["T1"], T_ISO)
            nxt = os.path.join(td, "next.csv")
            self.assertEqual(S.main(["next", "--model", model, "--objective", "eta", "--k", "4",
                                     "--constraint", "PR>=4.25", "--fixed", "T1=288.15", "--out", nxt]), 0)
            Xn, _, _, _ = S.read_table(nxt, SPACE.names, [])
            self.assertEqual(Xn.shape, (4, 3))
            sur = S.Surrogate.load(model)
            self.assertEqual(sur.outputs, OUT)


if __name__ == "__main__":
    unittest.main(verbosity=2)
