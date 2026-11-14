#!/usr/bin/env python3
"""Модель-заменитель (суррогат) дорогого расчёта и оптимизация по ней.

Узел платформы трека «Простая литая ГТУ-заменитель». Дорогой расчёт (SU2, OpenFOAM,
CalculiX) даёт таблицу «параметры -> показатели»; модуль:

  1) строит план расчётов в заданных границах (Sobol / латинский гиперкуб, scipy.stats.qmc);
  2) обучает две модели: ансамбль MLP (scikit-learn) и гауссовский процесс (базовая линия);
  3) проверяет их перекрёстно по блокам (пространственные кластеры или порции прогонов)
     и на отложенной выборке, с ошибкой по каждому показателю;
  4) ищет оптимум (scipy differential_evolution) только внутри изученной области:
     штраф за подход к краю, штраф за «дыры» в данных, осторожная оценка mu -/+ kappa*sigma,
     флаг неопределённости (дисперсия GP, разброс ансамбля MLP, расхождение MLP и GP);
  5) выбирает следующую порцию точек для дорасчёта по ожидаемому улучшению (EI)
     с учётом вероятности выполнить ограничения; порция — методом «Kriging believer».

Командная строка:
  python surrogate.py plan     --problem problem.json --n 64 --method sobol --out plan.csv
  python surrogate.py fit      --problem problem.json --data runs.csv [--holdout h.csv] --out model.joblib
  python surrogate.py optimize --model model.joblib --objective eta --max --constraint "PR>=4.25" --fixed T1=288.15
  python surrogate.py next     --model model.joblib --k 8 --objective eta --max --constraint "PR>=4.25"

problem.json: {"inputs": {"U2": [350, 650], "m": [1.5, 25]}, "outputs": ["eta", "PR"]}
Подробности и проверка — README.md рядом, тест — test_surrogate.py.
"""
from __future__ import annotations

import argparse
import csv
import functools
import inspect
import json
import math
import os
import re
import sys
import time
import warnings

import numpy as np
from scipy.optimize import differential_evolution
from scipy.spatial import cKDTree
from scipy.stats import norm, qmc

import joblib
import sklearn
from sklearn.cluster import KMeans
from sklearn.exceptions import ConvergenceWarning
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.model_selection import GroupKFold
from sklearn.neural_network import MLPRegressor
from threadpoolctl import ThreadpoolController

__all__ = [
    "Space", "design", "Surrogate", "GPModel", "MLPEnsemble", "metrics", "block_cv",
    "validate", "optimize", "propose", "expected_improvement", "parse_constraint",
    "read_table", "write_table",
]


# Матрицы малы (сотни точек): многопоточный BLAS только мешает, а на занятой машине (соседние
# расчёты) замедляет обучение в 4-14 раз (замер в README). Число потоков BLAS — SURROGATE_BLAS_THREADS.
_BLAS = int(os.environ.get("SURROGATE_BLAS_THREADS", "1"))
_TPC = None


def _blas_limited(fn):
    @functools.wraps(fn)
    def wrap(*a, **k):
        global _TPC
        if _BLAS <= 0:
            return fn(*a, **k)
        if _TPC is None:
            _TPC = ThreadpoolController()
        with _TPC.limit(limits=_BLAS, user_api="blas"):
            return fn(*a, **k)
    return wrap


def _rng_kw(fn, seed):
    """scipy >= 1.15 принимает rng=, старые версии — seed="""
    try:
        params = inspect.signature(fn).parameters
    except (TypeError, ValueError):
        params = {}
    return {"rng": seed} if "rng" in params else {"seed": seed}


# ======================================================================= пространство параметров
class Space:
    """границы параметров; внутри модуля всё считается в единичном кубе [0, 1]^d"""

    def __init__(self, bounds):
        if not bounds:
            raise ValueError("пустой список параметров")
        self.names = list(bounds)
        self.lo = np.array([float(bounds[k][0]) for k in self.names])
        self.hi = np.array([float(bounds[k][1]) for k in self.names])
        bad = [k for k, a, b in zip(self.names, self.lo, self.hi) if not b > a]
        if bad:
            raise ValueError(f"верхняя граница не больше нижней: {bad}")

    @property
    def d(self):
        return len(self.names)

    def bounds(self):
        return {k: [float(a), float(b)] for k, a, b in zip(self.names, self.lo, self.hi)}

    def to_unit(self, X):
        return (np.asarray(X, float) - self.lo) / (self.hi - self.lo)

    def from_unit(self, U):
        return self.lo + np.asarray(U, float) * (self.hi - self.lo)


def design(space, n, method="sobol", seed=0):
    """План из n точек в границах space.
    sobol  — скремблированная последовательность Соболя; равномерность гарантирована при n = 2^k,
             иначе берутся первые n точек из 2^k (предупреждение);
    lhs    — латинский гиперкуб с улучшением центрированной L2-неравномерности (random-cd);
    random — равномерные случайные точки (для отложенной выборки).
    Возвращает (X, info), info — метод, неравномерность (centered L2 discrepancy) плана."""
    n = int(n)
    if n < 2:
        raise ValueError("n >= 2")
    d = space.d
    if method == "sobol":
        eng = qmc.Sobol(d, scramble=True, **_rng_kw(qmc.Sobol, seed))
        k = int(math.ceil(math.log2(n)))
        if 2 ** k != n:
            warnings.warn(f"Sobol: n={n} не степень двойки, взяты первые {n} из {2 ** k} точек")
        U = eng.random_base2(k)[:n]
    elif method == "lhs":
        eng = qmc.LatinHypercube(d, optimization="random-cd", **_rng_kw(qmc.LatinHypercube, seed))
        U = eng.random(n)
    elif method == "random":
        U = np.random.default_rng(seed).random((n, d))
    else:
        raise ValueError(f"неизвестный метод плана: {method}")
    info = dict(method=method, n=n, seed=seed, discrepancy_cd=float(qmc.discrepancy(U)))
    return space.from_unit(U), info


# ======================================================================= модели
def _gp_kernel(d, nu):
    return (ConstantKernel(1.0, (1e-3, 1e3))
            * Matern(length_scale=np.full(d, 0.5), length_scale_bounds=(1e-2, 1e2), nu=nu)
            + WhiteKernel(1e-3, (1e-9, 1e0)))


class GPModel:
    """по гауссовскому процессу на каждый показатель: Matern(nu) с отдельной длиной по каждой оси
    + белый шум (оценка шума расчёта). Работает со стандартизованными показателями."""

    name = "gp"

    def __init__(self, nu=2.5, restarts=3, seed=0):
        self.nu, self.restarts, self.seed = nu, restarts, seed
        self.gps, self.noise_ = [], np.zeros(0)
        self.white = True       # шум входит в ядро (WhiteKernel); у копии condition() — в alpha

    @_blas_limited
    def fit(self, U, Ys):
        self.gps = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            for j in range(Ys.shape[1]):
                gp = GaussianProcessRegressor(_gp_kernel(U.shape[1], self.nu), normalize_y=False,
                                              n_restarts_optimizer=self.restarts, random_state=self.seed + j)
                gp.fit(U, Ys[:, j])
                self.gps.append(gp)
        self.noise_ = np.array([gp.kernel_.k2.noise_level for gp in self.gps])
        return self

    def noise_var(self):
        """оценка дисперсии шума по каждому показателю (в стандартизованных единицах)"""
        return self.noise_

    @_blas_limited
    def predict(self, U, latent=True):
        """среднее и СКО; latent=True — СКО самой функции без шума расчёта"""
        mu, sd = [], []
        for gp, nv in zip(self.gps, self.noise_):
            m, s = gp.predict(U, return_std=True)
            if self.white:
                v = s ** 2 - (nv if latent else 0.0)
            else:
                v = s ** 2 + (0.0 if latent else nv)
            mu.append(m)
            sd.append(np.sqrt(np.clip(v, 0.0, None)))
        return np.column_stack(mu), np.column_stack(sd)

    @_blas_limited
    def condition(self, U_add, Ys_add):
        """копия модели с добавленными точками при неизменных гиперпараметрах (Kriging believer).
        Добавленные точки считаются точными (дорогой расчёт детерминирован: повтор той же точки
        ничего не даёт), исходные — с оценённым шумом; шум переносится из ядра в alpha."""
        new = GPModel(self.nu, 0, self.seed)
        new.white, new.noise_ = False, self.noise_.copy()
        for j, gp in enumerate(self.gps):
            k = gp.kernel_.k1 if self.white else gp.kernel_
            n0 = len(gp.X_train_)
            a0 = (np.full(n0, self.noise_[j]) if self.white else np.asarray(gp.alpha, float) * np.ones(n0))
            alpha = np.concatenate([a0 + 1e-10, np.full(len(U_add), 1e-8)])
            g = GaussianProcessRegressor(k, alpha=alpha, optimizer=None, normalize_y=False)
            g.fit(np.vstack([gp.X_train_, U_add]), np.concatenate([gp.y_train_, Ys_add[:, j]]))
            new.gps.append(g)
        return new

    def describe(self):
        return [str(gp.kernel_) for gp in self.gps]


@_blas_limited
def _fit_member(Ub, Yb, hidden, alpha, activation, max_iter, seed):
    net = MLPRegressor(hidden_layer_sizes=hidden, activation=activation, solver="lbfgs", alpha=alpha,
                       max_iter=max_iter, tol=1e-8, random_state=seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        net.fit(2.0 * Ub - 1.0, Yb)
    return net


def _alpha_job(U, Yc, tr, te, hidden, alpha, activation, max_iter, seed):
    net = _fit_member(U[tr], Yc[tr], hidden, alpha, activation, max_iter, seed)
    e = np.asarray(net.predict(2.0 * U[te] - 1.0)).reshape(len(te), -1) - Yc[te].reshape(len(te), -1)
    return float(np.sum(e ** 2))


class MLPEnsemble:
    """ансамбль многослойных перцептронов (lbfgs); у каждого члена своя начальная инициализация и,
    при bootstrap=True, своя бутстреп-выборка. Разброс членов — мера неопределённости.
    alpha='auto' — L2-регуляризация подбирается по каждому показателю внутренней 4-кратной
    перекрёстной проверкой из alpha_grid (шумному показателю нужна сильнее, гладкому — слабее)."""

    name = "mlp"

    def __init__(self, members=8, hidden=(32, 32), alpha="auto", activation="tanh", bootstrap=True,
                 per_output=True, alpha_grid=(1e-4, 1e-3, 1e-2, 1e-1, 1.0), max_iter=3000, seed=0, n_jobs=1):
        self.members, self.hidden, self.alpha = members, tuple(hidden), alpha
        self.activation, self.bootstrap, self.max_iter = activation, bootstrap, max_iter
        self.per_output, self.alpha_grid, self.seed, self.n_jobs = per_output, tuple(alpha_grid), seed, n_jobs
        self.nets, self.alpha_ = [], []

    def _par(self, fn, jobs):
        if self.n_jobs == 1:
            return [fn(*j) for j in jobs]
        return joblib.Parallel(n_jobs=self.n_jobs)(joblib.delayed(fn)(*j) for j in jobs)

    def _select_alpha(self, U, Ys, cols):
        from sklearn.model_selection import KFold
        if len(U) < 4:
            raise ValueError(f"alpha='auto' (4-кратная проверка) требует не меньше 4 точек, а их {len(U)}; "
                             "задайте alpha числом")
        splits = list(KFold(4, shuffle=True, random_state=self.seed).split(U))
        jobs, keys = [], []
        for ci, c in enumerate(cols):
            Yc = Ys[:, c[0]] if len(c) == 1 else Ys[:, c]
            for a in self.alpha_grid:
                for tr, te in splits:
                    jobs.append((U, Yc, tr, te, self.hidden, a, self.activation, self.max_iter, self.seed))
                    keys.append((ci, a))
        sse = self._par(_alpha_job, jobs)
        table = {}
        for (ci, a), e in zip(keys, sse):
            table[(ci, a)] = table.get((ci, a), 0.0) + e
        self.alpha_cv_ = {f"{ci}:{a:g}": v / len(U) for (ci, a), v in table.items()}
        return [min(self.alpha_grid, key=lambda a: table[(ci, a)]) for ci in range(len(cols))]

    def fit(self, U, Ys):
        """per_output=True — своя сеть на каждый показатель (у показателей разная гладкость и шум),
        иначе одна сеть с несколькими выходами"""
        n, p = Ys.shape
        cols = [[j] for j in range(p)] if self.per_output else [list(range(p))]
        self.alpha_ = (self._select_alpha(U, Ys, cols) if self.alpha == "auto"
                       else [float(self.alpha)] * len(cols))
        jobs = []
        for i in range(self.members):
            rng = np.random.default_rng(self.seed + 1000 * i)
            idx = rng.integers(0, n, n) if self.bootstrap else np.arange(n)
            for ci, c in enumerate(cols):
                yb = Ys[np.ix_(idx, c)]
                jobs.append((U[idx], yb[:, 0] if len(c) == 1 else yb, self.hidden, self.alpha_[ci],
                             self.activation, self.max_iter, self.seed + i))
        nets = self._par(_fit_member, jobs)
        self.nets = [nets[i * len(cols):(i + 1) * len(cols)] for i in range(self.members)]
        return self

    @_blas_limited
    def predict(self, U, latent=True):
        Z = 2.0 * U - 1.0
        P = np.stack([np.column_stack([np.asarray(net.predict(Z)).reshape(len(U), -1) for net in member])
                      for member in self.nets])
        sd = P.std(axis=0, ddof=1) if len(self.nets) > 1 else np.zeros(P.shape[1:])
        return P.mean(axis=0), sd

    def describe(self):
        return dict(members=self.members, hidden=list(self.hidden), alpha=self.alpha, alpha_used=self.alpha_,
                    activation=self.activation, bootstrap=self.bootstrap, per_output=self.per_output,
                    iters_max=max(int(net.n_iter_) for member in self.nets for net in member))


class Surrogate:
    """модель-заменитель: входы -> показатели. Обучает выбранные виды моделей ('mlp', 'gp')
    на одной таблице; хранит обучающие данные (изученная область, плотность точек).
    log — показатели, которые учатся в логарифме (положительные величины с разбросом в разы:
    диаметр, расход, напряжение); прогноз — exp(среднего логарифма), СКО — exp(mu)*СКО логарифма."""

    def __init__(self, space, outputs, kinds=("mlp", "gp"), mlp=None, gp=None, log=(), seed=0, n_jobs=1):
        self.space, self.outputs = space, list(outputs)
        self.kinds = tuple(kinds)
        bad = [o for o in log if o not in self.outputs]
        if bad:
            raise ValueError(f"log: нет таких показателей {bad}")
        self.log = np.array([o in set(log) for o in self.outputs])
        self.mlp_kw, self.gp_kw = dict(mlp or {}), dict(gp or {})
        self.seed, self.n_jobs = seed, n_jobs
        self.models, self.report = {}, {}

    # ----------------------------------------------------------------- обучение
    def fit(self, X, Y):
        X, Y = np.asarray(X, float), np.asarray(Y, float).reshape(len(X), -1)
        if Y.shape[1] != len(self.outputs):
            raise ValueError("число столбцов Y не совпадает с числом показателей")
        if not (np.isfinite(X).all() and np.isfinite(Y).all()):
            raise ValueError("в данных NaN/inf — отбросьте неудачные прогоны до обучения (read_table это делает)")
        if (Y[:, self.log] <= 0).any():
            raise ValueError("показатель из списка log должен быть положительным")
        self.X, self.Y = X, Y
        self.y_scale = Y.std(0)                      # масштаб показателя для штрафов и допусков
        self.y_scale[self.y_scale == 0] = 1.0
        Yt = Y.copy()
        Yt[:, self.log] = np.log(Y[:, self.log])
        self.y_mu, self.y_sd = Yt.mean(0), Yt.std(0)
        self.y_sd[self.y_sd == 0] = 1.0
        U, Ys = self.space.to_unit(X), (Yt - self.y_mu) / self.y_sd
        self.models = {}
        for kind in self.kinds:
            if kind == "mlp":
                kw = dict(seed=self.seed, n_jobs=self.n_jobs)
                kw.update(self.mlp_kw)
                self.models[kind] = MLPEnsemble(**kw).fit(U, Ys)
            elif kind == "gp":
                kw = dict(seed=self.seed)
                kw.update(self.gp_kw)
                self.models[kind] = GPModel(**kw).fit(U, Ys)
            else:
                raise ValueError(f"неизвестный вид модели: {kind}")
        # изученная область: прямоугольник по данным и типичное расстояние между соседями.
        # Расстояние — по различным точкам: повторы (пересчёт той же точки) дают нулевого соседа,
        # и при половине повторов медиана стала бы 0 — штраф за дыры приковал бы оптимум к данным.
        self.box_lo, self.box_hi = X.min(0), X.max(0)
        Ub = self._ub(X)
        self.tree = cKDTree(Ub)
        Uu = np.unique(np.round(Ub, 12), axis=0)
        self.n_unique = int(len(Uu))
        if self.n_unique >= 2:
            self.nn_med = float(np.median(cKDTree(Uu).query(Uu, k=2)[0][:, 1]))
        else:
            self.nn_med = 1.0
        return self

    def _ub(self, X):
        """координаты в единичном кубе изученной области"""
        span = np.where(self.box_hi > self.box_lo, self.box_hi - self.box_lo, 1.0)
        return (np.asarray(X, float) - self.box_lo) / span

    # ----------------------------------------------------------------- прогноз
    def predict(self, X, kind=None, latent=True):
        """среднее и СКО (в единицах показателей), массивы (n, p)"""
        kind = kind or self.kinds[0]
        X = np.atleast_2d(np.asarray(X, float))
        return self._inverse(*self.models[kind].predict(self.space.to_unit(X), latent=latent))

    def _inverse(self, mu_s, sd_s):
        """из стандартизованных (и, для log, логарифмированных) единиц модели — в единицы показателя"""
        mu, sd = mu_s * self.y_sd + self.y_mu, sd_s * self.y_sd
        if self.log.any():
            mu, sd = mu.copy(), sd.copy()
            mu[:, self.log] = np.exp(mu[:, self.log])
            sd[:, self.log] = mu[:, self.log] * sd[:, self.log]
        return mu, sd

    def noise_sd(self):
        """оценка СКО шума расчёта по GP, в единицах показателей"""
        if "gp" not in self.models:
            return None
        s = np.sqrt(self.models["gp"].noise_var()) * self.y_sd
        s = np.where(self.log, s * np.exp(self.y_mu), s)   # для log — относительный шум, умноженный на типичное значение
        return dict(zip(self.outputs, s.tolist()))

    def nn_ratio(self, X):
        """расстояние до ближайшей обучающей точки / медианное расстояние между соседями"""
        return self.tree.query(self._ub(np.atleast_2d(X)), k=1)[0] / max(self.nn_med, 1e-12)

    def save(self, path):
        joblib.dump(self, path)

    @staticmethod
    def load(path):
        return joblib.load(path)

    def describe(self):
        return dict(inputs=self.space.bounds(), outputs=self.outputs, n=len(self.X),
                    n_unique=getattr(self, "n_unique", None), log=[o for o, f in zip(self.outputs, self.log) if f],
                    studied=dict(zip(self.space.names, zip(self.box_lo.tolist(), self.box_hi.tolist()))),
                    nn_median_unit=self.nn_med, noise_sd_gp=self.noise_sd(),
                    models={k: m.describe() for k, m in self.models.items()},
                    versions=dict(sklearn=sklearn.__version__, numpy=np.__version__,
                                  scipy=__import__("scipy").__version__))


# ======================================================================= проверка
def metrics(Y_true, Y_pred, Y_std=None, names=None):
    """ошибка по каждому показателю: RMSE, MAE, максимум модуля, MAPE и максимум относительной
    ошибки (%), R2; при Y_std — доля точек, где |ошибка| <= 2 СКО (калибровка неопределённости)"""
    Y_true = np.asarray(Y_true, float).reshape(len(Y_true), -1)
    Y_pred = np.asarray(Y_pred, float).reshape(Y_true.shape)
    names = names or [f"y{j}" for j in range(Y_true.shape[1])]
    out = {}
    for j, nm in enumerate(names):
        e = Y_pred[:, j] - Y_true[:, j]
        denom = np.where(np.abs(Y_true[:, j]) > 1e-12, np.abs(Y_true[:, j]), np.nan)
        rel = np.abs(e) / denom * 100.0
        ss = float(np.sum((Y_true[:, j] - Y_true[:, j].mean()) ** 2))
        r = dict(n=int(len(e)), rmse=float(np.sqrt(np.mean(e ** 2))), mae=float(np.mean(np.abs(e))),
                 max_abs=float(np.max(np.abs(e))), mape_pct=float(np.nanmean(rel)),
                 p95_rel_pct=float(np.nanpercentile(rel, 95)), max_rel_pct=float(np.nanmax(rel)),
                 r2=float(1.0 - np.sum(e ** 2) / ss) if ss > 0 else float("nan"))
        if Y_std is not None:
            s = np.asarray(Y_std, float).reshape(Y_true.shape)[:, j]
            r["cover_2sd"] = float(np.mean(np.abs(e) <= 2.0 * s))
        out[nm] = r
    return out


def validate(sur, X, Y, latent=False):
    """ошибка каждой модели на отложенной выборке (X, Y); latent=False — в СКО входит шум"""
    return {k: metrics(Y, *sur.predict(X, kind=k, latent=latent), names=sur.outputs) for k in sur.kinds}


def _fold(space, outputs, kinds, mlp, gp, log, seed, X, Y, tr, te):
    s = Surrogate(space, outputs, kinds=kinds, mlp=dict(mlp or {}, n_jobs=1), gp=gp, log=log, seed=seed, n_jobs=1)
    s.fit(X[tr], Y[tr])
    return te, {k: s.predict(X[te], kind=k, latent=False) for k in kinds}


def block_groups(X, space, n_blocks=20, mode="kmeans", seed=0):
    """метки блоков: kmeans — пространственные кластеры (выбрасывается целый район);
    order — подряд идущие порции строк (порции прогонов)"""
    n = len(X)
    if mode == "kmeans":
        return KMeans(n_clusters=n_blocks, n_init=10, random_state=seed).fit_predict(space.to_unit(X))
    if mode == "order":
        return np.arange(n) * n_blocks // n
    raise ValueError(f"неизвестный способ блоков: {mode}")


def block_cv(X, Y, space, outputs, n_folds=5, groups="kmeans", n_blocks=20, kinds=("mlp", "gp"),
             mlp=None, gp=None, log=(), seed=0, n_jobs=1):
    """перекрёстная проверка по блокам: точки делятся на блоки (groups — 'kmeans', 'order' или
    массив меток), блоки — на n_folds частей; модель учится без части и предсказывает её.
    Возвращает {вид: {'metrics': ..., 'pred': (mu, sd)}} и метки блоков."""
    X, Y = np.asarray(X, float), np.asarray(Y, float).reshape(len(X), -1)
    g = block_groups(X, space, n_blocks, groups, seed) if isinstance(groups, str) else np.asarray(groups)
    if len(np.unique(g)) < n_folds:
        raise ValueError(f"блоков {len(np.unique(g))} меньше числа частей {n_folds}")
    splits = list(GroupKFold(n_splits=n_folds).split(X, Y, g))
    res = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(_fold)(space, outputs, kinds, mlp, gp, log, seed, X, Y, tr, te) for tr, te in splits)
    out = {}
    for k in kinds:
        mu, sd = np.zeros_like(Y), np.zeros_like(Y)
        for te, pr in res:
            mu[te], sd[te] = pr[k]
        out[k] = dict(metrics=metrics(Y, mu, sd, outputs), pred=(mu, sd))
    return out, g


# ======================================================================= оптимизация
_CON = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*(>=|<=)\s*"
                  r"([-+]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][-+]?[0-9]+)?)\s*$")


def parse_constraint(c):
    """'PR>=4.25' -> ('PR', '>=', 4.25); кортеж пропускается как есть"""
    if isinstance(c, (tuple, list)):
        return str(c[0]), str(c[1]), float(c[2])
    m = _CON.match(c)
    if not m:
        raise ValueError(f"ограничение не разобрано: {c!r} (формат: имя>=число или имя<=число)")
    return m.group(1), m.group(2), float(m.group(3))


def _check_task(sur, objective, kind, constraints):
    """проверка цели, вида модели и имён в ограничениях; возвращает разобранные ограничения"""
    if objective not in sur.outputs:
        raise ValueError(f"цель {objective!r} не среди показателей {sur.outputs}")
    if kind not in sur.models:
        raise ValueError(f"модель {kind!r} не обучена")
    cons = [parse_constraint(c) for c in constraints]
    for name, _, _ in cons:
        if name not in sur.outputs and name not in sur.space.names:
            raise ValueError(f"ограничение по неизвестной величине {name!r}")
    return cons


def _full(sur, fixed):
    fixed = dict(fixed or {})
    names = sur.space.names
    unknown = [k for k in fixed if k not in names]
    if unknown:
        raise ValueError(f"нет таких входов: {unknown}")
    for k, v in fixed.items():
        j = names.index(k)
        if not (sur.box_lo[j] <= v <= sur.box_hi[j]):
            raise ValueError(f"{k}={v} вне изученной области [{sur.box_lo[j]:g}, {sur.box_hi[j]:g}]")
    free = [j for j, k in enumerate(names) if k not in fixed]
    base = np.array([fixed.get(k, np.nan) for k in names], float)
    return free, base


def _evaluate(sur, X, kind, objective, sign, cons, kappa):
    """осторожная оценка цели и нарушения ограничений (в долях СКО данных) для массива точек X"""
    mu, sd = sur.predict(X, kind=kind)
    io = sur.outputs.index(objective)
    val = mu[:, io] - sign * kappa * sd[:, io]
    viol = np.zeros(len(X))
    for name, op, c in cons:
        if name in sur.outputs:
            j = sur.outputs.index(name)
            g, s, scale = mu[:, j], sd[:, j], sur.y_scale[j]
        else:
            j = sur.space.names.index(name)
            g, s, scale = X[:, j], 0.0, sur.space.hi[j] - sur.space.lo[j]
        gc = g - kappa * s if op == ">=" else g + kappa * s
        v = (c - gc) if op == ">=" else (gc - c)
        viol += np.clip(v / scale, 0.0, None)
    return val, viol, mu, sd


def optimize(sur, objective, maximize=True, constraints=(), fixed=None, kind="mlp", kappa=1.0,
             edge=0.05, w_edge=1.0, gap_ratio=2.0, w_gap=1.0, w_con=1e3, tol=None,
             seed=0, maxiter=300, popsize=30, polish=True):
    """Оптимум по модели (differential_evolution) только внутри изученной области.

    Поиск идёт в прямоугольнике обучающих данных (не номинальных границ). Штрафуемая функция:
      -sign*(mu - sign*kappa*sigma)/s_obj
      + w_con * сумма нарушений ограничений (ограничения тоже осторожные: mu -/+ kappa*sigma)
      + w_edge * сумма по свободным осям (1 - расстояние_до_края/edge)^2 в полосе edge у края
      + w_gap  * (nn/nn_med - gap_ratio)^2, если ближайшая обучающая точка дальше gap_ratio
               типичных расстояний (дыра в данных).
    s_obj — СКО цели по данным. fixed — входы, закреплённые на значении (например, T1 = ISO).
    Флаги: sigma — СКО модели в оптимуме больше tol; disagree — MLP и GP расходятся больше tol;
    gap — дыра в данных; constraint — ограничения не выполнены по среднему прогнозу. Отдельно
    at_edge — оси, по которым без штрафа за край модель ушла бы дальше (оптимум упирается в край
    изученной области; истинный может лежать снаружи). tol — допуск цели (см. _tol).
    Возвращает словарь: точка, прогноз обеих моделей, флаги, сведения о поиске."""
    cons = _check_task(sur, objective, kind, constraints)
    sign = 1.0 if maximize else -1.0
    free, base = _full(sur, fixed)
    if not free:
        raise ValueError("все входы закреплены — оптимизировать нечего")
    io = sur.outputs.index(objective)
    s_obj = float(sur.y_scale[io])

    def assemble(Xf):  # Xf: (S, n_free)
        X = np.tile(base, (len(Xf), 1))
        X[:, free] = Xf
        return X

    def penal(X, with_edge=True):
        val, viol, _, _ = _evaluate(sur, X, kind, objective, sign, cons, kappa)
        J = -sign * val / s_obj + w_con * (viol + viol ** 2)
        if with_edge:
            Ub = sur._ub(X)[:, free]
            t = np.minimum(Ub, 1.0 - Ub) / edge
            J += w_edge * np.sum(np.clip(1.0 - t, 0.0, None) ** 2, axis=1)
        r = sur.nn_ratio(X)
        J += w_gap * np.clip(r - gap_ratio, 0.0, None) ** 2
        return J

    def fobj(x):
        x = np.asarray(x, float)
        if x.ndim == 1:
            return float(penal(assemble(x[None, :]))[0])
        return penal(assemble(x.T))

    bounds = [(sur.box_lo[j], sur.box_hi[j]) for j in free]
    t0 = time.time()
    r = differential_evolution(fobj, bounds, maxiter=maxiter, popsize=popsize, tol=1e-8, polish=polish,
                               init="sobol", vectorized=True, updating="deferred",
                               **_rng_kw(differential_evolution, seed))
    x = assemble(np.asarray(r.x)[None, :])
    tol = _tol(sur, objective, tol)
    # упирается ли оптимум в край: без штрафа за край стало бы лучше, если сдвинуть точку на край оси
    J0, at_edge = float(penal(x, False)[0]), []
    for j in free:
        u = sur._ub(x)[0, j]
        xe = x.copy()
        xe[0, j] = sur.box_lo[j] if u < 0.5 else sur.box_hi[j]
        if float(penal(xe, False)[0]) < J0 - 0.1 * tol / s_obj:
            at_edge.append(f"{sur.space.names[j]}:{'низ' if u < 0.5 else 'верх'}")
    return _report_point(sur, x, objective, sign, cons, kind, kappa, edge, gap_ratio, tol, free, at_edge,
                         de=dict(success=bool(r.success), message=str(r.message), nfev=int(r.nfev),
                                 nit=int(r.nit), fun=float(r.fun), seconds=round(time.time() - t0, 2)))


def _tol(sur, objective, tol):
    """допуск цели: заданный; иначе 2 RMSE перекрёстной проверки (если модель обучалась через CLI);
    иначе 0,1 СКО цели по данным"""
    if tol is None:
        tol = (sur.report.get("tol") or {}).get(objective)
    if tol is None:
        tol = 0.1 * float(sur.y_scale[sur.outputs.index(objective)])
    return float(tol)


def _report_point(sur, x, objective, sign, cons, kind, kappa, edge, gap_ratio, tol, free, at_edge, de=None):
    pred = {}
    for k in sur.kinds:
        mu, sd = sur.predict(x, kind=k)
        pred[k] = {o: dict(mean=float(mu[0, j]), sd=float(sd[0, j])) for j, o in enumerate(sur.outputs)}
    flags, notes = [], []
    s_main = pred[kind][objective]["sd"]
    if s_main > tol:
        flags.append("sigma")
        notes.append(f"СКО модели {kind} по {objective} = {s_main:.4g} больше допуска {tol:.4g}")
    others = [k for k in sur.kinds if k != kind]
    for k in others:
        dmu = abs(pred[k][objective]["mean"] - pred[kind][objective]["mean"])
        if dmu > tol:
            flags.append("disagree")
            notes.append(f"{kind} и {k} расходятся по {objective} на {dmu:.4g} (> {tol:.4g})")
    nn = float(sur.nn_ratio(x)[0])
    if nn > gap_ratio:
        flags.append("gap")
        notes.append(f"ближайшая обучающая точка в {nn:.2f} типичных расстояниях (> {gap_ratio})")
    _, viol, _, _ = _evaluate(sur, x, kind, objective, sign, cons, 0.0)
    if viol[0] > 1e-9:
        flags.append("constraint")
        notes.append("ограничения не выполнены даже по среднему прогнозу")
    Ub = sur._ub(x)[0]
    if at_edge:
        notes.append("оптимум упирается в край изученной области (" + ", ".join(at_edge)
                     + "): без штрафа за край модель ушла бы дальше; истинный оптимум может лежать "
                       "снаружи — расширить план в эту сторону")
    return dict(x=dict(zip(sur.space.names, x[0].tolist())), objective=objective,
                maximize=sign > 0, kind=kind, kappa=kappa, tol=tol, constraints=[list(c) for c in cons],
                pred=pred, nn_ratio=nn, edge_unit=dict(zip(sur.space.names, Ub.tolist())),
                at_edge=at_edge, uncertain=bool(flags), flags=flags, notes=notes, de=de)


# ======================================================================= следующая порция точек
def expected_improvement(mu, sd, best, xi=0.0, maximize=True):
    """ожидаемое улучшение относительно best (для максимума или минимума)"""
    mu, sd = np.asarray(mu, float), np.asarray(sd, float)
    imp = (mu - best - xi) if maximize else (best - mu - xi)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(sd > 0, imp / sd, 0.0)
        ei = np.where(sd > 0, imp * norm.cdf(z) + sd * norm.pdf(z), np.clip(imp, 0.0, None))
    return np.clip(ei, 0.0, None)


def _pof(sur, X, mu, sd, cons):
    p = np.ones(len(X))
    for name, op, c in cons:
        if name in sur.outputs:
            j = sur.outputs.index(name)
            s = np.maximum(sd[:, j], 1e-12)
            z = (mu[:, j] - c) / s if op == ">=" else (c - mu[:, j]) / s
            p *= norm.cdf(z)
        else:
            j = sur.space.names.index(name)
            p *= (X[:, j] >= c) if op == ">=" else (X[:, j] <= c)
    return p


def propose(sur, k, objective, maximize=True, constraints=(), fixed=None, kind="gp", domain="bounds",
            n_cand=4096, local_frac=0.25, local_sd=0.05, xi=None, exclude=0.5, beta=2.0, seed=0):
    """k новых точек по ожидаемому улучшению, умноженному на вероятность выполнить ограничения.

    Кандидаты: точки Соболя в номинальных границах (domain='bounds') или в изученной области
    ('studied') + облако вокруг текущего лучшего прогноза (local_frac кандидатов, СКО local_sd в
    долях диапазона). Лучшее значение — лучший средний прогноз среди кандидатов, выполняющих
    ограничения по среднему (подходит для зашумлённых расчётов). Порция: kind='gp' — Kriging
    believer (точка добавляется в GP как точная со своим прогнозом, гиперпараметры прежние, EI
    пересчитывается); kind='mlp' — разброс ансамбля как СКО. При обоих видах кандидаты ближе
    exclude*nn_med (типичных расстояний между соседями) к уже выбранным точкам исключаются.
    Когда EI исчезает везде (максимум меньше 1e-6 СКО цели — уточнённый believer-точками район
    модели уже нечего улучшать), точка берётся по верхней доверительной границе
    sign*mu + beta*sigma среди кандидатов с вероятностью выполнить ограничения >= 0,5
    (режим 'ucb' в строке результата), а не по шуму в EI.
    В строке результата: mode — ei | ucb | pof; acq — значение критерия, по которому выбрана точка
    (EI*PoF, sign*mu + beta*sigma или PoF — единицы разные); ei — EI*PoF в единицах цели при любом
    режиме (0, если допустимых кандидатов нет): по нему и судят, исчерпан ли дорасчёт."""
    cons = _check_task(sur, objective, kind, constraints)
    if domain not in ("bounds", "studied"):
        raise ValueError(f"domain: bounds или studied, а не {domain!r}")
    sign = 1.0 if maximize else -1.0
    io = sur.outputs.index(objective)
    free, base = _full(sur, fixed) if domain == "studied" else _free_nominal(sur, fixed)
    if not free:
        raise ValueError("все входы закреплены — выбирать нечего")
    lo = sur.box_lo if domain == "studied" else sur.space.lo
    hi = sur.box_hi if domain == "studied" else sur.space.hi
    rng = np.random.default_rng(seed)
    n_loc = int(n_cand * local_frac)
    eng = qmc.Sobol(len(free), scramble=True, **_rng_kw(qmc.Sobol, seed))
    Ug = eng.random_base2(int(math.ceil(math.log2(max(2, n_cand - n_loc)))))[: n_cand - n_loc]
    Xc = np.tile(base, (len(Ug), 1))
    Xc[:, free] = lo[free] + Ug * (hi[free] - lo[free])
    if n_loc > 0:
        mu, sd = sur.predict(Xc, kind=kind)
        ok = _feasible_mean(sur, Xc, mu, cons)
        if ok.any():
            centre = Xc[np.argmax(np.where(ok, sign * mu[:, io], -np.inf))]
        else:
            centre = Xc[np.argmax(_pof(sur, Xc, mu, sd, cons))]
        Xl = np.tile(centre, (n_loc, 1))
        Xl[:, free] += rng.normal(0.0, local_sd, (n_loc, len(free))) * (hi[free] - lo[free])
        Xl[:, free] = np.clip(Xl[:, free], lo[free], hi[free])
        Xc = np.vstack([Xc, Xl])
    if xi is None:
        xi = 0.01 * float(sur.y_scale[io])
    model = sur.models[kind]
    Uc = sur.space.to_unit(Xc)
    chosen, rows = [], []
    for step in range(int(k)):
        mu, sd = sur._inverse(*model.predict(Uc))
        ok = _feasible_mean(sur, Xc, mu, cons)
        pof = _pof(sur, Xc, mu, sd, cons)
        if ok.any():
            best = float(np.max(sign * mu[ok, io])) * sign
            acq = expected_improvement(mu[:, io], sd[:, io], best, xi, maximize) * pof
        else:
            best, acq = None, pof.copy()
        eipof = acq.copy() if best is not None else np.zeros(len(Xc))   # EI x PoF, в единицах цели
        banned = np.zeros(len(Xc), bool)
        if chosen:  # не ближе exclude типичных расстояний к уже выбранным точкам порции
            d = cKDTree(sur._ub(np.array(chosen))).query(sur._ub(Xc), k=1)[0]
            banned = d < exclude * sur.nn_med
        acq[banned] = -1.0
        mode = "ei" if best is not None else "pof"
        if best is not None and acq.max() < 1e-6 * float(sur.y_scale[io]):
            mode = "ucb"
            ucb = sign * mu[:, io] + beta * sd[:, io]
            good = (pof >= 0.5) & ~banned
            acq = np.where(good if good.any() else ~banned, ucb, -np.inf)
        i = int(np.argmax(acq))
        x = Xc[i].copy()
        chosen.append(x)
        rows.append(dict(x=dict(zip(sur.space.names, x.tolist())), mode=mode, acq=float(acq[i]),
                         ei=float(eipof[i]), pof=float(pof[i]),
                         best=best, pred={o: dict(mean=float(mu[i, j]), sd=float(sd[i, j]))
                                          for j, o in enumerate(sur.outputs)}))
        if kind == "gp":  # Kriging believer: точка входит в модель со своим прогнозом
            u = sur.space.to_unit(x[None, :])
            m_u, _ = model.predict(u)
            model = model.condition(u, m_u)
    return np.array(chosen), rows


def _free_nominal(sur, fixed):
    fixed = dict(fixed or {})
    names = sur.space.names
    unknown = [k for k in fixed if k not in names]
    if unknown:
        raise ValueError(f"нет таких входов: {unknown}")
    for k, v in fixed.items():
        j = names.index(k)
        if not (sur.space.lo[j] <= v <= sur.space.hi[j]):
            raise ValueError(f"{k}={v} вне границ задачи")
    free = [j for j, k in enumerate(names) if k not in fixed]
    return free, np.array([fixed.get(k, np.nan) for k in names], float)


def _feasible_mean(sur, X, mu, cons):
    ok = np.ones(len(X), bool)
    for name, op, c in cons:
        g = mu[:, sur.outputs.index(name)] if name in sur.outputs else X[:, sur.space.names.index(name)]
        ok &= (g >= c) if op == ">=" else (g <= c)
    return ok


# ======================================================================= таблицы
def read_table(path, inputs, outputs, group_col=None):
    """CSV с заголовком -> X, Y, метки блоков (или None), число отброшенных строк.
    Строки с пустым/нечисловым показателем (неудачный прогон) отбрасываются и считаются."""
    X, Y, G, dropped = [], [], [], 0
    with open(path, newline="", encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        need = list(inputs) + list(outputs) + ([group_col] if group_col else [])
        miss = [c for c in need if c not in (rd.fieldnames or [])]
        if miss:
            raise ValueError(f"{path}: нет столбцов {miss}")
        for row in rd:
            try:
                x = [float(row[c]) for c in inputs]
                y = [float(row[c]) for c in outputs]
            except (TypeError, ValueError):
                dropped += 1
                continue
            if not (np.isfinite(x).all() and np.isfinite(y).all()):
                dropped += 1
                continue
            X.append(x)
            Y.append(y)
            G.append(row.get(group_col) if group_col else None)
    return np.array(X), np.array(Y), (np.array(G) if group_col else None), dropped


def write_table(path, X, names, extra=None):
    """X (n, d) с именами столбцов names + необязательные столбцы extra {имя: массив};
    нечисловые значения (например, режим выбора точки) пишутся как текст"""
    extra = extra or {}
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(list(names) + list(extra))
        for i, row in enumerate(np.asarray(X, float)):
            w.writerow([f"{v:.10g}" for v in row] + [_cell(extra[c][i]) for c in extra])


def _cell(v):
    if isinstance(v, str):
        return v
    return f"{float(v):.10g}"


# ======================================================================= командная строка
def _load_problem(path):
    with open(path, encoding="utf-8") as fh:
        p = json.load(fh)
    return Space(p["inputs"]), list(p["outputs"]), p


def _fmt_metrics(title, m):
    lines = [title, f"  {'показатель':12s} {'RMSE':>10s} {'макс|e|':>10s} {'MAPE,%':>8s} {'макс,%':>8s} {'R2':>7s}"]
    for o, r in m.items():
        lines.append(f"  {o:12s} {r['rmse']:10.4g} {r['max_abs']:10.4g} {r['mape_pct']:8.3f} "
                     f"{r['max_rel_pct']:8.3f} {r['r2']:7.4f}")
    return "\n".join(lines)


def _strip(obj):
    """убрать массивы предсказаний из отчёта перед записью в JSON"""
    if isinstance(obj, dict):
        return {k: _strip(v) for k, v in obj.items() if k != "pred"}
    return obj


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("plan", help="план расчётов")
    a.add_argument("--problem", required=True)
    a.add_argument("--n", type=int, required=True)
    a.add_argument("--method", default="sobol", choices=["sobol", "lhs", "random"])
    a.add_argument("--seed", type=int, default=0)
    a.add_argument("--out", required=True)

    a = sub.add_parser("fit", help="обучение и проверка")
    a.add_argument("--problem", required=True)
    a.add_argument("--data", required=True)
    a.add_argument("--holdout")
    a.add_argument("--kinds", default="mlp,gp")
    a.add_argument("--folds", type=int, default=5)
    a.add_argument("--blocks", default="kmeans", help="kmeans | order | имя столбца с метками порций")
    a.add_argument("--n-blocks", type=int, default=20)
    a.add_argument("--log", default="", help="показатели через запятую, которые учить в логарифме")
    a.add_argument("--members", type=int, default=8)
    a.add_argument("--hidden", default="32,32")
    a.add_argument("--alpha", default="auto", help="auto или число (L2-регуляризация MLP)")
    a.add_argument("--jobs", type=int, default=4)
    a.add_argument("--seed", type=int, default=0)
    a.add_argument("--out", required=True)
    a.add_argument("--report")

    for nm in ("optimize", "next"):
        a = sub.add_parser(nm, help="оптимум по модели" if nm == "optimize" else "следующая порция точек (EI)")
        a.add_argument("--model", required=True)
        a.add_argument("--objective", required=True)
        g = a.add_mutually_exclusive_group()
        g.add_argument("--max", dest="maximize", action="store_true", default=True)
        g.add_argument("--min", dest="maximize", action="store_false")
        a.add_argument("--constraint", action="append", default=[])
        a.add_argument("--fixed", action="append", default=[], help="имя=значение")
        a.add_argument("--seed", type=int, default=0)
        a.add_argument("--out")
        if nm == "optimize":
            a.add_argument("--kind", default="mlp")
            a.add_argument("--kappa", type=float, default=1.0)
            a.add_argument("--edge", type=float, default=0.05)
            a.add_argument("--tol", type=float)
        else:
            a.add_argument("--kind", default="gp")
            a.add_argument("--k", type=int, default=8)
            a.add_argument("--domain", default="bounds", choices=["bounds", "studied"])

    args = ap.parse_args(argv)

    if args.cmd == "plan":
        space, outputs, _ = _load_problem(args.problem)
        X, info = design(space, args.n, args.method, args.seed)
        write_table(args.out, X, space.names)
        print(json.dumps(info, ensure_ascii=False))
        return 0

    if args.cmd == "fit":
        space, outputs, _ = _load_problem(args.problem)
        blocks_col = args.blocks if args.blocks not in ("kmeans", "order") else None
        X, Y, G, dropped = read_table(args.data, space.names, outputs, blocks_col)
        kinds = tuple(args.kinds.split(","))
        mlp = dict(members=args.members, hidden=tuple(int(v) for v in args.hidden.split(",")),
                   alpha=args.alpha if args.alpha == "auto" else float(args.alpha))
        t0 = time.time()
        log = tuple(v for v in args.log.split(",") if v)
        cv, _ = block_cv(X, Y, space, outputs, n_folds=args.folds, groups=G if blocks_col else args.blocks,
                         n_blocks=args.n_blocks, kinds=kinds, mlp=mlp, log=log, seed=args.seed, n_jobs=args.jobs)
        sur = Surrogate(space, outputs, kinds=kinds, mlp=mlp, log=log, seed=args.seed, n_jobs=args.jobs).fit(X, Y)
        rep = dict(data=args.data, n=len(X), dropped=dropped, cv=_strip(cv),
                   tol={o: 2.0 * cv[kinds[0]]["metrics"][o]["rmse"] for o in outputs})
        print(_fmt_metrics(f"перекрёстная по блокам ({args.blocks}, {args.folds} частей), {kinds[0]}:",
                           cv[kinds[0]]["metrics"]))
        if args.holdout:
            Xh, Yh, _, dh = read_table(args.holdout, space.names, outputs)
            rep["holdout"] = dict(data=args.holdout, n=len(Xh), dropped=dh, metrics=validate(sur, Xh, Yh))
            for k in kinds:
                print(_fmt_metrics(f"отложенная выборка, {k}:", rep["holdout"]["metrics"][k]))
        rep["model"] = sur.describe()
        rep["seconds"] = round(time.time() - t0, 1)
        sur.report = rep
        sur.save(args.out)
        if args.report:
            with open(args.report, "w", encoding="utf-8") as fh:
                json.dump(rep, fh, ensure_ascii=False, indent=1)
        print(f"модель: {args.out}; строк {len(X)}, отброшено {dropped}; {rep['seconds']} с")
        return 0

    sur = Surrogate.load(args.model)
    fixed = {}
    for s in args.fixed:
        k, v = s.split("=", 1)
        fixed[k.strip()] = float(v)
    if args.cmd == "optimize":
        r = optimize(sur, args.objective, args.maximize, args.constraint, fixed, kind=args.kind,
                     kappa=args.kappa, edge=args.edge, tol=args.tol, seed=args.seed)
        txt = json.dumps(r, ensure_ascii=False, indent=1)
        print(txt)
        if args.out:
            with open(args.out, "w", encoding="utf-8") as fh:
                fh.write(txt)
        return 0
    Xn, rows = propose(sur, args.k, args.objective, args.maximize, args.constraint, fixed,
                       kind=args.kind, domain=args.domain, seed=args.seed)
    if args.out:
        io = args.objective
        write_table(args.out, Xn, sur.space.names,
                    {"mode": [r["mode"] for r in rows], "acq": [r["acq"] for r in rows],
                     "ei": [r["ei"] for r in rows], "pof": [r["pof"] for r in rows],
                     f"pred_{io}": [r["pred"][io]["mean"] for r in rows],
                     f"sd_{io}": [r["pred"][io]["sd"] for r in rows]})
    print(json.dumps(rows, ensure_ascii=False, indent=1))
    return 0


if __name__ == "__main__":
    # Запуск скриптом: классы иначе попали бы в model.joblib как __main__.Surrogate, и модель,
    # сохранённая `python surrogate.py fit`, не грузилась бы через `import surrogate`.
    # Поэтому работаем через модуль, импортированный под своим именем.
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import surrogate as _mod
    sys.exit(_mod.main())
