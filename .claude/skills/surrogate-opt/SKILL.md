---
name: surrogate-opt
description: Модель-заменитель (суррогат) по выборке дорогих расчётов SU2, OpenFOAM и CalculiX трека «Простая литая ГТУ-заменитель» и оптимизация по ней скриптом routes/gtu/platform/surrogate/surrogate.py — план расчётов (Соболь, латинский гиперкуб), нормировка, ансамбль MLP и гауссовский процесс, перекрёстная проверка по блокам и отложенная выборка, оценка неопределённости, активный дорасчёт по ожидаемому улучшению, оптимум только внутри изученной области с проверкой прямым расчётом, хранение выборки и модели под ключом прогона. Используй, когда планируешь серию дорогих расчётов, строишь модель по их итогам, ищешь по ней оптимум колеса, камеры или диска, выбираешь следующую порцию точек или проверяешь найденный оптимум. Триггеры — «суррогат», «модель-заменитель», «нейросеть по расчётам», «план расчётов», «Соболь», «латинский гиперкуб», «гауссовский процесс», «ожидаемое улучшение», «активный дорасчёт», «оптимизация по модели», «surrogate.py».
---

# Модель-заменитель и оптимизация по выборке дорогих расчётов

Узел платформы трека (issue #59): дорогой расчёт (SU2 — ступени, OpenFOAM — горение и теплообмен, CalculiX —
прочность) даёт таблицу «параметры → показатели»; модель по ней ищет оптимум и выбирает, что считать дальше.
На 14.11.2026 `surrogate.py` проверен только на заменителе `gtu_explore.comp_stage` (3 входа); реальных
выборок SU2 / OpenFOAM / CalculiX нет. Пометки: «(код)» — по исходнику, «(док.)» — по документации, «(опыт)» —
проверено запуском 14.11.2026, «(поиск)» — по выдержке поиска, полный текст не открыт.

## 1. Что где лежит

| что | где |
|---|---|
| план, обучение, проверка, оптимум, дорасчёт | `routes/gtu/platform/surrogate/surrogate.py`; описание и итоги — `README.md` там же |
| тест на заменителе и его итоги | `test_surrogate.py`, `test_surrogate_results.json` (запись — `SURROGATE_RECORD=1`) |
| дешёвый прямой расчёт для обкатки | `routes/gtu/gtu_explore.py` → `comp_stage` (навык `gtu-cycle`) |
| дорогие расчёты точек | SU2 — навык `su2-turbo`; OpenFOAM — `openfoam-combustor`; CalculiX — `routes/gtu/platform/fea/disk.py` |
| облачные прогоны, ключ прогона | `routes/gtu/platform/aws/runner.py` (`new_run_id`, `key_prefix`), README там же |
| выборки и модели прогонов | `routes/gtu/platform/surrogate/runs/<run_id>/` — папки ещё нет, договорённость этого навыка |

## 2. Установка и ресурсы

```bash
pip install scikit-learn scipy numpy      # joblib и threadpoolctl — зависимости scikit-learn
pip install cantera pytest                # только для test_surrogate.py
```

Стоят (опыт): Python 3.13.18, scikit-learn 1.9.1 (= последняя на PyPI 14.11.2026; зависимости по PyPI — numpy,
scipy, joblib, narwhals, threadpoolctl), scipy 1.18.1, numpy 2.5.3, joblib 1.6.0. `fit` на 32 точках — 16 с
(опыт), на 128 — 73 с (README); на 4 ядрах `--jobs 4`, BLAS в один поток (`SURROGATE_BLAS_THREADS=1`, код).
Местный расчёт — не дольше 20 мин; дорогие — на AWS через `runner.py`; здесь ключи AWS не действуют.

## 3. Порядок работы

**Шаг 0. Начало сессии** — по CLAUDE.md: `.agent/manifest.json`, `data/catalog.json`, `git ls-remote --heads origin`.

**Шаг 1. Постановка — `problem.json`.**

```json
{"inputs": {"U2": [350, 650], "m": [1.5, 25], "T1": [243.15, 318.15]}, "outputs": ["eta", "PR", "D2"]}
```

- Границы — физически допустимые и считаемые решателем: титан до 620 м/с, ограничения литья — `LIM` в
  `gtu_explore.py` (навык `gtu-cycle`). Где решатель не сходится — пустая строка в данных, а не граница на глаз.
- Входов — минимум: объём выборки растёт с размерностью. Дискретные входы (число лопаток) модуль считает
  непрерывными (код): найденное значение округлить и пересчитать прямо.
- Цель и ограничения назвать до плана: `eta` — цель, `PR>=4.25`, `D2<=0.40` — ограничения.

**Шаг 2. План и отложенная выборка.**

```bash
python surrogate.py plan --problem problem.json --n 64 --method sobol --seed 1 --out plan.csv
python surrogate.py plan --problem problem.json --n 16 --method random --seed 7 --out hold_in.csv
```

- Соболь (`qmc.Sobol`, скремблирование) — только n = 2^k: «lose their balance properties if one uses a sample
  size that is not a power of 2» (док. scipy). При другом n модуль берёт первые n из 2^k и печатает
  предупреждение (опыт). Латинский гиперкуб — `--method lhs` = `qmc.LatinHypercube(optimization="random-cd")`,
  n любое (код).
- Начальный объём — предварительно n ≈ 10·d (Loeppky, Sacks, Welch 2009, поиск); хватает ли — решает
  проверка по блокам (шаг 4), а не правило.
- Отложенная выборка — отдельный план с другим зерном (`hold_in.csv` → тот же решатель, те же настройки →
  `hold.csv`); в обучении и в выборе настроек модели не участвует.
- `discrepancy_cd` в выводе — центрированная L2-неравномерность (`qmc.discrepancy`, `method="CD"`, код).

**Шаг 3. Дорогие расчёты под ключом прогона.**

- Ключ `run_id` — формат `runner.new_run_id()`: `YYYYMMDDTHHMMSSZ-xxxx` (UTC и 4 шестнадцатеричных знака, код).
  В облаке результаты задачи — `s3://<бакет>/<префикс>/gtu/<run_id>/<job>/` (`key_prefix`, код); бакет и
  префикс — из `ALPHA_O_S3`, запись только под ним.
- `python3 runner.py run jobs.json --branch B` — проверка прав (DryRun); `--live` — живой спот в us-east-1
  с потолком `--ceiling 300` (код). Живой запуск — только с разрешения владельца.
- Каждая точка — строка `runs.csv`: входы, показатели, `batch` (номер порции), `run_id`, `job`, решатель и
  версия, уровень сетки; текстовые столбцы сверх входов и показателей не мешают (опыт). Не сошлось —
  показатели пустые: `read_table` отбрасывает строку и пишет число в `dropped` (код, опыт).

**Шаг 4. Обучение и проверка.**

```bash
python surrogate.py fit --problem problem.json --data runs.csv --holdout hold.csv --blocks batch \
       --log D2 --kinds mlp,gp --jobs 4 --out model.joblib --report report.json
```

- Нормировка (код): входы → [0, 1] по границам `problem.json`; показатели стандартизуются по обучающей части
  каждой складки; `--log` — положительные показатели с разбросом в разы (диаметр, расход, напряжение).
- MLP — ансамбль `--members 8` сетей `MLPRegressor` (tanh, `--hidden 32,32`, `solver="lbfgs"`), у каждой своё
  зерно и бутстреп; `--alpha auto` — L2 по внутренней 4-кратной проверке (код). lbfgs: «For small datasets,
  however, 'lbfgs' can converge faster and perform better»; ансамбль: у MLP невыпуклые потери, «different random
  weight initializations can lead to different validation accuracy» (док. scikit-learn). СКО — разброс членов.
- GP — `GaussianProcessRegressor`, ядро `ConstantKernel × Matern(nu=2.5)` со своей длиной по каждой оси
  `+ WhiteKernel` (оценка шума), 3 перезапуска оптимизатора (код; `n_restarts_optimizer`, `WhiteKernel` — док.).
- Блоки: `--blocks batch` (столбец порций) — после дорасчёта; `kmeans` (`--n-blocks 20` пространственных
  кластеров) — для одного плана; `order` — подряд идущие порции. Блоки делятся на `--folds 5` через
  `GroupKFold`: группа не бывает сразу в обучении и в проверке, групп не меньше, чем складок (док.).
  Случайное деление структурированных данных занижает ошибку (Roberts et al. 2017, поиск).
- Приёмка по `report.json` → `cv` и `holdout`, по каждому показателю: `rmse` и `max_abs` в единицах показателя
  ниже точности, нужной решению (для КПД — доли п. п.); `cover_2sd` (доля точек с |ошибкой| ≤ 2 СКО) близка к
  0,95 при честной оценке неопределённости; MLP и GP согласны. Не прошла — больше точек (шаг 6,
  `--domain bounds`) или меньше входов, а не подгонка настроек под отложенную выборку.

**Шаг 5. Оптимум — только внутри изученной области.**

```bash
python surrogate.py optimize --model model.joblib --objective eta --max --kind mlp \
       --constraint "PR>=4.25" --constraint "D2<=0.40" --fixed T1=288.15 --out opt_mlp.json
```

- Повторить с `--kind gp --out opt_gp.json`: две модели — две независимые оценки одной точки.
- `scipy.optimize.differential_evolution` ищет только в прямоугольнике обучающих точек (`studied` в
  `report.json` → `model`), по осторожной оценке μ − κσ (`--kappa 1`, ограничения тоже осторожные), со штрафом
  в полосе `--edge 0.05` у края и за удаление от данных больше двух типичных расстояний (код). `--fixed` вне
  изученной области — `ValueError … вне изученной области` (опыт).
- Флаги `opt_*.json`: `sigma`, `disagree` (MLP и GP разошлись больше `tol`), `gap`, `constraint`; отдельно
  `at_edge` — оптимум упёрся в край, истинный может лежать снаружи (код). `tol` по умолчанию — 2 RMSE
  перекрёстной проверки (код); лучше задать `--tol` из требования решения.
- Есть флаг или `at_edge` — оптимум не публикуется: шаг 6 (расширить план) и повтор.

**Шаг 6. Активный дорасчёт.**

```bash
python surrogate.py next --model model.joblib --objective eta --max --k 6 --kind gp --domain bounds \
       --constraint "PR>=4.25" --fixed T1=288.15 --out next.csv
```

- Критерий — ожидаемое улучшение EI = (μ − f*)·Φ(z) + σ·φ(z), z = (μ − f*)/σ (Jones, Schonlau, Welch 1998 —
  поиск; формула — код `expected_improvement`), умноженное на вероятность выполнить ограничения. f* — лучшее
  среднее прогноза среди допустимых кандидатов, а не лучший зашумлённый расчёт (код).
- Порция из k точек — «Kriging believer» (Ginsbourger, Le Riche, Carraro 2010, поиск): выбранная точка входит в
  GP с прогнозом вместо ответа, EI пересчитывается. EI исчез везде — точки по μ + 2σ (`ucb` в выводе, код).
- `k` — по числу параллельных машин (`runner.py --max-parallel`, по умолчанию 6, код). `next.csv` → расчёты →
  строки в `runs.csv` с новым `batch` → шаг 4 с `--blocks batch`.
- Стоп: максимум `acq` в `next.csv` меньше допуска цели две порции подряд или исчерпан бюджет (договорённость
  навыка, на реальных данных не проверено).

**Шаг 7. Проверка оптимума прямым расчётом.** В точке `x` из `opt_*.json` — тот же решатель, та же сетка, те же
настройки. Принято, если |прогноз − прямой| ≤ `tol` и прямой расчёт выполняет ограничения. Нет — точку в
`runs.csv`, шаги 4–5 заново. Для SU2 — ещё сетка coarse → medium (навык `su2-turbo`). Смок 14.11.2026 (опыт):
прогноз 0,04657 ± 0,00098, прямой 0,04563, ограничение g ≤ 1,4 выполнено (1,3994).

**Шаг 8. Хранение под ключом прогона.** `runs/<run_id>/`: `problem.json`, `plan.csv`, `hold.csv`, `runs.csv`,
`report.json`, `opt_*.json`, `next.csv`, `direct.json` (прямые расчёты оптимума), `manifest.json` — `run_id`,
коммит (`git rev-parse HEAD`), команды и зёрна, решатель, версия, сетка, sha256 `runs.csv`, версии библиотек из
`report.json` → `model.versions`. Это набор, который scikit-learn рекомендует хранить рядом с моделью: данные,
код, версии, оценка перекрёстной проверки (док.). `model.joblib` — в S3 под тем же ключом
(`<префикс>/gtu/<run_id>/surrogate/`); в репозиторий по умолчанию не кладётся — модель воспроизводится `fit`
из `manifest.json`. Число в дневник — из `direct.json` и results-JSON; каталог `node tools/build_catalog.mjs` —
в чистой копии HEAD (`git worktree add --detach`), вторым коммитом (CLAUDE.md, урок 1).

## 4. Обязательные правила

1. Опубликованное число — из прямого расчёта; прогноз модели — только с пометкой «прогноз модели» и СКО.
2. Отложенная выборка не участвует ни в обучении, ни в выборе настроек; использована для настройки — новая.
3. Оптимум — только внутри изученной области. Флаг не убирают ростом `--tol` или `--kappa 0`.
4. Неудачный расчёт — пустой показатель, не ноль и не подстановка; число отброшенных — в отчёт.
5. Одна выборка — один решатель, одна версия, один уровень сетки. Иначе — столбец признака и отдельные блоки.
6. Ключ прогона — в каждой строке и каждом файле; файлы прогона не переписываются, порция дописывается.
7. `model.joblib` грузится той же версией scikit-learn и только свой: загрузка модели другой версии не
   поддерживается, pickle из чужого источника — выполнение чужого кода (док.).
8. AWS: живой запуск — с разрешения владельца; us-east-1, спот, потолок 300 долл., запись только под
   `ALPHA_O_S3` в `gtu/`; ключей и секретов в файлах нет.

## 5. Типовые ошибки

- Соболь с n не степенью двойки — равномерность потеряна (док.); нужно 2^k или `lhs`.
- `qmc.discrepancy` по точкам в физических единицах — `ValueError: Sample is not in unit hypercube` (код scipy).
- Случайная k-кратная проверка по точкам порций дорасчёта: соседи порции и в обучении, и в проверке — ошибка
  занижена. Нужно `--blocks batch`.
- Нормировка по всем данным вместе с проверочными — утечка (док. scikit-learn, «Data leakage»).
- Относительная ошибка у показателя около нуля: в смоке f = −0,87…0,04, MAPE 22 % при R² 0,96 (опыт). Мерить
  `rmse` и `max_abs` в единицах показателя.
- Допуск по умолчанию от грубой модели: в смоке (32 точки, 4 порции) `tol` = 0,10 = 2 × RMSE 0,050, то есть
  11 % размаха f — `disagree` молчит, пока MLP и GP расходятся меньше 0,10 (опыт). Задавать `--tol` явно.
- СКО GP с `WhiteKernel` в scikit-learn включает шум: `predict(return_std=True)` дал 0,099 при оценке шума 0,093
  (опыт, 1.9.1). Такая σ в EI завышает разведку; `GPModel` модуля отделяет СКО функции от шума (код).
- Одна сеть MLP вместо ансамбля — итог зависит от начальных весов (док.). Рост разброса ансамбля вне данных
  показан на задачах статьи (Lakshminarayanan et al. 2017, поиск), не гарантирован — вне `studied` не искать.
- Оптимум с `at_edge` объявлен оптимумом: на заменителе при данных U2 ≥ 521 м/с найдено 527,3, прямой КПД
  0,8544 против 0,8586 у истинного оптимума при U2 = 505,5 (README).
- GP на тысячах точек — O(n³) по времени, O(n²) по памяти; больше 208 точек не проверено (README).

## 6. Проверка результата

1. `python -m pytest -q -s test_surrogate.py` — 8 тестов пройдены 14.11.2026 за 203 с на 4 ядрах (опыт).
2. `report.json`: `cv` и `holdout` — `rmse` ниже нужной точности, `cover_2sd` около 0,95, MLP и GP согласны.
3. `opt_mlp.json` и `opt_gp.json`: флагов нет, `at_edge` пуст, цель в найденных точках совпадает в пределах `tol`.
4. Прямой расчёт в оптимуме: |прогноз − прямой| ≤ `tol`, ограничения выполнены, сетка уточнена.
5. `manifest.json` полон; повторный `fit` из него даёт те же метрики при тех же зёрнах.
6. `node tools/validate.mjs`, `node tools/build_catalog.mjs --check` (в чистой копии HEAD), `node tools/smoke.mjs`.

## 7. Скрипты трека

- `routes/gtu/platform/surrogate/` (`surrogate.py`, тест, итоги, README); `platform/aws/runner.py` — прогоны и ключ.
- `routes/gtu/gtu_explore.py` (`comp_stage`; схема 130 в `gtu_explore_results.json`: колесо 0,435 м, PR 4,25),
  `gtu_dropin.py` (Taurus 60, одно колесо, PR 6,75), `platform/fea/disk.py`, `platform/geom/impeller.py`.

## Источники (сверено 14.11.2026)

- scipy v1.18.1, исходник с docstring (сайт docs.scipy.org из контейнера не открылся — DNS):
  https://raw.githubusercontent.com/scipy/scipy/v1.18.1/scipy/stats/_qmc.py (`Sobol`, `random_base2`,
  `LatinHypercube`, `discrepancy`), …/v1.18.1/scipy/optimize/_differentialevolution.py ; страница:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.qmc.Sobol.html
- scikit-learn 1.9.1, исходник документации (сайт scikit-learn.org не открылся — DNS):
  https://raw.githubusercontent.com/scikit-learn/scikit-learn/1.9.1/doc/modules/neural_networks_supervised.rst ,
  …/1.9.1/doc/modules/gaussian_process.rst , …/1.9.1/doc/modules/cross_validation.rst ,
  …/1.9.1/doc/common_pitfalls.rst , …/1.9.1/doc/model_persistence.rst , …/1.9.1/sklearn/model_selection/_split.py ,
  …/1.9.1/sklearn/gaussian_process/_gpr.py , …/1.9.1/sklearn/neural_network/_multilayer_perceptron.py ;
  страницы: https://scikit-learn.org/stable/modules/neural_networks_supervised.html ,
  https://scikit-learn.org/stable/model_persistence.html
- PyPI (версии): https://pypi.org/pypi/scikit-learn/json , https://pypi.org/pypi/scipy/json , https://pypi.org/pypi/numpy/json
- Только выдержки поиска: Jones, Schonlau, Welch 1998, EGO — https://doi.org/10.1023/A:1008306431147 ;
  Ginsbourger, Le Riche, Carraro 2010 — https://hal-emse.ccsd.cnrs.fr/emse-00436126 ; Loeppky, Sacks, Welch 2009 —
  https://doi.org/10.1198/TECH.2009.08040 ; Lakshminarayanan et al. 2017 — https://arxiv.org/abs/1612.01474 ;
  Roberts et al. 2017 — https://doi.org/10.1111/ecog.02881
- Опыт: все четыре команды `surrogate.py` на аналитической функции (32 точки Соболя, 16 отложенных), проверка
  σ GP с `WhiteKernel`, `test_surrogate.py` — в контейнере 14.11.2026.
