---
name: gtu-cycle
description: Модель цикла трека «Простая литая ГТУ-заменитель» — routes/gtu/gtu_cycle.py (cycle, cycle_r, recup на Cantera) и перебор схем routes/gtu/gtu_explore.py (comp_stage, evaluate, finalize, optimize, concepts) — расчёт точки и схемы, размер центробежного колеса, калибровка по SGT-400 и Taurus 60, сверка с Mercury 50 и OPRA OP16, принятые ограничения литья и окружных скоростей, погрешность ±2 п. п., новая схема в перебор без порчи сверки; каждое опубликованное число — из results-JSON. Используй, когда считаешь цикл или схему двигателя-заменителя, правишь модель или перебор, обновляешь results-JSON трека или цитируешь его числа. Триггеры — «цикл ГТУ», «gtu_cycle», «gtu_explore», «перебор схем», «новая схема», «калибровка модели», «сверка с паспортом», «размер колеса», «КПД двигателя-заменителя».
---

# Модель цикла ГТУ-заменителя: расчёт, калибровка, перебор

Модель — уровень облика (0D-цикл и 1D-размеры), не расчёт течения. Её числа отбирают схемы и дают
вход для SU2, OpenFOAM и CalculiX. Навык — как считать на ней, не сломав сверку, и как публиковать
числа так, чтобы любой повторил их по репозиторию.

## 1. Что где лежит (`routes/gtu/`)

| что | скрипт → результат | функции и ключи |
|---|---|---|
| цикл на Cantera (GRI-Mech 3.0, метан, ISO) | `gtu_cycle.py` → `gtu_cycle_results.json` | `cycle`, `cycle_r`, `recup` |
| калибровка КПД узлов | `gtu_study.py` → `gtu_study_results.json` | `calib`; ключи `SGT400`, `T60`, `grid` |
| точки базового осевого варианта | `gtu_point.py` (читает `gtu_study_results.json`) | `A:…`, `B:…` |
| простые схемы, прецедент OP16 | `gtu_simple.py` → `gtu_simple_results.json` | `OP16-подобный …`, `простой …` |
| перебор 249 схем, размеры колёс | `gtu_explore.py` → `gtu_explore_results.json` | `lim`, `res[]` (поле `id`) |
| замена в пакетах Taurus 60 и SGT-400 | `gtu_dropin.py` → `gtu_dropin_results.json` | `ARCH`, `COOL`, `SLOT`, `res[]` |
| семейство по паспортам | `gtu_family.py` | паспорта `fleet/` — навык `gtu-fleet` |

Цель, сверка, принятые величины — `README.md`; платформа — `platform/` (`aws`, `fea`, `geom`, `surrogate`).

## 2. Как устроена модель

**`cycle(PR, T41_C, eta_pc, eta_pt_gg, eta_pt_pt, x_ngv, x_rot, x_leak, dp_comb=0.04, dp_duct=0.01,
p_exh_ratio=1.02, eta_b=0.995, eta_m=0.99, W1=1.0)`** — двухвальная ГТУ со свободной силовой
турбиной. Сжатие и расширение — 40 политропных шагов при замороженном составе; горение полное
(CH₄ + 2O₂ → CO₂ + 2H₂O), недожог `eta_b` — потеря тепла. `T41_C` — газ перед рабочим колесом 1-й
ступени после подмешивания охлаждения соплового аппарата, а не T4 камеры и не температура из
паспорта. Доли от `W1`: `x_ngv` возвращается до ротора 1, `x_rot` — после ТВД, `x_leak` теряется.
Возврат: `P_MW`, `eff` (на валу, по низшей теплоте 50,03 МДж/кг), `sp_work_kJkg`, `EGT_C`, `W_exh`,
`T3_C`, `T4_C`, `T45_C`, `GG_ER`, `PT_ER`, `wc_kJkg`, `Wf`. Смешение — `ct.Quantity(g, constant="HP")`
(по умолчанию Cantera смешивает при постоянных UV); корни — `scipy.optimize.brentq`.
`cycle_r(…, T3r=None, dp_rec_air=0.0)` — копия `cycle` с подогретым воздухом перед камерой;
`recup(PR, T41_C, eps, dp_air=0.025, dp_gas=0.035, **kw)` — рекуператор: 6 итераций по выхлопу,
добавляет `T3r_C` и `Tstack_C` (оценка через 0,93). На клеммы — `eff * 0.985 * 0.972` (редуктор и
генератор; `gtu_simple.py`, `gtu_point.py`, `ETA_GEN` в `gtu_family.py`; оценка).

**`gtu_explore.py`** — своя быстрая модель: свойства из таблиц Cantera с шагом 2 К, смесь по доле
продуктов, утечка 4 %, мощность 6 МВт.
- `comp_stage(T1, p1, PRs, m, kind, dbl=False)`, `kind` ∈ `ti`/`al`/`mono`: КПД 0,880 −
  0,010·(πₛ − 2,5) (−0,02 для монотора) − 0,3·зазор/b₂; U₂ = √(w/0,70); D₂ из m = ρ₁·U₂·D₂²·0,085;
  n = 60·U₂/(π·D₂); `ok` — U₂ в пределе, для алюминия ещё T₂ ≤ 473 К; `dbl` — двусторонний вход.
- `evaluate(c, PR, TIT_C, P=6.0e6)` — цикл схемы `c`; `None` при нарушении предела.
- `finalize(c, r, P=6.0e6)` — размеры, «литьё» (`да`/`на пределе`/отказ), индекс сложности `K`,
  `eff_el`; `None` — отливка больше предела. Аргумент `P` внутри не используется.
- `optimize(c)` — TIT от 800 °C до предела шагом 25, π 2,5…16 шагом 0,25, у монотора `kD`
  1,0/1,1/1,2; лучшая — по наибольшему `eff`, не по простоте; `evaluate` зовёт без `P`.
- `concepts()` — список схем, семейства A, B, C, E; `id` = позиция в списке + 1.

## 3. Порядок работы

**Шаг 0. Среда.** `pip install cantera scipy numpy`. Сверено 14.11.2026: cantera 3.2.0 (PyPI:
Python ≥ 3.10 и < 3.15), scipy 1.18.1, numpy 2.5.3, Python 3.13; `gri30.yaml` входит в пакет Cantera.

**Шаг 1. Копия.** Скрипт делает `chdir` в свою папку и перезаписывает свой results-JSON (с новой
`date`); импорт `gtu_explore` тоже меняет текущую папку. Пробы — только в копии:

```bash
T=$(mktemp -d); cp -r routes/gtu "$T"/ && cd "$T"/gtu && rm -rf platform __pycache__
python3 -c "
from gtu_cycle import cycle, recup
import gtu_explore as E
r = cycle(6.75, 1050, eta_pc=0.855, eta_pt_gg=0.865, eta_pt_pt=0.885, x_ngv=0.03, x_rot=0.0, x_leak=0.015)
print(round(r['eff']*100, 1), round(r['sp_work_kJkg']), round(r['EGT_C']))
st = E.comp_stage(E.T0, E.P0, 6.75, 22.5, 'ti'); print(round(st['D2'], 2), round(st['rpm']), round(st['U2']), st['ok'])
c = E.concepts()[129]; print(c['name'], E.optimize(c))   # схема id 130
"
```

Ожидается (14.11.2026, 2 с): `27.1 256 627` — как `gtu_dropin_results.json → res[arch=1Ц, cool=N]`;
`0.59 19951 618 True` — в JSON 19 967 об/мин, там расход не округлён до 22,5; схема 130 — как
`gtu_explore_results.json → res[id=130]`. Разница в последнем знаке — довод брать число из JSON.

**Шаг 2. Длительность.** Замер 14.11.2026, 4 ядра: `gtu_cycle.py` 12 с, `gtu_study.py` 99 с,
`gtu_point.py` 17 с, `gtu_simple.py` 50 с, `gtu_dropin.py` 312 с, `gtu_explore.py` 79 с;
все — в одном процессе (параллелит только `gtu_family.py`, 4 процесса). Больше 20 минут — делить.

**Шаг 3. Сверка** — разделы 4 и 9: после правки модели пересчитать калибровку и прецеденты.

**Шаг 4. Публикация.** Скрипт запускается на месте; скрипт и results-JSON — первым коммитом;
каталог (`node tools/build_catalog.mjs`, трек `gas-turbine` по префиксу `gtu_`) — в чистой копии
HEAD (`git worktree add --detach`) вторым (`.agent/LESSONS.md`, п. 1); затем проверки CLAUDE.md.
На вкладки сайта трек пока не выносится.

## 4. Калибровка и сверка

| объект | вход в коде | модель (results-JSON) | замечание |
|---|---|---|---|
| SGT-400 13 МВт, привод | 13,4 МВт; 36,2 %; π 16,8; выхлоп 39,4 кг/с; 555 °C | `gtu_study_results.json → SGT400`: КПД узлов 0,888/0,868/0,888, T41 1219 °C, выхлоп 549 °C | страница Siemens 14.11.2026, MD 13 МВт: 13,4 МВт; 36,4 %; 9886 кДж/(кВт·ч); π 17,1; 39,1 кг/с; 543 °C. 555 °C — колонка PG |
| Taurus 60, привод | 5,74 МВт; 31,96 %; π 12,2; 21,63 кг/с; 510 °C | `→ T60`: 0,874/0,864/0,884, T41 1082 °C, выхлоп 524 °C | Solar, предварительно (выдержка поиска): 5740 кВт; 11 265 кДж/(кВт·ч); 77 880 кг/ч; 510 °C. π 12,2 — только в `fleet/*_prelim.json` |
| OPRA OP16, радиальная | π 6,7 и T41 приняты, 8,8 кг/с; паспорт 1,876 МВт эл., 25,1 % | `gtu_simple_results.json → OP16-подобный одновальный T41=1000`: `eff_el` 0,237, `P_el` 1,795 МВт | T41 950/1050 °C: 23,2/24,1 %. Паспорт — предварительно (выдержка поиска) |
| Mercury 50, рекуператор | — | **нет** ни скрипта, ни results-JSON | README: «расчёт 37,6–39,6 %» при паспорте 38,5 % — из репозитория не воспроизводится |

- КПД SGT-400 и Taurus 60 подогнан: `calib` ищет общую поправку `d` ко всем трём КПД, T41 — под
  мощность. Совпадение КПД — не проверка. Проверка — выхлоп (−6 и +14 °C к входу в коде) и
  прецеденты OP16 и Mercury 50.
- Две модели: `gtu_explore.py` печатает `сверка: PR 6.5, 1000 C -> 24.5 % ; уд. работа 209 кДж/кг;
  T5 586` и `сверка с рекуператором -> 34.8 %`; в `gtu_simple_results.json → простой π=6.5 T41=1000` —
  24,5 % и 599 °C (наборы потерь разные). Печать в JSON не попадает — не цитировать, пока не сохранена.
- ±2 п. п. по КПД — принятая погрешность отбора схем (README), а не вычисленная: схемы, которые
  различаются меньше чем на 2 п. п., моделью не различаются.

## 5. Принятые ограничения (`LIM` в `gtu_explore.py`, оценка)

Окружная скорость колеса компрессора: титан 620, алюминий 520 (T за колесом ≤ 473 К), монотор
заодно с турбиной 480 м/с. Литое радиальное колесо турбины: 500 м/с при газе до 900 °C, выше —
минус 0,7 м/с на градус; газ до 950 °C. Осевая неохлаждаемая лопатка (монокристалл): 420 м/с, газ
до 1050 °C. Газ на входе в рекуператор ≤ 700 °C. Цельная отливка: монотор 0,45 м («да») и до 0,60
(«на пределе»), радиальное колесо 0,50/0,70 м, осевое 0,50 м (больше — диск и лопатки); титановое
колесо ≤ 0,90 м. Классы охлаждения `gtu_dropin.py` и `gtu_family.py`: U — T41 ≤ 1000 °C; N — ≤ 1050 °C,
охлаждаемый сопловой аппарат 3 %; K — ≤ 1200 °C, 5,5 % + 5,5 %. Прочность (`platform/fea/README.md`,
условные профили): 620 м/с для титана держится только при ступице без сквозного отверстия.

## 6. Новая схема в перебор

1. Строка `add(fam=…, name=…, N=…, spools=…, cstages=…, comp=…, rotor=…, ggt=…, shaft=…, pt=…,
   eps=…, use=…, drive=…)` — **в конец** `concepts()`, после семейства E. Вставка в середину
   сдвигает `id` всех следующих: «схема 130» в README, `platform/geom/impeller.py --case core130`, папки
   `platform/fea/runs/*130*` станут указывать на чужую схему.
2. Только значения, которые разбирают и `evaluate`, и `finalize`: `comp` ∈ ti/al/mono, `rotor` ∈
   two/mono, `ggt` ∈ radial/ax1/ax2, `shaft` ∈ single/free, `pt` ∈ None/radial/ax1/ax2, `drive` ∈
   direct/fast, `ic` — промохлаждение. `pt="radial"` в коде есть, но ни одной схемой не опробован —
   при первом применении сверить с `gtu_cycle.cycle` на той же точке.
3. Новый узел или материал — новый ключ `LIM` с пометкой «принято (оценка)» и строка в README.
   Старые пределы не трогать: от них зависят все 249 строк.
4. Прогон в копии и сравнение (раздел 9): старые `id` совпадают полностью, строка `сверка:` та же,
   добавились только новые `id` в хвосте `res`.
5. Новая схема `feasible: false` — сначала `evaluate` и `finalize` напрямую: `optimize` глотает
   любое исключение (`except Exception: r = None`), опечатка в ключе выглядит как «невыполнима».
6. Затем запуск на месте, «249» в README → новое число, каталог, PR со ссылкой на issue #59.

## 7. Обязательные правила

- **Каждое опубликованное число — из results-JSON.** README, дневник, описание PR: число стоит в
  `*_results.json` рядом со скриптом, с тем же округлением, рядом — «файл → ключ». Печать в
  терминал, журнал, пересказ — не источник.
- Чего нет в JSON (Mercury 50), не публиковать, пока нет скрипта; до тех пор — запись в
  `data/known_gaps.json` (CLAUDE.md, «Как добавить этап», п. 4).
- КПД узлов, потери, пределы — «принято (оценка)»; паспортные числа — со ссылкой (навык `gtu-fleet`);
  число из выдержки поиска — «предварительно» и не в заголовок.
- Изменился вход калибровки — отдельный PR: пересчёт `gtu_study.py`, `gtu_point.py` и зависимых,
  таблица «было/стало». Тихо заменять число в коде нельзя.
- Расчёт дольше 20 минут — не здесь; облачные прогоны — `platform/aws`, живой запуск только флагом.

## 8. Типовые ошибки

- Скрипт запущен в репозитории ради пробы: перезаписан results-JSON, сменилась `date`. Вернуть —
  `git checkout -- routes/gtu/<скрипт>_results.json`.
- Правка в `cycle` без такой же в `cycle_r`: это две копии одного кода, `recup` идёт через `cycle_r`.
- КПД на валу сравнили с клеммами: паспорт GS/PG — клеммы генератора, MD/CS — вал.
- T41 модели приравняли к температуре из паспорта: определения разные (раздел 2).
- Ждали, что калибровка дойдёт до перебора: `gtu_dropin.py` и `gtu_family.py` берут литеральные
  `eta_pt_gg=0.865`, `eta_pt_pt=0.885`, а не вывод `gtu_study.py`; у `gtu_explore.py` своя модель.
- Вход SGT-400 в коде расходится со страницей Siemens: 36,2 против 36,4 %, π 16,8 против 17,1, 39,4
  против 39,1 кг/с, 555 °C — колонка PG (у MD 543 °C). Выхлоп модели 549 °C — −6 °C к входу в коде,
  но +6 °C к странице. Решать отдельным PR (раздел 7), не по ходу.
- Схема 130 — алюминиевое колесо, U₂ 505 м/с при пределе 520 (`res[id=130]`); титановый двойник с
  теми же размерами — `id` 139; прочность в `platform/fea` считала титан. Номер — вместе с материалом.
- `brentq` бросает `ValueError` («f(a) and f(b) must have different signs»): доля топлива ищется в
  0,005…0,05, давление за ТВД — в 1,05·P0…0,95·p4; при малой π или крайней T41 корня нет.
  `gtu_dropin.py` такие точки молча пропускает (`except Exception: continue`) — проверять пустоты.
- Размеры перебора — при 6 МВт: `optimize` не передаёт `P`. Другая мощность — `evaluate(c, PR, TIT, P=…)`
  напрямую, затем `finalize` (он берёт расход из результата `evaluate`).

## 9. Проверка результата

Полный повтор в копии и сравнение с репозиторием (из корня репозитория):

```bash
R=$PWD/routes/gtu; T=$(mktemp -d); cp -r "$R" "$T"/ && cd "$T"/gtu && rm -rf platform __pycache__ *_results.json
L="gtu_cycle gtu_study gtu_point gtu_simple gtu_dropin gtu_explore"
for f in $L; do python3 $f.py > $f.log 2>&1 || echo "СБОЙ $f"; done; grep -E "сверка|выполнимых" gtu_explore.log
for f in $L; do python3 - ${f}_results.json "$R/${f}_results.json" <<'EOF'
import json, sys
def walk(a, b, p, out):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted((set(a) | set(b)) - {"date", "software"}):
            walk(a.get(k), b.get(k), f"{p}/{k}", out) if k in a and k in b else out.append(f"{p}/{k}: нет пары")
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): out.append(f"{p}: длина {len(a)} против {len(b)}")
        for i, (x, y) in enumerate(zip(a, b)): walk(x, y, f"{p}[{i}]", out)
    elif isinstance(a, float) and isinstance(b, (int, float)):
        if abs(a - b) > 1e-6 * max(1.0, abs(b)): out.append(f"{p}: {a} против {b}")
    elif a != b: out.append(f"{p}: {a!r} против {b!r}")
out = []; walk(json.load(open(sys.argv[1])), json.load(open(sys.argv[2])), "", out)
print(sys.argv[1], "расхождений:", len(out)); print("\n".join(out[:10]))
EOF
done
```

До правки — 0 расхождений во всех шести файлах, печать `выполнимых: 223 из 249` и две строки
`сверка` из раздела 4 (повтор 14.11.2026, около 9,5 мин). После правки расходятся только задуманные
ключи, новые схемы — только новыми `id` в хвосте `res`. Затем из корня:
`node tools/build_catalog.mjs --check`, `node tools/validate.mjs`, `node tools/smoke.mjs`.

## Источники (сверено 14.11.2026)

- Cantera 3.2.0: docstring установленного пакета (`help(ct.Quantity)`, `TPX`, `TPY`, `HPY`, `HPX`,
  `SP`, `enthalpy_mass`, `entropy_mass`, `mean_molecular_weight`); исходник
  https://github.com/Cantera/cantera/blob/v3.2.0/interfaces/cython/cantera/composite.py (стр. 152,
  210, 218: `constant='UV'` по умолчанию, «combine at constant HP»); PyPI —
  https://pypi.org/pypi/cantera/3.2.0/json. https://cantera.org/3.2/python/index.html из контейнера
  не открылся (хост вне разрешённой сети).
- GRI-Mech 3.0: 53 вещества и 325 реакций — проверено на `gri30.yaml` (`n_species`, `n_reactions`);
  http://combustion.berkeley.edu/gri-mech/version30/text30.html не открывалась.
- SciPy 1.18.1 `brentq`: docstring установленного пакета («f(a) and f(b) must have opposite signs»,
  `xtol`); https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html — хост вне
  разрешённой сети.
- Siemens SGT-400, открыта curl:
  https://www.siemens-energy.com/global/en/home/products-services/product/sgt-400.html
- Предварительно, только выдержки поиска (страницы не открывались): Taurus 60 MD —
  https://www.solarturbines.com/en_US/products/mechanical-drive-packages/taurus-60.html; Mercury 50
  (4600 кВт, 38,5 % на клеммах, 2004) —
  https://www.modernpowersystems.com/features/featuremercury-50-recuperated/; OPRA OP16 (1876 кВт эл.,
  25,1 %, π 6,7:1) — https://pdf.directindustry.com/pdf/opra-turbines/op16-datasheet-metric/38704-643648.html
