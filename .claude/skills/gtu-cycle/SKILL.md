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
| цикл на Cantera (GRI-Mech 3.0, метан, ISO) | `gtu_cycle.py` → `gtu_cycle_results.json` | `cycle`, `cycle_r`, `recup`; ключи `SGT400`, `T60` — базовые КПД узлов **без** калибровки |
| калибровка КПД узлов | `gtu_study.py` → `gtu_study_results.json` | `calib`; ключи `SGT400`, `T60`, `grid` |
| точки базового осевого варианта | `gtu_point.py` → `gtu_point_results.json` (читает `gtu_study_results.json`) | `A:…`, `B:…` |
| простые схемы, прецедент OP16 | `gtu_simple.py` → `gtu_simple_results.json` | `OP16-подобный …`, `простой …` |
| перебор 249 схем, размеры колёс | `gtu_explore.py` → `gtu_explore_results.json` | `lim`, `res[]` (поле `id`) |
| Тесла, обращённый Брайтон | `gtu_exotic.py` → `gtu_exotic_results.json` | функции `gtu_explore`; `res[]` |
| добыча: топливный газ, обвод рекуператора | `gtu_upstream.py` → `gtu_upstream_results.json` | `E.evaluate`; `fuel`, `modes` |
| замена в пакетах Taurus 60 и SGT-400 | `gtu_dropin.py` → `gtu_dropin_results.json` | словари `ARCH`, `COOL`, `SLOT`; в JSON `arch`, `cool`, `slot`, `res[]` |
| семейство по паспортам | `gtu_family.py` → `gtu_family_results.json` | вход `fleet/gtu_fleet_<модель>.json` (16 моделей) — навык `gtu-fleet` |

Цель, сверка, принятые величины — `README.md`; платформа — `platform/` (`aws`, `fea`, `geom`, `surrogate`).
Зависимости: `gtu_cycle` нужен `gtu_study`, `gtu_point`, `gtu_simple`, `gtu_dropin`, `gtu_family`;
`gtu_explore` — `gtu_exotic`, `gtu_upstream`, `gtu_dropin` (`comp_stage`), `gtu_family`. Правка модели
меняет результаты всех зависимых.

Импорт любого модуля делает `chdir` в папку модуля. `gtu_cycle`, `gtu_explore`, `gtu_family` пишут
JSON только под `if __name__ == "__main__"`; у остальных шести защиты нет — `import gtu_dropin` (и study,
point, simple, exotic, upstream) прогоняет весь расчёт и перезаписывает свой results-JSON.

## 2. Как устроена модель

**`cycle(PR, T41_C, eta_pc, eta_pt_gg, eta_pt_pt, x_ngv, x_rot, x_leak, dp_comb=0.04, dp_duct=0.01,
p_exh_ratio=1.02, eta_b=0.995, eta_m=0.99, W1=1.0)`** — двухвальная ГТУ со свободной силовой
турбиной. Сжатие и расширение — 40 малых шагов с КПД (политропный процесс) при замороженном
составе; горение полное (CH₄ + 2O₂ → CO₂ + 2H₂O), недожог `eta_b` — потеря тепла. `T41_C` — газ перед
рабочим колесом 1-й ступени после подмешивания охлаждения соплового аппарата, а не T4 камеры и не
температура из паспорта. Доли от `W1`: `x_ngv` возвращается до ротора 1, `x_rot` — после ТВД,
`x_leak` теряется. Возврат: `P_MW`, `eff` (на валу, по низшей теплоте 50,03 МДж/кг), `sp_work_kJkg`,
`EGT_C`, `W_exh`, `T3_C`, `T4_C`, `T45_C`, `GG_ER`, `PT_ER`, `wc_kJkg`, `Wf`, а также `PR`, `T41_C`.
Смешение — `ct.Quantity(g, constant="HP")` (по умолчанию `Quantity` смешивает при постоянных UV;
при HP давления потоков обязаны совпадать, иначе `ValueError`); корни — `scipy.optimize.brentq`.
`cycle_r(…, T3r=None, dp_rec_air=0.0)` — копия `cycle` с подогретым воздухом перед камерой (воздух
охлаждения — по-прежнему при T3); `recup(PR, T41_C, eps, dp_air=0.025, dp_gas=0.035, **kw)` —
рекуператор: первый расчёт и 6 итераций по выхлопу, добавляет `T3r_C` и `Tstack_C` (оценка через
0,93). На клеммы — `eff * 0.985 * 0.972` (редуктор и генератор; `gtu_simple.py`, `gtu_point.py`,
`ETA_GEN` в `gtu_family.py`; оценка).

**`gtu_explore.py`** — своя быстрая модель: свойства из таблиц Cantera с шагом 2 К, смесь по доле
продуктов, утечка 4 %, мощность 6 МВт.
- `comp_stage(T1, p1, PRs, m, kind, dbl=False)`, `kind` ∈ `ti`/`al`/`mono`: КПД 0,880 −
  0,010·(πₛ − 2,5) (−0,02 для монотора) − 0,3·зазор/b₂; U₂ = √(w/0,70); D₂ из m = ρ₁·U₂·D₂²·0,085;
  n = 60·U₂/(π·D₂); `ok` — U₂ в пределе, для алюминия ещё T₂ ≤ 473 К; `dbl` — двусторонний вход.
- `evaluate(c, PR, TIT_C, P=6.0e6)` — цикл схемы `c`; `None` при нарушении предела, `w_net ≤ 0` или
  степени расширения СТ < 1,15.
- `finalize(c, r, P=6.0e6)` — размеры, «литьё» (`да`/`на пределе`), индекс сложности `K`, `eff_el`;
  `None` — отливка больше предела, титановое колесо больше 0,90 м или СТ на 3000 об/мин не
  выходит по размерам (h/D < 0,02 или D > 2,4 м). Аргумент `P` внутри не используется. `eff_el`:
  через редуктор `eff·0,985·0,972`, без редуктора `eff·0,972`, кластер без СТ `eff·0,955`.
- `optimize(c)` — TIT от 800 °C до предела (950 для радиальной ТГГ, 1050 для осевой) шагом 25,
  π 2,5…16 шагом 0,25, у монотора `kD` 1,0/1,1/1,2; лучшая — по наибольшему `eff`, не по простоте;
  `evaluate` зовёт без `P`. После вызова в `c` остаётся `kD` последнего перебора.
- `concepts()` — список схем, семейства A (63), B (54), C (108), E (24); `id` = позиция в списке + 1.

## 3. Порядок работы

**Шаг 0. Среда.** `pip install cantera scipy numpy`. Сверено 14.11.2026: cantera 3.2.0 (PyPI:
Python ≥ 3.10 и < 3.15), scipy 1.18.1, numpy 2.5.3, Python 3.13; `gri30.yaml` входит в пакет Cantera.

**Шаг 1. Копия.** Скрипт пишет свой results-JSON в свою папку (с новой `date`), импорт меняет
текущую папку (раздел 1). Пробы — только в копии; `platform/` (300 МБ) копировать не нужно, `fleet/` —
только для `gtu_family.py`. Из корня репозитория:

```bash
T=$(mktemp -d); mkdir "$T"/gtu && cp routes/gtu/*.py "$T"/gtu/ && cd "$T"/gtu
python3 -c "
from gtu_cycle import cycle, recup
import gtu_explore as E
r = cycle(6.75, 1050, eta_pc=0.855, eta_pt_gg=0.865, eta_pt_pt=0.885, x_ngv=0.03, x_rot=0.0, x_leak=0.015)
print(round(r['eff']*100, 1), round(r['sp_work_kJkg']), round(r['EGT_C']))
st = E.comp_stage(E.T0, E.P0, 6.75, 22.5, 'ti'); print(round(st['D2'], 2), round(st['rpm']), round(st['U2']), st['ok'])
c = E.concepts()[129]; print(c['name'], E.optimize(c))   # схема id 130
"
```

Ожидается (повторено 14.11.2026, около 3 с): `27.1 256 627` — как
`gtu_dropin_results.json → res[arch=1Ц, cool=N]`; `0.59 19951 618 True` — в JSON 19 967 об/мин:
там расход не округлён (22,464 кг/с); схема 130 совпадает с `gtu_explore_results.json → res[id=130]`
полностью (печать с `np.float64(…)` — это numpy). Разница в последнем знаке — довод брать число из JSON.

**Шаг 2. Длительность.** Замер автора навыка 14.11.2026, по одному: `gtu_cycle.py` 12 с,
`gtu_study.py` 99 с, `gtu_point.py` 17 с, `gtu_simple.py` 50 с, `gtu_dropin.py` 312 с,
`gtu_explore.py` 79 с; все — в одном процессе (параллелит только `gtu_family.py`, 4 процесса).
Повтор 14.11.2026: `gtu_family.py` 308 с, `gtu_exotic.py` и `gtu_upstream.py` — по 1 с. Параллель не помогает: четыре группы скриптов
одновременно на 4 ядрах заняли 613 с против 569 с по одному, каждый скрипт шёл в 2–5 раз медленнее.
Больше 20 минут — делить.

**Шаг 3. Сверка** — разделы 4 и 9: после правки модели пересчитать калибровку, прецеденты и все
зависимые скрипты (раздел 1).

**Шаг 4. Публикация.** Скрипт запускается на месте; скрипт и results-JSON — первым коммитом;
каталог (`node tools/build_catalog.mjs`, трек `gas-turbine` по префиксу `gtu_`) — в чистой копии
HEAD (`git worktree add --detach`) вторым (`.agent/LESSONS.md`, п. 1); затем проверки CLAUDE.md.
На вкладки сайта трек пока не выносится.

## 4. Калибровка и сверка

| объект | вход в коде | модель (results-JSON) | замечание |
|---|---|---|---|
| SGT-400 13 МВт, привод | 13,4 МВт; 36,2 %; π 16,8; выхлоп 39,4 кг/с; 555 °C | `gtu_study_results.json → SGT400`: КПД узлов 0,888/0,868/0,888, T41 1219 °C, выхлоп 549 °C | страница Siemens 14.11.2026, MD 13 МВт: 13,4 МВт; 36,4 %; 9886 кДж/(кВт·ч); π 17,1; 39,1 кг/с; 543 °C. 555 °C и 40,0 кг/с — колонка PG (12,9 МВт эл., 35,00 %). Сверенный паспорт — `fleet/gtu_fleet_sgt400.json → records[variant="MD 13 MW"]`, расхождение с кодом — там же, `common.code_note` |
| Taurus 60, привод | 5,74 МВт; 31,96 %; π 12,2; 21,63 кг/с; 510 °C; воздух `W1` = 21,63/1,004 = 21,54 кг/с | `→ T60`: 0,874/0,864/0,884, T41 1082 °C, выхлоп 524 °C | Solar DS60MD/713: 5740 кВт; 11 265 кДж/(кВт·ч) (= 31,96 %); 77 880 кг/ч; 510 °C; π 12,2; воздух на входе 21,3 кг/с. Паспорт — `fleet/gtu_fleet_taurus60.json → records[variant=MD]` |
| OPRA OP16, радиальная | π 6,7 и T41 приняты, 8,8 кг/с; паспорт 1,876 МВт эл., 25,1 % | `gtu_simple_results.json → OP16-подобный одновальный T41=1000`: `eff_el` 0,237, `P_el` 1,795 МВт | T41 950/1050 °C: 23,2/24,1 %. Паспорт — предварительно (выдержка поиска) |
| Mercury 50, рекуператор | — | **нет** ни скрипта, ни results-JSON | README: «расчёт 37,6–39,6 %» — из репозитория не воспроизводится. Паспорт Solar DSM50PG/1118: 4600 кВт эл.; 9350 кДж/(кВт·ч); 38,5 % на клеммах; π 9,9; воздух 17,9 кг/с; 365 °C — `fleet/gtu_fleet_mercury50.json` |

- КПД SGT-400 и Taurus 60 подогнан: `calib` ищет общую поправку `d` ко всем трём КПД (−0,0066 и
  −0,0007), T41 — под мощность. Совпадение КПД — не проверка. Проверка — выхлоп (−6 и +14 °C к входу
  в коде) и прецеденты OP16 и Mercury 50.
- `gtu_cycle_results.json → SGT400/T60` — не калибровка: там базовые КПД без `d`, T41 подобрана под
  мощность, КПД 36,84 и 32,02 %, выхлоп 535,6 и 522,7 °C. README («549 °C при 555 °C») — из
  `gtu_study_results.json`.
- На странице Siemens градусы Фаренгейта и фунты не сходятся с соседними числами: MD 543 °C
  (1,020 °F = 549 °C), PG 555 °C (1,019 °F = 548 °C), PG 40,0 кг/с (86,3 lb/s = 39,1 кг/с). Выхлоп
  модели 549 °C — от −6 до +6 °C к этим числам; точнее проверка не различает.
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
колесо ≤ 0,90 м (то же `D2_MAX` в `gtu_family.py`). Классы охлаждения `gtu_dropin.py` и
`gtu_family.py`: U — T41 ≤ 1000 °C; N — ≤ 1050 °C, охлаждаемый сопловой аппарат 3 %; K — ≤ 1200 °C,
5,5 % + 5,5 %; `gtu_dropin.py` считает ровно на верхней границе класса. Прочность
(`platform/fea/README.md`, условные профили): 620 м/с для титана держится только при ступице без
сквозного отверстия.

## 6. Новая схема в перебор

1. Строка `add(fam=…, name=…, N=…, spools=…, cstages=…, comp=…, rotor=…, ggt=…, shaft=…, pt=…,
   eps=…, use=…, drive=…)` (при нужде `dbl=…`, `ic=…`) — **в конец** `concepts()`, после семейства E.
   Вставка в середину сдвигает `id` всех следующих: «схема 130» в README,
   `platform/geom/impeller.py --case core130`, папки `platform/fea/runs/*130*` станут указывать на
   чужую схему.
2. Только значения, которые разбирают и `evaluate`, и `finalize`: `comp` ∈ ti/al/mono, `rotor` ∈
   two/mono, `ggt` ∈ radial/ax1/ax2, `shaft` ∈ single/free, `pt` ∈ None/radial/ax1/ax2, `drive` ∈
   direct/fast (обязателен при `pt`), `use` ∈ gen/both (`gen` с `fast` — редуктор), `dbl` —
   двусторонний вход, `ic` — промохлаждение. Неизвестное значение `ggt`/`pt` код молча считает двумя
   осевыми ступенями. `pt="radial"` в коде есть, но ни одной схемой не опробован — при первом
   применении сверить с `gtu_cycle.cycle` на той же точке.
3. Новый узел или материал — новый ключ `LIM` с пометкой «принято (оценка)» и строка в README.
   Старые пределы не трогать: от них зависят все 249 строк.
4. Новая проверка с `return None` внутри цикла `evaluate`: служебные `c["_T3"]`, `c["_m"]` снимаются
   только при нормальном выходе — после досрочного `None` они перейдут в следующий вызов той же схемы.
5. Прогон в копии и сравнение (раздел 9): старые `id` совпадают полностью, строка `сверка:` та же,
   добавились только новые `id` в хвосте `res`.
6. Новая схема `feasible: false` — сначала `evaluate` и `finalize` напрямую на свежем словаре из
   `concepts()`: `optimize` глотает любое исключение `evaluate` (`except Exception: r = None`), опечатка
   в ключе выглядит как «невыполнима»; исключение в `finalize` не ловится и роняет весь перебор.
7. Затем запуск на месте, «249» в README → новое число, каталог, PR со ссылкой на issue #59.

## 7. Обязательные правила

- **Каждое опубликованное число — из results-JSON.** README, дневник, описание PR: число стоит в
  `*_results.json` рядом со скриптом, с тем же округлением, рядом — «файл → ключ». Печать в
  терминал, журнал, пересказ — не источник.
- Чего нет в JSON (Mercury 50), не публиковать, пока нет скрипта; до тех пор — запись в
  `data/known_gaps.json` (CLAUDE.md, «Как добавить этап», п. 4; поля записи `id`, `reference`,
  `cited_in`, `what_happened`, `impact`, `resolution`, `noted`). На 14.11.2026 записи о Mercury 50 нет.
- КПД узлов, потери, пределы — «принято (оценка)»; паспортные числа — со ссылкой (навык `gtu-fleet`);
  число из выдержки поиска — «предварительно» и не в заголовок.
- Изменился вход калибровки — отдельный PR: пересчёт `gtu_study.py`, `gtu_point.py` и зависимых,
  таблица «было/стало». Вход повторён литералами в четырёх файлах: `gtu_cycle.py` (блок `__main__`),
  `gtu_study.py` (`calib`), `gtu_dropin.py` (`SLOT`), `gtu_upstream.py` (`cases`) — менять все сразу.
  Тихо заменять число в коде нельзя.
- Расчёт дольше 20 минут — не здесь; облачные прогоны — `platform/aws`, живой запуск только флагом.

## 8. Типовые ошибки

- Скрипт запущен в репозитории ради пробы (или импортирован, раздел 1): перезаписан results-JSON,
  сменилась `date`. Вернуть — `git checkout -- routes/gtu/<скрипт>_results.json`.
- Правка в `cycle` без такой же в `cycle_r`: это две копии одного кода, `recup` идёт через `cycle_r`.
- КПД на валу сравнили с клеммами: паспорт GS/PG — клеммы генератора, MD/CS — вал.
- T41 модели приравняли к температуре из паспорта: определения разные (раздел 2).
- Процитировали `gtu_cycle_results.json → SGT400` как калибровку: там КПД 36,84 % и выхлоп 535,6 °C
  без поправки `d` (раздел 4).
- Ждали, что калибровка дойдёт до перебора: `gtu_dropin.py` и `gtu_family.py` берут литеральные
  `eta_pt_gg=0.865`, `eta_pt_pt=0.885`, а не вывод `gtu_study.py`; у `gtu_explore.py` своя модель.
- Вход SGT-400 в коде расходится со страницей Siemens: 36,2 против 36,4 %, π 16,8 против 17,1, 39,4
  против 39,1 кг/с, 555 °C — колонка PG (у MD 543 °C). Решать отдельным PR (раздел 7), не по ходу.
- Схема 130 — кластер из 4 ядер с рекуператором ε 0,85 и СТ на 3000 об/мин: `m` 39,8 кг/с — на все
  ядра (на ядро около 10), `Dc` 0,435 м и 22 197 об/мин — одного ядра. Колесо алюминиевое, U₂ 505 м/с
  при пределе 520 (`res[id=130]`); титановый двойник с теми же размерами — `id` 139; прочность в
  `platform/fea` считала титан. Номер — вместе с материалом.
- `gtu_upstream.py` держит числа схемы 130 литералами (π 4,25, 950 °C, 33,7 %, 39,8 кг/с) — после
  правки перебора сверить их с `res[id=130]` вручную.
- `brentq` бросает `ValueError` («f(a) and f(b) must have different signs»): доля топлива ищется в
  0,005…0,05, давление за ТВД — в 1,05·P0…0,95·p4; при малой π или крайней T41 корня нет.
  `gtu_dropin.py` такие точки молча пропускает (`except Exception: continue`), `gtu_family.py` —
  тоже (`_run` возвращает `None`) — проверять пустоты.
- Размеры перебора — при 6 МВт: `optimize` не передаёт `P`. Другая мощность — `evaluate(c, PR, TIT, P=…)`
  напрямую, затем `finalize` (он берёт расход из результата `evaluate`). Словарь — свежий из
  `concepts()`, для монотора `c["kD"]` = `res[id]["kD"]`: после `optimize` в `c` лежит последний `kD`.

## 9. Проверка результата

Полный повтор в копии и сравнение с репозиторием (из корня репозитория):

```bash
L="gtu_cycle gtu_study gtu_point gtu_simple gtu_dropin gtu_explore gtu_exotic gtu_upstream gtu_family"
R=$PWD/routes/gtu; T=$(mktemp -d); mkdir "$T"/gtu "$T"/ref && cp "$R"/*.py "$T"/gtu/ && cp -r "$R"/fleet "$T"/gtu/
for f in $L; do git show HEAD:routes/gtu/${f}_results.json > "$T"/ref/${f}_results.json; done
cd "$T"/gtu
for f in $L; do python3 $f.py > $f.log 2>&1 || echo "СБОЙ $f"; done; grep -E "сверка|выполнимых" gtu_explore.log
for f in $L; do python3 - ${f}_results.json ../ref/${f}_results.json <<'EOF'
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

Порядок в `L` важен: `gtu_point.py` читает `gtu_study_results.json` этого же прогона. Эталон — JSON
из HEAD, а не из рабочей папки: 14.11.2026 во время прогона чужой коммит (83355f2) сменил паспорта
`fleet/` и `gtu_family_results.json`, и сравнение с рабочей папкой дало 2537 ложных расхождений.
Скрипты и `fleet/` берутся из рабочей папки — так проверяется ваша правка. Повтор 14.11.2026: без
правок 0 расхождений во всех девяти файлах, печать `выполнимых: 223 из 249` и две строки `сверка` из
раздела 4; весь прогон — около 15 мин (шаг 2).
После правки расходятся только задуманные ключи, новые схемы — только новыми `id` в хвосте `res`.

Затем каталог — в чистой копии HEAD после коммита всех правок (`.agent/LESSONS.md`, п. 1; в рабочей
папке лежат сгенерированные `platform/geom/out/*.json`, и проверка там ничего не доказывает):

```bash
wt=$(mktemp -d) && git worktree add --detach "$wt" HEAD && (cd "$wt" && node tools/build_catalog.mjs --check); git worktree remove --force "$wt"
node tools/validate.mjs && node tools/smoke.mjs
```

## Источники (сверено 14.11.2026)

- Cantera 3.2.0: исходник
  https://github.com/Cantera/cantera/blob/v3.2.0/interfaces/cython/cantera/composite.py — файл с
  raw.githubusercontent.com побайтно совпал с установленным пакетом; стр. 152 — класс `Quantity`,
  196 и 210 — примеры «combine at constant UV» / «combine at constant HP», 218 —
  `__init__(…, constant='UV')`, 332–345 — сложение (при HP — `ValueError`, если давления не равны).
  Docstring установленного пакета (`TPX`, `TPY`, `HPY`, `HPX`, `SP`, `enthalpy_mass`, `entropy_mass`,
  `mean_molecular_weight`); PyPI — https://pypi.org/pypi/cantera/3.2.0/json (`requires_python`
  `<3.15,>=3.10`). https://cantera.org/3.2/python/index.html из контейнера не открывается (хост вне
  разрешённой сети).
- GRI-Mech 3.0: 53 вещества и 325 реакций — проверено на `gri30.yaml` (`n_species`, `n_reactions`);
  http://combustion.berkeley.edu/gri-mech/version30/text30.html не открывается.
- SciPy 1.18.1 `brentq`: текст ошибки «f(a) and f(b) must have different signs» —
  https://github.com/scipy/scipy/blob/v1.18.1/scipy/optimize/_zerosmodule.c (стр. 110) и вызов в
  установленном пакете; docstring — «f(a) and f(b) must have opposite signs», `xtol`.
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html — хост вне
  разрешённой сети.
- Siemens SGT-400, открыта curl, две таблицы (Simple cycle power generation, Mechanical drive):
  https://www.siemens-energy.com/global/en/home/products-services/product/sgt-400.html
- Solar, листы данных открыты curl и `pdftotext`: Taurus 60 MD, DS60MD/713 —
  https://s7d2.scene7.com/is/content/Caterpillar/C10550209 («Output Power 5740 kW», «Heat Rate 11 265
  kJ/kW-hr», «Exhaust Flow 77 880 kg/hr», «Exhaust Temp. 510°C», «Pressure Ratio: 12.2:1», «Inlet
  Airflow: 21.3 kg/sec», «Engine efficiency: 32%»); Mercury 50 PG, DSM50PG/1118 —
  https://s7d2.scene7.com/is/content/Caterpillar/CM20150710-52396-21070 («Output Power 4600 kWe»,
  «Heat Rate 9350 kJ/kWe-hr», «Compression Ratio: 9.9:1», «17.9 kg/sec», «Exhaust Temperature 365°C»,
  «Engine efficiency: 38.5%»). solarturbines.com из контейнера не открывается.
- Предварительно, только выдержка поиска (страница из контейнера не открывается): OPRA OP16
  (1876 кВт эл., 25,1 %, π 6,7:1) — https://pdf.directindustry.com/pdf/opra-turbines/op16-datasheet-metric/38704-643648.html
