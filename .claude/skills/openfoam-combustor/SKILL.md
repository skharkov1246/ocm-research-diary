---
name: openfoam-combustor
description: OpenFOAM для камеры сгорания и теплообмена трека «Простая литая ГТУ-заменитель» — выбор версии (openfoam.com v2606 или openfoam.org v14) и установка, решатели reactingFoam / rhoReactingFoam / XiFoam / sprayFoam / chtMultiRegionFoam и их замены foamRun / foamMultiRun в v14, кинетика GRI-Mech 3.0 и сокращённые механизмы, диффузионное горение и горение с предварительным смешением, распыл жидкого топлива для перехода газ — дизель, сопряжённый теплообмен со стенкой жаровой трубы, проверка на открытом эталоне Sandia Flame D, оценка NOx и CO с приведением к 15 % O2, прогон на AWS через routes/gtu/platform/aws/runner.py. Используй, когда ставишь OpenFOAM, собираешь или считаешь кейс камеры сгорания, выбираешь механизм, считаешь температуру стенки или выбросы двигателя-заменителя. Триггеры — «камера сгорания», «OpenFOAM», «reactingFoam», «chtMultiRegionFoam», «Sandia Flame D», «GRI-Mech», «NOx», «CO на выхлопе», «распыл дизеля», «температура жаровой трубы».
---

# OpenFOAM: камера сгорания и теплообмен двигателя-заменителя

Платформа трека: SU2 — ступени (навык `su2-turbo`), OpenFOAM — горение и теплообмен, CalculiX — прочность
(`routes/gtu/platform/fea/`), Gmsh / Open CASCADE — геометрия. Камера получает воздух за центробежным колесом
и отдаёт газ на литую турбину без охлаждения или с охлаждаемым первым сопловым аппаратом, поэтому от расчёта
нужны три вещи: поле температуры на выходе, температура стенки жаровой трубы, NOx и CO. Цикл и размеры ядра — навык
`gtu-cycle`, интерфейсы пакета — `gtu-dropin`.

Состояние на 14.11.2026: кейса камеры в репозитории нет, боевых прогонов OpenFOAM не было. Проверено в контейнере
(4 ядра, Ubuntu 24.04): apt-пакет `openfoam` 1912 ставится; учебный `counterFlowFlame2D` из исходников v2606
прошёл на нём до 0,5 с за 39 с на одном ядре (T max 2016 K; блок `functions` из controlDict удалён); chemkinToFoam
переводит GRI-Mech 3.0. Официальные репозитории `dl.openfoam.org` и `dl.openfoam.com` из контейнера закрыты
прокси (403), поэтому v2606 и v14 здесь не ставились: имена решателей, словарей и учебников сверены по исходным
кодам (клоны GitHub OpenFOAM-14 и зеркала gitlab.com/openfoam/core/openfoam, api=2606).

## 1. Что где лежит

| что | где |
|---|---|
| вход камеры: `PR`, `T3_C`, `T4_C`, `T41_C`, `Wf` (топливо), `W_exh` | `routes/gtu/gtu_cycle_results.json` (`SGT400`, `T60`) |
| вход по схемам замены: `PR`, `T3`, `T41`, расход `G` в пакете | `routes/gtu/gtu_dropin_results.json` → `res[]` |
| ядро схемы 130 (колесо 0,435 м, расход около 10 кг/с, PR 4,25) | `routes/gtu/gtu_explore_results.json` |
| прогоны на AWS (спот, потолок 300 долл., выгрузка в `$ALPHA_O_S3/gtu/`) | `routes/gtu/platform/aws/runner.py`, `README.md` |
| кейсы камеры | `routes/gtu/platform/combustor/` — предложение, папки ещё нет |
| образцы | учебники версии: `SandiaD_LTS`, `counterFlowFlame2D_GRI_TDAC`, `aachenBomb`, `reverseBurner` (раздел 3) |

Степень сжатия в results-JSON трека — от 4,25 до 16,8, значит давление в камере порядка 4–17 бар. Граничные
условия берутся только из этих файлов, с указанием файла и ключа в README кейса; своих чисел не придумывать.

## 2. Версия и установка

Две ветви с несовместимыми словарями: openfoam.com (OpenCFD), выпуск v2606, и openfoam.org (OpenFOAM
Foundation), выпуск 14 (даты 26.06.2026 и 14.07.2026 — предварительно, по выдаче поиска; номера сверены по
исходникам). **Основная — v2606**: в ней раздельные решатели, названные в техзадании (`reactingFoam`,
`rhoReactingFoam`, `XiFoam`, `sprayFoam`, `reactingParcelFoam`, `chtMultiRegionFoam` с горением в газовой
области — в `createFluidFields.H` есть `CombustionModel<rhoReactionThermo>`), учебник
`XiDyMFoam/annularCombustorTurbine` (в учебниках v14 его нет), и это тот же род, что apt-1912 для местной обкатки.
Четыре образца из раздела 1 и модель LISA есть в обеих ветвях. v14 — допустимая замена по таблице раздела 3.
Один кейс — одна версия; версия (`api`/шапка `Version:` в словарях) пишется в results-JSON.

```bash
# openfoam.com v2606 (по README зеркала SourceForge, там пример 2112 — имя 2606 по шаблону, не проверено)
curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash
sudo apt-get update && sudo apt-get install -y openfoam2606-default
source /usr/lib/openfoam/openfoam2606/etc/bashrc      # или оболочка: openfoam2606
# openfoam.org v14 (по выдаче поиска, страница из контейнера не открылась; сверить перед запуском)
sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc"
sudo add-apt-repository http://dl.openfoam.org/ubuntu && sudo apt update && sudo apt -y install openfoam14
source /opt/openfoam14/etc/bashrc
# местная обкатка в контейнере: Ubuntu 24.04 universe, openfoam 1912.200626 (проверено)
sudo apt-get install -y openfoam && source /usr/share/openfoam/etc/bashrc
```

## 3. Соответствие имён: v2606 → v14

| назначение | openfoam.com v2606 | openfoam.org v14 |
|---|---|---|
| реагирующий газ | `reactingFoam` (psiReactionThermo), `rhoReactingFoam` (rhoReactionThermo) | `foamRun -solver multicomponentFluid`; в controlDict `solver multicomponentFluid;` |
| предварительное смешение, модель Xi | `XiFoam` (psiuReactionThermo), `XiDyMFoam` | `foamRun -solver XiFluid` (переписан в v14) |
| распыл | `sprayFoam`; `reactingParcelFoam` (есть LTS-учебник), `simpleReactingParcelFoam` | `multicomponentFluid` + в `fvModels` `clouds { type clouds; }`, `constant/cloudProperties` (`type sprayCloud;`) |
| сопряжённый теплообмен | `chtMultiRegionFoam`; `constant/regionProperties`: `regions (fluid (gas) solid (solid));` | `foamMultiRun`; controlDict: `regionSolvers { gas multicomponentFluid; solid solid; }` |
| стык газ — стенка | `compressible::turbulentTemperatureRadCoupledMixed` | `coupledTemperature` |
| модель горения | `combustionProperties`: `combustionModel EDC;` + `EDCCoeffs {version v2005;}` | `reactionProperties`: `reactionModel EDC;` + `EDC {version v2005;}` |
| смесь | `thermophysicalProperties`: `mixture reactingMixture;`, `inertSpecie N2;` | `physicalProperties`: `mixture multicomponentMixture;`, `defaultSpecie N2;` |
| прочее | `turbulenceProperties`, `fvOptions`, `odeCoeffs`, `atomizationModel`, `postProcess` | `momentumTransport`, `fvModels`, `ode`, `atomisationModel`, `foamPostProcess` |
| 0D кинетика, перевод CHEMKIN | `chemFoam`, `chemkinToFoam` | `chemFoam`, `chemkinToFoam` |

В v14 `reactingFoam`, `XiFoam`, `chtMultiRegionFoam` остались сценариями в `bin/`, которые печатают замену и
запускают `foamRun`/`foamMultiRun`; `sprayFoam` и `rhoReactingFoam` в `bin/` v14 нет. Модели горения в обеих
ветвях: `laminar`, `EDC`, `PaSR` (`Cmix`), `infinitelyFastChemistry`, `diffusion`, `FSD`, `singleStep…`, `zone…`.

## 4. Порядок работы

- **Шаг 0.** Начало сессии по CLAUDE.md: манифест, каталог, хвосты веток. Один проход — один кейс.
- **Шаг 1. Вход.** Из results-JSON (раздел 1): давление и температура за компрессором, расход воздуха, топлива,
  целевая температура на выходе. Ссылку на файл и ключ — в README кейса.
- **Шаг 2. Механизм** (раздел 5). Природный газ — полный GRI-Mech 3.0; дизель — отдельный механизм.
- **Шаг 3. Эталон.** Sandia Flame D (раздел 6) — до любого числа по своей камере, на той же версии и тех же
  моделях. Без пройденного эталона числа камеры идут с пометкой «не подтверждено эталоном».
- **Шаг 4. Геометрия и сетка.** Gmsh → `gmshToFoam`, затем `checkMesh`; для стенки — `splitMeshRegions -cellZones`.
  Не меньше двух сеток; разница итоговых величин между ними — в results-JSON.
- **Шаг 5. Горение.** Диффузионное и частично перемешанное — `reactingFoam` + `EDC` (version v2005) или `PaSR`.
  Стационар — псевдовремя LTS: `ddtSchemes { default localEuler; }`, как в `SandiaD_LTS`. Как в образце: сначала
  без химии (`chemistry off`, до 1500 итераций), затем с химией. Бедное предварительное смешение: `XiFoam` даёт
  фронт и тепловыделение, но не состав газа — NOx и CO по нему не считаются; для них — `reactingFoam` с механизмом.
  Ускорение химии — TDAC: `chemistryType { solver ode; method TDAC; }`, блоки `reduction` (DAC) и `tabulation` (ISAT).
- **Шаг 6. Распыл (переход газ — дизель).** Образец — `aachenBomb`: `coneNozzleInjection`, распад `ReitzDiwakar`,
  испарение `liquidEvaporationBoil`, теплообмен капли `RanzMarshall`; для центробежной форсунки есть `LISAAtomization`.
  Жидкости в библиотеке обеих ветвей: `C12H26`, `C16H34`, `C7H16`, `IDEA` (70 % н-декана + 30 % α-метилнафталина).
  `sprayFoam` нестационарен; установившийся режим — `reactingParcelFoam` с LTS (учебник `verticalChannelLTS`).
- **Шаг 7. Стенка.** `chtMultiRegionFoam` с областями `gas` (горение) и `solid` (жаровая труба), образец
  `reverseBurner`; излучение — `radiationModel P1` (как в образцах). Свойства металла брать из
  `routes/gtu/platform/fea/materials.json`, если там есть нужный сплав; иначе — источник со ссылкой.
- **Шаг 8. Выбросы** (раздел 7), затем второй способ — Cantera.
- **Шаг 9. Запись.** Кейс, `Allrun`, журнал и results-JSON с версией, механизмом, сетками, невязками, статусом
  эталона; чего нет — в `data/known_gaps.json`; затем `node tools/build_catalog.mjs` и проверки CLAUDE.md.

## 5. Кинетика

- **Полный GRI-Mech 3.0** — 53 вещества, 325 реакций, с химией азота (NO, NO2, N2O, HCN, NH3). Готов в формате
  OpenFOAM: `counterFlowFlame2D_GRI_TDAC/constant/reactionsGRI` и `thermo.compressibleGasGRI` (проверено подсчётом;
  в Cantera 3.2 `gri30.yaml` — те же 53/325). Разработчики велят брать механизм как есть: замена констант или
  выброс веществ ухудшает его; поведение проверено для метана, этана, водорода и CO.
- **Диапазон по давлению** для камеры 4–17 бар не сверен: README30 с сайта Berkeley из контейнера не открылся, а
  вторичные источники расходятся. До сверки результаты выше 10 атм — с пометкой.
- **Сокращённый механизм:** не вырезать вещества руками, а включать DAC/ISAT (TDAC) на полном GRI.
- **Дизель:** в учебнике `aachenBomb` — одна глобальная реакция н-гептана (5 веществ); для CO и NOx не годится.
  Механизм для C12H26 в учебниках нет — выбирать отдельно, со ссылкой (не выбран, не проверен).
- **Перевод Cantera YAML → OpenFOAM:** `yaml2ck mech.yaml --mechanism m.inp --thermo t.dat --transport tr.dat`, затем
  `chemkinToFoam m.inp t.dat <transportProperties> reactions thermo` (третий аргумент — словарь OpenFOAM, образец
  `SandiaD_LTS/chemkin/transportProperties`). Файлы yaml2ck (Cantera 3.2) chemkinToFoam 1912 не читает; помогли три
  правки (проверено, 53/325): `THERMO` → `THERMO ALL` и строка диапазона температур в колонках `%10.3f` без пустых
  строк после неё; `Ar` → `AR` в `ELEM`; `REACTIONS CAL/MOLE MOLE` → `… MOLES`. На v2606/v14 не проверялось.

## 6. Эталон: Sandia Flame D

Пилотируемая струя (предварительно, по выдаче поиска — ERCOFTAC AC2-09; перед сравнением сверить с архивом TNF):
топливо 25 % CH4 + 75 % воздуха по объёму, сопло D = 7,2 мм, среднерасходная скорость 49,6 м/с, Re = 22 400;
пилот — кольцо 7,7/18,2 мм, 11,4 м/с, бедная смесь при 1880 K; местного срыва почти нет. Измерены T, доля
смешения, N2, O2, H2O, H2, CH4, CO, CO2, OH и NO (OH и NO — лазерная флуоресценция). Первоисточник — Barlow, Frank,
Proc. Combust. Inst. 27 (1998) 1087–1095. Учебник `SandiaD_LTS` повторяет геометрию (проверено по `blockMeshDict`:
радиус сопла 3,6 мм, пилот 3,85–9,1 мм, `scale 0.001`) и считает RANS `kEpsilon` + EDC v2005 + P1 + TDAC.

1. Прогнать учебник как есть: `Allrun` (chemkinToFoam → blockMesh → setFields → без химии до 1500 → с химией до 5000).
2. Сравнить с данными TNF осевую T и долю смешения и радиальные профили на z/D = 7,5; 15; 30; 45; 60; 75 (станции
   `system/sampleDict` в v2606). Критерий согласия записать до прогона.
3. Повторить с полным GRI (раздел 5) и сравнить CO и NO — только после этого считать выбросы своей камеры.
4. Чувствительность к `kEpsilonCoeffs { C1 …; }`: по умолчанию 1,44 (проверено в исходнике); для круглой струи в
   литературе берут 1,6 (Pope, 1978). Прогон с 1,6 в этом проекте не делался.

Время прогона эталона на 4 ядрах не замерено. Местный прогон — не дольше 20 минут; длиннее — на AWS (раздел 8).

## 7. NOx и CO

1. Мольные доли: функция `moleFractions` (`postProcess -func moleFractions` в v2606, `foamPostProcess` в v14).
2. Средние по выходу с весом расхода: `surfaceFieldValue` на патче выхода, `weightField phi;`, операция
   `weightedAverage` (v2606); в списке операций v14 весовых нет, вес задаётся `weightField` — проверить на прогоне.
3. Сухой газ: X_сух = X / (1 − X_H2O). NOx = NO + NO2 в ppm сухого газа.
4. Приведение к 15 % O2: C₁₅ = C_сух · (20,9 − 15) / (20,9 − O2,сух[%]). Нормы 40 CFR 60 подраздел KKKK даны
   «at 15 percent O2»; дословно формулу в тексте правила не сверял.
5. Второй способ: Cantera 3.2 (`ct.Solution('gri30.yaml')`, `ct.IdealGasConstPressureReactor` в `ct.ReactorNet`) при
   тех же p, T, коэффициенте избытка воздуха и времени пребывания из CFD. Расхождение — остановка и разбор.
6. Отдельной модели NOx (постобработки по Зельдовичу) в исходниках обеих ветвей нет: NO берётся только из механизма.

## 8. Прогон на AWS

Ключи AWS в этом контейнере не действуют — реальных вызовов не делать, проверять задачу офлайн
(`runner.load_jobs`). Поле `"apt": ["openfoam"]` поставит на Ubuntu 24.04 версию 1912 — нужную версию ставить в
`command` (команда идёт от root, из корня клона). Пример прошёл `load_jobs`; боевого запуска не было:

```json
{"name": "of-sandiaD", "max_min": 120,
 "command": "curl -s https://dl.openfoam.com/add-debian-repo.sh | bash\napt-get install -y openfoam2606-default\n. /usr/lib/openfoam/openfoam2606/etc/bashrc\ncp -r \"$WM_PROJECT_DIR/tutorials/combustion/reactingFoam/RAS/SandiaD_LTS\" \"$GTU_OUT/SandiaD_LTS\"\ncd \"$GTU_OUT/SandiaD_LTS\" && ./Allrun"}
```

Не проверено: кладёт ли `openfoam2606-default` учебники в `$WM_PROJECT_DIR/tutorials`. Параллельно — `decomposePar`
и `runParallel <решатель>` (как в `Allrun` учебника `reverseBurner`). Потолок 300 долл., запись только под префиксом
`ALPHA_O_S3`, подпапка `gtu/`; ключей и секретов в файлах нет.

## 9. Типовые ошибки

- **NOx = 0 на образце Sandia.** `SandiaD_LTS/chemkin/grimech30.dat` в обеих ветвях — урезанный GRI без азота:
  36 веществ, 219 реакций (проверено двумя способами: подсчёт в файле и вывод chemkinToFoam). Для NO — полный GRI.
- **apt-пакет вместо выпуска.** `apt-get install openfoam` на Ubuntu 24.04 даёт 1912 (2019 г.), не v2606.
- **Пустые учебники.** `openfoam-examples` из Ubuntu кладёт учебники в `/usr/share/doc/…`, а в урезанном образе
  контейнера там остаются пустые папки (проверено). Учебники брать из исходников выпуска.
- **Словари другой ветви.** Перенос файлов v2606 ↔ v14 без переименования по таблице раздела 3 (combustion ↔ reaction,
  thermophysical ↔ physical, turbulence ↔ momentumTransport, fvOptions ↔ fvModels, `atomization` ↔ `atomisation`).
- **Версия EDC по умолчанию.** Справка `EDC.H` в v2606 называет умолчание «v2015» (в v14 — v2005): `version` писать явно.
- **Опечатка в образце.** v2606 `SandiaD_LTS/system/sampleDict`, набор `Radial_15`: `end (0.024 0.108);` — две
  компоненты вместо трёх; исправить на `(0.024 0 0.108)` до выборки.
- **Химия с холодного старта.** Без предварительного прогона `chemistry off` поле течения не развито — так образец не делает.
- **Состав из XiFoam / глобальной реакции.** Ни модель Xi, ни одноступенчатая реакция н-гептана не дают NOx и CO.

## 10. Проверка результата

- Журнал решателя кончается `End`; невязки и `min/max(T)` выписаны; T не выходит за физические пределы.
- Баланс массы вход/выход и баланс энергии: средняя T выхода против `T4_C`/`T41` того же режима в results-JSON.
- Две сетки, разница итогов записана; эталон Sandia D пройден на той же версии и моделях.
- NOx и CO сведены вторым способом (Cantera); способ приведения к 15 % O2 записан рядом с числом.
- В results-JSON: версия, механизм (число веществ и реакций), модели, сетки, ссылки на входные ключи; непроверенное
  помечено, пробелы — в `data/known_gaps.json`.

## 11. Скрипты трека

`routes/gtu/README.md`; `gtu_cycle.py` (`cycle`, `cycle_r`, `recup`) → `gtu_cycle_results.json`; `gtu_explore.py`
(`comp_stage`, 249 схем) → `gtu_explore_results.json`; `gtu_dropin.py` (пакеты Taurus 60 и SGT-400) →
`gtu_dropin_results.json`; `platform/aws/runner.py`; навыки `gtu-cycle`, `gtu-dropin`, `gtu-fleet`, `su2-turbo`.

## Источники (сверено 14.11.2026)

- Исходники OpenFOAM 14 (клон, коммит 64cfdf54 от 11.11.2026): `bin/reactingFoam`, `bin/XiFoam`,
  `bin/chtMultiRegionFoam`, учебники `multicomponentFluid/*`, `multiRegion/CHT/reverseBurner`, `EDC.H`, `kEpsilon.C`,
  `surfaceFieldValue.H` — https://github.com/OpenFOAM/OpenFOAM-14
- Исходники openfoam.com v2606 (зеркало, коммит 6e259574 от 28.07.2026, `META-INFO/api-info`: api=2606): решатели,
  `createFields.H`, учебники, `EDC.H`, `surfaceFieldValue.C` — https://gitlab.com/openfoam/core/openfoam
- Выпуск 14 и пакеты Ubuntu (по выдаче поиска, страницы из контейнера не открылись) — https://openfoam.org/release/14/ ,
  https://openfoam.org/news/v14-patch/ , https://openfoam.org/download/14-ubuntu/
- Выпуск v2606 (по выдаче поиска) — https://www.openfoam.com/news/main-news/openfoam-v2606 , https://www.openfoam.com/download
- Установка openfoam.com, README репозитория (прочитан) — https://sourceforge.net/projects/openfoam/files/repos/deb/
- GRI-Mech 3.0 (сайт не открылся; README — по выдаче поиска) — http://combustion.berkeley.edu/gri-mech/version30/text30.html ,
  https://ansyshelp.ansys.com/public/views/secured/corp/v251/en/chemkin_gs/i23788.html
- Sandia Flame D (по выдаче поиска) — https://kbwiki.ercoftac.org/w/index.php/Description_AC2-09 ,
  https://www.kbwiki.ercoftac.org/w/index.php/Test_Data_AC2-09 ; архив TNF (не открывался) —
  https://tnfworkshop.org/data-archives/pilotedjet/ch4-air/
- C1 = 1,6 для круглой струи — Reinhardt, Kleiser, PAMM 2014, doi:10.1002/pamm.201410310 (по выдаче поиска)
- Нормы «at 15 percent O2» — 40 CFR 60 подраздел KKKK, таблица 1 —
  https://www.ecfr.gov/current/title-40/chapter-I/subchapter-C/part-60/subpart-KKKK
