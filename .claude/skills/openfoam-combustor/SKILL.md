---
name: openfoam-combustor
description: OpenFOAM для камеры сгорания и теплообмена трека «Простая литая ГТУ-заменитель» — выбор версии (openfoam.com v2606 или openfoam.org v14) и установка, решатели reactingFoam / rhoReactingFoam / XiFoam / sprayFoam / chtMultiRegionFoam и их замены foamRun / foamMultiRun в v14, кинетика GRI-Mech 3.0 и сокращённые механизмы, диффузионное горение и горение с предварительным смешением, распыл жидкого топлива для перехода газ — дизель, сопряжённый теплообмен со стенкой жаровой трубы, проверка на открытом эталоне Sandia Flame D, оценка NOx и CO с приведением к 15 % O2, прогон на AWS через routes/gtu/platform/aws/runner.py. Используй, когда ставишь OpenFOAM, собираешь или считаешь кейс камеры сгорания, выбираешь механизм, считаешь температуру стенки или выбросы двигателя-заменителя. Триггеры — «камера сгорания», «OpenFOAM», «reactingFoam», «chtMultiRegionFoam», «Sandia Flame D», «GRI-Mech», «NOx», «CO на выхлопе», «распыл дизеля», «температура жаровой трубы».
---

# OpenFOAM: камера сгорания и теплообмен двигателя-заменителя

Платформа трека: SU2 — ступени (навык `su2-turbo`), OpenFOAM — горение и теплообмен, CalculiX — прочность
(`routes/gtu/platform/fea/`, навык `calculix-rotor`), Gmsh / Open CASCADE — геометрия. Камера получает воздух за
центробежным колесом и отдаёт газ на литую турбину без охлаждения или с охлаждаемым первым сопловым аппаратом,
поэтому от расчёта нужны три вещи: поле температуры на выходе, температура стенки жаровой трубы, NOx и CO. Цикл и
размеры ядра — навык `gtu-cycle`, интерфейсы пакета — `gtu-dropin`, облако — `aws-gtu-runs`.

Состояние на 14.11.2026: кейса камеры в репозитории нет, боевых прогонов OpenFOAM не было. Проверено в контейнере
(4 ядра, Ubuntu 24.04): apt-пакет `openfoam` 1912 ставится; учебный `counterFlowFlame2D` из исходников v2606
прошёл на нём до 0,5 с за 39 с на одном ядре (T max 2016 K; блок `functions` из controlDict удалён — на 1912 любой
function object падает, раздел 9); chemkinToFoam 1912 переводит GRI-Mech 3.0. Официальные репозитории
`dl.openfoam.org` и `dl.openfoam.com` из контейнера закрыты прокси (403), сайты openfoam.org, develop.openfoam.com,
combustion.berkeley.edu не разрешаются (их надо разрешить в сетевых настройках среды), поэтому v2606 и v14 здесь не
ставились: имена решателей, словарей и учебников сверены по исходным кодам (клоны GitHub OpenFOAM-14, OpenFOAM-13 и
зеркала gitlab.com/openfoam/core/openfoam, api=2606).

## 1. Что где лежит

| что | где |
|---|---|
| вход камеры: `PR`, `T3_C`, `T4_C`, `T41_C`, `Wf` (топливо), `W_exh` | `routes/gtu/gtu_cycle_results.json` (`SGT400`, `T60`) |
| вход по схемам замены: `PR`, `T3`, `T41`; расход в пакете — `T60.G`, `S400.G` | `routes/gtu/gtu_dropin_results.json` → `res[]` |
| ядро схемы 130: `Dc` 0,435 м, `PR` 4,25, `m` 39,8 кг/с на 4 ядра (`N`), около 10 кг/с на ядро | `routes/gtu/gtu_explore_results.json` → `res[]`, `id` 130 |
| прогоны на AWS (спот, потолок 300 долл., выгрузка в `$ALPHA_O_S3/gtu/`) | `routes/gtu/platform/aws/runner.py`, `README.md` |
| кейсы камеры | `routes/gtu/platform/combustor/` — предложение, папки ещё нет |
| образцы | учебники версии: `SandiaD_LTS`, `counterFlowFlame2D_GRI_TDAC`, `aachenBomb`, `reverseBurner` (раздел 3) |

Степень сжатия в этих трёх файлах — от 3,0 (`gtu_explore`) до 16,8 (`gtu_cycle`, `SGT400`), значит давление в камере
порядка 3–17 бар. Граничные условия берутся только из этих файлов, с указанием файла и ключа в README кейса; своих
чисел не придумывать.

## 2. Версия и установка

Две ветви с несовместимыми словарями: openfoam.com (OpenCFD), выпуск v2606, и openfoam.org (OpenFOAM
Foundation), выпуск 14 (даты 26.06.2026 и 14.07.2026 — предварительно, по выдаче поиска; номера сверены по
исходникам). **Основная — v2606**: в ней раздельные решатели, названные в техзадании (`reactingFoam`,
`rhoReactingFoam`, `XiFoam`, `sprayFoam`, `simpleSprayFoam`, `reactingParcelFoam`, `chtMultiRegionFoam` с горением в
газовой области — в `createFluidFields.H` есть `CombustionModel<rhoReactionThermo>`), учебник
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
| предварительное смешение, модель Xi | `XiFoam` (psiuReactionThermo), `XiDyMFoam` | `foamRun -solver XiFluid`; модуль переработан (148 файлов против 67 в v13), словари XiFoam один к одному не переносятся |
| распыл | `sprayFoam` (нестационарный), `simpleSprayFoam` (стационарный, SIMPLE); оба с облаком `sprayCloud` | `multicomponentFluid` + в `fvModels` `clouds { type clouds; libs ("liblagrangianParcel.so"); }`, `constant/cloudProperties` (`type sprayCloud;`) |
| облако без распада капель | `reactingParcelFoam`, `simpleReactingParcelFoam` (облако `reactingCloud1`) | `cloudProperties`: `type reactingCloud;` или `reactingMultiphaseCloud;` |
| впрыск, атомизация | `coneNozzleInjection`, `atomizationModel`, `LISAAtomization` | `coneInjection`, `atomisationModel`, `LISAAtomisation` |
| сопряжённый теплообмен | `chtMultiRegionFoam`; `constant/regionProperties`: `regions (fluid (gas) solid (solid));` | `foamMultiRun`; controlDict: `regionSolvers { gas multicomponentFluid; solid solid; }` |
| стык газ — стенка | `compressible::turbulentTemperatureRadCoupledMixed` | `coupledTemperature` |
| модель горения | `combustionProperties`: `combustionModel EDC;` + `EDCCoeffs {version v2005;}` | `reactionProperties`: `reactionModel EDC;` + `EDC {version v2005;}` |
| смесь | `thermophysicalProperties`: `mixture reactingMixture;`, `inertSpecie N2;` | `physicalProperties`: `mixture multicomponentMixture;`, `defaultSpecie N2;` |
| турбулентность | `turbulenceProperties`: `RAS { RASModel kEpsilon; }` (ключ `model` тоже читается) | `momentumTransport`: `RAS { model kEpsilon; }` (ключ `RASModel` тоже читается) |
| TDAC (DAC + ISAT) | `chemistryProperties`: `chemistryType { solver ode; method TDAC; }`, `reduction { method DAC; }`, `tabulation { method ISAT; }`, `odeCoeffs` | `#includeEtc "caseDicts/solvers/chemistry/TDAC/chemistryProperties.cfg"` (там `reduction`, `tabulation`), без `chemistryType`; `ode` |
| мольные доли | function object `psiReactionThermoMoleFractions` / `rhoReactionThermoMoleFractions` | `foamPostProcess -solver multicomponentFluid -func moleFractions` |
| доля смешения | function object `BilgerMixtureFraction` (образец `counterFlowFlame2D/system/FOBilgerMixtureFraction`) | готового нет (по именам файлов исходников) — считать по составу |
| прочее | `fvOptions`, `postProcess` | `fvModels` (+ `fvConstraints`), `foamPostProcess` |
| 0D кинетика, перевод CHEMKIN | `chemFoam`, `chemkinToFoam` | `chemFoam`, `chemkinToFoam` |

В v14 `reactingFoam`, `XiFoam`, `chtMultiRegionFoam` остались сценариями в `bin/`, которые печатают замену и
запускают `foamRun`/`foamMultiRun`; `sprayFoam` и `rhoReactingFoam` в `bin/` v14 нет. Модели горения в обеих
ветвях: `laminar`, `EDC`, `PaSR` (`Cmix`), `infinitelyFastChemistry`, `diffusion`, `FSD`, `singleStep…`, `zone…`;
только в v2606 — `diffusionMulticomponent`, `eddyDissipationDiffusionModel`. Методы сокращения в обеих: `DAC`, `DRG`,
`DRGEP`, `EFA`, `PFA`; табуляция — `ISAT`.

## 4. Порядок работы

- **Шаг 0.** Начало сессии по CLAUDE.md: манифест, каталог, хвосты веток. Один проход — один кейс.
- **Шаг 1. Вход.** Из results-JSON (раздел 1): давление и температура за компрессором, расход воздуха, топлива,
  целевая температура на выходе. Ссылку на файл и ключ — в README кейса.
- **Шаг 2. Механизм** (раздел 5). Природный газ — полный GRI-Mech 3.0; дизель — отдельный механизм.
- **Шаг 3. Эталон.** Sandia Flame D (раздел 6) — до любого числа по своей камере, на той же версии и тех же
  моделях. Без пройденного эталона числа камеры идут с пометкой «не подтверждено эталоном».
- **Шаг 4. Геометрия и сетка.** Gmsh → `gmshToFoam`, затем `checkMesh`. Области для стенки — как в `reverseBurner`:
  v2606 `splitMeshRegions -cellZonesOnly -overwrite`, затем `decomposePar -allRegions`; v14 `createZones`, затем
  `splitMeshRegions -cellZones "gas solid"`. Не меньше двух сеток; разница итоговых величин между ними — в results-JSON.
- **Шаг 5. Горение.** Диффузионное и частично перемешанное — `reactingFoam` + `EDC` (version v2005) или `PaSR`.
  Стационар — псевдовремя LTS: `ddtSchemes { default localEuler; }`, как в `SandiaD_LTS`. Как в образце: сначала
  без химии (`chemistry off`, до 1500 итераций), затем с химией. Бедное предварительное смешение: `XiFoam` даёт
  фронт и тепловыделение, но не состав газа — NOx и CO по нему не считаются; для них — `reactingFoam` с механизмом.
  Ускорение химии — TDAC (синтаксис по версиям — таблица раздела 3).
- **Шаг 6. Распыл (переход газ — дизель).** Образец — `aachenBomb`: впрыск `coneNozzleInjection` (v14 —
  `coneInjection`), распад `ReitzDiwakar`, испарение `phaseChangeModel liquidEvaporationBoil`, теплообмен капли
  `RanzMarshall`, `atomizationModel none` (v14 — `atomisationModel`); для центробежной форсунки — `LISAAtomization`
  (v14 — `LISAAtomisation`). Жидкости в библиотеке обеих ветвей: `C10H22`, `C12H26`, `C16H34`, `C7H16`, `IDEA`
  (70 % н-декана + 30 % α-метилнафталина). Установившийся режим в v2606 — `simpleSprayFoam` (облако распыла,
  учебника нет); `reactingParcelFoam` с LTS (учебник `verticalChannelLTS`) годится только для облака без распада
  и атомизации. В v14 `multicomponentFluid` считает облака и с LTS (`verticalChannelLTS`), и стационарно
  (`verticalChannelSteady`), но там облака `reactingCloud`/`reactingMultiphaseCloud`; `sprayCloud` с LTS не проверялся.
- **Шаг 7. Стенка.** `chtMultiRegionFoam` / `foamMultiRun` с областями `gas` (горение) и `solid` (жаровая труба),
  образец `reverseBurner` (в обеих ветвях горение `PaSR`). Излучение: в v2606 `reverseBurner` оно выключено
  (`radiationModel none`), в v14 — `P1` в газе и `opaqueSolid` в стенке; для камеры брать `P1`, как в `SandiaD_LTS`.
  Свойства металла: в `routes/gtu/platform/fea/materials.json` есть только `ti64` и `in713c`, и только плотность,
  упругость, КТР, прочность — теплопроводности и теплоёмкости нет; брать из источника со ссылкой и дописывать.
- **Шаг 8. Выбросы** (раздел 7), затем второй способ — Cantera.
- **Шаг 9. Запись.** Кейс, `Allrun`, журнал и results-JSON с версией, механизмом, сетками, невязками, статусом
  эталона; чего нет — в `data/known_gaps.json`; затем `node tools/build_catalog.mjs` и проверки CLAUDE.md.

## 5. Кинетика

- **Полный GRI-Mech 3.0** — 53 вещества, 325 реакций, с химией азота (NO, NO2, N2O, HCN, NH3). Готов в формате
  OpenFOAM в учебнике `counterFlowFlame2D_GRI_TDAC`: v2606 — `constant/reactionsGRI` и `thermo.compressibleGasGRI`,
  v14 — `constant/reactionsGRI` и `constant/thermo` (проверено подсчётом; в Cantera 3.2 `gri30.yaml` — те же
  53/325). Указание разработчиков брать механизм как есть, без замены констант и выброса веществ, — предварительно,
  по выдаче поиска (сайт GRI-Mech из контейнера не открылся).
- **Диапазон по давлению** для камеры 3–17 бар не сверен: README30 с сайта Berkeley из контейнера не открылся, а
  вторичные источники расходятся. До сверки результаты выше 10 атм — с пометкой.
- **Сокращённый механизм:** не вырезать вещества руками, а включать DAC/ISAT (TDAC) на полном GRI.
- **Полный GRI в образце Sandia** (v2606): положить в `SandiaD_LTS/constant` оба файла из
  `counterFlowFlame2D_GRI_TDAC/constant` (тот же `thermoType`: hePsiThermo, sutherland, janaf) и убрать из `Allrun`
  строку chemkinToFoam; в v14 файл `thermo` положить под именем `thermo.compressibleGasGRI`. Поле `0/Ydefault` уже
  есть. Прогоном не проверялось.
- **Дизель:** в учебнике `aachenBomb` — одна глобальная реакция н-гептана (5 веществ); для CO и NOx не годится. В
  исходниках v14 есть детальный механизм н-гептана `test/chemistry/nc7h16/chemkin/chem.inp` (544 вещества — два
  подсчёта; около 2450 реакций; CO есть, NO нет), его берёт 0D-учебник `multicomponentFluid/nc7h16`. Механизма для
  C12H26 в учебниках нет — выбирать отдельно, со ссылкой (не выбран, не проверен).
- **Перевод Cantera YAML → OpenFOAM:** `yaml2ck mech.yaml --mechanism m.inp --thermo t.dat --transport tr.dat`, затем
  `chemkinToFoam m.inp t.dat <transportProperties> reactions thermo` (третий аргумент — словарь OpenFOAM, образец
  `SandiaD_LTS/chemkin/transportProperties`). Файлы yaml2ck (Cantera 3.2) chemkinToFoam 1912 не читает; нужны все
  три правки (проверено дважды: с тремя — 53/325, без любой одной — отказ): `THERMO` → `THERMO ALL` и строка
  диапазона температур в колонках `%10.3f` без пустых строк после неё; `Ar` → `AR` в `ELEM`;
  `REACTIONS CAL/MOLE MOLE` → `… MOLES`. На v2606/v14 не проверялось.

## 6. Эталон: Sandia Flame D

Пилотируемая струя (предварительно, по выдаче поиска — ERCOFTAC AC2-09; перед сравнением сверить с архивом TNF):
топливо 25 % CH4 + 75 % воздуха по объёму, сопло D = 7,2 мм, среднерасходная скорость 49,6 м/с, Re = 22 400;
пилот — кольцо 7,7/18,2 мм, 11,4 м/с, бедная смесь при 1880 K; местного срыва почти нет. Измерены T, доля
смешения, N2, O2, H2O, H2, CH4, CO, CO2, OH и NO (OH и NO — лазерная флуоресценция). Первоисточник — Barlow, Frank,
Proc. Combust. Inst. 27 (1998) 1087–1095. Учебник `SandiaD_LTS` повторяет геометрию (проверено по `blockMeshDict`:
радиус сопла 3,6 мм, пилот 3,85–9,1 мм; v2606 `scale 0.001`, v14 `units [mm]`) и считает RANS `kEpsilon` +
EDC v2005 + P1 + TDAC.

1. Прогнать учебник как есть: `Allrun` (v2606: chemkinToFoam → blockMesh → setFields → без химии до 1500 → с химией
   до 5000; v14 — то же без setFields, решатель `foamRun`).
2. Сравнить с данными TNF осевую T и долю смешения и радиальные профили на z/D = 7,5; 15; 30; 45; 60; 75 (станции
   `system/sampleDict` в v2606, в v14 файла нет). В `sampleDict` нет NO и доли смешения: NO дописать в `fields`, долю
   смешения в v2606 — function object `BilgerMixtureFraction`. Критерий согласия записать до прогона.
3. Повторить с полным GRI (раздел 5) и сравнить CO и NO — только после этого считать выбросы своей камеры.
4. Чувствительность к C1: по умолчанию 1,44 (проверено в исходниках обеих ветвей); для круглой струи в литературе
   берут 1,6 (предварительно, по выдаче поиска — раздел «Источники»). Где задать: v2606 `turbulenceProperties` →
   `RAS { kEpsilonCoeffs { C1 1.6; } }`; v14 `momentumTransport` → `RAS { kEpsilon { C1 1.6; } }` (читается и
   `kEpsilonCoeffs`). Прогон с 1,6 в этом проекте не делался.

Время прогона эталона на 4 ядрах не замерено. Местный прогон — не дольше 20 минут; длиннее — на AWS (раздел 8).

## 7. NOx и CO

1. Мольные доли (поля `X_<вещество>`). v2606: в `controlDict` → `functions` запись
   `{ type psiReactionThermoMoleFractions; libs (fieldFunctionObjects); }` (для `rhoReactingFoam` и
   `chtMultiRegionFoam` — `rhoReactionThermoMoleFractions`); справка `moleFractions.H`: утилитой `postProcess` не
   работает, `postProcess -func moleFractions` отвечает «Cannot find functionObject file moleFractions» (проверено на
   1912). v14: `foamPostProcess -solver multicomponentFluid -func moleFractions` (объект ищет термодинамику решателя,
   без `-solver` её нет).
2. Средние по выходу с весом расхода: `surfaceFieldValue` на патче выхода, `weightField phi;`. v2606 — операция
   `weightedAverage`; v14 — `operation average;` — с `weightField` это Σ(w·X)/Σw (`surfaceFieldValueTemplates.C`).
   Оба способа сверены по исходникам, прогоном не проверены.
3. Сухой газ: X_сух = X / (1 − X_H2O). NOx = NO + NO2 в ppm сухого газа.
4. Приведение к 15 % O2: C₁₅ = C_сух · (20,9 − 15) / (20,9 − O2,сух[%]). Нормы 40 CFR 60 подраздел KKKK даны
   «at 15 percent O2»; дословно формулу в тексте правила не сверял.
5. Второй способ: Cantera 3.2 (`ct.Solution('gri30.yaml')`, `ct.IdealGasConstPressureReactor` в `ct.ReactorNet`) при
   тех же p, T, коэффициенте избытка воздуха и времени пребывания из CFD. Расхождение — остановка и разбор.
6. Отдельной модели NOx (постобработки по Зельдовичу) в исходниках обеих ветвей нет (поиск по именам файлов
   полного дерева): NO берётся только из механизма.

## 8. Прогон на AWS

Порядок запуска, бюджет и гашение машин — навык `aws-gtu-runs`. Ключи AWS в этом контейнере не действуют —
реальных вызовов не делать, проверять задачу офлайн (`runner.load_jobs`). Поле `"apt": ["openfoam"]` поставит на
Ubuntu 24.04 версию 1912 — нужную версию ставить в `command` (команда идёт от root, из корня клона). Пример прошёл
`load_jobs`; боевого запуска не было:

```json
{"name": "of-sandiaD", "max_min": 120,
 "command": "curl -s https://dl.openfoam.com/add-debian-repo.sh | bash\napt-get update\napt-get install -y openfoam2606-default\n. /usr/lib/openfoam/openfoam2606/etc/bashrc\ncp -r \"$WM_PROJECT_DIR/tutorials/combustion/reactingFoam/RAS/SandiaD_LTS\" \"$GTU_OUT/SandiaD_LTS\"\ncd \"$GTU_OUT/SandiaD_LTS\" && ./Allrun"}
```

Не проверено: кладёт ли `openfoam2606-default` учебники в `$WM_PROJECT_DIR/tutorials`. Параллельно —
`decomposePar` (для нескольких областей `-allRegions`) и `runParallel <решатель>`, как в `Allrun` учебника
`reverseBurner`; так как команда идёт от root, перед ней
`export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1` (раздел 9). Потолок 300 долл., запись только под
префиксом `ALPHA_O_S3`, подпапка `gtu/`; ключей и секретов в файлах нет.

## 9. Типовые ошибки

- **NOx = 0 на образце Sandia.** `SandiaD_LTS/chemkin/grimech30.dat` в обеих ветвях — урезанный GRI без азота:
  36 веществ, 219 реакций (проверено подсчётом в файле и выводом chemkinToFoam). Для NO — полный GRI.
- **apt-пакет вместо выпуска.** `apt-get install openfoam` на Ubuntu 24.04 даёт 1912 (2019 г.), не v2606.
- **Function objects на apt-1912.** Любой блок `functions` в controlDict и любой `postProcess -func …` (даже
  `mag(U)`, `fieldMinMax`) падают с `FOAM FATAL IO ERROR: error in IOstream "sha1"` (проверено 14.11.2026). На
  1912 блок `functions` удалять, постобработку делать на v2606/v14.
- **Пустые учебники.** `openfoam-examples` из Ubuntu кладёт учебники в `/usr/share/doc/…`, а в образе контейнера
  `path-exclude=/usr/share/doc/*` (в `/etc/dpkg/dpkg.cfg.d`) — остаются пустые папки (проверено). Учебники брать
  из исходников выпуска.
- **mpirun от root.** OpenMPI 4.1.6 отказывается: «mpirun has detected an attempt to run as root» (проверено);
  на AWS команда идёт от root — задать `OMPI_ALLOW_RUN_AS_ROOT=1` и `OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1` или
  `--allow-run-as-root`.
- **Словари другой ветви.** Перенос файлов v2606 ↔ v14 без переименования по таблице раздела 3 (combustion ↔ reaction,
  thermophysical ↔ physical, turbulence ↔ momentumTransport, fvOptions ↔ fvModels, `atomization` ↔ `atomisation`,
  `coneNozzleInjection` ↔ `coneInjection`; v2606 читает только `odeCoeffs`, не `ode`).
- **Распыл без распада.** `reactingParcelFoam` / `simpleReactingParcelFoam` строят облако `reactingCloud1` — моделей
  распада и атомизации в нём нет; для дизеля — `sprayFoam` / `simpleSprayFoam` (v2606) или `type sprayCloud;` (v14).
- **Версия EDC по умолчанию.** В коде обеих ветвей умолчание v2005 (v2606 — `EDCs.C`), но справка `EDC.H` в v2606
  пишет «v2015» — такой версии нет (есть v1981, v1996, v2005, v2016). `version` писать явно.
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
`gtu_dropin_results.json`; `platform/aws/runner.py`; навыки `gtu-cycle`, `gtu-dropin`, `gtu-fleet`, `su2-turbo`,
`calculix-rotor`, `aws-gtu-runs`.

## Источники (сверено 14.11.2026)

- Исходники OpenFOAM 14 (клон, коммит 64cfdf54 от 11.11.2026): `bin/reactingFoam`, `bin/XiFoam`,
  `bin/chtMultiRegionFoam`, учебники `multicomponentFluid/*`, `multiRegion/CHT/reverseBurner`, `EDC.C`, `kEpsilon.C`,
  `RASModel.C`, `surfaceFieldValueTemplates.C`, `moleFractions.C`, `foamPostProcess.C`,
  `etc/caseDicts/solvers/chemistry/TDAC/`, `test/chemistry/nc7h16` — https://github.com/OpenFOAM/OpenFOAM-14
- Исходники OpenFOAM 13 (клон, коммит 18870c2 от 23.06.2026) — только для сравнения `applications/modules/XiFluid` —
  https://github.com/OpenFOAM/OpenFOAM-13
- Исходники openfoam.com v2606 (зеркало, коммит 6e259574 от 28.07.2026, `META-INFO/api-info`: api=2606): решатели,
  `createFields.H`, учебники, `EDC.H`, `EDCs.C`, `surfaceFieldValue.H`, `moleFractions.H`,
  `moleFractionsFunctionObjects.C`, `simpleSprayFoam.C` — https://gitlab.com/openfoam/core/openfoam
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
