---
name: gtu-fleet
description: Паспорта и посадочные размеры установок-образцов трека «Простая литая ГТУ-заменитель» — Solar (Saturn 20 … Titan 130) и Siemens (SGT-100 … SGT-500) до 20 МВт — мощность ISO генератора и привода, КПД, расход и температура выхлопа, степень сжатия, ступени, валы, частота выходного вала, габариты и масса пакета и двигателя, присоединительные размеры; каждое число со ссылкой и цитатой, сверка двумя способами, файл routes/gtu/fleet/gtu_fleet_<модель>.json — вход gtu_family.py. Используй, когда собираешь или уточняешь паспорт установки, разбираешь расхождение источников или готовишь вход перебора семейства. Триггеры — «паспорт установки», «посадочные размеры», «линейка Solar», «линейка Siemens», «Turbomachinery Package Specification», «данные Taurus», «данные SGT», «fleet».
---

# Паспорта установок-образцов: Solar и Siemens до 20 МВт

Заменитель встаёт в чужой пакет, поэтому паспорт исходной установки — его техзадание: мощность,
частота и фланец выходного вала, расход через вход и выхлоп пакета, габарит и посадка. Навык — как
собрать паспорт одной модели с проверяемым по ссылке каждым числом в формате для `gtu_family.py`.

## 1. Что где лежит

| что | где |
|---|---|
| файл модели | `routes/gtu/fleet/gtu_fleet_<модель>.json` (раздел 5) |
| предварительные паспорта | `routes/gtu/fleet/gtu_fleet_solar_prelim.json` (18 записей), `gtu_fleet_siemens_prelim.json` (17): собраны по выдержкам поиска, страницы не открывались (см. их `meta` / `_meta`) |
| потребитель | `routes/gtu/gtu_family.py` — без аргументов читает все `fleet/*.json` |
| паспорта, вписанные в код | `gtu_study.py` (калибровка по SGT-400 и Taurus 60), `gtu_dropin.py` (словарь `SLOT`) |

Имя: `gtu_fleet_` + модель строчными, не буквы и не цифры → `_` (`gtu_fleet_sgt_400.json`). Префикс
обязателен: `tools/build_catalog.mjs` относит файл к треку `gas-turbine` по `^gtu_`, иначе «не определён».

## 2. Перечень

- **Solar:** Saturn 20, Centaur 40, Centaur 50, Mercury 50, Taurus 60, Taurus 65, Taurus 70, Mars 90,
  Mars 100, Titan 130. GS — генераторная установка, MD и CS — привод. Mercury 50 и Taurus 65 в
  предварительном файле — только GS (не сверено); Mercury 50 — с рекуператором, `gtu_family.py` его
  не подбирает (`RECUP`).
- **Siemens:** SGT-100 … SGT-500; PG — генератор, MD — привод. Рейтинги одной модели (SGT-100 5,1 и
  5,4 МВт; SGT-400 13 и 15 МВт) — разные записи `models[]` одного файла.
- **Вне трека:** всё свыше 20 МВт (Titan 250 и выше, SGT-600 и выше). SGT-A05 есть в предварительном
  файле, SGT-50 (2 МВт) — в брошюре Siemens, но в перечне трека они не названы — без решения владельца не заводить.

## 3. Порядок работы

**Шаг 0.** Начало сессии по CLAUDE.md. Предварительный файл — адреса для проверки, не данные; его не править.

**Шаг 1. Документы.** Одна модель за проход. Порядок источников:
1. Производитель. Solar — лист данных (код в последней строке колонтитула, например `DS40GS/1112/EO`;
   стр. 1 — General Specifications: ступени, π, расход воздуха, частоты, редуктор; стр. 2 — Performance и
   размеры). У генераторных листов серии Power Generation (`DS20PG/0113`, `DS60PG/0113`) длина и ширина
   есть только на чертеже «Enclosure Access and Maintenance Space», в текстовом слое PDF их нет: страницу
   рендерить (`pdftoppm`, ниже) и читать картинку; на том же чертеже — зона выкатки двигателя. Siemens —
   таблицы «Performance data for simple cycle power generation», «Performance data for mechanical drive»
   и «Physical dimensions and weight» на странице модели; плакат `Poster-SGT…` и брошюра
   `GT-Portfolio-Brochure` (ссылки — в HTML страниц SGT-100, -300, -400).
2. Solar Turbomachinery Package Specification (код в колонтитуле: `TPS70CSMD/715` — Taurus 70
   Compressor Set and Mechanical Drive; код генераторных TPS не проверен): ступени, частоты, масса узла
   двигателя, салазки. В TPS — запрет распространения
   за пределы организации прямого заказчика и пометка «Caterpillar Confidential Green». Репозиторий
   открытый: число, которого нет в листах данных и на страницах, из TPS не вносить до решения
   владельца — строка в `gaps` («есть в TPS…, ждёт решения»); текст и таблицы TPS не копировать.
3. Третьи стороны (ETN Global, Forecast International, дилеры, копии на scribd, manualmachine,
   DirectIndustry, archive.org) — только второй способ сверки, в `doc` — «третья сторона».

Доступ из облачного контейнера (проверено 14.11.2026): `solarturbines.com` — curl 403, WebFetch —
хост вне разрешённой сети; те же PDF лежат на CDN `s7d2.scene7.com/is/content/Caterpillar/<id>`,
искать — WebSearch с `allowed_domains: ["s7d2.scene7.com"]`, качать — curl. На `siemens-energy.com`
страницы `…/products-services/product/sgt-100.html`, `sgt-300.html`, `sgt-400.html` отдают таблицы;
`sgt-200.html` и `sgt-500.html` — 404, по SGT-500 есть только `…/service/sgt-500-gas-turbine.html`.
Брошюра `GT-Portfolio-Brochure` (78 МБ, «© Siemens Energy 2025», PDF от 21.01.2026) качается curl; в её
текстовом слое — сводная таблица PG (SGT-100, -300, -400; SGT-200 и SGT-500 в ней нет), таблицы MD текстом
не извлекаются. SGT-500 — только service-страница (блок «Specifications»: PG 19,1 МВт, КПД MD и PG, HR,
частота СТ 3600, π 13, выхлоп 97,9 кг/с и 376 °C, ступени; габаритов и массы нет). SGT-200 — ни страницы, ни
строки в брошюре: что не нашлось у производителя — `null` и `gaps`, третья сторона — только сверка.
etn.global, pdf.directindustry.com, scribd, manualmachine, archive.org, soar.hk, iteh, nist.gov — не
открываются: хост разрешают в настройках сети среды, обходить запрет нельзя.

```bash
D=$(mktemp -d); cd "$D"      # скачанное — во временную папку, не в репозиторий
curl -sL -A "Mozilla/5.0" -o doc.pdf "https://s7d2.scene7.com/is/content/Caterpillar/C10550201"
pdfinfo doc.pdf | grep -E "^(Title|Pages)"   # Title бывает чужим или пустым; документ — по колонтитулу
pdftotext -layout doc.pdf doc.txt    # -layout сохраняет строки таблиц
grep -n -E "Output Power|Heat Rate|Exhaust|Speed|Stage|Ratio|Airflow|Rotation|Length|Width|Height|Weight" doc.txt
pdftoppm -f 2 -l 2 -r 100 -png doc.pdf p   # p-2.png: чертёж с длиной и шириной пакета, читать глазами
curl -sL -A "Mozilla/5.0" -o p.html "https://www.siemens-energy.com/global/en/home/products-services/product/sgt-300.html"
python3 -I -c "import re,html,sys;t=open(sys.argv[1],encoding='utf-8',errors='ignore').read();t=re.sub(r'<(script|style).*?</\1>','',t,flags=re.S);print('\n'.join(l.strip() for l in html.unescape(re.sub(r'<[^>]+>','\n',t)).split('\n') if l.strip()))" p.html \
  | grep -n -A2 -E "^(Power output|Gross efficiency|Efficiency|Heat rate|Turbine speed|Drive shaft speed|Pressure ratio|Exhaust|Approx. weight|Length|Width|Height|.*exclude)"
```

**Шаг 2. Запись.** Число пишется, только если найдено в скачанном тексте: `sources[поле] =
{url, doc, quote}`; `doc` — код документа и таблица («DS40GS/1112, Performance»; с чертежа — «DS60PG/0113,
чертёж Enclosure Access»), `quote` — дословно, с единицами источника. Пересчитанное — `derived[поле]`:
формула и цитата исходного числа. Выдержка
поиска — не источник; оценки в паспорт не пишутся; не нашлось — `null` и строка в `gaps`.
Единицы: 1 lb = 0,45359237 кг; 1 ft = 0,3048 м; 1 hp = 550 ft·lbf/s = 0,74569987 кВт; 1 Btu(IT) =
1,05505585 кДж; °C = (°F − 32)/1,8; кг/ч ÷ 3600. КПД: η = 3600 / HR [кДж/(кВт·ч)], по низшей теплоте.

**Шаг 3. Сверка двумя способами** — мощность, КПД, расход и температура выхлопа, степень сжатия,
частота вала, габариты и масса пакета. Итог — в `check[поле]`.
1. Второй документ другого типа или редакции: лист данных ↔ страница ↔ плакат ↔ справочник.
2. Внутренняя согласованность: η = 3600/HR с точностью 0,15 п. п.; пары единиц в строке (кг↔lb,
   м↔ft, °C↔°F, кВт↔hp); выхлоп больше воздуха на долю топлива (`FAR = 0.017` в `gtu_family.py`);
   энергобаланс q = G·1,10·(T − 15) / (P·(100/η − 1)), 1,10 кДж/(кг·К) — принятая средняя теплоёмкость
   выхлопа (Cantera, gri30: продукты метана в сухом воздухе при доле топлива 0,017, 15→500 °C — 1,085;
   влажность ISO её немного поднимает). На 34 записях предварительных файлов q = 0,925…1,029 (GS/PG в среднем ниже:
   потери генератора и редуктора уходят не в выхлоп); вне 0,88…1,06 — ошибка в числах.

Не совпало — `conflicts[]`: `{field, values, source, taken, why}`; без второго подтверждения
`taken: null`. Сомнительное на листе в файле не «чинить» — только `conflicts`.

**Шаг 4.** Проверка — раздел 7. **Шаг 5.** По CLAUDE.md: ветка; файл модели — первым коммитом;
каталог (`node tools/build_catalog.mjs`) — в чистой рабочей копии HEAD (`git worktree add --detach`)
вторым коммитом; проверки перед PR; PR ссылается на issue трека (#59).

## 4. Расхождения источников: ISO, клеммы, вал

- **Условия.** Лист Solar: «Nominal Rating – per ISO», 15 °C, уровень моря, без потерь на входе и
  выхлопе, влажность 60 %, природный газ LHV = 35 МДж/нм³, без потерь на вспомогательные; на листе MD
  вместо последнего — «Optimum power turbine speed» и «AC-driven accessories» (рейтинг — при оптимальной,
  а не максимальной частоте СТ). Это условия ISO 3977-2 (101,3 кПа, 15 °C, 60 %; ГОСТ Р 52200-2004 —
  аутентичный текст ISO 3977-2:1997) — по выдержкам поиска, тексты стандартов из контейнера не открылись.
  На страницах Siemens условия не указаны: в `basis` — «не указаны; ISO — допущение». Площадочные
  рейтинги, старые редакции и данные третьих сторон с ISO не смешивать.
- **Клеммы и вал.** GS/PG — мощность на клеммах генератора (Solar: «Engine efficiency … (measured at
  generator terminals)»; Siemens PG: «Gross efficiency», MD — просто «Efficiency»); MD/CS — на валу
  силовой турбины. В файле —
  как в источнике, с `power_at`; в вал переводит `gtu_family.py` (`ETA_GEN = 0.972 * 0.985` —
  принятая оценка в коде). GS и MD одной модели не сравнивать как одно число.
- **Частота вала.** Solar одновальный GS — «Turbine … Max. Speed» на входе эпициклического редуктора
  на 1500/1800 об/мин (на DS60PG/0113 частоты турбины нет вовсе — `null`); Solar MD — «Power Turbine … Max. Speed»; Siemens PG — «Turbine speed» (у
  двухвального SGT-400 PG это СТ перед редуктором генератора); Siemens MD — «Turbine speed» (SGT-300,
  SGT-400) или «Drive shaft speed» (SGT-100), три числа без подписи (SGT-400: «4,750 - 9,500 – 9,975 rpm» =
  50 – 100 – 105 % среднего числа; толкование «мин – ном – макс» — наше). В `out_rpm` — одно число, в
  `out_rpm_at` — что это, диапазон — строкой в `out_rpm_range`, частота генератора или компрессора —
  `driven_rpm`. `gtu_family.py` берёт из строки первое число от 1000 до 40000: «4,750 - 9,500 –
  9,975» дало бы 4750.
- **Пакет.** Solar: «Typical Weight» (листы Oil & Gas) или «Package Weight» (листы Power Generation);
  с кожухом или без и с маслом или без — на листе обычно не сказано, тогда так и писать; GS — с
  генератором, MD — салазки двигателя. Siemens «Approx. weight»; сноска страницы: «Dimensions exclude
  inlet filter housing and exhaust stack. For power generation, the AC generator is included. For
  mechanical drive, the driven equipment is excluded.» Основание — в `pkg_basis`; числа разных оснований
  и редакций не усреднять (Taurus 60 GS: DS60PG/0113 — 39 100 кг, предварительный файл — 37 920).

## 5. Формат файла

```json
{"model": "Taurus 70", "maker": "Solar Turbines", "collected": "2026-11-14", "basis": "ISO: 15 °C, …", "gaps": ["…"],
 "models": [{"variant": "MD", "rating": "…", "power_kW": 0, "power_at": "вал силовой турбины",
   "heat_rate_kJ_kWh": 0, "eff_pct": 0, "exh_flow_kg_s": 0, "exh_temp_C": 0, "inlet_flow_kg_s": 0, "PR": 0,
   "comp_stages": 0, "turb_stages": "ТГГ 2 + СТ 2", "shafts": "two", "out_rpm": 0, "out_rpm_at": "макс. частота СТ",
   "out_rpm_range": null, "gg_rpm": 0, "driven_rpm": null, "pkg_L_m": 0, "pkg_W_m": 0, "pkg_H_m": 0, "pkg_kg": 0,
   "pkg_basis": "салазки без кожуха, без масла", "engine_L_m": null, "engine_D_m": null, "engine_kg": null,
   "interface": {"flange": null, "rotation": null, "mounting": null, "coupling": null, "inlet": null, "exhaust": null},
   "sources": {"power_kW": {"url": "https://…", "doc": "DS…, Performance", "quote": "Output Power … kW"}},
   "derived": {"eff_pct": "3600 / heat_rate_kJ_kWh"}, "check": {"eff_pct": "энергобаланс q = 0,96"},
   "conflicts": [{"field": "…", "values": ["…", "…"], "source": "…", "taken": null, "why": "…"}]}]}
```

Поля образца обязательны (неизвестное — `null`); `rating`, `out_rpm_range`, `driven_rpm`, `check` —
по возможности. Числа — числами: `gtu_family.py` (`num`) из строки выкидывает запятые и берёт первое
слово, так что «35,00 %» стало бы 3500, а «5,1» — 51. `variant` — `GS`, `PG`, `MD` или `CS` (`GS*`/`PG*` →
генератор, прочее — привод). `interface`: фланец выходного вала (тип, центрирующий поясок, болты),
направление вращения и вид, опоры на раме, муфта, присоединение к шахтам входа и выхлопа. В открытых
листах из этого бывает только направление вращения без указания, откуда смотреть (DS60PG/0113: «Clockwise
Rotation») — писать дословно, вид — `null`. Остальное, по тексту TPS, — в установочных чертежах для
заказчика. Пока нет — `null` и `gaps`.

## 6. Типовые ошибки (найдены 14.11.2026)

- Пары единиц врут: SGT-400 PG «40.0 kg/s (86.3 lb/s)», а 86,3 lb/s = 39,1 кг/с (у MD той же
  страницы — «39.1 kg/s (86.3 lb/s)»); SGT-300 PG «58.967 kg (130,844 lb)», а 130 844 lb ≈ 59 350 кг;
  SGT-100 PG «2.7 m (9.4 ft)», а 9,4 ft = 2,87 м.
- Единицы переставлены: Taurus 60 GS (DS60PG/0113) — «Inlet Airflow: 47.6 kg/sec (21.6 lb/sec)»; 21,6 кг/с =
  47,6 lb/s. В файл — исправленное только через `conflicts` и второй источник.
- Не та единица: у MD SGT-300 и SGT-400 в скобках «Btu/kWh», а число — Btu/(hp·ч) (SGT-400: 9886 кДж/(кВт·ч)
  = 6987 Btu/(hp·ч); на плакате SGT-300 — «Btu / hph»). У SGT-100 так не выходит: MD «10,803 kJ/kWh (10,195
  Btu/kWh)» при 10 239 по пересчёту, PG 5,4 МВт «11,613 kJ/kWh ( 11,292 Btu/kWh)» при 11 007. Каждую пару — пересчётом.
- КПД не равен 3600/HR: SGT-100 MD — 33,1 % при 10 803 кДж/(кВт·ч) (это 33,3 %); SGT-100 PG 5,4 МВт —
  30,3 % при 11 613 (31,0 %); SGT-500 (service-страница) — PG 33,6 % при 10 651 (33,8 %).
- Десятичная запятая: SGT-400 PG 13 МВт — «Gross efficiency 35,00 %» (в соседнем столбце «35.6%»).
- Брошюра и страница расходятся: SGT-400 PG π — брошюра GT-Portfolio 16,9 и 18,5, страницы global и US
  14.11.2026 — 17,1 и 19,5; SGT-300 PG HR — 11 704 и 11 594. В предварительном файле 16,9 приписано
  US-странице — откуда выдержка, не установлено. Дата снятия — в `collected`.
- Плакат и страница расходятся: SGT-300 PG — плакат 7,9 МВт и «543°C (993°F)» (543 °C = 1009 °F),
  страница 8,0 МВт и «534 °C (993° F)». SGT-500: предварительный файл — 369 °C, service-страница — 376 °C.
- Воздух и выхлоп не сходятся на одном листе: Centaur 40 GS (DS40GS/1112) — воздух 18,4 кг/с, выхлоп
  68 365 кг/ч = 18,99 кг/с (+3,2 %), а топливо по мощности и HR листа — около 1,4 % воздуха (LHV метана
  50 МДж/кг — допущение). Причина не установлена; в `gtu_family.py` идёт выхлоп.
- Чужие метаданные PDF: у `TPS70CSMD/715` заголовок «Taurus 60 Generator Set»; документ — по колонтитулу.
- Подозрительное на листе: Saturn 20 MD (DS20MD/1109) — газогенератор «Max. Speed: 15,000 rpm»,
  силовая турбина «22,300 rpm», а у одновального GS (DS20PG/0113) турбина «Max. Speed: 22,300 rpm».
  Возможна перестановка — `conflicts` до второго источника.
- Модель посчитана дважды: `gtu_family.py` без аргументов читает и `*_prelim.json` — запускать со списком (раздел 7).

## 7. Проверка результата

Проверщик формата (из корня репозитория; код возврата 0 — чисто):

```bash
python3 - routes/gtu/fleet/gtu_fleet_<модель>.json <<'EOF'
import json, re, sys, os
NUM = ["power_kW", "heat_rate_kJ_kWh", "eff_pct", "exh_flow_kg_s", "exh_temp_C", "inlet_flow_kg_s", "PR", "comp_stages",
       "out_rpm", "gg_rpm", "pkg_L_m", "pkg_W_m", "pkg_H_m", "pkg_kg", "engine_L_m", "engine_D_m", "engine_kg"]
REQ = NUM + ["variant", "power_at", "turb_stages", "shafts", "out_rpm_at", "pkg_basis", "interface", "sources", "derived", "conflicts"]
bad = []; err = lambda f, m: (bad.append(m), print(f"ОШИБКА {f}: {m}"))
for f in sys.argv[1:]:
    txt = open(f, encoding="utf-8").read(); d = json.loads(txt)
    for k in ("model", "maker", "collected", "basis", "models", "gaps"):
        if k not in d: err(f, f"нет поля {k}")
    if os.path.basename(f) != "gtu_fleet_" + re.sub(r"[^a-z0-9]+", "_", d.get("model", "").lower()).strip("_") + ".json": err(f, "имя файла не gtu_fleet_<модель>.json")
    if re.search(r"/home/|/tmp/|[A-Za-z]:\\\\|@[\w-]+\.[a-z]{2,}|\+\d[\d ()-]{8,}", txt): err(f, "локальный путь, почта или телефон")
    for r in d.get("models", []):
        tag = f"{f} {r.get('variant')} {r.get('rating', '')}"
        if r.get("variant") not in ("GS", "PG", "MD", "CS"): err(tag, "variant не из GS/PG/MD/CS")
        for k in REQ:
            if k not in r: err(tag, f"нет поля {k} (неизвестное — null)")
        src, der, conf = r.get("sources") or {}, r.get("derived") or {}, {c.get("field") for c in r.get("conflicts") or []}
        isnum = lambda v: isinstance(v, (int, float)) and not isinstance(v, bool)
        for k in (k for k in NUM if r.get(k) is not None):
            if not isnum(r[k]): err(tag, f"{k} не число: {r[k]!r}")
            s = src.get(k)
            if not (k in der or (isinstance(s, dict) and str(s.get("url", "")).startswith("http") and s.get("quote"))):
                err(tag, f"{k} без ссылки с цитатой и без формулы в derived")
        P, e, G, T, hr = (r.get(k) if isnum(r.get(k)) else None for k in ("power_kW", "eff_pct", "exh_flow_kg_s", "exh_temp_C", "heat_rate_kJ_kWh"))
        if e and hr and abs(3600 / hr * 100 - e) > 0.15 and "eff_pct" not in conf:
            err(tag, f"КПД {e} не равен 3600/heat_rate = {3600 / hr * 100:.2f}; разбор — в conflicts")
        if P and e and G and T and not 0.88 <= G * 1.10 * (T - 15) / (P * (100 / e - 1)) <= 1.06:
            err(tag, f"энергобаланс выхлопа {G * 1.10 * (T - 15) / (P * (100 / e - 1)):.2f} вне 0,88–1,06")
print("итог: ошибок", len(bad)); sys.exit(1 if bad else 0)
EOF
```

Вход перебора — в копии папки: `gtu_family.py` перезаписывает `gtu_family_results.json` рядом с собой.

```bash
T=$(mktemp -d); cp -r routes/gtu "$T"/ && cd "$T"/gtu
F=$(ls fleet/gtu_fleet_*.json | grep -v _prelim)
[ -n "$F" ] && python3 gtu_family.py $F   # пустой список = чтение всех fleet/*.json, вместе с _prelim
```

Время: около 4 мин на 4 свободных ядрах (сетка циклов одна на любой набор файлов); на машине, занятой
другими расчётами, 14.11.2026 — 6,5 мин (7,5 мин процессорного времени).
Ожидается строка модели с исходными `P`, `G`, `n` и подбором схем; запись без `power_kW`, `exh_flow_kg_s`
или `exh_temp_C` пропускается молча. Затем — из корня репозитория — `node tools/validate.mjs`, а
`node tools/build_catalog.mjs --check` — в чистой рабочей копии HEAD после коммита (CLAUDE.md, урок 1).

## 8. Скрипты трека (`routes/gtu/`)

`README.md` — цель, сверка модели, принятые величины; `gtu_cycle.py` — цикл на Cantera (`cycle`, `cycle_r`,
`recup`); `gtu_explore.py` — перебор 249 схем, `comp_stage` — размер центробежного колеса; `gtu_dropin.py` —
замена в пакетах Taurus 60 и SGT-400 (`SLOT` в коде); `gtu_family.py` — семейство по `fleet/`; `gtu_study.py` —
калибровка по паспортам SGT-400 и Taurus 60 (числа в коде): уточнённый паспорт сверить с ними. Уже видно:
`gtu_study.py` берёт для SGT-400 13 МВт MD π 16,8, выхлоп 39,4 кг/с, η 36,2 %, а страница 14.11.2026 — 17,1;
39,1; 36,4 %; в `gtu_dropin.SLOT["S400"]` мощность MD (13,4 МВт) стоит вместе с температурой выхлопа PG (555 °C;
у MD — 543 °C). Расхождение — в `check` паспорта и строкой владельцу; код этим навыком не правится.

## Источники (открыты и прочитаны 14.11.2026, если не сказано иное)

- Solar на CDN: https://s7d2.scene7.com/is/content/Caterpillar/C10550201 (DS40GS/1112, Centaur 40 GS),
  …/C10550212 (DS20MD/1109, Saturn 20 MD), …/C10550179 (DS20PG/0113, Saturn 20 GS), …/C10550195 (DS60PG/0113,
  Taurus 60 GS), …/C10660338 (TPS70CSMD/715 — только структура и примечания).
- Siemens: https://www.siemens-energy.com/global/en/home/products-services/product/sgt-100.html, …/sgt-300.html,
  …/sgt-400.html, https://www.siemens-energy.com/us/en/home/products-services/product/sgt-400.html,
  https://www.siemens-energy.com/global/en/home/products-services/service/sgt-500-gas-turbine.html,
  https://assets.siemens-energy.com/dam/c71bf020-1e39-48e7-8e56-b0840000a12a/Poster-SGT-300-pdf_Original%20file.pdf,
  https://assets.siemens-energy.com/dam/4b92fab0-08bf-4ffe-8eb9-b3da0096cddd/GT-Portfolio-Brochure__HQ-pdf_Original%20file.pdf
  (сводная таблица PG). Плакаты SGT-100 и SGT-400 — ссылки в HTML страниц, не скачивались.
- `pdftotext -layout`, `pdfinfo`, `pdftoppm -f -l -r -png`: справка `pdftotext -h`, `pdftoppm -h` (poppler 24.02.0)
- Не открывались, по выдержкам поиска: ISO 3977-2 — https://standards.iteh.ai/catalog/standards/iso/0bf7c7bf-0eef-4e83-82b9-ab1d2abc2935/iso-3977-2-2023;
  ГОСТ Р 52200-2004 — https://docs.cntd.ru/document/1200028959; Btu(IT) = 1055,05585262 Дж —
  https://www.knowledgedoor.com/2/units_and_constants_handbook/british-thermal-unit_group.html. Множитель
  hp выведен из точных ft и lb (550 · 0,3048 · 0,45359237 · 9,80665 = 745,6999 Вт); hp и Btu сверены
  по парам значений Siemens (SGT-400 MD, SGT-300 MD). NIST SP 811 из контейнера не открылся.
