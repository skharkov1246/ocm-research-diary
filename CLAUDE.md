# MatterForge / OCM research diary — правила проекта

Открытый билингвальный (RU/EN, уровни simple/tech) дневник квантово-химических
экспериментов по дизайну катализаторов. Статический сайт, GitHub Pages собирает
его **из ветки `main`** — что не в `main`, того нет на сайте.

Live: https://skharkov1246.github.io/ocm-research-diary/

## Главный принцип — достоверность

Каждое число — из реального расчёта. Что нельзя проверить — помечаем явно.
Ничего не приукрашиваем; неудачи и артефакты показываем. Контент берём из
`assets/matterforge/matterforge_tabs_content.json` и results-JSON в `routes/`
и `calc/` — не выдумываем свои цифры.

## Структура

| Путь | Что это |
|---|---|
| `index.html` + `entries.js` + `app.js` | OCM-дневник (вкладка 0), записи — массив `ENTRIES` |
| `femoco.html` + `femoco_entries.js` | FeMoco/аммиак-дневник, отдельная страница |
| `matterforge.js` + `matterforge.css` | Рендер вкладок продуктов (tab 1..N) |
| `assets/matterforge/matterforge_tabs_content.json` | Единый источник контента всех продукт-вкладок (`tabs[]`, у каждой — массив `updates`) + `validated_track` (Дорожка B) |
| `assets/`, `assets/femoco/`, `assets/matterforge/` | Картинки (PNG) |
| `routes/` | Расчётные скрипты и results-JSON по продукт-трекам |
| `calc/` | Расчёты Fe-цеолит / alpha-O (CH₄→метанол), AWS-скрипты |
| `MATTERFORGE_TRACKB_RESEARCH.md` | Research-блок Дорожки B (валидация метода) |

## Как добавить этап

1. OCM-дневник: дописать объект в `ENTRIES` (`entries.js`); FeMoco: в `femoco_entries.js`;
   продукт-вкладка: дописать запись в `updates` соответствующего таба в
   `matterforge_tabs_content.json`.
2. Картинки — в `assets/…`. При правке `matterforge.js`/`matterforge.css`
   обновить cache-bust штамп `?v=YYYYMMDD…` в `index.html`.
3. Проверить: `python3 -m json.tool assets/matterforge/matterforge_tabs_content.json`
   и `node --check` для изменённых `.js`.

## Git-процесс (обязательный — работа уже терялась в незамерженных ветках)

- В `main` напрямую не коммитить. Работа — в фичевой ветке.
- **Пушить ветку перед любым завершением сессии** — облачные контейнеры
  эфемерны, незапушенное исчезает.
- Сразу после пуша — открыть PR; PR вливать в `main` как только этап готов.
  **Не оставлять ветки незамерженными** — сайт собирается только из `main`.
- В начале сессии проверить хвосты: `git ls-remote --heads origin` — если висят
  старые `claude/*`-ветки с работой, сказать об этом пользователю.
- Каждый продукт-трек имеет свой GitHub Issue («Track: …») с чеклистом стадий;
  PR ссылается на issue своего трека.

## Треки (id вкладок в tabs JSON)

`h2-electrocatalysis`, `ocm`, `ammonia` (FeMoco), `polyolefins`,
`battery-cathode`, `co2-to-fuels`, `h2o2-direct`, `cc-activation-recycling`,
`oer-scaling` (группа «К тестированию»), `ch4-to-methanol` (группа
«К тестированию», расчёты в `calc/`).
