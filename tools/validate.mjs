#!/usr/bin/env node
/**
 * Проверка содержимого сайта перед публикацией.
 *
 * Сайт собирается GitHub Pages прямо из main и рендерится только на клиенте:
 * любая опечатка в entries.js или в JSON вкладок даёт пустую страницу, и узнать
 * об этом можно лишь открыв сайт. Этот скрипт ловит такие поломки до мержа.
 *
 * Что проверяем:
 *   • синтаксис всех .js (node --check);
 *   • записи дневников: id, date, двуязычные stage/title/simple/tech, уникальность id;
 *   • JSON вкладок: обязательные поля таба и записей, длина stage_progress,
 *     треки flagship ссылаются на существующие вкладки;
 *   • все картинки, на которые ссылается контент, есть на диске;
 *   • ссылки на расчёты (routes/…, calc/…) ведут на существующие файлы.
 *
 *   node tools/validate.mjs            проверить
 *   node tools/validate.mjs --quiet    только ошибки
 */
import { readFileSync, existsSync, readdirSync, statSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { createContext, runInContext } from "node:vm";
import path from "node:path";

const ROOT = path.resolve(import.meta.dirname, "..");
const quiet = process.argv.includes("--quiet");
const errors = [];
const warns = [];
const err = (m) => errors.push(m);
const warn = (m) => warns.push(m);
const ok = (m) => { if (!quiet) console.log(`  ✓ ${m}`); };

const rel = (p) => path.relative(ROOT, p);
const listFiles = (dir, ext) => {
  const out = [];
  const walk = (d) => {
    for (const name of readdirSync(d)) {
      if (name === "node_modules" || name === ".git") continue;
      const full = path.join(d, name);
      if (statSync(full).isDirectory()) walk(full);
      else if (!ext || full.endsWith(ext)) out.push(full);
    }
  };
  walk(dir);
  return out;
};

// ---------------------------------------------------------------- 1. синтаксис JS
const jsFiles = listFiles(ROOT, ".js").filter((f) => !f.includes("/assets/"));
for (const f of jsFiles) {
  try {
    execFileSync(process.execPath, ["--check", f], { stdio: "pipe" });
  } catch (e) {
    err(`${rel(f)}: синтаксическая ошибка — ${String(e.stderr || e).split("\n")[0]}`);
  }
}
ok(`синтаксис ${jsFiles.length} файлов .js`);

// ---------------------------------------------------------------- 2. дневники
function loadEntries(file, varName) {
  const ctx = createContext({});
  runInContext(readFileSync(path.join(ROOT, file), "utf8") + `\n;globalThis.__out=${varName};`, ctx);
  return ctx.__out;
}

const BILINGUAL = ["stage", "title", "simple", "tech"];
function checkEntries(file, varName, label) {
  let entries;
  try {
    entries = loadEntries(file, varName);
  } catch (e) {
    err(`${file}: не загружается — ${e.message}`);
    return [];
  }
  if (!Array.isArray(entries)) {
    err(`${file}: ${varName} не массив`);
    return [];
  }
  const seen = new Set();
  entries.forEach((e, i) => {
    const who = `${file}[${i}]`;
    if (!e.id) err(`${who}: нет id — якорь записи станет #undefined`);
    else if (seen.has(e.id)) err(`${who}: id «${e.id}» уже используется`);
    else seen.add(e.id);
    if (!e.date) err(`${who}: нет date`);
    else if (!/^\d{4}-\d{2}(-\d{2})?$/.test(e.date)) err(`${who}: date «${e.date}» не в формате YYYY-MM-DD`);
    for (const key of BILINGUAL) {
      if (!e[key]) { err(`${who} (${e.id || "?"}): нет блока ${key}`); continue; }
      for (const lang of ["ru", "en"]) {
        if (!e[key][lang]) err(`${who} (${e.id || "?"}): ${key}.${lang} пуст — билингвальность нарушена`);
      }
    }
    for (const fig of e.figures || []) {
      if (!fig.src) { err(`${who}: у картинки нет src`); continue; }
      if (!existsSync(path.join(ROOT, fig.src))) err(`${who}: нет файла картинки ${fig.src}`);
    }
  });
  ok(`${label}: ${entries.length} записей`);
  return entries;
}

const ocm = checkEntries("entries.js", "ENTRIES", "дневник OCM");
const femoco = checkEntries("femoco_entries.js", "ENTRIES", "дневник FeMoco");

// ---------------------------------------------------------------- 3. вкладки MatterForge
const TABS_PATH = "assets/matterforge/matterforge_tabs_content.json";
let tabs = [];
try {
  const doc = JSON.parse(readFileSync(path.join(ROOT, TABS_PATH), "utf8"));
  tabs = doc.tabs || [];
  const ids = new Set();
  for (const t of tabs) {
    const who = `${TABS_PATH} → ${t.id || "без id"}`;
    if (!t.id) err(`${who}: у вкладки нет id`);
    else if (ids.has(t.id)) err(`${who}: id повторяется`);
    else ids.add(t.id);
    for (const key of ["title", "subtitle"]) {
      if (!t[key]) warn(`${who}: нет поля ${key}`);
    }
    // matterforge.js резолвит имена относительно assets/matterforge/ (IMG_DIR) и на
    // ошибку загрузки сам подставляет заглушку «иллюстрация готовится» — поэтому
    // отсутствующая иллюстрация вкладки страницу не ломает и это предупреждение
    for (const png of [t.molecules_png, t.descriptor_png].filter(Boolean)) {
      if (!existsSync(path.join(ROOT, "assets/matterforge", png))) warn(`${who}: нет иллюстрации ${png} (на сайте будет заглушка)`);
    }
    (t.updates || []).forEach((u, i) => {
      if (!u.date) err(`${who}: update[${i}] без date`);
      if (!u.text) err(`${who}: update[${i}] без text`);
    });
  }
  for (const f of doc.flagship || []) {
    for (const tr of f.tracks || []) {
      if (!ids.has(tr)) err(`${TABS_PATH}: flagship ссылается на несуществующую вкладку «${tr}»`);
    }
  }
  ok(`вкладки MatterForge: ${tabs.length}, записей в них ${tabs.reduce((s, t) => s + (t.updates || []).length, 0)}`);
} catch (e) {
  err(`${TABS_PATH}: ${e.message}`);
}

// ---------------------------------------------------------------- 4. ссылки на расчёты
const textOf = (o) => JSON.stringify(o);
const REF = /\b((?:routes|calc|payloads)\/[A-Za-z0-9_./-]+\.(?:json|py|sh|xyz|npz))/g;
const sources = [
  ["entries.js", textOf(ocm)],
  ["femoco_entries.js", textOf(femoco)],
  [TABS_PATH, existsSync(path.join(ROOT, TABS_PATH)) ? readFileSync(path.join(ROOT, TABS_PATH), "utf8") : ""],
];
// известные пробелы провенанса перечислены явно и не роняют сборку,
// всё остальное битое — ошибка
const known = new Set();
const GAPS = "data/known_gaps.json";
if (existsSync(path.join(ROOT, GAPS))) {
  for (const g of JSON.parse(readFileSync(path.join(ROOT, GAPS), "utf8")).gaps || []) known.add(g.reference);
}
let refTotal = 0, refBad = 0, refKnown = 0;
for (const [name, text] of sources) {
  for (const m of new Set((text.match(REF) || []))) {
    refTotal++;
    if (existsSync(path.join(ROOT, m))) continue;
    if (known.has(m)) { refKnown++; warn(`${name}: ${m} отсутствует — учтено в ${GAPS}`); }
    else { refBad++; err(`${name}: ссылка на несуществующий файл расчёта ${m}`); }
  }
}
ok(`ссылок на расчёты: ${refTotal}, битых ${refBad}, известных пробелов ${refKnown}`);

// ---------------------------------------------------------------- 5. картинки-сироты
const assetFiles = new Set(listFiles(path.join(ROOT, "assets"), ".png").map(rel));
const used = new Set();
for (const [, text] of sources) {
  for (const m of text.match(/assets\/[A-Za-z0-9_./-]+\.png/g) || []) used.add(m);
  for (const m of text.match(/"[A-Za-z0-9_.-]+\.png"/g) || []) used.add("assets/matterforge/" + m.replace(/"/g, ""));
}
const orphans = [...assetFiles].filter((a) => !used.has(a));
if (orphans.length) warn(`картинок не используется в контенте: ${orphans.length} (${orphans.slice(0, 3).join(", ")}…)`);

// ---------------------------------------------------------------- итог
if (!quiet) for (const w of warns) console.log(`  ⚠ ${w}`);
if (errors.length) {
  console.error("\n✗ проверка не пройдена, публиковать нельзя:");
  for (const e of errors) console.error(`   - ${e}`);
  process.exit(1);
}
console.log(`\n✓ контент прошёл проверку (предупреждений: ${warns.length})`);
