#!/usr/bin/env node
/**
 * Каталог расчётов и данных: data/catalog.json.
 *
 * В репозитории 391 файл в routes/ и 37 в calc/ одной кучей, и понять «что уже
 * посчитано, каким скриптом, попало ли это в дневник» можно только чтением всего
 * подряд. Каталог собирается сканированием и отвечает на это машинно.
 *
 * Для каждого результата: путь, трек, объём, ключи верхнего уровня, порождающий
 * скрипт (если найден), упоминается ли в контенте сайта, когда и кем менялся.
 *
 *   node tools/build_catalog.mjs           пересобрать
 *   node tools/build_catalog.mjs --check   проверить актуальность (для CI)
 */
import { readFileSync, writeFileSync, existsSync, readdirSync, statSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { createContext, runInContext } from "node:vm";
import path from "node:path";

const ROOT = path.resolve(import.meta.dirname, "..");
const OUT = path.join(ROOT, "data", "catalog.json");
const check = process.argv.includes("--check");

// трек определяется по префиксу имени файла; порядок важен — первое совпадение выигрывает
const TRACKS = [
  [/^ocm_field|^ocm_mnw|^ocm_stage|^ocm_aws|^spin_field/, "ocm"],
  [/^nh3_|^femoco|^fe4s4/, "ammonia"],
  [/^h2_|^oer/, "h2-electrocatalysis"],
  [/^polyolefins|^lao_/, "polyolefins"],
  [/^co2_to_fuels/, "co2-to-fuels"],
  [/^h2o2/, "h2o2-direct"],
  [/^cc_|^pyrolysis|^chemloop/, "cc-activation-recycling"],
  [/^fe_zeolite|^alpha_o|^feo_|^fe_ferryl/, "ch4-to-methanol"],
  [/^battery|^cathode/, "battery-cathode"],
  [/^nacn|^msa|^pom_|^urea|^acetic|^dife|^antibep/, "прочее"],
];
const trackOf = (name) => (TRACKS.find(([rx]) => rx.test(name)) || [null, "не определён"])[1];

const sh = (...args) => {
  try { return execFileSync(args[0], args.slice(1), { cwd: ROOT, encoding: "utf8", timeout: 60000 }).trim(); }
  catch { return ""; }
};

const gitMeta = (rel) => {
  const line = sh("git", "log", "-1", "--format=%ad|%an|%h", "--date=short", "--", rel);
  if (!line.includes("|")) return { last_change: null, last_author: null, last_commit: null };
  const [date, author, sha] = line.split("|");
  return { last_change: date, last_author: author, last_commit: sha };
};

const listFiles = (dir) => {
  const out = [];
  const walk = (d) => {
    if (!existsSync(d)) return;
    for (const n of readdirSync(d)) {
      const f = path.join(d, n);
      if (statSync(f).isDirectory()) walk(f); else out.push(f);
    }
  };
  walk(dir);
  return out;
};

// ------------------------------------------------------------------ контент сайта
const loadEntries = (file, varName) => {
  try {
    const ctx = createContext({});
    runInContext(readFileSync(path.join(ROOT, file), "utf8") + `\n;globalThis.__o=${varName};`, ctx);
    return ctx.__o || [];
  } catch { return []; }
};
const ocmEntries = loadEntries("entries.js", "ENTRIES");
const femocoEntries = loadEntries("femoco_entries.js", "ENTRIES");
const tabsPath = "assets/matterforge/matterforge_tabs_content.json";
const tabsRaw = existsSync(path.join(ROOT, tabsPath)) ? readFileSync(path.join(ROOT, tabsPath), "utf8") : "{}";
const contentBlob = JSON.stringify(ocmEntries) + JSON.stringify(femocoEntries) + tabsRaw;

// ------------------------------------------------------------------ наборы данных
const scriptStems = new Set(
  [...listFiles(path.join(ROOT, "routes")), ...listFiles(path.join(ROOT, "calc"))]
    .filter((f) => f.endsWith(".py"))
    .map((f) => path.basename(f, ".py")));

const datasets = [];
for (const dir of ["routes", "calc"]) {
  for (const f of listFiles(path.join(ROOT, dir))) {
    const rel = path.relative(ROOT, f);
    const ext = path.extname(f);
    if (![".json", ".npz", ".xyz"].includes(ext)) continue;
    const size = statSync(f).size;
    const base = path.basename(f, ext);
    const entry = {
      path: rel,
      track: trackOf(path.basename(f)),
      kind: base.endsWith("_results") ? "результат расчёта" : ext === ".xyz" ? "геометрия" : ext === ".npz" ? "кэш волновой функции" : "данные",
      bytes: size,
      referenced_in_site: contentBlob.includes(rel),
    };
    // порождающий скрипт: точное совпадение имени или имя без суффикса _results
    const stem = base.replace(/_results$/, "");
    for (const cand of [base, stem]) {
      if (scriptStems.has(cand)) { entry.produced_by = `${dir}/${cand}.py`; break; }
    }
    if (ext === ".json") {
      try {
        const d = JSON.parse(readFileSync(f, "utf8"));
        if (Array.isArray(d)) { entry.shape = "список"; entry.records = d.length; }
        else if (d && typeof d === "object") { entry.shape = "объект"; entry.top_keys = Object.keys(d).slice(0, 25); }
        // метаданные расчёта, если автор их записал
        const meta = ["method", "basis", "functional", "software", "date", "commit"].filter((k) => d && k in d);
        if (meta.length) entry.provenance_fields = meta;
      } catch (e) { entry.parse_error = String(e.message).slice(0, 120); }
    }
    Object.assign(entry, gitMeta(rel));
    datasets.push(entry);
  }
}
datasets.sort((a, b) => (a.track === b.track ? a.path.localeCompare(b.path) : a.track.localeCompare(b.track)));

const byTrack = {};
for (const d of datasets) byTrack[d.track] = (byTrack[d.track] || 0) + 1;
const results = datasets.filter((d) => d.kind === "результат расчёта");

const catalog = {
  schema: "matterforge.data-catalog/1",
  repo: "skharkov1246/ocm-research-diary",
  generated_by: "tools/build_catalog.mjs",
  how_to_refresh: "node tools/build_catalog.mjs",
  summary: {
    datasets: datasets.length,
    total_bytes: datasets.reduce((s, d) => s + d.bytes, 0),
    results_files: results.length,
    results_with_script: results.filter((d) => d.produced_by).length,
    results_referenced_in_site: results.filter((d) => d.referenced_in_site).length,
    results_with_provenance_fields: results.filter((d) => d.provenance_fields).length,
    by_track: byTrack,
  },
  content: {
    ocm_entries: ocmEntries.length,
    femoco_entries: femocoEntries.length,
    tabs: (JSON.parse(tabsRaw).tabs || []).length,
    updates: (JSON.parse(tabsRaw).tabs || []).reduce((s, t) => s + (t.updates || []).length, 0),
  },
  known_gaps: existsSync(path.join(ROOT, "data/known_gaps.json")) ? "data/known_gaps.json" : null,
  datasets,
};

const text = JSON.stringify(catalog, null, 2) + "\n";
if (check) {
  if (!existsSync(OUT)) { console.error("✗ нет data/catalog.json — выполните: node tools/build_catalog.mjs"); process.exit(1); }
  // Проверяем СОСТАВ и ФОРМУ каталога: не появился ли набор мимо каталога и не
  // исчез ли описанный. Всё, что меняется само по себе, из сравнения исключаем —
  // поля из git и объём файла: на событии pull_request проверяется merge-коммит
  // с актуальным main, где данные могли уже обновиться.
  const strip = (c) => c.datasets.map(({ last_change, last_author, last_commit, bytes, records, ...rest }) => rest);
  const old = JSON.parse(readFileSync(OUT, "utf8"));
  if (JSON.stringify(strip(old)) !== JSON.stringify(strip(catalog))) {
    console.error("✗ каталог устарел — выполните: node tools/build_catalog.mjs");
    process.exit(1);
  }
  console.log(`✓ каталог актуален: ${catalog.datasets.length} наборов`);
} else {
  writeFileSync(OUT, text);
  const s = catalog.summary;
  console.log(`✓ data/catalog.json: ${s.datasets} наборов (${(s.total_bytes / 1e6).toFixed(1)} МБ), ` +
    `результатов ${s.results_files}, из них со скриптом ${s.results_with_script}, ` +
    `упомянуто на сайте ${s.results_referenced_in_site}, с метаданными расчёта ${s.results_with_provenance_fields}`);
}
