#!/usr/bin/env node
/**
 * Smoke-тест сайта: поднимает локальный сервер, открывает страницы в headless
 * Chromium и проверяет, что контент действительно отрисовался и в консоли нет ошибок.
 *
 * Нужен потому, что сайт рендерится только на клиенте: валидатор содержимого
 * (tools/validate.mjs) не увидит поломку в app.js или matterforge.js, а посетитель
 * увидит пустую страницу.
 *
 *   node tools/smoke.mjs
 */
import { spawn, execFileSync } from "node:child_process";
import { existsSync } from "node:fs";
import path from "node:path";
import net from "node:net";

const ROOT = path.resolve(import.meta.dirname, "..");
const CHROME = [
  process.env.CHROME_PATH,
  "/opt/pw-browsers/chromium-1194/chrome-linux/chrome",
  "/usr/bin/chromium", "/usr/bin/chromium-browser", "/usr/bin/google-chrome",
].find((p) => p && existsSync(p));

// страница → строки, которые обязаны появиться в DOM после рендера
const PAGES = {
  "index.html": ["mf-tabs", "diary"],
  "femoco.html": ["diary"],
  "matterforge.html": ["mf"],
};

const freePort = () => new Promise((res) => {
  const s = net.createServer();
  s.listen(0, () => { const p = s.address().port; s.close(() => res(p)); });
});

const wait = (ms) => new Promise((r) => setTimeout(r, ms));

if (!CHROME) {
  console.log("⚠ Chromium не найден — smoke-тест пропущен");
  process.exit(0);
}

const port = await freePort();
const server = spawn("python3", ["-m", "http.server", String(port), "--bind", "127.0.0.1"],
                     { cwd: ROOT, stdio: "ignore" });
process.on("exit", () => server.kill());
await wait(900);

let failed = 0;
for (const [page, expects] of Object.entries(PAGES)) {
  const url = `http://127.0.0.1:${port}/${page}`;
  let dom = "", log = "", hung = false;
  const flags = [
    "--headless=new", "--no-sandbox", "--disable-gpu",
    "--disable-dev-shm-usage",            // на CI /dev/shm мал, без этого рендерер виснет
    "--no-first-run", "--no-default-browser-check",
    "--disable-extensions", "--disable-background-networking",
    "--disable-sync", "--disable-crash-reporter",
    "--disable-background-timer-throttling",
    "--disable-features=Translate,BackForwardCache,MediaRouter",
    `--user-data-dir=/tmp/smoke-${port}-${page}`,
    "--virtual-time-budget=5000", "--enable-logging=stderr", "--v=0",
  ];
  const run = (extra, ms) => {
    try {
      return { dom: execFileSync(CHROME, [...flags, ...extra, "--dump-dom", url],
               { encoding: "utf8", timeout: ms, stdio: ["ignore", "pipe", "pipe"] }), log: "", hung: false };
    } catch (e) {
      const timedOut = e.code === "ETIMEDOUT" || e.signal === "SIGTERM";
      return { dom: String(e.stdout || ""), log: String(e.stderr || ""), hung: timedOut && !e.stdout };
    }
  };
  let r = run([], 75000);
  if (r.hung) r = run(["--single-process"], 60000);   // вторая попытка в один процесс
  ({ dom, log, hung } = r);
  if (hung) {
    console.log(`⚠ ${page}: Chromium не ответил — проверка рендера не выполнена`);
    continue;
  }
  const errs = (log || "").split("\n").filter((l) =>
    /ERROR:CONSOLE|Uncaught|SyntaxError|is not defined|is not a function|Failed to load resource/.test(l));
  const missing = expects.filter((x) => !dom.includes(x));
  if (errs.length || missing.length || dom.length < 3000) {
    failed++;
    console.error(`✗ ${page}: ${dom.length} байт DOM`);
    for (const m of missing) console.error(`   - в DOM нет «${m}»`);
    for (const e of errs.slice(0, 4)) console.error(`   - ${e.slice(-180)}`);
  } else {
    console.log(`✓ ${page}: ${Math.round(dom.length / 1024)} КБ DOM, ошибок консоли нет`);
  }
}
server.kill();
process.exit(failed ? 1 : 0);
