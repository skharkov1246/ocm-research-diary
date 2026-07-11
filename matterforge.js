/* ============================================================================
   MatterForge — site renderer (two-level, "grown-up" structure).
   Content: assets/matterforge/matterforge_tabs_content.json (single source).

   Navigation:
     • HOME (landing) = whole-project dashboard / overview. Reached by clicking
       the brand/logo. It is NOT a tab.
     • Top tab bar = the catalyst PROJECTS (h2, ocm, ammonia, polyolefins, battery).
     • Inside a project, a sub-nav with two views:
         – "Дашборд"  = the project's detailed, movement-first dashboard (landing).
         – "Дневник операций" = the project's operations log.
       For OCM the operations log is the existing bilingual diary (entries.js,
       rendered into .hero/#toc/#diary/#glossary); other projects use an optional
       per-tab `updates[]` array in the JSON (empty → honest placeholder).

   Nothing on screen is invented — every value comes from the JSON / entries.js.
   IIFE-scoped so helpers never collide with app.js.
   ============================================================================ */
(() => {
  "use strict";
  if ("scrollRestoration" in history) history.scrollRestoration = "manual";

  const CONTENT_URL = "assets/matterforge/matterforge_tabs_content.json";
  const IMG_DIR = "assets/matterforge/";
  const HOME = "__home__";
  const PROVEN = "proven";   // Track-B validation-benchmark tab ("Обкатанные технологии")
  const TESTING = "__testing__";
  const JOURNAL = "__journal__";   // единый архив: все операции всех треков + OCM

  /* author-supplied figure captions (not in the JSON) */
  const CAPTIONS = {
    molecules:  { _default: "Молекулы и интермедиаты", "h2-electrocatalysis": "Лестница OER-интермедиатов",
                  "ocm": "Развилка метильного радикала: C2 vs CO2" },
    descriptor: { _default: "Descriptor landscape — куда целимся",
                  "ocm": "ΔΔE‡ — дескриптор селективности" },
  };

  const $ = (sel, root = document) => root.querySelector(sel);
  const el = (tag, cls, html) => {
    const n = document.createElement(tag);
    if (cls) n.className = cls;
    if (html != null) n.innerHTML = html;
    return n;
  };
  const diarySections = () => [...document.querySelectorAll(".hero, #toc, #diary, #glossary")];
  const lang = () => document.body.dataset.lang || localStorage.getItem("ocm-lang") || "ru";
  const tt = (o) => (o ? (o[lang()] ?? o.ru ?? o.en ?? "") : "");
  const diaryEntries = () => (typeof ENTRIES !== "undefined" && Array.isArray(ENTRIES)) ? ENTRIES : [];

  /* ---- dates / "movement" helpers ---- */
  function parseDate(s) {
    const m = String(s || "").match(/^(\d{4})-(\d{2})(?:-(\d{2}))?/);
    return m ? new Date(+m[1], +m[2] - 1, +(m[3] || 1)) : null;
  }
  function relTime(s) {
    const d = parseDate(s); if (!d) return "";
    const days = Math.floor((Date.now() - d.getTime()) / 86400000);
    if (days <= 0) return "сегодня";
    if (days === 1) return "вчера";
    if (days < 30) return days + " дн. назад";
    const mo = Math.round(days / 30); return mo + " мес. назад";
  }
  /* Unified "operations" list for a project: OCM = the bilingual diary; others = JSON updates[]. */
  function updatesFor(t) {
    let items;
    if (t.id === "ocm") {
      items = diaryEntries().map((e) => ({ date: e.date, tag: tt(e.stage), title: tt(e.title) }));
    } else {
      items = (t.updates || []).map((u) => ({ date: u.date, tag: u.stage_label || "", title: u.text || u.title || "" }));
    }
    return items.filter((x) => x.date && x.title).sort((a, b) => (parseDate(b.date) - parseDate(a.date)));
  }

  let DATA = null, TABS = [], PROD = [], TEST = [], HAS_DIARY = false;
  let active = HOME;   // HOME | projectId
  let sub = "dash";    // "dash" | "log"

  async function boot() {
    const nav = $("#mf-tabs"), view = $("#mf-view");
    if (!nav || !view) return;
    HAS_DIARY = diarySections().length > 0;

    try {
      const res = await fetch(CONTENT_URL, { cache: "no-store" });
      if (!res.ok) throw new Error("HTTP " + res.status);
      DATA = await res.json();
    } catch (e) {
      console.warn("MatterForge disabled — couldn't load " + CONTENT_URL + ":", e.message);
      return;
    }
    TABS = [...(DATA.tabs || [])].sort((a, b) => (a.rank ?? 99) - (b.rank ?? 99));
    if (!TABS.length) return;
    PROD = TABS.filter((t) => (t.group || "production") !== "testing");
    TEST = TABS.filter((t) => (t.group || "production") === "testing");

    routeFromHash();

    // brand / logo → home
    const brand = $(".brand");
    if (brand) {
      brand.style.cursor = "pointer";
      brand.setAttribute("title", "На главную — обзор проекта");
      brand.addEventListener("click", goHome);
    }
    // standalone-page hooks
    const s1 = $("#mf-subtitle"); if (s1) s1.textContent = DATA._meta?.project || "";
    const s2 = $("#mf-purpose"); if (s2) s2.textContent = DATA._meta?.purpose || "";

    renderTabs();
    applyView();
  }

  let journalFilter = null;   // track id to pre-filter the journal, or null

  function routeFromHash() {
    const [id, s] = location.hash.replace(/^#/, "").split("/");
    if (id === PROVEN && DATA && DATA.validated_track) { active = PROVEN; sub = "dash"; }
    else if (id === "testing") { active = TESTING; sub = "dash"; }
    else if (id === "journal") { active = JOURNAL; sub = "dash"; journalFilter = s || null; }
    else if (TABS.some((t) => t.id === id)) { active = id; sub = s === "log" ? "log" : "dash"; }
    else { active = HOME; sub = "dash"; }
  }
  function syncHash() {
    const target = active === HOME ? location.pathname + location.search
      : active === TESTING ? "#testing"
      : active === JOURNAL ? "#journal" + (journalFilter ? "/" + journalFilter : "")
      : "#" + active + (sub === "log" ? "/log" : "");
    history.replaceState(null, "", target);
  }
  function goJournal(filter) { active = JOURNAL; sub = "dash"; journalFilter = filter || null; renderTabs(); applyView(); syncHash(); window.scrollTo({ top: 0, behavior: "smooth" }); }

  function goHome() { active = HOME; sub = "dash"; renderTabs(); applyView(); syncHash(); window.scrollTo({ top: 0 }); }
  function goTesting() { active = TESTING; sub = "dash"; renderTabs(); applyView(); syncHash(); window.scrollTo({ top: 0, behavior: "smooth" }); }
  function selectProject(id) { active = id; sub = "dash"; renderTabs(); applyView(); syncHash(); window.scrollTo({ top: 0, behavior: "smooth" }); }
  function selectSub(s) { if (sub === s) return; sub = s; applyView(); syncHash(); window.scrollTo({ top: 0, behavior: "smooth" }); }

  function tabButton(id, color, label) {
    const b = el("button", "mf-tab" + (id === active ? " active" : ""), `<span class="mf-tab-dot"></span>${label}`);
    b.style.setProperty("--tab", color);
    b.dataset.id = id;
    b.setAttribute("role", "tab");
    b.setAttribute("aria-selected", String(id === active));
    b.addEventListener("click", () => selectProject(id));
    return b;
  }
  function renderTabs() {
    const nav = $("#mf-tabs");
    nav.innerHTML = "";
    const home = el("button", "mf-tab mf-tab-home" + (active === HOME ? " active" : ""), "★ Обзор применений");
    home.type = "button"; home.addEventListener("click", goHome);
    nav.appendChild(home);
    // единый журнал вместо пер-трековых вкладок
    const onJ = active === JOURNAL || TABS.some((t) => t.id === active);
    const allUps = mergeAllUpdates();
    const jb = el("button", "mf-tab mf-tab-journal" + (onJ ? " active" : ""),
      `📓 Журнал <span class="mf-tab-count">${allUps.length}</span>`);
    jb.type = "button";
    jb.setAttribute("title", "Единый лабораторный журнал — все операции всех направлений и OCM");
    jb.addEventListener("click", () => goJournal(null));
    nav.appendChild(jb);
    if (DATA && DATA.validated_track) {
      const vt = DATA.validated_track;
      nav.appendChild(tabButton(PROVEN, vt.color || "#48bb78", tt(vt.title)));
    }
  }

  function applyView() {
    const view = $("#mf-view");
    const showDiary = HAS_DIARY && active === "ocm" && sub === "log";
    if (HAS_DIARY) diarySections().forEach((s) => { s.hidden = !showDiary; });
    view.hidden = false;
    if (active === HOME) renderHome(view);
    else if (active === PROVEN) renderValidated(view);
    else if (active === JOURNAL) renderJournal(view);
    else if (active === TESTING) renderJournal(view);
    else if (TABS.some((t) => t.id === active)) { journalFilter = active; renderJournal(view); }
    else renderHome(view);
  }

  /* ---------- shared content helpers ---------- */
  function currentStageIndex(progress) {
    const i = progress.findIndex((p) => p < 1);
    return i === -1 ? progress.length - 1 : i;
  }
  function block(title, body, kind) {
    const b = el("div", "mf-block mf-block-" + kind);
    b.appendChild(el("div", "mf-block-title", title));
    b.appendChild(el("p", "mf-block-body", body));
    return b;
  }
  function appendFigures(view, t) {
    const figs = el("div", "mf-figs");
    figs.appendChild(figure(t.molecules_png, "molecules", t.id));
    figs.appendChild(figure(t.descriptor_png, "descriptor", t.id));
    const qp = t.molecules_png && t.molecules_png.replace(/_molecules\.png$/i, "_quantum_path.png");
    if (qp && qp !== t.molecules_png) {
      const qf = el("figure", "mf-figure");
      const qi = el("img"); qi.src = IMG_DIR + qp; qi.loading = "lazy"; qi.alt = "Путь к квантовому преимуществу";
      qi.addEventListener("error", () => qf.remove());
      qf.appendChild(qi);
      qf.appendChild(el("figcaption", null, "Путь к квантовому преимуществу: предел DFT+U и рост активного пространства (иллюстративно)"));
      figs.appendChild(qf);
    }
    view.appendChild(figs);
  }
  function figure(file, kind, tabId) {
    const caption = CAPTIONS[kind][tabId] || CAPTIONS[kind]._default;
    const fig = el("figure", "mf-figure");
    const img = el("img"); img.src = IMG_DIR + file; img.loading = "lazy"; img.alt = caption;
    img.addEventListener("error", () => {
      img.remove();
      fig.prepend(el("div", "mf-fig-placeholder", `<span>иллюстрация готовится</span><small>${IMG_DIR}${file}</small>`));
    });
    fig.appendChild(img);
    fig.appendChild(el("figcaption", null, caption));
    return fig;
  }
  function roadmap(stages, progress) {
    const wrap = el("div", "mf-section");
    wrap.appendChild(el("h3", "mf-h3", "Дорожная карта до партнёрства"));
    const cur = currentStageIndex(progress);
    const road = el("div", "mf-roadmap");
    stages.forEach((name, i) => {
      const p = Math.max(0, Math.min(1, progress[i] ?? 0));
      const pct = Math.round(p * 100);
      const step = el("div", "mf-stage" + (i === cur ? " current" : "") + (p >= 1 ? " done" : ""));
      step.appendChild(el("div", "mf-stage-top",
        `<span class="mf-stage-num">${p >= 1 ? "✓" : i + 1}</span>` + (i === cur ? `<span class="mf-stage-now">сейчас</span>` : "")));
      step.appendChild(el("div", "mf-stage-name", name));
      const bar = el("div", "mf-stage-bar"); const fill = el("div", "mf-stage-fill"); fill.style.width = pct + "%";
      bar.appendChild(fill); step.appendChild(bar);
      step.appendChild(el("div", "mf-stage-pct", pct + "%"));
      road.appendChild(step);
    });
    wrap.appendChild(road);
    return wrap;
  }
  function card(label, value, sub2, { meter = null, isText = false } = {}) {
    const c = el("div", "mf-card");
    c.appendChild(el("div", "mf-card-label", label));
    const v = el("div", "mf-card-value" + (isText ? " is-text" : "")); v.innerHTML = value; c.appendChild(v);
    if (meter != null) {
      const m = el("div", "mf-meter"); const f = el("div", "mf-meter-fill");
      f.style.width = Math.max(0, Math.min(10, meter)) * 10 + "%"; m.appendChild(f); c.appendChild(m);
    }
    if (sub2) c.appendChild(el("div", "mf-card-sub", sub2));
    return c;
  }
  function economics(t) {
    const e = t.economics || {};
    const wrap = el("div", "mf-section");
    wrap.appendChild(el("h3", "mf-h3", "Экономика и контекст"));
    const grid = el("div", "mf-econ");
    grid.appendChild(card("TAM 2035", e.tam_2035_busd != null ? "$" + e.tam_2035_busd + " млрд" : "—", "целевой рынок"));
    const sv = e.saving_busd_year;
    const svDisp = sv == null ? "—" : (typeof sv === "string" && !/^\s*\d/.test(sv)) ? sv : ("$" + sv + " млрд/год");
    grid.appendChild(card("Потенциальная экономия", svDisp, "порядок величины, не гарантия"));
    grid.appendChild(card("Квантовый потенциал",
      e.quantum_potential != null ? e.quantum_potential + "<span class='mf-card-unit'> / 10</span>" : "—",
      "применимость квантовых вычислений", e.quantum_potential != null ? { meter: e.quantum_potential } : {}));
    grid.appendChild(card("Таймлайн quantum advantage", t.timeline || "—", "оценка порядка", { isText: true }));
    grid.appendChild(card("Первые клиенты",
      (t.first_clients || []).map((c) => `<span class="mf-chip">${c}</span>`).join("") || "—", "", { isText: true }));
    if (e.catalyst_cost_share_pct != null)
      grid.appendChild(card("Доля катализатора", "≈" + e.catalyst_cost_share_pct + "%", "в стоимости передела"));
    wrap.appendChild(grid);
    if (e.caveat) {
      const cav = el("div", "mf-caveat");
      cav.appendChild(el("div", "mf-caveat-label", "⚠ Ограничение метода"));
      cav.appendChild(el("p", "mf-caveat-body", e.caveat));
      wrap.appendChild(cav);
    }
    return wrap;
  }

  /* movement-first elements */
  function statusBar(t) {
    const meta = DATA._meta || {};
    const stages = meta.stage_model || [];
    const prog = t.stage_progress || [];
    const cur = currentStageIndex(prog.length ? prog : [0]);
    const overall = prog.length ? Math.round(prog.reduce((a, b) => a + (b || 0), 0) / prog.length * 100) : 0;
    const ups = updatesFor(t);
    const s = el("div", "mf-status");
    s.innerHTML =
      `<span class="mf-status-stage">Сейчас: <b>${stages[cur] || "—"}</b></span>` +
      `<span class="mf-status-prog">готовность <b>${overall}%</b></span>` +
      (ups[0] ? `<span class="mf-status-upd live">● обновлено ${relTime(ups[0].date)}</span>`
              : `<span class="mf-status-upd">журнал пуст</span>`);
    return s;
  }
  /* one operations entry — expandable (long) or a flat row (short). Interactive. */
  function opItem(u) {
    const full = u.title || "";
    const long = full.length > 96;
    const lead = long ? full.slice(0, 92).trim() + "…" : full;
    const headHtml = `<span class="mf-op-date">${u.date}</span>` +
      (u.tag ? `<span class="mf-op-tag">${u.tag}</span>` : "") +
      `<span class="mf-op-lead">${lead}</span>`;
    if (!long) {
      const row = el("div", "mf-op mf-op-flat", headHtml);
      if (u.tag) row.dataset.stage = u.tag;
      return row;
    }
    const d = el("details", "mf-op");
    if (u.tag) d.dataset.stage = u.tag;
    d.appendChild(el("summary", "mf-op-head", headHtml));
    d.appendChild(el("div", "mf-op-body", full));
    return d;
  }
  /* stage filter chips for the full log */
  function filterChips(view, ups) {
    const stages = [...new Set(ups.map((u) => u.tag).filter(Boolean))];
    if (stages.length < 2) return null;
    const bar = el("div", "mf-filter");
    const mk = (label, val) => {
      const b = el("button", "mf-chip-btn" + (val === null ? " active" : ""), label);
      b.type = "button";
      b.addEventListener("click", () => {
        bar.querySelectorAll(".mf-chip-btn").forEach((x) => x.classList.toggle("active", x === b));
        view.querySelectorAll(".mf-feed-full .mf-op").forEach((op) => {
          op.style.display = (val === null || op.dataset.stage === val) ? "" : "none";
        });
      });
      return b;
    };
    bar.appendChild(mk("Все", null));
    stages.forEach((s) => bar.appendChild(mk(s, s)));
    return bar;
  }
  function feedSection(t) {
    const ups = updatesFor(t);
    const sec = el("div", "mf-section");
    const head = el("div", "mf-section-head");
    head.appendChild(el("h3", "mf-h3", "Последние операции"));
    if (ups.length) {
      const more = el("button", "mf-link", "Все записи →"); more.type = "button";
      more.addEventListener("click", () => selectSub("log"));
      head.appendChild(more);
    }
    sec.appendChild(head);
    if (ups.length) {
      const list = el("div", "mf-feed");
      ups.slice(0, 4).forEach((u) => list.appendChild(opItem(u)));
      sec.appendChild(list);
    } else {
      sec.appendChild(el("p", "mf-empty", "Записи операций появятся здесь по мере работы — открой «Дневник операций»."));
    }
    return sec;
  }

  /* ---------- HOME: whole-project overview ---------- */
  function homeCard(t) {
    const meta = DATA._meta || {};
    const prog = t.stage_progress || [];
    const cur = currentStageIndex(prog.length ? prog : [0]);
    const ups = updatesFor(t);
    const c = el("button", "mf-dash-card"); c.type = "button"; c.style.setProperty("--tab", t.color);
    c.addEventListener("click", () => selectProject(t.id));
    const strip = prog.map((p, i) => {
      const cls = p >= 1 ? "done" : i === cur ? "cur" : "";
      return `<span class="mf-dash-seg ${cls}"><span style="width:${Math.round((p || 0) * 100)}%"></span></span>`;
    }).join("");
    const e = t.economics || {};
    c.innerHTML =
      `<div class="mf-dash-top"><span class="mf-dash-dot"></span><span class="mf-dash-name">${t.title}</span>` +
      (t.composite != null ? `<span class="mf-dash-score">${t.composite}</span>` : "") + `</div>` +
      `<div class="mf-dash-substage">Сейчас: <b>${(meta.stage_model || [])[cur] || "—"}</b></div>` +
      `<div class="mf-dash-strip">${strip}</div>` +
      `<div class="mf-dash-foot"><span class="mf-dash-upd ${ups[0] ? "live" : ""}">${ups[0] ? "● обновлено " + relTime(ups[0].date) : "журнал пуст"}</span>` +
      `<span class="mf-dash-tam">TAM $${e.tam_2035_busd ?? "—"} млрд</span></div>`;
    return c;
  }
  const READY = ["идея", "модель", "расчёт", "проверено", "к железу"];
  function readyWord(v) { return READY[Math.max(0, Math.min(4, Math.round((v / 10) * 4)))]; }
  function meterRow(label, val, hint) {
    const pct = Math.max(0, Math.min(10, val)) * 10;
    const r = el("div", "mf-flag-meter");
    r.innerHTML =
      `<div class="mf-flag-meter-top"><span class="mf-flag-meter-l">${label}</span>` +
      `<span class="mf-flag-meter-v">${hint}</span></div>` +
      `<div class="mf-meter"><div class="mf-meter-fill" style="width:${pct}%"></div></div>`;
    return r;
  }
  function flagCard(app) {
    const c = el("button", "mf-flag");
    c.type = "button";
    c.style.setProperty("--tab", app.color || "#4c9be8");
    const ups = (app.tracks || []).reduce((n, id) => {
      const t = TABS.find((x) => x.id === id); return n + (t ? updatesFor(t).length : 0);
    }, 0);
    c.innerHTML =
      `<div class="mf-flag-head"><span class="mf-flag-dot"></span><span class="mf-flag-title">${app.title}</span></div>` +
      `<div class="mf-flag-flip"><span class="mf-flag-kicker">Дерзкий ход</span>${app.flip}</div>` +
      `<div class="mf-flag-impact"><span class="mf-flag-kicker">Что меняет</span>${app.impact}</div>` +
      `<div class="mf-flag-result"><span class="mf-flag-kicker">Уже посчитано</span>${app.result}</div>`;
    const meters = el("div", "mf-flag-meters");
    meters.appendChild(meterRow("Готовность к железу", app.hw, readyWord(app.hw)));
    meters.appendChild(meterRow("Проработанность идеи", app.mat, app.mat + "/10"));
    c.appendChild(meters);
    c.appendChild(el("div", "mf-flag-foot", `${ups} записей в журнале →`));
    c.addEventListener("click", () => goJournal((app.tracks || [])[0] || null));
    return c;
  }
  function renderHome(view) {
    const meta = DATA._meta || {};
    const apps = DATA.flagship || [];
    view.className = "mf-panel mf-dash mf-home2";
    view.style.removeProperty("--tab");
    view.innerHTML = "";

    const head = el("div", "mf-head");
    head.appendChild(el("h2", "mf-title-h", "Технологии, меняющие мир"));
    head.appendChild(el("p", "mf-sub",
      "Наш необычный подход — переворачивать способ, которым меняют состояние вещества: поле вместо жара, сила вместо тепла, воздух и вода вместо ископаемого сырья. Ниже — дерзкие применения, что из этого уже посчитано и насколько они готовы к проверке на железе. Всё детальное — в едином журнале."));
    view.appendChild(head);

    const withResult = apps.length;
    const agg = el("div", "mf-dash-agg");
    agg.innerHTML =
      `<div class="mf-agg-item"><span class="mf-agg-n">${apps.length}</span><span class="mf-agg-l">дерзких применений</span></div>` +
      `<div class="mf-agg-item"><span class="mf-agg-n">${withResult}</span><span class="mf-agg-l">с реальными расчётами</span></div>` +
      `<div class="mf-agg-item"><span class="mf-agg-n">${mergeAllUpdates().length}</span><span class="mf-agg-l">операций в журнале</span></div>`;
    view.appendChild(agg);

    const grid = el("div", "mf-flag-grid");
    apps.forEach((a) => grid.appendChild(flagCard(a)));
    view.appendChild(grid);

    const jcta = el("button", "mf-testing-promo");
    jcta.type = "button";
    jcta.innerHTML =
      `<span class="mf-testing-promo-l">📓 Единый журнал</span>` +
      `<span class="mf-testing-promo-r">${mergeAllUpdates().length} операций всех направлений и OCM — хронологический архив →</span>`;
    jcta.addEventListener("click", () => goJournal(null));
    view.appendChild(jcta);

    if (DATA.validated_track) {
      const vt = DATA.validated_track;
      const cta = el("div", "mf-proven-cta");
      cta.style.setProperty("--tab", vt.color || "#48bb78");
      cta.appendChild(el("div", "mf-proven-cta-tag", tt({ ru: "Дорожка B · валидация метода", en: "Track B · validation" })));
      cta.appendChild(el("div", "mf-proven-cta-body", tt(vt.subtitle)));
      const go = el("button", "mf-link", tt({ ru: "Открыть «Обкатанные технологии» →", en: "Open \"Proven systems\" →" }));
      go.type = "button"; go.addEventListener("click", () => selectProject(PROVEN));
      cta.appendChild(go);
      view.appendChild(cta);
    }
    if (meta.honesty_note) view.appendChild(el("p", "mf-disclaimer", `<strong>Достоверность:</strong> ${meta.honesty_note}`));
  }

  /* ---------- JOURNAL: единый архив всех операций ---------- */
  function mergeAllUpdates() {
    const out = [];
    (TABS || []).forEach((t) => {
      updatesFor(t).forEach((u) => out.push({ ...u, track: t.id, trackTitle: t.title, color: t.color }));
    });
    return out.sort((a, b) => (parseDate(b.date) - parseDate(a.date)));
  }
  function journalItem(u) {
    const full = u.title || "";
    const long = full.length > 96;
    const lead = long ? full.slice(0, 92).trim() + "…" : full;
    const head =
      `<span class="mf-op-date">${u.date}</span>` +
      `<span class="mf-op-track" style="--tab:${u.color || "#888"}">${u.trackTitle}</span>` +
      (u.tag ? `<span class="mf-op-tag">${u.tag}</span>` : "") +
      `<span class="mf-op-lead">${lead}</span>`;
    if (!long) { const row = el("div", "mf-op mf-op-flat", head); row.dataset.track = u.track; return row; }
    const d = el("details", "mf-op"); d.dataset.track = u.track;
    d.appendChild(el("summary", "mf-op-head", head));
    d.appendChild(el("div", "mf-op-body", full));
    return d;
  }
  function renderJournal(view) {
    const meta = DATA._meta || {};
    view.className = "mf-panel mf-journal";
    view.style.removeProperty("--tab");
    view.innerHTML = "";
    const all = mergeAllUpdates();

    const head = el("div", "mf-head");
    head.appendChild(el("h2", "mf-title-h", "Лабораторный журнал"));
    head.appendChild(el("p", "mf-sub",
      "Единый хронологический архив всех операций по всем направлениям и OCM: расчёты, проверки, тупики и выводы. Нажми запись, чтобы развернуть; фильтруй по направлению."));
    view.appendChild(head);

    // фильтр по направлению
    const tracks = [...new Set(all.map((u) => u.track))]
      .map((id) => TABS.find((t) => t.id === id)).filter(Boolean);
    const bar = el("div", "mf-filter");
    const mkBtn = (label, val, color) => {
      const b = el("button", "mf-chip-btn" + (val === journalFilter ? " active" : ""), label);
      b.type = "button";
      if (color) b.style.setProperty("--tab", color);
      b.addEventListener("click", () => {
        journalFilter = val;
        bar.querySelectorAll(".mf-chip-btn").forEach((x) => x.classList.toggle("active", x === b));
        view.querySelectorAll(".mf-feed-full .mf-op").forEach((op) => {
          op.style.display = (val === null || op.dataset.track === val) ? "" : "none";
        });
        syncHash();
      });
      return b;
    };
    bar.appendChild(mkBtn("Все направления", null));
    tracks.forEach((t) => bar.appendChild(mkBtn(t.title, t.id, t.color)));
    view.appendChild(bar);

    view.appendChild(el("p", "mf-op-count", "Операций в журнале: " + all.length + " · нажми запись, чтобы развернуть"));
    const list = el("div", "mf-feed mf-feed-full");
    all.forEach((u) => {
      const it = journalItem(u);
      if (journalFilter && u.track !== journalFilter) it.style.display = "none";
      list.appendChild(it);
    });
    view.appendChild(list);
    if (meta.honesty_note) view.appendChild(el("p", "mf-disclaimer", `<strong>Достоверность:</strong> ${meta.honesty_note}`));
  }
  /* ---------- TESTING: early-screening group (inverse-design hypotheses) ---------- */
  function renderTesting(view) {
    const meta = DATA._meta || {};
    view.className = "mf-panel mf-dash mf-testing";
    view.style.removeProperty("--tab");
    view.innerHTML = "";
    const head = el("div", "mf-head");
    head.appendChild(el("h2", "mf-title-h", "🧪 К тестированию"));
    head.appendChild(el("p", "mf-sub",
      "Ранний скрининг — направления «просчёт от обратного». Это рабочие гипотезы: рамка и дескриптор заданы, а прогресс у каждого направления свой — от пустого журнала до реальных расчётов (см. дашборд и журнал внутри). Держим их отдельно от основных направлений, пока не подтвердим квантовый edge расчётом."));
    view.appendChild(head);
    const grid = el("div", "mf-dash-grid");
    TEST.forEach((t) => grid.appendChild(homeCard(t)));
    view.appendChild(grid);
    if (meta.honesty_note) view.appendChild(el("p", "mf-disclaimer", `<strong>Достоверность:</strong> ${meta.honesty_note}`));
  }

  /* ---------- PROJECT: sub-nav + (dashboard | operations log) ---------- */
  function renderProject(view, t) {
    view.className = "mf-panel mf-project";
    view.style.setProperty("--tab", t.color);
    view.innerHTML = "";

    const subnav = el("div", "mf-subnav");
    const mk = (key, label) => {
      const b = el("button", "mf-subtab" + (sub === key ? " active" : ""), label);
      b.type = "button"; b.addEventListener("click", () => selectSub(key));
      return b;
    };
    subnav.appendChild(mk("dash", "Дашборд"));
    subnav.appendChild(mk("log", "Дневник операций"));
    view.appendChild(subnav);

    if (sub === "log") appendLog(view, t);
    else appendDashboard(view, t);
  }

  function appendDashboard(view, t) {
    const meta = DATA._meta || {};
    const head = el("div", "mf-head");
    head.appendChild(el("div", "mf-accent-bar"));
    head.appendChild(el("h2", "mf-title-h", t.title));
    head.appendChild(el("p", "mf-sub", t.subtitle));
    head.appendChild(statusBar(t));
    head.appendChild(el("div", "mf-badge", `<span class="mf-badge-label">Нарративная роль</span>${t.narrative_role || ""}`));
    view.appendChild(head);

    view.appendChild(feedSection(t)); // movement first

    if (Array.isArray(meta.stage_model) && Array.isArray(t.stage_progress))
      view.appendChild(roadmap(meta.stage_model, t.stage_progress));

    const two = el("div", "mf-twocol");
    two.appendChild(block("Что оптимизируют квантовые вычисления", t.what_quantum_optimizes, "q"));
    two.appendChild(block("Почему классика не справляется", t.why_classical_fails, "c"));
    view.appendChild(two);

    view.appendChild(economics(t));
    appendFigures(view, t);

    if (meta.honesty_note) view.appendChild(el("p", "mf-disclaimer", `<strong>Достоверность:</strong> ${meta.honesty_note}`));
  }

  function appendLog(view, t) {
    const head = el("div", "mf-loghead");
    head.appendChild(el("h2", "mf-title-h", "Дневник операций"));
    if (t.id === "ocm") {
      head.appendChild(el("p", "mf-sub", "Хронологический журнал расчётов, проверок и выводов по OCM. Полные записи — ниже."));
      view.appendChild(head);
      // the full bilingual diary (.hero/#toc/#diary/#glossary) is revealed below #mf-view by applyView()
    } else {
      view.appendChild(head);
      const ups = updatesFor(t);
      if (ups.length) {
        const chips = filterChips(view, ups);
        if (chips) view.appendChild(chips);
        view.appendChild(el("p", "mf-op-count", "Операций в журнале: " + ups.length + " · нажми запись, чтобы развернуть"));
        const list = el("div", "mf-feed mf-feed-full");
        ups.forEach((u) => list.appendChild(opItem(u)));
        view.appendChild(list);
      } else {
        view.appendChild(el("p", "mf-empty",
          "Записей операций пока нет — они появятся здесь по мере работы над направлением. (Источник: поле updates[] этой вкладки в matterforge_tabs_content.json.)"));
      }
    }
  }

  /* ---------- PROVEN: Track-B validation benchmark (bilingual) ---------- */
  function langMini() {
    const wrap = el("div", "mf-lang-mini");
    ["ru", "en"].forEach((lg) => {
      const btn = el("button", "mf-lang-b" + (lang() === lg ? " active" : ""), lg.toUpperCase());
      btn.type = "button";
      btn.addEventListener("click", () => {
        document.body.dataset.lang = lg;
        try { localStorage.setItem("ocm-lang", lg); } catch (e) {}
        renderTabs(); applyView();
      });
      wrap.appendChild(btn);
    });
    return wrap;
  }
  function systemCard(s) {
    const nm = typeof s.name === "string" ? s.name : tt(s.name);
    const card = el("div", "mf-sys mf-sys-" + (s.status || "target"));
    const head = el("div", "mf-sys-head");
    head.appendChild(el("span", "mf-sys-name", nm));
    head.appendChild(el("span", "mf-sys-badge mf-sys-badge-" + (s.status || "target"), tt(s.status_label)));
    card.appendChild(head);
    if (s.result) card.appendChild(el("div", "mf-sys-result", tt(s.result)));
    const grid = el("div", "mf-sys-grid");
    const row = (k, v) => {
      const r = el("div", "mf-sys-row");
      r.appendChild(el("div", "mf-sys-k", k));
      const vv = el("div", "mf-sys-v"); vv.innerHTML = v; r.appendChild(vv);
      grid.appendChild(r);
    };
    row(tt({ ru: "Активный центр · реакция", en: "Active site · reaction" }), tt(s.site));
    row(tt({ ru: "Известный ответ (эксп./лит.)", en: "Known answer (exp./lit.)" }), tt(s.known));
    const dois = (s.dois || []).map((d) =>
      `<a class="mf-doi" href="https://doi.org/${d}" target="_blank" rel="noopener">${d}</a>`).join(" · ");
    row("DOI", dois || "—");
    row("Active space", s.active_space || "—");
    row(tt({ ru: "Кубиты (JW)", en: "Qubits (JW)" }), s.qubits || "—");
    row(tt({ ru: "Где DFT врёт", en: "Where DFT fails" }), tt(s.why_dft));
    row(tt({ ru: "Статус", en: "Status" }), tt(s.current));
    card.appendChild(grid);
    return card;
  }
  function renderValidated(view) {
    const vt = DATA.validated_track || {};
    view.className = "mf-panel mf-proven";
    view.style.setProperty("--tab", vt.color || "#48bb78");
    view.innerHTML = "";

    const head = el("div", "mf-head");
    head.appendChild(el("div", "mf-accent-bar"));
    const titlerow = el("div", "mf-proven-titlerow");
    titlerow.appendChild(el("h2", "mf-title-h", tt(vt.title)));
    titlerow.appendChild(langMini());
    head.appendChild(titlerow);
    head.appendChild(el("p", "mf-sub", tt(vt.subtitle)));
    head.appendChild(el("div", "mf-badge",
      `<span class="mf-badge-label">${tt({ ru: "Нарративная роль", en: "Narrative role" })}</span>${tt(vt.role)}`));
    view.appendChild(head);

    const split = el("div", "mf-track-split");
    (vt.tracks || []).forEach((tr) => {
      const c = el("div", "mf-track-card mf-track-" + (tr.kind || "rd"));
      c.appendChild(el("div", "mf-track-tag", tt(tr.tag)));
      c.appendChild(el("p", "mf-track-body", tt(tr.body)));
      split.appendChild(c);
    });
    view.appendChild(split);

    const pipe = el("div", "mf-section");
    pipe.appendChild(el("h3", "mf-h3", tt({ ru: "Единый движок — стадии", en: "Shared engine — stages" })));
    const flow = el("div", "mf-pipeline");
    (vt.stages || []).forEach((st, i) => {
      const step = el("div", "mf-pipe-step" + (st.trackB ? " mf-pipe-b" : ""));
      step.appendChild(el("span", "mf-pipe-n", String(i + 1)));
      step.appendChild(el("span", "mf-pipe-t", tt(st)));
      if (st.trackB) step.appendChild(el("span", "mf-pipe-tag", tt({ ru: "только B", en: "B only" })));
      flow.appendChild(step);
    });
    pipe.appendChild(flow);
    view.appendChild(pipe);

    const b = vt.badge || {};
    const bd = el("div", "mf-badge-def");
    bd.appendChild(el("div", "mf-badge-def-chip", tt(b.label)));
    bd.appendChild(el("p", "mf-badge-means", tt(b.means)));
    const notList = (b.not && (b.not[lang()] || b.not.ru || b.not.en)) || [];
    if (notList.length) {
      const nw = el("div", "mf-badge-not");
      nw.appendChild(el("div", "mf-badge-not-h", tt({ ru: "Плашка НЕ означает:", en: "The badge does NOT mean:" })));
      const ul = el("ul", "mf-badge-list");
      notList.forEach((x) => ul.appendChild(el("li", null, x)));
      nw.appendChild(ul);
      bd.appendChild(nw);
    }
    bd.appendChild(el("div", "mf-badge-status",
      `<strong>${tt({ ru: "Статус сейчас", en: "Status now" })}:</strong> ${tt(b.status_now)}`));
    view.appendChild(bd);

    const sysSec = el("div", "mf-section");
    sysSec.appendChild(el("h3", "mf-h3", tt({ ru: "Системы Дорожки B", en: "Track B systems" })));
    (vt.systems || []).forEach((s) => sysSec.appendChild(systemCard(s)));
    view.appendChild(sysSec);

    if (vt.caveat) {
      const cav = el("div", "mf-proven-caveat");
      cav.appendChild(el("div", "mf-caveat-label", "⚠ " + tt({ ru: "Оговорка", en: "Caveat" })));
      cav.appendChild(el("p", "mf-caveat-body", tt(vt.caveat)));
      view.appendChild(cav);
    }
  }

  /* re-render on language toggle (the OCM feed/log are bilingual) */
  document.addEventListener("click", (ev) => {
    const b = ev.target.closest('.toggle[data-set="lang"]');
    if (b && DATA) { renderTabs(); applyView(); }
  });
  window.addEventListener("hashchange", () => {
    const prev = active + "/" + sub;
    routeFromHash();
    if (prev !== active + "/" + sub) { renderTabs(); applyView(); }
  });

  boot();
})();
