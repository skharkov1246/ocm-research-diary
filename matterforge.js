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

  let DATA = null, TABS = [], HAS_DIARY = false;
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

  function routeFromHash() {
    const [id, s] = location.hash.replace(/^#/, "").split("/");
    if (TABS.some((t) => t.id === id)) { active = id; sub = s === "log" ? "log" : "dash"; }
    else { active = HOME; sub = "dash"; }
  }
  function syncHash() {
    const target = active === HOME ? location.pathname + location.search
      : "#" + active + (sub === "log" ? "/log" : "");
    history.replaceState(null, "", target);
  }

  function goHome() { active = HOME; sub = "dash"; renderTabs(); applyView(); syncHash(); window.scrollTo({ top: 0 }); }
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
    const home = el("button", "mf-tab mf-tab-home" + (active === HOME ? " active" : ""), "← Обзор");
    home.type = "button"; home.addEventListener("click", goHome);
    nav.appendChild(home);
    TABS.forEach((t) => nav.appendChild(tabButton(t.id, t.color, t.title)));
  }

  function applyView() {
    const view = $("#mf-view");
    const showDiary = HAS_DIARY && active === "ocm" && sub === "log";
    if (HAS_DIARY) diarySections().forEach((s) => { s.hidden = !showDiary; });
    view.hidden = false;
    if (active === HOME) renderHome(view);
    else renderProject(view, TABS.find((t) => t.id === active));
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
    grid.appendChild(card("Таймлайн quantum advantage", t.timeline || "—", "честная оценка", { isText: true }));
    grid.appendChild(card("Первые клиенты",
      (t.first_clients || []).map((c) => `<span class="mf-chip">${c}</span>`).join("") || "—", "", { isText: true }));
    if (e.catalyst_cost_share_pct != null)
      grid.appendChild(card("Доля катализатора", "≈" + e.catalyst_cost_share_pct + "%", "в стоимости передела"));
    wrap.appendChild(grid);
    if (e.caveat) {
      const cav = el("div", "mf-caveat");
      cav.appendChild(el("div", "mf-caveat-label", "⚠ Честная оговорка"));
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
      ups.slice(0, 4).forEach((u) => {
        const it = el("div", "mf-feed-item");
        it.innerHTML = `<span class="mf-feed-date">${u.date}</span>` +
          (u.tag ? `<span class="mf-feed-tag">${u.tag}</span>` : "") +
          `<span class="mf-feed-title">${u.title}</span>`;
        list.appendChild(it);
      });
      sec.appendChild(list);
    } else {
      sec.appendChild(el("p", "mf-empty", "Записи операций появятся здесь по мере работы — открой «Дневник операций»."));
    }
    return sec;
  }

  /* ---------- HOME: whole-project overview ---------- */
  function renderHome(view) {
    const meta = DATA._meta || {};
    const stages = meta.stage_model || [];
    view.className = "mf-panel mf-dash";
    view.style.removeProperty("--tab");
    view.innerHTML = "";
    const totalTam = TABS.reduce((s, t) => s + (Number(t.economics && t.economics.tam_2035_busd) || 0), 0);

    const head = el("div", "mf-head");
    head.appendChild(el("h2", "mf-title-h", "Обзор проекта"));
    head.appendChild(el("p", "mf-sub",
      "Пять направлений квантовой разработки катализаторов и наш прогресс по этапам — от литературы до пилота с партнёром. Открой направление, чтобы увидеть его дашборд и дневник операций."));
    view.appendChild(head);

    const agg = el("div", "mf-dash-agg");
    agg.innerHTML =
      `<div class="mf-agg-item"><span class="mf-agg-n">${TABS.length}</span><span class="mf-agg-l">направления</span></div>` +
      `<div class="mf-agg-item"><span class="mf-agg-n">$${totalTam} млрд</span><span class="mf-agg-l">суммарный TAM 2035</span></div>` +
      `<div class="mf-agg-item"><span class="mf-agg-n">${stages.length}</span><span class="mf-agg-l">этапов до пилота</span></div>`;
    view.appendChild(agg);

    const grid = el("div", "mf-dash-grid");
    TABS.forEach((t) => {
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
      grid.appendChild(c);
    });
    view.appendChild(grid);
    if (meta.honesty_note) view.appendChild(el("p", "mf-disclaimer", `<strong>Честно:</strong> ${meta.honesty_note}`));
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

    if (meta.honesty_note) view.appendChild(el("p", "mf-disclaimer", `<strong>Честно:</strong> ${meta.honesty_note}`));
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
        const list = el("div", "mf-feed mf-feed-full");
        ups.forEach((u) => {
          const it = el("div", "mf-feed-item");
          it.innerHTML = `<span class="mf-feed-date">${u.date}</span>` +
            (u.tag ? `<span class="mf-feed-tag">${u.tag}</span>` : "") +
            `<span class="mf-feed-title">${u.title}</span>`;
          list.appendChild(it);
        });
        view.appendChild(list);
      } else {
        view.appendChild(el("p", "mf-empty",
          "Записей операций пока нет — они появятся здесь по мере работы над направлением. (Источник: поле updates[] этой вкладки в matterforge_tabs_content.json.)"));
      }
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
