/* ============================================================================
   MatterForge journal — tab renderer.
   Drives any page that provides #mf-tabs (the bar) and #mf-view (the panel).
   Content comes from assets/matterforge/matterforge_tabs_content.json — all five
   tabs share one schema, so one generic renderer fills the whole bar. Nothing on
   screen is invented; every value is read from the JSON.

   Structure:
     • tab 0 = "Обзор" dashboard (project overview + per-direction stage progress) — the landing.
     • tabs 1..N = the catalyst projects; each project keeps its own diary INSIDE its tab.
       The OCM research diary (.hero/#toc/#diary/#glossary) is embedded in the "ocm" tab.

   Wrapped in an IIFE so its helpers ($/el/…) never collide with app.js, which
   declares its own globals of the same name in the shared classic-script scope.
   ============================================================================ */
(() => {
  "use strict";

  // a deep-linked tab should open at the top, not at a restored diary scroll position
  if ("scrollRestoration" in history) history.scrollRestoration = "manual";

  const CONTENT_URL = "assets/matterforge/matterforge_tabs_content.json";
  const IMG_DIR = "assets/matterforge/";
  const DIARY_ID = "__diary__";
  const DASH = "__dashboard__"; // landing view: whole-project overview + progress

  /* figure captions are author-supplied (not in the JSON); generic default with
     a per-tab override for the H₂ OER ladder requested in the spec. */
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

  let DATA = null;
  let TABS = [];
  let HAS_DIARY = false;
  let active = null;

  async function boot() {
    const nav = $("#mf-tabs");
    const view = $("#mf-view");
    if (!nav || !view) return; // host page didn't provide the hooks
    HAS_DIARY = diarySections().length > 0;

    try {
      const res = await fetch(CONTENT_URL, { cache: "no-store" });
      if (!res.ok) throw new Error("HTTP " + res.status);
      DATA = await res.json();
    } catch (e) {
      // graceful: leave the host page fully working, just no MatterForge tabs
      console.warn("MatterForge tabs disabled — couldn't load " + CONTENT_URL + ":", e.message);
      return;
    }

    TABS = [...(DATA.tabs || [])].sort((a, b) => (a.rank ?? 99) - (b.rank ?? 99));
    if (!TABS.length) return;

    const fromHash = location.hash.replace(/^#/, "");
    active = TABS.some((t) => t.id === fromHash) ? fromHash : DASH; // landing = dashboard

    // optional header/footer hooks on a standalone page
    const sub = $("#mf-subtitle"); if (sub) sub.textContent = DATA._meta?.project || "";
    const purpose = $("#mf-purpose"); if (purpose) purpose.textContent = DATA._meta?.purpose || "";

    renderTabs();
    applyView();
  }

  function tabButton(id, color, label) {
    const b = el("button", "mf-tab" + (id === active ? " active" : ""),
      `<span class="mf-tab-dot"></span>${label}`);
    b.style.setProperty("--tab", color);
    b.dataset.id = id;
    b.setAttribute("role", "tab");
    b.setAttribute("aria-selected", String(id === active));
    b.addEventListener("click", () => selectTab(id));
    return b;
  }

  function renderTabs() {
    const nav = $("#mf-tabs");
    nav.innerHTML = "";
    nav.appendChild(tabButton(DASH, "var(--accent)", "📊 Обзор"));
    TABS.forEach((t) => nav.appendChild(tabButton(t.id, t.color, t.title)));
  }

  function selectTab(id) {
    if (active === id) return;
    active = id;
    history.replaceState(null, "", id === DASH ? location.pathname + location.search : "#" + id);
    renderTabs();
    applyView();
    window.scrollTo({ top: 0, behavior: "smooth" });
  }

  function applyView() {
    const view = $("#mf-view");
    // Dashboard is the landing; each project's diary lives inside its tab (OCM = .hero/#toc/#diary/#glossary).
    const showDiary = HAS_DIARY && active === "ocm";
    if (HAS_DIARY) diarySections().forEach((s) => { s.hidden = !showDiary; });
    view.hidden = false;
    if (active === DASH) renderDashboard(view);
    else renderPanel(view, TABS.find((t) => t.id === active));
  }

  /* current stage = first one not yet complete (else the last). */
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

  function figure(file, kind, tabId) {
    const caption = CAPTIONS[kind][tabId] || CAPTIONS[kind]._default;
    const fig = el("figure", "mf-figure");
    const img = el("img");
    img.src = IMG_DIR + file;
    img.loading = "lazy";
    img.alt = caption;
    img.addEventListener("error", () => {
      img.remove();
      fig.prepend(el("div", "mf-fig-placeholder",
        `<span>иллюстрация готовится</span><small>${IMG_DIR}${file}</small>`));
    });
    fig.appendChild(img);
    fig.appendChild(el("figcaption", null, caption));
    return fig;
  }

  function roadmap(stages, progress) {
    const wrap = el("div", "mf-section");
    wrap.appendChild(el("h3", "mf-h3", "Дорожная карта"));
    const cur = currentStageIndex(progress);
    const road = el("div", "mf-roadmap");
    stages.forEach((name, i) => {
      const p = Math.max(0, Math.min(1, progress[i] ?? 0));
      const pct = Math.round(p * 100);
      const step = el("div", "mf-stage" + (i === cur ? " current" : "") + (p >= 1 ? " done" : ""));
      step.appendChild(el("div", "mf-stage-top",
        `<span class="mf-stage-num">${p >= 1 ? "✓" : i + 1}</span>` +
        (i === cur ? `<span class="mf-stage-now">сейчас</span>` : "")));
      step.appendChild(el("div", "mf-stage-name", name));
      const bar = el("div", "mf-stage-bar");
      const fill = el("div", "mf-stage-fill");
      fill.style.width = pct + "%";
      bar.appendChild(fill);
      step.appendChild(bar);
      step.appendChild(el("div", "mf-stage-pct", pct + "%"));
      road.appendChild(step);
    });
    wrap.appendChild(road);
    return wrap;
  }

  function card(label, value, sub, { meter = null, isText = false } = {}) {
    const c = el("div", "mf-card");
    c.appendChild(el("div", "mf-card-label", label));
    const v = el("div", "mf-card-value" + (isText ? " is-text" : ""));
    v.innerHTML = value;
    c.appendChild(v);
    if (meter != null) {
      const m = el("div", "mf-meter");
      const f = el("div", "mf-meter-fill");
      f.style.width = Math.max(0, Math.min(10, meter)) * 10 + "%";
      m.appendChild(f);
      c.appendChild(m);
    }
    if (sub) c.appendChild(el("div", "mf-card-sub", sub));
    return c;
  }

  function economics(t) {
    const e = t.economics || {};
    const wrap = el("div", "mf-section");
    wrap.appendChild(el("h3", "mf-h3", "Экономика и контекст"));
    const grid = el("div", "mf-econ");
    grid.appendChild(card("TAM 2035", "$" + e.tam_2035_busd + " млрд", "целевой рынок"));
    // saving may be a plain range ("10–20") or a contextual phrase (OCM) — don't $-wrap the phrase
    const _sv = e.saving_busd_year;
    const _svDisp = (typeof _sv === "string" && !/^\s*\d/.test(_sv)) ? _sv : ("$" + _sv + " млрд/год");
    grid.appendChild(card("Потенциальная экономия", _svDisp, "порядок величины, не гарантия"));
    grid.appendChild(card("Квантовый потенциал",
      e.quantum_potential + "<span class='mf-card-unit'> / 10</span>",
      "применимость квантовых вычислений", { meter: e.quantum_potential }));
    grid.appendChild(card("Таймлайн quantum advantage", t.timeline, "честная оценка", { isText: true }));
    grid.appendChild(card("Первые клиенты",
      (t.first_clients || []).map((c) => `<span class="mf-chip">${c}</span>`).join(""), "", { isText: true }));
    if (e.catalyst_cost_share_pct != null)
      grid.appendChild(card("Доля катализатора", "≈" + e.catalyst_cost_share_pct + "%", "в стоимости передела"));
    wrap.appendChild(grid);

    // honest economics caveat — shown prominently as a warning block, never hidden (OCM, polyolefins).
    if (e.caveat) {
      const cav = el("div", "mf-caveat");
      cav.appendChild(el("div", "mf-caveat-label", "⚠ Честная оговорка"));
      cav.appendChild(el("p", "mf-caveat-body", e.caveat));
      wrap.appendChild(cav);
    }
    return wrap;
  }

  function renderDashboard(view) {
    const meta = DATA._meta || {};
    const stages = meta.stage_model || [];
    view.className = "mf-panel mf-dash";
    view.style.removeProperty("--tab");
    const totalTam = TABS.reduce((s, t) => s + (Number(t.economics && t.economics.tam_2035_busd) || 0), 0);

    const head = el("div", "mf-head");
    head.appendChild(el("h2", "mf-title-h", "Обзор проекта"));
    head.appendChild(el("p", "mf-sub",
      "Пять направлений квантовой разработки катализаторов и наш прогресс по этапам — от литературы до пилота с партнёром. Открой направление, чтобы увидеть детали и его дневник."));

    const agg = el("div", "mf-dash-agg");
    agg.innerHTML =
      `<div class="mf-agg-item"><span class="mf-agg-n">${TABS.length}</span><span class="mf-agg-l">направления</span></div>` +
      `<div class="mf-agg-item"><span class="mf-agg-n">$${totalTam} млрд</span><span class="mf-agg-l">суммарный TAM 2035</span></div>` +
      `<div class="mf-agg-item"><span class="mf-agg-n">${stages.length}</span><span class="mf-agg-l">этапов до пилота</span></div>`;

    const grid = el("div", "mf-dash-grid");
    TABS.forEach((t) => {
      const prog = t.stage_progress || [];
      const cur = currentStageIndex(prog.length ? prog : [0]);
      const c = el("button", "mf-dash-card");
      c.type = "button";
      c.style.setProperty("--tab", t.color);
      c.addEventListener("click", () => selectTab(t.id));
      const strip = prog.map((p, i) => {
        const cls = p >= 1 ? "done" : i === cur ? "cur" : "";
        return `<span class="mf-dash-seg ${cls}"><span style="width:${Math.round((p || 0) * 100)}%"></span></span>`;
      }).join("");
      const e = t.economics || {};
      c.innerHTML =
        `<div class="mf-dash-top"><span class="mf-dash-dot"></span><span class="mf-dash-name">${t.title}</span>` +
        (t.composite != null ? `<span class="mf-dash-score">${t.composite}</span>` : "") + `</div>` +
        `<div class="mf-dash-substage">Сейчас: <b>${stages[cur] || "—"}</b></div>` +
        `<div class="mf-dash-strip">${strip}</div>` +
        `<div class="mf-dash-meta"><span>TAM $${e.tam_2035_busd ?? "—"} млрд</span>` +
        `<span>квант ${e.quantum_potential ?? "—"}/10</span><span>${t.timeline || ""}</span></div>`;
      grid.appendChild(c);
    });

    view.innerHTML = "";
    view.appendChild(head);
    view.appendChild(agg);
    view.appendChild(grid);
    if (meta.honesty_note)
      view.appendChild(el("p", "mf-disclaimer", `<strong>Честно:</strong> ${meta.honesty_note}`));
  }

  function renderPanel(view, t) {
    const meta = DATA._meta || {};
    view.className = "mf-panel";
    view.style.setProperty("--tab", t.color);
    view.innerHTML = "";

    // 1. title + subtitle + accent
    const head = el("div", "mf-head");
    head.appendChild(el("div", "mf-accent-bar"));
    head.appendChild(el("h2", "mf-title-h", t.title));
    head.appendChild(el("p", "mf-sub", t.subtitle));
    // 2. narrative_role badge — the per-tab positioning (e.g. OCM = "почему мы, а не IBM")
    head.appendChild(el("div", "mf-badge",
      `<span class="mf-badge-label">Нарративная роль</span>${t.narrative_role}`));
    view.appendChild(head);

    // 3 + 4. reasoning blocks
    const two = el("div", "mf-twocol");
    two.appendChild(block("Что оптимизируют квантовые вычисления", t.what_quantum_optimizes, "q"));
    two.appendChild(block("Почему классика не справляется", t.why_classical_fails, "c"));
    view.appendChild(two);

    // 5 + 6. figures
    const figs = el("div", "mf-figs");
    figs.appendChild(figure(t.molecules_png, "molecules", t.id));
    figs.appendChild(figure(t.descriptor_png, "descriptor", t.id));
    // optional 3rd figure by convention: NN_<name>_quantum_path.png — appended only
    // if the file actually loads (most tabs have none; it self-removes on 404).
    const qpFile = t.molecules_png && t.molecules_png.replace(/_molecules\.png$/i, "_quantum_path.png");
    if (qpFile && qpFile !== t.molecules_png) {
      const qf = el("figure", "mf-figure");
      const qi = el("img");
      qi.src = IMG_DIR + qpFile;
      qi.loading = "lazy";
      qi.alt = "Путь к квантовому преимуществу";
      qi.addEventListener("error", () => qf.remove());
      qf.appendChild(qi);
      qf.appendChild(el("figcaption", null,
        "Путь к квантовому преимуществу: предел DFT+U и рост активного пространства (иллюстративно)"));
      figs.appendChild(qf);
    }
    view.appendChild(figs);

    // 7. roadmap
    if (Array.isArray(meta.stage_model) && Array.isArray(t.stage_progress))
      view.appendChild(roadmap(meta.stage_model, t.stage_progress));

    // 8. economics
    view.appendChild(economics(t));

    // 9. honesty disclaimer
    view.appendChild(el("p", "mf-disclaimer",
      `<strong>Честно:</strong> ${meta.honesty_note || ""}`));
  }

  window.addEventListener("hashchange", () => {
    const id = location.hash.replace(/^#/, "");
    const next = TABS.some((t) => t.id === id) ? id : DASH;
    if (next !== active) { active = next; renderTabs(); applyView(); }
  });

  boot();
})();
