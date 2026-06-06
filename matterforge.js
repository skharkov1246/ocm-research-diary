/* ============================================================================
   MatterForge journal — tab renderer.
   Drives any page that provides #mf-tabs (the bar) and #mf-view (the panel).
   Content comes from assets/matterforge/matterforge_tabs_content.json — all five
   tabs share one schema, so one generic renderer fills the whole bar. Nothing on
   screen is invented; every value is read from the JSON.

   Context-adaptive:
     • standalone journal page (no diary)  → just the tabs, first one (H₂) active.
     • mounted inside the OCM diary index   → adds a "Дневник" tab (tab 0) and
       hides .hero/#toc/#diary/#glossary while a MatterForge tab is shown.

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
    active = TABS.some((t) => t.id === fromHash) ? fromHash
           : HAS_DIARY ? DIARY_ID : TABS[0].id;

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
    if (HAS_DIARY) nav.appendChild(tabButton(DIARY_ID, "var(--accent)", "📔 Дневник OCM"));
    TABS.forEach((t) => nav.appendChild(tabButton(t.id, t.color, t.title)));
  }

  function selectTab(id) {
    if (active === id) return;
    active = id;
    history.replaceState(null, "", id === DIARY_ID ? location.pathname + location.search : "#" + id);
    renderTabs();
    applyView();
    window.scrollTo({ top: 0, behavior: "smooth" });
  }

  function applyView() {
    const view = $("#mf-view");
    const showDiary = HAS_DIARY && active === DIARY_ID;
    if (HAS_DIARY) diarySections().forEach((s) => { s.hidden = !showDiary; });
    if (showDiary) {
      view.hidden = true;
      view.innerHTML = "";
      return;
    }
    view.hidden = false;
    renderPanel(view, TABS.find((t) => t.id === active));
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
    const next = TABS.some((t) => t.id === id) ? id : (HAS_DIARY ? DIARY_ID : active);
    if (next !== active) { active = next; renderTabs(); applyView(); }
  });

  boot();
})();
