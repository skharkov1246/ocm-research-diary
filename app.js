/* OCM diary — renderer. State = {lang: 'ru'|'en', level: 'simple'|'tech'}. */
const REPO_URL = "https://github.com/skharkov1246/ocm-research-diary";

const state = {
  lang: localStorage.getItem("ocm-lang") || "ru",
  level: localStorage.getItem("ocm-level") || "simple",
};

const $ = (sel) => document.querySelector(sel);
const t = (obj) => (obj ? (obj[state.lang] ?? obj.en ?? "") : "");
const el = (tag, cls, html) => {
  const n = document.createElement(tag);
  if (cls) n.className = cls;
  if (html != null) n.innerHTML = html;
  return n;
};

function renderHeader() {
  document.documentElement.lang = state.lang;
  document.body.dataset.lang = state.lang;
  document.body.dataset.level = state.level;

  $("#site-title").textContent = t(SITE.title);
  $("#site-subtitle").textContent = t(SITE.subtitle);
  $("#lang-label").textContent = t(UI.langLabel);
  $("#level-label").textContent = t(UI.levelLabel);
  $("#lvl-simple").textContent = t(UI.simple);
  $("#lvl-tech").textContent = t(UI.tech);

  $("#hero-text").innerHTML = SITE.hero[state.level][state.lang];
  $("#honesty-note").textContent = t(UI.honesty);
  $("#footer-honesty").textContent = t(UI.honesty);
  $("#repo-link").href = REPO_URL;

  document.querySelectorAll(".toggle").forEach((b) => {
    b.classList.toggle("active", state[b.dataset.set] === b.dataset.value);
  });
}

function renderToc() {
  const toc = $("#toc");
  toc.innerHTML = "";
  ENTRIES.forEach((e) => {
    const a = el("a", null, `<span class="dot">●</span>${t(e.stage)} — ${t(e.title)}`);
    a.href = "#" + e.id;
    toc.appendChild(a);
  });
}

function renderTable(tbl) {
  const wrap = el("div", "data-table");
  wrap.appendChild(el("div", "tt", t(tbl.title)));
  const table = el("table");
  const thead = el("thead");
  const htr = el("tr");
  t(tbl.head).forEach((h) => htr.appendChild(el("th", null, h)));
  thead.appendChild(htr);
  table.appendChild(thead);
  const tbody = el("tbody");
  tbl.rows.forEach((row) => {
    const tr = el("tr");
    row.forEach((c) => tr.appendChild(el("td", null, c)));
    tbody.appendChild(tr);
  });
  table.appendChild(tbody);
  wrap.appendChild(table);
  return wrap;
}

let figNum = 0;

function renderEntry(e, idx) {
  const card = el("section", "entry");
  card.id = e.id;

  // Optional per-entry accent: recolour everything that reads var(--accent*) inside this card.
  if (e.accent) {
    card.classList.add("accented");
    card.style.setProperty("--accent", e.accent);
    card.style.setProperty("--accent-2", e.accent);
  }

  const meta = el("div", "entry-meta");
  meta.appendChild(el("span", "stage-tag", t(e.stage)));
  meta.appendChild(el("span", "entry-date", e.date));
  card.appendChild(meta);

  card.appendChild(el("h2", null, t(e.title)));

  const body = el("div", "entry-body", e[state.level][state.lang]);
  card.appendChild(body);

  // Optional honest economics caveat — shown, never hidden.
  if (e.economics && e.economics.caveat) {
    const ec = el("div", "economics-caveat");
    ec.innerHTML =
      `<span class="ec-label">${t(UI.economicsLabel)}</span>` + t(e.economics.caveat);
    card.appendChild(ec);
  }

  if (e.table) card.appendChild(renderTable(e.table));

  if (e.figures && e.figures.length) {
    const figs = el("div", "figures");
    e.figures.forEach((f, i) => {
      const fig = el("figure", "figure");
      const img = el("img");
      img.src = f.src;
      img.loading = "lazy";
      img.alt = t(f.caption);
      figNum++;
      const cap = el("figcaption", null,
        `<span class="fignum">${t(UI.figurePrefix)} ${figNum}</span>${t(f.caption)}`);
      fig.appendChild(img);
      fig.appendChild(cap);
      figs.appendChild(fig);
    });
    card.appendChild(figs);
  }
  return card;
}

function renderDiary() {
  const d = $("#diary");
  d.innerHTML = "";
  figNum = 0;
  ENTRIES.forEach((e, i) => d.appendChild(renderEntry(e, i)));
}

function renderGlossary() {
  $("#glossary-title").textContent = t(UI.glossaryTitle);
  const list = $("#glossary-list");
  list.innerHTML = "";
  GLOSSARY.forEach((g) => {
    list.appendChild(el("dt", null, t(g.term)));
    list.appendChild(el("dd", null, t(g.def)));
  });
}

function renderAll() {
  renderHeader();
  renderToc();
  renderDiary();
  renderGlossary();
}

document.addEventListener("click", (ev) => {
  const btn = ev.target.closest(".toggle");
  if (!btn) return;
  state[btn.dataset.set] = btn.dataset.value;
  localStorage.setItem("ocm-" + btn.dataset.set, btn.dataset.value);
  renderAll();
});

renderAll();
