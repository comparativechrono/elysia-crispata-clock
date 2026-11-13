/* Elysia clocks: comparison page. Data arrive as window.GENUS from data/genus.js,
   written by scripts/genus/20_build_genus_site.py. */
(function () {
  "use strict";
  const D = window.GENUS;
  if (!D) return;
  const SVGNS = "http://www.w3.org/2000/svg";
  const byId = (id) => document.getElementById(id);
  function el(tag, attrs, kids) {
    const n = document.createElement(tag);
    if (attrs) for (const k in attrs) {
      if (k === "class") n.className = attrs[k];
      else if (k === "text") n.textContent = attrs[k];
      else if (k === "html") n.innerHTML = attrs[k];
      else if (k.startsWith("on")) n.addEventListener(k.slice(2), attrs[k]);
      else if (attrs[k] !== null && attrs[k] !== undefined) n.setAttribute(k, attrs[k]);
    }
    (kids || []).forEach((c) => c && n.append(c));
    return n;
  }
  function sv(tag, attrs, kids) {
    const n = document.createElementNS(SVGNS, tag);
    if (attrs) for (const k in attrs) {
      if (k === "text") n.textContent = attrs[k];
      else if (attrs[k] !== null && attrs[k] !== undefined) n.setAttribute(k, attrs[k]);
    }
    (kids || []).forEach((c) => c && n.append(c));
    return n;
  }
  const STATUS = {
    complete: { cls: "ok", word: "intact" }, complete_discrepancy: { cls: "warn", word: "whole, sequence does not translate cleanly" },
    fragmented: { cls: "part", word: "in pieces" }, partial: { cls: "part", word: "partial" },
    divergent: { cls: "miss", word: "divergent copy only" }, absent: { cls: "miss", word: "not found" },
  };
  const SP = D.species, G = {};
  D.genes.forEach((g) => (G[g.symbol] = g));
  const NAME = {};
  SP.forEach((s) => (NAME[s.tag] = s.name));
  Object.keys(D.outgroups).forEach((t) => (NAME[t] = D.outgroups[t]));
  const pair = (a, b) => [a, b].sort().join("|");
  const pct = (x) => (x === null || x === undefined ? "–" : Math.round(x * 100) + "%");
  const glyph = (cls) => el("i", { class: "glyph " + cls, "aria-hidden": "true" });

  /* ---------- gene by species */
  function renderMatrix() {
    const t = byId("matrix");
    if (!t) return;
    const head = el("tr", {}, [el("th", { class: "g", text: "Gene" })].concat(SP.map((s) => el("th", {}, [el("a", { href: s.slug + "/index.html", html: "<i>" + s.short + "</i>" })]))));
    const body = el("tbody");
    Object.keys(D.categories).forEach((cat) => {
      const gs = D.genes.filter((g) => g.category === cat);
      if (!gs.length) return;
      body.append(el("tr", { class: "grp" }, [el("td", { colspan: SP.length + 1, text: D.categories[cat] })]));
      gs.forEach((g) => {
        const cells = SP.map((s) => {
          const st = STATUS[g.status[s.slug]] || STATUS.absent;
          const id = s.tag === "ELYCR" ? null : g.identity[pair(s.tag, "ELYCR")];
          const comp = g.completeness[s.slug], copies = g.copies[s.slug];
          const tip = `${s.name}: ${st.word}` + (comp ? `, fills ${Math.round(comp * 100)}% of the shared positions` : "") + (copies && copies > 1 ? `, ${copies} places in the assembly` : "");
          return el("td", {}, [el("a", { class: "cell", href: `${s.slug}/index.html#gene/${encodeURIComponent(g.symbol)}`, title: tip, "aria-label": tip },
            [glyph(st.cls), document.createTextNode(s.tag === "ELYCR" ? "" : (g.status[s.slug] === "absent" ? "" : pct(id)))])]);
        });
        body.append(el("tr", {}, [el("td", { class: "g" }, [el("button", { type: "button", text: g.label, title: g.name, onclick: () => pick(g.symbol, true) })])].concat(cells)));
      });
    });
    t.replaceChildren(el("thead", {}, [head]), body);
  }

  /* ---------- conservation bars */
  function renderCons() {
    const host = byId("cons");
    if (!host) return;
    const gs = D.genes.filter((g) => g.mean_identity && g.n_compared >= 2).sort((a, b) => b.mean_identity - a.mean_identity);
    const W = 1040, X0 = 150, PW = W - X0 - 70, ROW = 19, TOP = 22, lo = 0.4;
    const H = TOP + gs.length * ROW + 8;
    const x = (v) => X0 + ((Math.max(lo, v) - lo) / (1 - lo)) * PW;
    const svg = sv("svg", { viewBox: `0 0 ${W} ${H}`, class: "cons", role: "img", "aria-label": "Mean pairwise identity of each clock protein among Elysia species" });
    [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1].forEach((t) => {
      svg.append(sv("line", { x1: x(t), x2: x(t), y1: TOP - 6, y2: H - 6, stroke: "var(--grid)" }));
      svg.append(sv("text", { x: x(t), y: 12, "text-anchor": "middle", text: Math.round(t * 100) + "%" }));
    });
    gs.forEach((g, i) => {
      const y = TOP + i * ROW;
      const row = sv("g", { class: "row", tabindex: 0, role: "button", "aria-label": `${g.label}: ${pct(g.mean_identity)} mean identity, lowest pair ${pct(g.min_identity)}` }, [
        sv("title", { text: `${g.label} (${g.name}): mean ${pct(g.mean_identity)}, lowest pair ${pct(g.min_identity)}, ${g.n_compared} species compared` }),
        sv("rect", { x: 0, y: y - 2, width: W, height: ROW, fill: "transparent" }),
        sv("text", { x: X0 - 10, y: y + 11, "text-anchor": "end", text: g.label }),
        sv("rect", { class: "bar", x: X0, y: y + 1, width: Math.max(2, x(g.mean_identity) - X0), height: 12, rx: 3 }),
        sv("text", { class: "val", x: x(g.mean_identity) + 7, y: y + 11, text: pct(g.mean_identity) }),
      ]);
      row.addEventListener("click", () => pick(g.symbol, true));
      row.addEventListener("keydown", (e) => { if (e.key === "Enter") pick(g.symbol, true); });
      svg.append(row);
    });
    host.replaceChildren(svg);
  }

  /* ---------- one protein, all species */
  function renderPick() {
    const sel = byId("pick");
    if (!sel) return;
    Object.keys(D.categories).forEach((cat) => {
      const gs = D.genes.filter((g) => g.category === cat);
      if (!gs.length) return;
      sel.append(el("optgroup", { label: D.categories[cat] }, gs.map((g) => el("option", { value: g.symbol, text: `${g.label}: ${g.name}` }))));
    });
    sel.addEventListener("change", () => pick(sel.value, false));
  }
  function pick(sym, scroll) {
    const g = G[sym];
    if (!g) return;
    const sel = byId("pick");
    if (sel) sel.value = sym;
    const tags = SP.map((s) => s.tag).filter((t) => g.alignment[t]);
    const outs = Object.keys(D.outgroups).filter((t) => g.alignment[t]);
    const note = byId("pick-note");
    if (note) note.textContent = g.note || "";
    // identity table
    const idm = byId("idmatrix");
    if (idm) {
      const all = tags.concat(outs);
      const head = el("tr", {}, [el("th", { text: "Identity" })].concat(all.map((t) => el("th", { class: "num", html: "<i>" + NAME[t].replace(/^(\w)\w+ /, "$1. ") + "</i>" }))));
      const rows = all.map((a, i) => el("tr", {}, [el("td", { html: "<i>" + NAME[a] + "</i>" + (outs.includes(a) ? ' <span class="muted small">other gastropod</span>' : "") })].concat(
        all.map((b, j) => el("td", { class: "num", text: i === j ? "" : (j > i ? "" : pct(g.identity[pair(a, b)])) })))));
      idm.replaceChildren(el("div", { class: "tablewrap" }, [el("table", { class: "data" }, [el("thead", {}, [head]), el("tbody", {}, rows)])]));
    }
    // alignment
    const host = byId("aln");
    if (host) {
      const L = g.alignment[tags[0]].length;
      const cons = [];
      for (let i = 0; i < L; i++) {
        const c = {};
        tags.forEach((t) => { const r = g.alignment[t][i]; if (r !== "-") c[r] = (c[r] || 0) + 1; });
        let best = "-", n = 0;
        for (const r in c) if (c[r] > n) { best = r; n = c[r]; }
        cons.push(best);
      }
      const esc = (s) => s.replace(/&/g, "&amp;").replace(/</g, "&lt;");
      let ruler = "";
      for (let i = 1; i <= L; i += 10) ruler += String(i).padEnd(10, " ");
      const line = (t, og) => {
        const s = g.alignment[t];
        let out = "", run = "", cls = null;
        const flush = () => { if (run) out += cls ? `<span class="${cls}">${run}</span>` : run; run = ""; };
        for (let i = 0; i < L; i++) {
          const r = s[i], k = r === "-" ? "gp" : (r !== cons[i] ? "v" : "");
          if (k !== cls) { flush(); cls = k; }
          run += r;
        }
        flush();
        return `<span class="nm${og ? " og" : ""}" title="${esc(NAME[t])}">${esc(NAME[t].replace(/^(\w)\w+ /, "$1. "))}</span>${og ? '<span class="og">' + out + "</span>" : out}`;
      };
      host.innerHTML = `<span class="nm"> </span><span class="ruler">${ruler.slice(0, L)}</span>\n` + tags.map((t) => line(t, false)).join("\n") +
        (outs.length ? "\n" + `<span class="nm"> </span><span class="gp">${"·".repeat(Math.min(L, 60))}</span>\n` + outs.map((t) => line(t, true)).join("\n") : "");
    }
    if (scroll) {
      const target = byId("pick");
      if (target) target.scrollIntoView({ behavior: window.matchMedia("(prefers-reduced-motion: reduce)").matches ? "auto" : "smooth", block: "start" });
    }
    history.replaceState(null, "", "#gene/" + encodeURIComponent(sym));
  }

  /* ---------- species tree */
  function parseNewick(s) {
    let i = 0;
    function node() {
      const n = { kids: [], name: "", len: 0, sup: null };
      if (s[i] === "(") {
        i++;
        for (;;) { n.kids.push(node()); if (s[i] === ",") { i++; continue; } break; }
        i++; // )
        let lab = "";
        while (i < s.length && !":,);".includes(s[i])) lab += s[i++];
        if (lab) n.sup = parseFloat(lab);
      } else {
        while (i < s.length && !":,);".includes(s[i])) n.name += s[i++];
      }
      if (s[i] === ":") { i++; let v = ""; while (i < s.length && !",);".includes(s[i])) v += s[i++]; n.len = parseFloat(v); }
      return n;
    }
    return node();
  }
  function renderTree() {
    const host = byId("tree");
    if (!host || !D.tree) return;
    const root = parseNewick(D.tree.trim());
    const leaves = [];
    (function walk(n, d) { n.depth = d + n.len; if (!n.kids.length) { n.y = leaves.length; leaves.push(n); } else { n.kids.forEach((k) => walk(k, n.depth)); n.y = (n.kids[0].y + n.kids[n.kids.length - 1].y) / 2; } })(root, 0);
    root.depth = 0;
    const maxD = Math.max(...leaves.map((l) => l.depth));
    const W = 880, X0 = 16, PW = 520, ROW = 30, TOP = 14, H = TOP + leaves.length * ROW + 34;
    const x = (d) => X0 + (d / maxD) * PW, y = (n) => TOP + n.y * ROW + 12;
    const svg = sv("svg", { viewBox: `0 0 ${W} ${H}`, class: "tree", role: "img", "aria-label": "Tree of the species from concatenated clock proteins" });
    (function draw(n, px) {
      svg.append(sv("path", { class: "br", d: `M${px},${y(n)}H${x(n.depth)}` }));
      if (n.kids.length) {
        svg.append(sv("path", { class: "br", d: `M${x(n.depth)},${y(n.kids[0])}V${y(n.kids[n.kids.length - 1])}` }));
        if (n.sup !== null && !isNaN(n.sup)) svg.append(sv("text", { class: "sup", x: x(n.depth) - 4, y: y(n) - 4, "text-anchor": "end", text: Math.round(n.sup) }));
        n.kids.forEach((k) => draw(k, x(n.depth)));
      } else {
        const og = !!D.outgroups[n.name];
        const t = sv("text", { x: x(n.depth) + 8, y: y(n) + 4, "font-style": "italic", class: og ? "og" : "", text: NAME[n.name] || n.name });
        svg.append(t);
        if (og) svg.append(sv("text", { class: "sup", x: x(n.depth) + 8 + (NAME[n.name] || n.name).length * 7.2 + 8, y: y(n) + 4, text: "other gastropod" }));
      }
    })(root, X0);
    const bar = 0.1;
    svg.append(sv("path", { class: "br", d: `M${X0},${H - 14}H${X0 + (bar / maxD) * PW}` }));
    svg.append(sv("text", { class: "sup", x: X0 + (bar / maxD) * PW + 8, y: H - 10, text: "0.1 substitutions per site" }));
    host.replaceChildren(svg);
  }

  /* ---------- linkage */
  function renderLinkage() {
    const t = byId("linkage");
    if (!t) return;
    const cs = D.chromosome_species.map((slug) => SP.find((s) => s.slug === slug));
    const rows = D.genes.filter((g) => cs.every((s) => g.chrom[s.slug])).sort((a, b) => (+a.chrom[cs[0].slug]) - (+b.chrom[cs[0].slug]) || a.label.localeCompare(b.label));
    const head = el("tr", {}, [el("th", { text: "Gene" })].concat(cs.map((s) => el("th", { class: "num", html: "Chromosome in <i>" + s.short + "</i>" }))));
    let last = null;
    const body = el("tbody");
    rows.forEach((g) => {
      const c0 = g.chrom[cs[0].slug];
      const tr = el("tr", {}, [el("td", {}, [el("strong", { text: g.label })])].concat(cs.map((s) => el("td", { class: "num", text: g.chrom[s.slug] }))));
      if (c0 !== last && last !== null) tr.style.borderTop = "2px solid var(--rule-strong)";
      last = c0;
      body.append(tr);
    });
    t.replaceChildren(el("thead", {}, [head]), body);
  }

  function init() {
    renderMatrix(); renderCons(); renderPick(); renderTree(); renderLinkage();
    const m = location.hash.match(/^#gene\/(.+)$/);
    pick(m && G[decodeURIComponent(m[1])] ? decodeURIComponent(m[1]) : "Period", false);
    if (!m) history.replaceState(null, "", location.pathname + location.search);
  }
  (document.fonts && document.fonts.ready ? document.fonts.ready : Promise.resolve()).then(init);
})();
