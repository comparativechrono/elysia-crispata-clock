/* Elysia crispata clock — gene browser (index page).
   Data arrive as window.CLOCK from data/clock.js, which scripts/30_build_site.py writes. */
(function () {
  "use strict";
  const D = window.CLOCK;
  const SVGNS = "http://www.w3.org/2000/svg";
  const byId = (id) => document.getElementById(id);
  const G = {};
  D.genes.forEach((g) => (G[g.symbol] = g));

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
      else if (k.startsWith("on")) n.addEventListener(k.slice(2), attrs[k]);
      else if (attrs[k] !== null && attrs[k] !== undefined) n.setAttribute(k, attrs[k]);
    }
    (kids || []).forEach((c) => c && n.append(c));
    return n;
  }
  const fmt = (n) => Number(n).toLocaleString("en-GB");
  const mb = (bp) => (bp / 1e6).toFixed(2) + " Mb";

  /* ---------- status vocabulary: one word and one glyph per state, used everywhere */
  const STATUS = {
    complete: { cls: "ok", word: "intact" },
    complete_discrepancy: { cls: "warn", word: "whole, sequence differs" },
    fragmented: { cls: "part", word: "in pieces" },
    partial: { cls: "part", word: "partial" },
    divergent: { cls: "miss", word: "divergent copy only" },
    absent: { cls: "miss", word: "missing" },
  };
  const stOf = (g, asm) => STATUS[g.assemblies[asm || "sanger"].status];
  function glyph(cls) { return el("i", { class: "glyph " + cls, "aria-hidden": "true" }); }
  function statusHtml(g, asm, withWord) {
    const s = stOf(g, asm);
    const a = g.assemblies[asm || "sanger"];
    let word = s.word;
    if (a.status === "partial") word = "partial, " + Math.round(a.coverage * 100) + "%";
    return el("span", { class: "st" }, [glyph(s.cls), withWord === false ? el("span", { class: "sr", text: word }) : document.createTextNode(word)]);
  }
  function svgGlyph(cls, cx, cy, r) {
    r = r || 5;
    const css = getComputedStyle(document.documentElement);
    const col = { ok: "--st-ok", warn: "--st-warn", part: "--st-serious", miss: "--st-crit" }[cls];
    const c = "var(" + col + ")";
    const g = sv("g", { "aria-hidden": "true" });
    if (cls === "ok") g.append(sv("circle", { cx, cy, r, fill: c }));
    else if (cls === "warn") g.append(sv("circle", { cx, cy, r, fill: c, stroke: "rgba(0,0,0,.3)", "stroke-width": 1 }));
    else if (cls === "part") {
      g.append(sv("path", { d: `M${cx} ${cy - r}A${r} ${r} 0 0 0 ${cx} ${cy + r}Z`, fill: c }));
      g.append(sv("circle", { cx, cy, r: r - 0.8, fill: "none", stroke: c, "stroke-width": 1.6 }));
    } else g.append(sv("circle", { cx, cy, r: r - 1, fill: "none", stroke: c, "stroke-width": 2 }));
    return g;
  }

  /* ---------- tooltip */
  const tip = el("div", { class: "tip-pop", role: "tooltip" });
  document.body.append(tip);
  function showTip(html, x, y) {
    tip.innerHTML = html; tip.classList.add("on");
    const r = tip.getBoundingClientRect();
    let left = x + 14, top = y + 14;
    if (left + r.width > innerWidth - 8) left = x - r.width - 14;
    if (top + r.height > innerHeight - 8) top = y - r.height - 14;
    tip.style.left = Math.max(6, left) + "px"; tip.style.top = Math.max(6, top) + "px";
  }
  function hideTip() { tip.classList.remove("on"); }
  function tipFor(node, htmlFn) {
    node.addEventListener("mousemove", (e) => showTip(htmlFn(), e.clientX, e.clientY));
    node.addEventListener("mouseleave", hideTip);
    node.addEventListener("focus", () => { const r = node.getBoundingClientRect(); showTip(htmlFn(), r.right, r.top); });
    node.addEventListener("blur", hideTip);
  }
  function geneTip(g) {
    const a = g.assemblies.sanger;
    const loc = a.pieces.length ? a.pieces.map((p) => `${chrName(p.chrom)}:${fmt(p.start)}`).join(", ") : "not in the assembly";
    return `<b>${g.label}</b><br>${g.name}<br>${fmt(g.length)} aa, ${stOf(g).word} in xgElyCris1.1<br>${loc}`;
  }
  const chrName = (c) => (/^\d+$/.test(c) ? "chr" + c : c);

  /* ---------- schematic */
  const LAYOUT = {
    zones: [
      { id: "light", x: 14, y: 14, w: 252, h: 562, label: "Light input" },
      { id: "loop2", x: 282, y: 14, w: 616, h: 136, label: "Second loop", note: "sets when Clock and Bmal1 are transcribed" },
      { id: "core", x: 282, y: 164, w: 616, h: 190, label: "Core loop" },
      { id: "pace", x: 282, y: 368, w: 616, h: 208, label: "Pace-setters", note: "phosphorylation and degradation of Period and Timeless" },
      { id: "out", x: 914, y: 14, w: 252, h: 562, label: "Output", note: "what the clock drives" },
    ],
    flows: [
      { zone: "light", x: 30, y: 62, w: 224, genes: ["Cry1"], sub: null },
      { zone: "light", x: 30, y: 150, w: 224, genes: ["Ropsin1", "Ropsin2", "RopsinL1", "RopsinL2", "Xenopsin", "GoOpsin1", "GoOpsin2", "Neuropsin", "Retinochrome"], sub: "Opsins" },
      { zone: "light", x: 30, y: 408, w: 224, genes: ["Phr6-4", "PhrCPD", "CryDASH1", "CryDASH2"], sub: "Photolyases (Cry relatives, DNA repair)" },
      { zone: "pace", x: 392, y: 400, w: 500, genes: ["Ck1de", "Ck1a", "Ck2a", "Ck2b", "Gsk3", "Nemo"], side: "Kinases", stack: "pace" },
      { zone: "pace", x: 392, y: 0, w: 500, genes: ["Pp1a", "Pp1b", "Pp2aC", "Pp2aB55", "Pp2aB56e", "Pp2aB56d"], side: "Phosphatases", stack: "pace" },
      { zone: "pace", x: 392, y: 0, w: 500, genes: ["Slimb", "Jetlag", "Sirt1"], side: "Degradation", stack: "pace" },
      { zone: "out", x: 930, y: 62, w: 224, genes: ["Pdf", "Pdfr"], sub: null },
      { zone: "out", x: 930, y: 150, w: 224, genes: ["Lark", "Atx2", "Nocturnin", "NonA"], sub: "RNA-binding and clock-controlled" },
      { zone: "out", x: 930, y: 276, w: 224, genes: ["Tipin"], sub: "Partner of Timeout" },
    ],
    fixed: [
      { gene: "Pdp1", x: 330, y: 82 }, { gene: "Vrille", x: 410, y: 82 },
      { gene: "E75", x: 520, y: 82 }, { gene: "Hr3", x: 592, y: 82 },
      { gene: "Cwo", x: 760, y: 82 },
      { gene: "Clock", x: 372, y: 250 }, { gene: "Bmal1", x: 452, y: 250 },
      { gene: "Period", x: 664, y: 226 }, { gene: "Timeless", x: 752, y: 226 },
      { gene: "Cry2", x: 664, y: 270 }, { gene: "Timeout", x: 752, y: 270 },
    ],
  };
  const PILL_H = 28;
  let measureCtx = null;
  function textW(s) {
    if (!measureCtx) { measureCtx = document.createElement("canvas").getContext("2d"); }
    measureCtx.font = '550 13px "Bricolage Grotesque", "Segoe UI", system-ui, sans-serif';
    return measureCtx.measureText(s).width;
  }
  const pillW = (g) => Math.ceil(textW(g.label)) + 36;
  function pill(g, x, y) {
    const w = pillW(g);
    const n = sv("g", { class: "pill", tabindex: 0, role: "button", "aria-label": `${g.label}: ${g.name}. ${stOf(g).word} in xgElyCris1.1`, transform: `translate(${x} ${y})` }, [
      sv("rect", { x: 0, y: 0, width: w, height: PILL_H, rx: PILL_H / 2 }),
      svgGlyph(stOf(g).cls, 14, PILL_H / 2, 5),
      sv("text", { x: 25, y: 18.5, text: g.label }),
    ]);
    n.addEventListener("click", () => openGene(g.symbol));
    n.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); openGene(g.symbol); } });
    tipFor(n, () => geneTip(g));
    return { node: n, w };
  }
  function renderSchematic() {
    const host = byId("schematic");
    if (!host) return;
    const svg = sv("svg", { class: "schematic", viewBox: "0 0 1180 590", role: "group", "aria-label": "Schematic of the clock genes found in Elysia crispata" });
    const defs = sv("defs");
    defs.innerHTML =
      '<marker id="arr" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M1 1L9 5L1 9" fill="none" stroke="context-stroke" stroke-width="1.6"/></marker>' +
      '<marker id="bar" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="9" markerHeight="9" orient="auto"><path d="M5 0.5V9.5" stroke="context-stroke" stroke-width="2"/></marker>';
    svg.append(defs);
    const pos = {};
    LAYOUT.zones.forEach((z) => {
      svg.append(sv("rect", { class: "zone", x: z.x, y: z.y, width: z.w, height: z.h, rx: 14 }));
      svg.append(sv("text", { class: "zone-label", x: z.x + 16, y: z.y + 26, text: z.label }));
      if (z.note) svg.append(sv("text", { class: "zone-note", x: z.x + 16 + textW(z.label) + 16, y: z.y + 26, text: z.note }));
    });
    const wires = sv("g");
    svg.append(wires);
    const stackY = {};
    LAYOUT.flows.forEach((f) => {
      let x = f.x, y = f.stack && stackY[f.stack] ? stackY[f.stack] : f.y;
      if (f.sub) svg.append(sv("text", { class: "zone-note", x: f.x, y: y - 10, text: f.sub }));
      if (f.side) svg.append(sv("text", { class: "zone-note", x: 298, y: y + 18, text: f.side }));
      f.genes.forEach((s) => {
        const g = G[s]; if (!g) return;
        const w = pillW(g);
        if (x + w > f.x + f.w) { x = f.x; y += PILL_H + 6; }
        const p = pill(g, x, y); svg.append(p.node);
        pos[s] = { x, y, w }; x += w + 8;
      });
      if (f.stack) stackY[f.stack] = y + PILL_H + 12;
    });
    LAYOUT.fixed.forEach((f) => { const g = G[f.gene]; if (!g) return; const p = pill(g, f.x, f.y); svg.append(p.node); pos[f.gene] = { x: f.x, y: f.y, w: p.w }; });
    const cx = (s) => pos[s].x + pos[s].w / 2, top = (s) => pos[s].y, bot = (s) => pos[s].y + PILL_H, right = (s) => pos[s].x + pos[s].w, left = (s) => pos[s].x;
    function wire(d, cls, end) { wires.append(sv("path", { class: "wire" + (cls ? " " + cls : ""), d, "marker-end": `url(#${end || "arr"})` })); }
    function label(x, y, t, anchor) { svg.append(sv("text", { class: "wire-label", x, y, "text-anchor": anchor || "start", text: t })); }
    if (pos.Clock && pos.Bmal1) {
      // heterodimer bracket
      const bx = left("Clock") - 10, bw = right("Bmal1") - left("Clock") + 20;
      wires.append(sv("rect", { x: bx, y: top("Clock") - 10, width: bw, height: PILL_H + 20, rx: 22, fill: "none", stroke: "var(--rule-strong)", "stroke-dasharray": "2 4" }));
      label(bx + bw / 2, top("Clock") - 18, "heterodimer", "middle");
      // activation of E-box targets
      const ax = bx + bw + 6, ay = top("Clock") + PILL_H / 2, tx = left("Period") - 14;
      wire(`M${ax} ${ay}H${tx}`);
      label((ax + tx) / 2, ay - 8, "switch on", "middle");
      // repression feedback
      const fy = bot("Cry2") + 30;
      wire(`M${cx("Cry2") + 40} ${bot("Cry2") + 6}V${fy}H${bx + bw / 2}V${bot("Clock") + 14}`, "", "bar");
      label((cx("Cry2") + bx + bw / 2) / 2 + 20, fy + 16, "proteins return and switch the pair off", "middle");
      // second loop inputs
      wire(`M${(cx("Pdp1") + cx("Vrille")) / 2} ${bot("Pdp1") + 4}V${top("Clock") - 34}H${cx("Clock")}V${top("Clock") - 12}`, "soft");
      wire(`M${(cx("E75") + cx("Hr3")) / 2} ${bot("E75") + 4}V${top("Bmal1") - 34}H${cx("Bmal1")}V${top("Bmal1") - 12}`, "soft");
      label((cx("Pdp1") + cx("Vrille")) / 2 + 8, bot("Pdp1") + 22, "on / off: Clock");
      label((cx("E75") + cx("Hr3")) / 2 + 8, bot("E75") + 22, "off / on: Bmal1");
      wire(`M${cx("Cwo")} ${bot("Cwo") + 4}V${top("Period") - 8}`, "soft", "bar");
      label(cx("Cwo") - 8, bot("Cwo") + 22, "damps E-box genes", "end");
      // pace-setters act on Period / Timeless
      const px = (right("Period") + left("Timeless")) / 2;
      wire(`M${px} ${400 - 8}V${bot("Cry2") + 8}`, "soft");
      label(px + 8, 384, "tune repressor lifetime");
      // light: Cry1 feeds the core loop
      if (pos.Cry1) {
        const ly = top("Clock") + PILL_H / 2 + 44;
        wire(`M${right("Cry1") + 4} ${top("Cry1") + PILL_H / 2}H${274}V${ly}H${296}`, "soft");
        label(30, top("Cry1") + PILL_H + 18, "light resets the core loop");
      }
      // output
      wire(`M${right("Timeless") + 14} ${top("Timeless") + PILL_H / 2}H${906}`, "soft");
      if (pos.Pdf && pos.Pdfr) wire(`M${right("Pdf") + 3} ${top("Pdf") + PILL_H / 2}H${left("Pdfr") - 4}`, "");
    }
    host.replaceChildren(svg);
    // phone-width fallback: plain grouped lists
    const list = byId("schematic-list");
    if (list) {
      const groups = [
        ["Core loop", ["Clock", "Bmal1", "Period", "Timeless", "Cry2", "Timeout"]],
        ["Second loop", ["Pdp1", "Vrille", "E75", "Hr3", "Cwo"]],
        ["Pace-setters", ["Ck1de", "Ck1a", "Ck2a", "Ck2b", "Gsk3", "Nemo", "Pp1a", "Pp1b", "Pp2aC", "Pp2aB55", "Pp2aB56e", "Pp2aB56d", "Slimb", "Jetlag", "Sirt1"]],
        ["Light input", ["Cry1", "Ropsin1", "Ropsin2", "RopsinL1", "RopsinL2", "Xenopsin", "GoOpsin1", "GoOpsin2", "Neuropsin", "Retinochrome", "Phr6-4", "PhrCPD", "CryDASH1", "CryDASH2"]],
        ["Output", ["Pdf", "Pdfr", "Lark", "Atx2", "Nocturnin", "NonA", "Tipin"]],
      ];
      list.replaceChildren(...groups.flatMap(([h, ss]) => [el("h3", { text: h }), el("div", { class: "chips" }, ss.filter((s) => G[s]).map(chip))]));
    }
  }
  function chip(s) {
    const g = G[s];
    return el("button", { class: "chip", type: "button", onclick: () => openGene(s), "aria-label": `${g.label}: ${stOf(g).word} in xgElyCris1.1` }, [glyph(stOf(g).cls), document.createTextNode(g.label)]);
  }

  /* ---------- karyotype */
  function renderKaryo() {
    const host = byId("karyo");
    if (!host) return;
    const chrs = D.chromosomes, maxLen = Math.max(...chrs.map((c) => c.length));
    const X0 = 46, W = 1100, PW = W - X0 - 30, ROW = 56, TOP = 26;
    const H = TOP + chrs.length * ROW + 26;
    const svg = sv("svg", { viewBox: `0 0 ${W} ${H}`, role: "group", "aria-label": "Positions of clock genes on the 15 chromosomes of xgElyCris1.1" });
    const x = (bp) => X0 + (bp / maxLen) * PW;
    const axis = sv("g", { class: "axis" });
    for (let t = 0; t <= maxLen; t += 10e6) {
      axis.append(sv("line", { x1: x(t), x2: x(t), y1: TOP - 12, y2: H - 22 }));
      axis.append(sv("text", { x: x(t), y: H - 8, "text-anchor": "middle", text: t / 1e6 + (t === 0 ? " Mb" : "") }));
    }
    svg.append(axis);
    chrs.forEach((c, i) => {
      const y = TOP + i * ROW + 20;
      svg.append(sv("text", { class: "chr-label", x: X0 - 10, y: y + 8, "text-anchor": "end", text: c.name }));
      svg.append(sv("rect", { class: "chr", x: X0, y, width: Math.max(4, x(c.length) - X0), height: 8, rx: 4 }));
      const marks = [];
      D.genes.forEach((g) => g.assemblies.sanger.pieces.forEach((p, k) => { if (p.chrom === c.name && k === 0) marks.push({ g, bp: (p.start + p.end) / 2 }); }));
      marks.sort((a, b) => a.bp - b.bp);
      const used = { up: -1e9, down: -1e9 };
      marks.forEach((m) => {
        const mx = x(m.bp), w = textW(m.g.label) * 0.86 + 18;
        let side = "up";
        if (mx - w / 2 < used.up + 4) side = mx - w / 2 < used.down + 4 ? "up2" : "down";
        let lx = mx;
        if (side === "up") used.up = mx + w / 2;
        else if (side === "down") used.down = mx + w / 2;
        else { lx = Math.max(mx, used.up + w / 2 + 4); used.up = lx + w / 2; }
        const up = side !== "down";
        const ty = up ? y - 8 : y + 25;
        const gnode = sv("g", { class: "mark", tabindex: 0, role: "button", "aria-label": `${m.g.label} on chromosome ${c.name} at ${mb(m.bp)}` }, [
          sv("line", { x1: mx, x2: mx, y1: y - 3, y2: y + 11 }),
          lx !== mx ? sv("line", { x1: mx, x2: lx, y1: y - 3, y2: y - 5, "stroke-width": 1 }) : null,
          svgGlyph(stOf(m.g).cls, lx - w / 2 + 6, ty - 4, 4),
          sv("text", { x: lx - w / 2 + 14, y: ty, text: m.g.label }),
        ]);
        gnode.addEventListener("click", () => openGene(m.g.symbol));
        gnode.addEventListener("keydown", (e) => { if (e.key === "Enter") openGene(m.g.symbol); });
        tipFor(gnode, () => geneTip(m.g));
        svg.append(gnode);
      });
    });
    host.replaceChildren(svg);
    const off = byId("offmap");
    if (off) {
      const unplaced = D.genes.filter((g) => g.assemblies.sanger.pieces.length && !g.assemblies.sanger.pieces.some((p) => /^\d+$/.test(p.chrom)));
      const none = D.genes.filter((g) => !g.assemblies.sanger.pieces.length);
      off.replaceChildren(
        el("div", {}, [el("h4", { text: `On scaffolds not yet placed on a chromosome (${unplaced.length})` }), el("div", { class: "chips" }, unplaced.map((g) => chip(g.symbol)))]),
        el("div", {}, [el("h4", { text: `Not in the assembly at all (${none.length})` }), el("div", { class: "chips" }, none.map((g) => chip(g.symbol)))])
      );
    }
  }

  /* ---------- table */
  const CATS = D.categories;
  const state = { cat: "all", status: "all", q: "" };
  function matches(g) {
    if (state.cat !== "all" && g.category !== state.cat) return false;
    const cls = stOf(g).cls;
    if (state.status === "intact" && cls !== "ok") return false;
    if (state.status === "problem" && cls === "ok") return false;
    if (state.q) {
      const hay = [g.symbol, g.label, g.name, g.ensembl_overlap.join(" "), g.ensembl_genes.join(" "), g.purdue_proteins.map((p) => p.accession + " " + p.locus_tag).join(" "), g.domains.map((d) => d.name).join(" "),
        g.orthologues.map((o) => o.gene).join(" ")].join(" ").toLowerCase();
      if (!state.q.toLowerCase().split(/\s+/).every((w) => hay.includes(w))) return false;
    }
    return true;
  }
  function renderFilters() {
    const host = byId("filters");
    if (!host) return;
    const mk = (group, value, text) => el("button", { class: "fbtn", type: "button", "aria-pressed": String(state[group] === value), text,
      onclick: () => { state[group] = value; renderFilters(); renderTable(); } });
    host.replaceChildren(
      el("input", { type: "search", placeholder: "Search gene, ID or domain", "aria-label": "Search genes", value: state.q, oninput: (e) => { state.q = e.target.value; renderTable(); } }),
      mk("status", "all", "All"), mk("status", "intact", "Intact in xgElyCris1.1"), mk("status", "problem", "Not intact"),
      el("span", { class: "muted small", text: "|" }),
      mk("cat", "all", "Every role"), ...Object.keys(CATS).map((c) => mk("cat", c, CATS[c])),
      el("span", { class: "count", id: "count" })
    );
  }
  function renderTable() {
    const host = byId("genes-body");
    if (!host) return;
    const rows = D.genes.filter(matches);
    const maxDepth = Math.max(...D.genes.map((g) => (g.reads_support ? g.reads_support.median_reads : 0)), 1);
    const maxTpm = Math.max(...D.genes.map((g) => (g.expression ? Math.log10(1 + g.expression.tpm) : 0)), 0.1);
    host.replaceChildren(...rows.map((g) => {
      const a = g.assemblies.sanger;
      const locus = a.pieces.length ? a.pieces.map((p) => `${chrName(p.chrom)}:${fmt(p.start)}–${fmt(p.end)}`).join("; ") : (a.best_other && a.status === "divergent" ? `closest: ${chrName(a.best_other.chrom)}:${fmt(a.best_other.start)}` : "");
      const tr = el("tr", { onclick: () => openGene(g.symbol) }, [
        el("td", { class: "sym" }, [el("button", { type: "button", text: g.label, "aria-label": `Open ${g.label}` }), el("span", { class: "gname", text: g.name })]),
        el("td", { text: CATS[g.category] || g.category }),
        el("td", { class: "num", text: fmt(g.length) }),
        el("td", {}, [statusHtml(g, "sanger"), locus ? el("span", { class: "locus", text: locus }) : null]),
        el("td", {}, [el("span", { class: "others" }, [["alt", "alt"], ["purdue", "Purdue"], ["tampa", "Tampa"]].map(([k, lab]) => {
          const s = el("span", {}, [glyph(stOf(g, k).cls), document.createTextNode(lab)]);
          s.title = `${D.assemblies[k]}: ${stOf(g, k).word}`;
          return s;
        }))]),
        el("td", { class: "num" }, g.reads_support ? [el("i", { class: "minibar", style: `width:${Math.round((g.reads_support.median_reads / maxDepth) * 54)}px` }), document.createTextNode(Math.round(g.reads_support.median_reads))] : [document.createTextNode("–")]),
        el("td", { class: "num" }, g.expression ? [el("i", { class: "minibar rna", style: `width:${Math.round((Math.log10(1 + g.expression.tpm) / maxTpm) * 54)}px` }), document.createTextNode(g.expression.tpm >= 10 ? g.expression.tpm.toFixed(0) : g.expression.tpm.toFixed(1))] : [document.createTextNode("–")]),
      ]);
      return tr;
    }));
    if (!rows.length) host.append(el("tr", {}, [el("td", { colspan: 7, class: "empty", text: "No gene matches. Clear the search or pick another filter." })]));
    const c = byId("count"); if (c) c.textContent = `${rows.length} of ${D.genes.length} genes`;
  }

  /* ---------- gene dialog */
  const dlg = byId("gene");
  function sentence(g) {
    const a = g.assemblies.sanger, s = a.status;
    const where = a.pieces.length ? a.pieces.map((p) => `${chrName(p.chrom)}:${fmt(p.start)}–${fmt(p.end)} (${p.strand} strand)`).join(" and ") : "";
    const elsewhere = ["purdue", "tampa", "alt"].filter((k) => g.assemblies[k].status === "complete").map((k) => ({ purdue: "the Purdue assembly", tampa: "the Tampa assembly", alt: "the alternate haplotype" }[k]));
    const reads = g.reads_support ? ` HiFi reads from the sequenced animal cover ${g.reads_support.covered} of ${g.reads_support.segments} exon segments (median ${Math.round(g.reads_support.median_reads)} reads each).` : "";
    if (s === "complete") return `Intact in xgElyCris1.1 at ${where}.`;
    if (s === "complete_discrepancy") {
      const p = a.pieces[0];
      return `The whole gene is in xgElyCris1.1 at ${where}, but the assembly sequence does not translate cleanly against the reference protein (${p.fs ? p.fs + " frameshift" + (p.fs > 1 ? "s" : "") : ""}${p.fs && p.stop ? ", " : ""}${p.stop ? p.stop + " in-frame stop" + (p.stop > 1 ? "s" : "") : ""}; ${(p.ident * 100).toFixed(1)}% identity). Treat the assembly sequence here with care; the protein below comes from ${g.model_kind === "ensembl" ? "the Ensembl model" : "another assembly"}.`;
    }
    if (s === "fragmented") return `All of the gene is in xgElyCris1.1 but in ${a.pieces.length} separate pieces (${where}) that are not in transcript order, so no single gene model can span it there. It is whole in ${elsewhere.join(" and ") || "no other assembly"}.`;
    if (s === "partial") return `Only ${Math.round(a.coverage * 100)}% of the protein is encoded in xgElyCris1.1 (${where}). The rest is missing from the assembly, not from the animal: it is whole in ${elsewhere.join(" and ") || "no other assembly"}.${reads}`;
    if (s === "divergent") return `xgElyCris1.1 has no clean copy. The closest locus (${chrName(a.best_other.chrom)}:${fmt(a.best_other.start)}–${fmt(a.best_other.end)}) is ${(a.best_other.ident * 100).toFixed(1)}% identical with ${a.best_other.fs} frameshift${a.best_other.fs === 1 ? "" : "s"}; the gene is whole in ${elsewhere.join(" and ") || "no other assembly"}.${reads}`;
    return `Missing from xgElyCris1.1${a.best_other && a.best_other.paralogue_of ? ` (the nearest match there is the paralogue ${G[a.best_other.paralogue_of] ? G[a.best_other.paralogue_of].label : a.best_other.paralogue_of})` : ""}. It is missing from the assembly, not from the animal: the gene is whole in ${elsewhere.join(" and ") || "no other assembly"}.${reads}`;
  }
  function exonFigure(g) {
    // draw the locus where the gene is most complete: Sanger if intact, otherwise the first complete assembly
    let key = "sanger";
    if (!["complete", "complete_discrepancy"].includes(g.assemblies.sanger.status)) key = ["purdue", "tampa", "alt"].find((k) => g.assemblies[k].status === "complete") || "sanger";
    const pcs = g.assemblies[key].pieces;
    if (!pcs.length) return null;
    const p = pcs[0], W = 760, X0 = 10, PW = W - 20, span = p.end - p.start + 1;
    const x = (bp) => X0 + ((p.strand === "-" ? p.end - bp : bp - p.start) / span) * PW;
    const svg = sv("svg", { viewBox: `0 0 ${W} 58`, class: "exons", role: "img", "aria-label": `Exon structure of ${g.label}: ${p.exons.length} coding exons over ${fmt(span)} bp` });
    svg.append(sv("line", { class: "gene-line", x1: X0, x2: X0 + PW, y1: 22, y2: 22 }));
    p.exons.forEach((e) => { const a = x(e[0]), b = x(e[1]); svg.append(sv("rect", { class: "exon", x: Math.min(a, b), y: 14, width: Math.max(1.5, Math.abs(b - a)), height: 16, rx: 1.5 })); });
    svg.append(sv("text", { class: "lab", x: X0, y: 50, text: `${chrName(p.chrom)}:${fmt(p.strand === "-" ? p.end : p.start)}` }));
    svg.append(sv("text", { class: "lab", x: X0 + PW, y: 50, "text-anchor": "end", text: fmt(p.strand === "-" ? p.start : p.end) }));
    svg.append(sv("text", { class: "lab", x: X0 + PW / 2, y: 50, "text-anchor": "middle", text: `${p.exons.length} coding exons over ${fmt(span)} bp, drawn 5′ to 3′` }));
    return el("figure", { class: "fig" }, [svg, el("figcaption", { text: `Coding exons in ${D.assemblies[key]}${key !== "sanger" ? " (the gene is not whole in xgElyCris1.1)" : ""}, drawn 5′ to 3′.` })]);
  }
  function domainFigure(g) {
    const W = 760, X0 = 10, PW = W - 20, L = g.length;
    const x = (aa) => X0 + (aa / L) * PW;
    const rows = g.domains.length ? 1 : 0;
    const svg = sv("svg", { viewBox: `0 0 ${W} 64`, class: "domains", role: "img", "aria-label": `Protein domains of ${g.label}: ${g.domains.map((d) => d.name).join(", ") || "none detected"}` });
    svg.append(sv("rect", { class: "backbone", x: X0, y: 22, width: PW, height: 8, rx: 4 }));
    let lastEnd = -1e9, up = true;
    g.domains.forEach((d, i) => {
      const a = x(d.start), b = x(d.end);
      const r = sv("rect", { class: "dom", x: a, y: 16, width: Math.max(3, b - a), height: 20, rx: 4, tabindex: 0 });
      tipFor(r, () => `<b>${d.name}</b> (${d.acc})<br>residues ${d.start}–${d.end}<br>i-E-value ${d.ievalue.toExponential(1)}`);
      svg.append(r);
      const tw = textW(d.name) * 0.9, mid = (a + b) / 2;
      up = mid - tw / 2 < lastEnd + 6 ? !up : true;
      svg.append(sv("text", { x: Math.min(Math.max(mid, X0 + tw / 2), X0 + PW - tw / 2), y: up ? 11 : 50, "text-anchor": "middle", text: d.name }));
      lastEnd = mid + tw / 2;
    });
    svg.append(sv("text", { class: "tick", x: X0, y: 62, text: "1" }));
    svg.append(sv("text", { class: "tick", x: X0 + PW, y: 62, "text-anchor": "end", text: fmt(L) + " aa" }));
    return el("figure", { class: "fig" }, [svg, el("figcaption", { text: g.domains.length ? "Pfam domains (hmmsearch, gathering thresholds)." : "No Pfam domain above the gathering threshold." })]);
  }
  function depthFigure(g) {
    const rs = g.reads_support; if (!rs) return null;
    const n = rs.per_segment.length, W = 760, X0 = 34, PW = W - X0 - 10, H = 120, base = 96;
    const maxv = Math.max(D.meta.typical_exon_depth * 1.6, ...rs.per_segment.map((s) => s.reads));
    const bw = Math.min(26, PW / n - 3);
    const svg = sv("svg", { viewBox: `0 0 ${W} ${H}`, class: "depthbars", role: "img", "aria-label": `Reads covering each of ${n} exon segments of ${g.label}` });
    const yy = (v) => base - (v / maxv) * 78;
    svg.append(sv("line", { class: "ref", x1: X0, x2: X0 + PW, y1: yy(D.meta.typical_exon_depth), y2: yy(D.meta.typical_exon_depth) }));
    svg.append(sv("text", { x: X0 + PW, y: yy(D.meta.typical_exon_depth) - 4, "text-anchor": "end", text: `typical clock-gene exon: ${Math.round(D.meta.typical_exon_depth)} reads` }));
    svg.append(sv("text", { x: X0 - 6, y: base + 3, "text-anchor": "end", text: "0" }));
    svg.append(sv("text", { x: X0 - 6, y: yy(maxv) + 9, "text-anchor": "end", text: Math.round(maxv) }));
    rs.per_segment.forEach((s, i) => {
      const x = X0 + (i + 0.5) * (PW / n) - bw / 2, h = Math.max(1.5, base - yy(s.reads));
      const r = sv("rect", { class: s.reads === 0 ? "zero" : "", x, y: s.reads === 0 ? base - 3 : base - h, width: bw, height: s.reads === 0 ? 3 : h, rx: 2, tabindex: 0 });
      tipFor(r, () => `<b>${s.reads} reads</b><br>exon segment ${i + 1} of ${n}<br>${s.segment}`);
      svg.append(r);
    });
    svg.append(sv("text", { x: X0, y: H - 6, text: "exon segments, 5′ → 3′ of the Purdue-assembly gene model" }));
    return el("figure", { class: "fig" }, [svg, el("figcaption", { text: `HiFi reads from the sequenced animal (xgElyCris1) overlapping each exon by at least 150 bp. ${rs.covered} of ${rs.segments} segments have three or more reads.` })]);
  }
  function fasta(name, desc, seq) { return `>${name} ${desc}\n${seq.replace(/(.{60})/g, "$1\n")}\n`; }
  function seqBlock(title, name, desc, seq) {
    const text = fasta(name, desc, seq);
    const copy = el("button", { class: "btn", type: "button", text: "Copy FASTA" });
    copy.addEventListener("click", () => navigator.clipboard.writeText(text).then(() => { copy.textContent = "Copied"; setTimeout(() => (copy.textContent = "Copy FASTA"), 1400); }));
    const dl = el("a", { class: "btn", download: `${name}.${title.startsWith("Protein") ? "faa" : "fna"}`, href: "data:text/plain;charset=utf-8," + encodeURIComponent(text), text: "Download" });
    return el("div", {}, [el("div", { class: "seqtools" }, [el("h4", { text: title }), copy, dl]), el("div", { class: "seqbox", tabindex: 0, text: seq })]);
  }
  function openGene(sym, push) {
    const g = G[sym]; if (!g || !dlg) return;
    hideTip();
    const A = g.assemblies;
    const asmRows = ["sanger", "alt", "purdue", "tampa"].map((k) => {
      const a = A[k];
      const where = a.pieces.length ? a.pieces.map((p) => `${p.chrom}:${fmt(p.start)}–${fmt(p.end)} (${p.strand}) aa ${p.qstart + 1}–${p.qend}, ${(p.ident * 100).toFixed(1)}% id${p.fs ? `, ${p.fs} frameshift${p.fs > 1 ? "s" : ""}` : ""}${p.stop ? `, ${p.stop} stop` : ""}`).join("\n")
        : (a.best_other ? `closest match ${a.best_other.chrom}:${fmt(a.best_other.start)}–${fmt(a.best_other.end)}, ${(a.best_other.ident * 100).toFixed(1)}% id${a.best_other.paralogue_of ? " (paralogue " + a.best_other.paralogue_of + ")" : ""}` : "no match");
      return el("tr", {}, [el("td", { text: D.assemblies[k] }), el("td", {}, [statusHtml(g, k)]), el("td", { class: "mono", style: "white-space:pre-line", text: where })]);
    });
    const orth = g.orthologues.map((o) => el("tr", {}, [el("td", { html: `<i>${D.species[o.species] || o.species}</i>` }), el("td", { class: "mono", text: o.protein }), el("td", { text: o.gene || "" }), el("td", { class: "num", text: fmt(o.length) }), el("td", { class: "num", text: o.pident.toFixed(0) + "%" }), el("td", { text: o.reciprocal ? "yes" : "no" })]));
    const cands = g.candidates.slice(0, 8).map((c) => el("tr", {}, [el("td", { class: "mono", text: c.candidate }), el("td", { class: "num", text: fmt(c.length) }), el("td", { class: "num", text: c.score.toFixed(0) }), el("td", { text: c.chosen ? "used" : (c.qc || "") })]));
    const tp = g.tree_placement;
    const treeText = tp ? `In the ${tp.tree} family tree ${g.label} groups with ${tp.seeds_in_clade.join(", ")}${tp.support !== null && tp.support !== undefined ? ` (clade of ${tp.clade_size} sequences, support ${tp.support > 1 ? Math.round(tp.support) + "%" : tp.support.toFixed(2)})` : ` (clade of ${tp.clade_size} sequences)`}.${tp.support !== null && tp.support !== undefined && (tp.support > 1 ? tp.support / 100 : tp.support) < 0.7 ? " That support is weak, so the name rests mainly on reciprocal best hits." : ""}` : "";
    const ids = [];
    if (g.ensembl_overlap.length) ids.push("Ensembl gene" + (g.ensembl_overlap.length > 1 ? "s " : " ") + g.ensembl_overlap.join(", "));
    if (g.purdue_proteins.length) ids.push("Purdue " + g.purdue_proteins.map((p) => `${p.locus_tag} (${p.accession})`).join(", "));
    const body = el("div", { class: "gene-body" }, [
      el("section", {}, [el("p", { class: "verdict", text: sentence(g) }), g.note ? el("p", { class: "muted", text: g.note.charAt(0).toUpperCase() + g.note.slice(1) + "." }) : null, ids.length ? el("p", { class: "small muted", text: ids.join("; ") }) : null]),
      el("section", {}, [el("h3", { text: "Where it is" }), exonFigure(g),
        el("table", { class: "kv" }, [el("thead", {}, [el("tr", {}, [el("th", { text: "Assembly" }), el("th", { text: "State" }), el("th", { text: "Location of the protein-coding alignment" })])]), el("tbody", {}, asmRows)])]),
      el("section", {}, [el("h3", { text: "The protein" }), el("p", { class: "small muted", text: g.model_description + "." }), domainFigure(g),
        seqBlock(`Protein, ${fmt(g.length)} aa`, "Ecr_" + g.symbol, g.name + " [Elysia crispata] " + g.model_source, g.protein),
        g.cds ? seqBlock(`Coding sequence, ${fmt(g.cds.length)} nt`, "Ecr_" + g.symbol + "_cds", g.name + " [Elysia crispata] " + g.model_source, g.cds) : el("p", { class: "small muted", text: "No coding sequence is given for this model." })]),
      el("section", {}, [el("h3", { text: "Why it is called " + g.label }), el("p", { text: treeText }),
        D.trees[g.symbol] ? el("figure", { class: "fig" }, [el("div", { html: D.trees[g.symbol] }), el("figcaption", { text: "Blue: this gene. Squares: reference proteins that define the family. Numbers are branch supports of 0.7 or more. The dashed branch is the neighbouring lineage, collapsed." })]) : null,
        el("h4", { text: "Best reciprocal hits in other species", style: "margin:16px 0 8px" }),
        el("table", { class: "kv" }, [el("thead", {}, [el("tr", {}, [el("th", { text: "Species" }), el("th", { text: "Protein" }), el("th", { text: "Gene" }), el("th", { class: "num", text: "Length" }), el("th", { class: "num", text: "Identity" }), el("th", { text: "Reciprocal" })])]), el("tbody", {}, orth)])]),
      el("section", {}, [el("h3", { text: "Evidence in the sequenced animal" }), depthFigure(g),
        g.expression ? el("p", { text: `Head RNA-seq (one sample, one time point): ${g.expression.tpm.toFixed(1)} TPM${g.expression.per_locus.length > 1 ? " summed over " + g.expression.per_locus.map((x) => x.locus_tag + " " + x.tpm.toFixed(1)).join(", ") : ""}. A single library says the gene is transcribed; it says nothing about rhythm.` }) : el("p", { class: "muted", text: "No expression estimate." })]),
      el("section", {}, [el("h3", { text: "Gene models compared" }), el("p", { class: "small muted", text: "Score is the mean DIAMOND bit score against orthologues in four outgroup molluscs; higher means more of the conserved protein is captured." }),
        el("table", { class: "kv" }, [el("thead", {}, [el("tr", {}, [el("th", { text: "Candidate model" }), el("th", { class: "num", text: "aa" }), el("th", { class: "num", text: "Score" }), el("th", { text: "" })])]), el("tbody", {}, cands)])]),
    ]);
    const head = el("div", { class: "gene-head" }, [el("div", {}, [el("h2", { text: g.label }), el("div", { class: "gname", text: g.name }), el("div", { class: "gcat", text: CATS[g.category] || g.category })]),
      el("button", { class: "close", type: "button", text: "Close", onclick: () => dlg.close() })]);
    dlg.replaceChildren(head, body);
    if (!dlg.open) dlg.showModal();
    dlg.scrollTop = 0;
    if (push !== false) history.replaceState(null, "", "#gene/" + encodeURIComponent(sym));
  }
  if (dlg) {
    dlg.addEventListener("click", (e) => { if (e.target === dlg) dlg.close(); });
    dlg.addEventListener("close", () => { if (location.hash.startsWith("#gene/")) history.replaceState(null, "", location.pathname + location.search); });
  }
  window.openGene = openGene;

  function init() {
    renderSchematic(); renderKaryo(); renderFilters(); renderTable();
    const m = location.hash.match(/^#gene\/(.+)$/);
    if (m) openGene(decodeURIComponent(m[1]), false);
  }
  (document.fonts && document.fonts.ready ? document.fonts.ready : Promise.resolve()).then(init);
})();
