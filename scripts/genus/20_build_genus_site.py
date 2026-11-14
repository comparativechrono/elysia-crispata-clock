#!/usr/bin/env python3
"""Build the genus-level pages of the site and the clock-gene page of every species except
E. crispata (whose pages are written by scripts/30_build_site.py).

usage: 20_build_genus_site.py repo_dir
Reads  config/species.tsv, config/species_catalogue.tsv, results/clock_genes.json,
       results/genus/<slug>.clock.json, results/genus/genus.json, docs/assets/photos/credits.json
       and, when present, results/genus/<slug>.kleptoplast.json / <slug>.transfer.json
Writes docs/index.html (species), docs/compare.html, docs/methods.html, docs/<slug>/index.html,
       docs/data/<slug>.js, docs/data/genus.js, docs/downloads/*, and forwarding pages for the
       addresses the E. crispata pages had before the site covered more than one species.
"""
import sys, os, json, html, datetime, shutil, csv, collections, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chrome
from chrome import LABEL, CATS, STATUS_CLASS, STATUS_WORD

repo = sys.argv[1]
R = os.path.join(repo, "results"); RG = os.path.join(R, "genus"); SITE = os.path.join(repo, "docs")
os.makedirs(os.path.join(SITE, "data"), exist_ok=True); os.makedirs(os.path.join(SITE, "downloads"), exist_ok=True)
e = html.escape; today = datetime.date.today().isoformat()
def it(s): return f"<i>{e(s)}</i>"
def num(n): return f"{n:,}"
def jload(p, default=None): return json.load(open(p)) if os.path.exists(p) else default

SP = list(csv.DictReader(open(os.path.join(repo, "config", "species.tsv")), delimiter="\t"))
BYSLUG = {s["slug"]: s for s in SP}
PRE = {"elysia-crispata": "Ecr", "elysia-chlorotica": "Ech", "elysia-marginata": "Ema", "elysia-timida": "Eti", "elysia-leucolegnote": "Ele",
       "elysia-atroviridis": "Eat", "elysia-papillosa": "Epa"}
SHORT = {s["slug"]: "E. " + s["species"].split()[1] for s in SP}
crisp = json.load(open(os.path.join(R, "clock_genes.json"))); CG = {g["symbol"]: g for g in crisp["genes"]}
symbols = [g["symbol"] for g in crisp["genes"]]
for g in crisp["genes"]: g["label"] = LABEL.get(g["symbol"], g["symbol"])
GN = json.load(open(os.path.join(RG, "genus.json"))); GG = GN["genes"]
CAT = {s["slug"]: jload(os.path.join(RG, s["slug"] + ".clock.json")) for s in SP if s["slug"] != "elysia-crispata"}
CAT = {k: v for k, v in CAT.items() if v}
credits = jload(os.path.join(SITE, "assets", "photos", "credits.json"), {})
NOTES = {}
np_ = os.path.join(repo, "config", "genus_notes.tsv")
if os.path.exists(np_):
    for l in open(np_):
        if l.strip() and not l.startswith("#"):
            k, v = l.rstrip("\n").split("\t", 1); NOTES[k] = v
KLEPTO = {s["slug"]: jload(os.path.join(RG, s["slug"] + ".kleptoplast.json")) for s in SP}
TRANSFER = {s["slug"]: jload(os.path.join(RG, s["slug"] + ".transfer.json")) for s in SP}
ck = jload(os.path.join(R, "kleptoplast.json"))
STEXT = jload(os.path.join(repo, "config", "species_text.json"), {})

def pair(a, b): return "|".join(sorted((a, b)))
def status_of(slug, sym):
    if slug == "elysia-crispata": return CG[sym]["assemblies"]["sanger"]["status"]
    return next(g for g in CAT[slug]["genes"] if g["symbol"] == sym)["status"]
def gene_of(slug, sym): return next(g for g in CAT[slug]["genes"] if g["symbol"] == sym)
def counts(slug):
    c = collections.Counter(STATUS_CLASS[status_of(slug, s)] for s in symbols)
    return c
def statusbar(slug):
    c = counts(slug); n = len(symbols)
    lab = f'{c["ok"]} intact, {c["warn"]} whole but with a break in the sequence, {c["part"]} partial or in pieces, {c["miss"]} not found'
    return f'<div class="statusbar" role="img" aria-label="{lab}" title="{lab}">' + "".join(f'<i class="{k}" style="flex:{c[k]}"></i>' for k in ("ok", "warn", "part", "miss") if c[k]) + "</div>"
def tabs_for(slug):
    t = [("index.html", "Clock genes")]
    # the pages of the other species are written by 21_build_species_pages.py, which runs first
    if slug == "elysia-crispata" or os.path.exists(os.path.join(SITE, slug, "kleptoplasts.html")): t.append(("kleptoplasts.html", "Whose chloroplasts"))
    if slug == "elysia-crispata" or os.path.exists(os.path.join(SITE, slug, "gene-transfer.html")): t.append(("gene-transfer.html", "Gene transfer"))
    if slug == "elysia-crispata": t.append(("methods.html", "Methods for this species"))
    return t

# ------------------------------------------------------------------ per-species clock pages
LEVEL_NOTE = {"chromosome": "chromosome-level", "scaffold": "scaffold-level", "contig": "contig-level"}
CAVEAT = {
 "elysia-chlorotica": "The assembly is in 9,989 scaffolds built from short reads and early PacBio reads. Genes that come out partial or with a broken reading frame here are more likely casualties of the assembly than of evolution.",
 "elysia-marginata": "The assembly was built from short reads alone, in 14,149 scaffolds with 82 Mb of gaps. Expect genes to be split or to end at a gap.",
 "elysia-timida": "A chromosome-level assembly. Its gene annotation is not in GenBank, so every model here is the <i>E. crispata</i> protein aligned to the genome.",
 "elysia-leucolegnote": "The public sequence of this assembly is 35% N: it looks hard-masked, with repeats blanked out. An exon that falls in a masked stretch cannot be found, so partial genes here may be whole in the animal.",
 "elysia-atroviridis": "A chromosome-level assembly from accurate long reads. Its gene annotation is not public, so every model here is the <i>E. crispata</i> protein aligned to the genome.",
 "elysia-papillosa": "A contig-level assembly from accurate long reads, not yet scaffolded into chromosomes and without a gene annotation. Every model here is the <i>E. crispata</i> protein aligned to the genome.",
}
def species_data(slug):
    c = CAT[slug]; sp = BYSLUG[slug]; A = c["assembly"]
    genes = []
    for g in c["genes"]:
        m = CG[g["symbol"]]; sym = g["symbol"]; lab = LABEL.get(sym, sym)
        ident = g.get("identity_to_crispata")
        asm = dict(status=g["status"], coverage=g.get("consensus_coverage", 0), pieces=g.get("pieces", []))
        if g["status"] == "absent" and g.get("nearest"):
            n = g["nearest"]; asm["best_other"] = dict(chrom=n["chrom"], start=n["start"], end=n["end"], ident=n["identity"], fs=0,
                                                       paralogue_of=n["looks_like"] if n.get("looks_like") in CG and n["looks_like"] != sym else None)
        if g["status"] == "absent":
            ident_txt = "No place in the genome passes the reciprocal test for this gene."
        elif g["model_kind"] == "annotation":
            ident_txt = (f"This is the species’ own annotated protein, GenBank {g['protein_id']}. It passes the reciprocal test: among all E. crispata proteins its best match is {lab}. "
                         f"It is {round(100 * ident) if ident else '?'}% identical to the E. crispata protein and fills {round(100 * g['consensus_coverage'])}% of the positions that most species share.")
        else:
            ident_txt = (f"Found by aligning the E. crispata {lab} protein to the genome. The place it landed passes the reciprocal test: translated, its best match among all E. crispata proteins is {lab} again. "
                         f"The protein is {round(100 * ident) if ident else '?'}% identical to the E. crispata one and fills {round(100 * g['consensus_coverage'])}% of the positions that most species share."
                         + (" The species’ annotation has a gene here too, but the aligned model covers more of the protein." if g.get("annotation_protein") else ""))
        if g.get("extra_loci"):
            ident_txt += f" {len(g['extra_loci'])} more place{'s' if len(g['extra_loci']) > 1 else ''} in the assembly also pass{'' if len(g['extra_loci']) > 1 else 'es'} the test ({'; '.join(x['chrom'] + ':' + num(x['start']) + '–' + num(x['end']) for x in g['extra_loci'][:3])}): a second copy, or the other haplotype left in the assembly."
        if sym in NOTES: ident_txt += " " + NOTES[sym]
        genes.append(dict(symbol=sym, label=lab, name=m["name"], category=m["category"], note=m.get("note", ""), protein=g.get("protein", ""), length=g.get("length", 0), cds=g.get("cds", ""),
                          model_kind=g.get("model_kind"), model_source=g.get("model_source", ""), model_description=(g.get("model_source") or "No model"),
                          domains=g.get("domains", []), assemblies={"primary": asm}, identity_to_crispata=ident, identified=ident_txt, long_intron=g.get("long_intron", False),
                          ids=[f"GenBank protein {g['protein_id']}"] if g.get("protein_id") else []))
    across = {}
    for sym in symbols:
        rows = []
        for o in SP:
            if o["slug"] == slug or (o["slug"] != "elysia-crispata" and o["slug"] not in CAT): continue
            st = status_of(o["slug"], sym)
            ln = len(CG[sym]["protein"]) if o["slug"] == "elysia-crispata" else gene_of(o["slug"], sym).get("length", 0)
            rows.append(dict(name=o["species"], href=f"../{o['slug']}/index.html#gene/{sym}", status=st, length=ln, identity=GG[sym]["identity"].get(pair(o["tag"], sp["tag"]))))
        across[sym] = rows
    return dict(genes=genes, not_found=crisp["not_found"], assemblies={"primary": f"{A['name']} ({A['accession']})"}, primary="primary", primary_name=A["name"], others=[],
                species_name=sp["species"], prefix=PRE[slug], columns=["gene", "role", "aa", "status", "ident", "model"], categories=CATS, species={},
                chromosomes=[dict(name=x["name"], length=x["length"]) for x in A["chromosomes"]], trees={}, meta=dict(built=today, typical_exon_depth=0),
                across=across, compare_href="../compare.html")

WORDS = ["no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve"]
def word(k): return WORDS[k] if 0 <= k < len(WORDS) else str(k)
def names(gs): return ", ".join(LABEL.get(g["symbol"], g["symbol"]) for g in gs) or "none"
for slug, c in CAT.items():
    sp = BYSLUG[slug]; A = c["assembly"]; data = species_data(slug)
    with open(os.path.join(SITE, "data", slug + ".js"), "w") as fo: fo.write("window.CLOCK = " + json.dumps(data, separators=(",", ":")) + ";\n")
    G = c["genes"]; n = len(G)
    ok = [g for g in G if g["status"] == "complete"]; disc = [g for g in G if g["status"] == "complete_discrepancy"]
    part = [g for g in G if g["status"] in ("partial", "fragmented")]; miss = [g for g in G if g["status"] == "absent"]
    n_ann = sum(1 for g in G if g.get("model_kind") == "annotation"); n_proj = sum(1 for g in G if g.get("model_kind") == "projection")
    ids = [g["identity_to_crispata"] for g in G if g.get("identity_to_crispata")]
    nf_found = [x for x in c.get("not_found", []) if x["found"]]
    chrom = bool(A["chromosomes"])
    how = (f"{n_ann} models are the species’ own annotated genes and {n_proj} are the <i>E. crispata</i> protein aligned to the genome, used where the annotation has no gene or a shorter one."
           if n_ann else f"The assembly has no public gene annotation, so all {n_proj} models are the <i>E. crispata</i> protein aligned to the genome.")
    body = f"""
<section class="intro"><div class="wrap">
<h1 style="max-width:18ch">The clock genes of {it(sp["species"])}</h1>
<p class="lede">Where each of the 47 clock genes catalogued in <i>E. crispata</i> sits in the {LEVEL_NOTE[sp["level"]]} assembly {e(A["name"])} ({A["span"] / 1e6:,.0f} Mb), whether the assembly has it whole, and the protein it encodes. The sequenced animal came from {e(sp["origin"])}.</p>
<div class="clockwork">
<div id="schematic"></div>
<div class="clockwork-list" id="schematic-list"></div>
<div class="clockwork-legend" aria-label="Key to gene status in {e(A["name"])}">
<span><i class="glyph ok"></i>intact in {e(A["name"])}</span>
<span><i class="glyph warn"></i>whole, but the sequence does not translate cleanly</span>
<span><i class="glyph part"></i>partial or in pieces</span>
<span><i class="glyph miss"></i>not found</span>
</div>
</div>
<p class="clockwork-caveat">Select a gene to open it. The wiring is the textbook fly and mouse clock; which gene does what in <i>Elysia</i> has not been tested.</p>
</div></section>

<section class="block"><div class="wrap">
<div class="cols">
<div><h3>What was found</h3>
<p>{len(ok)} of the {n} genes are intact{f" and {word(len(disc))} more {'is' if len(disc) == 1 else 'are'} whole but {'needs' if len(disc) == 1 else 'need'} a frameshift or stop codon to read through ({names(disc)})" if disc else ""}. {f"{word(len(part)).capitalize()} {'is' if len(part) == 1 else 'are'} partial or in pieces ({names(part)})" if part else "None is partial"}{f", and {word(len(miss))} {'has' if len(miss) == 1 else 'have'} no counterpart ({names(miss)})" if miss else ""}.</p></div>
<div><h3>How the genes were found</h3>
<p>{how} Every model had to pass a reciprocal test against all <i>E. crispata</i> proteins, which keeps paralogues apart. The proteins are {round(100 * min(ids))} to {round(100 * max(ids))}% identical to their <i>E. crispata</i> counterparts (median {round(100 * sorted(ids)[len(ids) // 2])}%).</p></div>
<div><h3>What to keep in mind</h3>
<p>{CAVEAT.get(slug, "")}</p></div>
</div>
</div></section>
""" + (f"""
<section class="block"><div class="wrap">
<h2>Where the genes sit</h2>
<p>The {len(A["chromosomes"])} chromosomes of {e(A["name"])}, drawn to scale and numbered as in the assembly. Chromosome numbers are not comparable between species. Each mark is a clock gene; select one to open it.</p>
<div class="karyo-card"><div class="karyo" id="karyo"></div><div class="offmap" id="offmap"></div></div>
</div></section>
""" if chrom else "") + f"""
<section class="block"><div class="wrap">
<h2>All {n} genes</h2>
<div class="filters" id="filters"></div>
<div class="tablewrap"><table class="genes">
<thead><tr><th>Gene</th><th>Role</th><th class="num">Protein, aa</th><th>In {e(A["name"])}</th><th class="num" title="Identical residues over the aligned part of the two proteins">Identity to <i>E. crispata</i></th><th>Model</th></tr></thead>
<tbody id="genes-body"></tbody></table></div>
<p class="small muted" style="margin-top:12px;max-width:90ch">A model is called intact when it fills at least 90% of the alignment positions that most species share and the genome reads through without a frameshift or stop. Models aligned from <i>E. crispata</i> can miss poorly conserved ends of a protein.</p>

<h3 style="margin-top:34px">Looked for and not found</h3>
<p style="max-width:80ch">The seven clock-associated genes that <i>E. crispata</i> lacks were searched for here too, with the same fly and mammal proteins. {"None has a counterpart in this genome either." if not nf_found else "Found here: " + ", ".join(x["symbol"] for x in nf_found) + "."}</p>
</div></section>
<dialog class="gene" id="gene" aria-label="Gene details"></dialog>
"""
    os.makedirs(os.path.join(SITE, slug), exist_ok=True)
    src = f"{it(sp['species'])}: genome assembly {e(A['name'])} ({e(A['accession'])}). Sources and papers for every species are listed on the <a href=\"../methods.html\">methods page</a>."
    open(os.path.join(SITE, slug, "index.html"), "w").write(chrome.page(
        f"{sp['species']} clock genes", body, root="../", species=sp["species"], tabs=tabs_for(slug), tab_current="index.html", today=today, sources=src,
        description=f"The circadian clock genes of the sea slug {sp['species']}, located in the genome assembly {A['name']}.",
        scripts=f'<script src="../data/{slug}.js"></script>\n<script src="../assets/app.js"></script>'))
    # downloads
    with open(os.path.join(SITE, "downloads", f"{slug}_clock_proteins.faa"), "w") as fo:
        for g in G:
            if g.get("protein"): fo.write(f">{PRE[slug]}_{g['symbol']} {CG[g['symbol']]['name']} [{sp['species']}] {g['status']} | {g['model_source']}\n" + re.sub(r"(.{60})", r"\1\n", g["protein"]).rstrip("\n") + "\n")
    with open(os.path.join(SITE, "downloads", f"{slug}_clock_cds.fna"), "w") as fo:
        for g in G:
            if g.get("cds"): fo.write(f">{PRE[slug]}_{g['symbol']}_cds {CG[g['symbol']]['name']} [{sp['species']}] {g['status']}\n" + re.sub(r"(.{60})", r"\1\n", g["cds"]).rstrip("\n") + "\n")

# ------------------------------------------------------------------ data for the comparison page
order = [s for s in SP if s["slug"] == "elysia-crispata" or s["slug"] in CAT]
def chrom_of(slug, sym):
    if slug == "elysia-crispata":
        p = CG[sym]["assemblies"]["sanger"]["pieces"]; return p[0]["chrom"] if p and re.fullmatch(r"\d+", p[0]["chrom"]) else None
    p = gene_of(slug, sym).get("pieces") or []
    return p[0]["chrom"] if p and re.fullmatch(r"\d+", str(p[0]["chrom"])) else None
chromsp = [s for s in order if s["level"] == "chromosome"]
cmp_genes = []
for sym in symbols:
    gg = GG[sym]
    st = {s["slug"]: status_of(s["slug"], sym) for s in order}
    comp = {s["slug"]: gg["completeness"].get(s["tag"]) for s in order}
    copies = {s["slug"]: (1 + len(gene_of(s["slug"], sym).get("extra_loci", [])) if s["slug"] in CAT and st[s["slug"]] != "absent" else None) for s in order}
    cmp_genes.append(dict(symbol=sym, label=LABEL.get(sym, sym), name=CG[sym]["name"], category=CG[sym]["category"], status=st, completeness=comp, copies=copies,
                          identity=gg["identity"], mean_identity=gg["mean_identity_elysia"], min_identity=gg["min_identity_elysia"], n_compared=gg["n_elysia_compared"],
                          alignment=gg["alignment"], chrom={s["slug"]: chrom_of(s["slug"], sym) for s in chromsp}, note=NOTES.get(sym, "")))
OUTNAMES = {"PLAOC": "Plakobranchus ocellatus", "APLCA": "Aplysia californica", "BIOGL": "Biomphalaria glabrata"}
genus = dict(species=[dict(slug=s["slug"], name=s["species"], short=SHORT[s["slug"]], tag=s["tag"], level=s["level"]) for s in order], outgroups=OUTNAMES, categories=CATS,
             genes=cmp_genes, tree=GN["tree"], tree_genes=GN["tree_genes"], tree_columns=GN["tree_columns"], chromosome_species=[s["slug"] for s in chromsp])
with open(os.path.join(SITE, "data", "genus.js"), "w") as fo: fo.write("window.GENUS = " + json.dumps(genus, separators=(",", ":")) + ";\n")

# linkage: are genes that share a chromosome in one species together in the others?
pairs_total = pairs_kept = 0
if len(chromsp) > 1:
    ref = chromsp[0]["slug"]
    loc = {sym: {s["slug"]: chrom_of(s["slug"], sym) for s in chromsp} for sym in symbols}
    full = [sym for sym in symbols if all(loc[sym][s["slug"]] for s in chromsp)]
    for i, a in enumerate(full):
        for b in full[i + 1:]:
            if loc[a][ref] == loc[b][ref]:
                pairs_total += 1; pairs_kept += all(loc[a][s["slug"]] == loc[b][s["slug"]] for s in chromsp)
else: full = []

# ------------------------------------------------------------------ comparison page
all_intact = [sym for sym in symbols if all(STATUS_CLASS[status_of(s["slug"], sym)] in ("ok", "warn") for s in order)]
mi = sorted((g for g in cmp_genes if g["mean_identity"]), key=lambda g: g["mean_identity"])
fast = mi[:5]; slow = mi[-5:][::-1]
core = [g for g in cmp_genes if g["category"] == "core" and g["mean_identity"]]
CORE = ["Clock", "Bmal1", "Period", "Timeless", "Cry1", "Cry2"]
core_whole = all(STATUS_CLASS[status_of(s["slug"], sym)] in ("ok", "warn") for s in order for sym in CORE if s["slug"] != "elysia-crispata")
CORE_TXT = ("Every one of the other species has both halves of the activator (Clock and Bmal1), the repressors Period, Timeless and a mammalian-type Cryptochrome, and a light-sensitive Cryptochrome, each as a whole gene."
            if core_whole else "The core loop (Clock, Bmal1, Period, Timeless and two Cryptochromes) is present in every species, though not every gene is whole in every assembly.")
bad = sorted(((sum(1 for sym in symbols if STATUS_CLASS[status_of(s["slug"], sym)] != "ok"), s) for s in order if s["slug"] != "elysia-crispata"), key=lambda x: -x[0])
WORST_TXT = "" if {bad[0][1]["slug"], bad[1][1]["slug"]} != {"elysia-chlorotica", "elysia-marginata"} else (f"The two assemblies built wholly or mainly from short reads have the most genes that are not intact ({it(bad[0][1]['species'])} {bad[0][0]}, {it(bad[1][1]['species'])} {bad[1][0]}); "
             f"the four built from accurate long reads have {min(b[0] for b in bad[2:])} to {max(b[0] for b in bad[2:])} each.")
# the core loop in E. crispata: whole in the reference, or only in another assembly of the species
cr_ref_gap = [sym for sym in CORE if STATUS_CLASS[status_of("elysia-crispata", sym)] not in ("ok", "warn")]
cr_any = all(any(v.get("status") in ("complete", "complete_discrepancy") for v in CG[sym]["assemblies"].values()) for sym in CORE)
HOME_CORE = ""
if core_whole and cr_any:
    HOME_CORE = "Every species has the whole core loop: Clock, Bmal1, Period, Timeless and two Cryptochromes"
    HOME_CORE += (f" (in <i>E. crispata</i>, {' and '.join(LABEL.get(x, x) for x in cr_ref_gap)} {'is' if len(cr_ref_gap) == 1 else 'are'} whole in another assembly of the species, not the reference)." if cr_ref_gap else ".")
compare_body = f"""
<section class="page-head"><div class="wrap">
<h1>The same clock in {word(len(order))} sea slugs</h1>
<p class="lede">The 47 clock genes of <i>E. crispata</i>, looked for in the genome of every <i>Elysia</i> species that has one, and their proteins laid side by side.</p>
<div class="answer">
<p><strong>The clock is the same throughout the genus.</strong> {len(all_intact)} of the 47 genes are whole in the reference assembly of all {word(len(order))} species. {CORE_TXT}</p>
<p><strong>Where a gene is missing or partial, the assembly is usually the reason.</strong> {WORST_TXT} Two genes are exceptions: the Go-opsins, whose number differs from species to species, and r-opsin-like 1, whose model is incomplete in <i>E. crispata</i> itself.</p>
<p><strong>The proteins differ in how fast they change.</strong> The most alike across the genus are {", ".join(g["label"] + " " + str(round(100 * g["mean_identity"])) + "%" for g in slow[:3])}; the least alike are {", ".join(g["label"] + " " + str(round(100 * g["mean_identity"])) + "%" for g in fast[:3])}. Among the core loop, {", ".join(g["label"] + " " + str(round(100 * g["mean_identity"])) + "%" for g in sorted(core, key=lambda g: -g["mean_identity"]))}.</p>
</div>
</div></section>

<section class="block"><div class="wrap">
<h2>Gene by species</h2>
<p>Each cell shows the state of the gene in that species’ reference assembly, and how identical its protein is to the <i>E. crispata</i> one. Select a gene name for the alignment, or a cell to open the gene in that species.</p>
<div class="legend"><span><i class="glyph ok"></i>intact</span><span><i class="glyph warn"></i>whole, sequence does not translate cleanly</span><span><i class="glyph part"></i>partial or in pieces</span><span><i class="glyph miss"></i>not found</span></div>
<div class="tablewrap"><table class="matrix" id="matrix"></table></div>
<p class="small muted" style="margin-top:12px;max-width:92ch">For <i>E. crispata</i> the state is that of the Sanger assembly xgElyCris1.1; all 47 genes are whole in at least one other assembly of that species. A model counts as intact when it fills at least 90% of the alignment positions that most sequences share.</p>
</div></section>

<section class="block"><div class="wrap">
<h2>How alike are the proteins?</h2>
<p>Mean identity between every pair of <i>Elysia</i> species, for each protein that is complete in at least two of them. Select a bar for the alignment.</p>
<div class="fig"><div id="cons"></div><p class="figcap">Identity is counted over the positions where both proteins have a residue. Opsins and photolyases are shown with the rest; low values there partly reflect paralogues that are hard to tell apart.</p></div>
</div></section>

<section class="block"><div class="wrap">
<h2>One protein, all species</h2>
<div class="gene-pick"><label for="pick">Protein</label><select id="pick"></select><span class="muted small" id="pick-note"></span></div>
<div id="idmatrix"></div>
<div class="aln-wrap" style="margin-top:14px"><div class="aln" id="aln"></div></div>
<p class="small muted" style="margin-top:10px;max-width:92ch">Shaded residues differ from the most common residue among the <i>Elysia</i> sequences at that position. The three sequences below the line are other gastropods, for orientation. Alignments are in <code>results/genus/alignments/</code> in the repository.</p>
</div></section>

<section class="block"><div class="wrap">
<h2>How the species are related</h2>
<p>A tree from {len(GN["tree_genes"])} clock proteins that are complete in every species, joined end to end ({num(GN["tree_columns"])} aligned positions), rooted on the sea hare <i>Aplysia</i>.</p>
<div class="fig"><div id="tree"></div><p class="figcap">Maximum likelihood (IQ-TREE, LG+F+G4). Numbers are ultrafast bootstrap supports. Branch lengths are substitutions per site. A tree from clock proteins alone is a guide to the species on this site, not a phylogeny of the genus.</p></div>
</div></section>
""" + (f"""
<section class="block"><div class="wrap">
<h2>Do the genes stay on the same chromosomes?</h2>
<p>{word(len(chromsp)).capitalize()} of the assemblies are chromosome-level, each with 15 chromosomes numbered by size, so the numbers themselves do not correspond. What can be compared is company: of the {pairs_total} pairs of clock genes that share a chromosome in <i>E. crispata</i>, {pairs_kept} share one in {"all " + word(len(chromsp) - 1) + " other species" if len(chromsp) > 2 else "the other species"} too.</p>
<div class="tablewrap"><table class="data" id="linkage"></table></div>
</div></section>
""" if len(chromsp) > 1 else "") + """
<script>window.COMPARE_ROOT = "";</script>
"""
open(os.path.join(SITE, "compare.html"), "w").write(chrome.page(
    "Clock genes compared across Elysia species", compare_body, root="", site_current="compare.html", today=today,
    scripts='<script src="data/genus.js"></script>\n<script src="assets/pages.js"></script>\n<script src="assets/compare.js"></script>'))

# ------------------------------------------------------------------ home page
def card(s):
    slug = s["slug"]; c = counts(slug); cr = credits.get(slug)
    A = CAT[slug]["assembly"] if slug in CAT else dict(name="xgElyCris1.1", accession="GCA_963854125.1", span=542.6e6)
    photo = (f'<a class="ph" href="{slug}/index.html"><img src="assets/photos/{cr["file"]}" alt="{e(s["species"])}" width="{cr["width"]}" height="{cr["height"]}" loading="lazy">'
             f'</a>' if cr else '<div class="ph"></div>')
    credit = f'<span class="credit">Photo: <a href="{cr["source"]}">{e(cr["photographer"])}</a>, {e(cr["licence"])}</span>' if cr else ""
    ann = {"genbank": "in GenBank", "ensembl": "Ensembl", "senckenberg": "not in GenBank", "none": "none public"}[s["annotation"]]
    k = KLEPTO.get(slug); kl = ""
    if slug == "elysia-crispata" and ck:
        top = ck["assignment"][:2]; kl = ", ".join(f'{it(a["taxon"])} {round(100 * a["fraction"])}%' for a in top)
    elif k and STEXT.get(slug, {}).get("card"): kl = STEXT[slug]["card"]
    links = "".join(f'<a href="{slug}/{h}">{t}</a>' for h, t in tabs_for(slug) if h != "methods.html")
    return f"""<article class="sp-card">
<div style="position:relative">{photo}{credit}</div>
<div class="tx">
<h3><a href="{slug}/index.html">{e(s["species"])}</a></h3>
<p class="where">Sequenced animal: {e(s["origin"])}</p>
<div><div class="small" style="margin-bottom:3px"><strong>{c["ok"]}</strong> of 47 clock genes intact in the assembly</div>{statusbar(slug)}</div>
<dl>
<dt>Assembly</dt><dd>{LEVEL_NOTE[s["level"]]}, {A["span"] / 1e6:,.0f} Mb <span class="muted small">{e(A["accession"])}</span></dd>
<dt>Gene annotation</dt><dd>{ann}</dd>
<dt>Long reads</dt><dd>{e(s["long_reads"])}</dd>
{f"<dt>Chloroplasts from</dt><dd>{kl}</dd>" if kl else ""}
</dl>
<div class="links">{links}</div>
</div></article>"""
others = []
for r in csv.DictReader(open(os.path.join(repo, "config", "species_catalogue.tsv")), delimiter="\t"):
    if r["Nuclear genome assembly"].lower().startswith("none"):
        def cl(x):
            x = re.sub(r"^x(\d+) ", lambda m: m.group(1) + " libraries, ", x.strip())
            return "–" if x.lower() in ("none", "", "none found") else e(x)
        dna = r["Genomic reads"]
        if dna.strip().lower() in ("none", "") and "metagenome" in r["Other omics / associated resources"].lower(): dna = r["Other omics / associated resources"]
        others.append(f'<tr><td>{it(r["Species"])}</td><td>{cl(dna)}</td><td>{cl(r["RNA-seq and transcriptome assemblies"])}</td><td>{cl(r["Sequencing-project status"])}</td></tr>')
n_sp = len(order)
# ---- what each sequenced animal carried, and the gene-transfer result, across species
def pct_small(x):
    v = 100 * x
    return f"{v:.1f}%" if v >= 1 else (f"{v:.2f}%" if v >= 0.1 else f"{v:.3f}%")
carried_rows = ""
all_transfer = all(os.path.exists(os.path.join(SITE, s_["slug"], "gene-transfer.html")) for s_ in order)
for s_ in order:
    slug = s_["slug"]
    if slug == "elysia-crispata" and ck:
        cp = ck["classes"]["cp"]["bp"] / ck["library"]["bp"]
        dp = sum(a["bp"] / pl for a, pl in zip(ck["assignment"][:2], (84040, 89458))); dn = ck["library"]["bp"] / 892.65e6
        carried_rows += (f'<tr><td><a href="{slug}/kleptoplasts.html">{it(s_["species"])}</a></td><td>{it("Penicillus capitatus")} and {it("P. lamourouxii")}</td><td>the animal assembled</td>'
                         f'<td class="num">{pct_small(cp)}</td><td class="num">{dp / dn:,.0f}</td></tr>')
        continue
    K = KLEPTO.get(slug); tx = STEXT.get(slug, {})
    if not K or not os.path.exists(os.path.join(SITE, slug, "kleptoplasts.html")): continue
    cells_f = []; cells_d = []
    span = CAT[slug]["assembly"]["span"]; pl = tx.get("plastome_length")
    for L in K["libraries"]:
        if L["plastid"]["bp"] < 5000: continue
        cells_f.append(pct_small(L["plastid"]["fraction"]))
        if pl:
            r = (L["plastid"]["bp"] / pl) / (L["bases"] / span); cells_d.append(f"{r:,.0f}" if r >= 10 else f"{r:.1f}")
    who = {"elysia-leucolegnote": "two other animals", "elysia-chlorotica": "two adults", "elysia-timida": "the animals assembled"}.get(slug, "the animal assembled")
    floor = " or more" if slug in ("elysia-timida", "elysia-leucolegnote") else ""
    carried_rows += (f'<tr><td><a href="{slug}/kleptoplasts.html">{it(s_["species"])}</a></td><td>{re.sub(r" [(]in two other animals[)]", "", tx.get("card", ""))}</td><td>{who}</td>'
                     f'<td class="num">{" and ".join(cells_f)}{floor}</td><td class="num">{" and ".join(cells_d)}</td></tr>')
carried_html = "" if not carried_rows else f"""
<section class="block"><div class="wrap">
<h2>What the sequenced animals were carrying</h2>
<p>The plastid DNA among the reads of each animal names the alga it had been eating, and shows how much chloroplast DNA it held when it was sequenced.</p>
<div class="tablewrap"><table class="data"><thead><tr><th>Species</th><th>Chloroplasts from</th><th>Reads from</th><th class="num">Plastid share of sequenced DNA</th><th class="num">Plastid genomes per slug genome</th></tr></thead>
<tbody>{carried_rows}</tbody></table></div>
<p class="small muted" style="margin-top:10px;max-width:92ch">The last column divides how deeply the plastid genome is covered by how deeply the slug genome is covered. It depends on which tissue was sequenced and on when the animal last fed. For <i>E. timida</i> only 49 kb of the alga’s plastid genome is known, so the figures are for that part. For <i>E. leucolegnote</i> they count the plastid DNA that resembles a reference, which leaves out much of what the other algae contributed.</p>
{'<p style="margin-top:18px;max-width:80ch"><strong>No algal DNA was found inserted in any slug genome.</strong> In every species, reads that join chloroplast DNA to slug DNA are rare and scattered across the genome, and they turn up at rates similar to joins that cannot be real, such as the slug’s mitochondrial DNA joined to its nuclear DNA. Library preparation makes such joins. Where an assembly shows chloroplast DNA inside or at the end of a slug sequence (<i>E. crispata</i>, <i>E. chlorotica</i>, <i>E. leucolegnote</i>), reads do not bear the join out. The assemblies also hold chloroplast DNA as sequences of its own (<i>E. chlorotica</i>, <i>E. marginata</i>), a piece of an algal mitochondrial genome (<i>E. papillosa</i>) and bacterial DNA, as much as 2 Mb of it (<i>E. chlorotica</i>, <i>E. marginata</i>, <i>E. leucolegnote</i>). Each species has a gene-transfer page with the evidence and its limits.</p>' if all_transfer else ''}
</div></section>
"""
home = f"""
<section class="home-head"><div class="wrap">
<h1>The circadian clock in {word(n_sp)} <i>Elysia</i> genomes</h1>
<p class="lede">Sea slugs of the genus <i>Elysia</i> keep chloroplasts taken from the algae they eat, some for months. This site catalogues the clock genes in every <i>Elysia</i> genome that has been sequenced, names the algae whose chloroplasts each sequenced animal carried, and tests whether any algal DNA has entered the slug genome.</p>
</div></section>

<section class="block"><div class="wrap">
<h2>The species</h2>
<p>{word(n_sp).capitalize()} species have a nuclear genome assembly. Each card opens that species’ pages.</p>
<div class="sp-grid">
{"".join(card(s) for s in order)}
</div>
</div></section>

{carried_html}
<section class="block"><div class="wrap">
<div class="cols">
<div><h3>One clock, {word(n_sp)} genomes</h3>
<p>{len(all_intact)} of the 47 clock genes are whole in the reference assembly of every species. {HOME_CORE} The <a href="compare.html">comparison page</a> has the gene-by-species table, the protein alignments and a tree of the species.</p></div>
<div><h3>What a gap means</h3>
<p>The assemblies differ in quality far more than the animals differ in genes. A gene shown as partial or missing is, in most cases, a hole in that assembly. Each species page says what its assembly can and cannot show.</p></div>
<div><h3>What is not here</h3>
<p>Nothing on this site measures a rhythm. Gene roles are inferred from flies and mammals, and each genome comes from a single animal.</p></div>
</div>
</div></section>

<section class="block"><div class="wrap">
<h2>Other <i>Elysia</i> with sequence data but no genome</h2>
<p>From a survey of the public archives on 13 November 2026. These species have RNA-seq or shotgun reads only, or a genome project that has not released data, so they have no pages here.</p>
<div class="tablewrap"><table class="data"><thead><tr><th>Species</th><th>Genomic reads</th><th>RNA-seq</th><th>Genome project</th></tr></thead>
<tbody>{"".join(others)}</tbody></table></div>
</div></section>
"""
open(os.path.join(SITE, "index.html"), "w").write(chrome.page("Elysia clocks: circadian clock genes in the sea slug genus Elysia", home, root="", site_current="index.html", today=today,
                                                           scripts='<script src="assets/pages.js"></script>'))

# ------------------------------------------------------------------ methods page (site-wide)
rows = "".join(f'<tr><td>{it(s["species"])}</td><td class="mono">{e(s["accession"])}</td><td>{e((CAT[s["slug"]]["assembly"]["name"] if s["slug"] in CAT else "xgElyCris1.1"))}</td><td>{LEVEL_NOTE[s["level"]]}</td>'
               f'<td>{ {"genbank": "GenBank", "ensembl": "Ensembl rapid release", "senckenberg": "authors’ site, not in GenBank", "none": "none public"}[s["annotation"]] }</td>'
               f'<td>{("annotated genes and aligned proteins" if (s["slug"] in CAT and CAT[s["slug"]]["annotated"]) else ("curated catalogue (four assemblies)" if s["slug"] == "elysia-crispata" else "aligned proteins"))}</td></tr>' for s in order)
photo_rows = "".join(f'<tr><td>{it(c["species"])}</td><td><a href="{c["source"]}">{e(c["photographer"])}</a></td><td><a href="{c["licence_url"]}">{e(c["licence"])}</a></td><td>{e(c["place"])}, {e(c["date"])}</td></tr>' for c in credits.values())
shutil.copy(os.path.join(RG, "clock_matrix.tsv"), os.path.join(SITE, "downloads", "clock_matrix.tsv"))
if GN.get("tree"): open(os.path.join(SITE, "downloads", "species_tree.nwk"), "w").write(GN["tree"] + "\n")
dl = [("clock_matrix.tsv", "Gene by species table (TSV)", "state, completeness, identity to E. crispata and location of every gene in every species"),
      ("species_tree.nwk", "Species tree (Newick)", "from the clock proteins complete in every species")]
for s in order:
    if s["slug"] in CAT:
        dl.append((f"{s['slug']}_clock_proteins.faa", f"{s['species']} clock proteins (FASTA)", "")); dl.append((f"{s['slug']}_clock_cds.fna", f"{s['species']} coding sequences (FASTA)", ""))
read_rows = ""
rp = os.path.join(repo, "config", "species_reads.tsv")
if os.path.exists(rp):
    for r in csv.DictReader(open(rp), delimiter="\t"):
        K = KLEPTO.get(r["slug"]) or {}
        L = next((x for x in K.get("libraries", []) if x["library"] == r["library"]), None)
        if not L: continue
        read_rows += (f'<tr><td>{it(BYSLUG[r["slug"]]["species"])}</td><td>{e(r["animal"][0].upper() + r["animal"][1:])}</td><td>{e(r["platform"])}</td><td class="mono">{e(r["runs"].replace(",", ", "))}</td>'
                      f'<td class="num">{L["bases"] / 1e9:.1f} Gb</td></tr>')
for s_ in order:
    for suffix, label in (("_kleptoplast.json", "chloroplast results (JSON)"), ("_gene_transfer.json", "gene-transfer results (JSON)")):
        if os.path.exists(os.path.join(SITE, "downloads", s_["slug"] + suffix)): dl.append((s_["slug"] + suffix, f"{s_['species']} {label}", ""))
    for f in sorted(os.listdir(os.path.join(SITE, "downloads"))):       # plastid contigs and barcode sequences (45_export_tables.py)
        m = re.match(re.escape(s_["slug"]) + r"(?:_(\w+?))?_(plastid_contigs|barcodes)\.fa$", f)
        if m: dl.append((f, f"{s_['species']} {'plastid contigs' if m.group(2) == 'plastid_contigs' else 'rbcL and tufA sequences'} (FASTA)",
                         {"adult": "second adult, short reads", "hainan": "Hainan animal", None: ("Shenzhen animal" if s_["slug"] == "elysia-leucolegnote" else "")}.get(m.group(1), m.group(1) or "")))
reads_methods = "" if not read_rows else f"""
<h2 id="reads">Whose chloroplasts: reading the reads</h2>
<p>For <i>E. crispata</i> all 47 Gb of reads were examined (see <a href="elysia-crispata/methods.html">its methods page</a>). For the other species the read sets below were used.</p>
<div class="tablewrap"><table class="data"><thead><tr><th>Species</th><th>DNA from</th><th>Sequencing</th><th>Runs</th><th class="num">Examined</th></tr></thead><tbody>{read_rows}</tbody></table></div>
<details class="method" open><summary>1. Reference panel</summary>
<p>1,271 algal plastid sequences (241 whole plastid genomes of green algae; plastid records of Dasycladales, Xanthophyceae and Cladophorales; a land plant, a diatom, a red alga and a cyanobacterium to catch anything else), 117 algal mitochondrial genomes, 7,352 nuclear ribosomal DNA sequences of algae, six slug ribosomal DNA sequences and the mitochondrial genomes of nine <i>Elysia</i> species. Low-complexity sequence in the panel is masked for the first pass.</p></details>
<details class="method"><summary>2. First pass: every read against the panel</summary>
<p>Each run was fetched from the SRA, converted and mapped to the panel as one stream (minimap2; k 15, w 10 for long reads, the short-read preset for read pairs). Reads with a match were kept, both mates of a pair if either matched. Two kinds of read were thinned to a sample: reads whose only match is a short local one to a panel window that collects such matches from very many reads (slug repeats that resemble a stretch of a panel sequence), and reads that match nothing but ribosomal DNA at slug-like identity.</p></details>
<details class="method"><summary>3. Second pass: base by base</summary>
<p>The kept reads were aligned to the species’ own assembly together with the panel. Each read is described by its best set of non-overlapping alignments. A stretch counts as plastid, mitochondrial or slug when it matches at 97% or better (82% for the older PacBio RS reads of <i>E. chlorotica</i>). Sequences of an assembly that are organelle DNA from end to end are counted as organelle, and plastid-like stretches inside other sequences are hidden for this step, so that a plastid read is not explained away by a plastid fragment left in the assembly.</p></details>
<details class="method"><summary>4. Naming the alga</summary>
<p>Plastid reads are assigned to the reference they match best. Two barcode genes, rbcL and tufA, are then compared with 11,600 GenBank records of Ulvophyceae and Xanthophyceae. With accurate long reads, every read that spans a barcode contributes its own sequence, identical sequences are counted, and each distinct sequence is compared. With short reads, the plastid pairs are first assembled into contigs (MEGAHIT) and the barcodes are cut from the contigs.</p></details>

<h2 id="transfer">Gene transfer: where algal DNA was looked for</h2>
<details class="method" open><summary>1. In the assembly, as DNA</summary>
<p>Each assembly was cut into 10 kb windows (1 kb overlap) and compared with the panel. Sequences covered 50% or more by organelle matches are reported as organelle DNA. Elsewhere, every match of 200 bases or more at 80% identity or better is listed, after removing matches to low-complexity sequence, matches to panel windows hit from 15 or more places in the assembly, and the slug’s own ribosomal DNA.</p></details>
<details class="method"><summary>2. In the assembly, as protein</summary>
<p>471 proteins were aligned to each assembly with miniprot, keeping weak and partial matches: the plastid-encoded proteins of <i>Udotea flabellum</i>, photosynthesis proteins of <i>Chlamydomonas</i> and <i>Arabidopsis</i> (nuclear and plastid), UniProt proteins of siphonous green algae and all proteins of <i>Vaucheria litorea</i>. Each matching place was translated and compared with animal, algal, cyanobacterial, bacterial and fungal proteins (DIAMOND). A place is the slug’s own gene when the best match among animals other than sea slugs is at least as good as the best algal or microbial match, or within 10% of it while a protein of another sea slug beats both (very conserved proteins such as GAPDH). The gene set published with an assembly is never counted as animal evidence for that assembly: the published gene sets of <i>E. chlorotica</i> and <i>E. marginata</i> include genes on bacterial and plastid sequences, and such a gene would otherwise vouch for itself.</p>
<p>Every sequence on which a place of another kind was found was then translated as a whole in all six frames (DIAMOND blastx with frameshifts and range culling) and each protein-coding stretch classified in the same way. A sequence is called bacterial when three or more stretches, and at least half of all, are of bacterial type. The nearest bacterial genera come from the Swiss-Prot entries matched best; they name relatives, not species.</p></details>
<details class="method"><summary>3. In the reads</summary>
<p>A long read is chimeric when two neighbouring stretches of 300 bases or more are of different kinds, for example plastid and slug. A read pair is mixed when its mates are of different kinds. The rate for plastid DNA is compared with the rate for the slug’s own mitochondrial DNA, which measures how often the library joins unrelated molecules. A real insertion would be seen in many independent molecules at one place; joins are grouped by position and molecules sharing both ends are counted once.</p></details>
<details class="method"><summary>4. Testing a place in an assembly with reads</summary>
<p>Where an assembly holds an organelle-like stretch inside another sequence, reads were aligned to that place with 3 kb either side. At each edge the reads that run through are counted against the reads that stop there and continue into something else. For <i>E. chlorotica</i> the reads tested were those that matched the panel, to which the places had been added. For the one place in <i>E. leucolegnote</i>, whose plastid part is too unlike the panel for that, every read of both libraries was compared with the place itself.</p></details>
"""
methods = f"""
<section class="page-head"><div class="wrap">
<h1>Methods and data</h1>
<p class="lede">How the genes were found in each genome and compared, how the reads and assemblies were searched for algal DNA, and where the data came from. The <i>E. crispata</i> work is described in more detail on <a href="elysia-crispata/methods.html">its own methods page</a>.</p>
</div></section>
<section class="block article"><div class="wrap">
<h2>Genomes used</h2>
<div class="tablewrap"><table class="data"><thead><tr><th>Species</th><th>Assembly</th><th>Name</th><th>Level</th><th>Gene annotation</th><th>Clock gene models from</th></tr></thead><tbody>{rows}</tbody></table></div>
<p class="small muted" style="margin-top:10px">Species and data were taken from a survey of the public archives made on 13 November 2026.</p>

<h2>Finding the clock genes</h2>
<details class="method" open><summary>1. Starting point</summary>
<p>The 47 genes of the <i>E. crispata</i> catalogue, each identified there by reciprocal best hits and a family tree against fly, mammal and mollusc proteins.</p></details>
<details class="method"><summary>2. Candidates in each genome</summary>
<p>The 47 proteins were aligned to every genome with miniprot, allowing introns up to 500 kb and keeping weaker second alignments. Where GenBank has a gene annotation (<i>E. chlorotica</i>, <i>E. marginata</i>), the annotated proteins most similar to each catalogue protein were added as candidates.</p></details>
<details class="method"><summary>3. Telling orthologues from paralogues</summary>
<p>Each candidate was searched (DIAMOND) against every <i>E. crispata</i> protein of both annotations plus the 47 catalogue proteins. It was kept for a gene only if that gene was its best match among the catalogue and no other <i>E. crispata</i> protein fitted clearly better. The model for a gene is the kept candidate with the most residues identical to the <i>E. crispata</i> protein; an annotated gene is preferred when it is within 2% of the best alignment. Pieces of one protein that lie in order on one strand are joined as a gene with a long intron; pieces on different sequences count as a split gene.</p></details>
<details class="method"><summary>4. Completeness</summary>
<p>The models of all species were aligned (MAFFT) with the orthologues of three other gastropods. Completeness is the share of alignment positions present in at least half of the sequences that a model fills. Intact means 90% or more with no frameshift or stop in the genome alignment. This yardstick does not treat the <i>E. crispata</i> model as a gold standard, but a model aligned from <i>E. crispata</i> still cannot contain what the <i>E. crispata</i> model lacks.</p></details>
<details class="method"><summary>5. Comparison and tree</summary>
<p>Identity between two proteins is counted over the positions where both have a residue. The tree uses the proteins complete in every species, positions with at most 20% gaps, IQ-TREE with LG+F+G4 and 1,000 ultrafast bootstraps.</p></details>
<details class="method"><summary>6. Domains</summary><p>Pfam-A, hmmsearch with gathering thresholds.</p></details>
<p class="small muted">Scripts: <code>scripts/genus/</code> in the repository, numbered in the order they run.</p>
{reads_methods}

<h2>Downloads</h2>
<ul class="downloads">{"".join(f'<li><a href="downloads/{n}" download>{e(l)}</a><span>{e(t)}</span></li>' for n, l, t in dl)}</ul>

<h2 id="photos">Photographs</h2>
<p>The photographs are not of the sequenced animals. They come from iNaturalist observations and are used under the licence shown. Two species had no photo under a fully open licence, so theirs carry a non-commercial licence.</p>
<div class="tablewrap"><table class="data"><thead><tr><th>Species</th><th>Photographer</th><th>Licence</th><th>Where and when</th></tr></thead><tbody>{photo_rows}</tbody></table></div>
</div></section>
"""
open(os.path.join(SITE, "methods.html"), "w").write(chrome.page("Methods and data. Elysia clocks", methods, root="", site_current="methods.html", today=today, scripts='<script src="assets/pages.js"></script>'))

# ------------------------------------------------------------------ old addresses
for old in ("kleptoplasts.html", "gene-transfer.html"):
    open(os.path.join(SITE, old), "w").write(chrome.redirect("elysia-crispata/" + old))
old_data = os.path.join(SITE, "data", "clock.js")
if os.path.exists(old_data): os.remove(old_data)
print("home, compare, methods and", len(CAT), "species pages written;", len(all_intact), "genes whole in every genome; linkage pairs", pairs_kept, "of", pairs_total)
