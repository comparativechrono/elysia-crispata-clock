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
    lab = f'{c["ok"]} intact, {c["warn"]} complete with a frameshift or stop codon, {c["part"]} partial or fragmented, {c["miss"]} not found'
    return f'<div class="statusbar" role="img" aria-label="{lab}" title="{lab}">' + "".join(f'<i class="{k}" style="flex:{c[k]}"></i>' for k in ("ok", "warn", "part", "miss") if c[k]) + "</div>"
def tabs_for(slug):
    t = [("index.html", "Clock genes")]
    # the pages of the other species are written by 21_build_species_pages.py, which runs first
    if slug == "elysia-crispata" or os.path.exists(os.path.join(SITE, slug, "kleptoplasts.html")): t.append(("kleptoplasts.html", "Kleptoplasts"))
    if slug == "elysia-crispata" or os.path.exists(os.path.join(SITE, slug, "gene-transfer.html")): t.append(("gene-transfer.html", "Gene transfer"))
    if slug == "elysia-crispata": t.append(("methods.html", "Methods"))
    return t

# ------------------------------------------------------------------ per-species clock pages
LEVEL_NOTE = {"chromosome": "chromosome-level", "scaffold": "scaffold-level", "contig": "contig-level"}
CAVEAT = {
 "elysia-chlorotica": "The assembly consists of 9,989 scaffolds built from short reads and early PacBio reads. Genes that appear partial here, or that have an interrupted reading frame, are more likely to reflect errors in the assembly than real changes in the genes.",
 "elysia-marginata": "The assembly was built from short reads only and consists of 14,149 scaffolds with 82 Mb of gaps. Some genes are therefore split between scaffolds or truncated at a gap.",
 "elysia-timida": "This is a chromosome-level assembly. Because all gene models were derived from the <i>E. crispata</i> proteins, regions that are poorly conserved or that are absent from the <i>E. crispata</i> model may be missing from the models shown here.",
 "elysia-leucolegnote": "In the public version of this assembly 35% of the bases are N; the sequence appears to have been hard-masked for repeats. An exon that lies in a masked region cannot be detected, so genes that appear partial here may be complete in the genome.",
 "elysia-atroviridis": "This is a chromosome-level assembly built from PacBio HiFi reads. Because all gene models were derived from the <i>E. crispata</i> proteins, regions that are poorly conserved or that are absent from the <i>E. crispata</i> model may be missing from the models shown here.",
 "elysia-papillosa": "This is a contig-level assembly built from PacBio HiFi reads, which has not yet been scaffolded into chromosomes. Because all gene models were derived from the <i>E. crispata</i> proteins, regions that are poorly conserved or that are absent from the <i>E. crispata</i> model may be missing from the models shown here.",
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
            ident_txt = "No locus in the genome passed the reciprocal best-hit test for this gene."
        elif g["model_kind"] == "annotation":
            ident_txt = (f"The gene model is the annotated protein of this species (GenBank {g['protein_id']}). It passed the reciprocal best-hit test: its best match among all E. crispata proteins is {lab}. "
                         f"The protein is {round(100 * ident) if ident else '?'}% identical to the E. crispata protein and covers {round(100 * g['consensus_coverage'])}% of the alignment positions that are present in most species.")
        else:
            ident_txt = (f"The gene was identified by aligning the E. crispata {lab} protein to the genome. The locus passed the reciprocal best-hit test: the best match of its translation among all E. crispata proteins is {lab}. "
                         f"The protein is {round(100 * ident) if ident else '?'}% identical to the E. crispata protein and covers {round(100 * g['consensus_coverage'])}% of the alignment positions that are present in most species."
                         + (" The annotation of this species also contains a gene at this locus, but the aligned model covers more of the protein." if g.get("annotation_protein") else ""))
        if g.get("extra_loci"):
            nx = len(g["extra_loci"])
            ident_txt += (f" {word(nx).capitalize()} further {'loci' if nx > 1 else 'locus'} in the assembly also passed the test ({'; '.join(x['chrom'] + ':' + num(x['start']) + '–' + num(x['end']) for x in g['extra_loci'][:3])}). "
                          f"{'These may be additional copies' if nx > 1 else 'This may be a second copy'} of the gene, or the other haplotype retained in the assembly.")
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

word = chrome.word
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
    how = (f"Of the {n_ann + n_proj} gene models, {n_ann} are annotated genes of this species. The other {n_proj} were obtained by aligning the <i>E. crispata</i> protein to the genome, at loci where the annotation has no gene or only a shorter one."
           if n_ann else f"No gene annotation is available in GenBank for this assembly, so all {n_proj} gene models were obtained by aligning the <i>E. crispata</i> proteins to the genome.")
    specimen = "specimens" if slug == "elysia-timida" else "specimen"
    body = f"""
<section class="intro"><div class="wrap">
<h1 style="max-width:22ch">Circadian clock genes of {it(sp["species"]).replace(" ", "&nbsp;")}</h1>
<p class="lede">The 47 clock genes catalogued in <i>E. crispata</i> were searched for in {e(A["name"])}, the {LEVEL_NOTE[sp["level"]]} assembly of {it(sp["species"])} ({A["span"] / 1e6:,.0f} Mb). For each gene, this page gives its position in the assembly, whether it is complete, and the protein it encodes. The sequenced {specimen} came from {e(sp["origin"])}.</p>
<div class="clockwork">
<div id="schematic"></div>
<div class="clockwork-list" id="schematic-list"></div>
<div class="clockwork-legend" aria-label="Key to gene status in {e(A["name"])}">
<span><i class="glyph ok"></i>intact in {e(A["name"])}</span>
<span><i class="glyph warn"></i>complete, but with a frameshift or stop codon in the assembly sequence</span>
<span><i class="glyph part"></i>partial or fragmented</span>
<span><i class="glyph miss"></i>not found</span>
</div>
</div>
<p class="clockwork-caveat">Select a gene to see its details. The diagram shows the interactions established in <i>Drosophila</i> and mouse; the functions of these genes have not been tested in <i>Elysia</i>.</p>
</div></section>

<section class="block"><div class="wrap">
<div class="cols">
<div><h3>Genes found</h3>
<p>Of the {n} genes, {len(ok)} are intact{f", and {word(len(disc))} more {'is' if len(disc) == 1 else 'are'} complete but {'contains' if len(disc) == 1 else 'contain'} a frameshift or stop codon in the assembly sequence ({names(disc)})" if disc else ""}. {f"{word(len(part)).capitalize()} {'is' if len(part) == 1 else 'are'} partial or fragmented ({names(part)})" if part else "No gene is partial or fragmented"}{f", and {word(len(miss))} {'was' if len(miss) == 1 else 'were'} not found ({names(miss)})" if miss else ""}.</p></div>
<div><h3>Identification</h3>
<p>{how} Each gene model was accepted only if it passed a reciprocal best-hit test against all <i>E. crispata</i> proteins, which distinguishes orthologues from paralogues. The proteins are {round(100 * min(ids))} to {round(100 * max(ids))}% identical to their <i>E. crispata</i> orthologues (median {round(100 * sorted(ids)[len(ids) // 2])}%).</p></div>
<div><h3>Limitations of the assembly</h3>
<p>{CAVEAT.get(slug, "")}</p></div>
</div>
</div></section>
""" + (f"""
<section class="block"><div class="wrap">
<h2>Chromosomal locations</h2>
<p>The {len(A["chromosomes"])} chromosomes of {e(A["name"])} are drawn to scale and numbered as in the assembly; chromosome numbers do not correspond between species. Each mark is a clock gene; select a mark to see the details of that gene.</p>
<div class="karyo-card"><div class="karyo" id="karyo"></div><div class="offmap" id="offmap"></div></div>
</div></section>
""" if chrom else "") + f"""
<section class="block"><div class="wrap">
<h2>All {n} genes</h2>
<div class="filters" id="filters"></div>
<div class="tablewrap"><table class="genes">
<thead><tr><th>Gene</th><th>Category</th><th class="num">Protein length (aa)</th><th>Status in {e(A["name"])}</th><th class="num" title="Percentage of identical residues over the aligned part of the two proteins">Identity to <i>E. crispata</i></th><th>Gene model</th></tr></thead>
<tbody id="genes-body"></tbody></table></div>
<p class="small muted" style="margin-top:12px;max-width:90ch">A gene model is classed as intact when it covers at least 90% of the alignment positions that are present in most species and its coding sequence in the assembly contains no frameshift or stop codon. Models obtained by aligning the <i>E. crispata</i> protein may lack poorly conserved regions at the ends of a protein.</p>

<h3 style="margin-top:34px">Genes searched for but not found</h3>
<p style="max-width:80ch">The seven clock-associated genes that were not found in <i>E. crispata</i> were also searched for in this genome, with the same <i>Drosophila</i> and mammalian query proteins. {"None of them was found." if not nf_found else "The following were found: " + ", ".join(x["symbol"] for x in nf_found) + "."}</p>
</div></section>
<dialog class="gene" id="gene" aria-label="Gene details"></dialog>
"""
    os.makedirs(os.path.join(SITE, slug), exist_ok=True)
    src = f"Data for {it(sp['species'])}: genome assembly {e(A['name'])} ({e(A['accession'])}). The data sources and publications for all species are listed on the <a href=\"../methods.html\">methods page</a>."
    open(os.path.join(SITE, slug, "index.html"), "w").write(chrome.page(
        f"{sp['species']} clock genes", body, root="../", species=sp["species"], tabs=tabs_for(slug), tab_current="index.html", today=today, sources=src,
        description=f"Circadian clock genes of the sea slug {sp['species']}, identified in the genome assembly {A['name']}.",
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
CORE_TXT = ("In each of the other species the two core activators (Clock and Bmal1), the repressors Period, Timeless and a mammalian-type Cryptochrome, and a light-sensitive <i>Drosophila</i>-type Cryptochrome are all present as complete genes."
            if core_whole else "The genes of the core loop (Clock, Bmal1, Period, Timeless and two Cryptochromes) are present in every species, although not every gene is complete in every assembly.")
bad = sorted(((sum(1 for sym in symbols if STATUS_CLASS[status_of(s["slug"], sym)] != "ok"), s) for s in order if s["slug"] != "elysia-crispata"), key=lambda x: -x[0])
WORST_TXT = "" if {bad[0][1]["slug"], bad[1][1]["slug"]} != {"elysia-chlorotica", "elysia-marginata"} else (f"The two assemblies built entirely or mainly from short reads have the largest numbers of genes that are not intact ({bad[0][0]} in {it(SHORT[bad[0][1]['slug']])} and {bad[1][0]} in {it(SHORT[bad[1][1]['slug']])}), "
             f"whereas the four assemblies built from high-accuracy long reads have {min(b[0] for b in bad[2:])} to {max(b[0] for b in bad[2:])} each.")
# the core loop in E. crispata: whole in the reference, or only in another assembly of the species
cr_ref_gap = [sym for sym in CORE if STATUS_CLASS[status_of("elysia-crispata", sym)] not in ("ok", "warn")]
cr_any = all(any(v.get("status") in ("complete", "complete_discrepancy") for v in CG[sym]["assemblies"].values()) for sym in CORE)
HOME_CORE = ""
if core_whole and cr_any:
    HOME_CORE = "The genes of the core loop (Clock, Bmal1, Period, Timeless and two Cryptochromes) are complete in every species"
    HOME_CORE += (f" (in <i>E. crispata</i>, {' and '.join(LABEL.get(x, x) for x in cr_ref_gap)} {'is' if len(cr_ref_gap) == 1 else 'are'} complete in another assembly of the species but not in the reference assembly)." if cr_ref_gap else ".")
def idlist(gs):
    """'A (95%), B (94%) and C (93%)': mean identity of each protein among the Elysia species"""
    xs = [f'{g["label"]} ({round(100 * g["mean_identity"])}%)' for g in gs]
    return xs[0] if len(xs) == 1 else ", ".join(xs[:-1]) + " and " + xs[-1]
compare_body = f"""
<section class="page-head"><div class="wrap">
<h1>Comparison of clock genes across {word(len(order))} <i>Elysia</i> species</h1>
<p class="lede">The 47 clock genes catalogued in <i>E. crispata</i> were searched for in every <i>Elysia</i> species with a genome assembly, and the encoded proteins were aligned and compared.</p>
<div class="answer">
<p><strong>The set of clock genes is conserved across the genus.</strong> Of the 47 genes, {len(all_intact)} are complete in the reference assembly of all {word(len(order))} species. {CORE_TXT}</p>
<p><strong>Missing and partial genes mostly reflect the quality of the assemblies.</strong> {WORST_TXT} The exceptions are the Go-opsins, which vary in number between species, and r-opsin-like 1, for which the gene model is incomplete in <i>E. crispata</i> itself.</p>
<p><strong>The proteins differ in their degree of conservation.</strong> The most conserved across the genus, by mean pairwise identity, are {idlist(slow[:3])}; the least conserved are {idlist(fast[:3])}. The values for the core loop are {idlist(sorted(core, key=lambda g: -g["mean_identity"]))}.</p>
</div>
</div></section>

<section class="block"><div class="wrap">
<h2>Gene status by species</h2>
<p>Each cell shows the status of the gene in the reference assembly of that species and the identity of its protein to the <i>E. crispata</i> protein. Select a gene name to see the alignment, or a cell to see the gene in that species.</p>
<div class="legend"><span><i class="glyph ok"></i>intact</span><span><i class="glyph warn"></i>complete, with a frameshift or stop codon</span><span><i class="glyph part"></i>partial or fragmented</span><span><i class="glyph miss"></i>not found</span></div>
<div class="tablewrap"><table class="matrix" id="matrix"></table></div>
<p class="small muted" style="margin-top:12px;max-width:92ch">For <i>E. crispata</i> the status refers to the Sanger Institute assembly xgElyCris1.1; all 47 genes are complete in at least one other assembly of that species. A gene model is classed as intact when it covers at least 90% of the alignment positions that are present in most sequences and contains no frameshift or stop codon.</p>
</div></section>

<section class="block"><div class="wrap">
<h2>Protein conservation</h2>
<p>Mean pairwise identity among the <i>Elysia</i> species, for each protein that is complete in at least two of them. Select a bar to see the alignment.</p>
<div class="fig"><div id="cons"></div><p class="figcap">Identity is calculated over the positions at which both proteins have a residue. For the opsins and photolyases, low values partly reflect the difficulty of distinguishing paralogues.</p></div>
</div></section>

<section class="block"><div class="wrap">
<h2>Protein alignments</h2>
<div class="gene-pick"><label for="pick">Protein</label><select id="pick"></select><span class="muted small" id="pick-note"></span></div>
<div id="idmatrix"></div>
<div class="aln-wrap" style="margin-top:14px"><div class="aln" id="aln"></div></div>
<p class="small muted" style="margin-top:10px;max-width:92ch">Shaded residues differ from the most common residue at that position among the <i>Elysia</i> sequences. The three sequences below the dotted line are from other gastropods and are included as outgroups. The alignments are available in <code>results/genus/alignments/</code> in the repository.</p>
</div></section>

<section class="block"><div class="wrap">
<h2>Species tree</h2>
<p>The tree was inferred from a concatenated alignment of the {len(GN["tree_genes"])} clock proteins that are complete in every species ({num(GN["tree_columns"])} aligned positions) and is rooted on the sea hare <i>Aplysia californica</i>.</p>
<div class="fig"><div id="tree"></div><p class="figcap">Maximum-likelihood tree (IQ-TREE, LG+F+G4). Numbers are ultrafast bootstrap support values, and branch lengths are in substitutions per site. A tree based on clock proteins alone is intended as a guide to the relationships among the species on this site, not as a phylogeny of the genus.</p></div>
</div></section>
""" + (f"""
<section class="block"><div class="wrap">
<h2>Conservation of chromosomal linkage</h2>
<p>{word(len(chromsp)).capitalize()} of the assemblies are chromosome-level. Each has 15 chromosomes, which are numbered by size, so chromosome numbers do not correspond between species. Linkage can nevertheless be compared: of the {pairs_total} pairs of clock genes that lie on the same chromosome in <i>E. crispata</i>, {pairs_kept} also lie on the same chromosome in {"all " + word(len(chromsp) - 1) + " other species" if len(chromsp) > 2 else "the other species"}.</p>
<div class="tablewrap"><table class="data" id="linkage"></table></div>
</div></section>
""" if len(chromsp) > 1 else "") + """
<script>window.COMPARE_ROOT = "";</script>
"""
open(os.path.join(SITE, "compare.html"), "w").write(chrome.page(
    "Comparison of clock genes across seven Elysia species", compare_body, root="", site_current="compare.html", today=today,
    scripts='<script src="data/genus.js"></script>\n<script src="assets/pages.js"></script>\n<script src="assets/compare.js"></script>'))

# ------------------------------------------------------------------ home page
def card(s):
    slug = s["slug"]; c = counts(slug); cr = credits.get(slug)
    A = CAT[slug]["assembly"] if slug in CAT else dict(name="xgElyCris1.1", accession="GCA_963854125.1", span=542.6e6)
    photo = (f'<a class="ph" href="{slug}/index.html"><img src="assets/photos/{cr["file"]}" alt="{e(s["species"])}" width="{cr["width"]}" height="{cr["height"]}" loading="lazy">'
             f'</a>' if cr else '<div class="ph"></div>')
    credit = f'<span class="credit">Photo: <a href="{cr["source"]}">{e(cr["photographer"])}</a>, {e(cr["licence"])}</span>' if cr else ""
    ann = {"genbank": "GenBank", "ensembl": "Ensembl", "senckenberg": "not in GenBank", "none": "not publicly available"}[s["annotation"]]
    k = KLEPTO.get(slug); kl = ""
    if slug == "elysia-crispata" and ck:
        top = ck["assignment"][:2]; kl = ", ".join(f'{it(a["taxon"])} {round(100 * a["fraction"])}%' for a in top)
    elif k and STEXT.get(slug, {}).get("card"): kl = STEXT[slug]["card"]
    links = "".join(f'<a href="{slug}/{h}">{t}</a>' for h, t in tabs_for(slug) if h != "methods.html")
    return f"""<article class="sp-card">
<div style="position:relative">{photo}{credit}</div>
<div class="tx">
<h3><a href="{slug}/index.html">{e(s["species"])}</a></h3>
<p class="where">Sequenced specimen: {e(s["origin"])}</p>
<div><div class="small" style="margin-bottom:3px"><strong>{c["ok"]}</strong> of 47 clock genes intact in the assembly</div>{statusbar(slug)}</div>
<dl>
<dt>Assembly</dt><dd>{LEVEL_NOTE[s["level"]]}, {A["span"] / 1e6:,.0f} Mb <span class="muted small">{e(A["accession"])}</span></dd>
<dt>Gene annotation</dt><dd>{ann}</dd>
<dt>Long reads</dt><dd>{e(s["long_reads"])}</dd>
{f"<dt>Kleptoplasts from</dt><dd>{kl}</dd>" if kl else ""}
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
        carried_rows += (f'<tr><td><a href="{slug}/kleptoplasts.html">{it(s_["species"])}</a></td><td>{it("Penicillus capitatus")} and {it("P. lamourouxii")}</td><td>assembled individual</td>'
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
            # depth of the plastid genome against depth of the slug genome: the bases that are neither plastid nor mitochondrial, over the assembly
            r = (L["plastid"]["bp"] / pl) / ((L["bases"] - L["plastid"]["bp"] - L["mitochondrion"]["bp"]) / span); cells_d.append(f"{r:,.0f}" if r >= 10 else f"{r:.1f}")
    who = {"elysia-leucolegnote": "two other individuals", "elysia-chlorotica": "two adults", "elysia-timida": "assembled individuals", "elysia-papillosa": "assembled individual and a second individual"}.get(slug, "assembled individual")
    floor = " or more" if slug in ("elysia-timida", "elysia-leucolegnote") else ""
    carried_rows += (f'<tr><td><a href="{slug}/kleptoplasts.html">{it(s_["species"])}</a></td><td>{re.sub(r" [(]in two other individuals[)]", "", tx.get("card", ""))}</td><td>{who}</td>'
                     f'<td class="num">{" and ".join(cells_f)}{floor}</td><td class="num">{" and ".join(cells_d)}</td></tr>')
# ---- the gene-transfer question in one paragraph, from the collected results
def _long_libs():
    for slug_, T_ in TRANSFER.items():
        for R_ in (T_ or {}).get("reads", []):
            if R_.get("shapes"): yield slug_, R_
_ins = {}; _cells = []
for slug_, R_ in _long_libs():
    n_ = sum(1 for c in R_["shapes"].get("candidates", []) if c["shape"] == "insertion" and c["insert"] == "cp")
    if n_: _ins[slug_] = _ins.get(slug_, 0) + n_
    sd_ = ((R_["shapes"].get("background") or {}).get("depth") or {}).get("spanning_depth") or {}
    # libraries whose figure carries a caveat of its own (E. timida: most of the plastid genome is unknown) are left out of the range
    if not n_ and sd_.get("500") and R_["kind"] == "hifi" and ("could_show_" + R_["library"]) not in (STEXT.get(slug_) or {}): _cells.append(sd_["500"])
_cr = jload(os.path.join(RG, "elysia-crispata.shapes.json")) or {}          # E. crispata: its page is built by scripts/31_build_pages.py
_sd = (((_cr.get("background") or {}).get("depth") or {}).get("spanning_depth") or {}).get("500")
if _sd and not _cr.get("candidates"): _cells.append(_sd)
# every long-read library: cell equivalents of DNA in reads long enough to show an insertion (for the methods page)
_alld = [x for x in [(((R_["shapes"].get("background") or {}).get("depth") or {}).get("spanning_depth") or {}).get("500") for slug_, R_ in _long_libs()] + [_sd] if x]
_all_cells_txt = f"ranges from {min(_alld) / 2:.0f} to {max(_alld) / 2:.0f}" if _alld else "is a few"
_ins_txt = ("in none of the libraries" if not _ins else "in only one library, the PacBio RS reads of " + ", ".join(it(BYSLUG[k]["species"].replace("Elysia", "E.")) for k in _ins) + f" ({word(sum(_ins.values()))} reads), which also contains reads with the reciprocal structure, nuclear DNA flanked by chloroplast DNA"
            if len(_ins) == 1 else "in the libraries of " + ", ".join(it(BYSLUG[k]["species"].replace("Elysia", "E.")) for k in _ins))
def _rng(a, b): return f"{a:,}" if a == b else f"{a:,} to {b:,}"
# how often chloroplast DNA is joined to slug DNA, against mitochondrial DNA, by kind of library: libraries with enough
# joins of both kinds to give a ratio (30 mitochondrial ones examined; 30 chloroplast ones, 10 in long reads). Pairs at the
# places listed as recurrent (bacterial ribosomal RNA genes) are left out, as on the species pages.
_ratio = {"short": [], "long": []}
for slug_, T_ in TRANSFER.items():
    for R_ in (T_ or {}).get("reads", []):
        r_ = (R_.get("chimera") or {}).get("rates") or {}; cp_, mt_ = r_.get("cp") or {}, r_.get("slugmt") or {}
        if R_["kind"] == "short":
            n_cp, n_mt = cp_.get("pairs_with_mate"), mt_.get("pairs_with_mate"); j_mt = mt_.get("pairs_also_with_slug", 0)
            j_cp = cp_.get("pairs_also_with_slug", 0) - sum(x["pairs"] for x in R_["chimera"].get("recurrent", []) if x["kind"] == "cp")
            enough = (mt_.get("examined") or {}).get("pairs_also_with_slug", j_mt) >= 30 and (cp_.get("examined") or {}).get("pairs_also_with_slug", j_cp) >= 30
        else:
            n_cp, n_mt, j_cp, j_mt = cp_.get("reads_with_segment"), mt_.get("reads_with_segment"), cp_.get("reads_also_with_slug", 0), mt_.get("reads_also_with_slug", 0)
            enough = j_mt >= 30 and j_cp >= 10
        if n_cp and n_mt and enough: _ratio["short" if R_["kind"] == "short" else "long"].append((j_cp / n_cp) / (j_mt / n_mt))
def _r1(v): return f"{min(v):.1f}" if f"{min(v):.1f}" == f"{max(v):.1f}" else f"{min(v):.1f} to {max(v):.1f}"
_ratio_txt = ""
if _ratio["short"] and _ratio["long"]:
    _long_txt = "as frequently as mitochondrial DNA is" if all(0.9 <= x <= 1.1 for x in _ratio["long"]) else f'{_r1(_ratio["long"])} times as frequently as mitochondrial DNA is'
    # the piece at the end of the fragment in mates that read through a join, pooled over the short-read libraries with ten or more mitochondrial ones
    _ep = {"cp": [0, 0], "slugmt": [0, 0]}
    for slug_, T_ in TRANSFER.items():
        for R_ in (T_ or {}).get("reads", []):
            e_ = (R_.get("chimera") or {}).get("split_end_piece") or {}
            if R_["kind"] != "short" or sum((e_.get("slugmt") or {}).values()) < 10: continue
            for k_ in _ep:
                _ep[k_][0] += (e_.get(k_) or {}).get("organelle", 0); _ep[k_][1] += sum((e_.get(k_) or {}).values())
    _ep_txt = ""
    if _ep["slugmt"][1] >= 30 and _ep["cp"][1] >= 30 and _ep["slugmt"][0] / _ep["slugmt"][1] < 0.75 * _ep["cp"][0] / _ep["cp"][1]:
        _ep_txt = (f'Insertions in somatic cells would produce such an excess, but the reads that span a junction indicate that it arises during library preparation. A read of this kind starts with the short segment at the end of its fragment. '
                   f'This terminal segment is chloroplast DNA in {_ep["cp"][0]} of {_ep["cp"][1]} junctions with chloroplast DNA, but mitochondrial DNA in only {_ep["slugmt"][0]} of {_ep["slugmt"][1]} junctions with mitochondrial DNA. '
                   'If the junctions were insertions in chromosomes that had been fragmented at random, the proportion would be the same for both organelles. '
                   'The difference shows that the fragments of the two organelle DNAs differ in size in these libraries, so the two are not joined to other DNA at the same rate during library preparation, and the excess cannot be interpreted as evidence of insertions '
                   '(the counts and the controls are given on the <a href="elysia-leucolegnote/gene-transfer.html"><i>E. leucolegnote</i></a> page). ')
    else:
        _ep_txt = ('Insertions in somatic cells would produce such an excess, but so would a library preparation that joins short, degraded fragments of DNA more readily than intact ones. '
                   'The reads cannot distinguish between these explanations, but an excess that depends on the type of library is more likely to arise from the library. ')
    _ratio_txt = (f'In the long-read libraries, chloroplast DNA is joined to nuclear DNA {_long_txt}. '
                  f'In the {word(len(_ratio["short"]))} short-read libraries that contain enough junctions of both kinds for a comparison, it is joined {_r1(_ratio["short"])} times as frequently. ' + _ep_txt)
_cell_txt = ((" The PacBio HiFi libraries of the other species contain no such reads. This sets only a weak upper limit on the frequency of insertions. Each read derives from a single cell, so a library that covers the genome "
              f"{_rng(round(min(_cells)), round(max(_cells)))} times with reads long enough to reveal an insertion has sampled DNA equivalent to {_rng(round(min(_cells) / 2), round(max(_cells) / 2))} cells. "
              f"About three reads containing an insertion would be expected if there were one insertion per {_rng(round(min(_cells) / 6), round(max(_cells) / 6))} cells. "
              'In tobacco, where the rate has been measured, chloroplast DNA is transferred to the nucleus in about one in 16,000 pollen grains (<a href="https://doi.org/10.1038/nature01435">Huang et al. 2003</a>) '
              'and about one in five million leaf cells (<a href="https://doi.org/10.1073/pnas.1430924100">Stegemann et al. 2003</a>). Insertions at such frequencies are far below what can be detected by sequencing a whole animal.') if _cells else "")
transfer_summary = (
    '<h3 style="margin:30px 0 8px">Has algal DNA been transferred to the slug genome?</h3>'
    '<p style="max-width:80ch">There is no evidence of algal DNA in the germline genome of any of these species. No assembly contains algal DNA that is supported by its reads, and egg DNA of <i>E. chlorotica</i> contains none. '
    'The same egg library does contain a fragment of mitochondrial DNA inserted in a chromosome, which is spanned by two independent molecules, so the approach can detect an inherited insertion of organelle DNA.</p>'
    '<p style="max-width:80ch">Kleptoplasts are located in somatic cells, however, so an insertion of chloroplast DNA would be present in one cell or in a small clone of cells, and would be represented by a single molecule among many from the same locus that lack it. '
    'All reads of every library were therefore searched for single molecules in which chloroplast DNA is joined to slug nuclear DNA. Such chimeric molecules were found in every library. '
    'However, every library also contains chimeras that cannot reflect the genome, for example between the slug’s mitochondrial DNA and its chromosomes, or between two different chromosomes. '
    'Chimeras of this kind form during library preparation, and a genuine junction that is sequenced only once cannot be distinguished from them.</p>'
    f'<p style="max-width:80ch">{_ratio_txt}</p>'
    f'<p style="max-width:80ch">Reads in which chloroplast DNA is flanked on both sides by contiguous nuclear sequence, the one structure that the joining of two fragments does not produce, were found {_ins_txt}.{_cell_txt}</p>'
    '<p style="max-width:80ch">The assemblies themselves contain chloroplast DNA at the ends of nuclear contigs (<i>E. crispata</i>, <i>E. chlorotica</i>, <i>E. leucolegnote</i>), which can be traced to one or two chimeric molecules of this kind where the reads are public; '
    'chloroplast DNA as separate sequences (<i>E. chlorotica</i>, <i>E. marginata</i>); a fragment of an algal mitochondrial genome (<i>E. papillosa</i>); and bacterial DNA, up to 2 Mb of it (<i>E. chlorotica</i>, <i>E. marginata</i>, <i>E. leucolegnote</i>). '
    'The gene-transfer page of each species describes the candidate molecules, the comparisons and the limits of detection.</p>')
carried_html = "" if not carried_rows else f"""
<section class="block"><div class="wrap">
<h2>Kleptoplasts of the sequenced animals</h2>
<p>The plastid DNA among the reads of each animal identifies the algae on which it had fed, and shows how much plastid DNA the animal contained when it was sequenced.</p>
<div class="tablewrap"><table class="data"><thead><tr><th>Species</th><th>Algal source of the kleptoplasts</th><th>Source of reads</th><th class="num">Plastid DNA as a share of sequenced bases</th><th class="num">Plastid genome copies per nuclear genome copy</th></tr></thead>
<tbody>{carried_rows}</tbody></table></div>
<p class="small muted" style="margin-top:10px;max-width:92ch">The last column is the read depth of the plastid genome divided by the read depth of the slug nuclear genome. It depends on the tissue that was sequenced and on how recently the animal had fed. For <i>E. timida</i> only 49 kb of the algal plastid genome has been sequenced, so the values refer to that part of the genome. For <i>E. leucolegnote</i> the values include only plastid DNA that resembles a reference sequence, which excludes much of the DNA from the other algae that these animals had eaten.</p>
{transfer_summary if all_transfer else ''}
</div></section>
"""
home = f"""
<section class="home-head"><div class="wrap">
<h1>Circadian clock genes in the genomes of {word(n_sp)} <i>Elysia</i> species</h1>
<p class="lede">Sea slugs of the genus <i>Elysia</i> retain chloroplasts from the algae they feed on, in some species for months. This site catalogues the circadian clock genes of every <i>Elysia</i> species with a sequenced genome, identifies the algal source of the kleptoplasts in each sequenced animal, and tests whether algal DNA has been transferred to the slug genome.</p>
</div></section>

<section class="block"><div class="wrap">
<h2>Species</h2>
<p>A nuclear genome assembly is available for {word(n_sp)} species. Select a card to open the pages for that species.</p>
<div class="sp-grid">
{"".join(card(s) for s in order)}
</div>
</div></section>

{carried_html}
<section class="block"><div class="wrap">
<div class="cols">
<div><h3>Clock genes across the genus</h3>
<p>Of the 47 clock genes, {len(all_intact)} are complete in the reference assembly of every species. {HOME_CORE} The <a href="compare.html">comparison page</a> gives the status of every gene in every species, the protein alignments and a species tree.</p></div>
<div><h3>Interpreting missing genes</h3>
<p>The assemblies differ far more in quality than the species differ in gene content. In most cases, a gene shown as partial or not found lies in a gap or a poorly assembled region of that assembly. The clock-gene page of each species describes the limitations of its assembly.</p></div>
<div><h3>Scope and limitations</h3>
<p>No rhythms were measured for this site. Gene functions are inferred from homology with <i>Drosophila</i> and mammals, and each genome was sequenced from one or a few individuals.</p></div>
</div>
</div></section>

<section class="block"><div class="wrap">
<h2>Other <i>Elysia</i> species with sequence data but no genome assembly</h2>
<p>This list is based on a survey of the public sequence archives carried out on 13 November 2026. For these species only RNA-seq or shotgun reads are available, or a genome project is registered but has not yet released data, so they are not covered on this site.</p>
<div class="tablewrap"><table class="data"><thead><tr><th>Species</th><th>Genomic reads</th><th>RNA-seq</th><th>Genome project</th></tr></thead>
<tbody>{"".join(others)}</tbody></table></div>
</div></section>
"""
open(os.path.join(SITE, "index.html"), "w").write(chrome.page("Elysia clocks: circadian clock genes in the sea slug genus Elysia", home, root="", site_current="index.html", today=today,
                                                           scripts='<script src="assets/pages.js"></script>'))

# ------------------------------------------------------------------ methods page (site-wide)
rows = "".join(f'<tr><td>{it(s["species"])}</td><td class="mono">{e(s["accession"])}</td><td>{e((CAT[s["slug"]]["assembly"]["name"] if s["slug"] in CAT else "xgElyCris1.1"))}</td><td>{LEVEL_NOTE[s["level"]]}</td>'
               f'<td>{ {"genbank": "GenBank", "ensembl": "Ensembl rapid release", "senckenberg": "authors’ website, not in GenBank", "none": "not publicly available"}[s["annotation"]] }</td>'
               f'<td>{("annotated genes and protein alignments" if (s["slug"] in CAT and CAT[s["slug"]]["annotated"]) else ("curated catalogue (four assemblies)" if s["slug"] == "elysia-crispata" else "protein alignments"))}</td></tr>' for s in order)
# place names as iNaturalist records them, written out for the table
PLACE = {"9HMM+576, Vixville Dr, West Bay, Cayman Islands": "West Bay, Cayman Islands", "Hants County, NS, Canada": "Hants County, Nova Scotia, Canada",
         "Qingdao, CN-SD, CN": "Qingdao, Shandong, China", "Monroe County, FL, USA": "Monroe County, Florida, USA"}
photo_rows = "".join(f'<tr><td>{it(c["species"])}</td><td><a href="{c["source"]}">{e(c["photographer"])}</a></td><td><a href="{c["licence_url"]}">{e(c["licence"])}</a></td><td>{e(PLACE.get(c["place"], c["place"]))}, {e(chrome.long_date(c["date"]))}</td></tr>' for c in credits.values())
shutil.copy(os.path.join(RG, "clock_matrix.tsv"), os.path.join(SITE, "downloads", "clock_matrix.tsv"))
if GN.get("tree"): open(os.path.join(SITE, "downloads", "species_tree.nwk"), "w").write(GN["tree"] + "\n")
dl = [("clock_matrix.tsv", "Gene status by species (TSV)", "status, completeness, identity to E. crispata and location of every gene in every species"),
      ("species_tree.nwk", "Species tree (Newick)", "inferred from the clock proteins that are complete in every species")]
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
    for suffix, label in (("_kleptoplast.json", "kleptoplast results (JSON)"), ("_gene_transfer.json", "gene-transfer results (JSON)")):
        if os.path.exists(os.path.join(SITE, "downloads", s_["slug"] + suffix)): dl.append((s_["slug"] + suffix, f"{s_['species']} {label}", ""))
    for f in sorted(os.listdir(os.path.join(SITE, "downloads"))):       # plastid contigs and barcode sequences (45_export_tables.py)
        m = re.match(re.escape(s_["slug"]) + r"(?:_(\w+?))?_(plastid_contigs|barcodes)\.fa$", f)
        if m: dl.append((f, f"{s_['species']} {'plastid contigs' if m.group(2) == 'plastid_contigs' else 'rbcL and tufA sequences'} (FASTA)",
                         {"adult": "second adult, short reads", "hainan": "Hainan individual", "papi1": "second individual", None: ("Shenzhen individual" if s_["slug"] == "elysia-leucolegnote" else "")}.get(m.group(1), m.group(1) or "")))
reads_methods = "" if not read_rows else f"""
<h2 id="reads">Identification of the kleptoplasts from the reads</h2>
<p>For <i>E. crispata</i>, all 47 Gb of reads were analysed (see the <a href="elysia-crispata/methods.html">methods page for that species</a>). For the other species, the read sets listed below were analysed in full.</p>
<div class="tablewrap"><table class="data"><thead><tr><th>Species</th><th>Source of DNA</th><th>Platform</th><th>Runs</th><th class="num">Bases analysed</th></tr></thead><tbody>{read_rows}</tbody></table></div>
<details class="method" open><summary>1. Reference panel</summary>
<p>The reference panel consists of 1,271 algal plastid sequences (241 complete plastid genomes of green algae; plastid records of Dasycladales, Xanthophyceae and Cladophorales; and a land plant, a diatom, a red alga and a cyanobacterium to detect more distant matches), 117 algal mitochondrial genomes, 7,352 nuclear ribosomal DNA sequences of algae, six slug ribosomal DNA sequences and the mitochondrial genomes of nine <i>Elysia</i> species. Low-complexity sequence in the panel was masked for the first pass.</p></details>
<details class="method"><summary>2. First pass: mapping all reads to the panel</summary>
<p>Every run was processed in full, in batches, and each read was mapped to the panel with minimap2 (k = 15 and w = 10 for long reads; the short-read preset for read pairs). Runs of up to 11 GB were downloaded from the SRA; larger runs were streamed from the archive in several parallel streams. Long reads with a match to the panel were retained. Two classes of read were subsampled: reads whose only match is a short local alignment to a panel window that attracts such alignments from very many reads (slug repeats that resemble part of a panel sequence), and reads that match only ribosomal DNA, at an identity typical of the slug’s own ribosomal DNA.</p>
<p>Read pairs were mapped only to the organelle genomes of the panel and to the assembly loci under test (see below). Every pair with a match was retained, which includes every pair that could link organelle DNA to other DNA. The one exception is the library from the second adult of <i>E. chlorotica</i>, in which about one seventh of the DNA is chloroplast DNA, amounting to tens of millions of pairs. For this library every pair was mapped and counted, but only one in 50 of the pairs that match a single type of organelle genome and nothing else was retained, selected by its ordinal position in the run. All pairs that match a locus under test or two types of sequence were retained. In addition, one read in 100 to 1,000, selected by its ordinal position in the run, was retained irrespective of what it matched, to provide a random sample of each library.</p></details>
<details class="method"><summary>3. Second pass: base-level alignment</summary>
<p>The retained reads were aligned to the assembly of the species together with the panel. Each read was represented by its best set of non-overlapping alignments. A segment was classified as plastid, mitochondrial or nuclear when it matched at 97% identity or higher (82% for the older PacBio RS reads of <i>E. chlorotica</i>). Assembly sequences that consist entirely of organelle DNA were treated as organelle sequences, and plastid-like regions within other assembly sequences were masked for this step, so that a plastid read would not be assigned to the nuclear genome because of a plastid fragment retained in the assembly.</p></details>
<details class="method"><summary>4. Identification of the alga</summary>
<p>Plastid reads were assigned to the reference sequence that they match best. Two barcode genes, rbcL and tufA, were then compared with 11,600 GenBank records of Ulvophyceae and Xanthophyceae. For PacBio HiFi libraries, the barcode sequence was extracted from every read that spans the gene, identical sequences were counted, and each distinct sequence was compared with the GenBank records. For short-read libraries, the plastid read pairs were first assembled into contigs with MEGAHIT, and the barcode sequences were extracted from the contigs.</p></details>

<h2 id="transfer">Search for algal DNA in the slug genomes</h2>
<details class="method" open><summary>1. DNA-level search of the assemblies</summary>
<p>Each assembly was divided into 10 kb windows with 1 kb overlap, and the windows were aligned to the panel. Sequences in which organelle matches cover 50% or more of the length are reported as organelle DNA. For all other sequences, every match of at least 200 bp with at least 80% identity is listed, after removal of matches to low-complexity sequence, matches to panel windows that are hit from 15 or more loci in the assembly, and matches to the slug’s own ribosomal DNA.</p></details>
<details class="method"><summary>2. Protein-level search of the assemblies</summary>
<p>A set of 471 proteins was aligned to each assembly with miniprot, retaining weak and partial alignments. The set comprises the plastid-encoded proteins of <i>Udotea flabellum</i>, photosynthesis proteins of <i>Chlamydomonas</i> and <i>Arabidopsis</i> (both nuclear- and plastid-encoded), UniProt proteins of siphonous green algae and all available proteins of <i>Vaucheria litorea</i>. Each matching locus was translated and compared with animal, algal, cyanobacterial, bacterial and fungal proteins using DIAMOND. A locus was classified as a slug gene when its best match among animals other than sea slugs scored at least as high as its best algal or microbial match, or within 10% of that score when a protein from another sea slug scored higher than both (as happens for highly conserved proteins such as GAPDH). The gene set published with an assembly was never used as evidence of animal origin for that assembly. The published gene sets of <i>E. chlorotica</i> and <i>E. marginata</i> include genes located on bacterial and plastid sequences, and such genes would otherwise support their own classification as animal genes.</p>
<p>Every assembly sequence that contained a locus of non-animal type was then translated in full in all six reading frames (DIAMOND blastx with frameshift-aware alignment and range culling), and each protein-coding region was classified in the same way. A sequence was classified as bacterial when three or more of its coding regions, and at least half of the total, were of bacterial type. The closest bacterial genera were taken from the best-matching Swiss-Prot entries; they indicate related taxa, not the species present.</p></details>
<details class="method"><summary>3. Search of the reads</summary>
<p>An insertion of chloroplast DNA in a somatic cell would be represented by one molecule or a few. Every read in which organelle DNA is joined to slug nuclear DNA was therefore retained as a candidate and listed, and reads from the same locus that lack the junction were not taken as evidence against it. The question for each candidate is whether the read could have arisen in another way.</p>
<p><strong>Read structure.</strong> Each long read was represented by its best set of non-overlapping alignments to the assembly and the panel. A read was classified as having an insertion-like structure when it contained organelle DNA flanked on both sides by nuclear DNA, with the two flanks on the same assembly sequence, on the same strand, in the correct order and no more than 5 kb apart. An organelle segment too short to interrupt the alignment of the read to the nuclear sequence was detected as an insertion of at least 100 bases within that alignment (150 bases for the older reads) that the first pass had matched to an organelle sequence. Such a segment was then aligned to the organelle sequence at base level and retained if 80% of it aligned with at least 95% identity (75% for the older reads). This check is needed because a tandem repeat in the slug genome can share one seed per repeat unit with a tandem repeat in a plastid genome, which the first pass reports as a match. When the same organelle-like segment was found inserted at two different loci, it was classified as a mobile element that resembles organelle DNA. Reads with the same structure were also counted with mitochondrial DNA as the insert, with nuclear DNA from another locus as the insert, and with the reciprocal structure (nuclear DNA flanked by chloroplast DNA). A read with organelle DNA at one end only was classified as having a single junction.</p>
<p><strong>Frequency.</strong> The proportion of chloroplast-containing reads that are joined to nuclear DNA was compared with the corresponding proportion for mitochondrial DNA, and with the proportion of nuclear reads that are joined to nuclear DNA from a different assembly sequence or from more than 50 kb away. The last proportion was measured on the random sample of each library, aligned in the same way. For <i>E. atroviridis</i>, which was analysed before random samples were retained, it was measured on 40 evenly spaced blocks of the run, and for the long reads of <i>E. chlorotica</i> on all reads of 2 kb or longer.</p>
<p>For read pairs, two classes were counted: pairs with one read in organelle DNA and the other at a unique position in nuclear DNA (discordant pairs), and single reads that span a junction (split reads). These were compared with pairs from the random sample whose two reads map to different assembly sequences or more than 20 kb apart. Split reads were also classified by the type of DNA at the start of the read. A read starts at the end of its fragment, so the segment at the start of a split read is the short terminal segment of the fragment. If a chromosome carrying an insertion is fragmented at random, either type of DNA is equally likely to be the terminal segment; if two fragments are joined end to end, the terminal segment is more often the type of DNA that is more abundant among short fragments. Two controls were used to estimate the proportions expected under random fragmentation. In the first, 20,000 nuclear, mitochondrial and chloroplast reads were sampled and the first and last 60 bases of each (45 bases for 100-base reads) were tested separately, to confirm that a short segment is recognised as well at the end of a read as at its start. In the second, reads were simulated that start with a short segment of organelle DNA and continue as nuclear DNA, and these were mapped as in the first pass. Because such a fragment contains no other organelle DNA, the pair is lost unless the short segment is detected.</p>
<p><strong>Recurrence.</strong> Junctions were grouped by position. For long reads, the reads were also compared with each other directly, which detects junctions in repetitive DNA that cannot be assigned to a unique position. The 400 bases on either side of every junction were aligned between all pairs of reads. Two reads were considered to share a junction when the alignment extended across the junction for at least 100 bases on both sides, placed the junction of one read within 6 bases of the junction of the other (15 bases for the older reads), contained no gap of that size nearby, and was at least 98% identical on both sides (70% for the older reads, two of which agree at only about 75 to 80% of positions even when they derive from the same molecule). Reads from the same ZMW are passes over a single molecule and were counted once, as were reads that share an end, which are PCR duplicates of a single molecule. Pairs of reads that narrowly failed these criteria are reported as junctions that lie close together. Junctions in a repeat were also searched for in read pairs from the plastid side, by counting reads whose plastid alignment ends at the junction position and which continue into the repeat, and comparing this number with the number at control positions.</p>
<p>For read pairs, three lists were compiled: loci in the nuclear genome with three or more distinct pairs within 1 kb whose mates are organelle DNA; breakpoints that are spanned at the same base by split reads from two molecules; and any two molecules, whether discordant pairs or split reads, whose nuclear ends lie within 1 kb of each other and whose organelle ends lie within 1 kb of each other on the same organelle sequence. Duplicates of a single fragment, whether optical or PCR duplicates, were counted once. They were identified from the read sequences: two pairs were treated as duplicates when both reads were identical in one of their first three 40-base windows (allowing three mismatches and an offset of two bases), irrespective of their positions in the run. If junctions were distributed uniformly, fewer than one pair of molecules in ten libraries would lie this close together by chance. The egg library of <i>E. chlorotica</i> provides a positive control: two molecules span the same 88 bases of mitochondrial DNA inserted between two adjacent bases of a chromosome.</p>
<p><strong>PCR duplicates.</strong> An amplified library contains duplicates of its molecules, so the number of reads overestimates the number of molecules. In the random sample of each library, nuclear reads that share both ends with another read in the sample were counted. Because the sample is a known fraction of the library, this gives the proportion of reads that are duplicates.</p>
<p><strong>Detection limit.</strong> Each read derives from one chromosome copy in one cell. A library that covers the genome D times has therefore sampled DNA equivalent to about D/2 cells, with each locus sampled from D different cells. If each cell carries k insertions on average, each at a different locus, the library is expected to contain about k × D/2 reads that span an insertion. D was calculated from reads long enough to contain a 500-base insert with flanking sequence on both sides. When no such read is found, about three would have been expected for k = 6/D. This value is an average over all cells that contributed DNA to the library. For read pairs with reads of R bases and fragments of F bases, one insertion per cell adds about D × (R − 80) / R split reads with at least 40 bases on each side of the junction, and about D × (F − 160) / 2R discordant pairs. Here D is the depth of nuclear DNA: the number of bases in the library that can be placed on the assembly, divided by the number of non-N bases in the assembly.</p>
<p>For the long-read libraries analysed here, D/2 {_all_cells_txt} cells. The rate at which chloroplast DNA is transferred to the nucleus has been measured in tobacco plants carrying a marker gene in the chloroplast genome: the gene was transferred to the nucleus in about one in 16,000 pollen grains (<a href="https://doi.org/10.1038/nature01435">Huang et al. 2003</a>) and about one in five million leaf cells (<a href="https://doi.org/10.1073/pnas.1430924100">Stegemann et al. 2003</a>). Both estimates include only fragments large enough to carry a functional gene, so they are minimum values. Insertions at such frequencies could not be detected in any of these libraries. Sequencing a whole animal can detect somatic insertions only if they are frequent, of the order of one per few cells, or if an insertion is present in a clone of cells large enough to be sampled more than once.</p></details>
<details class="method"><summary>4. Testing loci in an assembly with reads</summary>
<p>Where an assembly contains an organelle-like region within another sequence, and at the locus of each <i>E. chlorotica</i> long read with an insertion-like structure, reads were aligned to the region together with 3 kb of flanking sequence on each side. At each boundary, the reads that span the boundary were counted, as were the reads whose alignment ends there and which continue into other sequence. A read that crosses a boundary but has a gap of 15 bases or more next to it (6 bases for HiFi reads) was counted separately, because it joins the two types of DNA at a different point from the one in the assembly. These loci had been added to the panel beforehand, so that their reads were retained in the first pass. At the loci of the long reads with an insertion-like structure, the short reads from the other animals were searched for pairs with one read within 1 kb of the locus and the other in organelle DNA, and for reads that extend from the locus into organelle DNA.</p></details>
"""
methods = f"""
<section class="page-head"><div class="wrap">
<h1>Methods and data</h1>
<p class="lede">This page describes how the clock genes were identified in each genome and compared, how the reads and assemblies were searched for algal DNA, and the sources of the data. The analysis of <i>E. crispata</i> is described in more detail on the <a href="elysia-crispata/methods.html">methods page for that species</a>.</p>
</div></section>
<section class="block article"><div class="wrap">
<h2>Genome assemblies</h2>
<div class="tablewrap"><table class="data"><thead><tr><th>Species</th><th>Accession</th><th>Assembly name</th><th>Level</th><th>Gene annotation</th><th>Source of clock gene models</th></tr></thead><tbody>{rows}</tbody></table></div>
<p class="small muted" style="margin-top:10px">The species and data sets were identified in a survey of the public sequence archives carried out on 13 November 2026.</p>

<h2>Identification of clock genes</h2>
<details class="method" open><summary>1. Reference gene set</summary>
<p>The starting point was the catalogue of 47 genes in <i>E. crispata</i>, each of which had been identified by reciprocal best hits and by phylogenetic analysis with <i>Drosophila</i>, mammalian and molluscan proteins.</p></details>
<details class="method"><summary>2. Candidate loci in each genome</summary>
<p>The 47 proteins were aligned to each genome with miniprot, allowing introns of up to 500 kb and retaining secondary alignments. For the two species with a gene annotation in GenBank (<i>E. chlorotica</i> and <i>E. marginata</i>), the annotated proteins most similar to each catalogue protein were added as candidates.</p></details>
<details class="method"><summary>3. Distinguishing orthologues from paralogues</summary>
<p>Each candidate was searched with DIAMOND against all <i>E. crispata</i> proteins from both annotations, together with the 47 catalogue proteins. A candidate was assigned to a gene only if that gene was its best match in the catalogue and no other <i>E. crispata</i> protein was a clearly better match. The gene model chosen for each gene was the retained candidate with the largest number of residues identical to the <i>E. crispata</i> protein; an annotated gene was preferred when it was within 2% of the best alignment. Partial alignments of one protein that lie in order on the same strand were merged into a single gene with a long intron, whereas partial alignments on different sequences were classed as a fragmented gene.</p></details>
<details class="method"><summary>4. Completeness</summary>
<p>The gene models from all species were aligned with MAFFT, together with the orthologues from three other gastropods. Completeness was defined as the proportion of conserved alignment positions (those occupied in at least half of the sequences) that are covered by a given model. A model was classed as intact if its completeness was at least 90% and its alignment to the genome contained no frameshift or stop codon. This measure does not treat the <i>E. crispata</i> model as a reference standard. However, a model obtained by aligning the <i>E. crispata</i> protein cannot include regions that the <i>E. crispata</i> model lacks.</p></details>
<details class="method"><summary>5. Comparison and species tree</summary>
<p>Identity between two proteins was calculated over the positions at which both have a residue. The species tree was inferred from the proteins that are complete in every species, using alignment positions with at most 20% gaps, with IQ-TREE under the LG+F+G4 model and 1,000 ultrafast bootstrap replicates.</p></details>
<details class="method"><summary>6. Protein domains</summary><p>Protein domains were identified by searching Pfam-A with hmmsearch, using the gathering thresholds.</p></details>
<p class="small muted">The scripts are in <code>scripts/genus/</code> in the repository and are numbered in the order in which they are run.</p>
{reads_methods}

<h2>Downloads</h2>
<ul class="downloads">{"".join(f'<li><a href="downloads/{n}" download>{e(l)}</a><span>{e(t)}</span></li>' for n, l, t in dl)}</ul>

<h2 id="photos">Photographs</h2>
<p>The photographs do not show the sequenced animals. They are taken from iNaturalist observations and are used under the licences shown. For two species no photograph was available under a fully open licence, so the photographs used carry a non-commercial licence.</p>
<div class="tablewrap"><table class="data"><thead><tr><th>Species</th><th>Photographer</th><th>Licence</th><th>Location and date</th></tr></thead><tbody>{photo_rows}</tbody></table></div>
</div></section>
"""
open(os.path.join(SITE, "methods.html"), "w").write(chrome.page("Methods and data – Elysia clocks", methods, root="", site_current="methods.html", today=today, scripts='<script src="assets/pages.js"></script>'))

# ------------------------------------------------------------------ old addresses
for old in ("kleptoplasts.html", "gene-transfer.html"):
    open(os.path.join(SITE, old), "w").write(chrome.redirect("elysia-crispata/" + old))
old_data = os.path.join(SITE, "data", "clock.js")
if os.path.exists(old_data): os.remove(old_data)
print("home, compare, methods and", len(CAT), "species pages written;", len(all_intact), "genes whole in every genome; linkage pairs", pairs_kept, "of", pairs_total)
