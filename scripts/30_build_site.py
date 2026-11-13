#!/usr/bin/env python3
"""Build the Elysia crispata pages of the site (docs/elysia-crispata/) from the result files.

usage: 30_build_site.py repo_dir [tree_svg_dir]
Reads  results/clock_genes.json, results/kleptoplast.json, results/hgt.json (the last two optional)
       and, when present, results/genus/genus.json (the same genes in the other species)
Writes docs/elysia-crispata/{index,kleptoplasts,gene-transfer,methods}.html,
       docs/data/elysia-crispata.js and copies the downloadable tables into
       docs/elysia-crispata/downloads/. The home page, the comparison page and the pages of
       the other species are written by scripts/genus/20_build_genus_site.py.
Every number on the pages is taken from the result files, so re-running the analysis and
then this script keeps the site in step with the data.
"""
import sys, os, json, shutil, html, datetime, glob

repo = sys.argv[1]
treedir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(repo, "results", "trees", "figures")
R = os.path.join(repo, "results"); SITE = os.path.join(repo, "docs"); DOCS = os.path.join(SITE, "elysia-crispata")
os.makedirs(os.path.join(SITE, "data"), exist_ok=True)
os.makedirs(os.path.join(DOCS, "downloads"), exist_ok=True)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "genus"))
import chrome

LABEL = {"Ck1de": "Ck1δ/ε", "Ck1a": "Ck1α", "Ck2a": "Ck2α", "Ck2b": "Ck2β", "Pp1a": "Pp1α", "Pp1b": "Pp1β", "Pp2aC": "Pp2a-C",
         "Pp2aB55": "Pp2a-B55", "Pp2aB56e": "Pp2a-B56ε", "Pp2aB56d": "Pp2a-B56δ", "Phr6-4": "Phr 6-4", "PhrCPD": "Phr CPD",
         "CryDASH1": "Cry-DASH 1", "CryDASH2": "Cry-DASH 2", "Ropsin1": "r-opsin 1", "Ropsin2": "r-opsin 2", "RopsinL1": "r-opsin-like 1",
         "RopsinL2": "r-opsin-like 2", "GoOpsin1": "Go-opsin 1", "GoOpsin2": "Go-opsin 2"}
CATS = {"core": "Core loop", "photolyase": "Photolyase family", "loop2": "Second loop", "kinase": "Kinase", "phosphatase": "Phosphatase",
        "degradation": "Degradation", "modulator": "Modulator", "output": "Output", "light": "Opsin"}
SPECIES = {"HUMAN": "Homo sapiens", "DROME": "Drosophila melanogaster", "MAGGI": "Magallana gigas", "APLCA": "Aplysia californica",
           "LOTGI": "Lottia gigantea", "BIOGL": "Biomphalaria glabrata", "POMCA": "Pomacea canaliculata", "MIZYE": "Mizuhopecten yessoensis",
           "ELYCH": "Elysia chlorotica", "ELYMA": "Elysia marginata", "PLAOC": "Plakobranchus ocellatus", "CAPTE": "Capitella teleta"}
ASM_SHORT = {"sanger": "xgElyCris1.1 primary (Sanger)", "alt": "xgElyCris1.1 alternate haplotype", "purdue": "Purdue ASM3367554v1", "tampa": "Tampa Elysia_crispata_v1"}
CHROMS = [("1", 90802811), ("2", 73756357), ("3", 56100361), ("4", 42579028), ("5", 30624988), ("6", 25209903), ("7", 23020430), ("8", 19671316),
          ("9", 19645690), ("10", 18345459), ("11", 17553211), ("12", 15649513), ("13", 14153947), ("14", 14121010), ("15", 13644184)]

def load(name):
    p = os.path.join(R, name)
    return json.load(open(p)) if os.path.exists(p) else None

clock = load("clock_genes.json")
klepto = load("kleptoplast.json")
hgt = load("hgt.json")
today = datetime.date.today().isoformat()
e = html.escape

# ------------------------------------------------------------------ data for the gene browser
import statistics
genes = clock["genes"]
trees = {}
for g in genes:
    g["label"] = LABEL.get(g["symbol"], g["symbol"])
    p = os.path.join(treedir, g["symbol"] + ".svg")
    if os.path.exists(p): trees[g["symbol"]] = open(p).read()
depths = [g["reads_support"]["median_reads"] for g in genes if g.get("reads_support")]
meta = dict(built=today, typical_exon_depth=statistics.median(depths) if depths else 0)
# the same genes in the other species (written by scripts/genus/12_compare_species.py)
across = {}
gp = os.path.join(R, "genus", "genus.json")
if os.path.exists(gp):
    GN = json.load(open(gp)); TAGS = {x["tag"]: x for x in GN["species"]}
    for sym, gg in GN["genes"].items():
        rows = []
        for sp in GN["species"]:
            if sp["slug"] == "elysia-crispata": continue
            cp = os.path.join(R, "genus", sp["slug"] + ".clock.json")
            if not os.path.exists(cp): continue
            og = next(x for x in json.load(open(cp))["genes"] if x["symbol"] == sym)
            rows.append(dict(name=sp["name"], href=f"../{sp['slug']}/index.html#gene/{sym}", status=og["status"], length=og.get("length", 0),
                             identity=gg["identity"].get("|".join(sorted((sp["tag"], "ELYCR"))))))
        across[sym] = rows
data = dict(genes=genes, not_found=clock["not_found"], assemblies=ASM_SHORT, categories=CATS, species=SPECIES,
            chromosomes=[dict(name=n, length=l) for n, l in CHROMS], trees=trees, meta=meta, across=across, compare_href="../compare.html")
with open(os.path.join(SITE, "data", "elysia-crispata.js"), "w") as fo:
    fo.write("window.CLOCK = " + json.dumps(data, separators=(",", ":")) + ";\n")

# ------------------------------------------------------------------ shared chrome
CRISPATA_SOURCES = ("<i>Elysia crispata</i>: genome assembly xgElyCris1.1 and reads of BioProject PRJEB67629 (Middlebrooks et al. 2026, Wellcome Open Research 11:148; "
                    "Aquatic Symbiosis Genomics, Wellcome Sanger Institute), released under the Tree of Life open data policy.")
def page(title, current, body, scripts=""):
    return chrome.page(title, body, root="../", species="Elysia crispata", tabs=chrome.SPECIES_TABS, tab_current=current, scripts=scripts, today=today, sources=CRISPATA_SOURCES,
                       description="Circadian clock genes, kleptoplast identity and a test for algal gene transfer in the lettuce sea slug Elysia crispata, from the xgElyCris1 genome data.")

def status_counts():
    c = {}
    for g in genes: c[g["assemblies"]["sanger"]["status"]] = c.get(g["assemblies"]["sanger"]["status"], 0) + 1
    return c

sc = status_counts()
n = len(genes)
n_ok = sc.get("complete", 0)
n_disc = sc.get("complete_discrepancy", 0)
n_bad = n - n_ok - n_disc
missing = [g for g in genes if g["assemblies"]["sanger"]["status"] in ("absent", "divergent")]
partial = [g for g in genes if g["assemblies"]["sanger"]["status"] in ("partial", "fragmented")]
other_complete = sum(1 for g in genes if g["assemblies"]["purdue"]["status"] == "complete" or g["assemblies"]["tampa"]["status"] == "complete")
low = sorted([g for g in missing + partial if g.get("reads_support")], key=lambda g: g["reads_support"]["median_reads"])
def names(gs): return ", ".join(g["label"] for g in gs)
typ = round(meta["typical_exon_depth"])
lowtxt = ""
if low:
    lo = [g for g in low if g["reads_support"]["median_reads"] < 0.4 * meta["typical_exon_depth"]]
    lowtxt = (f"Several of them sit in stretches the sequencing library barely sampled: {', '.join(g['label'] + ' (' + str(round(g['reads_support']['median_reads'])) + ' reads per exon)' for g in lo[:5])}, "
              f"against {typ} for a typical clock gene. The library was amplified by PCR from very little DNA, which leaves holes like these.") if lo else ""

index_body = f"""
<section class="intro"><div class="wrap">
<h1>The lettuce sea slug’s clock, gene by gene</h1>
<p class="lede"><i>Elysia crispata</i> keeps stolen chloroplasts photosynthesising for months. This is a catalogue of the circadian clock genes in its genome: where each one sits in the Sanger assembly xgElyCris1.1, whether the assembly has it whole, and the protein it encodes.</p>
<div class="clockwork">
<div id="schematic"></div>
<div class="clockwork-list" id="schematic-list"></div>
<div class="clockwork-legend" aria-label="Key to gene status in xgElyCris1.1">
<span><i class="glyph ok"></i>intact in xgElyCris1.1</span>
<span><i class="glyph warn"></i>whole, but the assembly sequence differs</span>
<span><i class="glyph part"></i>partial or in pieces</span>
<span><i class="glyph miss"></i>missing from the assembly</span>
</div>
</div>
<p class="clockwork-caveat">Select a gene to open it. The wiring is the textbook fly and mouse clock; which gene does what in <i>Elysia</i> has not been tested, and one RNA sample cannot show rhythm.</p>
</div></section>

<section class="block"><div class="wrap">
<div class="cols">
<div><h3>A complete clock, by homology</h3>
<p>{n} genes: both halves of the core activator (Clock and Bmal1), the repressors Period, Timeless and a mammalian-type Cryptochrome, a light-sensitive Cryptochrome, the second loop (Pdp1, Vrille, E75, Hr3, Clockwork orange), fifteen kinases, phosphatases and degradation factors, ten opsins and four photolyases. Each call rests on reciprocal best hits and a family tree.</p></div>
<div><h3>The assembly misses some of it</h3>
<p>{n_ok} genes are intact in xgElyCris1.1 and {n_disc} more are whole but with a frameshift or other mismatch in the assembly sequence. {len(partial)} are partial or scrambled ({names(partial)}) and {len(missing)} are absent ({names(missing)}).</p></div>
<div><h3>Missing from the assembly, not the animal</h3>
<p>All {n} genes are whole in at least one of the two other <i>E. crispata</i> assemblies, and reads from the sequenced animal itself cover their exons. {lowtxt}</p></div>
</div>
</div></section>

<section class="block"><div class="wrap">
<h2>Where the genes sit</h2>
<p>The 15 chromosome-scale scaffolds of xgElyCris1.1, drawn to scale. Each mark is a clock gene; select one to open it.</p>
<div class="karyo-card"><div class="karyo" id="karyo"></div><div class="offmap" id="offmap"></div></div>
</div></section>

<section class="block"><div class="wrap">
<h2>All {n} genes</h2>
<div class="filters" id="filters"></div>
<div class="tablewrap"><table class="genes">
<thead><tr><th>Gene</th><th>Role</th><th class="num">Protein, aa</th><th>In xgElyCris1.1</th><th>Other assemblies</th><th class="num" title="Median number of HiFi reads from the sequenced animal covering each exon">Reads per exon</th><th class="num" title="Transcripts per million in the head RNA-seq library">Head RNA, TPM</th></tr></thead>
<tbody id="genes-body"></tbody></table></div>
<p class="small muted" style="margin-top:12px;max-width:90ch">Other assemblies: the alternate haplotype of the same animal, and two independent assemblies of other individuals (Purdue, GCA_033675545.1; Tampa, GCA_984581095.1). Reads per exon counts PacBio HiFi reads from the sequenced animal. Head RNA is one library from the same animal.</p>

<h3 style="margin-top:34px">Looked for and not found</h3>
<div class="tablewrap" style="margin-top:10px"><table class="data">
<thead><tr><th>Gene</th><th>Searched with</th><th>Reading</th></tr></thead>
<tbody>{''.join(f'<tr><td><strong>{e(x["symbol"])}</strong><br><span class="muted small">{e(x["name"])}</span></td><td class="mono">{e(x["seeds"].replace(",", ", "))}</td><td>{e(x["interpretation"])}</td></tr>' for x in clock["not_found"])}</tbody></table></div>
</div></section>
<dialog class="gene" id="gene" aria-label="Gene details"></dialog>
"""
open(os.path.join(DOCS, "index.html"), "w").write(page("Elysia crispata clock genes", "index.html", index_body, '<script src="../data/elysia-crispata.js"></script>\n<script src="../assets/app.js"></script>'))

# downloads
for f in ("clock_genes.tsv", "clock_genes.json"):
    if os.path.exists(os.path.join(R, f)): shutil.copy(os.path.join(R, f), os.path.join(DOCS, "downloads", f))
for f in glob.glob(os.path.join(R, "sequences", "*")):
    shutil.copy(f, os.path.join(DOCS, "downloads", os.path.basename(f)))
print("index.html written:", n, "genes;", n_ok, "intact,", n_disc, "discrepant,", len(partial), "partial,", len(missing), "missing")

# other pages are written by 31_build_pages.py once the kleptoplast and gene-transfer results exist
if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)), "31_build_pages.py")):
    import runpy
    runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "31_build_pages.py"),
                   init_globals=dict(repo=repo, R=R, DOCS=DOCS, page=page, e=e, genes=genes, klepto=klepto, hgt=hgt, today=today, meta=meta))
