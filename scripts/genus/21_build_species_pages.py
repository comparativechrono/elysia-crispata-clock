#!/usr/bin/env python3
"""Build the "Whose chloroplasts" and "Gene transfer" pages of every species except E. crispata
(whose pages are written by scripts/31_build_pages.py).

usage: 21_build_species_pages.py REPO
Reads  results/genus/<slug>.kleptoplast.json and <slug>.transfer.json (39_collect_species.py),
       config/species.tsv and config/species_text.json (the conclusions, written by hand once the
       numbers were in; every number in them is filled in from the result files at build time).
Writes docs/<slug>/kleptoplasts.html, docs/<slug>/gene-transfer.html and copies the result files
       to docs/downloads/.
"""
import sys, os, json, csv, html, datetime, shutil, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chrome
repo = os.path.abspath(sys.argv[1]); RG = os.path.join(repo, "results", "genus"); SITE = os.path.join(repo, "docs")
e = html.escape; today = datetime.date.today().isoformat()
SP = list(csv.DictReader(open(os.path.join(repo, "config", "species.tsv")), delimiter="\t"))
TEXT = json.load(open(os.path.join(repo, "config", "species_text.json"))) if os.path.exists(os.path.join(repo, "config", "species_text.json")) else {}
def jload(p): return json.load(open(p)) if os.path.exists(p) else None
SPAN = {}
for _sp in SP:
    _c = jload(os.path.join(RG, _sp["slug"] + ".clock.json"))
    if _c: SPAN[_sp["slug"]] = _c["assembly"]["span"]
def it(s): return f"<i>{e(s)}</i>"
def num(x): return f"{x:,.0f}"
def pct(x, d=1): return f"{100 * x:.{d}f}%"
def pct_auto(x):
    """percentage with enough digits to show a small share"""
    v = 100 * x
    if 99.95 <= v < 100: return f"{v:.2f}%"
    return f"{v:.1f}%" if v >= 10 else (f"{v:.2f}%" if v >= 0.1 else (f"{v:.3f}%" if v >= 0.001 else ("0%" if v == 0 else f"{v:.4f}%")))
def gb(x): return f"{x / 1e9:.1f} Gb" if x >= 1e9 else f"{x / 1e6:.0f} Mb"
def plural(n, one, many=None): return one if n == 1 else (many or one + "s")
WORDS = ["no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve"]
def word(k): return WORDS[k] if 0 <= k < len(WORDS) else f"{k:,}"
def taxon_html(name):
    if re.match(r"^[A-Z][a-z]+ [a-z.]+$", name): return it(name)
    return e(name)
def hbar(rows, unit="%", width=760, label_w=250, digits=2):
    maxv = max(r[1] for r in rows) * 1.18 or 1
    rowh, top = 34, 8; H = top + rowh * len(rows) + 6; pw = width - label_w - 90
    out = [f'<svg viewBox="0 0 {width} {H}" role="img" aria-label="{e("; ".join(f"{r[0]}: {r[1]:.{digits}f}{unit}" for r in rows))}">',
           f'<line class="gridline" x1="{label_w}" x2="{label_w}" y1="{top}" y2="{H - 6}"/>']
    for i, r in enumerate(rows):
        y = top + i * rowh; w = max(2, r[1] / maxv * pw)
        out.append(f'<text class="lbl" x="{label_w - 10}" y="{y + 20}" text-anchor="end">{e(r[0])}</text>')
        out.append(f'<rect class="bar" x="{label_w}" y="{y + 7}" width="{w:.1f}" height="18" rx="4" data-tip="{e(r[2])}" tabindex="0"><title>{e(r[2])}</title></rect>')
        out.append(f'<text class="val" x="{label_w + w + 8:.1f}" y="{y + 20}">{pct_auto(r[1] / 100) if unit == "%" else f"{r[1]:.{digits}f}{unit}"}</text>')
    return "\n".join(out) + "</svg>"
TIP_JS = '<script src="../assets/pages.js"></script>'
KIND = {"hifi": "accurate long reads", "clr": "older long reads", "short": "short read pairs"}
CLS_WORD = {"cp": "plastid", "slugmt": "slug mitochondrial", "algmt": "algal mitochondrial", "algrDNA": "algal nuclear ribosomal", "bait": "bait"}

def fmt(text, ctx):
    """fill {name} placeholders of a hand-written sentence from the numbers of this species"""
    try: return text.format(**ctx)
    except KeyError as ex: raise SystemExit(f"config/species_text.json uses an unknown number {ex} in: {text[:80]}")

def has_page(slug, which):
    """a page is written only when its result file exists and its conclusion has been written"""
    f, key = {"klepto": (".kleptoplast.json", "klepto_answer"), "transfer": (".transfer.json", "transfer_answer")}[which]
    return os.path.exists(os.path.join(RG, slug + f)) and bool(TEXT.get(slug, {}).get(key))
def tabs_for(slug):
    t = [("index.html", "Clock genes")]
    if has_page(slug, "klepto"): t.append(("kleptoplasts.html", "Whose chloroplasts"))
    if has_page(slug, "transfer"): t.append(("gene-transfer.html", "Gene transfer"))
    return t

def context(sp, K, T):
    """numbers available to the hand-written text"""
    c = dict(species=sp["species"], short="E. " + sp["species"].split()[1], assembly=sp["asm_name"])
    if K:
        for L in K["libraries"]:
            p = "" if L["library"] == "main" else L["library"] + "_"
            c[p + "bases"] = gb(L.get("bases_screened") or L["bases"]); c[p + "reads"] = num(L["reads"])
            c[p + "cp_pct"] = pct_auto(L["plastid"]["fraction"]); c[p + "cp_reads"] = num(L["plastid"]["reads"] or 0); c[p + "cp_mb"] = f'{L["plastid"]["bp"] / 1e6:,.1f}'
            c[p + "mt_pct"] = pct_auto(L["mitochondrion"]["fraction"]); c[p + "mt_reads"] = num(L["mitochondrion"]["reads"] or 0)
            c[p + "rdna_reads"] = num(L["algal_rdna"]["reads"] or 0); c[p + "algmt_reads"] = num(L["algal_mitochondrion"]["reads"] or 0)
            c[p + "algmt_kb"] = f'{L["algal_mitochondrion"]["bp"] / 1e3:,.0f}'
            c[p + "cp_per_algmt"] = num(L["plastid"]["bp"] / L["algal_mitochondrion"]["bp"]) if L["algal_mitochondrion"]["bp"] else "–"
            # depth of the plastid genome against depth of the slug genome
            pl = TEXT.get(sp["slug"], {}).get("plastome_length"); span = SPAN.get(sp["slug"])
            if span:
                dn = L["bases"] / span; ds = (L.get("bases_screened") or L["bases"]) / span
                c[p + "nuc_depth"] = num(dn) if dn >= 10 else f"{dn:.1f}"; c[p + "screened_depth"] = num(ds) if ds >= 10 else f"{ds:.1f}"
            if pl and span and L["plastid"]["bp"]:
                dp = L["plastid"]["bp"] / pl
                c[p + "cp_depth"] = num(dp) if dp >= 10 else f"{dp:.1f}"; c[p + "depth_ratio"] = num(dp / dn) if dp / dn >= 10 else f"{dp / dn:.1f}"
            for i, d in enumerate(L["donors"][:4]):
                c[f"{p}donor{i + 1}"] = d["name"]; c[f"{p}donor{i + 1}_pct"] = pct_auto(d["fraction"]); c[f"{p}donor{i + 1}_id"] = f'{100 * d["median_identity"]:.1f}%'
                c[f"{p}donor{i + 1}_reads"] = num(d["reads"])
            for b in L["barcodes"]:
                k = f'{p}{b["barcode"]}'
                if k + "_id" not in c:
                    c[k + "_id"] = f'{b["identity"]:.1f}%' if b.get("identity") is not None else "–"; c[k + "_hit"] = (b.get("top_hit") or "").split(" [")[0]
                    c[k + "_len"] = num(b.get("alnlen") or 0); c[k + "_reads"] = num(b.get("reads") or 0)
                    ms = b.get("most_similar") or {}
                    c[k + "_best_id"] = f'{ms.get("identity", 0):.1f}%'; c[k + "_best_len"] = num(ms.get("longest_alnlen", 0)); c[k + "_best_species"] = ms.get("species", "")
            if L.get("contigs"): c[p + "contigs"] = num(L["contigs"]["n"]); c[p + "contig_kb"] = f'{L["contigs"]["bp"] / 1e3:,.0f}'; c[p + "contig_longest_kb"] = f'{L["contigs"]["longest"] / 1e3:,.1f}'
    if T:
        s = T["screen"]
        c["whole_n"] = word(len(s["whole"])); c["whole_bp"] = num(sum(w["length"] for w in s["whole"])); c["loci_n"] = word(len(s["loci"])); c["numt_n"] = word(len(s["numt_like"]))
        c["repeat_matches"] = num(s["removed"]["repeat_like_matches"]); c["own_rdna"] = num(s["removed"]["own_rdna_loci"])
        for R in T["reads"]:
            p = "" if R["library"] == "main" else R["library"] + "_"
            r = R["chimera"]["rates"]
            for k, v in r.items():
                tot = v.get("reads_with_segment", v.get("pairs_with_mate", 0)); ch = v.get("reads_also_with_slug", v.get("pairs_also_with_slug", 0))
                c[f"{p}{k}_n"] = num(tot); c[f"{p}{k}_chim"] = num(ch); c[f"{p}{k}_rate"] = pct_auto(ch / tot) if tot else "–"
            c[p + "recurrent"] = word(len(R["chimera"].get("recurrent", [])))
        if T.get("protein_screen"):
            ps = T["protein_screen"]["summary"]; c["prot_loci"] = num(ps["loci"]); c["prot_queries"] = num(T["protein_screen"]["n_queries"])
            c["prot_animal"] = num(ps.get("animal gene", 0)); c["prot_alga"] = word(ps.get("alga-like", 0)); c["prot_microbe"] = word(ps.get("microbe-like", 0))
            cen = [x for x in T["protein_screen"].get("non_animal_contigs", []) if "census" in x and x["census"]["stretches"].get("microbe", 0) >= 3
                   and x["census"]["stretches"].get("microbe", 0) >= 0.5 * sum(v for kk, v in x["census"]["stretches"].items() if kk != "own_only")]
            c["bact_n"] = word(len(cen)); c["Bact_n"] = word(len(cen)).capitalize(); c["bact_kb"] = num(round(sum(x["length"] for x in cen) / 1000)); c["bact_mb"] = f'{sum(x["length"] for x in cen) / 1e6:.1f}'
            c["bact_genes"] = num(sum(x["census"]["stretches"].get("microbe", 0) for x in cen))
            c["bact_listed"] = word(sum(x["listed_in_own_annotation"] for x in T["protein_screen"].get("non_animal_contigs", []))); c["Bact_listed"] = c["bact_listed"].capitalize()
            c["genome_gc"] = f'{100 * T["protein_screen"]["genome_gc"]:.0f}%' if T["protein_screen"].get("genome_gc") else ""
    return c

def lib_label(L, n_libs): return (L.get("label") or (L["animal"][0].upper() + L["animal"][1:])) if n_libs > 1 else "The reads"

# =====================================================================================
def klepto_page(sp, K, T):
    slug = sp["slug"]; tx = TEXT.get(slug, {}); ctx = context(sp, K, T); libs = K["libraries"]; main = libs[0]
    lede = fmt(tx.get("klepto_lede", "The plastid DNA among the reads of {species} says which alga the sequenced animal had been eating."), ctx)
    answer = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get("klepto_answer", []))
    rows_lib = "".join(f'<tr><td>{e(L["animal"][0].upper() + L["animal"][1:])}</td><td>{e(L["platform"])}</td><td class="mono">{e(", ".join(L["runs"]))}</td><td class="num">{gb(L.get("bases_screened") or L["bases"])}</td><td>{e(L["note"])}</td></tr>' for L in libs)
    charts = ""
    for L in libs:
        rows = [("Algal plastid DNA", 100 * L["plastid"]["fraction"], f'{num(L["plastid"]["reads"] or 0)} {"reads" if L["kind"] != "short" else "pairs"}, {L["plastid"]["bp"] / 1e6:,.2f} Mb'),
                ("Slug mitochondrion", 100 * L["mitochondrion"]["fraction"], f'{num(L["mitochondrion"]["reads"] or 0)} {"reads" if L["kind"] != "short" else "pairs"}, {L["mitochondrion"]["bp"] / 1e6:,.2f} Mb')]
        # short reads match algal mitochondrial genomes in conserved genes too; shown only when there is enough to mean something
        if L["algal_mitochondrion"]["bp"] and (L["kind"] != "short" or L["algal_mitochondrion"]["fraction"] >= 5e-6): rows.append(("Algal mitochondrial DNA", 100 * L["algal_mitochondrion"]["fraction"], f'{L["algal_mitochondrion"]["bp"] / 1e3:,.1f} kb'))
        if L["algal_rdna"]["bp"] and L["kind"] != "short": rows.append(("Algal nuclear ribosomal DNA", 100 * L["algal_rdna"]["fraction"], f'{L["algal_rdna"]["bp"] / 1e3:,.1f} kb'))
        digits = 2 if max(r[1] for r in rows) >= 0.1 else 4
        charts += f'<div class="chart"><h3>Share of sequenced bases{": " + e(L.get("label") or L["animal"]) if len(libs) > 1 else ""}</h3><p class="sub">{gb(L.get("bases_screened") or L["bases"])} of {e(KIND[L["kind"]])} examined. Counted here: bases that match the reference at {"97" if L["kind"] != "clr" else "82"}% or better. Everything else is slug nuclear DNA or matches nothing in the panel.</p>{hbar(rows, digits=digits)}</div>'
    donor_tables = ""
    for L in libs:
        if not L["donors"] or L["plastid"]["bp"] < 5000: continue
        rows = "".join(f'<tr><td>{taxon_html(d["name"])} <span class="mono muted">{e(d["accession"])}</span></td><td class="num">{pct_auto(d["fraction"])}</td><td class="num">{num(d["reads"])}</td>'
                       f'<td class="num">{100 * d["median_identity"]:.1f}%</td><td class="num">{100 * d["identity_q10"]:.1f} to {100 * d["identity_q90"]:.1f}%</td></tr>' for i, d in enumerate(L["donors"][:8]) if i == 0 or d["fraction"] >= 0.0002)
        donor_tables += (f'<h3>{e(lib_label(L, len(libs)))}</h3>' if len(libs) > 1 else "") + f'''<div class="tablewrap"><table class="data">
<thead><tr><th>Plastid reference matched best</th><th class="num">Share of plastid DNA</th><th class="num">{"Reads" if L["kind"] != "short" else "Mates"}</th><th class="num">Identity, median</th><th class="num">Middle 80%</th></tr></thead>
<tbody>{rows}</tbody></table></div>'''
    bc_tables = ""
    for L in libs:
        if not L["barcodes"]: continue
        def brow(b):
            ms = b.get("most_similar") or {}
            src = f'{num(b["reads"])} reads ({pct(b["fraction"], 0)})' if b.get("reads") is not None else f'contig of {num(b.get("contig_length") or 0)} bp'
            top = b.get("top_hit") or "no match"; sp_ = top.split(" [")[0]; acc = top.split("[")[-1].rstrip("]") if "[" in top else ""
            near = f'{ms.get("identity", 0):.1f}% over {num(ms.get("longest_alnlen", 0))} bp: {e(ms.get("species", ""))}' if ms else "–"
            return (f'<tr><td>{b["barcode"]}</td><td class="num">{num(b["length"])}</td><td>{src}</td><td>{taxon_html(sp_)} <span class="mono muted">{e(acc)}</span></td>'
                    f'<td class="num">{b["identity"]:.1f}% of {num(b["alnlen"])} bp</td><td>{near}</td></tr>') if b.get("identity") is not None else ""
        bs = sorted(L["barcodes"], key=lambda b: (b["barcode"], -(b.get("reads") or b.get("length") or 0)))
        rows = "".join(brow(b) for b in bs[:14])
        how = ("Each row is one distinct sequence (haplotype) read directly out of the reads that span the gene." if L["kind"] == "hifi" else
               f'Cut out of the plastid contigs assembled from these reads ({num(L["contigs"]["n"])} contigs, {L["contigs"]["bp"] / 1e3:,.0f} kb in all).' if L.get("contigs") else "")
        bc_tables += (f'<h3>{e(lib_label(L, len(libs)))}</h3>' if len(libs) > 1 else "") + f'''<p class="small muted" style="max-width:86ch">{how}</p><div class="tablewrap"><table class="data">
<thead><tr><th>Gene</th><th class="num">Length, bp</th><th>From</th><th>Longest good match in GenBank</th><th class="num">Identity</th><th>Most similar record of 400 bp or more</th></tr></thead>
<tbody>{rows}</tbody></table></div>'''
    nuc = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get("klepto_nuclei", []))
    notes = "".join(f"<li>{fmt(p, ctx)}</li>" for p in tx.get("klepto_caveats", []))
    which = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get("klepto_which", []))
    # downloads made by 45_export_tables.py
    dls = []
    for L in libs:
        tag = slug if L["library"] == "main" else f'{slug}_{L["library"]}'; lab = (" (" + lib_label(L, len(libs)) + ")") if len(libs) > 1 else ""
        if os.path.exists(os.path.join(SITE, "downloads", tag + "_plastid_contigs.fa")): dls.append(f'<a href="../downloads/{tag}_plastid_contigs.fa">plastid contigs{lab}</a>')
        if os.path.exists(os.path.join(SITE, "downloads", tag + "_barcodes.fa")): dls.append(f'<a href="../downloads/{tag}_barcodes.fa">barcode sequences{lab}</a>')
    rhythm = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get("klepto_rhythm", []))
    if rhythm and dls: rhythm += f'<p>Downloads: {", ".join(dls)}, and <a href="../downloads/{slug}_kleptoplast.json">all the numbers on this page</a>.</p>'
    body = f"""
<section class="page-head"><div class="wrap">
<h1>Whose chloroplasts was this slug carrying?</h1>
<p class="lede">{lede}</p>
<div class="answer">{answer}</div>
</div></section>

<section class="block article"><div class="wrap">
<h2>The reads</h2>
<div class="tablewrap"><table class="data">
<thead><tr><th>DNA from</th><th>Sequencing</th><th>Runs</th><th class="num">Examined</th><th>Note</th></tr></thead>
<tbody>{rows_lib}</tbody></table></div>
<p class="small muted" style="margin-top:10px;max-width:90ch">Every read was compared with a panel of 1,271 algal plastid sequences (241 whole plastid genomes of green algae, plastid records of Dasycladales, Xanthophyceae and Cladophorales, and distant plastids to catch anything else), 117 algal mitochondrial genomes, 7,352 nuclear ribosomal DNA sequences of algae and the mitochondrial genomes of nine <i>Elysia</i> species. Reads with a match were then aligned base by base to the panel and to this species’ own assembly.</p>

<h2>How much of the DNA is not slug</h2>
{charts}

<h2>Which alga</h2>
{which}
{donor_tables}
<p class="small muted" style="margin-top:10px;max-width:90ch">A read is assigned to the reference it matches best. Identity near 100% means the reference is the alga eaten, or one of its own species; identity well below that means only that the nearest sequenced relative is named.</p>

{"<h2>Two barcode genes</h2>" + bc_tables if bc_tables else ""}

{"<h2>Did algal nuclei come along?</h2>" + nuc if nuc else ""}

{"<h2>For rhythm work</h2>" + rhythm if rhythm else ""}

{"<h2>What to keep in mind</h2><ul>" + notes + "</ul>" if notes else ""}
</div></section>
"""
    src = f"{it(sp['species'])}: reads {e(', '.join(r for L in libs for r in L['runs']))}; assembly {e(sp['asm_name'])} ({e(sp['accession'])}). Methods and sources are on the <a href=\"../methods.html\">methods page</a>."
    return chrome.page(f"Whose chloroplasts? {sp['species']}", body, root="../", species=sp["species"], tabs=tabs_for(slug), tab_current="kleptoplasts.html", today=today, sources=src, scripts=TIP_JS,
                       description=f"Which alga's chloroplasts the sequenced {sp['species']} carried, worked out from the plastid DNA among its genome reads.")

# =====================================================================================
RNA_GENE = re.compile(r"ribosomal RNA|^rr[nls]|^rn[ls]$|tRNA|^trn", re.I)
def locus_reading(L, tx):
    key = f'{L["contig"]}:{L["start"]}'
    if key in tx.get("loci_readings", {}): return tx["loci_readings"][key]
    g = L.get("genes") or []
    if any(RNA_GENE.search(x) for x in g): return "A ribosomal or transfer RNA gene. These are conserved across bacteria, mitochondria and plastids, so a match at this identity does not point to an alga."
    if L["cls"] == "algal mitochondrial" and L["identity"] < 0.9: return "Matches an algal mitochondrial genome only loosely; most likely a ribosomal RNA gene of a bacterium or of the slug’s own mitochondrial ancestry."
    if L["identity"] < 0.9: return "A loose match to a plastid gene. See the text."
    return "Plastid DNA of the kind this species carries. See the text."

def transfer_page(sp, K, T):
    slug = sp["slug"]; tx = TEXT.get(slug, {}); ctx = context(sp, K, T); s = T["screen"]
    lede = fmt(tx.get("transfer_lede", "Two places to look: the assembly, and the reads it was built from."), ctx)
    answer = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get("transfer_answer", []))
    # ---- assembly
    whole = ""
    if s["whole"]:
        rows = "".join(f'<tr><td class="mono">{e(w["contig"])}</td><td class="num">{num(w["length"])}</td><td>{e(w["cls"])}</td><td>{taxon_html(" ".join(w["best_target_desc"].split()[:2]))} <span class="mono muted">{e(w["best_target"].split("|", 1)[1])}</span></td><td class="num">{100 * w["median_identity"]:.1f}%</td></tr>' for w in s["whole"][:40])
        whole = f'''<h3>Sequences that are organelle DNA from end to end</h3>
<div class="tablewrap"><table class="data"><thead><tr><th>Sequence</th><th class="num">Length, bp</th><th>What it is</th><th>Best match</th><th class="num">Identity</th></tr></thead><tbody>{rows}</tbody></table></div>'''
    loci = ""
    if s["loci"]:
        rows = "".join(f'<tr><td class="mono">{e(L["contig"])}:{num(L["start"])}–{num(L["end"])}</td><td class="num">{num(L["span"])}</td><td class="num">{100 * L["identity"]:.1f}%</td>'
                       f'<td>{taxon_html(" ".join(L["target_desc"].replace("UNVERIFIED: ", "").split()[:2]))} {e(L["cls"].replace("algal ", "").replace("plastid", "plastid genome"))}{": " + e(", ".join(g for g in L.get("genes", []) if len(g) < 14)[:60]) if L.get("genes") else ""}</td>'
                       f'<td>{locus_reading(L, tx)}</td></tr>' for L in sorted(s["loci"], key=lambda L: -L["identity"]))
        loci = f'''<h3>Organelle-like stretches inside other sequences</h3>
<div class="tablewrap"><table class="data"><thead><tr><th>Place in the assembly</th><th class="num">Length, bp</th><th class="num">Identity</th><th>Matches</th><th>Reading</th></tr></thead><tbody>{rows}</tbody></table></div>'''
    numt = ""
    if s["numt_like"]:
        rows = "".join(f'<tr><td class="mono">{e(N["contig"])}:{num(N["start"])}–{num(N["end"])}</td><td class="num">{num(N["end"] - N["start"])}</td><td class="num">{100 * N["identity"]:.1f}%</td><td class="mono">{e(N["target"].split("|", 1)[1])}:{num(N["tstart"])}–{num(N["tend"])}</td></tr>' for N in s["numt_like"][:30])
        numt = f'''<h3>Mitochondrial DNA in nuclear sequences</h3>
<p>{fmt(tx.get("numt_text", "Stretches that match the slug’s own mitochondrial genome inside nuclear sequences. Copies of mitochondrial DNA in the nucleus are common in animals; they are listed because the same screen that finds them would find plastid DNA."), ctx)}</p>
<div class="tablewrap"><table class="data"><thead><tr><th>Place in the assembly</th><th class="num">Length, bp</th><th class="num">Identity</th><th>Part of the mitochondrial genome</th></tr></thead><tbody>{rows}</tbody></table></div>'''
    asm_text = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get("transfer_assembly", []))
    removed = s["removed"]
    asm_method = (f'The assembly was cut into 10 kb windows and each window compared with the same panel as the reads. '
                  f'{num(removed["repeat_like_matches"])} matches were set aside because they fall on stretches of the panel that are hit from 15 or more places in the assembly: slug repeats that happen to resemble a short piece of a panel sequence. '
                  f'{num(removed["low_complexity_matches"])} more were to low-complexity sequence, and {num(removed["own_rdna_loci"])} were the slug’s own ribosomal DNA, which resembles algal ribosomal DNA at 80 to 90%.')
    # ---- reads
    reads_html = ""
    for R in T["reads"]:
        ch = R["chimera"]; rates = ch["rates"]; long_ = R["kind"] != "short"
        unit = "reads" if long_ else "pairs"
        order = [k for k in ("cp", "slugmt") if k in rates]
        rows = ""
        for k in order:
            v = rates[k]; tot = v.get("reads_with_segment", v.get("pairs_with_mate", 0)); c_ = v.get("reads_also_with_slug", v.get("pairs_also_with_slug", 0))
            if k == "bait" or not tot: continue
            rows += f'<tr><td>{CLS_WORD[k].capitalize()} DNA</td><td class="num">{num(tot)}</td><td class="num">{num(c_)}</td><td class="num">{pct_auto(c_ / tot) if tot else "–"}</td></tr>'
        rec = ch.get("recurrent", [])
        rec_cp = [x for x in rec if (x.get("pair", x.get("kind", "")).find("cp") >= 0)]
        rec_mt = [x for x in rec if (x.get("pair", x.get("kind", "")).find("slugmt") >= 0)]
        key = "reads_" + R["library"]
        para = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get(key, []))
        reads_html += f'''<h3>{e(R.get("label") or (R["animal"][0].upper() + R["animal"][1:]))}: {gb(R.get("bases_screened") or R["bases"])} of {e(KIND[R["kind"]])}</h3>
{f'<p class="small muted" style="margin:-4px 0 10px">{e(R["animal"][0].upper() + R["animal"][1:])}.</p>' if R.get("label") else ""}
<div class="tablewrap"><table class="data"><thead><tr><th>Kind of DNA</th><th class="num">{unit.capitalize()} that carry it</th><th class="num">Of these, also carrying slug DNA</th><th class="num">Rate</th></tr></thead><tbody>{rows}</tbody></table></div>
<p class="small muted" style="margin-top:8px;max-width:90ch">{"A read counts as carrying both when it has two neighbouring stretches, each of 300 bases or more, one matching slug DNA and one matching organelle DNA. Joins seen in two or more independent molecules" if long_ else "A pair counts as carrying both when one mate matches organelle DNA and the other matches slug DNA at one clear place. Places in the assembly that collect three or more such pairs"}: {len(rec_cp)} for plastid DNA, {len(rec_mt)} for mitochondrial DNA.</p>
{para}'''
    reads_text = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get("transfer_reads", []))
    # places in the assembly tested with reads
    tests = [(R, R["locus_tests"]) for R in T["reads"] if R.get("locus_tests")]
    locus_html = ""
    if tests:
        names = [L["name"] for L in tests[0][1]["loci"]]
        rows = ""
        for i, nm in enumerate(names):
            first = True
            for R, lt in tests:
                L = lt["loci"][i]
                def edge(k):
                    v = L["edges"][k]
                    if "note" in v: return '<td class="num muted" colspan="2">end of the sequence</td>'
                    return f'<td class="num"><strong>{num(v["through"])}</strong></td><td class="num">{num(v["stop_inside"] + v["stop_outside"])}</td>'
                span = len(tests)
                head = (f'<td rowspan="{span}"><span class="mono">{e(L["contig"])}:{num(L["start"])}–{num(L["end"])}</span><br><span class="muted small">{"plastid-like" if L["kind"] == "plastid" else "mitochondrial"}, {num(L["end"] - L["start"])} bp</span></td>') if first else ""
                rows += f'<tr>{head}<td>{e(R.get("label") or (R["animal"].split(" (")[0].capitalize() + ", " + KIND[R["kind"]]))}</td>{edge("left")}{edge("right")}</tr>'
                first = False
        locus_html = f'''<h3>Each place in the assembly, tested with reads</h3>
<p>{fmt(tx.get("locus_text", "For every organelle-like stretch inside an assembly sequence, reads were aligned to the stretch and 3 kb either side. At each edge, a read that runs through joins the stretch to its flank; a read that stops there continues into something else."), ctx)}</p>
<div class="tablewrap"><table class="data"><thead><tr><th rowspan="2">Place</th><th rowspan="2">Reads from</th><th colspan="2" class="num">Left edge</th><th colspan="2" class="num">Right edge</th></tr>
<tr><th class="num">run through</th><th class="num">stop</th><th class="num">run through</th><th class="num">stop</th></tr></thead><tbody>{rows}</tbody></table></div>'''
    reads_html += locus_html + "".join(f'<p{" style=\"margin-top:16px\"" if i == 0 else ""}>{fmt(p, ctx)}</p>' for i, p in enumerate(tx.get("locus_after", [])))
    # ---- proteins
    prot = ""
    P = T.get("protein_screen")
    if P:
        ps = P["summary"]
        def census_of(c):
            k = c.get("census")
            if not k: return None
            st = k["stretches"]; n = sum(v for kk, v in st.items() if kk != "own_only")
            return dict(n=n, microbe=st.get("microbe", 0), alga=st.get("alga", 0), animal=st.get("animal", 0) + st.get("slug_only", 0),
                        coding=k["coding_bp"] / max(1, k["examined_length"] - k["n_bases"]), gc=k["gc"])
        def creading(c):
            if c["contig"] in tx.get("protein_readings", {}): return tx["protein_readings"][c["contig"]]
            if c.get("whole") == "plastid" or "plastid-only contig" in c["classes"]: return "Plastid DNA from end to end (listed above)."
            if "plastid fragment (assembly artefact)" in c["classes"]: return "The plastid stretch described above."
            k = census_of(c)
            if k and k["microbe"] >= 3 and k["microbe"] >= 0.5 * k["n"]: return "Bacterial sequence."
            if k: return "Too few genes to say."
            if c["loci"] >= 2 and not c["animal_loci_on_contig"]: return "Bacterial-type genes side by side and no slug gene among the proteins searched: a bacterial sequence in the assembly."
            if c["loci"] >= 2: return "Several bacterial-type genes: a bacterial sequence, or a sequence joined to one."
            return "One gene of bacterial type."
        def ccell(c):
            k = census_of(c)
            if not k: return '<td class="num muted">–</td><td class="num muted">–</td><td class="num muted">–</td><td class="num muted">–</td>'
            return f'<td class="num">{100 * k["gc"]:.0f}%</td><td class="num">{k["n"]}</td><td class="num">{k["microbe"]}</td><td class="num">{round(100 * k["coding"])}%</td>'
        cs = P.get("non_animal_contigs", [])
        has_census = any("census" in c for c in cs)
        rows = "".join(f'<tr><td class="mono">{e(c["contig"])}</td><td class="num">{num(c["length"])}</td>{ccell(c) if has_census else ""}<td>{e(", ".join(c["genes"][:7]))}{" …" if len(c["genes"]) > 7 else ""}</td>'
                       f'<td>{creading(c)}</td></tr>' for c in cs[:40])
        ptext = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get("transfer_proteins", []))
        census_head = ('<th class="num">GC</th><th class="num" title="Protein-coding stretches found when the whole sequence is translated in all six frames">Coding stretches</th>'
                       '<th class="num">of bacterial type</th><th class="num" title="Share of the sequence taken up by those coding stretches">Coding share</th>') if has_census else ""
        own = {"elysia-chlorotica": "<i>E. chlorotica</i>", "elysia-marginata": "<i>E. marginata</i>"}.get(slug)
        own_note = (f" The published gene annotation of this assembly is not counted as evidence that a place is animal: it lists genes on bacterial and plastid sequences too, and such a gene would otherwise vouch for itself." if own else "")
        nm = P["no_match"]; n_nm = ps.get("no database match", 0)
        nomatch_text = ("short or weak alignments whose translation matches nothing" if not nm["low_complexity"] else
                        f'{num(nm["low_complexity"])} are low-complexity alignments' + ("" if nm["low_complexity"] >= n_nm else " and the rest are short or weak ones (25 to 40% identity, stitched together across several exons) whose translation matches nothing"))
        census_note = (f'<p>Each sequence that carries a place of another kind was then translated as a whole, in all six frames, and every protein-coding stretch found on it was compared with the same proteins. The table gives their number, how many are of bacterial type, and the share of the sequence they take up. Bacterial DNA is mostly coding; a slug sequence of this size would hold a gene or two. GC is the share of G and C among the bases; the genome as a whole has {100 * P["genome_gc"]:.0f}%.</p>' if has_census and P.get("genome_gc") else "")
        prot = f'''<h2>Algal genes, searched as proteins</h2>
<p>DNA comparison can miss a gene from an alga that has no close relative among the references. So {num(P["n_queries"])} proteins were aligned straight to the assembly: the plastid-encoded proteins of a green alga, photosynthesis proteins of <i>Chlamydomonas</i> and <i>Arabidopsis</i> (including those encoded in the nucleus), proteins of siphonous green algae, and every available protein of <i>Vaucheria litorea</i>. Each place that matched was translated and compared with animal, algal and bacterial proteins.</p>
<p>{num(ps["loci"])} places matched. Of these, {num(ps.get("animal gene", 0))} are the slug’s own genes: a protein of another animal is the best match, as when mitochondrial ATP synthase is found with the chloroplast one.{own_note} Another {num(ps.get("no database match", 0))} match no protein of any other organism: {nomatch_text}.</p>
{census_note}
{ptext}
{f'<div class="tablewrap"><table class="data"><thead><tr><th>Sequence</th><th class="num">Length, bp</th>{census_head}<th>Proteins of the search that matched</th><th>Reading</th></tr></thead><tbody>{rows}</tbody></table></div>' if rows else ""}'''
    caveats = "".join(f"<li>{fmt(p, ctx)}</li>" for p in tx.get("transfer_caveats", []))
    body = f"""
<section class="page-head"><div class="wrap">
<h1>Has algal DNA entered this slug’s genome?</h1>
<p class="lede">{lede}</p>
<div class="answer">{answer}</div>
</div></section>

<section class="block article"><div class="wrap">
<h2>In the assembly</h2>
<p>{asm_method}</p>
{asm_text}
{whole}
{loci}
{numt}

{"<h2>In the reads</h2>" + reads_text + reads_html if T["reads"] else ""}

{prot}

{"<h2>What this cannot rule out</h2><ul>" + caveats + "</ul>" if caveats else ""}
</div></section>
"""
    src = f"{it(sp['species'])}: assembly {e(sp['asm_name'])} ({e(sp['accession'])}){'; reads ' + e(', '.join(r for R in T['reads'] for r in R['runs'])) if T['reads'] else ''}. Methods and sources are on the <a href=\"../methods.html\">methods page</a>."
    return chrome.page(f"Gene transfer? {sp['species']}", body, root="../", species=sp["species"], tabs=tabs_for(slug), tab_current="gene-transfer.html", today=today, sources=src, scripts=TIP_JS,
                       description=f"A test for algal DNA in the genome of the sea slug {sp['species']}: its assembly and its reads searched for plastid and algal sequences.")

os.makedirs(os.path.join(SITE, "downloads"), exist_ok=True)
for sp in SP:
    slug = sp["slug"]
    if slug == "elysia-crispata": continue
    K = jload(os.path.join(RG, slug + ".kleptoplast.json")); T = jload(os.path.join(RG, slug + ".transfer.json"))
    os.makedirs(os.path.join(SITE, slug), exist_ok=True)
    done = []
    for which, fn, builder, dl in (("klepto", "kleptoplasts.html", klepto_page, "_kleptoplast.json"), ("transfer", "gene-transfer.html", transfer_page, "_gene_transfer.json")):
        path = os.path.join(SITE, slug, fn)
        if has_page(slug, which):
            open(path, "w").write(builder(sp, K, T)); done.append(fn)
            shutil.copy(os.path.join(RG, slug + {"klepto": ".kleptoplast.json", "transfer": ".transfer.json"}[which]), os.path.join(SITE, "downloads", slug + dl))
        elif os.path.exists(path): os.remove(path)
    print(slug, done or "no pages yet")
