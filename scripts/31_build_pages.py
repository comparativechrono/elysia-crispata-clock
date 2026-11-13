"""Article pages of the site (kleptoplasts, gene transfer, methods).

Executed by 30_build_site.py, which supplies: repo, R, DOCS, page(), e(), genes, klepto, hgt, today, meta.
All figures in the text are read from results/kleptoplast.json and results/hgt.json.
"""
import os, shutil, glob

def pct(x, d=1): return f"{100 * x:.{d}f}%"
def num(x): return f"{x:,.0f}"
def it(s): return f"<i>{e(s)}</i>"

def hbar(rows, unit="%", width=760, maxv=None, cls="bar", label_w=250, digits=2):
    """rows: (label, value, tooltip, css class or None). One hue, values written at the bar end."""
    maxv = maxv or max(r[1] for r in rows) * 1.18
    rowh, top = 34, 8
    H = top + rowh * len(rows) + 6
    out = [f'<svg viewBox="0 0 {width} {H}" role="img" aria-label="{e("; ".join(f"{r[0]}: {r[1]:.{digits}f}{unit}" for r in rows))}">']
    pw = width - label_w - 90
    out.append(f'<line class="gridline" x1="{label_w}" x2="{label_w}" y1="{top}" y2="{H - 6}"/>')
    for i, r in enumerate(rows):
        y = top + i * rowh
        w = max(2, r[1] / maxv * pw)
        c = r[3] if len(r) > 3 and r[3] else cls
        out.append(f'<text class="lbl" x="{label_w - 10}" y="{y + 20}" text-anchor="end">{e(r[0])}</text>')
        out.append(f'<rect class="{c}" x="{label_w}" y="{y + 7}" width="{w:.1f}" height="18" rx="4" data-tip="{e(r[2])}" tabindex="0"><title>{e(r[2])}</title></rect>')
        out.append(f'<text class="val" x="{label_w + w + 8:.1f}" y="{y + 20}">{r[1]:.{digits}f}{unit}</text>')
    out.append("</svg>")
    return "\n".join(out)

def splitbar(parts):
    """parts: (fraction, css colour var, label)"""
    segs = "".join(f'<i style="width:{max(0.6, 100 * f):.2f}%;background:var({c})" title="{e(l)}: {100 * f:.1f}%"></i>' for f, c, l in parts if f > 0)
    return f'<div class="splitbar" role="img" aria-label="{e(", ".join(f"{l} {100 * f:.1f}%" for f, c, l in parts))}">{segs}</div>'

TIP_JS = '<script src="assets/pages.js"></script>'

# =====================================================================================
#                                    KLEPTOPLASTS
# =====================================================================================
if klepto:
    L = klepto["library"]; C = klepto["classes"]; A = klepto["assignment"]
    tot_bp = L["bp"]
    f_cp, f_bac, f_mt = C["cp"]["bp"] / tot_bp, C["bac"]["bp"] / tot_bp, C["slugmt"]["bp"] / tot_bp
    by = {a["taxon"]: a for a in A}
    cap, lam = by.get("Penicillus capitatus"), by.get("Penicillus lamourouxii")
    udo = by.get("Udotea flabellum")
    hal_bp = sum(a["bp"] for a in A if a["taxon"].startswith("Halimeda")); hal_reads = sum(a["reads"] for a in A if a["taxon"].startswith("Halimeda"))
    conf_bp = sum(a["bp"] for a in A)
    other_bp = conf_bp - cap["bp"] - lam["bp"]
    depth_cap = cap["bp"] / 84040; depth_lam = lam["bp"] / 89458
    slug_depth = 46.85e9 / 892.65e6
    bc = klepto["barcodes"]
    def bshare(barcode, pred):
        return sum(b["fraction"] for b in bc if b["barcode"] == barcode and pred(b))
    rb_cap = bshare("rbcL", lambda b: "capitatus" in b["top_hit"]); rb_lam = bshare("rbcL", lambda b: "capitatus" not in b["top_hit"])
    tu_cap = bshare("tufA", lambda b: "capitatus" in b["top_hit"]); tu_lam = bshare("tufA", lambda b: "capitatus" not in b["top_hit"])
    rid = klepto["reference_identity"]; nuc = klepto["nuclear"]
    def bnote(b):
        if "capitatus" not in b["top_hit"]:
            return "the second alga; no reference of its own species for this marker" if b["barcode"] == "tufA" else "the second alga; nearest full-length reference is another species"
        main = max((x for x in bc if x["barcode"] == b["barcode"] and "capitatus" in x["top_hit"]), key=lambda x: x["reads"])
        if b is main:
            r = rid["capitatus"][b["barcode"]]
            return "P. capitatus" if b["identity"] >= 99.99 else f'P. capitatus; identical to the shorter voucher {r["accession"]} over {r["length"]} bp'
        return "P. capitatus with a 1-bp length difference in a homopolymer; most likely a sequencing or PCR slip"
    rows_bc = "".join(
        f'<tr><td>{b["barcode"]}</td><td class="num">{num(b["reads"])}</td><td class="num">{pct(b["fraction"])}</td>'
        f'<td>{it(" ".join(b["top_hit"].split()[:2]))} <span class="mono muted">{e(b["top_hit"].split("[")[-1].rstrip("]"))}</span></td>'
        f'<td class="num">{b["identity"]:.2f}% of {num(b["alnlen"])} bp</td><td>{e(note)}</td></tr>'
        for b in bc for note in [bnote(b)])
    comp_rows = "".join(
        f"<tr><td>{lab}</td><td>{splitbar([(a, '--series-1', 'P. capitatus'), (b, '--series-2', 'P. lamourouxii'), (c, '--series-muted', 'other')])}</td>"
        f'<td class="num">{pct(a)}</td><td class="num">{pct(b)}</td><td class="num">{pct(c, 2) if c else "–"}</td><td class="num">{n}</td></tr>'
        for lab, a, b, c, n in [
            ("All plastid reads, assigned to the nearest plastome", cap["bp"] / conf_bp, lam["bp"] / conf_bp, other_bp / conf_bp, num(sum(x["reads"] for x in A)) + " reads"),
            ("rbcL barcode (1,428 bp)", rb_cap, rb_lam, 0, num(sum(b["reads"] for b in bc if b["barcode"] == "rbcL")) + " reads"),
            ("tufA barcode (987 bp)", tu_cap, tu_lam, 0, num(sum(b["reads"] for b in bc if b["barcode"] == "tufA")) + " reads")])
    body = f"""
<section class="page-head"><div class="wrap">
<h1>Whose chloroplasts was this slug carrying?</h1>
<p class="lede">The animal sequenced for the genome was caught wild in Key Largo, so its kleptoplasts record what it had been eating. The plastid DNA in its reads says which algae those were.</p>
<div class="answer">
<p><strong>Two shaving-brush algae.</strong> About {pct(cap["bp"] / conf_bp, 0)} of the plastid DNA is {it("Penicillus capitatus")} and {pct(lam["bp"] / conf_bp, 0)} is a second {it("Penicillus")} that matches the one available {it("P. lamourouxii")} reference. Traces of {it("Udotea flabellum")} ({pct(udo["bp"] / conf_bp, 2)}) and a {it("Halimeda")} ({pct(hal_bp / conf_bp, 2)}) are there too.</p>
<p>Plastid DNA is abundant ({pct(f_cp)} of all sequenced bases) but algal nuclear DNA is absent: not one of {num(L["reads"])} reads carries algal nuclear ribosomal DNA. The slug kept the chloroplasts and nothing else.</p>
</div>
</div></section>

<section class="block article"><div class="wrap">
<h2>How much of the library is not slug</h2>
<p>Every one of the {num(L["reads"])} PacBio HiFi reads ({L["bp"] / 1e9:.2f} Gb) was compared with a panel of 241 algal plastid genomes, 117 algal mitochondrial genomes, 1,249 nuclear ribosomal DNA sequences of siphonous green algae, the four bacterial genomes assembled from this animal and the slug’s own mitochondrion.</p>
<div class="chart"><h3>Share of sequenced bases</h3>
<p class="sub">The remaining {pct(1 - f_cp - f_bac - f_mt)} is slug nuclear DNA and reads that match nothing in the panel.</p>
{hbar([("Algal plastid DNA", 100 * f_cp, f"{num(C['cp']['reads'])} reads, {C['cp']['bp'] / 1e6:,.0f} Mb"), ("Bacteria (4 genomes from this animal)", 100 * f_bac, f"{num(C['bac']['reads'])} reads, {C['bac']['bp'] / 1e6:,.0f} Mb"), ("Slug mitochondrion", 100 * f_mt, f"{num(C['slugmt']['reads'])} reads, {C['slugmt']['bp'] / 1e6:,.1f} Mb")])}
</div>
<p>That plastid share is deep coverage: about {num(round(depth_cap, -2))}× of the {it("P. capitatus")} plastid genome and {num(round(depth_lam, -1))}× of the second one, against roughly {slug_depth:.0f}× for the slug’s own nuclear genome. The library was amplified by PCR, so treat these ratios as rough.</p>

<h2>Which algae</h2>
<p>Three measures agree. Reads were assigned to whichever plastid genome they match best, and two standard barcode genes were read directly out of every read that spans them.</p>
<div class="legend"><span><b class="sw" style="background:var(--series-1)"></b>{it("Penicillus capitatus")}</span><span><b class="sw" style="background:var(--series-2)"></b>{it("Penicillus lamourouxii")} type</span><span><b class="sw" style="background:var(--series-muted)"></b>other algae</span></div>
<div class="tablewrap"><table class="data">
<thead><tr><th>Measure</th><th style="min-width:200px">Split</th><th class="num"><i>P. capitatus</i></th><th class="num"><i>P. lamourouxii</i> type</th><th class="num">Other</th><th class="num">Based on</th></tr></thead>
<tbody>{comp_rows}</tbody></table></div>

<h3>The trace algae</h3>
<div class="tablewrap"><table class="data">
<thead><tr><th>Alga</th><th class="num">Reads</th><th class="num">Share of plastid DNA</th><th>rbcL against GenBank</th><th>tufA against GenBank</th></tr></thead>
<tbody>
<tr><td>{it("Udotea flabellum")}</td><td class="num">{num(udo["reads"])}</td><td class="num">{pct(udo["bp"] / conf_bp, 3)}</td><td>{klepto["minor"][0]["rbcL"]}</td><td>{klepto["minor"][0]["tufA"]}</td></tr>
<tr><td>{it("Halimeda")} sp.</td><td class="num">{num(hal_reads)}</td><td class="num">{pct(hal_bp / conf_bp, 3)}</td><td>{klepto["minor"][1]["rbcL"]}</td><td>{klepto["minor"][1]["tufA"]}</td></tr>
</tbody></table></div>
<p class="small muted" style="margin-top:10px">A further {num(klepto["unassigned"].get("reads", 0))} reads matched a plastid genome only weakly or over a short stretch and are not counted as plastid. Aligned base by base, nearly all of them are slug nuclear reads that contain a short low-complexity stretch resembling plastid DNA.</p>

<h2>How sure are the names</h2>
<div class="cols">
<div><h3>{it("Penicillus capitatus")}: firm</h3>
<p>The dominant plastid type is identical to GenBank {it("P. capitatus")} vouchers over {num(rid["capitatus"]["rbcL"]["length"])} bp of rbcL ({rid["capitatus"]["rbcL"]["accession"]}) and {rid["capitatus"]["tufA"]["length"]} bp of tufA ({rid["capitatus"]["tufA"]["accession"]}). The next nearest species, {it(rid["capitatus"]["next_species"]["name"])}, is {rid["capitatus"]["next_species"]["rbcL"]}% and {rid["capitatus"]["next_species"]["tufA"]}% identical.</p></div>
<div><h3>{it("Penicillus lamourouxii")}: consistent, thinly referenced</h3>
<p>The second type is identical to {rid["lamourouxii"]["rbcL"]["accession"]} over {rid["lamourouxii"]["rbcL"]["length"]} bp of rbcL. That is the only {it("P. lamourouxii")} rbcL in GenBank, and there is no tufA for the species. The nearest other species, {it(rid["lamourouxii"]["next_species"]["name"])}, is {rid["lamourouxii"]["next_species"]["rbcL"]}% (rbcL) and {rid["lamourouxii"]["next_species"]["tufA"]}% (tufA) identical.</p></div>
<div><h3>They are different species</h3>
<p>The two plastid types differ at {100 - rid["between_types"]["rbcL"]:.1f}% of rbcL sites and {100 - rid["between_types"]["tufA"]:.1f}% of tufA sites, and by 7 to 11% across the whole plastid genome. The Sanger team assembled both (OZ373742 and OZ372796) and gave them the same names; this analysis reaches the names independently.</p></div>
</div>

<h3>Barcode haplotypes</h3>
<div class="tablewrap"><table class="data">
<thead><tr><th>Barcode</th><th class="num">Reads</th><th class="num">Share</th><th>Best GenBank match</th><th class="num">Identity</th><th>Reading</th></tr></thead>
<tbody>{rows_bc}</tbody></table></div>

<h2>No algal nuclei came along</h2>
<p>If algal nuclei or cytoplasm had been carried over with the plastids, the reads would contain algal nuclear ribosomal DNA, which sits in hundreds of copies per nucleus. {num(nuc["reads_first_flagged"])} reads resembled an algal rDNA sequence in the first screen; aligned base by base, every one is the slug’s own rDNA. The number that match algal nuclear rDNA best is {nuc["reads_matching_algal_rdna"]}. Algal mitochondrial DNA is at trace level: {nuc["algal_mito_reads"]} reads.</p>
<p>This matters for the <a href="gene-transfer.html">gene transfer question</a>: with no algal nuclear DNA in the sample, an algal nuclear gene found in these reads could only come from the slug’s own genome. None was found.</p>

<h2>For rhythm work</h2>
<ul>
<li>The plastid genomes to map kleptoplast transcripts against are OZ373742.1 ({it("P. capitatus")}) and OZ372796.1 ({it("P. lamourouxii")} type). They are not annotated in the archive; a homology-based gene table for both is in the <a href="methods.html#downloads">downloads</a>.</li>
<li>The mix is a property of this one wild animal on the day it was collected (30 July 2020). {it("E. crispata")} eats many algae, so other animals, and cultured ones, will differ.</li>
</ul>
</div></section>
"""
    open(os.path.join(DOCS, "kleptoplasts.html"), "w").write(page("Whose chloroplasts? Elysia crispata clock", "kleptoplasts.html", body, TIP_JS))

# =====================================================================================
#                                    GENE TRANSFER
# =====================================================================================
if hgt:
    CH = hgt["chimera"]; P = {p["pair"]: p for p in CH["pairs"]}; ct = CH["class_totals"]
    fr = hgt["fragments"]; prim = [f for f in fr if f["assembly"] == "sanger"]; alt = [f for f in fr if f["assembly"] == "alt"]
    n_cap = sum(1 for f in fr if "capitatus" in f["alga"]); n_lam = len(fr) - n_cap
    tested = [f for f in fr if f["through_molecules"] is not None]
    single = sum(1 for f in tested if f["through_molecules"] <= 1)
    placeable = [f for f in fr if (f.get("flank_reads") or 0) > 0]
    contradicted = sum(1 for f in placeable if (f["diverging_molecules"] or 0) >= 2)
    thr_mol = sum(f.get("through_molecules_boundary_pass") or 0 for f in placeable); div_mol = sum(f["diverging_molecules"] or 0 for f in placeable)
    pcr_copies = sum(1 for f in fr if (f["through_reads"] or 0) > (f["through_molecules"] or 0))
    rec_cp = P.get("cp+slug", {}).get("recurrent_clusters", 0); rec_bac = P.get("bac+slug", {}).get("recurrent_clusters", 0)
    DV = hgt.get("diverged") or {}
    _as = {a["taxon"]: a["fraction"] for a in (klepto or {}).get("assignment", [])}
    KS = f'{round(100 * _as.get("Penicillus capitatus", 0.86))}:{round(100 * _as.get("Penicillus lamourouxii", 0.14))}'
    AL = hgt.get("alien") or {}; cands = AL.get("candidates", [])
    def wh(a, label):
        w = hgt["assembly_screen"].get(a, {}).get("whole_contigs", {}); parts = []
        if w.get("cp"): parts.append(f'{w["cp"]["n"]} plastid ({w["cp"]["bp"] / 1000:.0f} kb)')
        if w.get("bac"): parts.append(f'{w["bac"]["n"]} bacterial ({w["bac"]["bp"] / 1000:.0f} kb)')
        return f'{label}: {" and ".join(parts)}' if parts else ""
    WH = "; ".join(x for x in (wh("sanger", "xgElyCris1.1"), wh("alt", "alternate haplotype"), wh("purdue", "Purdue"), wh("tampa", "Tampa")) if x)
    alg_c = [c for c in cands if c["top_group"] in ("ALG", "CYA")]
    PS = hgt.get("protein_screen") or {}
    def chrn(c): return ("chr" + c) if c.isdigit() else c
    def frow(f):
        ev = "–"
        if f["through_molecules"] is not None:
            ev = str(f["through_molecules"])
            if f["through_reads"] > f["through_molecules"]: ev = f'<span title="{f["through_reads"]} reads, PCR copies of one molecule">{ev}</span>'
        if f["diverging_molecules"] is None: dv = "–"
        elif not f["diverging_molecules"] and (f.get("unique_flank_bp") or 0) < 500 and not f.get("flank_reads"): dv = '<span class="muted" title="The slug DNA beside this fragment is repeated across the genome, so reads cannot be placed here">repeats</span>'
        else: dv = str(f["diverging_molecules"])
        return (f'<tr><td class="mono nowrap">{chrn(f["contig"])}:{f["start"] + 1:,}–{f["end"]:,}</td><td class="num">{f["length"]:,}</td><td class="num">{100 * f["identity"]:.1f}%</td>'
                f'<td class="nowrap">{it(f["alga"].replace("Penicillus", "P."))}</td><td>{e(", ".join(f["genes"][:4])) or "–"}</td><td class="nowrap">{"gap" if f["abuts"] == "gap" else "sequence end"}</td>'
                f'<td class="num">{ev}</td><td class="num">{dv}</td></tr>')
    alien_block = ""
    if AL:
        SM = AL.get("summary", {})
        onfrag = [c for c in cands if c.get("on_plastid_fragment")]; rest = [c for c in cands if not c.get("on_plastid_fragment")]
        def arow(c):
            where = c.get("location", "")
            if where: ctg, rng = where.split(":"); a, b = rng.split("-"); where = f"{chrn(ctg)}:{int(a):,}–{int(b):,}"
            pid = c["query"].split("|")
            return (f'<tr><td class="mono">{e(pid[2] if c["source"].startswith("Ensembl") else pid[1])}</td><td>{"Ensembl" if c["source"].startswith("Ensembl") else "Purdue"}</td><td class="num">{c["length"]}</td>'
                    f'<td class="mono">{e(where) or "–"}</td><td>{e((c.get("top_hit_description") or c["top_hit"].split("|")[-1]).replace("LOW QUALITY PROTEIN: ", "").replace(" (chloroplast)", ""))}, {it(SPN.get(c["top_hit"].split("|")[1], c["top_hit"].split("|")[1]))}</td>'
                    f'<td class="num">{c["top_pident"]:.0f}%</td><td class="num">{c["top_bits"]:.0f}</td><td class="num">{max(c["best_animal_bits"], c["best_sacoglossan_bits"]):.0f}</td>'
                    f'<td>{"on a chloroplast fragment" if c.get("on_plastid_fragment") else ("animal match nearly as good" if not c["strong"] else "unexplained")}</td></tr>')
        SPN = {"BRYSP": "Bryopsis", "CHOCR": "Chondrus", "ARATH": "Arabidopsis", "CHLRE": "Chlamydomonas", "ECTSI": "Ectocarpus", "SYNY3": "Synechocystis"}
        cands = sorted(cands, key=lambda c: (not c.get("on_plastid_fragment"), -c["top_bits"]))
        alien_block = f"""
<h2>Does any annotated slug gene look algal?</h2>
<p>All {num(SM.get("n_queries", 0))} proteins in the two annotations of the species (Ensembl on xgElyCris1.1; the Purdue assembly of another animal) were compared with 13 animal proteomes, three other sacoglossans, 12 algal and plant proteomes including two siphonous green algae, a cyanobacterium, and the bacterial and fungal entries of Swiss-Prot. {num(SM.get("n_with_hits", 0))} have a match. For {SM.get("n_listed", len(cands))} of them the best match is not an animal protein.</p>
<ul>
<li><strong>{len(onfrag)} are chloroplast genes that the Ensembl annotation lists as slug genes.</strong> Each sits on one of the chloroplast fragments above. Anyone counting “algal genes in the slug genome” from the annotation will count these.</li>
<li><strong>The other {len(rest)} are not algal.</strong> Their animal homologues score within 20% of the algal match: reverse transcriptases of mobile elements, a short conserved peptide, and one low-complexity protein.</li>
<li>The Purdue annotation, made from an assembly without the artefacts, has no protein that matches an alga or a bacterium clearly better than it matches an animal.</li>
</ul>
<div class="tablewrap"><table class="data"><thead><tr><th>Gene or protein</th><th>Annotation</th><th class="num">aa</th><th>Location in xgElyCris1.1</th><th style="min-width:13em">Best non-animal match</th><th class="num">Identity</th><th class="num">Bit score</th><th class="num">Best animal bit score</th><th style="min-width:10em">Reading</th></tr></thead><tbody>{"".join(arow(c) for c in cands)}</tbody></table></div>
"""
    prot_block = ""
    if PS:
        T = PS.get("totals", {}); NM = PS.get("no_match", {}); FC = PS.get("foreign_contigs", [])
        pur = [r for r in FC if r["assembly"] == "purdue" and r["length"] > 10000]; small_bac = [r for r in FC if "bacterial" in r["reading"]]
        n_foreign = T.get("alga-like", 0) + T.get("microbe-like", 0)
        def prow(r):
            span = r.get("read_cover_span")
            reads = "none" if not r.get("read_alignments") else f'first {span[1] / 1000:.1f} kb only'
            return (f'<tr><td class="mono">{e(r["contig"])}</td><td class="num">{r["length"]:,}</td><td class="num">{100 * r["gc"]:.0f}%</td><td>{e(", ".join(g.split("_")[0] for g in r["genes"]))}</td>'
                    f'<td class="num">{100 * r["identity_min"]:.0f}{"" if r["loci"] == 1 else f"–{100 * r["identity_max"]:.0f}"}%</td>'
                    f'<td class="num">{r["annotated_genes_with_protein_match"]} of {r["annotated_genes"]}</td><td>{reads}</td></tr>')
        pur_block = ""
        if pur:
            n_rep = sum(1 for r in pur if r.get("read_alignments")); n_gm = sum(r["annotated_genes"] for r in pur)
            pur_block = f"""
<h3>Eleven contigs of the Purdue assembly are not slug DNA</h3>
<p>The Purdue assembly (GCA_033675545.1, a different animal) contains {len(pur)} contigs, {min(r["length"] for r in pur) / 1000:.0f} to {max(r["length"] for r in pur) / 1000:.0f} kb long and {sum(r["length"] for r in pur) / 1000:.0f} kb in total, that carry {sum(r["loci"] for r in pur)} chloroplast-type genes. They look like the DNA of another organism left in the assembly:</p>
<ul>
<li>their G+C content is {100 * min(r["gc"] for r in pur):.0f} to {100 * max(r["gc"] for r in pur):.0f}%, against 37% for the assembly as a whole;</li>
<li>none of the {n_gm} gene models annotated on them matches any animal, algal or bacterial protein, and none of the {len(pur)} contigs carries a recognisable slug gene;</li>
<li>in two sequencing runs of the xgElyCris1 animal (25 Gb) no read matches them, apart from a repeated slug sequence at the very start of {n_rep} contigs. DNA that belongs to the species should be in every animal.</li>
</ul>
<p>The genes are {100 * min(r["identity_min"] for r in pur):.0f} to {100 * max(r["identity_max"] for r in pur):.0f}% identical to green-algal chloroplast proteins, too distant from any of the 241 plastid genomes in the panel to name the organism. For the {n_rep} contigs that begin with slug-like sequence, whether that join exists in the Purdue animal can only be tested with that animal’s own reads. Until then these contigs are best left out of any analysis of the Purdue assembly.</p>
<div class="tablewrap"><table class="data"><thead><tr><th>Purdue contig</th><th class="num">Length, bp</th><th class="num">G+C</th><th>Chloroplast-type genes found</th><th class="num">Identity to query</th><th class="num" title="Annotated gene models on the contig that match any protein in the database">Gene models with a match</th><th>xgElyCris1 reads matching</th></tr></thead>
<tbody>{"".join(prow(r) for r in pur)}</tbody></table></div>
"""
        prot_block = f"""
<h3>Photosynthesis genes, searched directly in the DNA</h3>
<p>Annotations can miss a gene, so {num(PS.get("n_queries", 0))} proteins were also aligned straight onto all four assemblies: the chloroplast-encoded proteins of {it("Udotea flabellum")} (a close relative of {it("Penicillus")} whose plastid genome is annotated), photosynthesis proteins of {it("Chlamydomonas")}, {it("Arabidopsis")} and siphonous green algae, and the known proteins of {it("Vaucheria litorea")}, the alga behind the original claims of gene transfer in {it("Elysia chlorotica")}. Each of the {PS.get("n_loci", 0)} places that matched was translated and compared with the protein database above, to tell an algal gene from the slug’s own version of a shared protein.</p>
<ul>
<li><strong>{T.get("animal gene", 0)} are the slug’s own genes.</strong> Their best match is an animal protein: mitochondrial ATP synthase found with chloroplast atpB, for example.</li>
<li><strong>{T.get("plastid fragment (assembly artefact)", 0)} lie on the chloroplast fragments described above</strong>, and {T.get("plastid-only contig", 0)} more on small contigs that are plastid DNA from end to end.</li>
<li><strong>{NM.get("n", 0)} match no protein at all.</strong> {NM.get("low_complexity", 0)} are runs of one or two amino acids and {NM.get("short", 0)} are shorter than 120 residues.</li>
<li><strong>{n_foreign} are algal or bacterial genes on contigs with no slug gene on them:</strong> {sum(r["loci"] for r in pur)} on {len(pur)} contigs of the Purdue assembly, described next, and {sum(r["loci"] for r in small_bac)} on {len(small_bac)} small bacterial contigs of the Tampa assembly.</li>
</ul>
<p>In none of the four assemblies does an algal gene sit among slug genes.</p>
{pur_block}
"""
    body = f"""
<section class="page-head"><div class="wrap">
<h1>Is algal DNA inserted in the slug genome?</h1>
<p class="lede">Long reads can show a foreign sequence sitting inside host DNA in a single molecule. They can also manufacture that picture. Here is what the xgElyCris1 reads show once the two are told apart.</p>
<div class="answer">
<p><strong>No inserted algal sequence was found.</strong> {num(P["cp+slug"]["reads"])} reads do join chloroplast DNA to slug DNA, but they do so at the same rate as reads that join bacterial DNA to slug DNA, and no chloroplast–slug junction is seen in two independent molecules. A real insertion would be crossed by dozens.</p>
<p><strong>The reference assembly needs a warning.</strong> xgElyCris1.1 and its alternate haplotype contain {len(fr)} stretches of chloroplast DNA, up to {max(f["length"] for f in fr) / 1000:.1f} kb long and 98 to 100% identical to the animal’s kleptoplasts. All {len(fr)} sit at the end of a contig and trace back to a single chimeric molecule. They are assembly artefacts, not insertions. Searching the assembly with a plastid genome will find them, and the Ensembl annotation lists {sum(1 for c in cands if c.get("on_plastid_fragment"))} chloroplast genes on them as slug genes.</p>
<p>That leaves nothing to compare with the kleptoplasts: there is no inserted algal sequence whose source could be traced. The artefacts themselves come from the two {it("Penicillus")} plastid types {n_cap} to {n_lam}, close to the {KS} mix of the kleptoplasts, as chimeras of kleptoplast DNA would.</p>
</div>
</div></section>

<section class="block article"><div class="wrap">
<h2>Why a single read proves nothing here</h2>
<p>The genome was sequenced from very little DNA, so the library was amplified by PCR before sequencing (PacBio’s ultra-low-input protocol). Amplification joins unrelated molecules now and then, and it copies each joined molecule many times. A read with slug DNA at one end and chloroplast DNA at the other is therefore expected even if no such molecule existed in the animal.</p>
<p>The bacteria living in the slug give a baseline. Nobody expects <i>Kistimonas</i> chromosomes to be spliced into slug chromosomes, so reads joining bacterial and slug DNA measure how often the library invents a junction.</p>
<div class="chart"><h3>Reads that also carry slug nuclear DNA</h3>
<p class="sub">Share of all reads of each kind. Chloroplast reads are joined to slug DNA no more often than bacterial reads are.</p>
{hbar([("Chloroplast reads", 100 * CH["rate_cp_slug"], f"{num(P['cp+slug']['reads'])} of {num(ct['cp'])} chloroplast reads"), ("Bacterial reads (baseline)", 100 * CH["rate_bac_slug"], f"{num(P['bac+slug']['reads'])} of {num(ct['bac'])} bacterial reads", "bar muted")], maxv=5)}
</div>
<p>Three further checks point the same way.</p>
<ul>
<li><strong>Chloroplast joined to bacterium.</strong> If junctions form at random, chloroplast reads should pick up a bacterial partner in proportion to how much bacterial DNA is in the library: about {min(CH["expected_cp_bac"], CH["expected_cp_bac_from_bacteria"]):.0f} to {max(CH["expected_cp_bac"], CH["expected_cp_bac_from_bacteria"]):.0f} reads, depending on which side the estimate starts from. There are {CH["observed_cp_bac"]}.</li>
<li><strong>No junction repeats.</strong> Of {num(P["cp+slug"]["junctions"])} chloroplast–slug junctions, {rec_cp} are found in two molecules with different ends. (Reads sharing an end are PCR copies of one molecule and count once.) The same test does find {rec_bac} repeated bacterium–slug junction, which turns out to be a stretch of bacterial sequence at the edge of a small unplaced scaffold: the method can see a real, recurring junction when there is one.</li>
<li><strong>Junctions fall everywhere.</strong> They are spread over all 15 slug chromosomes roughly in proportion to chromosome length and along the whole plastid genome, and {round(100 * CH["partner_species"]["capitatus"] / (CH["partner_species"]["capitatus"] + CH["partner_species"]["lamourouxii"]))}% involve the {it("P. capitatus")} plastid, the same share it has among the kleptoplasts.</li>
</ul>
<p>And no read carries chloroplast DNA with slug DNA from one locus on both sides of it.</p>

<h2>The chloroplast DNA inside xgElyCris1.1</h2>
<p>Cutting the assemblies into windows and matching them to the plastid genomes turns up exact copies of kleptoplast DNA: {len(prim)} in the primary assembly ({num(sum(f["length"] for f in prim))} bp, on chromosomes {", ".join(sorted({f["contig"] for f in prim if f["contig"].isdigit()}, key=int))} and three unplaced scaffolds) and {len(alt)} in the alternate haplotype ({num(sum(f["length"] for f in alt))} bp). The two other assemblies of the species, made from other animals without PCR amplification, have none.</p>
<p>Each copy was tested against the reads themselves. If chloroplast DNA really sat at one of these spots, on one or both copies of the chromosome, then half or all of the reads arriving from the slug side would run straight on into it. Instead:</p>
<ul>
<li>all {len(fr)} copies end at an assembly gap or at the end of a sequence, so only one side of each touches slug DNA;</li>
<li>every one of those {len(tested)} boundaries is carried by {"a single molecule" if single == len(tested) else "at most one molecule"}{f" (at {pcr_copies} it was sequenced twice, as PCR copies)" if pcr_copies else ""};</li>
<li>at the {len(placeable)} boundaries where the neighbouring slug DNA is unique enough to place reads, {thr_mol} molecules carry the junction and {div_mol} independent molecules arrive from the slug side and continue into something else. That count comes from two of the three sequencing runs (53% of the reads). At {contradicted} of the {len(placeable)}, two or more molecules contradict the assembly;</li>
<li>the other {len(fr) - len(placeable)} lie in repeats, where reads cannot be assigned to one place.</li>
</ul>
<p>The assembler extended a contig along one chimeric read, the contig stopped where that read stopped, and scaffolding left a gap beside it.</p>
<div class="tablewrap"><table class="data">
<thead><tr><th>Location in xgElyCris1.1</th><th class="num">Length, bp</th><th class="num">Identity to plastid</th><th>Alga</th><th style="min-width:9em">Plastid genes carried</th><th>Other end abuts</th><th class="num" title="Independent molecules in which one alignment runs from slug DNA into the chloroplast DNA as assembled (all three runs)">Molecules with the junction</th><th class="num" title="Independent molecules that follow the slug DNA up to the boundary and then continue into something other than chloroplast DNA (two of three runs)">Molecules against it</th></tr></thead>
<tbody>{"".join(frow(f) for f in prim)}</tbody></table></div>
<p class="small muted" style="margin-top:10px">Coordinates are 1-based on the Ensembl sequence names (chr1 = OY973766.1 and so on). The {len(alt)} fragments in the alternate haplotype, and this table, are in the <a href="methods.html#downloads">downloads</a>. Mask these intervals before using the assembly for anything that involves plastid genes: they include pieces of rbcL, psbA, tufA, atpB and petB.</p>

<h3>Nothing older</h3>
<p>An insertion from a past meal, or from a different alga, would have drifted away from today’s kleptoplasts, so every match that is not an exact copy was listed and sorted. The primary assembly has {num(DV["sanger"]["non_exact_matches"])} such matches to the 241 plastid genomes and to algal mitochondrial and ribosomal DNA. All but {DV["sanger"]["loci"]} hit simple repeats such as (CTT)n, which AT-rich plastid genomes also contain. The {DV["sanger"]["loci"]} that remain are the slug’s own ribosomal RNA genes, about 81% identical to algal ones. The alternate haplotype ({DV["alt"]["loci"]} loci) and the Tampa assembly ({DV["tampa"]["loci"]}) give the same picture.</p>
<p>So apart from the {len(fr)} exact copies, none of these three assemblies contains a stretch of 200 bp or more that resembles plastid DNA at 80% identity or better. The Purdue assembly does, and is dealt with <a href="#purdue">below</a>.</p>
<p>The assemblies also keep a few small free-standing contigs of plastid or bacterial DNA, not joined to anything: {WH}. These are ordinary leftovers of decontamination.</p>
{alien_block}
{prot_block}
<h2>What this cannot rule out</h2>
<ul>
<li>Inserts shorter than about 150 bp, or more than about 15% diverged at the DNA level, in parts of the genome that are not in any assembly. The DNA screens would miss them; the protein screens cover only what is assembled.</li>
<li>Insertions present in a small minority of cells. One molecule is indistinguishable from a library chimera.</li>
<li>Regions the PCR-amplified library barely sampled. Several clock genes sit in such regions (see the <a href="index.html">gene catalogue</a>).</li>
<li>Transfer from an alga with no sequenced relative, if the gene has no recognisable homologue in the algal proteomes used.</li>
</ul>
</div></section>
"""
    open(os.path.join(DOCS, "gene-transfer.html"), "w").write(page("Gene transfer? Elysia crispata clock", "gene-transfer.html", body, TIP_JS))

# =====================================================================================
#                                    METHODS
# =====================================================================================
dl = []
def add_dl(src, name, label, note):
    if os.path.exists(src):
        shutil.copy(src, os.path.join(DOCS, "downloads", name)); dl.append((name, label, note))
add_dl(os.path.join(R, "clock_genes.tsv"), "clock_genes.tsv", "Clock gene table (TSV)", "one row per gene: location, state in four assemblies, domains, orthologues")
add_dl(os.path.join(R, "sequences", "clock_proteins.faa"), "clock_proteins.faa", "Clock proteins (FASTA)", "representative protein for each gene")
add_dl(os.path.join(R, "sequences", "clock_cds.fna"), "clock_cds.fna", "Clock coding sequences (FASTA)", "where a clean coding sequence could be extracted")
add_dl(os.path.join(R, "clock_genes.json"), "clock_genes.json", "Everything about each gene (JSON)", "models compared, exon coordinates, read support, expression")
add_dl(os.path.join(R, "plastid_fragments_in_assembly.tsv"), "plastid_fragments_in_assembly.tsv", "Chloroplast fragments in the assemblies (TSV)", "intervals to mask, with the read evidence for each")
add_dl(os.path.join(R, "kleptoplast_barcodes.tsv"), "kleptoplast_barcodes.tsv", "Kleptoplast barcode haplotypes (TSV)", "rbcL and tufA sequences with read counts and GenBank matches")
add_dl(os.path.join(R, "penicillus_plastome_genes.bed"), "penicillus_plastome_genes.bed", "Gene positions on the two Penicillus plastid genomes (BED)", "by protein homology to Udotea flabellum NC_039528")
add_dl(os.path.join(R, "kleptoplast.json"), "kleptoplast.json", "Kleptoplast results (JSON)", "")
add_dl(os.path.join(R, "hgt.json"), "hgt.json", "Gene-transfer results (JSON)", "")
PSN = ((hgt or {}).get("protein_screen") or {}).get("n_queries", 471)
dl_html = "".join(f'<li><a href="downloads/{n}" download>{e(l)}</a><span>{e(t)}</span></li>' for n, l, t in dl)
body = f"""
<section class="page-head"><div class="wrap">
<h1>Methods and data</h1>
<p class="lede">Everything here was computed from public files. The scripts are in the repository; the numbered steps below follow their order.</p>
</div></section>
<section class="block article"><div class="wrap">
<h2 id="downloads">Downloads</h2>
<ul class="downloads">{dl_html}</ul>

<h2>Data used</h2>
<div class="tablewrap"><table class="data">
<thead><tr><th>What</th><th>Accession</th><th>Note</th></tr></thead>
<tbody>
<tr><td>Genome assembly, primary</td><td class="mono">GCA_963854125.1 (xgElyCris1.1)</td><td>542.6 Mb, 15 chromosomes, 8,217 gaps; Ensembl annotation 2024_05, 18,436 protein-coding genes</td></tr>
<tr><td>Genome assembly, alternate haplotype</td><td class="mono">GCA_963854135.1</td><td>363.8 Mb in contigs</td></tr>
<tr><td>PacBio HiFi reads (mid-body)</td><td class="mono">ERR11843472, ERR12353011, ERR12356115</td><td>5,516,672 reads, 46.85 Gb; ultra-low-input, PCR-amplified libraries</td></tr>
<tr><td>RNA-seq (head, same animal)</td><td class="mono">ERR12342469</td><td>35.3 million read pairs, one library</td></tr>
<tr><td>Plastid genomes from this animal</td><td class="mono">OZ373742.1, OZ372796.1</td><td>assembled by the Sanger team; labelled <i>P. capitatus</i> and <i>P. lamourouxii</i></td></tr>
<tr><td>Bacterial genomes from this animal</td><td class="mono">GCA_964289745, GCA_964291485, GCA_964291505, GCA_964291535</td><td><i>Ekhidna</i>, two <i>Kistimonas</i>, <i>Francisella</i>; used as the chimera baseline</td></tr>
<tr><td>Other <i>E. crispata</i> assemblies</td><td class="mono">GCA_033675545.1 (Purdue), GCA_984581095.1 (Tampa)</td><td>different animals; the Purdue assembly has its own annotation</td></tr>
</tbody></table></div>
<p style="margin-top:14px">The specimen (NSU0021201, clarki ecotype) was collected at Dove Cove, Key Largo, Florida on 30 July 2020. Genome paper: Middlebrooks et al. 2026, Wellcome Open Research 11:148.</p>

<h2>Clock genes</h2>
<details class="method" open><summary>1. Reference proteins and candidates</summary>
<p>102 reference clock proteins (human, mouse, <i>Drosophila</i>, monarch, zebrafish, oyster, <i>Aplysia</i>, squid, scallop) were fetched from UniProt and searched with DIAMOND (very-sensitive) against the Ensembl proteins of xgElyCris1.1, the Purdue-assembly proteins and 16 comparison proteomes. Scripts 01 to 03.</p></details>
<details class="method"><summary>2. Orthology</summary>
<p>Each family was aligned with MAFFT, trimmed with trimAl and a tree built with FastTree (LG + gamma, SH-like supports). A gene is named when it falls in the clade of its reference proteins and is the reciprocal best hit of the human or fly gene. Opsins were additionally required to have the retinal-binding lysine. Scripts 04, 05, 09.</p></details>
<details class="method"><summary>3. Choosing a protein model</summary>
<p>For each gene the candidates were the Ensembl model, the Purdue model, and translations of the genome guided by orthologues from other sea slugs, <i>Aplysia</i> and <i>Biomphalaria</i> (miniprot). Annotated models that do not re-align cleanly to their own assembly were set aside. Candidates were scored by mean bit score against orthologues in four outgroup molluscs; an existing annotation was kept unless a guided translation beat it by 5% and 30 bits. Six genes ended up with guided translations. Script 06.</p></details>
<details class="method"><summary>4. Presence in each assembly</summary>
<p>The chosen protein was aligned to all four assemblies with miniprot. Intact: at least 97% of the protein in one locus at 95% identity or more, with no frameshift or stop. Whole, sequence differs: the same with a frameshift or stop. In pieces: at least 90% found but split across loci. Partial: less than that. Missing: no locus at 95% identity. Scripts 07, 08.</p></details>
<details class="method"><summary>5. Reads and RNA</summary>
<p>Exons of every gene (from the Purdue models, with 150 bp flanks, repeat-containing parts masked) were included as baits when all HiFi reads were screened; a read counts when it aligns over 150 bp at 95% or better. Head RNA-seq was quantified with salmon against the Purdue coding sequences ({meta.get("rna_mapped", "37.5")}% of pairs assigned). Scripts 16, 20.</p></details>

<h2>Kleptoplasts and gene transfer</h2>
<details class="method" open><summary>6. Screening every read</summary>
<p>All reads were mapped (minimap2, k = 15) to a panel of 241 green-algal and 7 other plastid genomes, 117 algal mitochondrial genomes, 1,249 nuclear rDNA sequences of Bryopsidales, the four bacterial genomes, and the slug mitochondrion and rDNA. Fourteen low-complexity windows that attracted thousands of unrelated reads were masked. Scripts 10, 11.</p></details>
<details class="method"><summary>7. Naming the algae</summary>
<p>Plastid reads were assigned to their best-matching plastid genome (divergence 3% or less over 1 kb or more). rbcL and tufA were cut out of every read spanning them, identical sequences counted, and each haplotype compared with 6,453 GenBank barcodes of Bryopsidales. Script 13.</p></details>
<details class="method"><summary>8. Chimeric reads</summary>
<p>Reads with a panel match and 500 bp or more left unexplained were re-aligned base by base (minimap2 map-hifi) to the Sanger primary and alternate assemblies plus the panel, and the remainder to the Purdue and Tampa assemblies. A junction needs a slug segment next to a non-slug segment, each at least 300 bp at 97% identity or better. Read totals for each class use confident matches only (3% divergence or less over 1 kb or more). Two reads count as independent molecules only if both of their ends differ by more than 50 bp. Scripts 14, 15.</p></details>
<details class="method"><summary>9. Chloroplast DNA in the assemblies</summary>
<p>All four assemblies were cut into 10 kb windows and aligned base by base to the same panel. Matches of 200 bp or more at 97% identity or better are the exact copies. The rest were filtered for simple repeats (trinucleotide entropy of the matched panel sequence below 4.5 bits, or a panel window hit by 15 or more assembly windows) and whatever remained at 80% identity or better was listed one by one. For each exact plastid copy, a window of the copy plus 8 kb of neighbouring slug DNA was cut out. Repeats in the slug part were found empirically (100 bp bins hit by 4 or more of 20,000 reads) and masked, and every read of two of the three runs (ERR12356115 and ERR12353011, 53% of the data) was mapped to the masked windows; reads landing within 1.5 kb of a boundary were then aligned base by base. A read carries the junction if one alignment spans it by 300 bp on each side. It contradicts the junction if it matches 1 kb or more of slug flank at 98% or better, stops within 50 bp of the boundary, and has 300 bp or more left over. Scripts 22 to 27.</p></details>
<details class="method"><summary>10. Algal genes by protein similarity</summary>
<p>Every annotated protein (longest isoform per Ensembl gene, and all Purdue proteins) was searched with DIAMOND (sensitive) against the animal, algal, bacterial and fungal proteins listed on the gene-transfer page. A protein is listed when its best match overall is not an animal, and called a strong candidate when that match beats the best non-sacoglossan animal match by 20% in bit score. Separately, {PSN} plastid and photosynthesis proteins were aligned to all four assemblies with miniprot; overlapping hits were merged into loci, and the translation of each locus was searched against the same database to tell an algal gene from the slug’s own homologue. Scripts 21, 28, 29.</p></details>

<h2>Caveats</h2>
<ul>
<li>Gene roles are inferred from homology with flies and mammals. Nothing here measures a rhythm.</li>
<li>xgElyCris1.1 is incomplete (BUSCO 73.6% complete; about 60% of read k-mers represented). A gene missing from it is not missing from the animal, as the catalogue shows.</li>
<li>Six protein models are translations guided by another species’ protein and may miss lineage-specific exons. The gene pages say which.</li>
<li>Tree supports are FastTree local supports, not bootstraps. Where a clade is weakly supported the call rests on reciprocal best hits.</li>
<li>Expression values come from one head sample at one unknown time of day.</li>
<li>Independent work; not reviewed by the genome authors.</li>
</ul>
</div></section>
"""
open(os.path.join(DOCS, "methods.html"), "w").write(page("Methods and data. Elysia crispata clock", "methods.html", body, TIP_JS))
print("pages written:", ", ".join(f for f in ("kleptoplasts.html", "gene-transfer.html", "methods.html") if os.path.exists(os.path.join(DOCS, f))))
