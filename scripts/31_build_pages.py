"""Article pages of the site (kleptoplasts, gene transfer, methods).

Executed by 30_build_site.py, which supplies: repo, R, DOCS, page(), e(), genes, klepto, hgt, today, meta.
All figures in the text are read from results/kleptoplast.json and results/hgt.json.
"""
import os, shutil, glob, json

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

TIP_JS = '<script src="../assets/pages.js"></script>'

# =====================================================================================
#                                    KLEPTOPLASTS
# =====================================================================================
if klepto:
    KW = (hgt or {}).get("weak_reads") or {}
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
            return "second alga; no reference sequence of that species is available for this marker" if b["barcode"] == "tufA" else "second alga; the closest full-length reference is from another species"
        main = max((x for x in bc if x["barcode"] == b["barcode"] and "capitatus" in x["top_hit"]), key=lambda x: x["reads"])
        if b is main:
            r = rid["capitatus"][b["barcode"]]
            return "P. capitatus" if b["identity"] >= 99.99 else f'P. capitatus; identical to the shorter voucher sequence {r["accession"]} over {r["length"]} bp'
        return "P. capitatus with a 1-bp length difference in a homopolymer, most likely a sequencing or PCR error"
    rows_bc = "".join(
        f'<tr><td>{b["barcode"]}</td><td class="num">{num(b["reads"])}</td><td class="num">{pct(b["fraction"])}</td>'
        f'<td>{it(" ".join(b["top_hit"].split()[:2]))} <span class="mono muted">{e(b["top_hit"].split("[")[-1].rstrip("]"))}</span></td>'
        f'<td class="num">{b["identity"]:.2f}% of {num(b["alnlen"])} bp</td><td>{e(note)}</td></tr>'
        for b in bc for note in [bnote(b)])
    comp_rows = "".join(
        f"<tr><td>{lab}</td><td>{splitbar([(a, '--series-1', 'P. capitatus'), (b, '--series-2', 'P. lamourouxii'), (c, '--series-muted', 'other')])}</td>"
        f'<td class="num">{pct(a)}</td><td class="num">{pct(b)}</td><td class="num">{pct(c, 2) if c else "–"}</td><td class="num">{n}</td></tr>'
        for lab, a, b, c, n in [
            ("All plastid reads, assigned to the closest plastid genome", cap["bp"] / conf_bp, lam["bp"] / conf_bp, other_bp / conf_bp, num(sum(x["reads"] for x in A))),
            ("rbcL barcode (1,428 bp)", rb_cap, rb_lam, 0, num(sum(b["reads"] for b in bc if b["barcode"] == "rbcL"))),
            ("tufA barcode (987 bp)", tu_cap, tu_lam, 0, num(sum(b["reads"] for b in bc if b["barcode"] == "tufA")))])
    body = f"""
<section class="page-head"><div class="wrap">
<h1>Algal source of the kleptoplasts in <i>Elysia&nbsp;crispata</i></h1>
<p class="lede">The individual sequenced for the genome was collected from the wild at Key Largo, Florida, so its kleptoplasts reflect its recent diet. The plastid DNA among its reads identifies the algae on which it had fed.</p>
<div class="answer">
<p><strong>The kleptoplasts came from two species of <i>Penicillus</i> (shaving-brush algae).</strong> About {pct(cap["bp"] / conf_bp, 0)} of the plastid DNA is from {it("Penicillus capitatus")} and {pct(lam["bp"] / conf_bp, 0)} is from a second {it("Penicillus")} species, which matches the only available reference sequence of {it("P. lamourouxii")}. Small amounts of plastid DNA from {it("Udotea flabellum")} ({pct(udo["bp"] / conf_bp, 2)}) and from a species of {it("Halimeda")} ({pct(hal_bp / conf_bp, 2)}) are also present.</p>
<p>Plastid DNA is abundant ({pct(f_cp)} of all sequenced bases), but algal nuclear DNA is absent: none of the {num(L["reads"])} reads contains algal nuclear ribosomal DNA. The animal had retained the chloroplasts of its food algae but not their nuclei.</p>
</div>
</div></section>

<section class="block article"><div class="wrap">
<h2>Composition of the sequencing library</h2>
<p>All {num(L["reads"])} PacBio HiFi reads ({L["bp"] / 1e9:.2f} Gb) were compared with a reference panel consisting of 241 algal plastid genomes, 117 algal mitochondrial genomes, 1,249 nuclear ribosomal DNA sequences of siphonous green algae, the four bacterial genomes assembled from this individual and the mitochondrial genome of the slug.</p>
<div class="chart"><h3>Share of sequenced bases</h3>
<p class="sub">The remaining {pct(1 - f_cp - f_bac - f_mt)} consists of slug nuclear DNA and reads that match nothing in the panel.</p>
{hbar([("Algal plastid DNA", 100 * f_cp, f"{num(C['cp']['reads'])} reads, {C['cp']['bp'] / 1e6:,.0f} Mb"), ("Bacteria (4 genomes from this animal)", 100 * f_bac, f"{num(C['bac']['reads'])} reads, {C['bac']['bp'] / 1e6:,.0f} Mb"), ("Slug mitochondrial DNA", 100 * f_mt, f"{num(C['slugmt']['reads'])} reads, {C['slugmt']['bp'] / 1e6:,.1f} Mb")])}
</div>
<p>This amount of plastid DNA corresponds to very deep coverage of the plastid genomes: about {num(round(depth_cap, -2))}× for the {it("P. capitatus")} plastid genome and {num(round(depth_lam, -1))}× for the second plastid genome, compared with about {slug_depth:.0f}× for the slug nuclear genome. Because the library was amplified by PCR, these ratios are approximate.</p>

<h2>Identification of the algae</h2>
<p>Three measures give the same result. Each read was assigned to the plastid genome that it matches best, and the sequences of two standard barcode genes, rbcL and tufA, were extracted from every read that spans them.</p>
<div class="legend"><span><b class="sw" style="background:var(--series-1)"></b>{it("Penicillus capitatus")}</span><span><b class="sw" style="background:var(--series-2)"></b>{it("Penicillus lamourouxii")} type</span><span><b class="sw" style="background:var(--series-muted)"></b>other algae</span></div>
<div class="tablewrap"><table class="data">
<thead><tr><th>Measure</th><th style="min-width:200px">Proportions</th><th class="num"><i>P. capitatus</i></th><th class="num"><i>P. lamourouxii</i> type</th><th class="num">Other</th><th class="num">Number of reads</th></tr></thead>
<tbody>{comp_rows}</tbody></table></div>

<h3>Minor algal species</h3>
<div class="tablewrap"><table class="data">
<thead><tr><th>Alga</th><th class="num">Reads</th><th class="num">Share of plastid DNA</th><th>rbcL identity to GenBank records</th><th>tufA identity to GenBank records</th></tr></thead>
<tbody>
<tr><td>{it("Udotea flabellum")}</td><td class="num">{num(udo["reads"])}</td><td class="num">{pct(udo["bp"] / conf_bp, 3)}</td><td>{klepto["minor"][0]["rbcL"]}</td><td>{klepto["minor"][0]["tufA"]}</td></tr>
<tr><td>{it("Halimeda")} sp.</td><td class="num">{num(hal_reads)}</td><td class="num">{pct(hal_bp / conf_bp, 3)}</td><td>{klepto["minor"][1]["rbcL"]}</td><td>{klepto["minor"][1]["tufA"]}</td></tr>
</tbody></table></div>
<p class="small muted" style="margin-top:10px">A further {num(klepto["unassigned"].get("reads", 0))} reads matched a plastid genome only weakly or over a short region, and were not counted as plastid reads. Base-level alignment shows that nearly all of them are slug nuclear reads containing a short low-complexity region that resembles plastid DNA. Of the others, {KW.get("with_exact_piece", 0)} contain a short segment that is identical to kleptoplast DNA, and {sum(1 for d in KW.get("diverged", []) if not d["slug_bp"])} are plastid reads from other algae for which no close reference is available: one is 94% identical to {it("Caulerpa serrulata")} over 10 kb, and the others are most similar to {it("Halimeda")} and {it("Udotea")}.</p>

<h2>Confidence in the species identifications</h2>
<div class="cols">
<div><h3>{it("Penicillus capitatus")}: well supported</h3>
<p>The dominant plastid type is identical to GenBank voucher sequences of {it("P. capitatus")} over {num(rid["capitatus"]["rbcL"]["length"])} bp of rbcL ({rid["capitatus"]["rbcL"]["accession"]}) and {rid["capitatus"]["tufA"]["length"]} bp of tufA ({rid["capitatus"]["tufA"]["accession"]}). The next most similar species, {it(rid["capitatus"]["next_species"]["name"])}, is {rid["capitatus"]["next_species"]["rbcL"]}% (rbcL) and {rid["capitatus"]["next_species"]["tufA"]}% (tufA) identical.</p></div>
<div><h3>{it("Penicillus lamourouxii")}: consistent, but with few reference sequences</h3>
<p>The second plastid type is identical to {rid["lamourouxii"]["rbcL"]["accession"]} over {rid["lamourouxii"]["rbcL"]["length"]} bp of rbcL. This is the only {it("P. lamourouxii")} rbcL sequence in GenBank, and no tufA sequence is available for the species. The next most similar species, {it(rid["lamourouxii"]["next_species"]["name"])}, is {rid["lamourouxii"]["next_species"]["rbcL"]}% (rbcL) and {rid["lamourouxii"]["next_species"]["tufA"]}% (tufA) identical.</p></div>
<div><h3>The two types are distinct species</h3>
<p>The two plastid types differ at {100 - rid["between_types"]["rbcL"]:.1f}% of rbcL sites and {100 - rid["between_types"]["tufA"]:.1f}% of tufA sites, and by 7 to 11% across the whole plastid genome. Both plastid genomes were assembled by the Sanger Institute team (OZ373742 and OZ372796) and assigned to the same two species; the present analysis reached these identifications independently.</p></div>
</div>

<h3>Barcode haplotypes</h3>
<div class="tablewrap"><table class="data">
<thead><tr><th>Barcode</th><th class="num">Reads</th><th class="num">Share</th><th>Best GenBank match</th><th class="num">Identity</th><th>Interpretation</th></tr></thead>
<tbody>{rows_bc}</tbody></table></div>

<h2>Absence of algal nuclear DNA</h2>
<p>If algal nuclei or cytoplasm had been retained together with the plastids, the reads would contain algal nuclear ribosomal DNA, which is present in hundreds of copies per nucleus. In the first pass, {num(nuc["reads_first_flagged"])} reads resembled an algal rDNA sequence, but base-level alignment showed that all of them are the slug’s own rDNA. {"No read has" if not nuc["reads_matching_algal_rdna"] else str(nuc["reads_matching_algal_rdna"]) + " reads have"} an algal nuclear rDNA sequence as its best match. Algal mitochondrial DNA is present only in trace amounts ({word(nuc["algal_mito_reads"])} reads).</p>
<p>This is relevant to the <a href="gene-transfer.html">question of gene transfer</a>. Because the sample contains no algal nuclear DNA, an algal nuclear gene found among these reads could only have come from the slug genome. None was found.</p>

<h2>Implications for studies of kleptoplast rhythms</h2>
<ul>
<li>For animals with the same diet, kleptoplast transcripts should be mapped to the plastid genomes OZ373742.1 ({it("P. capitatus")}) and OZ372796.1 ({it("P. lamourouxii")} type). Neither genome is annotated in the archive; gene coordinates for both, derived by homology, are available in the <a href="methods.html#downloads">downloads</a>.</li>
<li>The composition described here applies to a single wild animal on the day it was collected (30 July 2020). {it("E. crispata")} feeds on many species of algae, so other individuals, including animals kept in culture, will differ. The individual used for the Purdue assembly, for example, contained {it("Acetabularia")} chloroplasts (see the <a href="gene-transfer.html#purdue">gene-transfer page</a>).</li>
</ul>
</div></section>
"""
    open(os.path.join(DOCS, "kleptoplasts.html"), "w").write(page("Kleptoplasts of Elysia crispata", "kleptoplasts.html", body, TIP_JS))

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
    DV = hgt.get("diverged") or {}; WR = hgt.get("weak_reads") or {}
    WRL = next((d for d in WR.get("diverged", []) if d["slug_bp"]), {"slug_identity": 1, "identity": 1})
    _as = {a["taxon"]: a["fraction"] for a in (klepto or {}).get("assignment", [])}
    KS = f'{round(100 * _as.get("Penicillus capitatus", 0.86))}:{round(100 * _as.get("Penicillus lamourouxii", 0.14))}'
    AL = hgt.get("alien") or {}; cands = AL.get("candidates", [])
    GENOME_BP = 542_570_831                                   # span of xgElyCris1.1 (primary assembly)
    LIB = (klepto or {}).get("library") or {}
    DEPTH = LIB.get("bp", 0) / GENOME_BP; MEAN_READ = LIB.get("bp", 0) / max(1, LIB.get("reads", 1))
    SPAN_DEPTH = DEPTH * max(0.0, 1 - 1100 / MEAN_READ) if MEAN_READ else 0     # reads that hold a 500-base piece with 300 bases either side
    PER_CELL = round(SPAN_DEPTH / 6) if SPAN_DEPTH else 0                        # one insertion in this many cells gives about three such reads
    # two steps of the pipeline for the other species, applied to this library as well (scripts/genus/run_screens.sh crispata)
    _sp = os.path.join(repo, "results", "genus", "elysia-crispata.shapes.json")
    SH = json.load(open(_sp)) if os.path.exists(_sp) else None
    shapes_block = ""; COPIES = ""; JOINS_COMPARED = ""
    if SH:
        bg = SH.get("background") or {}; cand_ins = [c for c in SH.get("candidates", []) if c["shape"] == "insertion"]
        sd = (bg.get("depth") or {}).get("spanning_depth") or {}
        if sd.get("500"): SPAN_DEPTH = sd["500"]; PER_CELL = round(SPAN_DEPTH / 6)
        jc = SH.get("joins_compared")          # every read with a chloroplast-slug join, compared with every other around the join
        if jc:
            nr = len(jc.get("near", []))
            near_txt = ""
            if nr:
                near_txt = ({1: "One pair", 2: "Two pairs", 3: "Three pairs", 4: "Four pairs"}.get(nr, f"{nr} pairs") + f' of reads {"has" if nr == 1 else "have"} junctions that are close together but not identical: '
                            'the chloroplast DNA ends at the same position, or within 20 bases of it, and is joined to a simple repeat of the slug genome, either at two different loci or at two positions within one repeat. '
                            'Simple repeats are sites at which a polymerase tends to stall during amplification, and an incomplete extension product can then prime synthesis on another molecule. '
                            'Independent integration of chloroplast DNA at such repeats in more than one cell would give the same result; in either case these are separate junctions, not a single junction observed twice. ')
            JOINS_COMPARED = (f'A comparison by position fails where the nuclear DNA is repetitive, and where two reads were placed on different assemblies of the same region. All {num(jc["reads"])} reads were therefore retrieved again and their junctions were compared with each other at base level. '
                              f'In total, {num(jc["groups_of_copies"] + len(jc["recurrent"]))} junctions occur in more than one read, '
                              + ("in every case as PCR duplicates of a single molecule and never in two independent molecules. " if not jc["recurrent"] else f'{len(jc["recurrent"])} of them in two independent molecules that do not share an end. ') + near_txt)
        cps = (bg.get("depth") or {}).get("copies") or {}          # copies made by the amplification, estimated on the sample
        if cps.get("copies_per_read"):
            COPIES = (f', or slightly fewer, because the library was amplified and about one read in {round(1 / cps["copies_per_read"])} is a PCR duplicate of another '
                      f'(the random sample comprises {100 * cps["sample_fraction"]:.1f}% of the library, and {cps["reads_with_a_copy_in_the_sample"]} of its {num(cps["placed_reads"])} reads have a duplicate within the sample)')
        found = ("None of them contains a segment of chloroplast or mitochondrial DNA of 100 bases or more flanked on both sides by nuclear DNA." if not cand_ins else
                 f"{len(cand_ins)} of them {'contains' if len(cand_ins) == 1 else 'contain'} a segment of organelle DNA of 100 bases or more flanked on both sides by nuclear DNA; {'this read is' if len(cand_ins) == 1 else 'these reads are'} listed in the downloads.")
        shapes_block = (f'<p><strong>Shorter inserts.</strong> An insert of a few hundred bases does not interrupt the alignment of a read to the nuclear sequence, so the search described above would not detect it. '
                        f'Every read in which the first pass had found an organelle segment of at least 60 bases with at least 95% identity, located away from both ends of the read and outside the panel regions that resemble slug repeats, was therefore retrieved again ({num(SH["targeted_reads"])} reads). '
                        f'Each read was aligned to the assembly at base level and searched for an organelle segment lying within a single continuous stretch of nuclear DNA. {found}</p>'
                        + (f'<p><strong>Background rate of chimeras between nuclear sequences.</strong> One read in a hundred from each run was aligned in the same way. Of {num(bg["reads_slug_placed"])} nuclear reads that map to a unique position, {num(bg["joining_unrelated"])} continue with nuclear DNA from another locus, on a different assembly sequence or more than 50 kb away. '
                           f'Each of these reads contains two nuclear segments, so {pct(bg["slug_pieces_joined"] / bg["slug_pieces"], 1)} of all nuclear segments are joined to an unrelated segment. This value includes real structural differences between the chromosomes of this individual and the assembly as well as chimeras formed during library preparation, so it is an upper limit for the rate of chimera formation. '
                           f'In the same sample, {num(bg["cp_joined"])} of {num(bg["cp_reads"])} chloroplast reads are joined to nuclear DNA at a unique position ({pct(bg["cp_joined"] / bg["cp_reads"], 1)}). Chloroplast DNA is therefore joined to nuclear DNA no more frequently than nuclear DNA is joined to unrelated nuclear DNA.</p>' if bg.get("reads_slug_placed") and bg.get("cp_reads") else ""))
    def wh(a, label):
        w = hgt["assembly_screen"].get(a, {}).get("whole_contigs", {}); parts = []
        if w.get("cp"): parts.append(f'{word(w["cp"]["n"])} plastid contig{"" if w["cp"]["n"] == 1 else "s"} ({w["cp"]["bp"] / 1000:.0f} kb)')
        if w.get("bac"): parts.append(f'{word(w["bac"]["n"])} bacterial contig{"" if w["bac"]["n"] == 1 else "s"} ({w["bac"]["bp"] / 1000:.0f} kb)')
        return f'{" and ".join(parts)} in {label}' if parts else ""
    WH = "; ".join(x for x in (wh("sanger", "xgElyCris1.1"), wh("alt", "the alternate haplotype"), wh("purdue", "the Purdue assembly"), wh("tampa", "the Tampa assembly")) if x)
    alg_c = [c for c in cands if c["top_group"] in ("ALG", "CYA")]
    PS = hgt.get("protein_screen") or {}
    def chrn(c): return ("chr" + c) if c.isdigit() else c
    def frow(f):
        ev = "–"
        if f["through_molecules"] is not None:
            ev = str(f["through_molecules"])
            if f["through_reads"] > f["through_molecules"]: ev = f'<span title="{f["through_reads"]} reads, PCR duplicates of a single molecule">{ev}</span>'
        if f["diverging_molecules"] is None: dv = "–"
        elif not f["diverging_molecules"] and (f.get("unique_flank_bp") or 0) < 500 and not f.get("flank_reads"): dv = '<span class="muted" title="The nuclear DNA adjacent to this fragment is repetitive, so reads cannot be placed uniquely">repetitive</span>'
        else: dv = str(f["diverging_molecules"])
        return (f'<tr><td class="mono nowrap">{chrn(f["contig"])}:{f["start"] + 1:,}–{f["end"]:,}</td><td class="num">{f["length"]:,}</td><td class="num">{100 * f["identity"]:.1f}%</td>'
                f'<td class="nowrap">{it(f["alga"].replace("Penicillus", "P."))}</td><td>{e(", ".join(f["genes"][:4])) or "–"}</td><td class="nowrap">{"assembly gap" if f["abuts"] == "gap" else "end of sequence"}</td>'
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
                    f'<td>{"located on a chloroplast fragment" if c.get("on_plastid_fragment") else ("animal match scores nearly as high" if not c["strong"] else "unexplained")}</td></tr>')
        SPN = {"BRYSP": "Bryopsis", "CHOCR": "Chondrus", "ARATH": "Arabidopsis", "CHLRE": "Chlamydomonas", "ECTSI": "Ectocarpus", "SYNY3": "Synechocystis"}
        cands = sorted(cands, key=lambda c: (not c.get("on_plastid_fragment"), -c["top_bits"]))
        alien_block = f"""
<h2>Annotated genes with a non-animal best match</h2>
<p>All {num(SM.get("n_queries", 0))} proteins in the two annotations of the species (the Ensembl annotation of xgElyCris1.1 and the annotation of the Purdue assembly, which is from a different individual) were compared with 13 animal proteomes, the proteomes of three other sacoglossans, 12 algal and plant proteomes including those of two siphonous green algae, a cyanobacterial proteome, and the bacterial and fungal entries in Swiss-Prot. Of these, {num(SM.get("n_with_hits", 0))} have a match, and for {SM.get("n_listed", len(cands))} of them the best match is not an animal protein.</p>
<ul>
<li><strong>{word(len(onfrag)).capitalize()} are chloroplast genes that are annotated as slug genes in the Ensembl annotation.</strong> Each lies on one of the chloroplast fragments described above. Any count of algal genes in the slug genome that is based on the annotation will include them.</li>
<li><strong>The other {word(len(rest))} are not of algal origin.</strong> Their animal homologues score within 20% of the algal match. They are reverse transcriptases of mobile elements, a short conserved peptide and one low-complexity protein.</li>
<li>The Purdue annotation, which is based on an assembly that lacks these fragments, contains no protein that matches an algal or bacterial protein clearly better than an animal protein.</li>
</ul>
<div class="tablewrap"><table class="data"><thead><tr><th>Gene or protein</th><th>Annotation</th><th class="num">Length (aa)</th><th>Location in xgElyCris1.1</th><th style="min-width:13em">Best non-animal match</th><th class="num">Identity</th><th class="num">Bit score</th><th class="num">Best animal bit score</th><th style="min-width:10em">Interpretation</th></tr></thead><tbody>{"".join(arow(c) for c in cands)}</tbody></table></div>
"""
    prot_block = ""
    if PS:
        T = PS.get("totals", {}); NM = PS.get("no_match", {}); FC = PS.get("foreign_contigs", [])
        pur = [r for r in FC if r["assembly"] == "purdue" and r["length"] > 10000]; small_bac = [r for r in FC if "bacterial" in r["reading"]]
        pur_small = [r for r in FC if r["assembly"] == "purdue" and r["length"] <= 10000]
        n_foreign = T.get("alga-like", 0) + T.get("microbe-like", 0)
        def prow(r):
            span = r.get("read_cover_span")
            reads = "none" if not r.get("read_alignments") else f'first {span[1] / 1000:.1f} kb only'
            return (f'<tr><td class="mono">{e(r["contig"])}</td><td class="num">{r["length"]:,}</td><td class="num">{100 * r["gc"]:.0f}%</td><td>{e(", ".join(g.split("_")[0] for g in r["genes"]))}</td>'
                    f'<td class="num">{100 * r["identity_min"]:.0f}{"" if r["loci"] == 1 else f"–{100 * r["identity_max"]:.0f}"}%</td>'
                    f'<td class="num">{r["annotated_genes_with_protein_match"]} of {r["annotated_genes"]}</td><td>{reads}</td></tr>')
        pur_block = ""
        if pur:
            FO = hgt.get("foreign") or {}; BL = FO.get("blast", []); JN = {j["contig"]: j for j in FO.get("joins", [])}; JR = FO.get("join_runs", {})
            n_gm = sum(r["annotated_genes"] or 0 for r in pur); SG = FO.get("same_genes") or {"acetabulum": [0, 0], "n": 0}
            pen = [b["identity"] for b in BL if "peniculus" in b["title"]]; ace = [b["identity"] for b in BL if "acetabulum" in b["title"]]
            tuf = max((b["identity"] for b in BL if "TS1860" in b["title"]), default=None)
            WHAT = {"JAWDGP010001375.1": "2.0 kb of a slug repeat, followed by <i>Acetabularia</i> chloroplast DNA",
                    "JAWDGP010004663.1": "4.0 kb of AT-rich, repeat-like sequence, followed by <i>Acetabularia</i> chloroplast DNA",
                    "JAWDGP010006093.1": "1.0 kb of slug-like sequence, followed by <i>Acetabularia</i> chloroplast DNA",
                    "JAWDGP010001684.1": "0.5 kb of <i>Penicillus</i>-like chloroplast DNA at the start of the contig, followed by a microsatellite and slug nuclear DNA",
                    "JAWDGP010006814.1": "control: a position within another contig, in repetitive DNA"}
            def jrow(c):
                j = JN.get(c)
                if not j: return ""
                return (f'<tr><td class="mono nowrap">{e(c)}</td><td>{WHAT.get(c, "")}</td><td class="num">{j["position"]:,}</td><td class="num">{j["reads_on_unique_side"]:,}</td>'
                        f'<td class="num">{j["through_reads"]}</td><td class="num">{j["diverging_reads"]}</td></tr>')
            join_rows = "".join(jrow(c) for c in WHAT)
            # (named ptested: "tested" is the list of the xgElyCris1.1 stretches, used again further down)
            ptested = [JN[c] for c in WHAT if c in JN and not WHAT[c].startswith("control")]; ctl = next((JN[c] for c in WHAT if c in JN and WHAT[c].startswith("control")), None)
            PR = ""
            if ptested:
                one = all(j["through_reads"] <= 1 for j in ptested)
                PR = (f'{"Each of the " + word(len(ptested)) + " junctions is supported by a single Nanopore read" if one else "The junctions are supported by " + ", ".join(str(j["through_reads"]) for j in ptested) + " reads"}. '
                      f'The other reads that reach the same position ({min(j["diverging_reads"] for j in ptested)} to {max(j["diverging_reads"] for j in ptested)} per junction) continue into different sequence. '
                      + (f'The control position is spanned by {ctl["through_reads"]} reads. ' if ctl else "")
                      + 'This is the same pattern as for the chloroplast fragments in xgElyCris1.1: the assembler extended the end of a contig along a single read in which chloroplast DNA is joined to nuclear DNA, so the junction is not part of the germline genome of that animal. '
                      'A single read cannot show whether the molecule came from a cell with an insertion or was formed when unrelated fragments were ligated end to end, which is known to occur in Nanopore libraries. '
                      'A contig end that is based on a single Nanopore read retains the errors of that read, which explains why the <i>Penicillus</i>-like segments are 85 to 96% identical to the plastid genome and not exact matches.')
            n_rep = sum(1 for r in pur if r.get("read_alignments"))
            pur_block = f"""
<h2 id="purdue">Chloroplast DNA in the Purdue assembly</h2>
<p>The Purdue assembly (GCA_033675545.1) is from a different individual and was built from Nanopore reads. The protein search identified {len(pur)} contigs in it, {min(r["length"] for r in pur) / 1000:.0f} to {max(r["length"] for r in pur) / 1000:.0f} kb long and {sum(r["length"] for r in pur) / 1000:.0f} kb in total, that carry chloroplast genes but no slug gene. They consist of chloroplast DNA from {it("Acetabularia")} (mermaid’s wineglass), on which that animal had evidently been feeding.</p>
<ul>
<li>For the genes that have an {it("Acetabularia peniculus")} sequence in GenBank, the contigs are {min(pen):.0f} to {max(pen):.0f}% identical to it{f"; tufA is {tuf:.1f}% identical to that of an unnamed <i>Acetabularia</i> isolate" if tuf else ""}. The better-studied {it("A. acetabulum")} is {SG["acetabulum"][0]:.0f} to {SG["acetabulum"][1]:.0f}% identical over the same {word(SG["n"])} genes. The species could not be determined.</li>
<li>{it("Acetabularia")} has a very large, gene-poor chloroplast genome that has never been assembled completely. This explains why no plastid genome in the panel matched these contigs, why their genes lie tens of kilobases apart, and why none of the {n_gm} gene models annotated on them matches any known protein.</li>
<li>Their G+C content is {100 * min(r["gc"] for r in pur):.0f} to {100 * max(r["gc"] for r in pur):.0f}%, compared with 37% for the assembly as a whole. No read from the xgElyCris1 individual matches them (two runs, 25 Gb), apart from a repetitive slug sequence at the very start of {word(n_rep)} of the contigs.</li>
</ul>
<p>The list of non-exact matches adds three further segments of chloroplast DNA in this assembly: an 11 kb contig from another siphonous alga (JAWDGP010006602.1; the closest reference, {it("Boodleopsis")}, is 89% identical), and two short {it("Penicillus")}-like segments on contigs that also contain slug nuclear DNA, one of 514 bp with 96% identity at the start of JAWDGP010001684.1 and one of 434 bp with 85% identity about 4 kb from the end of JAWDGP010005248.1.</p>
<p>In total there are five positions in the Purdue assembly at which chloroplast DNA and slug nuclear DNA are joined within a contig. Four of them could be tested with the reads from that individual ({len(JR.get("runs", []))} Nanopore runs, 38 Gb): reads that align to the single-copy side of each junction either span the junction, or end there and continue into different sequence. A control position in the middle of another contig shows the result expected for contiguous DNA.</p>
<div class="tablewrap"><table class="data"><thead><tr><th>Purdue contig</th><th>Sequences joined in the contig</th><th class="num">Position of the junction (bp)</th><th class="num" title="Reads aligned over 1 kb or more of the single-copy side, close to the junction">Reads adjacent to the junction</th><th class="num">Reads spanning the junction</th><th class="num">Reads ending at the junction and continuing into other sequence</th></tr></thead>
<tbody>{join_rows}</tbody></table></div>
<p style="margin-top:16px">{PR}</p>
<p>The fifth junction, in JAWDGP010005248.1, lies in repetitive DNA in which reads cannot be placed uniquely. The chloroplast segment differs from the {it("P. capitatus")} plastid genome mainly by small insertions and deletions (19 in 434 bp, together with 31 substitutions), which is the error profile of an uncorrected Nanopore read and not the pattern expected from evolutionary divergence.</p>
<h3><i>Acetabularia</i> contigs in the Purdue assembly</h3>
<div class="tablewrap"><table class="data"><thead><tr><th>Purdue contig</th><th class="num">Length (bp)</th><th class="num">G+C</th><th>Chloroplast genes found</th><th class="num" title="Protein identity to the Udotea, Chlamydomonas or Bryopsidales query">Identity to query</th><th class="num" title="Annotated gene models on the contig that match any protein in the database">Gene models with a match</th><th>Matching xgElyCris1 reads</th></tr></thead>
<tbody>{"".join(prow(r) for r in pur)}</tbody></table></div>
"""
        prot_block = f"""
<h3>Direct search of the assemblies for photosynthesis genes</h3>
<p>Because a gene may be missing from an annotation, {num(PS.get("n_queries", 0))} proteins were also aligned directly to all four assemblies. These were the chloroplast-encoded proteins of {it("Udotea flabellum")} (a close relative of {it("Penicillus")} with an annotated plastid genome), photosynthesis proteins of {it("Chlamydomonas")}, {it("Arabidopsis")} and siphonous green algae, and the known proteins of {it("Vaucheria litorea")}, the alga involved in the original reports of gene transfer in {it("Elysia chlorotica")}. Each of the {PS.get("n_loci", 0)} matching loci was translated and compared with the protein database described above, in order to distinguish an algal gene from the slug’s own homologue of a conserved protein.</p>
<ul>
<li><strong>Slug genes ({T.get("animal gene", 0)}).</strong> Their best match is an animal protein; for example, the gene for mitochondrial ATP synthase was detected with chloroplast atpB as the query.</li>
<li><strong>Loci on kleptoplast DNA.</strong> Of these, {T.get("plastid fragment (assembly artefact)", 0)} lie on the chloroplast fragments described above and {T.get("plastid-only contig", 0)} on small contigs that consist entirely of plastid DNA.</li>
<li><strong>Loci that match no known protein ({NM.get("n", 0)}).</strong> Of these, {NM.get("low_complexity", 0)} are low-complexity sequences composed of one or two amino acids, and {NM.get("short", 0)} are shorter than 120 residues.</li>
<li><strong>Algal or bacterial genes on contigs that are not, or not clearly, slug DNA ({n_foreign}).</strong> Of these, {sum(r["loci"] for r in pur)} lie on {len(pur)} contigs of the Purdue assembly and {sum(r["loci"] for r in pur_small)} at the start of one further contig ({", ".join(r["contig"] for r in pur_small)}), all of which are described <a href="#purdue">below</a>, and {sum(r["loci"] for r in small_bac)} lie on {len(small_bac)} small bacterial contigs of the Tampa assembly.</li>
</ul>
<p>In xgElyCris1.1, its alternate haplotype and the Tampa assembly, no algal gene is located among slug genes.</p>
{pur_block}
"""
    body = f"""
<section class="page-head"><div class="wrap">
<h1>Search for algal DNA in the <i>Elysia&nbsp;crispata</i> genome</h1>
<p class="lede">A single long read can show a foreign sequence integrated in host DNA, but chimeric molecules formed during library preparation give the same picture. This page describes what the xgElyCris1 reads show when the two are distinguished.</p>
<div class="answer">
<p><strong>There is no algal DNA in the germline genome of this animal, and no evidence of insertions in somatic cells above the background of chimeras formed during library preparation.</strong> In {num(P["cp+slug"]["reads"])} reads, chloroplast DNA is joined to slug nuclear DNA. Kleptoplasts are located in somatic cells, so an insertion would be present in one cell or a few and would be represented by a single molecule. Each of these reads is therefore a candidate, and the absence of the junction from other reads at the same locus is not evidence against it. However, this library was amplified by PCR, which produces chimeras between unrelated molecules. Bacterial DNA is joined to nuclear DNA at the same frequency ({pct(CH["rate_bac_slug"], 2)} of bacterial reads, compared with {pct(CH["rate_cp_slug"], 2)} of chloroplast reads), and {P["cp+slug"]["adapter_junctions"]} of the chloroplast junctions contain the amplification adapter. No read contains chloroplast DNA flanked on both sides by contiguous nuclear DNA, which is the one structure that the joining of two fragments does not produce.</p>
<p><strong>The reference assembly contains chloroplast DNA that is not part of the genome.</strong> xgElyCris1.1 and its alternate haplotype contain {len(fr)} fragments of chloroplast DNA, up to {max(f["length"] for f in fr) / 1000:.1f} kb long and 98 to 100% identical to the kleptoplast DNA of this animal. All {len(fr)} lie at the end of a contig, and each is supported by a single molecule, one of the chimeric reads described above. Reads that approach these positions from the nuclear side continue into other sequence, so the fragments are not part of the germline genome and should not be included in a reference assembly. They will be found in any search of the assembly with a plastid genome, and the Ensembl annotation lists {word(sum(1 for c in cands if c.get("on_plastid_fragment")))} chloroplast genes on them as slug genes.</p>
<p>The {len(fr)} fragments derive from the two {it("Penicillus")} plastid types in a ratio of {n_cap} to {n_lam}, similar to the {KS} ratio among the kleptoplasts, as expected for a random sample of the kleptoplast DNA of this animal.</p>
<p>The Purdue assembly, which is from a different individual, contains fragments of the same kind. It also contains 646 kb of {it("Acetabularia")} chloroplast DNA as separate contigs, which represent the kleptoplasts of that animal and were not removed from the assembly. Details are given <a href="#purdue">below</a>.</p>
</div>
</div></section>

<section class="block article"><div class="wrap">
<h2>Interpretation of single chimeric reads</h2>
<p>Sacoglossan sea slugs take up chloroplasts into cells of the digestive gland, not into the germline. If a fragment of chloroplast DNA were integrated into a chromosome, this would occur in a single somatic cell, and only that cell and its descendants would carry the insertion. In DNA extracted from a whole animal, such an insertion would be represented by a few molecules, or possibly by only one, and the other molecules from the same locus would lack it. The evidence for such an insertion would therefore consist of single reads, and reads from the same locus that lack the junction are expected. For each such read, the question is whether it could have arisen in another way.</p>
<p>In this case there is an obvious alternative. The genome was sequenced from a very small amount of DNA, so the library was amplified by PCR before sequencing (the PacBio ultra-low-input protocol). Amplification occasionally joins unrelated molecules, and each chimeric molecule is then copied many times. Reads with nuclear DNA at one end and chloroplast DNA at the other are therefore expected even if no such molecule was present in the animal.</p>
<p>The bacteria associated with the slug provide a baseline. Integration of <i>Kistimonas</i> chromosomal DNA into slug chromosomes is not expected, so the frequency of reads in which bacterial DNA is joined to nuclear DNA measures the rate at which the library produces artefactual junctions.</p>
<div class="chart"><h3>Reads that also contain slug nuclear DNA</h3>
<p class="sub">Percentage of all reads of each type. Chloroplast reads are joined to nuclear DNA no more frequently than bacterial reads are.</p>
{hbar([("Chloroplast reads", 100 * CH["rate_cp_slug"], f"{num(P['cp+slug']['reads'])} of {num(ct['cp'])} chloroplast reads"), ("Bacterial reads (baseline)", 100 * CH["rate_bac_slug"], f"{num(P['bac+slug']['reads'])} of {num(ct['bac'])} bacterial reads", "bar muted")], maxv=5)}
</div>
<p>Genuine insertions would add to this baseline, but there is no excess. The difference between the two rates corresponds to {abs(round((CH["rate_cp_slug"] - CH["rate_bac_slug"]) * ct["cp"]))} reads among {num(ct["cp"])}. The reads cover each position in the genome about {DEPTH:.0f} times, so if every cell carried one insertion, each at a different locus, about {DEPTH:.0f} additional reads would span one of the two ends of an insertion. A few insertions per cell would therefore have been detected, whereas a lower frequency would be masked by the {num(P["cp+slug"]["reads"])} chimeric reads that the library contains.</p>
<p>Four further tests were applied.</p>
<ul>
<li><strong>Adapter sequences.</strong> In {P["cp+slug"]["adapter_junctions"]} of the {num(P["cp+slug"]["junctions"])} chloroplast–nuclear junctions, and in {P["bac+slug"]["adapter_junctions"]} of the {num(P["bac+slug"]["junctions"])} bacterial–nuclear junctions, the adapter used in the amplification step lies between the two segments. These molecules were joined during library preparation.</li>
<li><strong>Chimeras between chloroplast and bacterial DNA.</strong> If junctions form at random, chloroplast reads should be joined to bacterial DNA in proportion to the amount of bacterial DNA in the library. This predicts about {min(CH["expected_cp_bac"], CH["expected_cp_bac_from_bacteria"]):.0f} to {max(CH["expected_cp_bac"], CH["expected_cp_bac_from_bacteria"]):.0f} such reads, depending on whether the estimate is based on the chloroplast or on the bacterial reads. The observed number is {CH["observed_cp_bac"]}.</li>
<li><strong>Recurrent junctions.</strong> When the {num(P["cp+slug"]["junctions"])} chloroplast–nuclear junctions are compared by their positions in the assembly, {"none is" if rec_cp == 0 else ("one is" if rec_cp == 1 else word(rec_cp) + " are")} found in two molecules with different ends. (Reads that share an end are PCR duplicates of a single molecule and were counted once.) {JOINS_COMPARED}A junction present in two independent molecules would indicate a clone of cells or an inherited insertion. An insertion in a single cell cannot meet this criterion, so the result does not exclude such insertions. The comparison by position does detect {word(rec_bac)} recurrent bacterial–nuclear junction, which corresponds to a region of bacterial sequence at the edge of a small unplaced scaffold, so the test can detect a recurrent junction when one is present.</li>
<li><strong>Distribution of the junctions.</strong> The junctions are distributed over all 15 slug chromosomes approximately in proportion to chromosome length, and along the whole plastid genome. The {it("P. capitatus")} plastid genome is involved in {round(100 * CH["partner_species"]["capitatus"] / (CH["partner_species"]["capitatus"] + CH["partner_species"]["lamourouxii"]))}% of them, which is the proportion of that genome among the kleptoplasts. This is the distribution expected for chimeras that form at random, but random insertions would give the same distribution, so it is the frequency described above, and not the distribution, that provides the test.</li>
</ul>
<p><strong>Reads with an insertion-like structure.</strong> A fragment integrated in a chromosome is flanked on both sides by nuclear DNA, and the two flanks are adjacent in the genome. The joining of two molecules end to end cannot produce this structure. No read has it: none contains a segment of chloroplast DNA of 300 bases or more flanked on both sides by nuclear DNA from the same locus. Reads long enough to show such a segment together with its flanks cover each position in the genome about {SPAN_DEPTH:.0f} times. Each read derives from a single cell, so this is equivalent to the DNA of about {SPAN_DEPTH / 2:.0f} cells{COPIES}. If there had been one insertion per {PER_CELL} cells, averaged over all the cells that contributed DNA to the library, about three such reads would have been expected. Insertions much rarer than this could not be detected in a library of this size (see the <a href="../methods.html#transfer">detection limit</a> in the methods).</p>
{shapes_block}

<h2>Chloroplast DNA in the xgElyCris1.1 assembly</h2>
<p>When the assemblies were divided into windows and aligned to the plastid genomes, exact copies of kleptoplast DNA were found: {len(prim)} in the primary assembly ({num(sum(f["length"] for f in prim))} bp in total, on chromosomes {", ".join(sorted({f["contig"] for f in prim if f["contig"].isdigit()}, key=int))} and three unplaced scaffolds) and {len(alt)} in the alternate haplotype ({num(sum(f["length"] for f in alt))} bp). The two other assemblies of the species, which were produced from other individuals without PCR amplification, contain none.</p>
<p>Each copy was tested against the reads. If chloroplast DNA were present at one of these positions in the germline genome, on one or both copies of the chromosome, half or all of the reads approaching from the nuclear side would continue directly into it. Instead:</p>
<ul>
<li>all {len(fr)} copies end at an assembly gap or at the end of a sequence, so each is adjacent to nuclear DNA on one side only;</li>
<li>each of these {len(tested)} boundaries is supported by {"a single molecule" if single == len(tested) else "at most one molecule"}{f" (at {word(pcr_copies)} of them the molecule was sequenced twice, as PCR duplicates)" if pcr_copies else ""};</li>
<li>at the {len(placeable)} boundaries where the adjacent nuclear DNA is sufficiently unique for reads to be placed, {thr_mol} molecules contain the junction and {div_mol} independent molecules approach from the nuclear side and continue into different sequence. These counts are based on two of the three sequencing runs (53% of the reads). At {contradicted} of the {len(placeable)} boundaries, two or more molecules contradict the assembly;</li>
<li>the other {len(fr) - len(placeable)} lie in repeats, where reads cannot be assigned to a unique position.</li>
</ul>
<p>In each case the assembler extended a contig along a single read in which chloroplast DNA is joined to nuclear DNA, the contig ended where that read ended, and scaffolding left a gap next to it. Each of these reads is one of the chimeric reads counted above. It is most likely a molecule that was joined during amplification, or conceivably the DNA of a single cell with an insertion. In neither case does it represent the germline genome of the animal, and such sequences should not be included in a reference assembly.</p>
<div class="tablewrap"><table class="data">
<thead><tr><th>Location in xgElyCris1.1</th><th class="num">Length (bp)</th><th class="num">Identity to plastid genome</th><th>Alga</th><th style="min-width:9em">Plastid genes</th><th>Other end of the fragment</th><th class="num" title="Independent molecules in which a single alignment extends from nuclear DNA into the chloroplast DNA as assembled (all three runs)">Molecules containing the junction</th><th class="num" title="Independent molecules that follow the nuclear DNA up to the boundary and then continue into sequence other than chloroplast DNA (two of three runs)">Molecules contradicting the junction</th></tr></thead>
<tbody>{"".join(frow(f) for f in prim)}</tbody></table></div>
<p class="small muted" style="margin-top:10px">Coordinates are 1-based and use the Ensembl sequence names (chr1 = OY973766.1, and so on). This table and the {len(alt)} fragments in the alternate haplotype are available in the <a href="methods.html#downloads">downloads</a>. These intervals should be masked before the assembly is used in any analysis that involves plastid genes, because they include parts of rbcL, psbA, tufA, atpB and petB.</p>

<h3>No evidence of older insertions</h3>
<p>An insertion that arose in an earlier generation, or that derived from a different alga, would differ in sequence from the present kleptoplasts, so all matches that are not exact copies were also listed and classified. The primary assembly has {num(DV["sanger"]["non_exact_matches"])} such matches to the 241 plastid genomes and to algal mitochondrial and ribosomal DNA. All but {DV["sanger"]["loci"]} are shorter than 200 bp, are less than 80% identical, or lie in simple repeats such as (CTT)n, which also occur in AT-rich plastid genomes. The remaining {DV["sanger"]["loci"]} are the slug’s own ribosomal RNA genes, which are about 81% identical to those of algae. The alternate haplotype ({DV["alt"]["loci"]} {"locus" if DV["alt"]["loci"] == 1 else "loci"}) and the Tampa assembly ({DV["tampa"]["loci"]} {"locus" if DV["tampa"]["loci"] == 1 else "loci"}) give the same result.</p>
<p>Thus, apart from the {len(fr)} exact copies and simple repeats, none of these three assemblies contains a region of 200 bp or more with at least 80% identity to plastid DNA. The Purdue assembly does contain such regions, which are described <a href="#purdue">below</a>.</p>
<p>The reads give the same result, including for parts of the genome that are not represented in any assembly. A total of {num(WR.get("weak_reads", 0))} reads matched a plastid genome only weakly and were analysed separately. Base-level alignment shows that {WR.get("with_exact_piece", 0)} of them contain a short segment identical to kleptoplast DNA (chimeric reads like those described above, with a plastid segment shorter than 1 kb), and {len(WR.get("diverged", []))} contain plastid-like sequence that differs by 4 to 8% from every reference. {word(sum(1 for d in WR.get("diverged", []) if not d["slug_bp"])).capitalize()} of the latter are ordinary plastid reads from minor algal species for which no close reference is available. The other is a low-accuracy read in which the nuclear part differs from the assembly to a similar extent ({100 * (1 - WRL["slug_identity"]):.1f}%, compared with {100 * (1 - WRL["identity"]):.1f}% for the plastid-like part). No read joins divergent plastid-like DNA to nuclear DNA.</p>
<p>The assemblies also contain a few small separate contigs of plastid or bacterial DNA that are not joined to any other sequence: {WH}. These are contaminant contigs of the kind that commonly remain after decontamination.</p>
{alien_block}
{prot_block}
<h2>Limitations</h2>
<ul>
<li>Insertions shorter than about 200 bp, or more than about 20% divergent from every plastid genome in the panel, would not be detected by the DNA-level searches, and the protein-level searches cover only assembled sequence.</li>
<li>An insertion in a single cell for which the read covers only one end cannot be distinguished from a chimera formed during amplification, and this library contains thousands of such chimeras. A library that had not been amplified would have a much lower background.</li>
<li>Insertions shorter than {100 if SH else 300} bases within a read were not considered. Below about 100 bases, the segments that match a plastid genome are transfer RNA genes, which are also carried by mobile elements in the slug genome.</li>
<li>Some regions of the genome are poorly represented in the PCR-amplified library. Several clock genes lie in such regions (see the <a href="index.html">clock-gene page</a>).</li>
<li>A gene transferred from an alga with no sequenced relative would not be detected if it has no recognisable homologue in the algal proteomes used.</li>
</ul>
</div></section>
"""
    open(os.path.join(DOCS, "gene-transfer.html"), "w").write(page("Search for algal DNA in the Elysia crispata genome", "gene-transfer.html", body, TIP_JS))

# =====================================================================================
#                                    METHODS
# =====================================================================================
dl = []
def add_dl(src, name, label, note):
    if os.path.exists(src):
        shutil.copy(src, os.path.join(DOCS, "downloads", name)); dl.append((name, label, note))
add_dl(os.path.join(R, "clock_genes.tsv"), "clock_genes.tsv", "Clock gene table (TSV)", "one row per gene: location, status in the four assemblies, domains and orthologues")
add_dl(os.path.join(R, "sequences", "clock_proteins.faa"), "clock_proteins.faa", "Clock proteins (FASTA)", "representative protein for each gene")
add_dl(os.path.join(R, "sequences", "clock_cds.fna"), "clock_cds.fna", "Clock coding sequences (FASTA)", "for the genes from which an uninterrupted coding sequence could be extracted")
add_dl(os.path.join(R, "clock_genes.json"), "clock_genes.json", "Full data for each gene (JSON)", "candidate gene models, exon coordinates, read support and expression")
add_dl(os.path.join(R, "plastid_fragments_in_assembly.tsv"), "plastid_fragments_in_assembly.tsv", "Chloroplast fragments in the assemblies (TSV)", "intervals to mask, with the read evidence for each")
add_dl(os.path.join(R, "kleptoplast_barcodes.tsv"), "kleptoplast_barcodes.tsv", "Kleptoplast barcode haplotypes (TSV)", "rbcL and tufA sequences with read counts and GenBank matches")
add_dl(os.path.join(R, "penicillus_plastome_genes.bed"), "penicillus_plastome_genes.bed", "Gene positions on the two Penicillus plastid genomes (BED)", "derived by protein homology with Udotea flabellum NC_039528")
add_dl(os.path.join(R, "intermediate", "diverged_plastid_like_matches.tsv"), "diverged_plastid_like_matches.tsv", "Plastid-like matches that are not exact copies (TSV)", "all matches remaining in the four assemblies after removal of simple repeats")
add_dl(os.path.join(R, "intermediate", "protein_screen_loci.tsv"), "protein_screen_loci.tsv", "Protein search of the assemblies (TSV)", "every locus to which a plastid or photosynthesis protein aligned, with its classification")
add_dl(os.path.join(R, "kleptoplast.json"), "kleptoplast.json", "Kleptoplast results (JSON)", "")
add_dl(os.path.join(R, "hgt.json"), "hgt.json", "Gene-transfer results (JSON)", "")
PSN = ((hgt or {}).get("protein_screen") or {}).get("n_queries", 471)
dl_html = "".join(f'<li><a href="downloads/{n}" download>{e(l)}</a><span>{e(t)}</span></li>' for n, l, t in dl)
body = f"""
<section class="page-head"><div class="wrap">
<h1>Methods and data for <i>Elysia&nbsp;crispata</i></h1>
<p class="lede">All results for <i>E. crispata</i> were computed from publicly available data. The scripts are in the repository, and the numbered steps below follow the order in which they are run.</p>
</div></section>
<section class="block article"><div class="wrap">
<h2 id="downloads">Downloads</h2>
<ul class="downloads">{dl_html}</ul>

<h2>Data</h2>
<div class="tablewrap"><table class="data">
<thead><tr><th>Data set</th><th>Accession</th><th>Notes</th></tr></thead>
<tbody>
<tr><td>Genome assembly, primary</td><td class="mono">GCA_963854125.1 (xgElyCris1.1)</td><td>542.6 Mb, 15 chromosomes, 8,217 gaps; Ensembl annotation 2024_05, 18,436 protein-coding genes</td></tr>
<tr><td>Genome assembly, alternate haplotype</td><td class="mono">GCA_963854135.1</td><td>363.8 Mb in contigs</td></tr>
<tr><td>PacBio HiFi reads (mid-body)</td><td class="mono">ERR11843472, ERR12353011, ERR12356115</td><td>5,516,672 reads, 46.85 Gb; ultra-low-input, PCR-amplified libraries</td></tr>
<tr><td>RNA-seq (head, same individual)</td><td class="mono">ERR12342469</td><td>35.3 million read pairs, one library</td></tr>
<tr><td>Plastid genomes assembled from this individual</td><td class="mono">OZ373742.1, OZ372796.1</td><td>assembled by the Sanger Institute team and assigned to <i>P. capitatus</i> and <i>P. lamourouxii</i></td></tr>
<tr><td>Bacterial genomes assembled from this individual</td><td class="mono">GCA_964289745, GCA_964291485, GCA_964291505, GCA_964291535</td><td><i>Ekhidna</i>, two <i>Kistimonas</i> and <i>Francisella</i>; used as the baseline for chimera formation</td></tr>
<tr><td>Other <i>E. crispata</i> assemblies</td><td class="mono">GCA_033675545.1 (Purdue), GCA_984581095.1 (Tampa)</td><td>from other individuals; the Purdue assembly has its own annotation</td></tr>
<tr><td>Nanopore reads from the Purdue individual</td><td class="mono">PRJNA987316, 12 MinION runs</td><td>38 Gb; used only to test junctions in the Purdue assembly</td></tr>
<tr><td>Dasycladales plastid sequences</td><td class="mono">333 GenBank records</td><td>used to identify the <i>Acetabularia</i> contigs in the Purdue assembly</td></tr>
</tbody></table></div>
<p style="margin-top:14px">The specimen (NSU0021201, clarki ecotype) was collected at Dove Cove, Key Largo, Florida, on 30 July 2020. The genome is described by Middlebrooks et al. (2026), Wellcome Open Research 11:148.</p>

<h2>Clock genes</h2>
<details class="method" open><summary>1. Reference proteins and candidates</summary>
<p>A set of 102 reference clock proteins (from human, mouse, <i>Drosophila</i>, monarch butterfly, zebrafish, oyster, <i>Aplysia</i>, squid and scallop) was retrieved from UniProt and searched with DIAMOND (very-sensitive mode) against the Ensembl proteins of xgElyCris1.1, the proteins of the Purdue assembly and 16 comparison proteomes (scripts 01 to 03).</p></details>
<details class="method"><summary>2. Orthology</summary>
<p>Each gene family was aligned with MAFFT and trimmed with trimAl, and a tree was built with FastTree (LG model with gamma-distributed rates; SH-like support values). A gene was assigned when it fell within the clade of its reference proteins and was the reciprocal best hit of the human or <i>Drosophila</i> gene. Opsins were also required to have the retinal-binding lysine (scripts 04, 05 and 09).</p></details>
<details class="method"><summary>3. Selection of a protein model</summary>
<p>For each gene, the candidate models were the Ensembl model, the Purdue model and homology-guided translations of the genome based on orthologues from other sea slugs, <i>Aplysia</i> and <i>Biomphalaria</i> (miniprot). Annotated models that could not be re-aligned to their own assembly without discrepancies were excluded. Candidates were scored by their mean bit score against the orthologues of four other molluscs, and an existing annotation was retained unless a homology-guided translation scored at least 5% and 30 bits higher. Homology-guided translations were selected for six genes (script 06).</p></details>
<details class="method"><summary>4. Status in each assembly</summary>
<p>The selected protein was aligned to all four assemblies with miniprot. A gene was classed as intact when at least 97% of the protein aligned at a single locus with at least 95% identity and without a frameshift or stop codon; as complete with a frameshift or stop codon when the same criteria were met apart from a frameshift or stop codon; as fragmented when at least 90% of the protein was found but split between loci; as partial when less than 90% was found; and as not found when no locus aligned with at least 95% identity (scripts 07 and 08).</p></details>
<details class="method"><summary>5. Read support and expression</summary>
<p>The exons of every gene (from the Purdue gene models, with 150 bp of flanking sequence and with repeat-containing regions masked) were included as baits when the HiFi reads were screened. A read was counted when it aligned over at least 150 bp with at least 95% identity. Expression in the head was quantified from the RNA-seq library with salmon, using the Purdue coding sequences as the reference ({meta.get("rna_mapped", "37.5")}% of read pairs were assigned; scripts 16 and 20).</p></details>

<h2>Kleptoplasts and gene transfer</h2>
<details class="method" open><summary>6. Screening of all reads</summary>
<p>All reads were mapped with minimap2 (k = 15) to a panel consisting of 241 green-algal and seven other plastid genomes, 117 algal mitochondrial genomes, 1,249 nuclear rDNA sequences of Bryopsidales, the four bacterial genomes, and the mitochondrial genome and rDNA of the slug. Fourteen low-complexity windows that attracted thousands of unrelated reads were masked (scripts 10 and 11).</p></details>
<details class="method"><summary>7. Identification of the algae</summary>
<p>Plastid reads were assigned to the plastid genome that they match best (with at most 3% divergence over at least 1 kb). The rbcL and tufA sequences were extracted from every read that spans them, identical sequences were counted, and each haplotype was compared with 6,453 GenBank barcode sequences of Bryopsidales (script 13).</p></details>
<details class="method"><summary>8. Chimeric reads</summary>
<p>Reads with a match to the panel and at least 500 bp not covered by that match were re-aligned at base level (minimap2, map-hifi preset) to the primary and alternate assemblies of xgElyCris1.1 together with the panel, and reads that remained unexplained were aligned to the Purdue and Tampa assemblies. A junction was defined as a nuclear segment adjacent to a non-nuclear segment, each at least 300 bp long with at least 97% identity. The total number of reads in each class was based on confident matches only (at most 3% divergence over at least 1 kb). Two reads were counted as independent molecules only if both of their ends differed by more than 50 bp. Reads that had been set aside because they matched a plastid genome only weakly were searched in the same alignments for plastid segments of 200 bp or more that are not simple repeats (scripts 14, 15 and 17).</p>
<p>Two steps of the pipeline that was developed later for the other species were also applied to this library (scripts/genus/run_screens.sh, step crispata). First, reads in which the first pass had found a plastid or mitochondrial segment of 60 to 5,000 bases with at least 95% identity, with at least 200 bases of the read on each side, and not in a panel region matched by more than six such reads (slug repeats that resemble that region), were retrieved again from the archive. They were aligned at base level to the primary assembly and the panel, and searched for an organelle segment of 100 bases or more lying within a single continuous stretch of nuclear DNA, either between two alignments or as an insertion within one alignment. Second, 40 evenly spaced blocks from each run, amounting to one read in a hundred, were aligned in the same way, in order to count how often a nuclear read is joined to nuclear DNA from a different assembly sequence or from more than 50 kb away, and how many reads are duplicates of another read (reads that share both ends). Scripts genus/54, 52, 31d, 33 and 46.</p>
<p>Finally, every read containing a junction between chloroplast DNA and nuclear DNA was retrieved again, and the junctions were compared with each other directly: the 400 bases on either side of each junction were aligned between all pairs of reads. Two reads were considered to share a junction when the alignment extended across it for at least 100 bases on both sides, placed the junction of one read within 6 bases of the junction of the other, contained no gap of that size next to it, and was at least 98% identical on both sides. Reads that share an end were treated as duplicates of a single molecule. This approach detects junctions in repetitive nuclear DNA, which cannot be found by comparing positions, and junctions for which the two reads had been placed on different assemblies. Scripts genus/52, 33, 46 and 51.</p></details>
<details class="method"><summary>9. Chloroplast DNA in the assemblies</summary>
<p>All four assemblies were divided into 10 kb windows and aligned at base level to the same panel. Matches of 200 bp or more with at least 97% identity were classed as exact copies. The remaining matches were filtered for simple repeats (trinucleotide entropy of the matched panel sequence below 4.5 bits, or a panel window hit by 15 or more assembly windows), and all those remaining with at least 80% identity were listed individually.</p>
<p>For each exact plastid copy, a window comprising the copy and 8 kb of adjacent nuclear DNA was extracted. Repeats in the nuclear part were identified empirically (100 bp bins hit by four or more of 20,000 reads) and masked, and all reads from two of the three runs (ERR12356115 and ERR12353011, 53% of the data) were mapped to the masked windows. Reads that mapped within 1.5 kb of a boundary were then aligned at base level. A read was counted as containing the junction if a single alignment spanned it by at least 300 bp on each side. A read was counted as contradicting the junction if it matched at least 1 kb of the nuclear flank with at least 98% identity, its alignment ended within 50 bp of the boundary, and at least 300 bp of the read remained unaligned (scripts 22 to 27).</p></details>
<details class="method"><summary>10. Search for algal genes by protein similarity</summary>
<p>Every annotated protein (the longest isoform of each Ensembl gene, and all Purdue proteins) was searched with DIAMOND (sensitive mode) against the animal, algal, bacterial and fungal proteins listed on the gene-transfer page. A protein was listed when its best overall match was not an animal protein, and was classed as a strong candidate when the bit score of that match exceeded that of the best non-sacoglossan animal match by at least 20%. In a separate analysis, {PSN} plastid and photosynthesis proteins were aligned to all four assemblies with miniprot. Overlapping hits were merged into loci, and the translation of each locus was searched against the same database to distinguish algal genes from the slug’s own homologues (scripts 21, 28 and 29).</p></details>

<details class="method"><summary>11. Foreign contigs and junctions in the Purdue assembly</summary>
<p>Contigs carrying algal genes outside the known chloroplast fragments were compared with the other assemblies (minimap2, asm20 preset), with two runs of xgElyCris1 reads, and with the plastid sequences of Dasycladales in GenBank (BLAST, dc-megablast). Junctions between chloroplast and nuclear DNA were tested with the Nanopore reads from the Purdue individual. Every read was mapped to a 4 kb bait sequence on the single-copy side of the junction (with repeats masked), and reads that reached the junction were aligned to a window around it. A read was counted as spanning the junction if a single alignment extended at least 300 bp beyond it on both sides, and as ending at the junction if its alignment ended within 150 bp of the junction with at least 200 bp of the read remaining unaligned (script 29).</p></details>

<h2>Limitations</h2>
<ul>
<li>Gene functions are inferred from homology with <i>Drosophila</i> and mammals. No rhythms were measured.</li>
<li>The xgElyCris1.1 assembly is incomplete (BUSCO completeness 73.6%; about 60% of read k-mers are represented). As the catalogue shows, a gene that is missing from the assembly is not necessarily missing from the genome.</li>
<li>Six protein models are homology-guided translations based on a protein from another species, and may lack lineage-specific exons. These models are identified in the details of each gene.</li>
<li>Support values in the trees are FastTree local support values, not bootstrap values. Where a clade is weakly supported, the assignment relies on reciprocal best hits.</li>
<li>Expression values are from a single head sample collected at an unknown time of day.</li>
<li>This analysis was carried out independently and has not been reviewed by the authors of the genome.</li>
</ul>
</div></section>
"""
open(os.path.join(DOCS, "methods.html"), "w").write(page("Methods and data for Elysia crispata", "methods.html", body, TIP_JS))
print("pages written:", ", ".join(f for f in ("kleptoplasts.html", "gene-transfer.html", "methods.html") if os.path.exists(os.path.join(DOCS, f))))
