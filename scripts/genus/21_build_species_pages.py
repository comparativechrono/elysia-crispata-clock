#!/usr/bin/env python3
"""Build the "Kleptoplasts" and "Gene transfer" pages of every species except E. crispata
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
word = chrome.word
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
KIND = {"hifi": "PacBio HiFi reads", "clr": "PacBio RS reads", "short": "paired-end short reads"}
CLS_WORD = {"cp": "plastid", "slugmt": "slug mitochondrial", "algmt": "algal mitochondrial", "algrDNA": "algal nuclear ribosomal", "bait": "bait"}

def per_cell(k):
    """words for a mean number of insertions per cell"""
    if k >= 1.5: return f"about {k:,.0f} insertions per cell" if k >= 10 else f"about {k:.1f} insertions per cell".replace(".0 ", " ")
    if k >= 0.75: return "about one insertion per cell"
    return f"one insertion per {round(1 / k):,} cells"

def fmt(text, ctx):
    """fill {name} placeholders of a hand-written sentence from the numbers of this species"""
    try: return text.format(**ctx)
    except KeyError as ex:
        if os.environ.get("DRAFT") == "1":                       # while results are still coming in: show what is missing
            class _D(dict):
                def __missing__(self, k): return "⟨" + k + "⟩"
            return text.format_map(_D(ctx))
        raise SystemExit(f"config/species_text.json uses an unknown number {ex} in: {text[:80]}")

def has_page(slug, which):
    """a page is written only when its result file exists and its conclusion has been written"""
    f, key = {"klepto": (".kleptoplast.json", "klepto_answer"), "transfer": (".transfer.json", "transfer_answer")}[which]
    return os.path.exists(os.path.join(RG, slug + f)) and bool(TEXT.get(slug, {}).get(key))
def tabs_for(slug):
    t = [("index.html", "Clock genes")]
    if has_page(slug, "klepto"): t.append(("kleptoplasts.html", "Kleptoplasts"))
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
                dp = L["plastid"]["bp"] / pl; dslug = (L["bases"] - L["plastid"]["bp"] - L["mitochondrion"]["bp"]) / span          # the bases that are neither plastid nor mitochondrial, over the assembly
                c[p + "cp_depth"] = num(dp) if dp >= 10 else f"{dp:.1f}"; c[p + "depth_ratio"] = num(dp / dslug) if dp / dslug >= 10 else f"{dp / dslug:.1f}"
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
            if L.get("contigs"):
                c[p + "contigs"] = num(L["contigs"]["n"]); c[p + "contig_kb"] = f'{L["contigs"]["bp"] / 1e3:,.0f}'; c[p + "contig_longest_kb"] = f'{L["contigs"]["longest"] / 1e3:,.1f}'
                for i, g in enumerate(L["contigs"].get("groups", []), 1):          # kinds of plastid contig, most deeply covered first
                    q = f"{p}cg{i}_"; c[q + "n"] = word(g["n"]); c[f"{p}Cg{i}_n"] = word(g["n"]).capitalize(); c[q + "kb"] = f'{g["bp"] / 1e3:,.1f}' if g["bp"] < 1e4 else f'{g["bp"] / 1e3:,.0f}'; c[q + "depth"] = f'{g["depth"]:,.0f}'
                    c[q + "id"] = f'{100 * g["identity_min"]:.0f} to {100 * g["identity_max"]:.0f}%'; c[q + "rel"] = f'{100 * g["relative"]:.0f}%' if g["relative"] >= 0.095 else f'{100 * g["relative"]:.1f}%'
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
                c[f"{p}{k}_chim_about"] = num(round(ch, 2 - len(str(int(ch))))) if ch >= 100 else num(ch)          # to two figures, for estimates
            for k, v in r.items():
                if "examined" in v: c[f"{p}{k}_chim_examined"] = num(v["examined"].get("pairs_also_with_slug", 0))
            for k, tag in (("cp", "cp"), ("slugmt", "mt")):          # which piece lies at the end of the fragment in the mates that read through a join
                ep = (R["chimera"].get("split_end_piece") or {}).get(k)
                if ep: c[f"{p}end_{tag}_org"] = num(ep["organelle"]); c[f"{p}end_{tag}_n"] = num(ep["organelle"] + ep["slug"])
                epc_ = R.get("end_piece_control") or {}          # 55_end_piece_control.py: what cutting at random would show, given how well each piece is found
                ex = (epc_.get("expected_share_organelle_first_found") or {}).get(k, (epc_.get("expected_share_organelle_first") or {}).get(k))
                if ex is not None: c[f"{p}end_{tag}_expected"] = f"{100 * ex:.0f}%"
                fp = (epc_.get("first_pass_finds_short_organelle_piece") or {}).get(k)
                if fp is not None: c[f"{p}fp_{tag}"] = f"{100 * fp:.0f}%"
            epc = (R.get("end_piece_control") or {}).get("mates") or {}
            for k, tag in (("slug", "slug"), ("slugmt", "mt"), ("cp", "cp")):
                if epc.get(k, {}).get("mates"):
                    c[f"{p}endpass_{tag}_start"] = f'{100 * epc[k]["start_passes"] / epc[k]["mates"]:.0f}%'; c[f"{p}endpass_{tag}_end"] = f'{100 * epc[k]["end_passes"] / epc[k]["mates"]:.0f}%'
            if R.get("end_piece_control"): c[p + "endpass_piece"] = str(R["end_piece_control"]["piece"])
            if r.get("cp", {}).get("rate") and r.get("slugmt", {}).get("rate") and "pairs_with_mate" in r["cp"]:
                # how much more often chloroplast DNA is joined to slug DNA than mitochondrial DNA is, and how many insertions in
                # every cell that difference would be if all of it were insertions (pairs at the recurrent places left out)
                at_places = sum(x["pairs"] for x in R["chimera"].get("recurrent", []) if x["kind"] == "cp")
                cp_pairs = r["cp"]["pairs_also_with_slug"] - at_places
                c[p + "cp_mt_ratio"] = f'{cp_pairs / r["cp"]["pairs_with_mate"] / r["slugmt"]["rate"]:.1f}'
                per = (R.get("sensitivity") or {}).get("mixed_pairs_per_insertion_per_cell")
                if per:
                    ex_ = (cp_pairs - r["slugmt"]["rate"] * r["cp"]["pairs_with_mate"]) / per
                    c[p + "excess_per_cell"] = num(round(ex_, 2 - len(str(int(ex_))))) if ex_ >= 10 else (f"{ex_:.0f}" if ex_ >= 1 else "less than one")
            c[p + "recurrent"] = word(len(R["chimera"].get("recurrent", [])))
            span = SPAN.get(sp["slug"])
            for i, L in enumerate((R.get("locus_tests") or {}).get("loci", [])):      # places of the assembly tested with reads
                for side in ("left", "right"):
                    v = L["edges"].get(side, {})
                    if "through" in v:
                        c[f"{p}locus{i + 1}_{side}_through"] = num(v.get("through_molecules", v["through"]) if R["kind"] != "short" else v["through"]); c[f"{p}locus{i + 1}_{side}_stop"] = num(v["stop_inside"] + v["stop_outside"])
                        c[f"{p}locus{i + 1}_{side}_stop_in"] = num(v["stop_inside"]); c[f"{p}locus{i + 1}_{side}_stop_out"] = num(v["stop_outside"])
                if L.get("depth_locus") is not None: c[f"{p}locus{i + 1}_depth"] = num(L["depth_locus"])
            cpl = (R.get("candidate_places") or {}).get("places")          # short reads of another animal at the places of the insertion-shaped long reads
            if cpl:
                for kind_, tag in (("cp", "cand"), ("slugmt", "candmt")):
                    xs = [x for x in cpl if x["insert"] == kind_ and x["reads_crossing"] is not None]
                    if xs:
                        lo, hi = min(x["reads_crossing"] for x in xs), max(x["reads_crossing"] for x in xs)
                        c[f"{p}{tag}_cross"] = num(lo) if lo == hi else f"{num(lo)} to {num(hi)}"; c[f"{p}{tag}_n"] = word(len(xs))
                        c[f"{p}{tag}_joins"] = word(sum(x["joins_here_same_organelle"] for x in xs)); c[f"{p}{tag}_joins_any"] = word(sum(x["joins_here"] for x in xs))
                        # pairs with one mate within 1 kb of a place and the other in this organelle, expected at the library's own rate of joining
                        sens_ = R.get("sensitivity") or {}; gsz = sens_.get("genome_acgt") or sens_.get("genome_bp"); mx = (r.get(kind_) or {}).get("pairs_also_with_slug", 0)
                        if gsz: ex_ = mx * 2000 * len(xs) / gsz; c[f"{p}{tag}_expected"] = f"{ex_:.0f}" if ex_ >= 1.5 else (f"{ex_:.1f}" if ex_ >= 0.05 else f"{ex_:.2g}")
            rjn = (R.get("repeat_joins") or {}).get("positions")          # short reads: mates that go from a plastid position into the repeat
            if rjn:
                for x in rjn:
                    q = f'{p}rj{x["plastid_position"]}_'; c[q + "on"] = num(x["at_the_join"]["mates_going_on"]); c[q + "rep"] = word(x["at_the_join"]["into_the_repeat"])
                    c[q + "ctl_mean"] = f'{x["controls_mean_into_the_repeat"]:.1f}'; c[q + "ctl_max"] = num(x["controls_max_into_the_repeat"])
            for P_ in R["chimera"].get("pairs", []):          # long reads: joins by kind, and how many hold the amplification adapter
                nm = {"cp+slug": "cp", "slug+slugmt": "mt"}.get(P_["pair"])
                if nm: c[f"{p}{nm}_joins"] = num(P_["junctions"]); c[f"{p}{nm}_adapter"] = num(P_.get("adapter_junctions", 0)); c[f"{p}{nm}_no_adapter"] = word(P_["junctions"] - P_.get("adapter_junctions", 0))
            if R.get("shapes"):          # long reads: shapes of the reads that join organelle DNA to slug DNA
                S = R["shapes"]; cands = S.get("candidates", []); bg = S.get("background") or {}
                for k, nm in (("cp", "cp"), ("slugmt", "mt")):
                    n_ins = sum(1 for x in cands if x["shape"] == "insertion" and x["insert"] == k)
                    c[f"{p}ins_{nm}"] = word(n_ins); c[f"{p}Ins_{nm}"] = word(n_ins).capitalize(); c[f"{p}ins_{nm}_n"] = n_ins
                    allr = bg.get("source") == "all reads"
                    tot = (bg.get(f"{nm}_reads") if allr else S["reads_with_piece"].get(k, 0)) or 0; j = (bg.get(f"{nm}_joined") if allr else S.get("reads_joined_to_placed_slug", {}).get(k, 0)) or 0
                    c[f"{p}join_{nm}_any"] = num(S["reads_joining_slug"].get(k, 0))
                    c[f"{p}join_{nm}"] = num(j); c[f"{p}join_{nm}_of"] = num(tot); c[f"{p}join_{nm}_rate"] = pct_auto(j / tot) if tot else "–"
                c[p + "mirror"] = word(sum(1 for x in cands if x["shape"] == "mirror"))
                if bg:
                    c[p + "bg_reads"] = num(bg["reads_slug_placed"]); c[p + "bg_joined"] = num(bg["joining_unrelated"]); c[p + "bg_read_rate"] = pct_auto(bg["joining_unrelated"] / bg["reads_slug_placed"]) if bg["reads_slug_placed"] else "–"
                    c[p + "bg_pieces"] = num(bg["slug_pieces"]); c[p + "bg_pieces_joined"] = num(bg["slug_pieces_joined"]); c[p + "bg_rate"] = pct_auto(bg["slug_pieces_joined"] / bg["slug_pieces"]) if bg["slug_pieces"] else "–"
                    c[p + "bg_fold_rate"] = pct_auto(bg["folding_back"] / bg["reads_slug_only"]) if bg["reads_slug_only"] else "–"
                    c[p + "bg_any_rate"] = pct_auto(bg["joining_unrelated_any"] / bg["reads_slug_only"]) if bg["reads_slug_only"] else "–"
                    c[p + "bg_ins"] = num(bg["insertion_slug_elsewhere"]["unique"]); c[p + "bg_ins_all"] = num(bg["insertion_slug_elsewhere"]["n"]); c[p + "bg_ins_alone"] = num(bg["insertion_slug_elsewhere"]["unique_alone_at_place"])
                    sd = (bg.get("depth") or {}).get("spanning_depth") or {}
                    if sd.get("500"):
                        c[p + "span500"] = f'{sd["500"]:,.1f}' if sd["500"] < 10 else num(sd["500"]); c[p + "per_cell"] = per_cell(6 / sd["500"])
                        if bg.get("depth", {}).get("sample_one_in"): c[p + "one_in"] = num(bg["depth"]["sample_one_in"])
                rj = S.get("recurrent_joins") or {}
                c[p + "rec_joins"] = word(len(rj.get("recurrent", []))); c[p + "rec_groups"] = num(rj.get("groups_same_well_only", 0) + len(rj.get("recurrent", [])))
                cps_ = (((S.get("background") or {}).get("depth") or {}).get("copies")) or {}          # copies of molecules, looked for among the reads of slug DNA
                if cps_.get("pairs_sharing_both_ends") is not None:
                    c[p + "copies_placed"] = num(cps_["placed_reads"]); c[p + "copies_pairs"] = num(cps_["pairs_sharing_both_ends"]); c[p + "copies_next"] = num(cps_["pairs_sharing_one_end_and_just_missing_the_other"])
                unc = sum(v for k, v in S["shapes"].items() if k.startswith("short piece not confirmed"))          # seed matches that are not organelle DNA base by base
                c[p + "unconfirmed"] = word(unc); c[p + "Unconfirmed"] = word(unc).capitalize(); c[p + "per_cell"] = c.get(p + "per_cell", "")
                c[p + "near_pairs"] = word(len(rj.get("near", []))); rp = rj.get("joins_into_repeat") or {}
                if rp:          # joins into one repeat of the slug genome (the satellite of E. chlorotica)
                    c[p + "rep_cp"] = num(rp["molecules"].get("cp", 0)); c[p + "rep_mt"] = num(rp["molecules"].get("slugmt", 0))
                    c[p + "rep_pairs"] = word(sum(1 for x in rp["pairs_within_40_bases"] if x["insert"] == "cp")); c[p + "rep_expected"] = f'{rp["pairs_expected_by_chance"].get("cp", 0):.1f}'
            else:                        # short reads
                for k, nm in (("cp", "cp"), ("slugmt", "mt")):
                    v = r.get(k) or {}
                    c[f"{p}split_{nm}"] = num(v.get("split_mates", 0)); c[f"{p}split_{nm}_rate"] = pct_auto(v.get("split_mates", 0) / v["pairs_with_mate"]) if v.get("pairs_with_mate") else "–"
                    c[f"{p}mix_{nm}_mol"] = num(v.get("molecules", v.get("pairs_also_with_slug", 0))); c[f"{p}split_{nm}_mol"] = num(v.get("split_molecules", v.get("split_mates", 0)))
                smp = R.get("sample") or {}; cc = smp.get("counts") or {}
                if cc.get("pairs_with_slug"):
                    c[p + "bg_pairs"] = num(cc["pairs_with_slug"]); c[p + "bg_apart"] = num(cc.get("slug_pairs_apart", 0)); c[p + "bg_rate"] = pct_auto(cc.get("slug_pairs_apart", 0) / cc["pairs_with_slug"])
                    c[p + "bg_split"] = num(smp.get("split_mates_joining_unrelated_slug", 0)); c[p + "bg_split_rate"] = pct_auto(smp.get("split_mates_joining_unrelated_slug", 0) / cc["pairs_with_slug"])
                    c[p + "insert"] = num(smp.get("insert_median") or 0)
                sens = R.get("sensitivity") or {}
                if sens.get("split_mates_per_insertion_per_cell"): c[p + "split_per_ins"] = num(sens["split_mates_per_insertion_per_cell"])
                if sens.get("depth"): c[p + "slug_depth"] = num(sens["depth"]) if sens["depth"] >= 10 else f'{sens["depth"]:.1f}'
                c[p + "same_break"] = word(len(R["chimera"].get("recurrent_breaks", []))); c[p + "same_join"] = word(len(same_joins(R["chimera"])))
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
def lc(label):
    """a library label within a sentence: lower case, unless it begins with a place name"""
    return label if label.split()[0] in ("Shenzhen", "Hainan") else label[0].lower() + label[1:]

# =====================================================================================
def klepto_page(sp, K, T):
    slug = sp["slug"]; tx = TEXT.get(slug, {}); ctx = context(sp, K, T); libs = K["libraries"]; main = libs[0]
    lede = fmt(tx.get("klepto_lede", "The plastid DNA among the reads of {species} identifies the alga on which the sequenced animal had fed."), ctx)
    answer = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get("klepto_answer", []))
    rows_lib = "".join(f'<tr><td>{e(L["animal"][0].upper() + L["animal"][1:])}</td><td>{e(L["platform"])}</td><td class="mono">{e(", ".join(L["runs"]))}</td><td class="num">{gb(L.get("bases_screened") or L["bases"])}</td><td>{e(L["note"])}</td></tr>' for L in libs)
    charts = ""
    for L in libs:
        rows = [("Algal plastid DNA", 100 * L["plastid"]["fraction"], f'{num(L["plastid"]["reads"] or 0)} {"reads" if L["kind"] != "short" else "read pairs"}, {L["plastid"]["bp"] / 1e6:,.2f} Mb'),
                ("Slug mitochondrial DNA", 100 * L["mitochondrion"]["fraction"], f'{num(L["mitochondrion"]["reads"] or 0)} {"reads" if L["kind"] != "short" else "read pairs"}, {L["mitochondrion"]["bp"] / 1e6:,.2f} Mb')]
        # short reads match algal mitochondrial genomes in conserved genes too; shown only when there is enough to mean something
        if L["algal_mitochondrion"]["bp"] and (L["kind"] != "short" or L["algal_mitochondrion"]["fraction"] >= 5e-6): rows.append(("Algal mitochondrial DNA", 100 * L["algal_mitochondrion"]["fraction"], f'{L["algal_mitochondrion"]["bp"] / 1e3:,.1f} kb'))
        if L["algal_rdna"]["bp"] and L["kind"] != "short": rows.append(("Algal nuclear ribosomal DNA", 100 * L["algal_rdna"]["fraction"], f'{L["algal_rdna"]["bp"] / 1e3:,.1f} kb'))
        digits = 2 if max(r[1] for r in rows) >= 0.1 else 4
        charts += f'<div class="chart"><h3>{e(lib_label(L, len(libs))) + ": share of sequenced bases" if len(libs) > 1 else "Share of sequenced bases"}</h3><p class="sub">In total, {gb(L.get("bases_screened") or L["bases"])} of {e(KIND[L["kind"]])} were analysed. The chart shows the bases that match a reference sequence with at least {"97" if L["kind"] != "clr" else "82"}% identity. All other bases are slug nuclear DNA or match nothing in the panel.</p>{hbar(rows, digits=digits)}</div>'
    donor_tables = ""
    for L in libs:
        if not L["donors"] or L["plastid"]["bp"] < 5000: continue
        rows = "".join(f'<tr><td>{taxon_html(d["name"])} <span class="mono muted">{e(d["accession"])}</span></td><td class="num">{pct_auto(d["fraction"])}</td><td class="num">{num(d["reads"])}</td>'
                       f'<td class="num">{100 * d["median_identity"]:.1f}%</td><td class="num">{100 * d["identity_q10"]:.1f} to {100 * d["identity_q90"]:.1f}%</td></tr>' for i, d in enumerate(L["donors"][:8]) if i == 0 or d["fraction"] >= 0.0002)
        donor_tables += (f'<h3>{e(lib_label(L, len(libs)))}</h3>' if len(libs) > 1 else "") + f'''<div class="tablewrap"><table class="data">
<thead><tr><th>Best-matching plastid reference</th><th class="num">Share of plastid DNA</th><th class="num">Reads</th><th class="num">Median identity</th><th class="num">10th to 90th percentile</th></tr></thead>
<tbody>{rows}</tbody></table></div>'''
    for L in libs:          # short reads: the plastid contigs, sorted into kinds
        G = (L.get("contigs") or {}).get("groups")
        if not G: continue
        def grow(g):
            if g["kind"] == "close": what = f'Within 5% of {taxon_html(g["label"])} <span class="mono muted">{e(", ".join(r["id"].split("|", 1)[-1] for r in g["references"][:2]))}</span>'
            else: what = ("No close reference; " + ("high read depth" if "plentiful" in g["kind"] else "low read depth") + '<br><span class="small muted">closest: ' + ", ".join(taxon_html(x) for x in list(dict.fromkeys(" ".join(r["name"].split()[:2]) for r in g["references"]))[:3]) + "</span>")
            ident = f'{100 * g["identity_min"]:.1f}%' if g["n"] == 1 else f'{100 * g["identity_min"]:.1f} to {100 * g["identity_max"]:.1f}%'
            return (f'<tr><td>{what}</td><td class="num">{num(g["n"])}</td><td class="num">{g["bp"] / 1e3:,.1f} kb</td><td class="num">{ident}</td>'
                    f'<td class="num">{g["depth"]:,.0f}×</td><td class="num">{100 * g["relative"]:.0f}%</td></tr>' if g["relative"] >= 0.095 else
                    f'<tr><td>{what}</td><td class="num">{num(g["n"])}</td><td class="num">{g["bp"] / 1e3:,.1f} kb</td><td class="num">{ident}</td><td class="num">{g["depth"]:,.1f}×</td><td class="num">{100 * g["relative"]:.1f}%</td></tr>')
        contig_tables_one = f'''<div class="tablewrap"><table class="data">
<thead><tr><th>Group of plastid contigs</th><th class="num">Contigs</th><th class="num">Total length</th><th class="num">Identity to the closest reference</th><th class="num">Read depth</th><th class="num">Depth relative to the first row</th></tr></thead>
<tbody>{"".join(grow(g) for g in G)}</tbody></table></div>
<p class="small muted" style="max-width:86ch">The plastid reads were assembled into contigs, and the contigs were grouped. A contig within 5% of a reference sequence is grouped with that reference. The remaining contigs have no close reference and are divided by read depth into those with at least one tenth of the depth of the first row and those with less. Depth is the mean read coverage of the contigs in a group, and indicates the relative abundance of each plastid genome in the animal. Only reads that resemble a reference sequence were assembled, so an alga without a close reference is represented only by its more conserved genes.</p>'''
        donor_tables += (f'<h3>{e(lib_label(L, len(libs)))}: plastid contigs</h3>' if len(libs) > 1 else "<h3>Plastid contigs</h3>") + contig_tables_one
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
        how = ("Each row is a distinct sequence (haplotype) extracted directly from the reads that span the gene." if L["kind"] == "hifi" else
               f'The sequences were extracted from the plastid contigs assembled from these reads ({num(L["contigs"]["n"])} contigs, {L["contigs"]["bp"] / 1e3:,.0f} kb in total).' if L.get("contigs") else "")
        bc_tables += (f'<h3>{e(lib_label(L, len(libs)))}</h3>' if len(libs) > 1 else "") + f'''<p class="small muted" style="max-width:86ch">{how}</p><div class="tablewrap"><table class="data">
<thead><tr><th>Gene</th><th class="num">Length (bp)</th><th>Source</th><th>Best match in GenBank</th><th class="num">Identity</th><th>Most similar GenBank record of 400 bp or more</th></tr></thead>
<tbody>{rows}</tbody></table></div>'''
    nuc = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get("klepto_nuclei", []))
    notes = "".join(f"<li>{fmt(p, ctx)}</li>" for p in tx.get("klepto_caveats", []))
    which = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get("klepto_which", []))
    # downloads made by 45_export_tables.py
    dls = []
    for L in libs:
        tag = slug if L["library"] == "main" else f'{slug}_{L["library"]}'; lab = (" (" + lc(lib_label(L, len(libs))) + ")") if len(libs) > 1 else ""
        if os.path.exists(os.path.join(SITE, "downloads", tag + "_plastid_contigs.fa")): dls.append(f'<a href="../downloads/{tag}_plastid_contigs.fa">plastid contigs{lab}</a>')
        if os.path.exists(os.path.join(SITE, "downloads", tag + "_barcodes.fa")): dls.append(f'<a href="../downloads/{tag}_barcodes.fa">barcode sequences{lab}</a>')
    rhythm = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get("klepto_rhythm", []))
    if rhythm and dls: rhythm += f'<p>Downloads: {", ".join(dls)}, and <a href="../downloads/{slug}_kleptoplast.json">all values shown on this page</a> (JSON).</p>'
    body = f"""
<section class="page-head"><div class="wrap">
<h1>Algal source of the kleptoplasts in {it(sp["species"]).replace(" ", "&nbsp;")}</h1>
<p class="lede">{lede}</p>
<div class="answer">{answer}</div>
</div></section>

<section class="block article"><div class="wrap">
<h2>Sequencing reads</h2>
<div class="tablewrap"><table class="data">
<thead><tr><th>Source of DNA</th><th>Platform</th><th>Runs</th><th class="num">Bases analysed</th><th>Notes</th></tr></thead>
<tbody>{rows_lib}</tbody></table></div>
<p class="small muted" style="margin-top:10px;max-width:90ch">Every read was compared with a reference panel of 1,271 algal plastid sequences (241 complete plastid genomes of green algae, plastid records of Dasycladales, Xanthophyceae and Cladophorales, and several distantly related plastid genomes to detect other matches), 117 algal mitochondrial genomes, 7,352 nuclear ribosomal DNA sequences of algae and the mitochondrial genomes of nine <i>Elysia</i> species. Reads with a match were then aligned at base level to the panel and to the assembly of this species.</p>

<h2>Composition of the sequenced DNA</h2>
{charts}

<h2>Identification of the alga</h2>
{which}
{donor_tables}
<p class="small muted" style="margin-top:10px;max-width:90ch">Each read is assigned to the reference sequence that it matches best. An identity close to 100% indicates that the reference sequence is from the same species as the alga on which the animal had fed. A substantially lower identity means only that the reference is the closest sequenced relative.</p>

{"<h2>Barcode genes (rbcL and tufA)</h2>" + bc_tables if bc_tables else ""}

{"<h2>Algal nuclear and mitochondrial DNA</h2>" + nuc if nuc else ""}

{"<h2>Implications for studies of kleptoplast rhythms</h2>" + rhythm if rhythm else ""}

{"<h2>Limitations</h2><ul>" + notes + "</ul>" if notes else ""}
</div></section>
"""
    src = f"Data for {it(sp['species'])}: reads {e(', '.join(r for L in libs for r in L['runs']))}; assembly {e(sp['asm_name'])} ({e(sp['accession'])}). The methods and data sources are described on the <a href=\"../methods.html\">methods page</a>."
    return chrome.page(f"Kleptoplasts of {sp['species']}", body, root="../", species=sp["species"], tabs=tabs_for(slug), tab_current="kleptoplasts.html", today=today, sources=src, scripts=TIP_JS,
                       description=f"Algal source of the kleptoplasts in the sea slug {sp['species']}, identified from the plastid DNA among the genome sequencing reads.")

# =====================================================================================
RNA_GENE = re.compile(r"ribosomal RNA|^rr[nls]|^rn[ls]$|tRNA|^trn", re.I)
def locus_reading(L, tx):
    key = f'{L["contig"]}:{L["start"]}'
    if key in tx.get("loci_readings", {}): return tx["loci_readings"][key]
    g = L.get("genes") or []
    if any(RNA_GENE.search(x) for x in g): return "A ribosomal or transfer RNA gene. These genes are conserved among bacteria, mitochondria and plastids, so a match with this identity does not indicate an algal origin."
    if L["cls"] == "algal mitochondrial" and L["identity"] < 0.9: return "A weak match to an algal mitochondrial genome; most likely a ribosomal RNA gene of bacterial origin or derived from the slug’s own mitochondrial genome."
    if L["identity"] < 0.9: return "A weak match to a plastid gene (see text)."
    return "Plastid DNA of the type present in this species (see text)."

READS_FRAMING = """<p>Sacoglossan sea slugs take up chloroplasts into cells of the digestive gland, not into the germline. If a fragment of chloroplast DNA were integrated into a chromosome, this would occur in a single somatic cell, and only that cell and its descendants would carry the insertion. In DNA extracted from a whole animal, such an insertion would be represented by a few molecules, or possibly by only one, and the other molecules from the same locus would lack it. A single read in which chloroplast DNA is joined to nuclear DNA is therefore the only evidence that can be expected, and other reads without the junction do not contradict it.</p>
<p>The difficulty is that library preparation also occasionally joins unrelated molecules, and in a single read such a chimera cannot be distinguished from a genuine junction. Three properties of the reads were examined.</p>
<ul>
<li><strong>Read structure.</strong> A fragment integrated in a chromosome is flanked on both sides by nuclear DNA, and the two flanks are adjacent in the genome. A chimera formed during library preparation consists of two molecules joined end to end. A read in which nuclear DNA is followed by organelle DNA and then by the continuation of the same nuclear sequence is described here as having an insertion-like structure. Two independent chimeric junctions that involve adjacent nuclear sequences are not expected by chance, so a read with this structure requires another explanation.</li>
<li><strong>Frequency.</strong> In principle, chimera formation during library preparation affects all DNA in the same way. Chloroplast DNA should then be joined to nuclear DNA as frequently as mitochondrial DNA is, and as frequently as one nuclear segment is joined to an unrelated nuclear segment. Genuine insertions would add to this background.</li>
<li><strong>Recurrence.</strong> A junction found in two molecules that are not duplicates of each other is not a library artefact. It indicates either a clone of cells or an inherited insertion. An insertion in a single cell cannot meet this criterion, so the absence of recurrent junctions does not exclude such insertions.</li>
</ul>"""
ORG_WORD = {"cp": "Chloroplast DNA", "slugmt": "Slug mitochondrial DNA", "algmt": "Algal mitochondrial DNA"}
def org_lc(k):
    """the kind of DNA as it reads in the middle of a sentence"""
    w = ORG_WORD.get(k, k).replace("Slug m", "M")
    return w[0].lower() + w[1:]
def gene_names(names, k=3):
    """short gene symbols out of the names of a feature table (rpoC1, not 'RNA polymerase beta subunit')"""
    out = []
    for g in names or []:
        g = {"rrl": "23S rRNA", "rrs": "16S rRNA", "rrf": "5S rRNA", "rrn16": "16S rRNA", "rrn23": "23S rRNA", "rrn5": "5S rRNA"}.get(g.split("_")[0], g.split("_")[0] if g.startswith("trn") else g)
        if " " in g and not g.endswith("rRNA"): continue
        if len(g) <= 12 and g not in out: out.append(g)
    return ", ".join(out[:k])
def rate_cell(n, tot): return f'<td class="num">{num(n)}</td><td class="num">{pct_auto(n / tot) if tot else "–"}</td>'
def long_reads_block(R, slug):
    S = R.get("shapes")
    if not S: return ""
    sh = S["shapes"]; bg = S.get("background") or {}
    def ins(k): return sh.get(f"insertion|{k}", 0)
    rows = ""
    allreads = bg.get("source") == "all reads"
    cands = S.get("candidates", [])
    n_ins = {k: sum(1 for c in cands if c["shape"] == "insertion" and c["insert"] == k) for k in ("cp", "slugmt", "algmt")}
    for k in ("cp", "slugmt"):
        tot = (bg.get("cp_reads" if k == "cp" else "mt_reads") if allreads else S["reads_with_piece"].get(k, 0)) or 0
        j = (bg.get("cp_joined" if k == "cp" else "mt_joined") if allreads else S.get("reads_joined_to_placed_slug", {}).get(k, 0)) or 0
        if not tot: continue
        rows += f'<tr><td>{ORG_WORD[k]}</td><td class="num">{num(tot)}</td>{rate_cell(j, tot)}<td class="num"><strong>{num(n_ins[k])}</strong></td></tr>'
    if bg:
        src = "all reads of 2 kb or longer" if allreads else f'random sample of one read in {bg["depth"]["sample_one_in"]:,}' if bg.get("depth") and bg["depth"].get("sample_one_in") else "random sample"
        e_ = bg["insertion_slug_elsewhere"]
        rows += (f'<tr><td>Nuclear DNA segments, for comparison <span class="muted small">({src})</span></td><td class="num">{num(bg["slug_pieces"])}</td>{rate_cell(bg["slug_pieces_joined"], bg["slug_pieces"])}'
                 f'<td class="num">{num(e_["unique"])}</td></tr>')
    note = ("The first three columns are based on all reads of 2 kb or longer, counted in the same way, so that the three frequencies are comparable. " if allreads else "")
    note += (f'A read is counted as containing a type of DNA when at least {S["min_segment"]} of its bases match that DNA with at least {100 * S["min_identity"]:.0f}% identity, and as joined to nuclear DNA when a nuclear segment of at least that length is adjacent to the organelle segment in the same read. '
             f'Only junctions with nuclear DNA that maps to a unique position (mapping quality of 30 or higher) are counted. The comparison row counts nuclear segments and not reads, because a read in which two nuclear segments are joined contains two of them. '
             f'The segments are those of reads that map to a unique position, and a segment is counted as joined when the read continues with nuclear DNA from a different assembly sequence or from more than 50 kb away (junctions at the very end of an assembly sequence are excluded). '
             f'This row includes real structural differences between the chromosomes of this individual and the assembly as well as chimeras formed during library preparation, so it is an upper limit for the frequency of chimera formation. '
             f'The last column gives the number of reads with an insertion-like structure; in the comparison row, these are reads in which the inserted segment maps to a unique position elsewhere in the genome.')
    out = f'''<div class="tablewrap"><table class="data"><thead><tr><th>Type of DNA</th><th class="num">Reads containing it (segments, for nuclear DNA)</th><th class="num">Joined to nuclear DNA at a unique position</th><th class="num">Frequency</th><th class="num">Flanked by nuclear DNA on both sides (insertion-like structure)</th></tr></thead><tbody>{rows}</tbody></table></div>
<p class="small muted" style="margin-top:8px;max-width:92ch">{note}</p>'''
    ins_rows = [c for c in cands if c["shape"] == "insertion"]; mir = [c for c in cands if c["shape"] == "mirror"]
    if ins_rows or mir:
        def crow(c):
            if c["shape"] == "insertion":
                place = c.get("place") or {}
                where = f'<span class="mono">{e(c["flank_sequence"])}:{num(place.get("position", c["left_to"]))}</span>'
                gn = gene_names(c.get("piece_genes"))
                what = f'{org_lc(c["insert"])}, <span class="mono">{num(c["piece_from"])}–{num(c["piece_to"])}</span>' + (f' ({e(gn)})' if gn else "")
            else:
                gn = gene_names(c.get("flank_genes"))
                where = f'none: the flanking sequence is chloroplast DNA, at position <span class="mono">{num(min(c["left_to"], c["right_to"]))}</span>' + (f' ({e(gn)})' if gn else "")
                what = f'nuclear DNA, <span class="mono">{e(c["piece_sequence"])}:{num(c["piece_from"])}–{num(c["piece_to"])}</span>'
            oc = c.get("other_reads_crossing")
            return (f'<tr><td class="mono">{e(c["read"])}</td><td>{where}</td><td>{what}</td><td class="num">{num(c["piece_bp"])}</td>'
                    f'<td class="num">{"–" if oc is None else num(oc)}</td><td>{e(tx_note(slug, c))}</td></tr>')
        out += f'''<p style="margin-top:14px"><strong>Reads with an insertion-like structure.</strong> Each row is a single molecule.</p>
<div class="tablewrap"><table class="data"><thead><tr><th>Read</th><th>Position in the slug genome</th><th>Inserted segment</th><th class="num">Length (bp)</th><th class="num" title="Other reads that span the same position in a single alignment, without the insert">Other reads spanning the position without the insert</th><th>Notes</th></tr></thead>
<tbody>{"".join(crow(c) for c in ins_rows + mir)}</tbody></table></div>'''
    d = (bg.get("depth") or {}).get("spanning_depth") if bg else None
    if d:
        s500 = d.get("500") or 0
        if s500 > 0:
            k = 6 / s500
            out += (f'<p class="small muted" style="margin-top:8px;max-width:92ch">Detection limit: reads long enough to contain a 500-base insert with {S["min_segment"]} bases of nuclear DNA on each side cover each position in the genome {s500:,.1f} times. '
                    f'Each read derives from a single cell, so these reads are equivalent to the DNA of about {(f"{s500 / 2:,.0f}" if s500 >= 6 else f"{s500 / 2:,.1f}")} cells. An insertion present in every cell on one of the two chromosome copies would be spanned by about that number of reads, '
                    f'and insertions at different loci in different cells contribute in the same way. At a frequency of {per_cell(k)}, averaged over all the cells that contributed DNA to the library, about three reads with an insertion-like structure would be expected. '
                    f'{copies_note(bg, R)}{TEXT.get(slug, {}).get("could_show_" + R["library"], "")}</p>')
    return out
def copies_note(bg, R):
    """amplified libraries: how many reads are copies of another read, estimated on the random sample"""
    cps = (bg.get("depth") or {}).get("copies") or {}
    if "amplified" not in R["platform"] or not cps.get("placed_reads"): return ""
    if cps["reads_with_a_copy_in_the_sample"]:
        return (f'The library was amplified, and about one read in {round(1 / cps["copies_per_read"])} is a PCR duplicate of another ({cps["reads_with_a_copy_in_the_sample"]} of {num(cps["placed_reads"])} reads in the random sample have a duplicate within the sample), so these values are slight overestimates. ')
    return (f'The library was amplified, which produces duplicates of molecules. None of the {num(cps["placed_reads"])} reads in the random sample has a duplicate within the sample, which indicates that duplicates make up less than {100 * cps["copies_per_read_upper"]:.0f}% of the reads. ')
def short_reads_block(R, slug):
    ch = R["chimera"]; rates = ch["rates"]; smp = R.get("sample") or {}; sens = R.get("sensitivity") or {}
    rows = ""
    for k in ("cp", "slugmt"):
        v = rates.get(k)
        if not v or not v.get("pairs_with_mate"): continue
        tot = v["pairs_with_mate"]
        rows += f'<tr><td>{ORG_WORD[k]}</td><td class="num">{num(tot)}</td>{rate_cell(v["pairs_also_with_slug"], tot)}{rate_cell(v.get("split_mates", 0), tot)}</tr>'
    if smp and smp.get("counts"):
        c = smp["counts"]; ns = c.get("pairs_with_slug", 0); one_in = next((r.get("sample_one_in") for r in R["chimera"].get("runs", [])), None)
        rows += (f'<tr><td>Nuclear DNA, for comparison <span class="muted small">(random sample of the library)</span></td><td class="num">{num(ns)}</td>{rate_cell(c.get("slug_pairs_apart", 0), ns)}'
                 f'{rate_cell(smp.get("split_mates_joining_unrelated_slug", 0), ns)}</tr>')
    out = f'''<div class="tablewrap"><table class="data"><thead><tr><th>Type of DNA</th><th class="num">Read pairs containing it</th><th class="num">Mate in nuclear DNA (discordant pairs)</th><th class="num">Frequency</th><th class="num">Junction within one read (split reads)</th><th class="num">Frequency</th></tr></thead><tbody>{rows}</tbody></table></div>
<p class="small muted" style="margin-top:8px;max-width:92ch">A pair is counted as containing a type of DNA when at least 80 bases of one of its reads match that DNA with at least 97% identity. In the third column the mate maps to a unique position in nuclear DNA; for nuclear DNA itself, the two reads map to different assembly sequences or more than 20 kb apart. In the fifth column a single read contains at least 40 bases of each type of DNA, so the junction itself is sequenced. All such pairs and reads {"that were analysed " if R.get("thinned", 1) > 1 else ""}are listed in the junction table available for download.</p>'''
    if R.get("thinned", 1) > 1:
        ex = (rates.get("cp") or {}).get("examined") or {}
        out += (f'<p class="small muted" style="margin-top:8px;max-width:92ch">Because so much of this library is organelle DNA, the pairs that match an organelle genome were subsampled: every read was mapped and counted, but only one in {R["thinned"]} of the matching pairs, selected by its ordinal position in the run, was analysed further. '
                f'The values in the first two rows are therefore estimates, obtained by multiplying the counts by {R["thinned"]}. In total, {num(ex.get("pairs_with_mate", 0))} chloroplast pairs were analysed, of which {num(ex.get("pairs_also_with_slug", 0))} have their mate in nuclear DNA and {num(ex.get("split_mates", 0))} contain a junction within one read. '
                f'All pairs that match one of the loci under test were retained.</p>')
    if sens:
        sp_ = sens.get("split_mates_per_insertion_per_cell"); mp = sens.get("mixed_pairs_per_insertion_per_cell")
        if sp_:
            cpv = rates.get("cp", {}); n_sp = cpv.get("split_mates", 0)
            out += (f'<p class="small muted" style="margin-top:8px;max-width:92ch">Detection limit: the fragments are about {num(sens["insert_median"])} bases long and each read covers {sens["read_length"]:.0f} bases, so a junction within a fragment is usually sequenced directly. '
                    f'One insertion per cell, at a different locus in each cell, would add about {sp_:,.0f} split reads' + (f' and about {mp:,.0f} discordant pairs' if mp and mp >= 1 else "") +
                    ('. The library contains no read that extends from chloroplast DNA into nuclear DNA.</p>' if not n_sp else
                     f'. For comparison, the library contains {"an estimated " if R.get("thinned", 1) > 1 else ""}{num(n_sp)} split reads that join chloroplast DNA to nuclear DNA, a number that includes the chimeras formed during library preparation.</p>') if sens.get("insert_median") else "")
    rb = [x for x in ch.get("recurrent_breaks", []) if x["kind"] in ORG_WORD]; rec = [x for x in ch.get("recurrent", []) if x["kind"] in ORG_WORD]
    sj = same_joins(ch)
    if rb or rec or sj:
        items = ""
        for x in sj:
            thru = [v for v in x.get("evidence", []) if v["evidence"] == "split mate"]
            pos = num(thru[0]["other_pos"]) if thru else (f'{num(x["other_from"])}–{num(x["other_to"])}' if x["other_to"] - x["other_from"] > 5 else num(x["other_from"]))
            items += (f'<li>{word(x["molecules"]).capitalize()} independent molecules join <span class="mono">{e(x["slug"])}:{num(thru[0]["slug_pos"] if thru else x["slug_from"])}</span> to {org_lc(x["kind"])} '
                      f'(<span class="mono">{e(x["other_t"].split("|", 1)[1])}:{pos}</span>). {tx_join(slug, x)}</li>')
        # breaks at one base of the slug sequence whose other sides lie at several places of a reference (the copies of a repeated gene) are one item
        merged = {}
        for x in rb:
            k_ = (x["kind"], x["slug"], x["slug_break"])
            if k_ in merged: merged[k_]["molecules"] += x["molecules"]; merged[k_]["others"].append(f'{x["other_t"].split("|", 1)[1]}:{num(x["other_break"])}')
            else: merged[k_] = dict(x, others=[f'{x["other_t"].split("|", 1)[1]}:{num(x["other_break"])}'])
        for x in merged.values():
            items += (f'<li>In {word(x["molecules"])} molecules, a read extends from <span class="mono">{e(x["slug"])}:{num(x["slug_break"])}</span> into sequence that matches {org_lc(x["kind"])} '
                      f'(<span class="mono">{e(" and ".join(x["others"]))}</span>). ' + (tx_break(slug, x) or ("This position is a gap in the assembly (a run of N). The reads continue into nuclear sequence that is missing from the assembly, and this sequence happens to resemble an organelle gene." if x.get("at_assembly_gap") else "")) + "</li>")
        for x in rec:
            items += (f'<li>In {word(x["distinct_pairs"])} distinct pairs, one read lies within 1 kb of <span class="mono">{e(x["slug"])}:{num(x["start"])}</span> and the other in {org_lc(x["kind"])}. {tx_place(slug, x)}</li>')
        out += f'<p style="margin-top:12px"><strong>Junctions found in more than one molecule.</strong></p><ul>{items}</ul>'
    return out
def tx_note(slug, c):
    n = TEXT.get(slug, {}).get("candidate_notes", {}).get(c["read"])
    if n: return n
    if c.get("assembly_locus"): return "At a locus where the assembly itself contains organelle-like DNA."
    if c["how"] == "inside one alignment": return "Apart from the inserted segment, the read is collinear with the nuclear sequence."
    return ""
def same_joins(ch):
    """joins carried by two molecules that the other two lists do not already hold (37_pair_analysis.py: same_join)"""
    rb = ch.get("recurrent_breaks", []); rec = ch.get("recurrent", []); out = []
    for x in ch.get("same_join", []):
        if x["kind"] not in ORG_WORD: continue
        if any(y["slug"] == x["slug"] and abs(y["slug_break"] - x["slug_from"]) <= 1000 for y in rb): continue
        if any(y["slug"] == x["slug"] and y["start"] - 1000 <= x["slug_from"] < y["start"] + 2000 for y in rec): continue
        if any(y["slug"] == x["slug"] and abs(y["slug_from"] - x["slug_from"]) <= 1000 for y in out): continue          # the same place on a second reference sequence
        out.append(x)
    return out
def tx_join(slug, x): return TEXT.get(slug, {}).get("join_notes", {}).get(f'{x["slug"]}:{x["slug_from"]}', "")
def tx_break(slug, x): return TEXT.get(slug, {}).get("break_notes", {}).get(f'{x["slug"]}:{x["slug_break"]}', "")
def tx_place(slug, x): return TEXT.get(slug, {}).get("place_notes", {}).get(f'{x["slug"]}:{x["start"]}', "")

def transfer_page(sp, K, T):
    slug = sp["slug"]; tx = TEXT.get(slug, {}); ctx = context(sp, K, T); s = T["screen"]
    lede = fmt(tx.get("transfer_lede", "Algal DNA was searched for in the assembly and in the reads from which it was built."), ctx)
    answer = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get("transfer_answer", []))
    # ---- assembly
    whole = ""
    if s["whole"]:
        rows = "".join(f'<tr><td class="mono">{e(w["contig"])}</td><td class="num">{num(w["length"])}</td><td>{e(w["cls"])}</td><td>{taxon_html(" ".join(w["best_target_desc"].split()[:2]))} <span class="mono muted">{e(w["best_target"].split("|", 1)[1])}</span></td><td class="num">{100 * w["median_identity"]:.1f}%</td></tr>' for w in s["whole"][:40])
        whole = f'''<h3>Assembly sequences that consist entirely of organelle DNA</h3>
<div class="tablewrap"><table class="data"><thead><tr><th>Sequence</th><th class="num">Length (bp)</th><th>Type</th><th>Best match</th><th class="num">Identity</th></tr></thead><tbody>{rows}</tbody></table></div>'''
    loci = ""
    if s["loci"]:
        rows = "".join(f'<tr><td class="mono">{e(L["contig"])}:{num(L["start"])}–{num(L["end"])}</td><td class="num">{num(L["span"])}</td><td class="num">{100 * L["identity"]:.1f}%</td>'
                       f'<td>{taxon_html(" ".join(L["target_desc"].replace("UNVERIFIED: ", "").split()[:2]))} {e(L["cls"].replace("algal ", "").replace("plastid", "plastid genome"))}{": " + e(", ".join(g for g in L.get("genes", []) if len(g) < 14)[:60]) if L.get("genes") else ""}</td>'
                       f'<td>{locus_reading(L, tx)}</td></tr>' for L in sorted(s["loci"], key=lambda L: -L["identity"]))
        loci = f'''<h3>Organelle-like regions within other assembly sequences</h3>
<div class="tablewrap"><table class="data"><thead><tr><th>Position in the assembly</th><th class="num">Length (bp)</th><th class="num">Identity</th><th>Best match</th><th>Interpretation</th></tr></thead><tbody>{rows}</tbody></table></div>'''
    numt = ""
    if s["numt_like"]:
        rows = "".join(f'<tr><td class="mono">{e(N["contig"])}:{num(N["start"])}–{num(N["end"])}</td><td class="num">{num(N["end"] - N["start"])}</td><td class="num">{100 * N["identity"]:.1f}%</td><td class="mono">{e(N["target"].split("|", 1)[1])}:{num(N["tstart"])}–{num(N["tend"])}</td></tr>' for N in s["numt_like"][:30])
        numt = f'''<h3>Mitochondrial DNA in nuclear sequences</h3>
<p>{fmt(tx.get("numt_text", "These regions within nuclear sequences match the slug’s own mitochondrial genome. Nuclear copies of mitochondrial DNA are common in animals; they are listed here because the search that detects them would also detect plastid DNA."), ctx)}</p>
<div class="tablewrap"><table class="data"><thead><tr><th>Position in the assembly</th><th class="num">Length (bp)</th><th class="num">Identity</th><th>Region of the mitochondrial genome</th></tr></thead><tbody>{rows}</tbody></table></div>'''
    asm_text = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get("transfer_assembly", []))
    removed = s["removed"]
    asm_method = (f'The assembly was divided into 10 kb windows, and each window was compared with the same reference panel as the reads. '
                  f'In total, {num(removed["repeat_like_matches"])} matches were excluded because they fall in panel regions that are hit from 15 or more loci in the assembly; these are slug repeats that happen to resemble a short region of a panel sequence. '
                  f'A further {num(removed["low_complexity_matches"])} matches were to low-complexity sequence, and {num(removed["own_rdna_loci"]) if removed["own_rdna_loci"] else "none"} were to the slug’s own ribosomal DNA, which is 80 to 90% identical to algal ribosomal DNA.')
    # ---- reads
    reads_html = ""
    for R in T["reads"]:
        key = "reads_" + R["library"]
        para = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get(key, []))
        span = SPAN.get(slug); depth = (R.get("bases_screened") or R["bases"]) / span if span else None
        sens_ = R.get("sensitivity") or {}
        if R["kind"] == "short" and sens_.get("slug_fraction_of_bases"):          # short reads: depth of the slug DNA that can be placed
            depth = sens_["depth"]
            sub = (e(R["animal"][0].upper() + R["animal"][1:]) + f'. Nuclear DNA that can be placed on the assembly makes up {pct(sens_["slug_fraction_of_bases"], 0)} of the bases, which corresponds to about '
                   + (f"{depth:,.0f}" if depth >= 10 else f"{depth:,.1f}") + "-fold coverage of the genome.")
        else:
            sub = e(R["animal"][0].upper() + R["animal"][1:]) + "." + (f" The reads correspond to about {depth:,.0f}-fold coverage of the genome." if depth and depth >= 10 else (f" The reads correspond to about {depth:,.1f}-fold coverage of the genome." if depth else ""))
        reads_html += f'''<h3>{e(R.get("label") or (R["animal"][0].upper() + R["animal"][1:]))}: {gb(R.get("bases_screened") or R["bases"])} of {e(KIND[R["kind"]])}</h3>
<p class="small muted" style="margin:-4px 0 10px">{sub}</p>
{long_reads_block(R, slug) if R["kind"] != "short" else short_reads_block(R, slug)}
{para}'''
    reads_text = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get("transfer_reads", []))
    # places in the assembly tested with reads
    tests = [(R, R["locus_tests"]) for R in T["reads"] if R.get("locus_tests")]
    locus_html = ""
    if tests:
        names = list(dict.fromkeys(L["name"] for _, lt in tests for L in lt["loci"]))
        rows = ""
        for i, nm in enumerate(names):
            first = True
            here = [(R, next((L for L in lt["loci"] if L["name"] == nm), None)) for R, lt in tests]
            # long reads: only the reads of the places that were in the panel of the first pass were kept
            here = [(R, L) for R, L in here if L and (R["kind"] == "short" or L.get("reads_aligned"))]
            for R, L in here:
                def edge(k):
                    v = L["edges"][k]
                    if "note" in v: return '<td class="num muted" colspan="2">end of the sequence</td>'
                    thr = v.get("through_molecules", v["through"]) if R["kind"] != "short" else v["through"]
                    gp = v.get("through_with_gap_molecules", 0) if R["kind"] != "short" else v.get("through_with_gap", 0)
                    off = f'<br><span class="small muted">and {num(gp)} with a different junction nearby</span>' if gp else ""
                    return f'<td class="num"><strong>{num(thr)}</strong>{off}</td><td class="num">{num(v["stop_inside"] + v["stop_outside"])}</td>'
                span = len(here)
                head = (f'<td rowspan="{span}"><span class="mono">{e(L["contig"])}:{num(L["start"])}–{num(L["end"])}</span><br><span class="muted small">{"plastid-like" if L["kind"] == "plastid" else "mitochondrial"}, {num(L["end"] - L["start"])} bp</span></td>') if first else ""
                rows += f'<tr>{head}<td>{e(R.get("label") or (R["animal"].split(" (")[0].capitalize() + ", " + KIND[R["kind"]]))}</td>{edge("left")}{edge("right")}</tr>'
                first = False
        locus_html = f'''<h3>Read-based tests of the organelle-like regions in the assembly</h3>
<p>{fmt(tx.get("locus_text", "For every organelle-like region within an assembly sequence, reads were aligned to the region together with 3 kb of flanking sequence on each side. At each boundary, a read that spans the boundary supports the junction between the region and its flank, whereas a read whose alignment ends there continues into different sequence."), ctx)}</p>
<div class="tablewrap"><table class="data"><thead><tr><th rowspan="2">Region</th><th rowspan="2">Source of reads</th><th colspan="2" class="num">Left boundary</th><th colspan="2" class="num">Right boundary</th></tr>
<tr><th class="num">reads spanning</th><th class="num">reads ending</th><th class="num">reads spanning</th><th class="num">reads ending</th></tr></thead><tbody>{rows}</tbody></table></div>'''
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
            if c.get("whole") == "plastid" or "plastid-only contig" in c["classes"]: return "Consists entirely of plastid DNA (listed above)."
            if "plastid fragment (assembly artefact)" in c["classes"]: return "The plastid region described above."
            k = census_of(c)
            if k and k["microbe"] >= 3 and k["microbe"] >= 0.5 * k["n"]: return "Bacterial sequence."
            if k: return "Too few genes to classify."
            if c["loci"] >= 2 and not c["animal_loci_on_contig"]: return "Adjacent genes of bacterial type and no slug gene among the proteins searched: a bacterial sequence in the assembly."
            if c["loci"] >= 2: return "Several genes of bacterial type: a bacterial sequence, or a sequence joined to one."
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
        census_head = ('<th class="num">G+C</th><th class="num" title="Protein-coding regions found when the whole sequence is translated in all six reading frames">Coding regions</th>'
                       '<th class="num">Of bacterial type</th><th class="num" title="Fraction of the sequence occupied by these coding regions">Coding fraction</th>') if has_census else ""
        own = {"elysia-chlorotica": "<i>E. chlorotica</i>", "elysia-marginata": "<i>E. marginata</i>"}.get(slug)
        own_note = (f" The published gene annotation of this assembly was not used as evidence that a locus is of animal origin, because it also includes genes on bacterial and plastid sequences, and such genes would otherwise support their own classification." if own else "")
        nm = P["no_match"]; n_nm = ps.get("no database match", 0)
        nomatch_text = ("these are short or weak alignments whose translations match no known protein" if not nm["low_complexity"] else
                        f'{num(nm["low_complexity"])} of these are low-complexity alignments' + ("" if nm["low_complexity"] >= n_nm else ", and the remainder are short or weak alignments (25 to 40% identity, spliced together across several exons) whose translations match no known protein"))
        census_note = (f'<p>Each assembly sequence that contains a locus of non-animal type was then translated in full, in all six reading frames, and every protein-coding region found on it was compared with the same protein database. The table gives the number of coding regions, the number that are of bacterial type, and the fraction of the sequence that they occupy. Bacterial DNA consists mostly of coding sequence, whereas a slug sequence of this size would contain only one or two genes. The G+C content of the genome as a whole is {100 * P["genome_gc"]:.0f}%.</p>' if has_census and P.get("genome_gc") else "")
        prot = f'''<h2>Protein-level search for algal genes</h2>
<p>A comparison at the DNA level can miss a gene from an alga that has no close relative among the reference sequences. A set of {num(P["n_queries"])} proteins was therefore aligned directly to the assembly: the plastid-encoded proteins of a green alga, photosynthesis proteins of <i>Chlamydomonas</i> and <i>Arabidopsis</i> (including those encoded in the nucleus), proteins of siphonous green algae, and all available proteins of <i>Vaucheria litorea</i>. Each matching locus was translated and compared with animal, algal and bacterial proteins.</p>
<p>In total, {num(ps["loci"])} loci matched. Of these, {num(ps.get("animal gene", 0))} are slug genes, for which the best match is a protein from another animal (for example, the gene for mitochondrial ATP synthase is detected with the chloroplast ATP synthase as the query).{own_note} A further {num(ps.get("no database match", 0))} match no protein from any other organism; {nomatch_text}.</p>
{census_note}
{ptext}
{f'<div class="tablewrap"><table class="data"><thead><tr><th>Sequence</th><th class="num">Length (bp)</th>{census_head}<th>Query proteins that matched</th><th>Interpretation</th></tr></thead><tbody>{rows}</tbody></table></div>' if rows else ""}'''
    caveats = "".join(f"<li>{fmt(p, ctx)}</li>" for p in tx.get("transfer_caveats", []))
    # tables written by 45_export_tables.py
    dls = []
    for R in T["reads"]:
        tag = slug if R["library"] == "main" else f'{slug}_{R["library"]}'; lab = (" (" + lc(R.get("label") or R["animal"]) + ")") if len(T["reads"]) > 1 else ""
        for suffix, what in (("_insertion_shapes.tsv", "reads with an insertion-like structure"), ("_junctions.tsv.gz", "all junctions within reads"), ("_junction_pairs.tsv.gz", "all discordant pairs and split reads that join organelle DNA to nuclear DNA")):
            if os.path.exists(os.path.join(SITE, "downloads", tag + suffix)): dls.append(f'<a href="../downloads/{tag}{suffix}">{what}{e(lab)}</a>')
    dl_html = f'<p class="small" style="margin-top:14px">Tables available for download: {"; ".join(dls)}; and <a href="../downloads/{slug}_gene_transfer.json">all values shown on this page</a> (JSON).</p>' if dls else ""
    body = f"""
<section class="page-head"><div class="wrap">
<h1>Search for algal DNA in the {it(sp["species"]).replace(" ", "&nbsp;")} genome</h1>
<p class="lede">{lede}</p>
<div class="answer">{answer}</div>
</div></section>

<section class="block article"><div class="wrap">
<h2>Search of the assembly</h2>
<p>{asm_method}</p>
{asm_text}
{whole}
{loci}
{numt}

{"<h2>Search of the reads</h2>" + READS_FRAMING + reads_text + reads_html + dl_html if T["reads"] else ""}

{prot}

{"<h2>Limitations</h2><ul>" + caveats + "</ul>" if caveats else ""}
</div></section>
"""
    src = f"Data for {it(sp['species'])}: assembly {e(sp['asm_name'])} ({e(sp['accession'])}){'; reads ' + e(', '.join(r for R in T['reads'] for r in R['runs'])) if T['reads'] else ''}. The methods and data sources are described on the <a href=\"../methods.html\">methods page</a>."
    return chrome.page(f"Search for algal DNA in the {sp['species']} genome", body, root="../", species=sp["species"], tabs=tabs_for(slug), tab_current="gene-transfer.html", today=today, sources=src, scripts=TIP_JS,
                       description=f"A search for algal DNA in the genome of the sea slug {sp['species']}, based on its assembly and its sequencing reads.")

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
