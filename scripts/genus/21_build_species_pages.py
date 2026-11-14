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

def per_cell(k):
    """words for a mean number of insertions per cell"""
    if k >= 1.5: return f"about {k:,.0f} insertions in every cell" if k >= 10 else f"about {k:.1f} insertions in every cell".replace(".0 ", " ")
    if k >= 0.75: return "about one insertion in every cell"
    return f"one insertion in every {round(1 / k):,} cells"

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
    for L in libs:          # short reads: the plastid contigs, sorted into kinds
        G = (L.get("contigs") or {}).get("groups")
        if not G: continue
        def grow(g):
            if g["kind"] == "close": what = f'Within 5% of {taxon_html(g["label"])} <span class="mono muted">{e(", ".join(r["id"].split("|", 1)[-1] for r in g["references"][:2]))}</span>'
            else: what = ("No close reference, " + ("plentiful" if "plentiful" in g["kind"] else "sparse") + '<br><span class="small muted">nearest: ' + ", ".join(taxon_html(x) for x in list(dict.fromkeys(" ".join(r["name"].split()[:2]) for r in g["references"]))[:3]) + "</span>")
            ident = f'{100 * g["identity_min"]:.1f}%' if g["n"] == 1 else f'{100 * g["identity_min"]:.1f} to {100 * g["identity_max"]:.1f}%'
            return (f'<tr><td>{what}</td><td class="num">{num(g["n"])}</td><td class="num">{g["bp"] / 1e3:,.1f} kb</td><td class="num">{ident}</td>'
                    f'<td class="num">{g["depth"]:,.0f}×</td><td class="num">{100 * g["relative"]:.0f}%</td></tr>' if g["relative"] >= 0.095 else
                    f'<tr><td>{what}</td><td class="num">{num(g["n"])}</td><td class="num">{g["bp"] / 1e3:,.1f} kb</td><td class="num">{ident}</td><td class="num">{g["depth"]:,.1f}×</td><td class="num">{100 * g["relative"]:.1f}%</td></tr>')
        contig_tables_one = f'''<div class="tablewrap"><table class="data">
<thead><tr><th>Kind of plastid contig</th><th class="num">Contigs</th><th class="num">Length together</th><th class="num">Identity to the nearest reference</th><th class="num">Depth</th><th class="num">Relative to the first row</th></tr></thead>
<tbody>{"".join(grow(g) for g in G)}</tbody></table></div>
<p class="small muted" style="max-width:86ch">The plastid reads were assembled into contigs, and the contigs sorted. A contig within 5% of a reference goes with that reference. The others have no close reference and are split by depth: a tenth or more of the depth of the first row, or less. Depth is how many times over the reads cover a contig, averaged over the contigs of a row; it shows how plentiful each plastid genome was in the animal. Only reads that resemble some reference were assembled, so of an alga with no close reference only the better conserved genes are here.</p>'''
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

READS_FRAMING = """<p>A slug takes chloroplasts into cells of its gut, not into its eggs or sperm. If a piece of chloroplast DNA were built into a chromosome, that would happen in one body cell, and only that cell and its daughters would carry it. In DNA from a whole animal such an insertion is carried by a few molecules, perhaps by one, and the other molecules from the same place do not have it. So a single read that joins chloroplast DNA to slug DNA is all the evidence there can be, and reads that disagree with it are expected.</p>
<p>The difficulty is that preparing DNA for sequencing also joins unrelated molecules to each other now and then, and in a single read such a join looks the same. Three things can be asked of the reads.</p>
<ul>
<li><strong>Shape.</strong> A fragment that sits in a chromosome has slug DNA on both sides, and the two sides are neighbours in the genome. A join made in the tube puts two molecules end to end. A read that goes slug DNA, organelle DNA, and then the same stretch of slug DNA carrying on is called an insertion shape here. Two accidental joins that happen to pick neighbours do not occur by chance, so this shape needs another explanation.</li>
<li><strong>Rate.</strong> Joining in the tube treats all DNA alike. Chloroplast DNA should then be joined to slug DNA as often as mitochondrial DNA is, and as often as one piece of slug DNA is joined to an unrelated one. Real insertions would come on top of that.</li>
<li><strong>The same join twice.</strong> A join found in two molecules that are not copies of each other is no accident of the library. It means a clone of cells, or an insertion that was inherited. An insertion in a single cell cannot pass this test and does not have to.</li>
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
        src = "every read of 2 kb or more" if allreads else f'a random sample, one read in {bg["depth"]["sample_one_in"]:,}' if bg.get("depth") and bg["depth"].get("sample_one_in") else "a random sample"
        e_ = bg["insertion_slug_elsewhere"]
        rows += (f'<tr><td>Pieces of slug DNA, for comparison <span class="muted small">({src})</span></td><td class="num">{num(bg["slug_pieces"])}</td>{rate_cell(bg["slug_pieces_joined"], bg["slug_pieces"])}'
                 f'<td class="num">{num(e_["unique"])}</td></tr>')
    note = ("The first three columns count every read of 2 kb or more in one and the same way, so that the three rates can be compared. " if allreads else "")
    note += (f'A read holds a kind of DNA when at least {S["min_segment"]} of its bases match it at {100 * S["min_identity"]:.0f}% or better. It is joined to slug DNA when a stretch of slug DNA of that length lies next to the piece in the same read. '
             f'Only joins to slug DNA that can be placed are counted (mapping quality 30 or more). The comparison row counts pieces of slug DNA, because a read in which two pieces of slug DNA are joined holds two of them: '
             f'the pieces are those of reads that can be placed, and a piece is joined when the read goes on with slug DNA from another sequence or more than 50 kb away (joins at the very end of an assembly sequence are left out). '
             f'That row takes in every real difference between the chromosomes of this animal and the assembly as well as joins made in the tube, so for the tube it is a ceiling. The last column counts insertion shapes: in the comparison row, those in which the piece in the middle comes from one clear place elsewhere in the genome.')
    out = f'''<div class="tablewrap"><table class="data"><thead><tr><th>Kind of DNA</th><th class="num">Reads that hold it (pieces, for slug DNA)</th><th class="num">Joined to slug DNA at one clear place</th><th class="num">Rate</th><th class="num">Inside slug DNA (insertion shape)</th></tr></thead><tbody>{rows}</tbody></table></div>
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
                where = f'none: the flanks are chloroplast DNA, at <span class="mono">{num(min(c["left_to"], c["right_to"]))}</span>' + (f' ({e(gn)})' if gn else "")
                what = f'slug DNA, <span class="mono">{e(c["piece_sequence"])}:{num(c["piece_from"])}–{num(c["piece_to"])}</span>'
            oc = c.get("other_reads_crossing")
            return (f'<tr><td class="mono">{e(c["read"])}</td><td>{where}</td><td>{what}</td><td class="num">{num(c["piece_bp"])}</td>'
                    f'<td class="num">{"–" if oc is None else num(oc)}</td><td>{e(tx_note(slug, c))}</td></tr>')
        out += f'''<p style="margin-top:14px"><strong>Reads with the shape of an insertion.</strong> Each row is one molecule.</p>
<div class="tablewrap"><table class="data"><thead><tr><th>Read</th><th>Place in the slug genome</th><th>Piece in the middle</th><th class="num">Its length, bp</th><th class="num" title="Other reads that cross the same place in one piece, without the insert">Other reads across the place, without it</th><th>Note</th></tr></thead>
<tbody>{"".join(crow(c) for c in ins_rows + mir)}</tbody></table></div>'''
    d = (bg.get("depth") or {}).get("spanning_depth") if bg else None
    if d:
        s500 = d.get("500") or 0
        if s500 > 0:
            k = 6 / s500
            out += (f'<p class="small muted" style="margin-top:8px;max-width:92ch">What this library could show: reads long enough to hold a 500-base insertion with {S["min_segment"]} bases of slug DNA on each side cover each place in the genome {s500:,.1f} times. '
                    f'Each read comes from one cell, so in such reads the library holds about {(f"{s500 / 2:,.0f}" if s500 >= 6 else f"{s500 / 2:,.1f}")} cells’ worth of DNA. An insertion carried by every cell on one of its two chromosome copies would be met by about that many reads, '
                    f'and insertions scattered over different cells add up in the same way: with {per_cell(k)}, averaged over all the cells that went into the tube, about three insertion-shaped reads would be expected. '
                    f'{copies_note(bg, R)}{TEXT.get(slug, {}).get("could_show_" + R["library"], "")}</p>')
    return out
def copies_note(bg, R):
    """amplified libraries: how many reads are copies of another read, estimated on the random sample"""
    cps = (bg.get("depth") or {}).get("copies") or {}
    if "amplified" not in R["platform"] or not cps.get("placed_reads"): return ""
    if cps["reads_with_a_copy_in_the_sample"]:
        return (f'The library was amplified, and about one read in {round(1 / cps["copies_per_read"])} is a copy of another ({cps["reads_with_a_copy_in_the_sample"]} of {num(cps["placed_reads"])} sample reads have their copy in the sample), so the figures are a little too high. ')
    return (f'The library was amplified, which makes copies of molecules. None of {num(cps["placed_reads"])} sample reads has its copy in the sample, which puts copies at under {100 * cps["copies_per_read_upper"]:.0f}% of the reads. ')
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
        rows += (f'<tr><td>Slug DNA, for comparison <span class="muted small">(a random sample of the library)</span></td><td class="num">{num(ns)}</td>{rate_cell(c.get("slug_pairs_apart", 0), ns)}'
                 f'{rate_cell(smp.get("split_mates_joining_unrelated_slug", 0), ns)}</tr>')
    out = f'''<div class="tablewrap"><table class="data"><thead><tr><th>Kind of DNA</th><th class="num">Pairs that hold it</th><th class="num">Other mate in slug DNA</th><th class="num">Rate</th><th class="num">Join read through in one mate</th><th class="num">Rate</th></tr></thead><tbody>{rows}</tbody></table></div>
<p class="small muted" style="margin-top:8px;max-width:92ch">A pair holds a kind of DNA when 80 bases or more of a mate match it at 97% or better. In the third column the other mate lies at one clear place in slug DNA; for slug DNA itself, the two mates lie on different sequences or more than 20 kb apart. In the fifth column a single mate carries at least 40 bases of each kind, so the join itself is read. Every such pair and mate {"that was examined " if R.get("thinned", 1) > 1 else ""}is in the junction table for download.</p>'''
    if R.get("thinned", 1) > 1:
        ex = (rates.get("cp") or {}).get("examined") or {}
        out += (f'<p class="small muted" style="margin-top:8px;max-width:92ch">So much of this library is organelle DNA that the pairs matching an organelle genome were thinned: every read was mapped and counted, and one pair in {R["thinned"]} of those that match was examined further, chosen by its number in the run. '
                f'The first two rows are therefore estimates, {R["thinned"]} times what was counted: {num(ex.get("pairs_with_mate", 0))} chloroplast pairs were examined, of which {num(ex.get("pairs_also_with_slug", 0))} have their other mate in slug DNA and {num(ex.get("split_mates", 0))} read through a join. '
                f'Pairs that match one of the places under test were all kept.</p>')
    if sens:
        sp_ = sens.get("split_mates_per_insertion_per_cell"); mp = sens.get("mixed_pairs_per_insertion_per_cell")
        if sp_:
            cpv = rates.get("cp", {}); n_sp = cpv.get("split_mates", 0)
            out += (f'<p class="small muted" style="margin-top:8px;max-width:92ch">What this library could show: the fragments are about {num(sens["insert_median"])} bases long and each mate reads {sens["read_length"]:.0f} of them, so a join inside a fragment is usually read through. '
                    f'One insertion in every cell, each cell with a different one, would add about {sp_:,.0f} mates that read through a join' + (f' and about {mp:,.0f} pairs that straddle one' if mp and mp >= 1 else "") +
                    ('. The library holds no mate that reads from chloroplast DNA into slug DNA.</p>' if not n_sp else
                     f'. The library holds {"an estimated " if R.get("thinned", 1) > 1 else ""}{num(n_sp)} such mates for chloroplast DNA already, at the rate at which it joins any two pieces of DNA.</p>') if sens.get("insert_median") else "")
    rb = [x for x in ch.get("recurrent_breaks", []) if x["kind"] in ORG_WORD]; rec = [x for x in ch.get("recurrent", []) if x["kind"] in ORG_WORD]
    sj = same_joins(ch)
    if rb or rec or sj:
        items = ""
        for x in sj:
            thru = [v for v in x.get("evidence", []) if v["evidence"] == "split mate"]
            pos = num(thru[0]["other_pos"]) if thru else (f'{num(x["other_from"])}–{num(x["other_to"])}' if x["other_to"] - x["other_from"] > 5 else num(x["other_from"]))
            items += (f'<li>{word(x["molecules"]).capitalize()} molecules that are not copies of each other join <span class="mono">{e(x["slug"])}:{num(thru[0]["slug_pos"] if thru else x["slug_from"])}</span> to {org_lc(x["kind"])} '
                      f'(<span class="mono">{e(x["other_t"].split("|", 1)[1])}:{pos}</span>). {tx_join(slug, x)}</li>')
        # breaks at one base of the slug sequence whose other sides lie at several places of a reference (the copies of a repeated gene) are one item
        merged = {}
        for x in rb:
            k_ = (x["kind"], x["slug"], x["slug_break"])
            if k_ in merged: merged[k_]["molecules"] += x["molecules"]; merged[k_]["others"].append(f'{x["other_t"].split("|", 1)[1]}:{num(x["other_break"])}')
            else: merged[k_] = dict(x, others=[f'{x["other_t"].split("|", 1)[1]}:{num(x["other_break"])}'])
        for x in merged.values():
            items += (f'<li>{word(x["molecules"]).capitalize()} molecules read from <span class="mono">{e(x["slug"])}:{num(x["slug_break"])}</span> into what matches {org_lc(x["kind"])} '
                      f'(<span class="mono">{e(" and ".join(x["others"]))}</span>). ' + (tx_break(slug, x) or ("The place is a gap in the assembly, a run of N: the reads carry on into slug DNA that the assembly lacks, and that DNA happens to resemble an organelle gene." if x.get("at_assembly_gap") else "")) + "</li>")
        for x in rec:
            items += (f'<li>{word(x["distinct_pairs"]).capitalize()} different pairs have one mate within 1 kb of <span class="mono">{e(x["slug"])}:{num(x["start"])}</span> and the other in {org_lc(x["kind"])}. {tx_place(slug, x)}</li>')
        out += f'<p style="margin-top:12px"><strong>The same join in more than one molecule.</strong></p><ul>{items}</ul>'
    return out
def tx_note(slug, c):
    n = TEXT.get(slug, {}).get("candidate_notes", {}).get(c["read"])
    if n: return n
    if c.get("assembly_locus"): return "At a place where the assembly itself holds organelle-like DNA."
    if c["how"] == "inside one alignment": return "The read follows the slug sequence without a break; the piece is extra."
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
        key = "reads_" + R["library"]
        para = "".join(f"<p>{fmt(p, ctx)}</p>" for p in tx.get(key, []))
        span = SPAN.get(slug); depth = (R.get("bases_screened") or R["bases"]) / span if span else None
        sens_ = R.get("sensitivity") or {}
        if R["kind"] == "short" and sens_.get("slug_fraction_of_bases"):          # short reads: depth of the slug DNA that can be placed
            depth = sens_["depth"]
            sub = (e(R["animal"][0].upper() + R["animal"][1:]) + f'. {pct(sens_["slug_fraction_of_bases"], 0)} of the bases are slug DNA that can be placed on the assembly, which covers each place in the genome about '
                   + (f"{depth:,.0f}" if depth >= 10 else f"{depth:,.1f}") + " times.")
        else:
            sub = e(R["animal"][0].upper() + R["animal"][1:]) + "." + (f" The reads cover each place in the genome about {depth:,.0f} times." if depth and depth >= 10 else (f" The reads cover each place in the genome about {depth:,.1f} times." if depth else ""))
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
                    off = f'<br><span class="small muted">and {num(gp)} with a join of {"its" if gp == 1 else "their"} own nearby</span>' if gp else ""
                    return f'<td class="num"><strong>{num(thr)}</strong>{off}</td><td class="num">{num(v["stop_inside"] + v["stop_outside"])}</td>'
                span = len(here)
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
    # tables written by 45_export_tables.py
    dls = []
    for R in T["reads"]:
        tag = slug if R["library"] == "main" else f'{slug}_{R["library"]}'; lab = (" (" + (R.get("label") or R["animal"]) + ")") if len(T["reads"]) > 1 else ""
        for suffix, what in (("_insertion_shapes.tsv", "reads with the shape of an insertion"), ("_junctions.tsv.gz", "every join inside a read"), ("_junction_pairs.tsv.gz", "every pair and mate that joins organelle DNA to slug DNA")):
            if os.path.exists(os.path.join(SITE, "downloads", tag + suffix)): dls.append(f'<a href="../downloads/{tag}{suffix}">{what}{e(lab)}</a>')
    dl_html = f'<p class="small" style="margin-top:14px">Tables for download: {"; ".join(dls)}; and <a href="../downloads/{slug}_gene_transfer.json">all the numbers on this page</a>.</p>' if dls else ""
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

{"<h2>In the reads</h2>" + READS_FRAMING + reads_text + reads_html + dl_html if T["reads"] else ""}

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
