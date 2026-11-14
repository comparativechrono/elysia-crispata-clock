#!/usr/bin/env python3
"""Collect, for every species except E. crispata, the results of the read screens, the assembly
screen and the protein screen into two files per species:
  results/genus/<slug>.kleptoplast.json   whose chloroplasts the sequenced animals carried
  results/genus/<slug>.transfer.json      organelle-like DNA in the assembly and in the reads

usage: 39_collect_species.py REPO READS_ROOT ASM_SCREEN_DIR [PROTSCREEN_DIR]
READS_ROOT/<slug>[.<library>]/ holds analysis.json (long reads, 34_read_analysis.py) or pairs.json
(short reads, 37_pair_analysis.py), barcodes/barcodes.json, kcp.contigs.tsv and locus_tests.json
where they exist, and for long reads candidates.json (46_candidates.py: the shape of every read
that joins organelle DNA to slug DNA), background.json (the same on the random sample of the
library) or allreads.json and control.json (47_library_control.sh, the same on every read), and
recurrent.json (51_recurrent_junctions.py). Libraries are listed in config/species_reads.tsv.
"""
import sys, os, json, csv, re, ssl, time, urllib.request, collections
repo, reads_root, screen_dir = sys.argv[1:4]
prot_dir = sys.argv[4] if len(sys.argv) > 4 else None
SP = list(csv.DictReader(open(os.path.join(repo, "config", "species.tsv")), delimiter="\t"))
LIB = collections.defaultdict(list)
for r in csv.DictReader(open(os.path.join(repo, "config", "species_reads.tsv")), delimiter="\t"): LIB[r["slug"]].append(r)
OUT = os.path.join(repo, "results", "genus"); os.makedirs(OUT, exist_ok=True)
CACHE = os.path.join(screen_dir, "feature_tables"); os.makedirs(CACHE, exist_ok=True)
ca = os.environ.get("SSL_CERT_FILE"); ctx = ssl.create_default_context(cafile=ca) if ca else ssl.create_default_context()

REF_NAME = {}
_pt = os.path.join(OUT, "panel.tsv.gz")
if os.path.exists(_pt):
    import gzip
    for _l in gzip.open(_pt, "rt"):
        _f = _l.rstrip("\n").split("\t")
        if len(_f) > 4: REF_NAME[_f[0]] = " ".join(re.sub(r" (chloroplast|plastid|genome assembly|DNA|strain|isolate|voucher|culture)\b.*|,.*", "", _f[4]).split()[:4])

def contig_groups(rows):
    """Plastid contigs of a short-read library sorted into kinds: contigs within 5% of a reference
    go with that reference; the others are split by how deeply they are covered, compared with the
    most plentiful kind (a tenth of its depth or more, or less). Depth is the assembler's own figure,
    averaged over the length of the contigs."""
    R = [dict(length=int(r["length"]), depth=float(r["depth"]), ref=r["best_panel_match"], identity=float(r["identity"])) for r in rows]
    if not R: return []
    def summary(ms, kind, label):
        bp = sum(m["length"] for m in ms); refs = collections.Counter()
        for m in ms: refs[m["ref"]] += m["length"]
        return dict(kind=kind, label=label, n=len(ms), bp=bp, depth=sum(m["depth"] * m["length"] for m in ms) / bp, identity_min=min(m["identity"] for m in ms), identity_max=max(m["identity"] for m in ms),
                    depth_min=min(m["depth"] for m in ms), depth_max=max(m["depth"] for m in ms), references=[dict(id=k, name=REF_NAME.get(k, k), bp=v) for k, v in refs.most_common(4)])
    close = collections.defaultdict(list); far = []
    for m in R: (close[REF_NAME.get(m["ref"], m["ref"])] if m["identity"] >= 0.95 else far).append(m)
    G = [summary(v, "close", k) for k, v in close.items()]
    top = max((g["depth"] for g in G), default=max(m["depth"] for m in R))
    for kind, ms in (("distant, plentiful", [m for m in far if m["depth"] >= 0.1 * top]), ("distant, sparse", [m for m in far if m["depth"] < 0.1 * top])):
        if ms: G.append(summary(ms, kind, ""))
    for g in G: g["relative"] = g["depth"] / top
    return sorted(G, key=lambda g: -g["depth"])

def features(acc):
    """gene and product names with their coordinates from the GenBank feature table of acc"""
    p = os.path.join(CACHE, acc + ".ft")
    if not os.path.exists(p):
        try:
            d = urllib.request.urlopen(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=ft&retmode=text", timeout=60, context=ctx).read().decode()
            time.sleep(0.4)
        except Exception: d = ""
        open(p, "w").write(d)
    feats = []; cur = None
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3 and f[0] and f[2]:
            a = int(re.sub(r"[<>^]", "", f[0])); b = int(re.sub(r"[<>^]", "", f[1]))
            cur = dict(start=min(a, b), end=max(a, b), kind=f[2], name=""); feats.append(cur)
        elif len(f) >= 5 and cur is not None and f[3] in ("gene", "product") and (not cur["name"] or f[3] == "gene"):
            cur["name"] = f[4]
    return [x for x in feats if x["kind"] in ("gene", "rRNA", "tRNA", "CDS") and x["name"]]
def genes_at(target, a, b):
    acc = target.split("|", 1)[1]
    names = []
    for x in features(acc):
        if x["start"] <= b and x["end"] >= a and x["end"] - x["start"] < 20000 and x["name"] not in names: names.append(x["name"])
    return names[:6]

def libdir(slug, lib): return os.path.join(reads_root, slug if lib == "main" else f"{slug}.{lib}")
def trim_candidate(c):
    """one read with the shape of an insertion (or its mirror image), for the page and the download"""
    mid = c["middle"]; L, R = c["left"], c["right"]
    org_t = mid[0]["t"]; a = min(m["ts"] for m in mid if m["t"] == org_t); b = max(m["te"] for m in mid if m["t"] == org_t)
    ids = [m.get("ident") for m in mid if m.get("ident") is not None]
    row = dict(read=c["read"], read_length=c["read_length"], shape=c["shape"], insert=c["insert"], how=c["how"],
               flank_sequence=L["t"], left_from=L["ts"], left_to=L["te"], right_from=R["ts"], right_to=R["te"], flank_strand=L["strand"], flank_mapq=min(L.get("mapq", 0), R.get("mapq", 0)),
               flank_identity=round(min(L["ident"], R["ident"]), 4), reference_gap=c["reference_gap"], read_gap=c["read_gap"],
               piece_sequence=org_t.split("|", 1)[-1], piece_from=a, piece_to=b, piece_bp=sum(m["qe"] - m["qs"] for m in mid), piece_pieces=len(mid),
               piece_identity=round(min(ids), 4) if ids else None, place=c.get("place"), other_reads_crossing=c.get("other_reads_crossing"),
               assembly_locus=c.get("assembly_locus"), second_pass=c.get("second_pass"), folds_back=c.get("folds_back_after", False))
    if c["shape"] == "insertion" and "|" in org_t: row["piece_genes"] = genes_at(org_t, a, b)
    elif c["shape"] == "mirror" and "|" in L["t"]: row["flank_genes"] = genes_at(L["t"], min(L["ts"], R["ts"]), max(L["te"], R["te"]))
    return row
def long_read_shapes(d, mode):
    """shapes of the reads that join organelle DNA to slug DNA, with what the library makes by itself"""
    C = load(os.path.join(d, "candidates.json"))
    if not C: return None
    out = dict(min_identity=C["min_identity"], min_segment=C["min_segment"], min_short_piece=C["min_short_piece"], reads_with_piece=C["reads_with_piece"],
               reads_joining_slug=C["reads_joining_slug"], reads_joined_to_placed_slug=C.get("reads_joined_to_placed_slug", {}), shapes=C["shapes"], shared_with_assembly=C.get("shared_with_assembly", {}), second_pass=C.get("second_pass", {}),
               depth=C.get("depth"))
    rows = {c["read"] + "|" + c["shape"] + "|" + c["insert"]: c for c in C["candidates"]}
    bg = None
    A = load(os.path.join(d, "allreads.json")); K = load(os.path.join(d, "control.json")); B = load(os.path.join(d, "background.json"))
    def rates(X):
        # slug DNA: reads that lie at one clear place, and those of them joined to a piece at another clear place;
        # organelle DNA: reads that hold it, and those of them joined to slug DNA at one clear place
        n = X["nuclear"]; rw = X["reads_with_piece"]; rj = X.get("reads_joined_to_placed_slug", {})
        # a read of slug DNA joined to a second piece of slug DNA holds two pieces, either of which could have been the one
        # that was joined; counted per piece, the figure compares with the share of organelle reads joined to slug DNA
        J = n.get("reads_placed_joined_to_placed_elsewhere", 0); P = n.get("reads_slug_placed", 0)
        return dict(reads_slug_only=n.get("reads_slug_only", 0), reads_slug_placed=P, joining_unrelated=J, slug_pieces=P + J, slug_pieces_joined=2 * J,
                    joining_unrelated_any=n.get("reads_joining_unrelated_slug_pieces", 0), folding_back=n.get("reads_folding_back", 0),
                    cp_reads=rw.get("cp", 0), cp_joined=rj.get("cp", 0), mt_reads=rw.get("slugmt", 0), mt_joined=rj.get("slugmt", 0))
    if A and K:          # every read of the library
        bg = dict(source="all reads", reads=A["reads_examined"], **rates(A), depth=A.get("depth"),
                  insertion_slug_elsewhere={k: K["insertion_slug_elsewhere"][k] for k in ("n", "unique", "unique_alone_at_place")}, shapes=K["shapes"], second_pass=K.get("second_pass", {}))
        for c in K["candidates"]: rows.setdefault(c["read"] + "|" + c["shape"] + "|" + c["insert"], c)
    elif B:              # a random sample
        bg = dict(source="sample", reads=B["reads_examined"], **rates(B), depth=B.get("depth"),
                  insertion_slug_elsewhere={k: B["insertion_slug_elsewhere"][k] for k in ("n", "unique", "unique_alone_at_place")}, shapes=B["shapes"], second_pass=B.get("second_pass", {}))
        for c in B["candidates"]: rows.setdefault(c["read"] + "|" + c["shape"] + "|" + c["insert"], c)
    out["background"] = bg
    out["candidates"] = sorted((trim_candidate(c) for c in rows.values()), key=lambda r: (r["shape"] != "insertion", r["insert"], -r["piece_bp"]))
    R = load(os.path.join(d, "recurrent.json"))
    if R: out["recurrent_joins"] = R
    return out
def load(p): return json.load(open(p)) if os.path.exists(p) else None
for sp in SP:
    slug = sp["slug"]
    if slug == "elysia-crispata":
        # E. crispata has its own pipeline (scripts/10 to 29). Two steps of this one were applied to it as well
        # (run_screens.sh crispata): the search for short organelle pieces inside reads, and the library's own
        # rate of joining, measured on a sample. They are collected here for scripts/31_build_pages.py.
        d = libdir(slug, "main"); sh = long_read_shapes(d, "hifi")
        if sh:
            tn = os.path.join(d, "targeted.names")
            sh["targeted_reads"] = sum(1 for l in open(tn) if l.strip()) if os.path.exists(tn) else None
            # every read with a join between chloroplast DNA and slug DNA, fetched again and compared with the others (run_screens.sh crispata)
            jd = libdir(slug, "joins"); rj = load(os.path.join(jd, "recurrent.json")); jn = os.path.join(jd, "joins.names")
            if rj:
                sh["joins_compared"] = dict(reads=sum(1 for l in open(jn) if l.strip()) if os.path.exists(jn) else None, joins=rj["joins_examined"], groups_of_copies=rj["groups_same_well_only"],
                                            recurrent=rj["recurrent"], near=rj.get("near", []), same_join_means=rj.get("same_join_means"))
            json.dump(sh, open(os.path.join(OUT, slug + ".shapes.json"), "w"), indent=1)
            print(slug, "shapes of the targeted reads and the sample collected")
        continue
    libs = []; tr_reads = []
    for L in LIB.get(slug, []):
        d = libdir(slug, L["library"])
        a = load(os.path.join(d, "analysis.json")) or load(os.path.join(d, "pairs.json"))
        if not a: continue
        cl = a["classes"]
        def part(k): return dict(reads=cl.get(k, {}).get("reads", cl.get(k, {}).get("pairs")), bp=cl.get(k, {}).get("bp", 0), fraction=cl.get(k, {}).get("fraction_of_bases", 0.0))
        contigs = None; ct = os.path.join(d, "kcp.contigs.tsv")
        if os.path.exists(ct):
            rows = list(csv.DictReader(open(ct), delimiter="\t"))
            contigs = dict(n=len(rows), bp=sum(int(r["length"]) for r in rows), longest=max((int(r["length"]) for r in rows), default=0), groups=contig_groups(rows),
                           top=[dict(contig=r["contig"], length=int(r["length"]), depth=float(r["depth"]), best_panel_match=r["best_panel_match"], identity=float(r["identity"])) for r in sorted(rows, key=lambda r: -int(r["length"]))[:8]])
        # a library of which only a sample of the matching pairs was analysed keeps the size of what was screened
        full = [os.path.join(d, r + ".count.full") for r in L["runs"].split(",")]
        screened = sum(int(open(f).read().split()[1]) for f in full if os.path.exists(f)) or a["bases"]
        lib = dict(library=L["library"], label=L.get("label", ""), runs=L["runs"].split(","), kind=L["kind"], platform=L["platform"], animal=L["animal"], note=L["note"], reads=a["reads"], bases=a["bases"], bases_screened=screened,
                   plastid=part("cp"), mitochondrion=part("slugmt"), algal_mitochondrion=part("algmt"), algal_rdna=part("algrDNA"),
                   plastid_any_identity=a.get("plastid_any_identity"), other_taxa=a.get("other_taxa", {}), donors=a["donors"], barcodes=load(os.path.join(d, "barcodes", "barcodes.json")) or [], contigs=contigs,
                   hot_windows=[h for r in a["runs"] for h in r.get("hot_windows", [])][:8])
        libs.append(lib)
        tr_reads.append(dict(library=L["library"], label=L.get("label", ""), runs=lib["runs"], kind=L["kind"], platform=L["platform"], animal=L["animal"], reads=a["reads"], bases=a["bases"], bases_screened=screened, chimera=a["chimera"],
                             min_identity=a.get("min_identity"), locus_tests=load(os.path.join(d, "locus_tests.json")),
                             shapes=long_read_shapes(d, L["kind"]) if L["kind"] != "short" else None, sample=a.get("sample"), sensitivity=a.get("sensitivity"), note=L["note"],
                             thinned=a.get("organelle_pairs_examined_one_in", 1), candidate_places=load(os.path.join(d, "candidate_places.json")), repeat_joins=load(os.path.join(d, "satellite_joins.json")),
                             end_piece_control=load(os.path.join(d, "end_piece_control.json"))))
    # the places of the insertion-shaped long reads, looked up in the short-read libraries of the species (other animals):
    # reads that cross the place in slug sequence (41_locus_test.py on the baits of the places), and pairs or mates
    # of that library that join organelle DNA to slug DNA there (its junction table)
    cands = [c for R in tr_reads if R.get("shapes") for c in R["shapes"]["candidates"] if c["shape"] == "insertion" and c.get("place")]
    for R in tr_reads:
        cp_ = R.get("candidate_places"); jt = os.path.join(libdir(slug, R["library"]), "junctions.tsv.gz")
        if R["kind"] != "short" or not cp_ or not cands or not os.path.exists(jt): continue
        import gzip
        J = [l.rstrip("\n").split("\t") for l in gzip.open(jt, "rt")][1:]
        rows = []
        for c in cands:
            t, pos = c["place"]["sequence"], c["place"]["position"]
            L = next((x for x in cp_["loci"] if x["contig"] == t and abs(x["start"] - pos) <= 150), None)
            thr = None if not L else max((v.get("through", 0) for v in L["edges"].values() if "through" in v), default=0)
            near = [j for j in J if j[3] == t and abs(int(j[4]) - pos) <= (1000 if j[1] == "mixed pair" else 30)]
            rows.append(dict(read=c["read"], insert=c["insert"], sequence=t, position=pos, reads_crossing=thr,
                             joins_here=len(near), joins_here_same_organelle=sum(1 for j in near if j[2] == c["insert"]), joins=[dict(molecule=j[0], evidence=j[1], organelle=j[2], organelle_sequence=j[7], organelle_position=int(j[8])) for j in near[:10]]))
        R["candidate_places"] = dict(places=rows)
    if libs:
        main = libs[0]
        head = [dict(name=d["name"], fraction=d["fraction"], median_identity=d["median_identity"]) for d in main["donors"] if d["fraction"] >= 0.01][:3]
        json.dump(dict(slug=slug, species=sp["species"], libraries=libs, donors=head, plastid_fraction=main["plastid"]["fraction"]), open(os.path.join(OUT, slug + ".kleptoplast.json"), "w"), indent=1)
    scr = load(os.path.join(screen_dir, slug + ".screen.json"))
    if scr:
        for L in scr["loci"]: L["genes"] = genes_at(L["target"], L["tstart"], L["tend"])
        for w in scr["whole"]: w.pop("covered", None)
        prot = None
        if prot_dir:
            p = load(os.path.join(prot_dir, "protein_screen_raw.json"))
            fs = load(os.path.join(prot_dir, "foreign_sequences.json")) or {}     # 43_foreign_sequences.py
            census = {r["contig"]: r for r in fs.get("sequences", {}).get(slug, [])}
            if p and slug in p["by_assembly"]:
                loci = [l for l in p["loci"] if l["assembly"] == slug]
                keep = [l for l in loci if l["class"] not in ("animal gene",)]
                whole_names = {w["contig"]: w["cls"] for w in scr["whole"]}
                by = collections.defaultdict(list)
                for l in keep:
                    if l["class"] != "no database match": by[l["contig"]].append(l)
                contigs = []
                for c, ls in by.items():
                    animal_here = sum(1 for l in loci if l["contig"] == c and l["class"] == "animal gene")
                    genes = []
                    for l in ls:
                        g = l["query"].split("|", 1)[1] if l["query"].startswith("UDOFL|") else ((l["query_description"] or l["query"]).split()[0])
                        g = g.split("_")[0]
                        if g not in genes: genes.append(g)
                    row = dict(contig=c, length=ls[0]["contig_length"], loci=len(ls), classes=sorted({l["class"] for l in ls}), genes=genes, animal_loci_on_contig=animal_here,
                               identity_min=round(min(l["identity"] for l in ls), 2), identity_max=round(max(l["identity"] for l in ls), 2), with_stops=sum(1 for l in ls if l["stops"]),
                               multi_exon=sum(1 for l in ls if l["exons"] > 1), whole=whole_names.get(c),
                               # places that the species' own published annotation lists as genes (its protein is the best match of all)
                               listed_in_own_annotation=sum(1 for l in ls if l.get("own_annotation_bits", 0) > max(l["alga_bits"], l["microbe_bits"], l.get("met_bits", 0))))
                    k = census.get(c)
                    if k: row["census"] = dict(gc=k["gc"], examined_length=k["examined_length"], n_bases=k["n_bases"], stretches=k["stretches"], coding_bp=k["coding_bp"],
                                               nearest_bacteria=k["nearest_bacteria"], bacterial_identity_median=k["bacterial_identity_median"])
                    contigs.append(row)
                nm = [l for l in loci if l["class"] == "no database match"]
                prot = dict(n_queries=p["n_queries"], summary=p["by_assembly"][slug], genome_gc=fs.get("genome_gc", {}).get(slug),
                            non_animal_contigs=sorted(contigs, key=lambda x: (-x["length"] if "census" in x else 0, -x["loci"])),
                            no_match=dict(n=len(nm), low_complexity=sum(1 for l in nm if l["peptide_entropy"] < 3.3), short=sum(1 for l in nm if l["peptide_entropy"] >= 3.3 and l["aligned_aa"] < 120),
                                          other=[dict(contig=l["contig"], start=l["start"], end=l["end"], query=(l["query_description"] or l["query"])[:50], identity=round(l["identity"], 2), aligned_aa=l["aligned_aa"]) for l in nm if l["peptide_entropy"] >= 3.3 and l["aligned_aa"] >= 120]),
                            non_animal=[{k: l.get(k) for k in ("contig", "contig_length", "start", "end", "class", "query", "query_description", "identity", "query_cover", "aligned_aa", "frameshifts", "stops", "animal_bits", "met_bits", "other_sacoglossan_bits", "own_annotation_bits", "alga_bits", "microbe_bits", "best_alga_hit", "peptide_entropy", "exons")} for l in keep])
        json.dump(dict(slug=slug, species=sp["species"], assembly=dict(name=sp["asm_name"], accession=sp["accession"], level=sp["level"]), screen=scr, reads=tr_reads, protein_screen=prot),
                  open(os.path.join(OUT, slug + ".transfer.json"), "w"), indent=1)
    print(slug, "libraries:", [l["library"] for l in libs], "| screen:", bool(scr), "| loci:", len(scr["loci"]) if scr else None)
