#!/usr/bin/env python3
"""Collect, for every species except E. crispata, the results of the read screens, the assembly
screen and the protein screen into two files per species:
  results/genus/<slug>.kleptoplast.json   whose chloroplasts the sequenced animals carried
  results/genus/<slug>.transfer.json      organelle-like DNA in the assembly and in the reads

usage: 39_collect_species.py REPO READS_ROOT ASM_SCREEN_DIR [PROTSCREEN_DIR]
READS_ROOT/<slug>[.<library>]/ holds analysis.json (long reads, 34_read_analysis.py) or pairs.json
(short reads, 37_pair_analysis.py), barcodes/barcodes.json, kcp.contigs.tsv and locus_tests.json
where they exist. Libraries are listed in config/species_reads.tsv.
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
def load(p): return json.load(open(p)) if os.path.exists(p) else None
for sp in SP:
    slug = sp["slug"]
    if slug == "elysia-crispata": continue
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
            contigs = dict(n=len(rows), bp=sum(int(r["length"]) for r in rows), longest=max((int(r["length"]) for r in rows), default=0),
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
                             min_identity=a.get("min_identity"), locus_tests=load(os.path.join(d, "locus_tests.json"))))
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
