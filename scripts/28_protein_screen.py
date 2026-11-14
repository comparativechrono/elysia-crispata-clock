#!/usr/bin/env python3
"""Step 2 of the protein screen: merge miniprot hits into loci and classify each locus.
usage: 28_protein_screen.py workdir alien_db.dmnd embedded_plastid.json queries.fa name=fasta ...

classes
  plastid fragment (assembly artefact)  the locus lies on one of the plastid fragments of script 24
  plastid-only contig / bacterial contig  the locus lies on a small contig that the DNA screen
                    (script 23, assembly_screen.json beside embedded_plastid.json) already found
                    to be plastid or bacterial DNA from end to end
  animal gene       the slug's own homologue. The best match among animals other than sea slugs
                    (MET) is at least as good as the best algal or microbial match; or it is within
                    10% of it and a protein of another sacoglossan species (SAC) beats both, as
                    happens for very conserved proteins such as GAPDH
  alga-like         otherwise, when the best non-animal match is an algal, plant or cyanobacterial protein
  microbe-like      otherwise, when it is a bacterial, archaeal or fungal protein
                    The published proteome of the species being screened is never counted as
                    animal evidence (OWN_PROTEOME below): gene annotations include genes on
                    bacterial and plastid sequences that were assembled along with the slug, and
                    such a gene would otherwise vouch for itself. Whether the locus is listed as a
                    gene in that annotation is recorded (own_annotation_bits).
  no database match the translated locus matches no protein at all (E <= 1e-3); these are short or
                    low-complexity alignments (Shannon entropy of the peptide is reported)
Writes protein_screen_raw.json (all loci) and, merging the follow-up checks of script 29 when
<assembly>_foreign_contigs.json exists, protein_screen.json for the site.
"""
import sys, os, json, subprocess, collections, re, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from miniprot_parse import parse
wd, db, emb, qfa = sys.argv[1:5]
asms = [a.split("=", 1)[0] for a in sys.argv[5:]]
D = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/diamond"
qdesc = {}
for line in open(qfa):
    if line[0] == ">":
        h = line[1:].rstrip(); qdesc[h.split()[0]] = h.split(" ", 1)[1] if " " in h else ""
frags = {a: [(r["contig"], r["start"], r["end"]) for r in rows] for a, rows in json.load(open(emb)).items()}
whole = {}
scr = os.path.join(os.path.dirname(os.path.abspath(emb)), "assembly_screen.json")
if os.path.exists(scr):
    for a, d in json.load(open(scr)).items():
        for x in d.get("whole_list", []): whole[(a, x["contig"])] = x["cls"]
loci = []
for asm in asms:
    hits = [h for hs in parse(f"{wd}/{asm}.gff").values() for h in hs if h["ident"] >= 0.25 and (h["qend"] - h["qstart"]) >= 40]
    hits.sort(key=lambda h: (h["chrom"], h["start"]))
    cur = None
    for h in hits:
        if cur and h["chrom"] == cur["contig"] and h["start"] < cur["end"]:
            cur["end"] = max(cur["end"], h["end"]); cur["hits"].append(h)
        else:
            cur = dict(assembly=asm, contig=h["chrom"], contig_length=h["chrom_len"], start=h["start"], end=h["end"], hits=[h]); loci.append(cur)
def entropy(p):
    c = collections.Counter(p); n = len(p)
    return round(-sum(v / n * math.log2(v / n) for v in c.values()), 2) if n else 0.0
with open(f"{wd}/loci.faa", "w") as fo:
    for i, l in enumerate(loci):
        b = max(l["hits"], key=lambda h: h["score"])
        pep = re.sub(r"[^A-Za-z]", "", b["trans"] or "").upper()
        l.update(id=f"L{i}", query=b["query"], query_description=qdesc.get(b["query"], ""), identity=b["ident"], query_cover=b["cov"], aligned_aa=b["qend"] - b["qstart"],
                 frameshifts=b["fs"], stops=b["stop"], n_queries=len({h["query"] for h in l["hits"]}), strand=b["strand"], exons=len(b["cds"]), pep_len=len(pep),
                 peptide_entropy=entropy(pep))
        if len(pep) >= 30: fo.write(f">L{i}\n{pep}\n")
        del l["hits"]
out = f"{wd}/loci_vs_db.tsv"
fresh = os.path.exists(out) and all(os.path.getmtime(out) > os.path.getmtime(f"{wd}/{a}.gff") for a in asms)
if not fresh: subprocess.run([D, "blastp", "-q", f"{wd}/loci.faa", "-d", db, "-o", out, "--sensitive", "-e", "1e-3", "-k", "50", "--threads", "2", "--quiet",
                "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore"], check=True)
# the proteome in the database that was annotated on the assembly being screened, if any
OWN_PROTEOME = {"elysia-chlorotica": "ELYCH", "elysia-marginata": "ELYMA"}
best = collections.defaultdict(dict); sac = collections.defaultdict(dict)
for line in open(out):
    q, s, pid, ln, ev, bs = line.rstrip().split("\t"); g = s.split("|")[0]; bs = float(bs)
    if g == "SAC":
        tag = s.split("|")[1]
        if bs > sac[q].get(tag, 0): sac[q][tag] = bs
        continue
    if g not in best[q] or bs > best[q][g][0]: best[q][g] = (bs, s, float(pid))
for l in loci:
    d = best.get(l["id"], {}); own_tag = OWN_PROTEOME.get(l["assembly"])
    met = d.get("MET", (0,))[0]; alga = max(d.get("ALG", (0,))[0], d.get("CYA", (0,))[0]); micro = max(d.get("BAC", (0,))[0], d.get("FUN", (0,))[0])
    slug = max((b for t, b in sac.get(l["id"], {}).items() if t != own_tag), default=0)      # other sacoglossan species
    own = sac.get(l["id"], {}).get(own_tag, 0)
    other = max(alga, micro)
    is_animal = (met > 0 or slug > 0) and (met >= other or (slug > other and met >= 0.9 * other))
    animal = max(met, slug) if is_animal else met
    l.update(animal_bits=max(met, slug), met_bits=met, other_sacoglossan_bits=slug, own_annotation_bits=own, alga_bits=alga, microbe_bits=micro,
             best_alga_hit=(d.get("ALG") or d.get("CYA") or (0, ""))[1])
    d = d or ({"SAC": 1} if slug else {})       # matched by nothing but its own annotation: no database match
    l["on_plastid_fragment"] = any(c == l["contig"] and l["start"] < b and l["end"] > a for c, a, b in frags.get(l["assembly"], []))
    w = whole.get((l["assembly"], l["contig"]))
    if l["on_plastid_fragment"]: l["class"] = "plastid fragment (assembly artefact)"
    elif w and animal < max(alga, micro, 1): l["class"] = "plastid-only contig" if w == "cp" else "bacterial contig"
    elif not d: l["class"] = "no database match"
    elif is_animal: l["class"] = "animal gene"
    elif alga >= micro: l["class"] = "alga-like"
    else: l["class"] = "microbe-like"
summary = {}
for asm in asms:
    c = collections.Counter(l["class"] for l in loci if l["assembly"] == asm)
    summary[asm] = dict(loci=sum(c.values()), **c)
n_q = len(qdesc)
json.dump(dict(n_queries=n_q, by_assembly=summary, loci=loci), open(f"{wd}/protein_screen_raw.json", "w"), indent=1)
print("queries", n_q)
for asm in asms: print(asm, summary[asm])
# ---- merge follow-up checks and write the summary used by the site
foreign = []
for asm in asms:
    p = f"{wd}/{asm}_foreign_contigs.json"
    if os.path.exists(p): foreign += json.load(open(p))
fmap = {(f["assembly"], f["contig"]): f for f in foreign}
contigs = collections.defaultdict(list)
for l in loci:
    if l["class"] in ("alga-like", "microbe-like"): contigs[(l["assembly"], l["contig"])].append(l)
rows = []
for (asm, c), ls in sorted(contigs.items()):
    f = fmap.get((asm, c), {})
    animal_here = sum(1 for l in loci if l["assembly"] == asm and l["contig"] == c and l["class"] == "animal gene")
    covered = sum(l["end"] - l["start"] for l in ls)
    row = dict(assembly=asm, contig=c, length=ls[0]["contig_length"], loci=len(ls), classes=sorted({l["class"] for l in ls}),
               genes=[l["query"].split("|", 1)[1] if l["query"].startswith("UDOFL|") else (l["query_description"].split()[0] if l["query_description"] else l["query"]) for l in ls],
               identity_min=min(l["identity"] for l in ls), identity_max=max(l["identity"] for l in ls), animal_loci_on_contig=animal_here,
               gc=f.get("gc"), annotated_genes=f.get("annotated_genes"), annotated_genes_with_protein_match=f.get("annotated_genes_with_protein_match"),
               read_alignments=f.get("read_alignments"), read_cover_span=f.get("read_cover_span"),
               aligned_to_other_assemblies={k[11:]: sum(a["end"] - a["start"] for a in v) for k, v in f.items() if k.startswith("aligned_to_")})
    in_others = sum(row["aligned_to_other_assemblies"].values()) if row["aligned_to_other_assemblies"] else 0
    if f and in_others > 0.2 * row["length"]: row["reading"] = "plastid-like sequence on a contig that is otherwise slug DNA: join tested with reads (script 29_join_test)"
    elif row["length"] <= 10000 and not animal_here: row["reading"] = "small contig with plastid genes only" if "alga-like" in row["classes"] else "small contig with bacterial genes only"
    elif f and not f.get("annotated_genes_with_protein_match") and (f.get("bp_covered_by_reads", 0) < 0.2 * row["length"]):
        row["reading"] = "DNA of another organism left in this assembly: absent from the reads of the xgElyCris1 animal"
    else: row["reading"] = "unexplained"
    rows.append(row)
nm = [l for l in loci if l["class"] == "no database match"]
tot = collections.Counter(l["class"] for l in loci)
json.dump(dict(n_queries=n_q, n_loci=len(loci), totals=dict(tot), by_assembly=summary, foreign_contigs=rows,
               no_match=dict(n=len(nm), low_complexity=sum(1 for l in nm if l["peptide_entropy"] < 3.3), short=sum(1 for l in nm if l["peptide_entropy"] >= 3.3 and l["aligned_aa"] < 120),
                             other=sum(1 for l in nm if l["peptide_entropy"] >= 3.3 and l["aligned_aa"] >= 120)),
               plastid_fragment_loci=[dict(assembly=l["assembly"], contig=l["contig"], start=l["start"], end=l["end"], query=l["query"], identity=l["identity"]) for l in loci if l["on_plastid_fragment"]],
               whole_contig_loci=[dict(assembly=l["assembly"], contig=l["contig"], length=l["contig_length"], query=l["query"], cls=l["class"]) for l in loci if l["class"] in ("plastid-only contig", "bacterial contig")],
               unexplained=sum(1 for r in rows if r["reading"] == "unexplained")), open(f"{wd}/protein_screen.json", "w"), indent=1)
print("totals", dict(tot)); print("no database match:", json.load(open(f"{wd}/protein_screen.json"))["no_match"])
for r in rows:
    print(f"{r['assembly']:7s} {r['contig']:20s} {r['length']:7d} bp loci={r['loci']} id {r['identity_min']:.2f}-{r['identity_max']:.2f} GC={r['gc']} genes={','.join(r['genes'])[:50]} | {r['reading']}")
