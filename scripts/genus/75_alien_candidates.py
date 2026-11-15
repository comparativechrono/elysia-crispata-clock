#!/usr/bin/env python3
"""Where the genes of non-animal origin found in four sacoglossan proteomes lie, and what kind of DNA carries them.

usage: 75_alien_candidates.py ALIEN_PROTEOMES.json PROTEOMES.fa ALIEN_DB.dmnd REGSCREEN_DIR WORKDIR OUT.json

ALIEN_PROTEOMES.json and PROTEOMES.fa come from 72_alien_proteomes.py, REGSCREEN_DIR from
71_regulator_screen.sh (its embedded_plastid.json: kleptoplast DNA inside assembled sequences).

1. The candidates are compared with the database again without a cap on the number of matches
   (DIAMOND --more-sensitive, E <= 1e-5, 1000 targets): with 25 targets, large families of algal and
   sea-slug proteins can crowd out the animal homologue. A candidate stays when its best algal, plant,
   cyanobacterial, bacterial or fungal match is at least 80 bits and better than its best match among
   the 13 animals that are not sea slugs (matches to other sea slugs are reported, not counted); it is
   strong when no animal matches or the margin is at least 1.2-fold.
2. Each is placed on its assembly: E. crispata from the Ensembl annotation (chromosome names turned
   into the INSDC sequence names by length), E. chlorotica and E. marginata through the locus tag of
   the UniProt entry (GN=) and the gene in the GenBank annotation. Plakobranchus ocellatus is not
   placed (its assembly is not part of this study).
3. Every sequence that carries one is translated and compared with the database
   (43_foreign_sequences.py), and judged by foreign_dna.py: kleptoplast DNA (the candidate overlaps a
   plastid stretch found by the DNA screen), a foreign sequence (bacterium, phage, virus or other
   contaminant), a foreign stretch at the end of a slug sequence, or slug DNA. Only candidates on slug DNA could be genes that an ancestor
   acquired; they are listed with their best matches, whether the other sea slugs have them, and
   whether their function concerns the chloroplast, light or signalling (by the words of their
   UniProt description and of their best algal match), and sorted into transposable elements (by
   the same words), conserved families whose algal or microbial match is less than 1.2 times the
   animal one, and candidate transfers.
"""
import collections
import gzip
import json
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import foreign_dna as fd                                         # noqa: E402

cand_json, prot_fa, db, regdir, wd, out = sys.argv[1:7]
BIN = fd.BIN
DATA = fd.DATA
HERE = os.path.dirname(os.path.abspath(__file__))
ASM = {"ELYCR": "elysia-crispata", "ELYCH": "elysia-chlorotica", "ELYMA": "elysia-marginata"}
UNIPROT = {t: os.path.join(DATA, "refprot", f"{t}.fa.gz") for t in ("ELYCH", "ELYMA", "PLAOC")}
ENSEMBL_CHROM = os.path.join(DATA, "ensembl", "Elysia_crispata-GCA_963854125.1-chromosomes.tsv.gz")
PLASTID_WORDS = re.compile(r"chloroplast|plastid|thylakoid|photosystem|photosynth|chlorophyll|light-harvesting|ferredoxin|thioredoxin|"
                           r"sigma factor|sigma-70|RNA polymerase sigma|phytochrome|cryptochrome|phototropin|aureochrome|LOV|blue.light|"
                           r"tetrapyrrole|porphyrin|protochlorophyllide|magnesium.chelatase|heme oxygenase|carotenoid|xanthophyll|"
                           r"Tic|Toc|translocon|Clp|FtsZ|FtsH|Deg|PsbP|PsbO|Rubisco|ribulose|plastoquinone|plastocyanin|"
                           r"kinase|phosphatase|histidine kinase|response regulator|two-component|calcium|calmodulin|redox|signal", re.I)
TE_WORDS = re.compile(r"transposase|transposon|reverse transcriptase|RNA-dependent DNA polymerase|integrase|Gag-Pol|retrotranspos|"
                      r"DDE Tnp|ALP1|Helitron|PIF1|ATP.dependa?e?nt DNA helicase|DNA helicase|CCHC-type zinc finger|Pao|jockey|piggyBac|"
                      r"Harbinger|Mariner|Tc1", re.I)
os.makedirs(wd, exist_ok=True)
J = json.load(open(cand_json))
cands = J["candidates"]


def fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, seq = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].rstrip(), []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


# ---------------------------------------------------------------- 1. matches without a cap
ids = {c["query"] for c in cands}
qfa, tsv = os.path.join(wd, "cands.faa"), os.path.join(wd, "cands_vs_db.tsv")
if not os.path.exists(tsv):
    with open(qfa, "w") as fo:
        for h, s in fasta(prot_fa):
            if h.split()[0] in ids:
                fo.write(f">{h}\n{s}\n")
    subprocess.run([BIN + "/diamond", "blastp", "-q", qfa, "-d", db, "--more-sensitive", "-e", "1e-5", "-k", "1000", "-b", "0.5", "-c", "2",
                    "-p", "2", "--quiet", "-o", tsv + ".tmp", "-f", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore"],
                   check=True)
    os.rename(tsv + ".tmp", tsv)
best = collections.defaultdict(dict)
for line in open(tsv):
    q, s, pid, ln, ql, sl, ev, bs = line.rstrip("\n").split("\t")
    g = s.split("|")[0]
    if g == "SAC" and s.split("|")[1] == q.split("|")[0]:
        continue
    bs = float(bs)
    if g not in best[q] or bs > best[q][g][0]:
        best[q][g] = (bs, s, float(pid))
kept = []
for c in cands:
    d = best[c["query"]]
    met = d.get("MET", (0, "", 0))
    alien = max(((d[g][0], g) for g in ("ALG", "CYA", "BAC", "FUN") if g in d), default=(0, None))
    if alien[0] < 80 or alien[0] <= met[0]:
        continue
    g = alien[1]
    c.update(group=g, top_hit=d[g][1], top_bits=d[g][0], top_identity=d[g][2], best_animal_bits=met[0], best_animal_hit=met[1],
             best_alga_bits=max(d.get("ALG", (0,))[0], d.get("CYA", (0,))[0]), best_bacterium_bits=d.get("BAC", (0,))[0],
             best_fungus_bits=d.get("FUN", (0,))[0], other_sacoglossan_bits=d.get("SAC", (0,))[0], other_sacoglossan_hit=d.get("SAC", (0, ""))[1],
             strong=bool(met[0] == 0 or alien[0] >= 1.2 * met[0]))
    if g in ("BAC", "FUN"):
        c["lineage"] = {"BAC": "bacterium or archaeon", "FUN": "fungus"}[g]
    kept.append(c)

# descriptions of the best matches: algal proteomes for algae, plants and the cyanobacterium, Swiss-Prot
# for bacteria and fungi
want = collections.defaultdict(list)
for c in kept:
    want[c["top_hit"].split("|")[2]].append(c)
algdir = os.path.join(DATA, "algae")
for p in sorted(os.listdir(algdir)):
    if not p.endswith((".fa.gz", ".fasta.gz")):
        continue
    with gzip.open(os.path.join(algdir, p), "rt") as fh:
        for line in fh:
            if line[0] != ">":
                continue
            first = line[1:].split()[0]
            acc = first.split("|")[1] if first.count("|") >= 2 else first
            if acc in want:
                dsc = line[1:].rstrip().split(" ", 1)[1].split(" OS=")[0].split(" [")[0] if " " in line.strip() else ""
                for c in want[acc]:
                    c["top_hit_description"] = dsc

# ---------------------------------------------------------------- 2. places
gn = {}
for t, p in UNIPROT.items():
    for h, _ in fasta(p):
        m = re.search(r" GN=(\S+)", h)
        gn[f"{t}|{h.split('|')[1]}"] = m.group(1) if m else None
gene = {}
for t in ("ELYCH", "ELYMA"):
    with gzip.open(os.path.join(DATA, "asm", ASM[t], "genomic.gff.gz"), "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] == "gene":
                m = re.search(r"locus_tag=([^;\n]+)", f[8])
                if m:
                    gene[(t, m.group(1))] = (f[0], int(f[3]), int(f[4]))
by_len = collections.defaultdict(list)
for c, n in fd.lengths(ASM["ELYCR"]).items():
    by_len[n].append(c)
chrom = {}
with gzip.open(ENSEMBL_CHROM, "rt") as fh:
    for line in fh:
        name, n = line.split()[:2]
        if len(by_len[int(n)]) == 1:
            chrom[name] = by_len[int(n)][0]
for c in kept:
    t = c["species"]
    if t == "ELYCR" and c.get("location"):
        ctg, span = c["location"].rsplit(":", 1)
        a, b = map(int, span.split("-"))
        c["place"] = (chrom.get(ctg, ctg), a, b)
    elif t in ("ELYCH", "ELYMA") and (t, gn.get(c["query"])) in gene:
        c["place"] = gene[(t, gn[c["query"]])]
        c["locus_tag"] = gn[c["query"]]

# ---------------------------------------------------------------- 3. what carries them
pdir = os.path.join(wd, "context")
os.makedirs(pdir, exist_ok=True)
raw = os.path.join(pdir, "protein_screen_raw.json")
loci = []
for c in kept:
    if "place" in c:
        asm = ASM[c["species"]]
        ctg, a, b = c["place"]
        loci.append(dict(id=c["query"], assembly=asm, contig=ctg, contig_length=fd.lengths(asm)[ctg], start=a, end=b,
                         **{"class": "alga-like" if c["group"] in ("ALG", "CYA") else "microbe-like"}))
new = json.dumps(dict(loci=loci), indent=1)
if not os.path.exists(raw) or open(raw).read() != new:
    open(raw, "w").write(new)
ctx_json = os.path.join(pdir, "foreign_sequences.json")
if not os.path.exists(ctx_json) or os.path.getmtime(ctx_json) < os.path.getmtime(raw):
    subprocess.run([sys.executable, os.path.join(HERE, "43_foreign_sequences.py"), pdir, db, ctx_json] +
                   [f"{a}={fd.genome(a)}" for a in sorted({l["assembly"] for l in loci})], check=True)
F = json.load(open(ctx_json))
regions = collections.defaultdict(list)
for asm, L in F["sequences"].items():
    for s in L:
        regions[(asm, s["contig"])].append(s)
emb = json.load(open(os.path.join(regdir, "embedded_plastid.json")))
bg = {asm: fd.background(asm) for asm in sorted({l["assembly"] for l in loci})}
judged = {}
for c in kept:
    if "place" not in c:
        c["verdict"] = "not placed (no assembly)" if c["species"] == "PLAOC" else "not placed"
        continue
    asm = ASM[c["species"]]
    ctg, a, b = c["place"]
    c["assembly"], c["sequence"], c["start"], c["end"] = asm, ctg, a, b
    frag = next((r for r in emb.get(asm, []) if r["contig"] == ctg and a < r["end"] and b > r["start"]), None)
    if frag:
        c["verdict"] = "kleptoplast DNA"
        c["plastid_stretch"] = f"{ctg}:{frag['start']}-{frag['end']}"
        continue
    s = next(s for s in regions[(asm, ctg)] if s["examined_start"] <= a < s["examined_start"] + s["examined_length"])
    key = (asm, ctg, s["examined_start"])
    if key not in judged:
        lo, hi = s["examined_start"], s["examined_start"] + s["examined_length"]
        mids = [(l["start"] + l["end"]) // 2 for l in loci if l["assembly"] == asm and l["contig"] == ctg and lo <= (l["start"] + l["end"]) // 2 < hi]
        region = None if s["examined_length"] >= s["contig_length"] else (s["examined_start"], s["examined_start"] + s["examined_length"])
        judged[key] = dict(fd.judge(asm, ctg, s, mids, bg[asm], region, F["genome_gc"].get(asm)), gc=s["gc"], stretches=s["stretches"],
                           nearest_bacteria=s["nearest_bacteria"], bacterial_identity_median=s["bacterial_identity_median"])
    j = judged[key]
    seg = j["foreign_stretch"]
    v = j["verdict"]
    if v == "foreign stretch at the end of a slug sequence" and not (seg["start"] <= (a + b) // 2 < seg["end"]):
        v = "slug DNA"
    c["verdict"] = v
    c["sequence_context"] = {k: j[k] for k in ("orf_share", "genetic_code", "gc", "stretches", "stretch_kinds", "nearest_bacteria", "bacterial_identity_median",
                                               "foreign_stretch")}

# ---------------------------------------------------------------- 4. summary
summary = {}
for t in ("ELYCR", "ELYCH", "ELYMA", "PLAOC"):
    ks = [c for c in kept if c["species"] == t]
    summary[t] = dict(proteins=J["summary"][t]["proteins"], first_pass=J["summary"][t].get("listed", 0), candidates=len(ks),
                      alga_like=sum(c["group"] in ("ALG", "CYA") for c in ks), microbe_like=sum(c["group"] in ("BAC", "FUN") for c in ks),
                      verdicts=dict(collections.Counter(c["verdict"] for c in ks)),
                      alga_like_verdicts=dict(collections.Counter(c["verdict"] for c in ks if c["group"] in ("ALG", "CYA"))))
on_slug = [c for c in kept if c["verdict"] == "slug DNA"]
for c in on_slug:
    words = " ".join([c.get("description") or "", c.get("top_hit_description") or ""])
    c["plastid_light_or_signalling_words"] = sorted({m.group(0).lower() for m in PLASTID_WORDS.finditer(words)})
    if TE_WORDS.search(words):
        c["kind"] = "transposable element"
    elif not c["strong"]:
        c["kind"] = "conserved family (algal or microbial match less than 1.2 times the animal one)"
    else:
        c["kind"] = "candidate transfer"
for t in summary:
    summary[t]["on_slug_dna"] = dict(collections.Counter(c["kind"] for c in on_slug if c["species"] == t))
FIELDS = ("query", "species", "locus_tag", "length", "description", "group", "lineage", "top_hit", "top_hit_description", "top_bits",
          "top_identity", "best_animal_bits", "best_animal_hit", "best_alga_bits", "best_bacterium_bits", "other_sacoglossan_bits",
          "other_sacoglossan_hit", "strong", "assembly", "sequence", "start", "end", "verdict", "plastid_stretch", "sequence_context",
          "kind", "plastid_light_or_signalling_words")
res = dict(summary=summary, background_orf_share=bg,
           on_slug_dna=[{k: c.get(k) for k in FIELDS if k in c} for c in sorted(on_slug, key=lambda c: (c["species"], -c["top_bits"]))],
           candidates=[{k: c.get(k) for k in FIELDS if k in c} for c in kept])
json.dump(res, open(out, "w"), indent=1, default=str)
print(json.dumps(summary, indent=1))
print(len(on_slug), "candidates on slug DNA")
