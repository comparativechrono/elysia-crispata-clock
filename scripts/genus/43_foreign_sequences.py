#!/usr/bin/env python3
"""What else is on a sequence that carries an alga-like or bacterial-type gene?

usage: 43_foreign_sequences.py PROTSCREEN_DIR alien_db.dmnd OUT.json slug=genome.fa.gz ...

The protein screen (scripts/28_protein_screen.py) searches with 471 proteins only, so it sees a
few genes of each sequence. Here every sequence on which it found a place whose best match is not
an animal protein is translated as a whole, in all six frames, and compared with the same protein
database (DIAMOND blastx with frameshifts and range culling, so that every gene along the
sequence is reported). Matches are merged into coding stretches, and each stretch is classified
by the rule of the protein screen:
  animal    the best match among animals other than sea slugs is at least as good as the best
            algal or microbial match (or within 10% of it while another sacoglossan beats both)
  alga      otherwise, best non-animal match algal, plant or cyanobacterial
  microbe   otherwise, best non-animal match bacterial, archaeal or fungal
  slug_only matched by sea-slug proteins only (the species' own annotation is not counted)
Sequences longer than 2 Mb are examined 100 kb either side of each place only.
Per sequence the output gives length, GC content, the number of coding stretches of each kind,
and the genera that the bacterial-type stretches match best (Swiss-Prot organism codes, turned
into names with UniProt's speclist.txt when SPECLIST points to a copy of
https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/speclist.txt).
"""
import sys, os, json, subprocess, collections
pdir, db, out = sys.argv[1:4]
genomes = dict(a.split("=", 1) for a in sys.argv[4:])
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin")
OWN_PROTEOME = {"elysia-chlorotica": "ELYCH", "elysia-marginata": "ELYMA"}
MAXLEN, FLANK = 2_000_000, 100_000
genus = {}
sl = os.environ.get("SPECLIST", "/home/claude/data/uniprot/speclist.txt")
if os.path.exists(sl):
    import re
    for line in open(sl, encoding="latin-1"):
        m = re.match(r"^(\w+)\s+[ABEVO]\s+\d+: N=(\S+)", line)
        if m: genus[m.group(1)] = m.group(2)
raw = json.load(open(os.path.join(pdir, "protein_screen_raw.json")))
flag = collections.defaultdict(lambda: collections.defaultdict(list))
for l in raw["loci"]:
    if l["class"] in ("alga-like", "microbe-like") and l["assembly"] in genomes: flag[l["assembly"]][l["contig"]].append(l)
wd = os.path.join(pdir, "foreign"); os.makedirs(wd, exist_ok=True)
fa = os.path.join(wd, "sequences.fa"); meta = {}
with open(fa, "w") as fo:
    for slug, cs in flag.items():
        for c, ls in cs.items():
            L = ls[0]["contig_length"]
            if L <= MAXLEN: regions = [(0, L)]
            else:
                regions = []
                for a, b in sorted((max(0, l["start"] - FLANK), min(L, l["end"] + FLANK)) for l in ls):
                    if regions and a <= regions[-1][1]: regions[-1] = (regions[-1][0], max(regions[-1][1], b))
                    else: regions.append((a, b))
            for a, b in regions:
                r = subprocess.run([BIN + "/samtools", "faidx", genomes[slug], f"{c}:{a + 1}-{b}"], capture_output=True, text=True, check=True).stdout
                seq = "".join(r.split("\n")[1:]).upper(); name = f"{slug}|{c}|{a}"
                fo.write(f">{name}\n{seq}\n")
                acgt = sum(seq.count(x) for x in "ACGT")
                meta[name] = dict(slug=slug, contig=c, offset=a, length=len(seq), contig_length=L, n=seq.count("N"), gc=round((seq.count("G") + seq.count("C")) / max(1, acgt), 4))
tsv = os.path.join(wd, "blastx.tsv")
if not (os.path.exists(tsv) and os.path.getmtime(tsv) > os.path.getmtime(os.path.join(pdir, "protein_screen_raw.json"))):
    subprocess.run([BIN + "/diamond", "blastx", "-q", fa, "-d", db, "-o", tsv, "--sensitive", "-e", "1e-5", "-F", "15", "--range-culling", "-k", "25",
                    "--threads", "2", "--quiet", "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qstart", "qend", "evalue", "bitscore"], check=True)
hits = collections.defaultdict(list)
for line in open(tsv):
    q, s, pid, ln, qs, qe, ev, bs = line.rstrip().split("\t"); a, b = sorted((int(qs), int(qe)))
    hits[q].append((float(bs), a, b, s, float(pid)))
res = collections.defaultdict(list)
for name, m in meta.items():
    own_tag = OWN_PROTEOME.get(m["slug"]); stretches = []
    for bs, a, b, s, pid in sorted(hits.get(name, []), reverse=True):
        g = s.split("|")[0]; tag = s.split("|")[1]
        key = "own" if (g == "SAC" and tag == own_tag) else {"MET": "met", "SAC": "slug", "ALG": "alga", "CYA": "alga", "BAC": "microbe", "FUN": "microbe"}[g]
        home = None
        for st in stretches:
            ov = min(b, st["end"]) - max(a, st["start"])
            if ov >= 0.5 * (b - a) or ov >= 0.5 * (st["end"] - st["start"]): home = st; break
        if home is None:
            if any(min(b, st["end"]) - max(a, st["start"]) > 0.2 * (b - a) for st in stretches): continue
            home = dict(start=a, end=b, best={}, top=s, top_identity=pid); stretches.append(home)
        if bs > home["best"].get(key, 0): home["best"][key] = bs
    rows = []
    for st in sorted(stretches, key=lambda x: x["start"]):
        d = st["best"]; met = d.get("met", 0); slug = d.get("slug", 0); other = max(d.get("alga", 0), d.get("microbe", 0))
        if not (met or slug or other): cls = "own_only"
        elif (met or slug) and (met >= other or (slug > other and met >= 0.9 * other)): cls = "animal" if met else "slug_only"
        elif d.get("alga", 0) >= d.get("microbe", 0): cls = "alga"
        else: cls = "microbe"
        rows.append(dict(start=st["start"] + m["offset"], end=st["end"] + m["offset"], cls=cls, best_hit=st["top"], best_identity=st["top_identity"],
                         bits={k: round(v) for k, v in d.items()}))
    c = collections.Counter(r["cls"] for r in rows)
    bac = [r for r in rows if r["cls"] == "microbe" and r["best_hit"].startswith("BAC|")]
    gen = collections.Counter(genus.get(r["best_hit"].split("|")[2].split("_")[-1], r["best_hit"].split("|")[2].split("_")[-1]) for r in bac)
    ids = sorted(r["best_identity"] for r in bac)
    res[m["slug"]].append(dict(contig=m["contig"], contig_length=m["contig_length"], examined_start=m["offset"], examined_length=m["length"], n_bases=m["n"], gc=m["gc"],
                               stretches=dict(c), coding_bp=sum(r["end"] - r["start"] for r in rows),
                               nearest_bacteria=gen.most_common(5), bacterial_identity_median=(ids[len(ids) // 2] if ids else None), list=rows))
# genome-wide GC for comparison (from a sample of each assembly: the first 20 Mb)
gc = {}
for slug, g in genomes.items():
    if slug not in res: continue
    p = subprocess.Popen(f"{BIN}/bgzip -dc {g} 2>/dev/null | head -c 20000000 | grep -v '>' | tr -d '\\n' | tr 'acgt' 'ACGT'", shell=True, stdout=subprocess.PIPE, text=True)
    seq = p.stdout.read(); acgt = sum(seq.count(x) for x in "ACGT"); gc[slug] = round((seq.count("G") + seq.count("C")) / max(1, acgt), 4)
json.dump(dict(genome_gc=gc, sequences=res), open(out, "w"), indent=1)
for slug, rows in res.items():
    print(slug, "genome GC", gc.get(slug))
    for r in sorted(rows, key=lambda r: -r["contig_length"]):
        print(f'  {r["contig"]:20s} {r["contig_length"]:8d} bp GC {r["gc"]:.3f} stretches {r["stretches"]} coding {r["coding_bp"]} bp | {", ".join(f"{g} {n}" for g, n in r["nearest_bacteria"][:4])} | median id {r["bacterial_identity_median"]}')
