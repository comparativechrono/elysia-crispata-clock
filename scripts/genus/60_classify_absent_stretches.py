#!/usr/bin/env python3
"""What are the stretches of an assembly that are absent from egg DNA (59_germline_coverage.py)?

usage: 60_classify_absent_stretches.py CENSUS.tsv GENOME.fa.gz PANEL.fa alien_db.dmnd OUT_PREFIX
                                       [--own TAG] [--extra panel_extra.tsv] [--other NAME=GENOME.fa.gz ...]

Each stretch is compared with
  the reference panel     minimap2 (k 15, w 10): share of the stretch that matches plastid sequences,
                          and whether these are the chloroplast contigs of this assembly that are
                          part of the panel (--extra: config/panel_extra.tsv)
  a protein database      DIAMOND blastx as in 43_foreign_sequences.py; every coding stretch is
                          classed as animal, slug only, alga or microbe by its best matches. --own
                          names a proteome of the database that was predicted from this assembly
                          and therefore proves nothing (it is not counted)
  the other assemblies    of the species (minimap2 -x asm20): share of the stretch found in any of
                          them at 85% identity or more
and put into one class:
  chloroplast   a chloroplast contig of the panel; or plastid sequences cover at least 10% of it;
                or it carries genes of algal type and none of animal type
  plastid-like  a record of the panel that panel_extra.tsv describes as plastid-like: a sequence that
                resembles part of a plastid genome but carries no plastid gene
  bacterial     it carries genes of bacterial type and none of animal type
  slug          it carries genes of animal type, or at least half of it is found in another
                assembly of the species: DNA of the slug that the parents of the eggs lack
  unassigned    none of these
For every stretch the share of its 12-base words that are distinct is given as well. A stretch in
which fewer than half of the words are distinct consists largely of tandem repeats; short reads
are often missing from such sequences for technical reasons, so their absence says little.
Output: OUT_PREFIX.tsv (the census lines with class and evidence), OUT_PREFIX.json (totals by class),
        OUT_PREFIX.hidden.tsv (sequence, start, end of the chloroplast and bacterial stretches, to
        be left out of the reference for the reads of the adults; a stretch that is not a contig
        of the panel and covers at least 80% of its sequence is extended to the whole sequence,
        because the census has a resolution of one window)
"""
import sys, os, json, subprocess, collections, argparse, csv
ap = argparse.ArgumentParser()
ap.add_argument("census"); ap.add_argument("genome"); ap.add_argument("panel"); ap.add_argument("db"); ap.add_argument("out")
ap.add_argument("--own", default=""); ap.add_argument("--extra"); ap.add_argument("--other", nargs="*", default=[])
A = ap.parse_args()
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/"
rows = list(csv.DictReader(open(A.census), delimiter="\t"))
for r in rows: r["id"] = f'{r["sequence"]}:{int(r["start"]) + 1}-{r["end"]}'
wd = A.out + ".tmp"; os.makedirs(wd, exist_ok=True)
fa = f"{wd}/stretches.fa"
with open(f"{wd}/regions.txt", "w") as fo: fo.write("".join(r["id"] + "\n" for r in rows))
with open(fa, "w") as fo: subprocess.run([BIN + "samtools", "faidx", A.genome, "-r", f"{wd}/regions.txt"], stdout=fo, check=True)
def union(iv):
    iv = sorted(iv); tot = 0; ce = -1
    for s, e in iv:
        if s > ce: tot += e - s; ce = e
        elif e > ce: tot += e - ce; ce = e
    return tot
# ---- the panel
extra = set(); plastid_like = set()
if A.extra:
    for l in open(A.extra):
        f = l.split("\t")
        if f[0] != "id":
            extra.add(f[0])
            if f[4].startswith("Plastid-like"): plastid_like.add(f[0])
subprocess.run(f"{BIN}minimap2 -x map-ont -k15 -w10 -c --secondary=no -t 2 {A.panel} {fa} 2>/dev/null | cut -f1-12 > {wd}/panel.paf", shell=True, check=True)
cp = collections.defaultdict(list); own_contig = collections.defaultdict(list); best = collections.defaultdict(collections.Counter)
for l in open(f"{wd}/panel.paf"):
    c = l.split("\t")
    if int(c[10]) < 200 or not c[5].startswith(("cp|", "asmcp|")): continue
    cp[c[0]].append((int(c[2]), int(c[3]))); best[c[0]][c[5]] += int(c[10])
    if c[5] in extra: own_contig[c[0]].append((int(c[2]), int(c[3])))
# ---- proteins
tsv = f"{wd}/blastx.tsv"
if not os.path.exists(tsv):
    subprocess.run(["flock", "/tmp/bigmem.lock", BIN + "diamond", "blastx", "-q", fa, "-d", A.db, "-o", tsv + ".tmp", "--sensitive", "-e", "1e-5", "-F", "15", "--range-culling", "-k", "25",
                    "--threads", "2", "--quiet", "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qstart", "qend", "evalue", "bitscore"], check=True)
    os.replace(tsv + ".tmp", tsv)
hits = collections.defaultdict(list)
for line in open(tsv):
    q, s, pid, ln, qs, qe, ev, bs = line.rstrip().split("\t"); a, b = sorted((int(qs), int(qe)))
    hits[q].append((float(bs), a, b, s, float(pid)))
genes = {}
for q, hs in hits.items():
    st_ = []
    for bs, a, b, s, pid in sorted(hs, reverse=True):
        g, tag = s.split("|")[0], s.split("|")[1]
        if g == "SAC" and tag == A.own: continue
        key = {"MET": "met", "SAC": "slug", "ALG": "alga", "CYA": "alga", "BAC": "microbe", "FUN": "microbe"}[g]
        home = None
        for st in st_:
            ov = min(b, st["end"]) - max(a, st["start"])
            if ov >= 0.5 * (b - a) or ov >= 0.5 * (st["end"] - st["start"]): home = st; break
        if home is None:
            if any(min(b, st["end"]) - max(a, st["start"]) > 0.2 * (b - a) for st in st_): continue
            home = dict(start=a, end=b, best={}); st_.append(home)
        if bs > home["best"].get(key, 0): home["best"][key] = bs
    c = collections.Counter()
    for st in st_:
        d = st["best"]; met = d.get("met", 0); slug = d.get("slug", 0); other = max(d.get("alga", 0), d.get("microbe", 0))
        if (met or slug) and (met >= other or (slug > other and met >= 0.9 * other)): c["animal" if met else "slug_only"] += 1
        elif d.get("alga", 0) >= d.get("microbe", 0): c["alga"] += 1
        else: c["microbe"] += 1
    genes[q] = c
# ---- the other assemblies of the species
found = collections.defaultdict(list)
for o in A.other:
    name, path = o.split("=", 1); paf = f"{wd}/other.{name}.paf"
    if not os.path.exists(paf):
        subprocess.run(f"flock /tmp/bigmem.lock {BIN}minimap2 -x asm20 -c --secondary=no -t 2 -I 400M {path} {fa} 2>/dev/null | cut -f1-12 > {paf}.tmp && mv {paf}.tmp {paf}", shell=True, check=True)
    for l in open(paf):
        c = l.split("\t")
        if int(c[10]) >= 300 and int(c[9]) / int(c[10]) >= 0.85: found[c[0]].append((int(c[2]), int(c[3])))
# ---- word diversity: tandem repeats
seqs = {}; nm = None
for l in open(fa):
    if l[0] == ">": nm = l[1:].strip(); seqs[nm] = []
    else: seqs[nm].append(l.strip().upper())
def distinct_words(q, k=12):
    s = "".join(seqs.get(q, [])); n = len(s) - k + 1
    return round(len({s[i:i + k] for i in range(n)}) / n, 3) if n > 0 else 1.0
# ---- classes
out_rows = []
for r in rows:
    q = r["id"]; L = int(r["end"]) - int(r["start"]); g = genes.get(q, collections.Counter())
    pc = union(cp[q]) / L; oc = union(own_contig[q]) / L; fo_ = union(found[q]) / L
    animal = g["animal"] + g["slug_only"]
    ref = best[q].most_common(1)[0][0] if best[q] else ""
    if oc >= 0.5 and ref in plastid_like: cls = "plastid-like"
    elif oc >= 0.5 or pc >= 0.10 or (g["alga"] and not animal): cls = "chloroplast"
    elif g["microbe"] and not animal: cls = "bacterial"
    elif animal or fo_ >= 0.5: cls = "slug"
    else: cls = "unassigned"
    ref = best[q].most_common(1)[0][0] if best[q] else ""
    out_rows.append(dict(r, cls=cls, plastid_share=round(pc, 3), panel_contig_share=round(oc, 3), best_plastid_match=ref, genes_animal=g["animal"], genes_slug_only=g["slug_only"], genes_alga=g["alga"], genes_microbe=g["microbe"], other_assemblies_share=round(fo_, 3), distinct_words=distinct_words(q)))
cols = ["sequence", "sequence_length", "start", "end", "bases", "gc", "reads_found", "reads_expected", "whole_sequence", "cls", "plastid_share", "panel_contig_share", "best_plastid_match", "genes_alga", "genes_microbe", "genes_animal", "genes_slug_only", "other_assemblies_share", "distinct_words"]
with open(A.out + ".tsv", "w") as fo:
    fo.write("\t".join(cols).replace("cls", "class") + "\n")
    for r in out_rows: fo.write("\t".join(str(r[c]) for c in cols) + "\n")
with open(A.out + ".hidden.tsv", "w") as fo:
    for r in out_rows:
        if r["cls"] not in ("chloroplast", "bacterial", "plastid-like"): continue
        a, b, L = int(r["start"]), int(r["end"]), int(r["sequence_length"])
        if r["panel_contig_share"] < 0.5 and b - a >= 0.8 * L: a, b = 0, L
        fo.write(f'{r["sequence"]}\t{a}\t{b}\t{r["cls"]}\n')
summ = {}
for k in ("chloroplast", "plastid-like", "bacterial", "slug", "unassigned"):
    xs = [r for r in out_rows if r["cls"] == k]
    summ[k] = dict(stretches=len(xs), bases=sum(int(r["bases"]) for r in xs), whole_sequences=sum(1 for r in xs if r["whole_sequence"] == "yes"),
                   in_panel=sum(1 for r in xs if r["panel_contig_share"] >= 0.5), bases_in_panel=sum(int(r["bases"]) for r in xs if r["panel_contig_share"] >= 0.5),
                   tandem_repeats=sum(1 for r in xs if r["distinct_words"] < 0.5), bases_tandem_repeats=sum(int(r["bases"]) for r in xs if r["distinct_words"] < 0.5),
                   with_animal_genes=sum(1 for r in xs if r["genes_animal"] + r["genes_slug_only"] > 0), in_other_assemblies=sum(1 for r in xs if r["other_assemblies_share"] >= 0.5))
json.dump(dict(classes=summ, stretches=len(out_rows), bases=sum(int(r["bases"]) for r in out_rows)), open(A.out + ".json", "w"), indent=1)
print(json.dumps(summ))
