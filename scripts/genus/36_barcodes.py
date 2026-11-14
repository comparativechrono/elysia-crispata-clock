#!/usr/bin/env python3
"""Name the plastids in a read set from two barcode genes, rbcL and tufA.

usage: 36_barcodes.py reads|contigs  SEQS.fa  barcodes.fa  OUTDIR
barcodes.fa: GenBank rbcL and tufA records of Ulvophyceae and Xanthophyceae (titles kept).

reads    accurate long reads (HiFi). For each gene the GenBank record that most reads span from end
         to end is taken as the frame; every read spanning it contributes that stretch of its own
         sequence, identical stretches are counted (scripts/13_barcode_census.py) and each resulting
         haplotype is compared with all GenBank records.
contigs  assembled plastid contigs (short-read species). The stretch of a contig that matches the
         best GenBank record of each gene is cut out and compared with all records.
Output: OUTDIR/barcodes.json
"""
import sys, os, subprocess, collections, json, re
mode, seqs, db, out = sys.argv[1:5]
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/"
HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(out, exist_ok=True)
titles = {}; length = {}; n = None
for line in open(db):
    if line[0] == ">":
        f = line[1:].rstrip().split(" ", 1); n = f[0]; titles[n] = f[1] if len(f) > 1 else ""; length[n] = 0
    else: length[n] += len(line.strip())
def gene(t):
    t = t.lower()
    if "rbcl" in t or "ribulose" in t: return "rbcL"
    if "tufa" in t or "elongation factor tu" in t: return "tufA"
    return None
def species(t):
    w = re.sub(r"UNVERIFIED(_ORG|_ASMBLY)?: ", "", t).split()
    return " ".join(w[:2]) if len(w) >= 2 else t
def blast_top(fa):
    r = subprocess.run(f"{BIN}blastn -query {fa} -subject {db} -outfmt '6 qseqid sseqid pident length bitscore' -max_hsps 1 -evalue 1e-20 2>/dev/null", shell=True, capture_output=True, text=True).stdout
    hits = collections.defaultdict(list)
    for line in r.splitlines():
        q, s, pid, ln, bs = line.split("\t"); hits[q].append((float(bs), float(pid), int(ln), s, titles.get(s, "")))
    res = {}
    for q, hh in hits.items():
        hh.sort(reverse=True); top = hh[0]
        same = collections.Counter(species(h[4]) for h in hh if h[1] >= top[1] - 0.5 and h[2] >= 0.9 * top[2])
        res[q] = dict(top_hit=species(top[4]) + " [" + top[3] + "]", identity=top[1], alnlen=top[2], species_near_best="; ".join(f"{k} ({v})" for k, v in same.most_common(6)))
    return res
rows = []
if mode == "reads":
    paf = f"{out}/reads_vs_db.paf"
    subprocess.run(f"{BIN}minimap2 -x asm20 -c --secondary=no -t 2 {db} {seqs} 2>/dev/null | cut -f1-12 > {paf}", shell=True, check=True)
    span = collections.Counter(); idsum = collections.Counter()
    for line in open(paf):
        c = line.split("\t"); t = c[5]; tl, ts, te = int(c[6]), int(c[7]), int(c[8])
        if ts <= 30 and te >= tl - 30: span[t] += 1; idsum[t] += int(c[9]) / int(c[10])
    refs = {}
    for g, minlen in (("rbcL", 1000), ("tufA", 700)):
        cand = [(span[t], idsum[t] / span[t], length[t], t) for t in span if gene(titles[t]) == g and length[t] >= minlen and length[t] <= 1600]
        if not cand:
            cand = [(span[t], idsum[t] / span[t], length[t], t) for t in span if gene(titles[t]) == g and length[t] >= 400 and length[t] <= 1600]
        if cand:
            cand.sort(reverse=True); refs[g] = cand[0][3]
            print(f"{g}: frame {cand[0][3]} ({length[cand[0][3]]} bp, {titles[cand[0][3]][:60]}) spanned by {cand[0][0]} reads at mean identity {cand[0][1]:.3f}")
    if refs:
        ids = " ".join(f"'{t}'" for t in refs.values())
        fa = subprocess.run(f"{BIN}seqkit grep -w 0 -n -r " + " ".join(f"-p '^{re.escape(t)}( |$)'" for t in refs.values()) + f" {db}", shell=True, capture_output=True, text=True).stdout
        with open(f"{out}/frames.fa", "w") as fo:
            name = None
            for line in fa.splitlines():
                if line.startswith(">"):
                    acc = line[1:].split()[0]; g = next(k for k, v in refs.items() if v == acc); fo.write(f">{g}|{acc}\n")
                else: fo.write(line + "\n")
        subprocess.run(f"python3 {HERE}/../13_barcode_census.py {seqs} {out}/frames.fa {db} {out}/census > {out}/census.log 2>&1", shell=True, check=True)
        for i, line in enumerate(open(f"{out}/census/barcode_haplotypes.tsv")):
            if i == 0: continue
            f = line.rstrip("\n").split("\t")
            rows.append(dict(barcode=f[0], haplotype=f[1], reads=int(f[2]), fraction=float(f[3]), length=int(f[4]), top_hit=f[5], identity=float(f[6]) if f[6] else None,
                             alnlen=int(f[7]) if f[7] else None, species_near_best=f[8], sequence=f[9], frame=refs.get(f[0])))
else:
    r = subprocess.run(f"{BIN}blastn -query {seqs} -subject {db} -outfmt '6 qseqid sseqid pident length bitscore qstart qend sstart send' -max_hsps 1 -evalue 1e-30 2>/dev/null", shell=True, capture_output=True, text=True).stdout
    best = {}
    for line in r.splitlines():
        q, s, pid, ln, bs, qs, qe, ss, se = line.split("\t"); g = gene(titles.get(s, ""))
        if not g: continue
        k = (g, q); v = (int(ln), float(bs), float(pid), s, min(int(qs), int(qe)), max(int(qs), int(qe)))
        if k not in best or v > best[k]: best[k] = v
    cl = {}
    for line in subprocess.run(f"{BIN}seqkit fx2tab -n -l {seqs}", shell=True, capture_output=True, text=True).stdout.splitlines():
        f = line.split("\t"); cl[f[0].split()[0]] = int(f[-1])
    with open(f"{out}/contig_barcodes.fa", "w") as fo:
        for (g, q), (ln, bs, pid, s, a, b) in sorted(best.items()):
            if ln < 300: continue
            sub = subprocess.run(f"{BIN}samtools faidx {seqs} '{q}:{a}-{b}'", shell=True, capture_output=True, text=True).stdout
            sq = "".join(sub.splitlines()[1:])
            fo.write(f">{g}|{q}:{a}-{b}\n{sq}\n")
    top = blast_top(f"{out}/contig_barcodes.fa")
    name = None
    for line in open(f"{out}/contig_barcodes.fa"):
        if line[0] == ">": name = line[1:].strip(); continue
        g, loc = name.split("|", 1); t = top.get(name, {})
        m = re.search(r"cov_([\d.]+)|multi=([\d.]+)", loc)
        rows.append(dict(barcode=g, haplotype=loc, length=len(line.strip()), contig_length=cl.get(loc.rsplit(":", 1)[0]), top_hit=t.get("top_hit"), identity=t.get("identity"), alnlen=t.get("alnlen"),
                         species_near_best=t.get("species_near_best", ""), sequence=line.strip()))
# GenBank records are often shorter than the frame. The longest good match is not always the most
# similar record, so the most similar record with >= 400 aligned bases is reported as well.
if rows:
    with open(f"{out}/final_haps.fa", "w") as fo:
        for i, r in enumerate(rows): fo.write(f">q{i}\n{r['sequence']}\n")
    rr = subprocess.run(f"{BIN}blastn -query {out}/final_haps.fa -subject {db} -outfmt '6 qseqid sseqid pident length bitscore' -max_hsps 1 -max_target_seqs 50000 -evalue 1e-20 2>/dev/null", shell=True, capture_output=True, text=True).stdout
    byq = collections.defaultdict(list)
    for line in rr.splitlines():
        q, sid, pid, ln, bs = line.split("\t")
        if int(ln) >= 400 and gene(titles.get(sid, "")) == rows[int(q[1:])]["barcode"]: byq[int(q[1:])].append((float(pid), int(ln), sid))
    for i, r in enumerate(rows):
        hh = sorted(byq.get(i, []), reverse=True)
        if not hh: continue
        top = hh[0][0]
        same = [h for h in hh if h[0] >= top - 0.1]
        sp = collections.Counter(species(titles.get(h[2], "")) for h in same)
        r["most_similar"] = dict(identity=top, records=len(same), longest_alnlen=max(h[1] for h in same), species="; ".join(f"{k} ({v})" for k, v in sp.most_common(6)),
                                 accession=max(same, key=lambda h: h[1])[2])
json.dump(rows, open(f"{out}/barcodes.json", "w"), indent=1)
for r in rows:
    print(f"  {r['barcode']:5s} {str(r.get('haplotype'))[:34]:34s} reads {r.get('reads', '-')!s:>6s} frac {r.get('fraction', '-')!s:>6s} len {r['length']:5d}  best {r.get('identity')}% over {r.get('alnlen')} bp  {r.get('top_hit')}  | near best: {r.get('species_near_best', '')[:70]} | most similar >=400 bp: {r.get('most_similar', {}).get('identity')}% ({r.get('most_similar', {}).get('longest_alnlen')} bp) {r.get('most_similar', {}).get('species', '')[:80]}")
