#!/usr/bin/env python3
"""Census of plastid barcode haplotypes (rbcL, tufA) carried by HiFi reads.

usage: 13_barcode_census.py reads.fa barcode_refs.fa barcode_db.fa outdir
barcode_refs.fa holds one reference per barcode (name = barcode id).  Reads are aligned
with minimap2 (asm20 preset, tolerant of ~20 % divergence); every read spanning a whole
barcode contributes its aligned sub-sequence.  Identical sub-sequences are counted, rare
variants within 3 edits of a more abundant haplotype are folded into it, and each
remaining haplotype is searched against the GenBank Bryopsidales barcode set.
"""
import sys, os, subprocess, collections, re
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/"
reads, refs, db, out = sys.argv[1:5]
os.makedirs(out, exist_ok=True)

def read_fa(p):
    d, n = {}, None
    for line in open(p):
        if line[0] == ">": n = line[1:].split()[0]; d[n] = []
        else: d[n].append(line.strip())
    return {k: "".join(v) for k, v in d.items()}
COMP = str.maketrans("ACGTacgt", "TGCAtgca")
ref = read_fa(refs)
titles = {}
for line in open(db):
    if line[0] == ">":
        f = line[1:].rstrip().split(" ", 1); titles[f[0]] = f[1] if len(f) > 1 else ""
paf = f"{out}/barcode.paf"
subprocess.run(f"{BIN}minimap2 -x asm20 -c --secondary=no -t 2 {refs} {reads} 2>/dev/null > {paf}", shell=True, check=True)
want = collections.defaultdict(dict)
SLACK = 30
for line in open(paf):
    c = line.split("\t")
    qlen, qs, qe = int(c[1]), int(c[2]), int(c[3])
    tlen, ts, te, nm = int(c[6]), int(c[7]), int(c[8]), int(c[9])
    if ts > SLACK or te < tlen - SLACK: continue        # must span (nearly) the whole barcode
    # extend the read interval to the reference ends that the aligner clipped
    left, right = ts, tlen - te
    if c[4] == "+": qs, qe = qs - left, qe + right
    else: qs, qe = qs - right, qe + left
    if qs < 0 or qe > qlen: continue
    bc = c[5].split("|")[0]
    if bc not in want[c[0]] or nm > want[c[0]][bc][3]:
        want[c[0]][bc] = (qs, qe, c[4], nm)
def fasta_iter(p):
    n, buf = None, []
    for line in open(p):
        if line[0] == ">":
            if n is not None: yield n, "".join(buf)
            n, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if n is not None: yield n, "".join(buf)

hap = collections.defaultdict(collections.Counter)
for n, s in fasta_iter(reads):
    if n not in want: continue
    for bc, (qs, qe, strand, _nm) in want[n].items():
        sub = s[qs:qe]
        if strand == "-": sub = sub.translate(COMP)[::-1]
        hap[bc][sub.upper()] += 1

def edits(a, b, maxd=3):
    """banded edit distance; returns maxd+1 when larger"""
    if abs(len(a) - len(b)) > maxd: return maxd + 1
    prev = list(range(len(b) + 1))
    for i in range(1, len(a) + 1):
        lo, hi = max(1, i - maxd), min(len(b), i + maxd)
        cur = [maxd + 1] * (len(b) + 1)
        if i <= maxd: cur[0] = i
        for j in range(lo, hi + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] != b[j - 1]))
        if min(cur[lo - 1:hi + 1]) > maxd: return maxd + 1
        prev = cur
    return prev[len(b)]

summary = open(f"{out}/barcode_haplotypes.tsv", "w")
summary.write("barcode\thaplotype\treads\tfraction\tlength\ttop_hit\ttop_hit_identity\ttop_hit_alnlen\tspecies_within_0.5pct_of_best\tsequence\n")
for bc in sorted(hap):
    cnt = hap[bc]; total = sum(cnt.values())
    order = [s for s, _ in cnt.most_common()]
    centres = []; merged = collections.Counter()
    for s in order:
        tgt = None
        if cnt[s] < 0.02 * total:        # only rare sequences can be folded into a neighbour
            for c in centres[:40]:
                if edits(s, c) <= 3: tgt = c; break
        if tgt is None: centres.append(s); tgt = s
        merged[tgt] += cnt[s]
    keep = [(s, k) for s, k in merged.most_common() if k >= max(5, 0.002 * total)]
    # second pass: every read goes to the nearest kept haplotype (<= 15 edits), so proportions use all reads
    cen = [s for s, _ in keep]
    assigned = collections.Counter(); unassigned = 0
    for s, k in cnt.items():
        best, bd = None, 16
        for c in cen:
            d = 0 if s == c else edits(s, c, 15)
            if d < bd: best, bd = c, d
            if d == 0: break
        if best is None: unassigned += k
        else: assigned[best] += k
    keep = sorted([(s, assigned[s]) for s in cen], key=lambda x: -x[1])
    total_assigned = sum(k for _, k in keep)
    fa = f"{out}/{bc}.haplotypes.fa"
    with open(fa, "w") as fo:
        for i, (s, k) in enumerate(keep, 1): fo.write(f">{bc}_h{i};reads={k}\n{s}\n")
    hits = collections.defaultdict(list)
    r = subprocess.run(f"{BIN}blastn -query {fa} -subject {db} -outfmt '6 qseqid sseqid pident length bitscore' -max_hsps 1 -evalue 1e-20 2>/dev/null",
                       shell=True, capture_output=True, text=True).stdout
    for line in r.splitlines():
        q, s, pid, ln, bs = line.split("\t")
        hits[q].append((float(bs), float(pid), int(ln), s, titles.get(s, "")))
    print(f"== {bc}: {total} spanning reads, {len(cnt)} distinct sequences, {len(keep)} haplotypes; {total_assigned} reads assigned, {unassigned} unassigned")
    for i, (s, k) in enumerate(keep, 1):
        hh = sorted(hits.get(f"{bc}_h{i};reads={k}", []), reverse=True)
        if hh:
            top = hh[0]
            def species(t):
                w = t.replace("UNVERIFIED: ", "").split()
                return " ".join(w[:2]) if len(w) >= 2 else t
            # species among hits within 0.5 % identity of the best and >= 90 % of its aligned length
            same = collections.Counter(species(h[4]) for h in hh if h[1] >= top[1] - 0.5 and h[2] >= 0.9 * top[2])
            same = [f"{k} ({v})" for k, v in same.most_common()]
            sp = species(top[4]) + " [" + top[3] + "]"
            summary.write(f"{bc}\th{i}\t{k}\t{k / max(1, total_assigned):.4f}\t{len(s)}\t{sp}\t{top[1]:.2f}\t{top[2]}\t{'; '.join(same[:6])}\t{s}\n")
            print(f"   h{i:<3d} {k:7d} reads ({100 * k / max(1, total_assigned):5.1f}%) len {len(s)}  best: {top[1]:.2f}% over {top[2]} bp  {sp}   | species within 0.5% of best: {'; '.join(same[:5])}")
        else:
            summary.write(f"{bc}\th{i}\t{k}\t{k / max(1, total_assigned):.4f}\t{len(s)}\t-\t\t\t\t{s}\n")
            print(f"   h{i:<3d} {k:7d} reads ({100 * k / max(1, total_assigned):5.1f}%) len {len(s)}  no hit in barcode set")
summary.close()
