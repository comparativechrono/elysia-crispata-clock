#!/usr/bin/env python3
"""Choose the reads of a pass-1 screen that go on to pass 2, and list hot panel windows.

usage: 31_select_hits.py RUN.p1.paf RUN.count OUT_PREFIX long|short
A panel window (200 bp) is "hot" when it collects short local matches (the aligned block is
under 600 bp of a read of 2 kb or more; for short reads, any match under 97% identity) from more
reads than a single-copy sequence could give at 20 times the sampled depth of an 800 Mb genome.
Such windows are slug repeats that resemble a stretch of a panel sequence. Reads whose only
hits are in hot windows are thinned to 60 per window. Reads that match nothing but ribosomal DNA
at slug-like identity (more than 5% from every algal record) are thinned to a sample of 300.
All other hit reads are kept.
"""
import sys, collections, json, random
paf, cnt, pre, kind = sys.argv[1:5]
n_reads, n_bases = (int(x) for x in open(cnt).read().split())
depth = n_bases / 8e8
hot_min = max(40, int(20 * depth))
random.seed(1)
hits = collections.defaultdict(list); shortc = collections.Counter(); rdna_only = {}
for line in open(paf):
    c = line.rstrip("\n").split("\t")
    qlen, qs, qe, t, ts, te, nm, al = int(c[1]), int(c[2]), int(c[3]), c[5], int(c[7]), int(c[8]), int(c[9]), int(c[10])
    local = (qe - qs < 600 and qlen >= 2000) if kind == "long" else (nm / max(1, al) < 0.97)
    b = (t, ((ts + te) // 2) // 200)
    hits[c[0]].append((b, local))
    if local: shortc[b] += 1
    # ribosomal DNA: the slug's own rDNA matches algal rDNA records at 80-90%. A read counts as
    # "own rDNA only" while every hit is to rDNA and none is within 5% of an algal record.
    cls = t.split("|")[0]
    dv = next((float(x[5:]) for x in c[12:] if x.startswith("dv:f:")), 1 - nm / max(1, al))
    own = cls == "slugrDNA" or (cls == "algrDNA" and dv > 0.05)
    rdna_only[c[0]] = rdna_only.get(c[0], True) and own
hot = {b for b, n in shortc.items() if n >= hot_min}
# neighbouring bins of a hot bin count as hot too
hot |= {(t, i + d) for (t, i) in list(hot) for d in (-1, 1)}
keep = []; per = collections.Counter(); only_hot = 0; n_rdna = 0; rdna_kept = 0
RDNA_SAMPLE = 300
for r, hh in hits.items():
    if rdna_only.get(r):
        n_rdna += 1
        if rdna_kept < RDNA_SAMPLE: rdna_kept += 1; keep.append(r)
    elif all(local and b in hot for b, local in hh):
        only_hot += 1
        b = hh[0][0]
        if per[b] < 60: per[b] += 1; keep.append(r)
    else:
        keep.append(r)
if kind == "short":      # both mates of every chosen pair
    spots = {r.rsplit(".", 1)[0] for r in keep}
    keep = [f"{s}.{i}" for s in sorted(spots) for i in (1, 2)]
open(pre + ".keep.names", "w").write("\n".join(keep) + "\n")
wins = collections.Counter()
for b in hot:
    if shortc[b]: wins[b] = shortc[b]
json.dump(dict(reads=n_reads, bases=n_bases, hit_reads=len(hits), reads_only_in_hot_windows=only_hot, own_rdna_reads=n_rdna, own_rdna_reads_kept=rdna_kept, kept=len(keep), hot_threshold=hot_min,
               hot_windows=[dict(target=t, start=i * 200, end=i * 200 + 200, reads=n) for (t, i), n in wins.most_common()]), open(pre + ".select.json", "w"), indent=1)
print(f"{len(hits)} reads with a hit; {n_rdna} match only ribosomal DNA at slug-like identity ({rdna_kept} kept as a sample); {only_hot} only in {len(wins)} hot windows (threshold {hot_min} reads); {len(keep)} kept")
