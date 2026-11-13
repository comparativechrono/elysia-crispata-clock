#!/usr/bin/env python3
"""Evidence for or against each exact plastid copy embedded in an assembly.

usage: 24_embedded_plastid_evidence.py assembly_screen.json pass2.paf genes.bed out.json name=fasta ...
For every embedded segment (>= 97 % identical to a plastome, >= 200 bp) this reports
  * whether each end abuts an assembly gap (run of N) or a sequence end,
  * how many HiFi reads cross each plastid/slug boundary (>= 200 bp on both sides) and how many of
    them are independent molecules (both read ends differ by > 50 bp),
  * which plastid genes the segment carries.
A genuine insertion present in the animal's genome should be crossed by tens of independent
reads at each boundary; a boundary that exists in one amplified molecule is a library chimera
that the assembler built into a contig.
"""
import sys, json, collections, gzip
scr, paf, bed, out = sys.argv[1:5]
fas = dict(a.split("=", 1) for a in sys.argv[5:])
D = json.load(open(scr))
genes = collections.defaultdict(list)
for line in open(bed):
    f = line.rstrip("\n").split("\t"); genes[f[0]].append((int(f[1]), int(f[2]), f[3]))
aln = collections.defaultdict(list)
for line in open(paf):
    c = line.split("\t")
    if "|" in c[5]: continue
    aln[c[5]].append((int(c[7]), int(c[8]), c[0], int(c[2]), int(c[3]), int(c[1])))
def load_contigs(path, want):
    op = gzip.open if path.endswith(".gz") else open
    seqs = {}; name = None; buf = []
    with op(path, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                if name in want: seqs[name] = "".join(buf)
                name = line[1:].split()[0]; buf = []
            elif name in want: buf.append(line.strip())
    if name in want: seqs[name] = "".join(buf)
    return seqs
def independent(reads):
    kept = []
    for r in sorted(reads):
        if all(abs(r[0] - k[0]) > 50 and abs(r[1] - k[1]) > 50 for k in kept): kept.append(r)
    return len(kept)
res = {}
for name, fa in fas.items():
    segs = [x for x in D[name]["embedded"] if x["ident"] >= 0.97 and x["length"] >= 200]
    seqs = load_contigs(fa, {x["contig"] for x in segs})
    rows = []
    for x in sorted(segs, key=lambda x: -x["length"]):
        s = seqs[x["contig"]]; L = len(s)
        def edge(pos, direction):
            """distance from pos to the nearest N-run or sequence end, looking outward"""
            if direction < 0:
                w = s[max(0, pos - 20000):pos]; i = max(w.rfind("N"), w.rfind("n"))
                d_gap = (len(w) - 1 - i) if i >= 0 else None
                d_end = pos
            else:
                w = s[pos:pos + 20000]; i = min([j for j in (w.find("N"), w.find("n")) if j >= 0], default=-1)
                d_gap = i if i >= 0 else None
                d_end = L - pos
            if d_end <= 20000 and (d_gap is None or d_end <= d_gap): return "sequence end", d_end
            if d_gap is not None: return "gap", d_gap
            return "none within 20 kb", None
        sides = {}
        for label, pos, direction in (("left", x["start"], -1), ("right", x["end"], +1)):
            kind, dist = edge(pos, direction)
            span = [(a, b) for a, b, r, qs, qe, ql in aln.get(x["contig"], []) if a <= pos - 200 and b >= pos + 200]
            sides[label] = dict(boundary=pos, nearest_break=kind, distance=dist, reads_crossing=len(span), independent_molecules=independent(span))
        tacc = x["target"].split("|")[1]
        lo, hi = min(x["tstart"], x["tend"]), max(x["tstart"], x["tend"])
        g = sorted({n for a, b, n in genes.get(tacc, []) if b > lo and a < hi})
        rows.append(dict(contig=x["contig"], start=x["start"], end=x["end"], length=x["length"], identity=x["ident"], plastome=tacc,
                         plastome_desc=x["target_desc"], plastome_start=lo, plastome_end=hi, genes=g, left=sides["left"], right=sides["right"]))
    res[name] = rows
    print(f"== {name}: {len(rows)} exact plastid segments, {sum(r['length'] for r in rows)} bp")
    for r in rows:
        l, rt = r["left"], r["right"]
        print(f"  {r['contig'][:20]:20s} {r['start']:>9d}-{r['end']:<9d} {r['length']:5d} bp {r['identity']*100:5.1f}% {r['plastome']:11s} | left: {l['nearest_break']} {l['distance']} bp; reads crossing {l['reads_crossing']} (indep {l['independent_molecules']}) | right: {rt['nearest_break']} {rt['distance']} bp; reads crossing {rt['reads_crossing']} (indep {rt['independent_molecules']}) | {','.join(r['genes'][:5])}")
json.dump(res, open(out, "w"), indent=1)
