#!/usr/bin/env python3
"""Count reads that support or contradict each assembled plastid/slug boundary.

usage: 27_boundary_evidence.py windows.tsv boundary.paf out.json
For each window (8 kb of slug flank + the plastid segment as assembled):
  through     reads aligned in one piece from >= 300 bp before to >= 300 bp after the boundary
              (they carry the junction exactly as assembled)
  diverge     reads aligned over >= 1 kb of the slug flank (>= 98 % identity) whose alignment
              stops within 50 bp of the boundary while >= 300 bp of the read continues beyond it:
              the animal's genome carries different sequence there, not plastid DNA
  flank_depth reads aligned over >= 1 kb of slug flank within 2 kb of the boundary
Independent molecules: reads that share an end are treated as PCR copies of one molecule, so a
read counts only if both of its ends differ by > 50 bp from every read already counted. The ends
are those of the whole read projected onto the window (alignment end plus any unaligned
overhang), so reads that stop aligning at the boundary or at the edge of the window are still
told apart by where the molecule itself begins and ends.
"""
import sys, json, collections
wtsv, paf, out = sys.argv[1:4]
W = {}
for line in list(open(wtsv))[1:]:
    f = line.rstrip("\n").split("\t")
    W[f[0]] = dict(assembly=f[1], contig=f[2], seg_start=int(f[3]), seg_end=int(f[4]), slug_side=f[5], win_start=int(f[6]), win_end=int(f[7]), bnd=int(f[8]))
A = collections.defaultdict(list)
for line in open(paf):
    c = line.split("\t")
    de = 0.0
    for t in c[12:]:
        if t.startswith("de:f:"): de = float(t[5:]); break
    a = dict(read=c[0], qlen=int(c[1]), qs=int(c[2]), qe=int(c[3]), strand=c[4], ts=int(c[7]), te=int(c[8]), ident=1 - de)
    lead, trail = (a["qs"], a["qlen"] - a["qe"]) if a["strand"] == "+" else (a["qlen"] - a["qe"], a["qs"])
    a["p0"] = a["ts"] - lead; a["p1"] = a["te"] + trail          # read ends projected onto the window
    A[c[5]].append(a)
def indep(al, ends=lambda a: (a["p0"], a["p1"])):
    kept = []
    for a in sorted(al, key=ends):
        x, y = ends(a)
        if all(abs(x - ends(k)[0]) > 50 and abs(y - ends(k)[1]) > 50 for k in kept): kept.append(a)
    return len(kept)
res = []
for w, info in W.items():
    b = info["bnd"]; left = info["slug_side"] == "left"
    through, diverge, depth = [], [], []
    for a in A.get(w, []):
        if a["ident"] < 0.98: continue
        slug_bp = (min(a["te"], b) - a["ts"]) if left else (a["te"] - max(a["ts"], b))
        if a["ts"] <= b - 300 and a["te"] >= b + 300:
            through.append(a); continue
        if slug_bp >= 1000:
            near_end = a["te"] if left else a["ts"]
            if abs(near_end - (b)) <= 2000: depth.append(a)
            # overhang of the read beyond the alignment end that faces the boundary
            if left: over = (a["qlen"] - a["qe"]) if a["strand"] == "+" else a["qs"]
            else: over = a["qs"] if a["strand"] == "+" else (a["qlen"] - a["qe"])
            if abs(near_end - b) <= 50 and over >= 300: diverge.append(a)
    dv = indep(diverge)
    res.append(dict(window=w, **info, through_reads=len(through), through_independent=indep(through), diverge_reads=len(diverge),
                    diverge_independent=dv, flank_reads=len(depth) + len(through)))
    print(f"{w[:44]:44s} slug on {info['slug_side']:5s} | through {len(through):3d} (indep {indep(through):2d}) | diverge at boundary {len(diverge):3d} (indep {dv:2d}) | slug-flank reads near boundary {len(depth) + len(through)}")
json.dump(res, open(out, "w"), indent=1)
t1 = sum(1 for r in res if r["through_independent"] <= 1); t2 = sum(1 for r in res if r["diverge_independent"] >= 2)
print(f"{len(res)} boundaries; {t1} carried by at most one independent molecule; {t2} contradicted by >= 2 independent reads")
print("totals: through molecules", sum(r["through_independent"] for r in res), "| diverging molecules", sum(r["diverge_independent"] for r in res), "| diverging reads", sum(r["diverge_reads"] for r in res))
