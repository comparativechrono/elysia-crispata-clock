#!/usr/bin/env python3
"""Summarise the pass-1 bait screen and choose reads for base-level re-alignment.

usage: 11_classify_bait_hits.py panel.tsv outdir pass1/*.bait.paf
Outputs (in outdir)
  read_classes.tsv     one line per read with a substantial hit: covered bp per class
  bait_summary.tsv     reads / bases per class and run
  cp_like.names        reads with any plastid / algal-mito / algal-rDNA hit  (-> barcode + chimera analysis)
  control.names        bacterial-MAG and slug-mito reads with >= 500 bp not explained by the panel
  clock.names          reads hitting clock-gene exon baits
"""
import sys, os, collections

panel_tsv, outdir = sys.argv[1:3]
pafs = sys.argv[3:]
os.makedirs(outdir, exist_ok=True)
MINLEN, MAXDV = 300, 0.25

def union(iv):
    iv = sorted(iv); tot = 0; cur = None
    for a, b in iv:
        if cur is None: cur = [a, b]
        elif a <= cur[1]: cur[1] = max(cur[1], b)
        else: tot += cur[1] - cur[0]; cur = [a, b]
    if cur: tot += cur[1] - cur[0]
    return tot

def grp(c):
    return "bac" if c.startswith("bac_") else c

summ = collections.defaultdict(lambda: collections.Counter())
fo = open(f"{outdir}/read_classes.tsv", "w")
fo.write("read\tlength\tmain_class\tcovered_bp_by_class\tunexplained_bp\tbest_target\tbest_dv\tbest_alnlen\n")
names = {k: open(f"{outdir}/{k}.names", "w") for k in ("cp_like", "control", "clock")}
for paf in pafs:
    run = os.path.basename(paf).split(".")[0]
    reads = collections.defaultdict(list)
    for line in open(paf):
        c = line.rstrip("\n").split("\t")
        al = int(c[10])
        if al < MINLEN: continue
        dv = 0.0
        for t in c[12:]:
            if t.startswith("dv:f:"): dv = float(t[5:]); break
        if dv > MAXDV: continue
        reads[c[0]].append((c[5].split("|")[0], c[5], int(c[1]), int(c[2]), int(c[3]), int(c[9]), al, dv))
    for r, hh in reads.items():
        L = hh[0][2]
        g = collections.defaultdict(list)
        for h in hh: g[grp(h[0])].append((h[3], h[4]))
        cov = {k: union(v) for k, v in g.items()}
        main = max(cov, key=cov.get)
        allcov = union([iv for v in g.values() for iv in v])
        best = max((h for h in hh if grp(h[0]) == main), key=lambda h: h[5])
        fo.write(f"{r}\t{L}\t{main}\t{';'.join(f'{k}:{v}' for k, v in sorted(cov.items()))}\t{L - allcov}\t{best[1]}\t{best[7]:.4f}\t{best[6]}\n")
        summ[run][main + "_reads"] += 1; summ[run][main + "_bp"] += cov[main]
        if any(k in g for k in ("cp", "cpx", "algmt", "algrDNA")) and not (set(g) <= {"algrDNA", "slugrDNA"} and main == "slugrDNA"):
            names["cp_like"].write(r + "\n")
        elif main in ("bac", "slugmt") and L - allcov >= 500:
            names["control"].write(r + "\n")
        if "clk" in g:
            names["clock"].write(r + "\n")
fo.close()
for f in names.values(): f.close()
classes = sorted({k.rsplit("_", 1)[0] for d in summ.values() for k in d})
with open(f"{outdir}/bait_summary.tsv", "w") as fo:
    fo.write("run\t" + "\t".join(f"{c}_reads\t{c}_bp" for c in classes) + "\n")
    tot = collections.Counter()
    for run, d in sorted(summ.items()):
        fo.write(run + "\t" + "\t".join(f"{d[c + '_reads']}\t{d[c + '_bp']}" for c in classes) + "\n")
        tot.update(d)
    fo.write("total\t" + "\t".join(f"{tot[c + '_reads']}\t{tot[c + '_bp']}" for c in classes) + "\n")
print(open(f"{outdir}/bait_summary.tsv").read())
