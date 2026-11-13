#!/usr/bin/env python3
"""Look inside the reads that matched a plastid genome only weakly in pass 1.

usage: 17_weak_plastid_reads.py read_classes.tsv pass2.paf panel.fa panel.tsv out.json
Reads whose best class is plastid but which are not confident plastid reads (<= 3 % divergence
over >= 1 kb) were set aside. Their base-level alignments (pass 2) are searched here for plastid
segments of >= 200 bp at >= 80 % identity whose matched panel sequence is not low-complexity
(trinucleotide entropy >= 4.5 bits). Such a segment is "exact" at >= 97 % identity (a short piece
of kleptoplast DNA) and "diverged" below that. For diverged segments the identity of the slug
part of the same read is reported: a read that is equally diverged everywhere is a low-accuracy
read, not diverged DNA.
"""
import sys, os, json, math, collections, subprocess
classes, paf, panel, ptsv, out = sys.argv[1:6]
ST = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/samtools"
desc = {}
for l in open(ptsv):
    f = l.rstrip("\n").split("\t"); desc[f[0]] = f[-1]
weak = set()
for l in list(open(classes))[1:]:
    f = l.rstrip("\n").split("\t")
    if f[2] in ("cp", "cpx") and not (float(f[6]) <= 0.03 and int(f[7]) >= 1000): weak.add(f[0])
A = collections.defaultdict(list)
for l in open(paf):
    f = l.rstrip("\n").split("\t")
    if f[0] not in weak: continue
    de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), 0.0); t = f[5]
    kind = "cp" if t.startswith(("cp|", "cpx|")) else ("slug" if "|" not in t else "other")
    A[f[0]].append(dict(qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), t=t, ts=int(f[7]), te=int(f[8]), ident=1 - de, ln=int(f[10]), kind=kind))
def entropy3(s):
    c = collections.Counter(s[i:i + 3] for i in range(len(s) - 2)); n = sum(c.values())
    return -sum(v / n * math.log2(v / n) for v in c.values()) if n else 0.0
cand = [(r, a) for r, al in A.items() for a in al if a["kind"] == "cp" and a["ln"] >= 200 and a["ident"] >= 0.80]
seqs = []
regs = [f'{a["t"]}:{a["ts"] + 1}-{a["te"]}' for r, a in cand]
for i in range(0, len(regs), 400):
    fa = subprocess.run([ST, "faidx", panel] + regs[i:i + 400], capture_output=True, text=True, check=True).stdout
    seqs += ["".join(b.split("\n")[1:]).upper() for b in fa.split(">")[1:]]
best = {}
for (r, a), s in zip(cand, seqs):
    if entropy3(s) < 4.5: continue
    if r not in best or a["ln"] > best[r]["ln"]: best[r] = a
exact = [r for r, a in best.items() if a["ident"] >= 0.97]
div = []
for r, a in best.items():
    if a["ident"] >= 0.97: continue
    slug = [x for x in A[r] if x["kind"] == "slug" and x["ln"] >= 300]
    sb = max(slug, key=lambda x: x["ln"]) if slug else None
    div.append(dict(read=r, read_length=a["qlen"], plastid_like_bp=a["qe"] - a["qs"], identity=round(a["ident"], 4), plastome=a["t"].split("|")[-1], plastome_species=" ".join(desc.get(a["t"], "").split()[:2]),
                    slug_bp=(sb["qe"] - sb["qs"]) if sb else 0, slug_identity=round(sb["ident"], 4) if sb else None))
div.sort(key=lambda d: -d["plastid_like_bp"])
res = dict(weak_reads=len(weak), realigned=len(A), with_exact_piece=len(exact),
           exact_with_slug=sum(1 for r in exact if any(x["kind"] == "slug" and x["ln"] >= 300 and x["ident"] >= 0.97 for x in A[r])),
           diverged=div, diverged_joined_to_slug=sum(1 for d in div if d["slug_bp"] and abs((1 - d["identity"]) - (1 - d["slug_identity"])) > 0.02))
json.dump(res, open(out, "w"), indent=1)
print({k: v for k, v in res.items() if k != "diverged"})
for d in div: print(" ", d)
