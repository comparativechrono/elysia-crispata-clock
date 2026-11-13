#!/usr/bin/env python3
"""Summarise the assembly-level organelle screen.

usage: 23_assembly_screen_summary.py panel.tsv outdir name=fai ... -- name=paf ...
For each assembly: contigs that are plastid / bacterial / mitochondrial as a whole, and
plastid- or algal-mitochondrion-like segments embedded in otherwise slug contigs.
"""
import sys, os, json, collections, re
args = sys.argv[1:]
panel_tsv, out = args[0], args[1]
sep = args.index("--")
fais = dict(a.split("=", 1) for a in args[2:sep]); pafs = dict(a.split("=", 1) for a in args[sep + 1:])
desc = {}
for line in open(panel_tsv):
    f = line.rstrip("\n").split("\t"); desc[f[0]] = f[4] if len(f) > 4 else ""
def grp(t):
    c = t.split("|")[0]
    return "bac" if c.startswith("bac_") else c
def union(iv):
    iv = sorted(iv); tot = 0; cur = None
    for a, b in iv:
        if cur is None: cur = [a, b]
        elif a <= cur[1]: cur[1] = max(cur[1], b)
        else: tot += cur[1] - cur[0]; cur = [a, b]
    if cur: tot += cur[1] - cur[0]
    return tot
res = {}
for name, paf in pafs.items():
    clen = {}
    for line in open(fais[name]):
        f = line.split("\t"); clen[f[0]] = int(f[1])
    hits = collections.defaultdict(list)
    for line in open(paf):
        c = line.rstrip("\n").split("\t")
        m = re.match(r"(.+)_sliding:(\d+)-(\d+)$", c[0])
        contig, off = m.group(1), int(m.group(2)) - 1
        al = int(c[10])
        de = 0.0
        for t in c[12:]:
            if t.startswith("de:f:"): de = float(t[5:]); break
        if al < 100: continue
        hits[contig].append(dict(s=off + int(c[2]), e=off + int(c[3]), cls=grp(c[5]), t=c[5], ts=int(c[7]), te=int(c[8]), aln=al, ident=round(1 - de, 4)))
    whole = collections.Counter(); whole_bp = collections.Counter(); embedded = []; whole_list = []
    for contig, hh in hits.items():
        L = clen.get(contig, 0)
        cov = {}
        for g in set(h["cls"] for h in hh):
            cov[g] = union([(h["s"], h["e"]) for h in hh if h["cls"] == g and h["ident"] >= 0.8])
        top = max(cov, key=cov.get)
        if L and cov[top] >= 0.5 * L and top in ("cp", "cpx", "bac", "slugmt", "algmt"):
            whole[top] += 1; whole_bp[top] += L
            whole_list.append(dict(contig=contig, length=L, cls=top, covered=cov[top], best=max((h for h in hh if h["cls"] == top), key=lambda h: h["aln"])["t"]))
            continue
        # embedded algal-organelle-like segments (merge windows overlaps)
        for g in ("cp", "cpx", "algmt"):
            segs = sorted([h for h in hh if h["cls"] == g], key=lambda h: h["s"])
            merged = []
            for h in segs:
                if merged and h["s"] <= merged[-1]["e"] + 50:
                    if h["aln"] > merged[-1]["aln"]: merged[-1].update(t=h["t"], ts=h["ts"], te=h["te"], aln=h["aln"], ident=h["ident"])
                    merged[-1]["e"] = max(merged[-1]["e"], h["e"])
                else: merged.append(dict(h))
            for h in merged:
                embedded.append(dict(contig=contig, contig_length=L, start=h["s"], end=h["e"], length=h["e"] - h["s"], cls=g, target=h["t"], target_desc=desc.get(h["t"], ""), tstart=h["ts"], tend=h["te"], ident=h["ident"]))
    res[name] = dict(whole_contigs={k: dict(n=whole[k], bp=whole_bp[k]) for k in whole}, whole_list=sorted(whole_list, key=lambda x: -x["length"])[:200], embedded=embedded)
    print(f"== {name}: whole-contig classes {dict(whole)} ({ {k: v for k, v in whole_bp.items()} } bp); embedded algal-organelle-like segments: {len(embedded)}")
    c = collections.Counter()
    for x in embedded:
        b = "<150" if x["length"] < 150 else ("150-300" if x["length"] < 300 else ("300-1000" if x["length"] < 1000 else ">=1000"))
        i = ">=95%" if x["ident"] >= 0.95 else (">=85%" if x["ident"] >= 0.85 else "<85%")
        c[(b, i)] += 1
    print("   embedded by length x identity:", dict(sorted(c.items())))
    for x in sorted(embedded, key=lambda x: -x["length"] * x["ident"])[:12]:
        print(f"   {x['contig'][:22]:22s} L={x['contig_length']:>9d} {x['start']:>9d}-{x['end']:<9d} len={x['length']:5d} id={x['ident']:.3f} {x['target'][:22]:22s} {x['tstart']}-{x['tend']} {x['target_desc'][:40]}")
os.makedirs(out, exist_ok=True)
json.dump(res, open(os.path.join(out, "assembly_screen.json"), "w"), indent=1)
