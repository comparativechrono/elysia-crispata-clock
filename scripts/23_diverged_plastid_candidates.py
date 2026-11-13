#!/usr/bin/env python3
"""List the matches of the window screen (script 22) that are neither exact plastid copies nor
simple repeats, so that every remaining match can be looked at one by one.

usage: 23_diverged_plastid_candidates.py asm_screen_dir panel.fa out.tsv out.json
A match is kept when it is >= 200 bp at >= 80 % identity to an algal plastid, mitochondrial or
rDNA sequence of the panel and
  - the matched panel sequence is not low-complexity (entropy of its trinucleotides >= 4.5 of a
    possible 6 bits; microsatellites such as (CTT)n score about 3),
  - the panel window is not one hit by >= 15 different assembly windows (repeat-like),
  - it is not on one of the exact plastid fragments or whole organelle contigs already listed.
Overlapping matches are merged into loci; the best (longest) match describes each locus.
"""
import sys, os, re, json, math, collections, subprocess
d, panel, out_tsv, out_json = sys.argv[1:5]
ST = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/samtools"
emb = json.load(open(f"{d}/embedded_plastid.json")); scr = json.load(open(f"{d}/assembly_screen.json"))
whole = {a: {x["contig"] for x in v["whole_list"]} for a, v in scr.items()}
frag = {a: [(r["contig"], r["start"], r["end"]) for r in rows] for a, rows in emb.items()}
def entropy3(s):
    c = collections.Counter(s[i:i + 3] for i in range(len(s) - 2)); t = sum(c.values())
    return -sum(v / t * math.log2(v / t) for v in c.values()) if t else 0.0
res = {}; rows = []
for asm in ("sanger", "alt", "purdue", "tampa"):
    p = f"{d}/{asm}.paf"
    if not os.path.exists(p): continue
    hits = []; tb = collections.Counter()
    for l in open(p):
        f = l.rstrip("\n").split("\t"); t = f[5]; cls = t.split("|")[0]
        if not (cls.startswith("cp") or cls in ("algmt", "algrDNA")): continue
        m = re.match(r"(.+)_sliding:(\d+)-(\d+)$", f[0]); off = int(m.group(2)) - 1
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), 0.0)
        h = dict(contig=m.group(1), start=off + int(f[2]), end=off + int(f[3]), target=t, tstart=int(f[7]), tend=int(f[8]), identity=round(1 - de, 4), length=int(f[10]),
                 cls={"algmt": "algal mitochondrial", "algrDNA": "algal nuclear rDNA"}.get(cls, "plastid"))
        hits.append(h)
        for b in range(h["tstart"] // 300, h["tend"] // 300 + 1): tb[(t, b)] += 1
    hot = {k for k, v in tb.items() if v >= 15}
    keep = []
    for h in hits:
        bins = [(h["target"], b) for b in range(h["tstart"] // 300, h["tend"] // 300 + 1)]
        if h["length"] < 200 or h["identity"] < 0.80 or sum(b in hot for b in bins) > 0.5 * len(bins): continue
        if h["contig"] in whole.get(asm, set()): continue
        if any(c == h["contig"] and h["start"] < b + 50 and h["end"] > a - 50 for c, a, b in frag.get(asm, [])): continue
        keep.append(h)
    if keep:
        regs = [f'{h["target"]}:{h["tstart"] + 1}-{h["tend"]}' for h in keep]
        fa = subprocess.run([ST, "faidx", panel] + regs, capture_output=True, text=True, check=True).stdout
        seqs = ["".join(b.split("\n")[1:]).upper() for b in fa.split(">")[1:]]
        for h, s in zip(keep, seqs): h["target_entropy"] = round(entropy3(s), 2)
    n_low = sum(1 for h in keep if h["target_entropy"] < 4.5)
    keep = sorted((h for h in keep if h["target_entropy"] >= 4.5), key=lambda h: (h["contig"], h["start"]))
    loci = []
    for h in keep:
        if loci and loci[-1]["contig"] == h["contig"] and h["start"] < loci[-1]["end"] + 300:
            L = loci[-1]; L["end"] = max(L["end"], h["end"])
            if h["length"] > L["length"]: L.update({k: h[k] for k in ("target", "tstart", "tend", "identity", "length", "cls", "target_entropy")})
        else: loci.append(dict(assembly=asm, **h))
    res[asm] = dict(non_exact_matches=len(hits), low_complexity_removed=n_low, loci=len(loci), by_class=dict(collections.Counter(L["cls"] for L in loci)),
                    best_plastid_identity=max((L["identity"] for L in loci if L["cls"] == "plastid"), default=None))
    rows += loci
    print(asm, res[asm])
with open(out_tsv, "w") as fo:
    cols = ["assembly", "contig", "start", "end", "length", "identity", "cls", "target", "tstart", "tend", "target_entropy"]
    fo.write("\t".join(cols) + "\n")
    for L in rows: fo.write("\t".join(str(L[c]) for c in cols) + "\n")
json.dump(dict(summary=res, loci=rows), open(out_json, "w"), indent=1)
