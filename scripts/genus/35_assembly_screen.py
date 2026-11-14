#!/usr/bin/env python3
"""Organelle-like DNA in an assembly, from the 10 kb window screen (scripts/22_assembly_organelle_screen.sh
run with the genus panel).

usage: 35_assembly_screen.py panel.tsv panel.fa genome.fa.gz.fai windows.paf out.json
Reports
  whole      sequences of the assembly that are organelle DNA as a whole (>= 50% covered by
             plastid, algal-mitochondrial or slug-mitochondrial matches)
  loci       every other place that matches an algal plastid, mitochondrial or nuclear rDNA
             sequence over >= 200 bp at >= 80% identity, after removing
               - matches to low-complexity panel sequence (trinucleotide entropy < 4.5 of 6 bits)
               - matches to panel windows hit from >= 15 places in the assembly (slug repeats that
                 resemble a stretch of a panel sequence)
             Overlapping matches are merged; the longest describes the locus.
  repeats    the panel windows removed as repeat-like, with the number of places that hit them
"""
import sys, os, re, json, math, collections, subprocess
panel_tsv, panel_fa, fai, paf, out = sys.argv[1:6]
ST = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/samtools"
desc = {}
for l in list(open(panel_tsv))[1:]:
    f = l.rstrip("\n").split("\t"); desc[f[0]] = f[4]
clen = {}
for l in open(fai):
    f = l.split("\t"); clen[f[0]] = int(f[1])
def union(iv):
    iv = sorted(iv); tot = 0; cur = None
    for a, b in iv:
        if cur is None: cur = [a, b]
        elif a <= cur[1]: cur[1] = max(cur[1], b)
        else: tot += cur[1] - cur[0]; cur = [a, b]
    if cur: tot += cur[1] - cur[0]
    return tot
def entropy3(s):
    c = collections.Counter(s[i:i + 3] for i in range(len(s) - 2)); t = sum(c.values())
    return -sum(v / t * math.log2(v / t) for v in c.values()) if t else 0.0
NAMES = {"cp": "plastid", "algmt": "algal mitochondrial", "algrDNA": "algal nuclear rDNA", "slugmt": "slug mitochondrial", "slugrDNA": "slug rDNA"}
hits = collections.defaultdict(list)
for line in open(paf):
    c = line.rstrip("\n").split("\t")
    m = re.match(r"(.+)_sliding:(\d+)-(\d+)$", c[0]); contig, off = m.group(1), int(m.group(2)) - 1
    al = int(c[10])
    if al < 100: continue
    de = next((float(t[5:]) for t in c[12:] if t.startswith("de:f:")), 0.0)
    hits[contig].append(dict(contig=contig, start=off + int(c[2]), end=off + int(c[3]), cls=c[5].split("|")[0], target=c[5], tstart=int(c[7]), tend=int(c[8]), length=al, identity=round(1 - de, 4)))
whole = []; rest = []
for contig, hh in hits.items():
    L = clen.get(contig, 0)
    cov = {g: union([(h["start"], h["end"]) for h in hh if h["cls"] == g and h["identity"] >= 0.8]) for g in set(h["cls"] for h in hh)}
    top = max(cov, key=cov.get)
    if L and cov[top] >= 0.5 * L and top in ("cp", "algmt", "slugmt"):
        best = max((h for h in hh if h["cls"] == top), key=lambda h: h["length"])
        ids = sorted(h["identity"] for h in hh if h["cls"] == top and h["length"] >= 500)
        whole.append(dict(contig=contig, length=L, cls=NAMES[top], covered=cov[top], best_target=best["target"], best_target_desc=desc.get(best["target"], ""),
                          median_identity=ids[len(ids) // 2] if ids else best["identity"]))
    else:
        rest += [h for h in hh if h["cls"] in ("cp", "algmt", "algrDNA", "slugmt")]
# repeat-like panel windows (algal classes only)
tb = collections.Counter()
alg = [h for h in rest if h["cls"] != "slugmt"]
for h in alg:
    for b in range(h["tstart"] // 300, h["tend"] // 300 + 1): tb[(h["target"], b)] += 1
hot = {k for k, v in tb.items() if v >= 15}
keep = []; n_hot = 0
for h in alg:
    bins = [(h["target"], b) for b in range(h["tstart"] // 300, h["tend"] // 300 + 1)]
    if h["length"] < 200 or h["identity"] < 0.80: continue
    if sum(b in hot for b in bins) > 0.5 * len(bins): n_hot += 1; continue
    keep.append(h)
n_low = 0
if keep:
    regs = [f'{h["target"]}:{h["tstart"] + 1}-{h["tend"]}' for h in keep]
    seqs = []
    for i in range(0, len(regs), 500):
        fa = subprocess.run([ST, "faidx", panel_fa] + regs[i:i + 500], capture_output=True, text=True, check=True).stdout
        seqs += ["".join(b.split("\n")[1:]).upper() for b in fa.split(">")[1:]]
    for h, s in zip(keep, seqs): h["target_entropy"] = round(entropy3(s), 2)
    n_low = sum(1 for h in keep if h["target_entropy"] < 4.5)
    keep = sorted((h for h in keep if h["target_entropy"] >= 4.5), key=lambda h: (h["contig"], h["start"]))
loci = []
for h in keep:
    if loci and loci[-1]["contig"] == h["contig"] and h["start"] < loci[-1]["end"] + 300:
        L = loci[-1]; L["end"] = max(L["end"], h["end"]); L["matches"] += 1
        if h["length"] > L["length"]: L.update({k: h[k] for k in ("target", "tstart", "tend", "identity", "length", "cls", "target_entropy")})
    else: loci.append(dict(h, matches=1))
for L in loci:
    L["cls"] = NAMES[L["cls"]]; L["span"] = L["end"] - L["start"]; L["contig_length"] = clen.get(L["contig"], 0); L["target_desc"] = desc.get(L["target"], "")
# ribosomal DNA is conserved across eukaryotes: the slug's own rDNA matches algal rDNA at 80-90%.
# Only matches of >= 92% identity are kept as possibly algal; the rest are counted.
own_rdna = [L for L in loci if L["cls"] == "algal nuclear rDNA" and L["identity"] < 0.92]
loci = [L for L in loci if not (L["cls"] == "algal nuclear rDNA" and L["identity"] < 0.92)]
# slug-mitochondrial matches outside the mitochondrial sequence itself (NUMT-like), merged
numt = []
for h in sorted((h for h in rest if h["cls"] == "slugmt" and h["length"] >= 200 and h["identity"] >= 0.80), key=lambda h: (h["contig"], h["start"])):
    if numt and numt[-1]["contig"] == h["contig"] and h["start"] < numt[-1]["end"] + 300:
        N = numt[-1]; N["end"] = max(N["end"], h["end"])
        if h["length"] > N["length"]: N.update({k: h[k] for k in ("target", "tstart", "tend", "identity", "length")})
    else: numt.append(dict(h))
reps = collections.Counter()
for (t, b), v in tb.items():
    if (t, b) in hot: reps[t] = max(reps[t], v)
res = dict(whole=sorted(whole, key=lambda x: -x["length"]), loci=loci, numt_like=numt,
           removed=dict(repeat_like_matches=n_hot, low_complexity_matches=n_low, repeat_like_panel_sequences=len(reps), own_rdna_loci=len(own_rdna)),
           repeats=[dict(target=t, target_desc=desc.get(t, ""), places=v) for t, v in reps.most_common(15)])
json.dump(res, open(out, "w"), indent=1)
print(f"whole organelle sequences: {collections.Counter(w['cls'] for w in whole)} ({sum(w['length'] for w in whole)} bp)")
print(f"algal-like loci kept: {len(loci)} ({collections.Counter(L['cls'] for L in loci)}); removed {n_hot} repeat-like and {n_low} low-complexity matches, {len(own_rdna)} loci of the slug's own rDNA; NUMT-like {len(numt)}")
for w in sorted(whole, key=lambda x: -x["length"])[:12]:
    print(f"  WHOLE {w['contig']:22s} {w['length']:>9d} {w['cls']:20s} id~{w['median_identity']:.3f} {w['best_target']:24s} {w['best_target_desc'][:50]}")
for L in sorted(loci, key=lambda L: -L["length"] * L["identity"])[:25]:
    print(f"  LOCUS {L['contig']:22s} {L['start']:>10d}-{L['end']:<10d} span {L['span']:6d} match {L['length']:5d} id {L['identity']:.3f} H {L['target_entropy']:.2f} {L['cls'][:12]:12s} {L['target']:24s} {L['tstart']}-{L['tend']} {L['target_desc'][:44]}")
