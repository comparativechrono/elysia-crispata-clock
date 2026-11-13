#!/usr/bin/env python3
"""Find HiFi reads that join slug nuclear DNA to plastid, bacterial or mitochondrial DNA.

usage: 14_chimera_analysis.py pass2.paf reads.fa read_classes.tsv outdir [extra_slug.paf]
pass2.paf : minimap2 -x map-hifi -c alignments of the selected reads against
            slug primary + alternate assemblies + the (unmasked) bait panel.
extra_slug.paf (optional): alignments of the same reads against the Purdue and Tampa
            assemblies, used only to recognise slug sequence missing from the Sanger assemblies.

A read is called chimeric when it carries a slug segment and a non-slug segment (each >= 300 bp
at >= 97 % identity) that overlap by <= 50 bp on the read.
Junctions are clustered (both breakpoints within 20 bp); reads sharing both outer ends
(within 20 bp) are treated as PCR duplicates of one molecule.
"""
import sys, os, collections, json

paf, reads_fa, classes_tsv, out = sys.argv[1:5]
extra = sys.argv[5] if len(sys.argv) > 5 else None
os.makedirs(out, exist_ok=True)
ADAPTER = "AAGCAGTGGTATCAACGCAGAGTACT"
COMP = str.maketrans("ACGT", "TGCA")
ADAPTER_RC = ADAPTER.translate(COMP)[::-1]

def cls_of(t):
    if "|" in t:
        c = t.split("|")[0]
        return "bac" if c.startswith("bac_") else c
    return "slug"

def load(p, force=None):
    segs = collections.defaultdict(list)
    for line in open(p):
        c = line.rstrip("\n").split("\t")
        al = int(c[10])
        if al < 200: continue
        de = 0.0
        for t in c[12:]:
            if t.startswith("de:f:"): de = float(t[5:]); break
        segs[c[0]].append(dict(qlen=int(c[1]), qs=int(c[2]), qe=int(c[3]), strand=c[4], t=c[5], ts=int(c[7]), te=int(c[8]),
                               aln=al, ident=1 - de, mapq=int(c[11]), cls=force or cls_of(c[5])))
    return segs

S = load(paf)
if extra:
    for r, ss in load(extra, force="slug_other").items(): S[r].extend(ss)
# class of each read from pass 1, kept only when the match is confident (<= 3 % divergence over >= 1 kb);
# these counts are the denominators of the chimera rates
main = {}
for line in list(open(classes_tsv))[1:]:
    f = line.rstrip("\n").split("\t")
    if float(f[6]) <= 0.03 and int(f[7]) >= 1000: main[f[0]] = f[2]

def has_adapter(seq):
    """adapter (<= 3 mismatches over 20 nt) anywhere in seq"""
    for a in (ADAPTER, ADAPTER_RC):
        core = a[3:23]
        for i in range(len(seq) - len(core) + 1):
            mm = 0
            for x, y in zip(seq[i:i + len(core)], core):
                if x != y:
                    mm += 1
                    if mm > 3: break
            if mm <= 3: return True
    return False

juncs = []; per_read = {}
for r, ss in S.items():
    good = [s for s in ss if s["aln"] >= 300 and ((s["cls"] in ("slug", "slug_other") and s["ident"] >= 0.97) or (s["cls"] not in ("slug", "slug_other") and s["ident"] >= 0.97))]
    # one best segment per read interval: greedy by aligned length, drop segments mostly covered by a better one
    good.sort(key=lambda s: -s["aln"] * s["ident"])
    kept = []
    for s in good:
        ok = True
        for k in kept:
            ov = min(s["qe"], k["qe"]) - max(s["qs"], k["qs"])
            if ov > 50 and ov > 0.2 * (s["qe"] - s["qs"]): ok = False; break
        if ok: kept.append(s)
    kept.sort(key=lambda s: s["qs"])
    per_read[r] = kept
    for a, b in zip(kept, kept[1:]):
        ca = "slug" if a["cls"].startswith("slug") else a["cls"]; cb = "slug" if b["cls"].startswith("slug") else b["cls"]
        if ca == cb: continue
        juncs.append((r, a, b))

need = {r for r, a, b in juncs}
seqs = {}; n = None; buf = []
for line in open(reads_fa):
    if line[0] == ">":
        if n in need: seqs[n] = "".join(buf).upper()
        n, buf = line[1:].split()[0], []
    elif n in need:
        buf.append(line.strip())
if n in need: seqs[n] = "".join(buf).upper()

def bp(seg, side):
    """breakpoint on the target for the read-side 'side' (left end = qs, right end = qe) and the outer end"""
    if side == "right":   # junction at qe
        return (seg["te"], seg["ts"]) if seg["strand"] == "+" else (seg["ts"], seg["te"])
    return (seg["ts"], seg["te"]) if seg["strand"] == "+" else (seg["te"], seg["ts"])

rows = []
for r, a, b in juncs:
    ja, oa = bp(a, "right"); jb, ob = bp(b, "left")
    gap = b["qs"] - a["qe"]
    js = seqs[r][max(0, a["qe"] - 40):b["qs"] + 40]
    pair = "+".join(sorted(["slug" if x["cls"].startswith("slug") else x["cls"] for x in (a, b)]))
    rows.append(dict(read=r, pair=pair, main=main.get(r, "?"), a_cls=a["cls"], a_t=a["t"], a_bp=ja, a_outer=oa, a_strand=a["strand"], a_len=a["aln"], a_ident=round(a["ident"], 4), a_mapq=a["mapq"],
                     b_cls=b["cls"], b_t=b["t"], b_bp=jb, b_outer=ob, b_strand=b["strand"], b_len=b["aln"], b_ident=round(b["ident"], 4), b_mapq=b["mapq"],
                     gap=gap, adapter=has_adapter(js), n_segments=len(per_read[r])))
# cluster junctions (both breakpoints within 25 bp) and count independent molecules:
# two reads count as independent only if BOTH of their outer ends differ by more than 50 bp
# (PCR duplicates of one amplified molecule share at least one end).
def near(a, b, tol): return abs(a - b) <= tol
def norm(row):
    x = sorted([(row["a_t"], row["a_bp"], row["a_outer"]), (row["b_t"], row["b_bp"], row["b_outer"])])
    return x
by_targets = collections.defaultdict(list)
for row in rows:
    n = norm(row); row["_n"] = n
    by_targets[(n[0][0], n[1][0])].append(row)
clusters = []
for key_, rr in by_targets.items():
    rr.sort(key=lambda r: (r["_n"][0][1], r["_n"][1][1]))
    used = [False] * len(rr)
    for i, r in enumerate(rr):
        if used[i]: continue
        grp = [r]; used[i] = True
        for j in range(i + 1, len(rr)):
            if used[j]: continue
            if rr[j]["_n"][0][1] - r["_n"][0][1] > 25: break
            if near(rr[j]["_n"][0][1], r["_n"][0][1], 25) and near(rr[j]["_n"][1][1], r["_n"][1][1], 25):
                grp.append(rr[j]); used[j] = True
        clusters.append(grp)
def independent(grp):
    """greedy count of reads whose two outer ends both differ from every read already counted"""
    kept = []
    for r in grp:
        if all(not near(r["_n"][0][2], k["_n"][0][2], 50) and not near(r["_n"][1][2], k["_n"][1][2], 50) for k in kept):
            kept.append(r)
    return len(kept)
for grp in clusters:
    k = independent(grp)
    for x in grp:
        x["cluster_reads"] = len(grp); x["cluster_molecules"] = k
def mol(row, w=20):
    return tuple(sorted([(row["a_t"], row["a_outer"] // w), (row["b_t"], row["b_outer"] // w)]))
cols = ["read", "pair", "main", "n_segments", "a_cls", "a_t", "a_bp", "a_outer", "a_strand", "a_len", "a_ident", "a_mapq", "b_cls", "b_t", "b_bp", "b_outer", "b_strand", "b_len", "b_ident", "b_mapq", "gap", "adapter", "cluster_reads", "cluster_molecules"]
with open(f"{out}/junctions.tsv", "w") as fo:
    fo.write("\t".join(cols) + "\n")
    for row in rows: fo.write("\t".join(str(row[c]) for c in cols) + "\n")
# flanked insertions: slug - X - slug with both slug flanks on the same sequence, same strand, colinear and close
flank = []
for r, kept in per_read.items():
    for i in range(1, len(kept) - 1):
        l, m, rr = kept[i - 1], kept[i], kept[i + 1]
        if l["cls"].startswith("slug") and rr["cls"].startswith("slug") and not m["cls"].startswith("slug") and l["t"] == rr["t"] and l["strand"] == rr["strand"]:
            jl, _ = bp(l, "right"); jr, _ = bp(rr, "left")
            if abs(jr - jl) < 20000:
                flank.append(dict(read=r, slug=l["t"], left_bp=jl, right_bp=jr, insert_cls=m["cls"], insert_t=m["t"], insert_ts=m["ts"], insert_te=m["te"], insert_len=m["qe"] - m["qs"], insert_ident=round(m["ident"], 4)))
json.dump(flank, open(f"{out}/flanked_inserts.json", "w"), indent=1)
# summaries
tot_main = collections.Counter(main.values())
summ = collections.OrderedDict()
pairs = collections.Counter(); pairs_reads = collections.defaultdict(set); pairs_mol = collections.defaultdict(set); adap = collections.Counter(); multi = collections.Counter()
for row in rows:
    pairs[row["pair"]] += 1; pairs_reads[row["pair"]].add(row["read"]); pairs_mol[row["pair"]].add(mol(row))
    adap[row["pair"]] += row["adapter"]
recurrent = []
for grp in clusters:
    p = grp[0]["pair"]
    if grp[0]["cluster_molecules"] >= 2:
        multi[p] += 1
        recurrent.append(dict(pair=p, reads=len(grp), independent_molecules=grp[0]["cluster_molecules"],
                              members=[{k: v for k, v in x.items() if not k.startswith("_")} for x in grp]))
print("class totals from pass 1 (reads by main class):", dict(tot_main))
print(f"{'pair':22s} {'junctions':>9s} {'reads':>7s} {'molecules':>9s} {'adapter@junction':>17s} {'junctions seen in >=2 independent molecules':>28s}")
res = []
for p in sorted(pairs, key=lambda p: -pairs[p]):
    print(f"{p:22s} {pairs[p]:9d} {len(pairs_reads[p]):7d} {len(pairs_mol[p]):9d} {adap[p]:10d} ({100 * adap[p] / pairs[p]:4.1f}%) {multi[p]:20d}")
    res.append(dict(pair=p, junctions=pairs[p], reads=len(pairs_reads[p]), molecules=len(pairs_mol[p]), adapter_junctions=adap[p], recurrent_clusters=multi[p]))
# how many examined reads of each class carry a slug segment at all
chim_reads = collections.defaultdict(set)
for row in rows:
    for side in ("a", "b"):
        c = row[side + "_cls"]
        if not c.startswith("slug"): chim_reads[(c, row["pair"])].add(row["read"])
json.dump(dict(pairs=res, class_totals=dict(tot_main), flanked=len(flank), recurrent=recurrent,
               reads_examined=len(S)), open(f"{out}/chimera_summary.json", "w"), indent=1)
print("flanked inserts (slug-X-slug, same locus):", collections.Counter(f["insert_cls"] for f in flank))
