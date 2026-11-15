#!/usr/bin/env python3
"""What the reads of one species contain, and whether any read joins slug DNA to organelle DNA.

usage: 34_read_analysis.py WORKDIR panel.tsv hifi|clr OUT.json [SCREEN.json]
WORKDIR holds pass2.paf and all_hits.fa (33_pass2.sh) and RUN.count / RUN.select.json files (pass 1).

Each read is described by its best set of non-overlapping alignments (segments). A segment is
  slug      on a sequence of the species' assembly
  cp        on a plastid sequence of the panel, or on an assembly sequence that is plastid DNA as a whole
  slugmt    on an Elysia mitochondrial genome;  algmt  algal mitochondrial;  algrDNA  algal nuclear rDNA
  circ      on a contig of the library's own reads that is a tandem array of a short unit (38c_tandem_contigs.py):
            a small circular element that the amplification of the library copied many times over
Confident segments: >= 300 bp at >= MINID identity (97% for HiFi, 82% for the older, noisier
PacBio RS reads). A read is chimeric when two neighbouring confident segments are of different
kinds and overlap by <= 50 bp on the read. Junctions are clustered (both break points within
25 bp); two reads count as independent molecules only if both of their outer ends differ by more
than 50 bp (PCR copies of one amplified molecule share their ends).
Where the first pass carried on only an even sample of the reads of one kind of organelle DNA
(31_select_hits.py, thin.tsv), a read of the sample counts K times in the composition of the
library, and the share of that kind joined to slug DNA is taken from the sample alone.
With LONG_READ_PLATFORM=ont (Nanopore reads, analysed as kind clr) a junction is tested for the
adapter of the ligation kit, which lies between two molecules that were read as one.
"""
import sys, os, json, glob, collections, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sam_segments import paf_segments, best_tiling, union, relabel_organelle_loci, thinned_reads
W, panel_tsv, kind, out = sys.argv[1:5]
SCREEN = sys.argv[5] if len(sys.argv) > 5 else None
MINID = 0.97 if kind == "hifi" else 0.82
MINSEG = 300 if kind == "hifi" else 500
CONF_LEN = 1000
ADAPTER = "AAGCAGTGGTATCAACGCAGAGTACT"
ONT = os.environ.get("LONG_READ_PLATFORM") == "ont"
ONT_TOP = "AATGTACTTCGTTCAGTTACGTATTGCT"[-22:]; ONT_BOTTOM = "GCAATACGTAACTGAACGAAGT"      # read start (its last 22 bases) and read end
THIN, THIN_K = thinned_reads(W)
def tk(c): return "cp" if c == "circ" else c          # the first pass knew the circular elements as plastid-like reads, and thinned them as such
def wgt(r, c): return THIN_K[tk(c)] if THIN.get(r) == tk(c) else 1
def counted(r, c): return tk(c) not in THIN_K or THIN.get(r) == tk(c)
COMP = str.maketrans("ACGT", "TGCA"); ADAPTER_RC = ADAPTER.translate(COMP)[::-1]
desc = {}
for l in list(open(panel_tsv))[1:]:
    f = l.rstrip("\n").split("\t"); desc[f[0]] = f[4]
def taxon(t):
    d = desc.get(t, "")
    w = d.replace("UNVERIFIED: ", "").replace("UNVERIFIED_ORG: ", "").replace("UNVERIFIED_ASMBLY: ", "").split()
    return " ".join(w[:2]) if len(w) >= 2 else (d or t)
def norm_cls(c): return {"asmcp": "cp", "asmalgmt": "algmt", "kcp": "cp"}.get(c, c)
def kcp_ref(t):
    """A plastid contig assembled from the reads of the library itself (38b_plastid_contigs_long.sh) is named
    kcp|NAME|ACCESSION|IDENTITY after the panel record that it matches best: (that record, identity)."""
    f = t.split("|")
    return ("cp|" + f[2], float(f[3])) if f[0] == "kcp" and len(f) >= 4 else None

runs = []
for c in sorted(glob.glob(os.path.join(W, "*.count"))):
    n, b = (int(x) for x in open(c).read().split()); run = os.path.basename(c)[:-6]
    sel = json.load(open(os.path.join(W, run + ".select.json"))) if os.path.exists(os.path.join(W, run + ".select.json")) else {}
    runs.append(dict(run=run, reads=n, bases=b, hit_reads=sel.get("hit_reads"), reads_only_in_hot_windows=sel.get("reads_only_in_hot_windows"), hot_windows=sel.get("hot_windows", [])[:10]))
tot_reads = sum(r["reads"] for r in runs); tot_bases = sum(r["bases"] for r in runs)

S = paf_segments(os.path.join(W, "pass2.paf"))
relabel_organelle_loci(S, SCREEN)
for ss in S.values():
    for s in ss: s["cls"] = norm_cls(s["cls"])
tiles = {r: best_tiling([s for s in ss if s["aln"] >= MINSEG and s["ident"] >= MINID]) for r, ss in S.items()}
loose = {r: best_tiling([s for s in ss if s["aln"] >= 300 and s["ident"] >= 0.80]) for r, ss in S.items()}

# ---- what the library contains
cls_reads = collections.Counter(); cls_bp = collections.Counter(); main = {}
by_taxon = {"algrDNA": collections.defaultdict(list), "algmt": collections.defaultdict(list)}
for r, kept in tiles.items():
    by = collections.Counter()
    for s in kept:
        if s["cls"] != "slug": by[s["cls"]] += s["qe"] - s["qs"]
        if s["cls"] in by_taxon and s["aln"] >= CONF_LEN: by_taxon[s["cls"]][taxon(s["t"])].append(s["ident"])
    for c, bp in by.items():
        cls_bp[c] += bp * wgt(r, c)
        if any(s["cls"] == c and s["aln"] >= CONF_LEN for s in kept): cls_reads[c] += wgt(r, c)
    conf = [s for s in kept if s["cls"] != "slug" and s["aln"] >= CONF_LEN]
    if conf: main[r] = max(conf, key=lambda s: s["aln"])["cls"]
# plastid DNA by the reference it matches best, at any identity >= 80% (a relative of the true donor still matches)
don = collections.defaultdict(lambda: dict(reads=0, bp=0, ids=[], acc=collections.Counter()))
cp_loose_bp = 0; cp_loose_reads = 0
don_run = collections.defaultdict(collections.Counter)          # the same by run, to compare the runs (or flow cells) of a library
for r, kept in loose.items():
    seen = set()
    for s in kept:
        if s["cls"] != "cp" or s["aln"] < 500: continue
        k = kcp_ref(s["t"])        # reads on a contig of the library's own plastid DNA count for the reference that the contig matches best, at the identity of the contig
        if s["t"].startswith("asmcp|"): tx = "plastid sequence in the assembly"
        elif k: tx = taxon(k[0])
        else: tx = taxon(s["t"])
        w_ = wgt(r, "cp")          # a read of an even sample of the plastid reads stands for K reads
        d = don[tx]; d["bp"] += (s["qe"] - s["qs"]) * w_; d["ids"].append((k[1] if k else s["ident"], s["aln"] * w_)); d["acc"][k[0] if k else s["t"]] += s["aln"] * w_
        don_run[r.split(".")[0]][tx] += (s["qe"] - s["qs"]) * w_
        cp_loose_bp += (s["qe"] - s["qs"]) * w_
        if tx not in seen: d["reads"] += w_; seen.add(tx)
    if seen: cp_loose_reads += wgt(r, "cp")
def wq(pairs, q):
    pairs = sorted(pairs); tot = sum(w for _, w in pairs); acc = 0
    for v, w in pairs:
        acc += w
        if acc >= q * tot: return v
    return pairs[-1][0] if pairs else None
donors = []
for tx, d in sorted(don.items(), key=lambda kv: -kv[1]["bp"]):
    if d["bp"] < 0.0005 * max(1, cp_loose_bp) and len(donors) >= 8: continue
    donors.append(dict(name=tx, reads=d["reads"], bp=d["bp"], fraction=round(d["bp"] / max(1, cp_loose_bp), 5), median_identity=round(wq(d["ids"], 0.5), 4),
                       identity_q10=round(wq(d["ids"], 0.1), 4), identity_q90=round(wq(d["ids"], 0.9), 4), accession=d["acc"].most_common(1)[0][0].split("|", 1)[1]))
donors = donors[:15]

# ---- chimeric reads
need = set(); juncs = []
for r, kept in tiles.items():
    for a, b in zip(kept, kept[1:]):
        if a["cls"] == b["cls"]: continue
        if b["qs"] - a["qe"] < -50: continue
        juncs.append((r, a, b)); need.add(r)
seqs = {}; n = None
for line in open(os.path.join(W, "all_hits.fa")):
    if line[0] == ">": n = line[1:].split()[0]
    elif n in need: seqs[n] = line.strip().upper()
def has_adapter(seq):
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
def edits(pat, text):
    """smallest number of differences between pat and any stretch of text"""
    prev = [0] * (len(text) + 1)
    for i in range(1, len(pat) + 1):
        cur = [i] + [0] * len(text)
        for j in range(1, len(text) + 1):
            cur[j] = min(prev[j - 1] + (pat[i - 1] != text[j - 1]), prev[j] + 1, cur[j - 1] + 1)
        prev = cur
    return min(prev)
def has_ont_adapter(seq): return bool(seq) and (edits(ONT_TOP, seq) <= 5 or edits(ONT_BOTTOM, seq) <= 5)
def bp(seg, side):
    if side == "right": return (seg["te"], seg["ts"]) if seg["strand"] == "+" else (seg["ts"], seg["te"])
    return (seg["ts"], seg["te"]) if seg["strand"] == "+" else (seg["te"], seg["ts"])
rows = []
for r, a, b in juncs:
    ja, oa = bp(a, "right"); jb, ob = bp(b, "left")
    js = seqs.get(r, "")[max(0, a["qe"] - (30 if ONT else 40)):b["qs"] + (30 if ONT else 40)]
    rows.append(dict(read=r, pair="+".join(sorted([a["cls"], b["cls"]])), a_cls=a["cls"], a_t=a["t"], a_bp=ja, a_outer=oa, a_len=a["aln"], a_ident=round(a["ident"], 4), a_mapq=a.get("mapq", 0),
                     b_cls=b["cls"], b_t=b["t"], b_bp=jb, b_outer=ob, b_len=b["aln"], b_ident=round(b["ident"], 4), b_mapq=b.get("mapq", 0), gap=b["qs"] - a["qe"],
                     adapter=has_adapter(js) if kind == "hifi" else has_ont_adapter(js) if ONT else False, n_segments=len(tiles[r])))
def near(a, b, tol): return abs(a - b) <= tol
TOL = 25 if kind == "hifi" else 150
by_t = collections.defaultdict(list)
for row in rows:
    row["_n"] = sorted([(row["a_t"], row["a_bp"], row["a_outer"]), (row["b_t"], row["b_bp"], row["b_outer"])])
    by_t[(row["_n"][0][0], row["_n"][1][0])].append(row)
clusters = []
for rr in by_t.values():
    rr.sort(key=lambda r: (r["_n"][0][1], r["_n"][1][1])); used = [False] * len(rr)
    for i, r in enumerate(rr):
        if used[i]: continue
        grp = [r]; used[i] = True
        for j in range(i + 1, len(rr)):
            if used[j]: continue
            if rr[j]["_n"][0][1] - r["_n"][0][1] > TOL: break
            if near(rr[j]["_n"][0][1], r["_n"][0][1], TOL) and near(rr[j]["_n"][1][1], r["_n"][1][1], TOL): grp.append(rr[j]); used[j] = True
        clusters.append(grp)
def spot(name):
    try: return int(name.rsplit(".", 1)[1])
    except Exception: return None
def same_well(a, b):
    """PacBio RS: the passes of one molecule are loaded as neighbouring reads (and an undetected
    adapter leaves both passes in one read), so reads numbered within 5 of each other are one molecule"""
    if kind != "clr": return a == b
    x, y = spot(a), spot(b)
    return a == b or (x is not None and y is not None and a.rsplit(".", 1)[0] == b.rsplit(".", 1)[0] and abs(x - y) <= 5)
def independent(grp):
    kept = []
    for r in grp:
        if all(not near(r["_n"][0][2], k["_n"][0][2], 50) and not near(r["_n"][1][2], k["_n"][1][2], 50) and not same_well(r["read"], k["read"]) for k in kept): kept.append(r)
    return len(kept)
for grp in clusters:
    k = independent(grp)
    for x in grp: x["cluster_reads"] = len(grp); x["cluster_molecules"] = k
pairs = collections.Counter(); pair_reads = collections.defaultdict(set); adap = collections.Counter(); rec = collections.Counter(); recurrent = []
for row in rows:
    pairs[row["pair"]] += 1; pair_reads[row["pair"]].add(row["read"]); adap[row["pair"]] += row["adapter"]
for grp in clusters:
    if grp[0]["cluster_molecules"] >= 2:
        rec[grp[0]["pair"]] += 1
        recurrent.append(dict(pair=grp[0]["pair"], reads=len(grp), independent_molecules=grp[0]["cluster_molecules"], members=[{k: v for k, v in x.items() if not k.startswith("_")} for x in grp[:12]]))
# reads of each organelle kind that also carry slug DNA (the kind's chimera rate)
chim_by = collections.defaultdict(set)
for row in rows:
    if "slug" in (row["a_cls"], row["b_cls"]):
        other = row["b_cls"] if row["a_cls"] == "slug" else row["a_cls"]
        chim_by[other].add(row["read"])
have = collections.Counter()
for r, kept in tiles.items():
    for c in set(s["cls"] for s in kept):
        if c != "slug" and counted(r, c): have[c] += 1
for c in list(chim_by):
    if tk(c) in THIN_K: chim_by[c] = {r for r in chim_by[c] if counted(r, c)}
rates = {c: dict(reads_with_segment=have[c], reads_also_with_slug=len(chim_by[c]), rate=round(len(chim_by[c]) / have[c], 5) if have[c] else None, **({"sample_one_in": THIN_K[tk(c)]} if tk(c) in THIN_K else {})) for c in have}
# slug - organelle - slug at one slug locus
flank = []
for r, kept in tiles.items():
    for i in range(1, len(kept) - 1):
        l, m, rr = kept[i - 1], kept[i], kept[i + 1]
        if l["cls"] == "slug" and rr["cls"] == "slug" and m["cls"] != "slug" and l["t"] == rr["t"] and l["strand"] == rr["strand"]:
            jl, _ = bp(l, "right"); jr, _ = bp(rr, "left")
            if abs(jr - jl) < 20000:
                flank.append(dict(read=r, slug=l["t"], left_bp=jl, right_bp=jr, insert_cls=m["cls"], insert_t=m["t"], insert_ts=m["ts"], insert_te=m["te"], insert_len=m["qe"] - m["qs"], insert_ident=round(m["ident"], 4)))
with open(os.path.join(W, "junctions.tsv"), "w") as fo:
    cols = ["read", "pair", "n_segments", "a_cls", "a_t", "a_bp", "a_outer", "a_len", "a_ident", "a_mapq", "b_cls", "b_t", "b_bp", "b_outer", "b_len", "b_ident", "b_mapq", "gap", "adapter", "cluster_reads", "cluster_molecules"]
    fo.write("\t".join(cols) + "\n")
    for row in rows: fo.write("\t".join(str(row[c]) for c in cols) + "\n")
other_taxa = {k: [dict(name=t, segments=len(v), mean_identity=round(sum(v) / len(v), 4)) for t, v in sorted(d.items(), key=lambda kv: -len(kv[1]))[:8]] for k, d in by_taxon.items()}
if "circ" in cls_bp and "cp" in THIN_K: THIN_K = dict(THIN_K, circ=THIN_K["cp"])
res = dict(kind=kind, min_identity=MINID, runs=runs, reads=tot_reads, bases=tot_bases, reads_examined=len(S), other_taxa=other_taxa, thinned=THIN_K, adapter_tested="nanopore ligation adapter" if ONT else "amplification adapter" if kind == "hifi" else None,
           classes={c: dict(reads=cls_reads[c], bp=cls_bp[c], fraction_of_bases=round(cls_bp[c] / max(1, tot_bases), 6)) for c in cls_bp},
           plastid_any_identity=dict(reads=cp_loose_reads, bp=cp_loose_bp, fraction_of_bases=round(cp_loose_bp / max(1, tot_bases), 6)),
           donors=donors, donors_by_run={run: dict(bp=sum(v.values()), donors=dict(v.most_common(12))) for run, v in sorted(don_run.items())},
           chimera=dict(pairs=[dict(pair=p, junctions=pairs[p], reads=len(pair_reads[p]), adapter_junctions=adap[p], recurrent_clusters=rec[p]) for p in sorted(pairs, key=lambda p: -pairs[p])],
                        rates=rates, recurrent=sorted(recurrent, key=lambda x: -x["independent_molecules"])[:60], flanked=flank[:60], n_flanked=len(flank)))
json.dump(res, open(out, "w"), indent=1)
print(f"{tot_reads} reads, {tot_bases / 1e9:.2f} Gb screened; {len(S)} reads examined in pass 2")
for c in sorted(cls_bp, key=lambda c: -cls_bp[c]): print(f"  {c:9s} confident reads {cls_reads[c]:7d}  bp {cls_bp[c]:12d}  ({100 * cls_bp[c] / max(1, tot_bases):.4f}% of bases)")
print(f"  plastid at any identity >= 80%: {cp_loose_reads} reads, {cp_loose_bp} bp ({100 * cp_loose_bp / max(1, tot_bases):.4f}%)")
for d in donors[:10]: print(f"    {d['name']:38s} {d['reads']:7d} reads {d['bp']:11d} bp {100 * d['fraction']:6.2f}%  identity median {d['median_identity']:.3f} (10-90%: {d['identity_q10']:.3f}-{d['identity_q90']:.3f}) {d['accession']}")
print("  chimeric junctions:", {p: pairs[p] for p in pairs}, "| recurrent clusters:", dict(rec), "| flanked inserts:", collections.Counter(f["insert_cls"] for f in flank))
print("  rates:", {c: (v["reads_also_with_slug"], v["reads_with_segment"], v["rate"]) for c, v in rates.items()})
