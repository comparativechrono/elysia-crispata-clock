#!/usr/bin/env python3
"""Helpers shared by the read-screen scripts: alignment segments from SAM text and from PAF."""
import re, subprocess, os
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
ST = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/samtools"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def sam_records(bam):
    """yield (name, flag, target, pos0, cigar, seq, tags) for mapped records of a BAM file"""
    p = subprocess.Popen([ST, "view", bam], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        flag = int(f[1])
        if flag & 4: continue
        yield f[0], flag, f[2], int(f[3]) - 1, f[5], f[9], f[11:]
    p.wait()

def segment(flag, target, pos0, cigar, tags):
    """alignment as a dict with read coordinates on the ORIGINAL read strand"""
    ops = [(int(n), o) for n, o in CIG.findall(cigar)]
    lead = ops[0][0] if ops and ops[0][1] in "SH" else 0
    trail = ops[-1][0] if len(ops) > 1 and ops[-1][1] in "SH" else 0
    qal = sum(n for n, o in ops if o in "MI=X"); tal = sum(n for n, o in ops if o in "MDN=X")
    qlen = lead + qal + trail
    rev = bool(flag & 16)
    qs, qe = (trail, trail + qal) if rev else (lead, lead + qal)
    de = None; nm = None
    for t in tags:
        if t.startswith("de:f:"): de = float(t[5:])
        elif t.startswith("NM:i:"): nm = int(t[5:])
    aln = sum(n for n, o in ops if o in "MID=X")
    if de is None: de = (nm / aln) if (nm is not None and aln) else 0.0
    return dict(qlen=qlen, qs=qs, qe=qe, strand="-" if rev else "+", t=target, ts=pos0, te=pos0 + tal, aln=aln, ident=1 - de,
                cls=target.split("|")[0] if "|" in target else "slug")

def paf_segments(path, minlen=200):
    import collections
    S = collections.defaultdict(list)
    for line in open(path):
        c = line.rstrip("\n").split("\t")
        al = int(c[10])
        if al < minlen: continue
        de = 0.0; cg = None
        for t in c[12:]:
            if t.startswith("de:f:"): de = float(t[5:])
            elif t.startswith("cg:Z:"): cg = t[5:]
        d = dict(qlen=int(c[1]), qs=int(c[2]), qe=int(c[3]), strand=c[4], t=c[5], tlen=int(c[6]), ts=int(c[7]), te=int(c[8]), aln=al, ident=1 - de,
                 mapq=int(c[11]), cls=c[5].split("|")[0] if "|" in c[5] else "slug")
        if cg: d["ins"] = cigar_insertions(cg, d["strand"], d["qs"], d["qe"], d["ts"])
        S[c[0]].append(d)
    return S

def cigar_insertions(cg, strand, qs, qe, ts, minlen=30):
    """read bases that an alignment skips as an insertion of >= minlen: (read_start, read_end, target_position).
    33_pass2.sh keeps the alignment string of an alignment only when it holds such an insertion."""
    out = []; q = 0; t = ts
    for n, o in CIG.findall(cg):
        n = int(n)
        if o in "M=X": q += n; t += n
        elif o == "I":
            if n >= minlen: out.append((qs + q, qs + q + n, t) if strand == "+" else (qe - q - n, qe - q, t))
            q += n
        elif o in "DN": t += n
    return out

def union(iv):
    iv = sorted(iv); tot = 0; cur = None
    for a, b in iv:
        if cur is None: cur = [a, b]
        elif a <= cur[1]: cur[1] = max(cur[1], b)
        else: tot += cur[1] - cur[0]; cur = [a, b]
    if cur: tot += cur[1] - cur[0]
    return tot

def best_tiling(segs, max_overlap=50):
    """greedy choice of non-overlapping segments along the read, best (aligned length x identity) first"""
    kept = []
    for s in sorted(segs, key=lambda s: -s["aln"] * s["ident"]):
        ok = True
        for k in kept:
            ov = min(s["qe"], k["qe"]) - max(s["qs"], k["qs"])
            if ov > max_overlap and ov > 0.2 * (s["qe"] - s["qs"]): ok = False; break
        if ok: kept.append(s)
    return sorted(kept, key=lambda s: s["qs"])

def revcomp(s): return s.translate(COMP)[::-1]


def relabel_organelle_loci(S, screen_json, extra=()):
    """Alignments that fall inside a mitochondrial-like stretch of an assembly sequence
    (35_assembly_screen.py, numt_like) are counted as mitochondrial, not slug: organelle reads
    often match such a stretch as well as they match the organelle genome itself.
    extra: further stretches given by hand as (sequence, start, end, kind), for organelle DNA in an
    assembly that the DNA screen could not recognise because no reference holds it."""
    import json, collections
    if not screen_json and not extra: return 0
    scr = json.load(open(screen_json)) if screen_json else {}; loci = collections.defaultdict(list)
    for c_, a_, b_, k_ in extra: loci[c_].append((a_, b_, k_))
    for N in scr.get("numt_like", []): loci[N["contig"]].append((N["start"], N["end"], "slugmt"))
    for L in scr.get("loci", []):
        if L["cls"] == "plastid" and L["identity"] >= 0.90: loci[L["contig"]].append((L["start"], L["end"], "cp"))
    n = 0
    for ss in S.values():
        for s in ss:
            if s["cls"] != "slug" or s["t"] not in loci: continue
            for a, b, c in loci[s["t"]]:
                ov = min(s["te"], b) - max(s["ts"], a)
                if ov >= 0.8 * (s["te"] - s["ts"]): s["cls"] = c; n += 1; break
    return n
