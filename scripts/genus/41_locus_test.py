#!/usr/bin/env python3
"""Is an organelle-like stretch inside an assembly sequence really joined to its flanks?
Test each one with reads.

usage: 41_locus_test.py baits.fa.tsv GENOME.fa.gz READS.fa clr|hifi|short OUT.json
baits.fa.tsv (40_locus_baits.py) lists the loci and a window of 3 kb either side. The reads are
aligned base by base to the windows (real sequence, nothing masked). For each of the two edges
of a locus the script counts
  through   alignments that run at least M bases past the edge on both sides: molecules in
            which the organelle-like stretch and the flank are one piece of DNA. An alignment that
            crosses the edge but has a gap of 15 bases or more (6 for accurate reads) within 80
            bases of it is counted apart (through_with_gap): that molecule joins the two kinds of
            DNA at a point of its own, not at the one in the assembly.
  stop      alignments on one side that end within 20 bases of the edge while the read goes on
            for another M bases or more: molecules that continue into something else there
            (inside = the read lies in the locus; outside = it lies in the flank)
and the mean depth of the locus and of 500 bases of flank on each side.
M is 100 for long reads and 25 for short reads. Alignments need >= 80% identity (long, noisy
reads) or >= 97% (accurate long reads, short reads).
"""
import sys, os, json, subprocess, collections
tsv, genome, reads, kind, out = sys.argv[1:6]
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/"
M = 25 if kind == "short" else 100
MINID = 0.80 if kind == "clr" else 0.97
PRESET = {"clr": "map-pb", "hifi": "map-hifi", "short": "sr"}[kind]
loci = []
for i, l in enumerate(open(tsv)):
    if i == 0: continue
    f = l.rstrip("\n").split("\t"); loci.append(dict(name=f[0], contig=f[1], ws=int(f[2]), we=int(f[3]), kind=f[4], start=int(f[5]), end=int(f[6])))
wd = out + ".tmp"; os.makedirs(wd, exist_ok=True)
with open(f"{wd}/windows.fa", "w") as fo:
    for i, L in enumerate(loci):
        fa = subprocess.run([BIN + "samtools", "faidx", genome, f'{L["contig"]}:{L["ws"] + 1}-{L["we"]}'], capture_output=True, text=True, check=True).stdout
        seq = "".join(fa.split("\n")[1:]).upper(); L["wlen"] = len(seq); fo.write(f">w{i}\n{seq}\n")
subprocess.run(f"{BIN}minimap2 -x {PRESET} -c --secondary=no -t 2 {wd}/windows.fa {reads} 2>/dev/null > {wd}/aln.paf", shell=True, check=True)
import re
GAP = 15 if kind == "clr" else 6
def gaps(cg, ts):
    """target positions and lengths of the insertions and deletions of GAP bases or more in an alignment"""
    out = []; t = ts
    for n, o in re.findall(r"(\d+)([MIDNSHP=X])", cg):
        n = int(n)
        if o in "M=X": t += n
        elif o in "DN":
            if n >= GAP: out.append((t, t + n, n))
            t += n
        elif o == "I" and n >= GAP: out.append((t, t, n))
    return out
A = collections.defaultdict(list)
for l in open(f"{wd}/aln.paf"):
    c = l.split("\t"); al = int(c[10])
    if al < (60 if kind == "short" else 300) or int(c[9]) / al < MINID: continue
    qlen, qs, qe = int(c[1]), int(c[2]), int(c[3])
    left_clip, right_clip = (qs, qlen - qe) if c[4] == "+" else (qlen - qe, qs)     # unaligned read bases beyond the target-left / target-right end
    cg = next((t[5:] for t in c[12:] if t.startswith("cg:Z:")), "")
    A[int(c[5][1:])].append(dict(read=c[0], ts=int(c[7]), te=int(c[8]), lclip=left_clip, rclip=right_clip, gaps=gaps(cg.strip(), int(c[7]))))
def molecules(names):
    """long reads: the several passes of one molecule are separate reads with neighbouring numbers (within 15); short reads: the two mates of a pair"""
    if kind == "short": return len({n.rsplit(".", 1)[0] for n in names})
    spots = sorted(int(n.rsplit(".", 1)[1]) for n in names if n.rsplit(".", 1)[1].isdigit()); k = len(names) - len(spots)
    for i, x in enumerate(spots):
        if i == 0 or x - spots[i - 1] > 15: k += 1
    return k
res = []
for i, L in enumerate(loci):
    a = L["start"] - L["ws"]; b = L["end"] - L["ws"]; al = A.get(i, [])
    def depth(x, y):
        if y <= x: return None
        return round(sum(max(0, min(r["te"], y) - max(r["ts"], x)) for r in al) / (y - x), 2)
    edges = {}
    for name, e, inside_is_right in (("left", a, True), ("right", b, False)):
        if e < M or e > L["wlen"] - M:
            edges[name] = dict(note="the locus reaches the end of the assembly sequence here; there is no flank to test"); continue
        crossing = [r for r in al if r["ts"] <= e - M and r["te"] >= e + M]
        def gap_at(r): return max((n for x, y, n in r["gaps"] if x <= e + 80 and y >= e - 80), default=0)
        gapped = {r["read"]: gap_at(r) for r in crossing if gap_at(r)}
        through = {r["read"] for r in crossing} - set(gapped)
        # reads that stop at the edge: coming from the right (start at e, clipped on their left) or from the left (end at e, clipped on their right)
        from_right = {r["read"] for r in al if abs(r["ts"] - e) <= 20 and r["lclip"] >= M} - through - set(gapped)
        from_left = {r["read"] for r in al if abs(r["te"] - e) <= 20 and r["rclip"] >= M} - through - set(gapped)
        inside, outside = (from_right, from_left) if inside_is_right else (from_left, from_right)
        edges[name] = dict(position=L["ws"] + e, through=len(through), through_reads=sorted(through)[:20], through_with_gap=len(gapped), through_with_gap_reads=dict(sorted(gapped.items())[:20]),
                           through_molecules=molecules(through), through_with_gap_molecules=molecules(gapped),
                           stop_inside=len(inside), stop_outside=len(outside))
    res.append(dict(name=L["name"], contig=L["contig"], start=L["start"], end=L["end"], kind=L["kind"], reads_aligned=len({r["read"] for r in al}),
                    depth_locus=depth(a, b), depth_left_flank=depth(max(0, a - 500), a), depth_right_flank=depth(b, min(L["wlen"], b + 500)), edges=edges))
    print(f'{L["contig"]}:{L["start"]}-{L["end"]} ({L["kind"]}) reads {res[-1]["reads_aligned"]} depth locus {res[-1]["depth_locus"]} flanks {res[-1]["depth_left_flank"]} / {res[-1]["depth_right_flank"]} | ' +
          " | ".join(f'{k}: ' + (v.get("note", "")[:28] or f'through {v["through"]} (+{v["through_with_gap"]} with a gap), stop inside {v["stop_inside"]}, stop outside {v["stop_outside"]}') for k, v in edges.items()))
json.dump(dict(kind=kind, margin=M, loci=res), open(out, "w"), indent=1)
