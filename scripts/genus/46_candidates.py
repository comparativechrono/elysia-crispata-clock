#!/usr/bin/env python3
"""Long reads: every molecule that joins organelle DNA to slug DNA, sorted by its shape.

usage: 46_candidates.py WORKDIR hifi|clr OUT.json [SCREEN.json] [--paf FILE] [--no-identity] [--no-p1]
                        [--genome-bp N] [--one-in S]
  --paf FILE     alignments to use instead of WORKDIR/pass2.paf: the random sample of the library
                 (sample.paf of 33_pass2.sh) or all its reads (47_library_control.sh)
  --genome-bp N  size of the assembly: gives the depth of the library and, on a sample, the depth of
                 reads long enough to show an insertion with both its flanks
  --one-in S     the reads of --paf are one in S of the library (default: from RUN.select.json)
  --panel FILE   the reference panel (indexed FASTA). A short piece that only the first pass matched
                 to an organelle sequence, by a chain of seeds, is then aligned to that sequence
                 base by base, and kept only if 80% of it aligns at 95% identity or better (75% for
                 the older reads). Without this, a tandem repeat of the slug that shares one seed
                 per unit with a tandem repeat of a plastid genome passes for a piece of that genome.
  --reads FILE   the reads (default WORKDIR/all_hits.fa), needed for --panel
  --depth-paf F  alignments of all reads of the library to the assembly (47_library_control.sh):
                 for every insertion-shaped read, the number of other reads that cross the same
                 place in one piece, without the insert

A somatic insertion of chloroplast DNA sits in the nuclei of some cells only, so it is carried by
few molecules, perhaps one, and the other reads of that place do not show it. One read is
therefore all the evidence there can be, and the question is what else could have made that read.
Each read is described by its best set of non-overlapping alignments to the species' assembly and
the reference panel (33_pass2.sh), and each organelle piece that touches slug DNA is given a shape:

  insertion   slug | organelle | slug, where the two slug pieces lie on the same assembly sequence,
              on the same strand and in the right order, no more than 5 kb apart: the organelle DNA
              sits inside one continuous stretch of slug DNA. Joining two unrelated molecules in
              the tube cannot give this; it would take two joins that happen to pick neighbours.
              Found in two ways:
                between alignments    the organelle piece is long enough to break the alignment of
                                      the read to the slug sequence in two (about 400 bases or more)
                inside one alignment  a shorter piece: the read aligns to the slug sequence across
                                      it, with the piece as an insertion of 100 bases or more (150 for
                                      the older reads) in that alignment, and the first pass (RUN.p1.paf)
                                      or the second matches the inserted bases to an organelle sequence.
                                      The same organelle stretch found inserted at two different places
                                      of the genome is a mobile element that happens to resemble it, and
                                      is counted apart ("element")
  fold-back   slug | organelle | slug where the second slug piece is the first one read backwards
              (a molecule read on both strands): one join, seen twice
  one join    slug | organelle, with nothing recognisable beyond the organelle piece, or with
              unrelated slug DNA beyond it. A real junction and a join made in the tube look the
              same in a single read.

The same shapes are counted with other DNA in the middle, as a measure of what else makes them:
  mitochondrial DNA between two neighbouring slug pieces
  slug DNA from elsewhere in the genome between two neighbouring slug pieces
  slug DNA between two neighbouring pieces of the plastid genome (the mirror image)
and reads made of two unrelated slug pieces give the rate at which the library joins molecules.
These counts mean most on all reads of a library or a random sample of it (47_library_control.sh),
not on the reads chosen by the first pass, nearly all of which hold organelle DNA.

An organelle piece that the read shares with the assembly at that place (the slug alignment covers
it without an insertion) is not an insertion in this animal: it is counted as "shared".
"""
import sys, os, json, glob, collections, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sam_segments import paf_segments, best_tiling, relabel_organelle_loci
ap = argparse.ArgumentParser()
ap.add_argument("workdir"); ap.add_argument("kind", choices=["hifi", "clr"]); ap.add_argument("out"); ap.add_argument("screen", nargs="?")
ap.add_argument("--paf"); ap.add_argument("--no-identity", action="store_true"); ap.add_argument("--no-p1", action="store_true")
ap.add_argument("--genome-bp", type=int, default=0); ap.add_argument("--one-in", type=int, default=0); ap.add_argument("--depth-paf")
ap.add_argument("--panel"); ap.add_argument("--reads")
A = ap.parse_args()
W, kind = A.workdir, A.kind
MINID = 0.0 if A.no_identity else (0.97 if kind == "hifi" else 0.82)
MINSEG = 300 if kind == "hifi" else 500
MINSHORT = 100 if kind == "hifi" else 150  # shortest organelle piece looked for inside an alignment. Below about 100 bases the
                                           # matches are transfer RNA genes, which mobile elements of the slug genome carry too
MAXGAP, MINGAP = 5000, -300          # reference distance allowed between the two flanks of an insertion
FAR = 50000                          # "elsewhere": another sequence, or further than this
UNIQ = 30                            # mapping quality from which a piece of slug DNA counts as placed
MINSEEDED = 0.3 if kind == "hifi" else 0.2   # first-pass match: share of its length covered by matching seeds
ORG = ("cp", "slugmt", "algmt")
def norm_cls(c): return {"asmcp": "cp", "asmalgmt": "algmt", "kcp": "cp"}.get(c, c)

S = paf_segments(A.paf or os.path.join(W, "pass2.paf"), minlen=MINSHORT)
relabel_organelle_loci(S, A.screen)
embedded = collections.defaultdict(list)          # places where the assembly itself holds organelle-like DNA
if A.screen:
    scr = json.load(open(A.screen))
    for L in scr.get("loci", []):
        if L["cls"] == "plastid": embedded[L["contig"]].append((L["start"], L["end"], "plastid-like stretch of the assembly"))
    for N in scr.get("numt_like", []): embedded[N["contig"]].append((N["start"], N["end"], "mitochondrial-like stretch of the assembly"))
def at_embedded(t, a, b):
    for s, e, what in embedded.get(t, []):
        if min(b, e + 100) - max(a, s - 100) > 0: return what
    return None

tiles = {}; short_org = collections.defaultdict(list)
for r, ss in S.items():
    for s in ss: s["cls"] = norm_cls(s["cls"])
    tiles[r] = best_tiling([s for s in ss if s["aln"] >= MINSEG and s["ident"] >= MINID])
    for s in ss:
        if s["cls"] in ORG and s["aln"] < MINSEG and s["ident"] >= MINID:
            short_org[r].append(dict(qs=s["qs"], qe=s["qe"], cls=s["cls"], t=s["t"], ts=s["ts"], te=s["te"], strand=s["strand"], ident=round(s["ident"], 4), source="second pass"))
if not A.no_p1 and not A.paf:
    for f in sorted(glob.glob(os.path.join(W, "*.p1.paf"))):
        for line in open(f):
            c = line.rstrip("\n").split("\t"); cls = norm_cls(c[5].split("|")[0])
            if cls not in ORG or c[0] not in tiles or int(c[3]) - int(c[2]) < MINSHORT: continue
            # a real piece is covered by matching seeds; a chain of a few seeds strung over hundreds of bases is a
            # repeat of the read that resembles a repeat of the panel sequence (a minisatellite, for one)
            if int(c[9]) < MINSEEDED * int(c[10]): continue
            dv = next((float(x[5:]) for x in c[12:] if x.startswith("dv:f:")), None)
            short_org[c[0]].append(dict(qs=int(c[2]), qe=int(c[3]), cls=cls, t=c[5], ts=int(c[7]), te=int(c[8]), strand=c[4],
                                        ident=None if dv is None else round(1 - dv, 4), source="first pass"))

# ---- base-level confirmation of the short pieces that come from the first pass
READSEQ = {}; unconfirmed = []
if A.panel and any(o["source"] == "first pass" for v in short_org.values() for o in v):
    need = {r for r, v in short_org.items() if any(o["source"] == "first pass" for o in v)}; name = None
    for line in open(A.reads or os.path.join(W, "all_hits.fa")):
        if line[0] == ">": name = line[1:].split()[0]
        elif name in need: READSEQ[name] = READSEQ.get(name, "") + line.strip().upper()
    import subprocess
    from Bio import Align
    _AL = Align.PairwiseAligner(); _AL.mode = "local"; _AL.match_score = 1; _AL.mismatch_score = -1; _AL.open_gap_score = -2; _AL.extend_gap_score = -1
    _COMP = str.maketrans("ACGTN", "TGCAN"); _BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin")
    MINCONF = 0.95 if kind == "hifi" else 0.75
def confirm(r, O):
    """(confirmed?, identity, share of the piece aligned) for a first-pass piece O of read r; None if it cannot be tested"""
    if not A.panel or r not in READSEQ: return None, None, None
    q = READSEQ[r][O["qs"]:O["qe"]]
    fa = subprocess.run([_BIN + "/samtools", "faidx", A.panel, f'{O["t"]}:{max(1, O["ts"] - 40)}-{O["te"] + 40}'], capture_output=True, text=True).stdout
    tseq = "".join(fa.split("\n")[1:]).upper()
    if O["strand"] == "-": tseq = tseq.translate(_COMP)[::-1]
    if not q or not tseq: return None, None, None
    a = _AL.align(q, tseq)[0]; qb, tb = a.aligned
    if not len(qb): return False, 0.0, 0.0
    m = sum(1 for (x0, x1), (y0, y1) in zip(qb, tb) for k in range(int(x1 - x0)) if q[int(x0) + k] == tseq[int(y0) + k])
    cols = sum(int(x1 - x0) for x0, x1 in qb) + sum(max(int(qb[i + 1][0] - qb[i][1]), int(tb[i + 1][0] - tb[i][1])) for i in range(len(qb) - 1))
    share = (int(qb[-1][1]) - int(qb[0][0])) / len(q); pid = m / cols if cols else 0.0
    return bool(share >= 0.8 and pid >= MINCONF), round(pid, 4), round(share, 3)

def ref_gap(L, R):
    """reference bases between two pieces that follow each other on the read; None if not on the same sequence and strand"""
    if L["t"] != R["t"] or L["strand"] != R["strand"]: return None
    return (R["ts"] - L["te"]) if L["strand"] == "+" else (L["ts"] - R["te"])
def neighbours(L, R):
    g = ref_gap(L, R)
    return g is not None and MINGAP <= g <= MAXGAP
def folded(L, R):
    """R is L read backwards: same sequence, opposite strand, overlapping or within 2 kb"""
    return L["t"] == R["t"] and L["strand"] != R["strand"] and min(L["te"], R["te"]) - max(L["ts"], R["ts"]) > -2000
def elsewhere(X, L):
    return X["t"] != L["t"] or min(abs(X["ts"] - L["te"]), abs(X["te"] - L["ts"])) > FAR
def seg_out(s): return dict(cls=s["cls"], t=s["t"], ts=s["ts"], te=s["te"], strand=s["strand"], qs=s["qs"], qe=s["qe"], aln=s["aln"], ident=round(s["ident"], 4), mapq=s.get("mapq", 0))
def place(L, R):
    """the stretch of the slug sequence between the two flanks"""
    a, b = (L["te"], R["ts"]) if L["strand"] == "+" else (R["te"], L["ts"])
    return min(a, b), max(a, b)

def second_pass(t, i, j):
    """A read that turns round and reads the same molecule on its other strand (a missed adapter, or a
    molecule that folds back) shows the place twice. Does the other strand hold the insert as well?
      'insert on both strands'   the middle piece is read again, on the other strand
      'insert on one strand'     the flank is read again and runs on, into the slug sequence, across
                                 the point where the first reading turned into the insert
    None when the read does not come back over the place."""
    L, R = t[i], t[j]; mid = t[i + 1:j]; out = None
    for flank, others, adjacent in ((R, t[j + 1:], R["qs"] - mid[-1]["qe"]), (L, t[:i], mid[0]["qs"] - L["qe"])):
        for x in others:
            if any(x["t"] == m["t"] and x["strand"] != m["strand"] and min(x["te"], m["te"]) - max(x["ts"], m["ts"]) > 0.5 * min(x["te"] - x["ts"], m["te"] - m["ts"]) for m in mid):
                return "insert on both strands"
            if x["t"] != flank["t"] or x["strand"] == flank["strand"] or adjacent > 60: continue
            # the end of the flank that touches the insert, on the reference
            inner_is_start = (flank is R) == (flank["strand"] == "+")
            beyond = (flank["ts"] - x["ts"]) if inner_is_start else (x["te"] - flank["te"])
            if beyond >= 80 and min(x["te"], flank["te"]) - max(x["ts"], flank["ts"]) > 100: out = "insert on one strand"
    return out

shape_n = collections.Counter(); cand = []; reads_with = collections.Counter(); nuclear = collections.Counter(); nuc_rows = []; shared = collections.Counter(); strands = collections.Counter()
for r, t in tiles.items():
    kinds = {s["cls"] for s in t}
    for k in kinds: reads_with[k] += 1
    if not t: continue
    n = len(t)
    # ---- insertion shapes between alignments: flanks i and j, everything between them of one other kind
    used = set(); found = []
    for i in range(n - 2):
        for j in range(i + 2, min(n, i + 6)):
            L, R = t[i], t[j]; mid = t[i + 1:j]
            if L["cls"] != R["cls"]: continue
            mk = {m["cls"] for m in mid}
            if L["cls"] == "slug" and neighbours(L, R):
                if mk <= set(ORG) and len(mk) == 1: found.append(("insertion", next(iter(mk)), i, j))
                elif mk == {"slug"} and all(elsewhere(m, L) and elsewhere(m, R) for m in mid): found.append(("insertion", "slug_elsewhere", i, j))
            elif L["cls"] == "cp" and neighbours(L, R) and mk == {"slug"}: found.append(("mirror", "slug_in_cp", i, j))
            elif L["cls"] == "slugmt" and neighbours(L, R) and mk == {"slug"}: found.append(("mirror", "slug_in_slugmt", i, j))
    for shape, what, i, j in found:
        L, R = t[i], t[j]; mid = t[i + 1:j]; used.update(range(i, j + 1))
        g = ref_gap(L, R); pa, pb = place(L, R)
        row = dict(read=r, read_length=t[0]["qlen"], shape=shape, insert=what, how="between alignments", left=seg_out(L), right=seg_out(R), middle=[seg_out(m) for m in mid], reference_gap=g,
                   read_gap=R["qs"] - L["qe"], insert_bp=sum(m["qe"] - m["qs"] for m in mid), unaligned_left=mid[0]["qs"] - L["qe"], unaligned_right=R["qs"] - mid[-1]["qe"],
                   n_segments=n, folds_back_after=any(folded(R, x) for x in t[j + 1:]) or any(folded(L, x) for x in t[:i]),
                   assembly_locus=at_embedded(L["t"], pa, pb) if L["cls"] == "slug" else None, second_pass=second_pass(t, i, j))
        shape_n[(shape, what)] += 1
        if row["second_pass"]: strands[(what, row["second_pass"])] += 1
        if what == "slug_elsewhere":
            row["insert_mapq_min"] = min(m.get("mapq", 0) for m in mid); nuc_rows.append(row)
        else: cand.append(row)
    # ---- shorter organelle pieces: inside an insertion of one slug alignment, or between two neighbouring ones
    if r in short_org:
        big = [(s["qs"], s["qe"]) for s in t if s["cls"] in ORG]
        slugs = [(k, s) for k, s in enumerate(t) if s["cls"] == "slug"]
        done_here = []
        for O in sorted(short_org[r], key=lambda o: -(o["qe"] - o["qs"])):
            ln = O["qe"] - O["qs"]
            if any(min(b, O["qe"]) - max(a, O["qs"]) > 0.5 * ln for a, b in big): continue          # already a piece of the tiling
            if any(min(b, O["qe"]) - max(a, O["qs"]) > 0.5 * min(ln, b - a) and c == O["cls"] for a, b, c in done_here): continue
            hit = None
            for k, s in slugs:
                for a, b, tp in s.get("ins", []):
                    ov = min(b, O["qe"]) - max(a, O["qs"])
                    if b - a >= MINSHORT and ov >= MINSHORT and ov >= 0.6 * ln and a - s["qs"] >= 100 and s["qe"] - b >= 100:
                        hit = dict(how="inside one alignment", left=seg_out(s), right=seg_out(s), reference_gap=0, read_gap=b - a, insert_bp=b - a, reference_position=tp,
                                   assembly_locus=at_embedded(s["t"], tp, tp)); break
                if hit: break
            if not hit:
                for (k1, L), (k2, R) in zip(slugs, slugs[1:]):
                    if k2 - k1 <= 3 and L["qe"] - 30 <= O["qs"] and O["qe"] <= R["qs"] + 30 and neighbours(L, R) and not any(t[m]["cls"] == "slug" for m in range(k1 + 1, k2)):
                        pa, pb = place(L, R)
                        hit = dict(how="between alignments (short piece)", left=seg_out(L), right=seg_out(R), reference_gap=ref_gap(L, R), read_gap=R["qs"] - L["qe"], insert_bp=ln,
                                   reference_position=pa, assembly_locus=at_embedded(L["t"], pa, pb)); break
            if hit and O["source"] == "first pass":
                ok_, pid_, share_ = confirm(r, O)
                if ok_ is False:          # the seeds matched, the bases do not: a repeat that resembles a repeat of the panel sequence
                    shape_n[("short piece not confirmed base by base", O["cls"])] += 1
                    unconfirmed.append(dict(read=r, piece=O, identity=pid_, share_aligned=share_, place=dict(sequence=hit["left"]["t"], position=hit["reference_position"])))
                    done_here.append((O["qs"], O["qe"], O["cls"])); continue
                if ok_: O = dict(O, ident=pid_, confirmed="base by base")
            if hit:
                between = hit["how"] != "inside one alignment"
                hit.update(read=r, read_length=t[0]["qlen"], shape="insertion", insert=O["cls"], middle=[O], n_segments=n, folds_back_after=False,
                           unaligned_left=O["qs"] - hit["left"]["qe"] if between else None, unaligned_right=hit["right"]["qs"] - O["qe"] if between else None)
                cand.append(hit); shape_n[("insertion", O["cls"])] += 1; shape_n[("insertion (short piece)", O["cls"])] += 1
                done_here.append((O["qs"], O["qe"], O["cls"]))
            else:
                for k, s in slugs:
                    cov = min(s["qe"], O["qe"]) - max(s["qs"], O["qs"]) - sum(max(0, min(b, O["qe"]) - max(a, O["qs"])) for a, b, _ in s.get("ins", []))
                    if cov >= 0.6 * ln:
                        shared[O["cls"]] += 1; done_here.append((O["qs"], O["qe"], O["cls"])); break
    # ---- joins between slug DNA and organelle DNA that are not part of an insertion shape
    for i in range(n - 1):
        a, b = t[i], t[i + 1]
        if i in used and i + 1 in used: continue
        pair = {a["cls"], b["cls"]}
        if "slug" in pair and len(pair) == 2 and (pair - {"slug"}) <= set(ORG) | {"algrDNA", "slugrDNA"}:
            o = (pair - {"slug"}).pop()
            if o not in ORG: continue
            sl = a if a["cls"] == "slug" else b
            # the same slug piece read backwards on the far side of the organelle piece: one join, read twice
            fold = (i + 2 < n and a["cls"] == "slug" and folded(a, t[i + 2])) or (i >= 1 and b["cls"] == "slug" and folded(b, t[i - 1]))
            shape = "fold-back" if fold else "one join"
            shape_n[(shape, o)] += 1
            cand.append(dict(read=r, read_length=a["qlen"], shape=shape, insert=o, left=seg_out(a), right=seg_out(b), gap=b["qs"] - a["qe"], n_segments=n, slug_mapq=sl.get("mapq", 0)))
    # ---- what the library does with slug DNA alone
    sl = [s for s in t if s["cls"] == "slug"]
    if sl and kinds == {"slug"}:
        nuclear["reads_slug_only"] += 1
        if len(sl) >= 2:
            pairs = list(zip(sl, sl[1:]))
            if any(not neighbours(a, b) and not folded(a, b) and elsewhere(b, a) for a, b in pairs): nuclear["reads_joining_unrelated_slug_pieces"] += 1
            if any(folded(a, b) for a, b in pairs): nuclear["reads_folding_back"] += 1
        # the comparable count: a read that lies at one clear place, joined to a piece that lies at another clear
        # place. Left out: joins at the end of an assembly sequence, where the read may simply run on into the
        # next contig, and pieces that cannot be placed (repeats).
        main = max(sl, key=lambda s: s["aln"])
        if main["mapq"] >= UNIQ:
            nuclear["reads_slug_placed"] += 1
            def at_end(a, b):
                # the two ends that meet at the join, on the reference
                ea = a["te"] if a["strand"] == "+" else a["ts"]; eb = b["ts"] if b["strand"] == "+" else b["te"]
                return min(ea, a["tlen"] - ea) < 2000 or min(eb, b["tlen"] - eb) < 2000
            for a, b in zip(sl, sl[1:]):
                if main not in (a, b): continue
                if a["mapq"] >= UNIQ and b["mapq"] >= UNIQ and not neighbours(a, b) and not folded(a, b) and elsewhere(b, a) and not at_end(a, b):
                    nuclear["reads_placed_joined_to_placed_elsewhere"] += 1; break
# a short piece of the same organelle stretch inserted at two or more different places: a mobile element
byw = collections.defaultdict(list)
for c in cand:
    if c["shape"] == "insertion" and c["how"] != "between alignments": byw[(c["middle"][0]["t"], c["middle"][0]["ts"] // 150)].append(c)
for grp in byw.values():
    if len({(c["left"]["t"], (c.get("reference_position") or 0) // 1000) for c in grp}) >= 2:
        for c in grp:
            shape_n[("insertion", c["insert"])] -= 1; shape_n[("insertion (short piece)", c["insert"])] -= 1; shape_n[("element", c["insert"])] += 1; c["shape"] = "element"
shape_n = collections.Counter({k: v for k, v in shape_n.items() if v})
reads_with_slug_join = collections.Counter()
seen = collections.defaultdict(set)
seen_placed = collections.defaultdict(set)
for c in cand:
    if c["shape"] == "element": continue
    seen[c["insert"]].add(c["read"])
    q = c.get("slug_mapq") if "slug_mapq" in c else min(c["left"].get("mapq", 0), c["right"].get("mapq", 0)) if c["shape"] == "insertion" else max(m.get("mapq", 0) for m in c.get("middle", [{}]))
    if (q or 0) >= UNIQ: seen_placed[c["insert"]].add(c["read"])
for k, v in seen.items(): reads_with_slug_join[k] = len(v)
reads_joined_to_placed_slug = {k: len(v) for k, v in seen_placed.items()}
n_reads = len(tiles)
# ---- the place of each insertion-shaped read on the slug sequence, and the other reads that cross it
def point(c):
    if c["how"] == "inside one alignment": return c["reference_position"]
    L, R = c["left"], c["right"]
    # the flank that the insert touches (fewer unexplained bases of the read between them) marks the place
    right_side = (c.get("unaligned_right") or 0) <= (c.get("unaligned_left") or 0)
    f = R if right_side else L
    return (f["ts"] if f["strand"] == "+" else f["te"]) if right_side else (f["te"] if f["strand"] == "+" else f["ts"])
for c in cand:
    if c["shape"] == "insertion": c["place"] = dict(sequence=c["left"]["t"], position=point(c))
if A.depth_paf:
    want = collections.defaultdict(list)
    for c in cand:
        if c["shape"] == "insertion": want[c["place"]["sequence"]].append(c); c["other_reads_crossing"] = 0; c["other_reads_crossing_names"] = []
    for line in open(A.depth_paf):
        f = line.split("\t", 12)
        if f[5] not in want: continue
        ts, te = int(f[7]), int(f[8])
        for c in want[f[5]]:
            p_ = c["place"]["position"]
            if ts <= p_ - 150 and te >= p_ + 150 and f[0] != c["read"] and int(f[11]) >= 30 and f[0] not in c["other_reads_crossing_names"]:
                c["other_reads_crossing"] += 1; c["other_reads_crossing_names"].append(f[0])
# ---- depth, and the depth of reads that could show an insertion of a given length with both flanks
depth = None
if A.genome_bp:
    counts = [open(f).read().split() for f in sorted(glob.glob(os.path.join(W, "*.count")))]
    lib_bases = sum(int(c[1]) for c in counts); lib_reads = sum(int(c[0]) for c in counts)
    depth = dict(genome_bp=A.genome_bp, library_bases=lib_bases, library_reads=lib_reads, depth=round(lib_bases / A.genome_bp, 2) if lib_bases else None)
    if A.paf:
        one_in = A.one_in
        if not one_in:
            sel = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(W, "*.select.json")))]
            one_in = max((x.get("sample_one_in", 0) for x in sel), default=0)
            thin = os.path.join(W, "sample.thin")            # 33_pass2.sh aligned every k-th read of the sample
            if os.path.exists(thin): one_in *= int(open(thin).read().split()[0])
        if one_in:
            lens = [ss[0]["qlen"] for ss in S.values() if ss]
            depth.update(sample_one_in=one_in, sample_reads_aligned=len(lens), sample_bases_aligned=sum(lens),
                         spanning_depth={str(L): round(one_in * sum(max(0, x - L - 2 * MINSEG) for x in lens) / A.genome_bp, 2) for L in (100, 500, 1000, 3000)})
            # copies: reads of the sample that lie at one clear place and share strand and both ends with another read
            # (within 10 bases; 30 for the older reads, whose ends are placed less exactly). With the sample being a
            # share f of the library, a read that has c copies in the library has one in the sample with probability
            # about c * f. Reads from one sequencing well of the older machines (numbers within 15 of each other) are
            # passes over one molecule and are not compared. As a measure of chance, the pairs that share one end and
            # miss the other by one to two times the tolerance are counted too.
            tol_ = 10 if kind == "hifi" else 30
            def well_(r):
                x = r.rsplit(".", 1)[-1]; return int(x) if x.isdigit() else None
            byts = collections.defaultdict(list)
            for r_, ss in S.items():
                m_ = max(ss, key=lambda x: x["qe"] - x["qs"]) if ss else None
                if m_ and m_["cls"] == "slug" and m_["mapq"] >= UNIQ and m_["qe"] - m_["qs"] >= 0.8 * m_["qlen"]:
                    lo_, hi_ = (m_["ts"] - m_["qs"], m_["te"] + m_["qlen"] - m_["qe"]) if m_["strand"] == "+" else (m_["ts"] - (m_["qlen"] - m_["qe"]), m_["te"] + m_["qs"])
                    byts[(m_["t"], m_["strand"])].append((lo_, hi_, r_))
            placed_ = sum(len(v) for v in byts.values()); withcopy = set(); pairs0 = pairs1 = 0
            for v in byts.values():
                v.sort()
                for i_ in range(len(v)):
                    j_ = i_ + 1
                    while j_ < len(v) and v[j_][0] - v[i_][0] <= tol_:
                        wi_, wj_ = well_(v[i_][2]), well_(v[j_][2]); d_ = abs(v[j_][1] - v[i_][1]); j_ += 1
                        if kind == "clr" and wi_ is not None and wj_ is not None and abs(wi_ - wj_) <= 15: continue
                        if d_ <= tol_: pairs0 += 1; withcopy.add(v[i_][2]); withcopy.add(v[j_ - 1][2])
                        elif d_ <= 2 * tol_: pairs1 += 1
            frac_ = len(lens) / lib_reads if lib_reads else 1 / one_in
            depth["copies"] = dict(placed_reads=placed_, reads_with_a_copy_in_the_sample=len(withcopy), pairs_sharing_both_ends=pairs0, pairs_sharing_one_end_and_just_missing_the_other=pairs1,
                                   tolerance=tol_, sample_fraction=round(frac_, 6),
                                   copies_per_read=round(len(withcopy) / placed_ / frac_, 4) if placed_ else None,
                                   copies_per_read_upper=round((len(withcopy) + 2 * (len(withcopy) ** 0.5) + 3) / placed_ / frac_, 4) if placed_ else None)
# insertion shapes with slug DNA from elsewhere: how many stand alone at their place
by_place = collections.Counter((x["left"]["t"], x["left"]["te"] // 2000) for x in nuc_rows)
for x in nuc_rows: x["reads_at_place"] = by_place[(x["left"]["t"], x["left"]["te"] // 2000)]
uniq = [x for x in nuc_rows if x["insert_mapq_min"] >= 30]
res = dict(kind=kind, min_identity=MINID, min_segment=MINSEG, min_short_piece=MINSHORT, reads_examined=n_reads, depth=depth, reads_with_piece=dict(reads_with),
           reads_joining_slug=dict(reads_with_slug_join), reads_joined_to_placed_slug=reads_joined_to_placed_slug,
           shapes={f"{a}|{b}": v for (a, b), v in sorted(shape_n.items())},
           shared_with_assembly=dict(shared), second_pass={f"{a}|{b}": v for (a, b), v in sorted(strands.items())},
           nuclear=dict(nuclear),
           insertion_slug_elsewhere=dict(n=len(nuc_rows), unique=len(uniq), unique_alone_at_place=sum(1 for x in uniq if x["reads_at_place"] == 1),
                                         unique_insert_bp=sorted(x["insert_bp"] for x in uniq), rows=sorted(nuc_rows, key=lambda x: -x["insert_mapq_min"])[:400]),
           candidates=[c for c in cand if c["shape"] in ("insertion", "mirror")], elements=[c for c in cand if c["shape"] == "element"][:200], short_pieces_not_confirmed=unconfirmed[:500],
           one_join=[c for c in cand if c["shape"] not in ("insertion", "mirror", "element")][:20000])
json.dump(res, open(A.out, "w"), indent=1)
print(f"{n_reads} reads; with a piece of:", dict(reads_with))
print("  shapes:", {f"{a}|{b}": v for (a, b), v in sorted(shape_n.items())}, "| organelle-like pieces shared with the assembly:", dict(shared))
print("  reads that come back over the place on the other strand:", {f"{a}|{b}": v for (a, b), v in sorted(strands.items())})
print("  joined to slug DNA at one clear place:", reads_joined_to_placed_slug)
print("  slug DNA alone:", dict(nuclear), "| insertion shape with slug DNA from elsewhere:", len(nuc_rows), "of which unique (MAPQ >= 30):", len(uniq),
      "of which alone at their place:", res["insertion_slug_elsewhere"]["unique_alone_at_place"])
if depth: print("  depth:", depth)
for c in res["candidates"][:40]:
    print(f'   {c["shape"]} of {c["insert"]} ({c["how"]}): {c["read"]} ({c["read_length"]} bp) {c["left"]["t"]}:{c["left"]["ts"]}-{c["left"]["te"]} | ' +
          ", ".join(f'{m["t"]}:{m["ts"]}-{m["te"]}' for m in c["middle"]) + f' | {c["right"]["ts"]}-{c["right"]["te"]}; flanks {c["reference_gap"]} bp apart, {c["read_gap"]} bp on the read' +
          (f' [{c["assembly_locus"]}]' if c.get("assembly_locus") else "") + (f' [{c["second_pass"]}]' if c.get("second_pass") else "") +
          (f' [{c["other_reads_crossing"]} other reads cross {c["place"]["sequence"]}:{c["place"]["position"]}]' if "other_reads_crossing" in c else ""))
