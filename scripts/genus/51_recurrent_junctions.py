#!/usr/bin/env python3
"""Long reads: is any join between organelle DNA and slug DNA carried by two molecules?

usage: 51_recurrent_junctions.py CANDIDATES.json READS.fa OUT.json [hifi|clr] [WORKDIR_TMP] [--repeat REPEAT.fa PANEL.fa.fai]
  --repeat REPEAT.fa PANEL.fa.fai
                       a repeat of the slug genome: the joins whose slug side is this repeat are
                       counted, with the number of pairs of molecules whose joins lie within 40
                       bases of each other on the organelle genome, and the number of such pairs
                       that chance gives (joins spread evenly over the organelle genome, whose
                       length is read from the index of the panel)

34_read_analysis.py looks for joins that share both break points on the reference sequences. That
misses a join into repeated slug DNA (a satellite, a transposon), because each read is then placed
on a different copy of the repeat. Here the reads themselves are compared. Around the join of every
read with one join (46_candidates.py; the two pieces must meet within 60 bases) a window of up to
400 bases either side is cut out, turned so that the organelle DNA lies on the left, and all
windows are aligned to each other (blastn). A pair of reads is a candidate when one alignment
covers at least 120 bases on both sides of the join in both reads. That alignment is lenient: it
bridges a stretch of some tens of bases that one read has and the other lacks. Each candidate pair
is therefore aligned again base by base, 200 bases either side of the join, and the two reads
share a join only if that alignment runs through the join for 100 bases on both sides, puts the
join of one read on the join of the other (within 6 bases for accurate reads, 15 for the older
ones), has no gap of that length or more within 80 bases of the join, and shows the two reads to
be the same sequence on each side of it (98% identical over those 100 bases for accurate reads,
70% for the older ones). Pairs that fail are joins that lie close to each other without being
the same, or the same organelle break point joined to two different members of a repeat family;
those from different molecules are listed as "near".
Reads from the same sequencing well (spot numbers within 15 of each other: the several passes of
one molecule) count as one molecule, and so do reads that share an end: the organelle piece or
the slug piece ends at the same point of its reference sequence, or the join lies at the same
distance from an end of the read (copies of one molecule made by PCR, in amplified libraries,
share its ends).
Output: groups of two or more molecules with the same join, and the near pairs.
"""
import sys, os, json, subprocess, collections
ARGV = list(sys.argv); REPEAT = None
if "--repeat" in ARGV:
    i_ = ARGV.index("--repeat"); REPEAT = ARGV[i_ + 1]; TLEN = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(ARGV[i_ + 2])}; del ARGV[i_:i_ + 3]
cand_json, reads_fa, out = ARGV[1:4]
kind = ARGV[4] if len(ARGV) > 4 else "clr"
tmp = ARGV[5] if len(ARGV) > 5 else out + ".tmp"
os.makedirs(tmp, exist_ok=True)
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin")
FL = 400; NEED = 120; ORG = ("cp", "slugmt", "algmt")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
d = json.load(open(cand_json))
J = {}
for c in d["one_join"]:
    a, b = c["left"], c["right"]
    if c["insert"] not in ORG or abs(c.get("gap", 0)) > 60: continue        # the two pieces must meet: no unexplained stretch between them
    org_left = a["cls"] != "slug"
    j = (a["qe"] + b["qs"]) // 2
    key = f'{c["read"]}|{j}'
    if key not in J: J[key] = dict(read=c["read"], j=j, org_left=org_left, insert=c["insert"], org=a if org_left else b, slug=b if org_left else a, wlen_read=c["read_length"])
want = collections.defaultdict(list)
for k, v in J.items(): want[v["read"]].append(k)
n = 0
with open(f"{tmp}/windows.fa", "w") as fo:
    name = None
    for line in open(reads_fa):
        if line[0] == ">": name = line[1:].split()[0]; continue
        if name in want:
            s = line.strip()
            for k in want[name]:
                v = J[k]; a = max(0, v["j"] - FL); b = min(len(s), v["j"] + FL); w = s[a:b]; jj = v["j"] - a
                if not v["org_left"]: w = w.translate(COMP)[::-1]; jj = len(w) - jj
                v["wj"] = jj; v["wlen"] = len(w); fo.write(f">{k}\n{w}\n"); n += 1
            name = None
print(n, "windows around joins")
opts = ["-task", "megablast", "-word_size", "20"] if kind == "hifi" else ["-task", "blastn", "-word_size", "11", "-reward", "1", "-penalty", "-1", "-gapopen", "2", "-gapextend", "1"]
subprocess.run([BIN + "/makeblastdb", "-in", f"{tmp}/windows.fa", "-dbtype", "nucl", "-out", f"{tmp}/w"], stdout=subprocess.DEVNULL, check=True)
with open(f"{tmp}/all.tsv", "w") as fo:
    subprocess.run([BIN + "/blastn"] + opts + ["-query", f"{tmp}/windows.fa", "-db", f"{tmp}/w", "-evalue", "1e-20", "-strand", "plus", "-num_threads", "1", "-max_target_seqs", "200",
                    "-outfmt", "6 qseqid sseqid qstart qend sstart send pident length"], stdout=fo, check=True)
def spot(r): return int(r.rsplit(".", 1)[1]) if r.rsplit(".", 1)[1].isdigit() else hash(r)
from Bio import Align
ALN = Align.PairwiseAligner(); ALN.mode = "local"; ALN.match_score = 1; ALN.mismatch_score = -1; ALN.open_gap_score = -2; ALN.extend_gap_score = -1
H = 200; NEEDV = 100; GAPMAX, TOL = (6, 6) if kind == "hifi" else (15, 15)
MINSAME = 0.98 if kind == "hifi" else 0.70          # two error-prone reads of one sequence agree at about 75 to 80% of positions
SEQ = {}; _n = None
for line in open(f"{tmp}/windows.fa"):
    if line[0] == ">": _n = line[1:].strip()
    else: SEQ[_n] = line.strip().upper()
def org_at(v):
    """position on the organelle sequence at which the organelle piece meets the slug DNA"""
    o = v["org"]; left = v["org_left"]
    return (o["te"] if o["strand"] == "+" else o["ts"]) if left else (o["ts"] if o["strand"] == "+" else o["te"])
_cmp = {}
def compare(q, s):
    """the joins of two reads, base by base: (same join?, what the alignment shows)"""
    key = (min(q, s), max(q, s))
    if key in _cmp: return _cmp[key]
    A, B = J[key[0]], J[key[1]]
    a0 = max(0, A["wj"] - H); a = SEQ[key[0]][a0:A["wj"] + H]; ja = A["wj"] - a0
    b0 = max(0, B["wj"] - H); b = SEQ[key[1]][b0:B["wj"] + H]; jb = B["wj"] - b0
    al = ALN.align(a, b)[0]; qb, sb = al.aligned
    if not len(qb): _cmp[key] = (False, dict(spans=False)); return _cmp[key]
    spans = bool(qb[0][0] <= ja - NEEDV and qb[-1][1] >= ja + NEEDV and sb[0][0] <= jb - NEEDV and sb[-1][1] >= jb + NEEDV)
    qb = [(int(x0), int(x1)) for x0, x1 in qb]; sb = [(int(y0), int(y1)) for y0, y1 in sb]
    gap = 0
    for i in range(1, len(qb)):
        if min(abs(qb[i][0] - ja), abs(qb[i - 1][1] - ja), abs(sb[i][0] - jb), abs(sb[i - 1][1] - jb)) <= 80:
            gap = max(gap, int(qb[i][0] - qb[i - 1][1]), int(sb[i][0] - sb[i - 1][1]))
    mapped = None
    for (x0, x1), (y0, y1) in zip(qb, sb):
        if x0 <= ja <= x1: mapped = int(y0 + (ja - x0)); break
        if x0 > ja: mapped = int(y0); break
    off = None if mapped is None else abs(mapped - jb)
    # identity of the two reads within 100 bases of the join, on each side: two copies of one sequence agree as far as
    # the reading errors allow; two different members of a repeat family do not
    def ident(lo, hi):
        m = n = 0
        for i, ((x0, x1), (y0, y1)) in enumerate(zip(qb, sb)):
            a_, b_ = max(x0, lo), min(x1, hi)
            if b_ > a_: m += sum(1 for k in range(a_, b_) if a[k] == b[y0 + (k - x0)]); n += b_ - a_
            if i + 1 < len(qb) and lo <= x1 < hi: n += int(max(qb[i + 1][0] - x1, sb[i + 1][0] - y1))
        return float(m / n) if n else 0.0
    idl, idr = ident(ja - NEEDV, ja), ident(ja, ja + NEEDV)
    same = bool(spans and gap < GAPMAX and off is not None and off <= TOL and min(idl, idr) >= MINSAME)
    _cmp[key] = (same, dict(runs_through_the_join=spans, longest_gap_near_join=gap, join_offset=off, identity_organelle_side=round(idl, 3), identity_slug_side=round(idr, 3))); return _cmp[key]
parent = {k: k for k in J}
def find(x):
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
links = {}; near = {}
for line in open(f"{tmp}/all.tsv"):
    q, s, qs, qe, ss, se, pid, ln = line.split("\t")
    if q == s or q not in J or s not in J: continue
    A, B = J[q], J[s]
    if A["insert"] != B["insert"] or "wj" not in A or "wj" not in B: continue
    if int(qs) <= A["wj"] - NEED and int(qe) >= A["wj"] + NEED and int(ss) <= B["wj"] - NEED and int(se) >= B["wj"] + NEED:
        # the join must sit at the same point of the alignment in both reads
        if abs((A["wj"] - int(qs)) - (B["wj"] - int(ss))) > (15 if kind == "hifi" else 40): continue
        same, how = compare(q, s)
        if same: links[(min(q, s), max(q, s))] = float(pid); parent[find(q)] = find(s)
        else: near[(min(q, s), max(q, s))] = how
def sides(k):
    """bases of the read on the organelle side and on the slug side of the join"""
    v = J[k]; return (v["j"], v["wlen_read"] - v["j"]) if v["org_left"] else (v["wlen_read"] - v["j"], v["j"])
def outer(v, side):
    """where on its reference sequence the organelle piece ("org") or the slug piece ("slug") ends, away from the join"""
    s = v[side]; left = (side == "org") == v["org_left"]
    return (s["ts"] if s["strand"] == "+" else s["te"]) if left else (s["te"] if s["strand"] == "+" else s["ts"])
def shares_end(k, m):
    """copies of one molecule made by PCR share an end: the same outer end on a reference sequence (within 10 bases),
    or the join at the same distance from an end of the read (within 30 bases, or 1% of that distance)"""
    A, B = J[k], J[m]
    for side in ("org", "slug"):
        if A[side]["t"] == B[side]["t"] and abs(outer(A, side) - outer(B, side)) <= 10: return True
    return any(abs(x - y) <= max(30, 0.01 * max(x, y)) for x, y in zip(sides(k), sides(m)))
groups = collections.defaultdict(list)
for k in J: groups[find(k)].append(k)
res = []
for g in groups.values():
    if len(g) < 2: continue
    mols = []
    for k in sorted(g, key=lambda k: spot(J[k]["read"])):
        # another molecule: not from the same well, and not a copy made by PCR. Copies of one molecule share its ends:
        # the join lies at the same distance from an end of the read, on the organelle side or on the slug side
        def other(m): return abs(spot(J[k]["read"]) - spot(J[m]["read"])) > 15 and not shares_end(k, m)
        if not mols or all(other(m) for m in mols): mols.append(k)
    v0 = J[g[0]]
    res.append(dict(insert=v0["insert"], reads=len({J[k]["read"] for k in g}), molecules=len(mols),
                    members=[dict(read=J[k]["read"], join_at=J[k]["j"], organelle=J[k]["org"]["t"], organelle_from=J[k]["org"]["ts"], organelle_to=J[k]["org"]["te"], organelle_strand=J[k]["org"]["strand"],
                                  slug=J[k]["slug"]["t"], slug_from=J[k]["slug"]["ts"], slug_to=J[k]["slug"]["te"], slug_mapq=J[k]["slug"]["mapq"]) for k in g]))
res.sort(key=lambda x: (-x["molecules"], -x["reads"]))
# joins of different molecules that lie close to each other without being the same join
def member(k): return dict(read=J[k]["read"], join_at=J[k]["j"], organelle=J[k]["org"]["t"], organelle_position=org_at(J[k]), organelle_from=J[k]["org"]["ts"], organelle_to=J[k]["org"]["te"],
                           organelle_strand=J[k]["org"]["strand"], slug=J[k]["slug"]["t"], slug_from=J[k]["slug"]["ts"], slug_to=J[k]["slug"]["te"], slug_mapq=J[k]["slug"]["mapq"])
near_rows = []; seen = set()
for (q, s), how in near.items():
    if find(q) == find(s): continue                                  # linked after all through a third read
    if abs(spot(J[q]["read"]) - spot(J[s]["read"])) <= 15: continue    # passes of one molecule that differ by a reading error
    if shares_end(q, s): continue                                      # copies of one molecule
    kk = tuple(sorted((find(q), find(s))))
    if kk in seen: continue
    seen.add(kk)
    near_rows.append(dict(insert=J[q]["insert"], organelle_positions_apart=abs(org_at(J[q]) - org_at(J[s])), alignment=how, members=[member(q), member(s)]))
near_rows.sort(key=lambda x: x["organelle_positions_apart"])
# ---- joins into one repeat of the slug genome: how close do they come to each other on the organelle genome?
repeat = None
if REPEAT:
    hit = collections.Counter()
    r = subprocess.run([BIN + "/blastn"] + opts + ["-query", f"{tmp}/windows.fa", "-subject", REPEAT, "-evalue", "1e-10", "-outfmt", "6 qseqid qstart qend"], capture_output=True, text=True, check=True)
    cov = collections.defaultdict(set)
    for line in r.stdout.splitlines():
        q, a, b = line.split("\t"); a, b = int(a), int(b)
        if q in J and "wj" in J[q]: cov[q].update(range(max(a, J[q]["wj"]), b))          # bases of the slug side that match the repeat
    mols = []
    for k in sorted((k for k in cov if len(cov[k]) >= 100), key=lambda k: spot(J[k]["read"])):
        v = J[k]; o = v["org"]
        below = (o["strand"] == "+") == v["org_left"]                # the organelle piece lies at positions below the join
        if mols and abs(spot(v["read"]) - spot(J[mols[-1]]["read"])) <= 15: continue      # another pass of the same molecule
        v["_below"] = below; mols.append(k)
    by = collections.defaultdict(list)
    for k in mols: by[(J[k]["insert"], J[k]["_below"])].append((org_at(J[k]), J[k]["read"]))
    close = []
    for (ins, below), v in by.items():
        v.sort()
        for i in range(len(v)):
            j = i + 1
            while j < len(v) and v[j][0] - v[i][0] <= 40: close.append(dict(insert=ins, positions=[v[i][0], v[j][0]], reads=[v[i][1], v[j][1]])); j += 1
    nn = collections.Counter(J[k]["insert"] for k in mols); olen = {J[k]["insert"]: TLEN.get(J[k]["org"]["t"]) for k in mols}
    expected = {ins: round(nn[ins] * (nn[ins] - 1) / 2 * (80 / olen[ins]) * 0.5, 2) for ins in nn if olen.get(ins)}
    repeat = dict(file=os.path.basename(REPEAT), molecules=dict(nn), pairs_within_40_bases=close, pairs_expected_by_chance=expected)
    print("joins into the repeat, by organelle:", dict(nn), "| pairs within 40 bases:", len(close), "| expected by chance:", expected)
json.dump(dict(joins_examined=n, groups_same_well_only=sum(1 for x in res if x["molecules"] < 2), recurrent=[x for x in res if x["molecules"] >= 2], near=near_rows, joins_into_repeat=repeat,
               same_join_means=dict(runs_through_join_bases=NEEDV, join_offset_max=TOL, gap_max=GAPMAX, identity_min=MINSAME)), open(out, "w"), indent=1)
print("groups of reads sharing a join:", len(res), "| of which from two or more molecules:", sum(1 for x in res if x["molecules"] >= 2), "| pairs of molecules with joins close to each other but not the same:", len(near_rows))
for x in near_rows[:10]:
    print(f'  near, {x["insert"]}: organelle positions {x["organelle_positions_apart"]} apart; {x["alignment"]}')
    for m in x["members"]: print(f'      {m["read"]:24s} {m["organelle"]}:{m["organelle_position"]} | {m["slug"]}:{m["slug_from"]}-{m["slug_to"]} mapq {m["slug_mapq"]}')
for x in res:
    if x["molecules"] >= 2:
        print(f'  {x["insert"]}: {x["reads"]} reads, {x["molecules"]} molecules')
        for m in x["members"]: print(f'      {m["read"]:24s} {m["organelle"]}:{m["organelle_from"]}-{m["organelle_to"]} {m["organelle_strand"]} | {m["slug"]}:{m["slug_from"]}-{m["slug_to"]} mapq {m["slug_mapq"]}')
