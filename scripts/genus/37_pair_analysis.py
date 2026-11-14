#!/usr/bin/env python3
"""Short reads: what the library contains, and every read pair that joins slug DNA to organelle DNA.

usage: 37_pair_analysis.py WORKDIR panel.tsv OUT.json [SCREEN.json] [--genome-bp N] [--genome-acgt N] [--genome GENOME.fa.gz] [--junctions OUT.tsv.gz]
                           [--organelle-stretch SEQUENCE:START-END=cp ...]
  --genome-acgt N   bases of the assembly that are not N (default: all of --genome-bp)
  --organelle-stretch   a stretch of the assembly that is organelle DNA although the DNA screen of the
                    assembly (35_assembly_screen.py) could not recognise it, because no reference holds
                    that part of the organelle genome: mates on it count as organelle DNA, not slug DNA

WORKDIR holds pass2.paf (each mate aligned on its own, minimap2 -x sr -c), the RUN.count,
RUN.select.json and RUN.p1.paf files of the first pass (31c_screen_run.py) and, if the first pass
kept a random sample of the library, sample.paf (33_pass2.sh).

A mate is  cp / slugmt / algmt / algrDNA  when >= 80 of its bases align to such a sequence at
>= 97% identity, and  slug  when >= 80 bases align to the assembly at >= 97% with mapping quality
>= 20. Two kinds of molecule join organelle DNA to slug DNA:
  mixed pair   one mate organelle, the other slug
  split mate   one mate carries >= 40 bases of each itself, so the join is read base by base
A chloroplast fragment that entered a chromosome in one cell gives at most one such molecule per
junction, so every one of them is a candidate and is written to the junction table. A join made in
the tube while the library was prepared looks exactly the same. What can be measured is whether
there are more of them than the library makes by itself:
  - the same count for mitochondrial DNA, which is also loose in the tube
  - in the random sample: pairs whose two mates lie on unrelated slug sequences
Under joining at random the three rates are equal. Joins found twice, in molecules that are not
copies of each other, are listed: a clone of cells, or an inherited insertion. Three lists are made:
places of the slug genome with three or more different mixed pairs within 1 kb (recurrent), breaks
read through at the same base by two molecules (recurrent_breaks), and any two molecules, of either
kind, whose slug sides and organelle sides both lie within 1 kb of each other (same_join).

Pairs that are organelle DNA at both ends ("plain") may have been thinned in the first pass
(--plain-keep K), and in a library rich in organelle DNA also the pairs that match one kind of
organelle genome and nothing else without being plain (--match-keep K). Both are identified again
from RUN.p1.paf by the rules of the first pass, and every count is weighted: a kept pair stands
for K pairs. The counts given for such a library are therefore estimates from one pair in K.
"""
import sys, os, json, glob, gzip, collections, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sam_segments import paf_segments, best_tiling, relabel_organelle_loci
ap = argparse.ArgumentParser()
ap.add_argument("workdir"); ap.add_argument("panel_tsv"); ap.add_argument("out"); ap.add_argument("screen", nargs="?")
ap.add_argument("--genome-bp", type=int, default=0); ap.add_argument("--genome-acgt", type=int, default=0); ap.add_argument("--junctions"); ap.add_argument("--genome")
ap.add_argument("--organelle-stretch", action="append", default=[])
A = ap.parse_args()
W, panel_tsv, out, SCREEN = A.workdir, A.panel_tsv, A.out, A.screen
EXTRA = []          # stretches of the assembly that are organelle DNA, given by hand: SEQUENCE:START-END=cp|slugmt
for x in A.organelle_stretch:
    loc, kind = x.rsplit("=", 1); c_, r_ = loc.rsplit(":", 1); a_, b_ = (int(v) for v in r_.split("-")); EXTRA.append((c_, a_, b_, kind))
FAR = 20000                     # two slug mates further apart than this, or on two sequences, are "apart"
desc = {}
for l in list(open(panel_tsv))[1:]:
    f = l.rstrip("\n").split("\t"); desc[f[0]] = f[4]
def taxon(t):
    if t.startswith("asmcp|") or "|" not in t: return "plastid sequence in the assembly"          # also a stretch of an assembly sequence counted as plastid DNA
    if t.startswith("kcp|"): return "plastid contigs assembled from these reads"
    d = desc.get(t, "")
    w = d.replace("UNVERIFIED: ", "").replace("UNVERIFIED_ORG: ", "").replace("UNVERIFIED_ASMBLY: ", "").split()
    return " ".join(w[:2]) if len(w) >= 2 else (d or t)
def norm(c): return {"asmcp": "cp", "kcp": "cp", "asmalgmt": "algmt"}.get(c, c)
# the panel sequences, to see what a mate matched: a run of T in a panel sequence matches any run of T
PANELSEQ = {}; _pf = os.path.join(os.path.dirname(os.path.abspath(panel_tsv)), "panel.fa")
def panel_stretch(t, a, b):
    if not PANELSEQ and os.path.exists(_pf):
        name = None; buf = []
        for line in open(_pf):
            if line[0] == ">":
                if name: PANELSEQ[name] = "".join(buf)
                name = line[1:].split()[0]; buf = []
            else: buf.append(line.strip().upper())
        if name: PANELSEQ[name] = "".join(buf)
    return PANELSEQ.get(t, "")[max(0, a):b]
def simple(t, a, b):
    """is the stretch of the panel sequence that a mate matched simple sequence (few different triplets)?"""
    q = panel_stretch(t, a, b)
    if len(q) < 30: return False
    return len({q[i:i + 3] for i in range(len(q) - 2)}) < 0.25 * (len(q) - 2)
ORGP1 = {"cp": "cp", "asmcp": "cp", "algmt": "algmt", "asmalgmt": "algmt", "slugmt": "slugmt"}

# ---------------------------------------------------------------- the first pass: counts and weights
# plain_spots, partial_spots: kept pair -> the kind of organelle genome it matched in the first pass, by which it was thinned
runs = []; weight = {}; plain_spots = {}; plain_tot = collections.Counter(); plain_kept = collections.Counter()
pweight = {}; partial_spots = {}; partial_tot = collections.Counter(); partial_kept = collections.Counter()
for c in sorted(glob.glob(os.path.join(W, "*.count"))):
    n, b = (int(x) for x in open(c).read().split()); run = os.path.basename(c)[:-6]
    sel = json.load(open(os.path.join(W, run + ".select.json"))) if os.path.exists(os.path.join(W, run + ".select.json")) else {}
    runs.append(dict(run=run, reads=n, bases=b, hit_reads=sel.get("hit_reads"), plain_keep_one_in=sel.get("plain_keep_one_in", 1), match_keep_one_in=sel.get("match_keep_one_in", 1), sample_one_in=sel.get("sample_one_in", 0),
                     plain_pairs={k[6:]: v for k, v in sel.items() if k.startswith("plain_") and not k.startswith("plain_kept_") and k != "plain_keep_one_in"},
                     other_pairs_with_a_match=sel.get("other_pairs_with_a_match"), hot_windows=sel.get("hot_windows", [])[:10]))
    for k, v in sel.items():
        if k.startswith("plain_kept_"): plain_kept[k[11:]] += v
        elif k.startswith("plain_") and k != "plain_keep_one_in": plain_tot[k[6:]] += v
        elif k.startswith("partial_kept_"): partial_kept[k[13:]] += v
        elif k.startswith("partial_"): partial_tot[k[8:]] += v
    if sel.get("method") == "pieces" and os.path.exists(os.path.join(W, run + ".p1.paf")):
        # which of the kept pairs are plain: the rule of 31c_screen_run.py
        best = {}
        for line in open(os.path.join(W, run + ".p1.paf")):
            c2 = line.split("\t", 10); r = c2[0]; al = int(c2[9])
            if r not in best or al > best[r][0]: best[r] = (al, ORGP1.get(c2[5].split("|")[0], "other"), (int(c2[3]) - int(c2[2])) >= 0.8 * int(c2[1]))
        thinned = sel.get("match_keep_one_in", 1) > 1
        for r, v in best.items():
            sp_ = r[:-2]; other = best.get(sp_ + (".2" if r.endswith(".1") else ".1"))
            if other and r.endswith(".2"): continue          # the pair is looked at once, from its first mate (or from its only matched mate)
            if other and v[2] and other[2] and v[1] == other[1] and v[1] != "other": plain_spots[sp_] = v[1]
            elif thinned and v[1] != "other" and (other is None or other[1] == v[1]): partial_spots[sp_] = v[1]
        del best
for k in plain_tot: weight[k] = plain_tot[k] / plain_kept[k] if plain_kept.get(k) else 1.0
for k in partial_tot: pweight[k] = partial_tot[k] / partial_kept[k] if partial_kept.get(k) else 1.0
def wt(sp, cls=None):
    """how many pairs of the library a kept pair stands for: as many as were counted for each kept one among
    the pairs of its kind (the kind is that of the first pass, whatever the second pass calls the mates)"""
    if sp in plain_spots: return weight.get(plain_spots[sp], 1.0)
    if sp in partial_spots: return pweight.get(partial_spots[sp], 1.0)
    return 1.0
wsp = wt
tot_reads = sum(r["reads"] for r in runs); tot_bases = sum(r["bases"] for r in runs)
read_len = tot_bases / max(1, tot_reads)

def classify(S):
    """confident kind of each mate, and the mates that carry two kinds"""
    mate = {}; split = []
    for r, ss in S.items():
        for s in ss: s["cls"] = norm(s["cls"])
        conf = best_tiling([s for s in ss if s["aln"] >= 40 and s["ident"] >= 0.97 and (s["cls"] != "slug" or s["mapq"] >= 20)], max_overlap=10)
        big = [s for s in conf if s["qe"] - s["qs"] >= 80]
        if big: mate[r] = max(big, key=lambda s: s["aln"])
        kinds = {s["cls"] for s in conf}
        if "slug" in kinds and len(kinds) > 1:
            a = max((s for s in conf if s["cls"] == "slug"), key=lambda s: s["aln"]); o = max((s for s in conf if s["cls"] != "slug"), key=lambda s: s["aln"])
            first, second = (a, o) if a["qs"] <= o["qs"] else (o, a)
            def right_end(s): return s["te"] if s["strand"] == "+" else s["ts"]
            def left_end(s): return s["ts"] if s["strand"] == "+" else s["te"]
            split.append(dict(read=r, other=o["cls"], slug=a["t"], slug_break=right_end(a) if first is a else left_end(a), slug_strand=a["strand"],
                              other_t=o["t"], other_break=left_end(o) if first is a else right_end(o), other_strand=o["strand"], other_ts=o["ts"], other_te=o["te"],
                              slug_len=a["qe"] - a["qs"], other_len=o["qe"] - o["qs"], bases_between=second["qs"] - first["qe"], slug_mapq=a["mapq"],
                              read_start=min(a["ts"], a["te"]), slug_far=a["ts"] if (right_end(a) if first is a else left_end(a)) == a["te"] else a["te"],
                              end_piece="slug" if first is a else "organelle"))          # which kind of DNA the mate begins with: the piece at the end of the fragment
        elif len([s for s in conf if s["cls"] == "slug"]) >= 2:
            sl = sorted((s for s in conf if s["cls"] == "slug"), key=lambda s: s["qs"])
            a, b = sl[0], sl[1]
            if a["t"] != b["t"] or abs(a["ts"] - b["ts"]) > FAR: split.append(dict(read=r, other="slug_elsewhere", slug=a["t"], slug_break=a["te"], other_t=b["t"], other_break=b["ts"]))
    return mate, split

# ---------------------------------------------------------------- the chosen reads
S = paf_segments(os.path.join(W, "pass2.paf"), minlen=40)
relabel_organelle_loci(S, SCREEN, EXTRA)
don = collections.defaultdict(lambda: dict(mates=0, bp=0.0, ids=[], acc=collections.Counter()))
for r, ss in S.items():
    loose = [s for s in ss if norm(s["cls"]) == "cp" and s["aln"] >= 80 and s["ident"] >= 0.80]
    if loose:
        b = max(loose, key=lambda s: s["aln"] * s["ident"]); d = don[taxon(b["t"])]
        w = wt(r.rsplit(".", 1)[0], "cp")
        d["mates"] += w; d["bp"] += w * (b["qe"] - b["qs"]); d["ids"].append(b["ident"]); d["acc"][b["t"]] += 1
mate, split = classify(S)
n_examined = len(S); del S
cls_bp = collections.Counter(); by_taxon = {"algrDNA": collections.defaultdict(list), "algmt": collections.defaultdict(list)}
spots = collections.defaultdict(dict)
for r, b in mate.items():
    sp, i = r.rsplit(".", 1); spots[sp][i] = b
    w = wt(sp, b["cls"])
    if b["cls"] != "slug": cls_bp[b["cls"]] += w * (b["qe"] - b["qs"])
    if b["cls"] in by_taxon: by_taxon[b["cls"]][taxon(b["t"])].append(b["ident"])
cp_loose = sum(d["bp"] for d in don.values())
donors = []
for tx, d in sorted(don.items(), key=lambda kv: -kv[1]["bp"])[:15]:
    ids = sorted(d["ids"])
    donors.append(dict(name=tx, reads=round(d["mates"]), bp=round(d["bp"]), fraction=round(d["bp"] / max(1, cp_loose), 5), median_identity=round(ids[len(ids) // 2], 4),
                       identity_q10=round(ids[len(ids) // 10], 4), identity_q90=round(ids[(9 * len(ids)) // 10], 4), accession=d["acc"].most_common(1)[0][0].split("|", 1)[-1]))
have = collections.Counter(); have_n = collections.Counter(); mixed = collections.defaultdict(list); plain_but_mixed = collections.Counter(); simple_dropped = collections.Counter()
for sp, m in spots.items():
    kinds = {b["cls"] for b in m.values()}; pl = sp in plain_spots
    for k in kinds:
        if k != "slug": have[k] += wt(sp, k); have_n[k] += 1
    if len(m) == 2 and "slug" in kinds and len(kinds) == 2:
        a = next(b for b in m.values() if b["cls"] == "slug"); o = next(b for b in m.values() if b["cls"] != "slug")
        if pl: plain_but_mixed[o["cls"]] += 1; continue
        if simple(o["t"], o["ts"], o["te"]): simple_dropped[o["cls"]] += 1; continue
        mixed[o["cls"]].append(dict(spot=sp, slug=a["t"], slug_pos=(a["ts"] + a["te"]) // 2, slug_strand=a["strand"], slug_mapq=a["mapq"], other_t=o["t"], other_pos=(o["ts"] + o["te"]) // 2, other_strand=o["strand"], w=wsp(sp), slug_far=a["ts"] if a["strand"] == "+" else a["te"]))
split_org = [x for x in split if x["other"] != "slug_elsewhere"]
_keep = []
for x in split_org:
    if simple(x["other_t"], x["other_ts"], x["other_te"]): simple_dropped[x["other"] + " (split mate)"] += 1
    else: _keep.append(x)
split_org = _keep
for x in split_org: x["w"] = wsp(x["read"].rsplit(".", 1)[0])
split_n = collections.Counter(x["other"] for x in split_org)
split_w = collections.Counter()
for x in split_org: split_w[x["other"]] += x["w"]
# Copies of one molecule (the same fragment read twice by the sequencer, or copied by PCR, wherever in the run they
# are) begin at the same base at both ends. They are recognised by the reads themselves: two pairs are copies
# when both mates agree base for base (mate 1 with mate 1 and 2 with 2, or the other way round) in one of the
# first three stretches of WIN bases, allowing MISM differences for errors of the sequencer and a shift of
# one or two bases. Three stretches, because a read whose first cycles failed still agrees further on, while a
# fragment that begins at another place agrees nowhere. Where the reads are not at hand, by position: evidence of the same kind whose slug side
# ends at the same base, give or take COPY bases.
COPY = 10; WIN = 40; PREFIX = 3 * WIN; MISM = 3
NEAR = 1000
_fai = os.path.join(os.path.dirname(os.path.abspath(panel_tsv)), "panel.fa.fai")
TLEN = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(_fai)} if os.path.exists(_fai) else {}
ev = collections.defaultdict(list)
for k, rows in mixed.items():
    for x in rows: ev[k].append(dict(molecule=x["spot"], evidence="mixed pair", slug=x["slug"], slug_pos=x["slug_pos"], slug_far=x["slug_far"], slug_mapq=x["slug_mapq"], other_t=x["other_t"], other_pos=x["other_pos"]))
for x in split_org:
    ev[x["other"]].append(dict(molecule=x["read"].rsplit(".", 1)[0], evidence="split mate", slug=x["slug"], slug_pos=x["slug_break"], slug_far=x["slug_far"], slug_mapq=x["slug_mapq"], other_t=x["other_t"], other_pos=x["other_break"]))
_want = {x["molecule"] for rows in ev.values() for x in rows}; HEAD = collections.defaultdict(dict)
_hits = os.path.join(W, "all_hits.fa")
if _want and os.path.exists(_hits):
    with open(_hits) as fh:
        name = None
        for line in fh:
            if line[0] == ">":
                r = line[1:].split()[0]; sp_, _, i_ = r.rpartition("."); name = (sp_, i_) if sp_ in _want else None
            elif name and len(HEAD[name[0]].get(name[1], "")) < PREFIX:
                HEAD[name[0]][name[1]] = (HEAD[name[0]].get(name[1], "") + line.strip().upper())[:PREFIX]
def _close(a, b):
    a = a or ""; b = b or ""
    for i_ in range(0, len(a) - WIN + 1, WIN):
        for d_ in (0, -1, 1, -2, 2):          # a copy made by PCR may have lost or gained a base in a run of one letter
            if i_ + d_ < 0 or i_ + d_ + WIN > len(b): continue
            if sum(1 for x_, y_ in zip(a[i_:i_ + WIN], b[i_ + d_:i_ + d_ + WIN]) if x_ != y_) <= MISM: return True
    return False
def copies(m1, m2):
    """are two pairs copies of one fragment? None when the reads are not at hand"""
    h1, h2 = HEAD.get(m1), HEAD.get(m2)
    if not h1 or not h2 or len(h1) < 2 or len(h2) < 2: return None
    return (_close(h1["1"], h2["1"]) and _close(h1["2"], h2["2"])) or (_close(h1["1"], h2["2"]) and _close(h1["2"], h2["1"]))
_par = {}
def _find(a):
    while _par.setdefault(a, a) != a:
        _par[a] = _par[_par[a]]; a = _par[a]
    return a
def runs_of(rows, key, near):
    """rows sorted by key, cut wherever two neighbours are more than near apart"""
    rows = sorted(rows, key=key); out_ = []; cur = []
    for x in rows:
        if cur and key(x) - key(cur[-1]) > near: out_.append(cur); cur = []
        cur.append(x)
    if cur: out_.append(cur)
    return out_
# copies can only be pairs whose slug sides lie close together: compare those
for k, rows in ev.items():
    by_place = collections.defaultdict(list)
    for x in rows: by_place[x["slug"]].append(x)
    for xs in by_place.values():
        for g in runs_of(xs, lambda x: x["slug_pos"], NEAR):
            ms = sorted({x["molecule"] for x in g})
            if len(ms) < 2 or len(ms) > 400: continue
            far = collections.defaultdict(list)
            for x in g: far[x["molecule"]].append((x["evidence"], x["slug_far"]))
            for i_, m1 in enumerate(ms):
                for m2 in ms[i_ + 1:]:
                    if _find(m1) == _find(m2): continue
                    c_ = copies(m1, m2)
                    if c_ is None: c_ = any(e1 == e2 and abs(f1 - f2) <= COPY for e1, f1 in far[m1] for e2, f2 in far[m2])
                    if c_: _par[_find(m2)] = _find(m1)
def molecules(rows): return len({_find(x["spot"]) for x in rows})
def split_molecules(k): return len({_find(x["read"].rsplit(".", 1)[0]) for x in split_org if x["other"] == k})
def mixed_w(k): return sum(x["w"] for x in mixed[k])
# counts are weighted (a kept pair of a thinned kind stands for K pairs); the pairs actually examined are given as well
rates = {k: dict(pairs_with_mate=round(have[k]), pairs_also_with_slug=round(mixed_w(k)), molecules=molecules(mixed[k]), rate=round(mixed_w(k) / have[k], 6) if have[k] else None, split_mates=round(split_w.get(k, 0)),
                 split_molecules=split_molecules(k), split_rate=round(split_w.get(k, 0) / have[k], 6) if have[k] else None,
                 examined=dict(pairs_with_mate=have_n[k], pairs_also_with_slug=len(mixed[k]), split_mates=split_n.get(k, 0))) for k in have}
thinned_one_in = max([r_["match_keep_one_in"] for r_ in runs] or [1])
# Which piece sits at the end of the fragment in the mates that read through a join. A mate begins at the end of
# its fragment, so the piece it begins with is the short one (40 bases to the length of the mate less 40), and the
# other piece runs on to the far end. Where two fragments were joined end to end in the tube, the short piece is a
# short fragment; a kind of DNA that comes in few short fragments is then seldom the end piece, and is joined less
# often altogether. A join inside a chromosome is cut at random and has either kind at the end equally often.
end_piece = {}
for k_ in split_n:
    seen_ = set(); c_ = collections.Counter()
    for x in split_org:
        root_ = _find(x["read"].rsplit(".", 1)[0])
        if x["other"] != k_ or root_ in seen_: continue
        seen_.add(root_); c_[x["end_piece"]] += 1
    end_piece[k_] = dict(organelle=c_["organelle"], slug=c_["slug"])

# ---------------------------------------------------------------- joins found more than once
recurrent = []
for k, rows in mixed.items():
    bins = collections.defaultdict(list)
    for x in rows: bins[(x["slug"], x["slug_pos"] // 1000)].append(x)
    for (c, b), xs in bins.items():
        distinct = {(x["slug_pos"] // 20, x["other_pos"] // 20) for x in xs}
        if len(distinct) >= 3:
            recurrent.append(dict(kind=k, slug=c, start=b * 1000, pairs=len(xs), distinct_pairs=len(distinct), other_targets=sorted({x["other_t"] for x in xs})[:4],
                                  other_positions=sorted(x["other_pos"] for x in xs)[:20]))
same_break = collections.defaultdict(list)
for x in split_org: same_break[(x["other"], x["slug"], x["slug_break"] // 6, x["other_t"].split("|")[0], x["other_break"] // 6)].append(x)
recurrent_breaks = []
for key, xs in same_break.items():
    starts = {_find(x["read"].rsplit(".", 1)[0]) for x in xs}
    if len(starts) >= 2:
        recurrent_breaks.append(dict(kind=key[0], slug=key[1], slug_break=xs[0]["slug_break"], other_t=xs[0]["other_t"], other_break=xs[0]["other_break"], reads=len(xs), molecules=len(starts)))

# The same join in two molecules, whatever shows it: the two mates of a pair lying either side of it, or one mate
# read through it. Two molecules count when their slug sides lie within NEAR bases of each other on one sequence
# and their organelle sides within NEAR bases on one organelle sequence, and when they are not copies of one
# fragment (see above). Chance: n(n-1)/2 pairs of joins, each pair close on both sides with probability
# (2 NEAR / genome) x (2 NEAR / organelle sequence), if joins fell evenly.
same_join = []; same_join_chance = {}
for k, rows in ev.items():
    n = len({_find(x["molecule"]) for x in rows}); lens = [TLEN.get(x["other_t"], 0) for x in rows]
    L = sorted(lens)[len(lens) // 2] if lens else 0
    if A.genome_bp and L: same_join_chance[k] = round(n * (n - 1) / 2 * min(1.0, 2 * NEAR / (A.genome_acgt or A.genome_bp)) * min(1.0, 2 * NEAR / L), 4)
    by_slug = collections.defaultdict(list)
    for x in rows: by_slug[(x["slug"], x["other_t"])].append(x)
    for (c_, t_), xs in by_slug.items():
        if len(xs) < 2: continue
        for g in runs_of(xs, lambda x: x["slug_pos"], NEAR):
            for h in runs_of(g, lambda x: x["other_pos"], NEAR):
                mol = {_find(x["molecule"]) for x in h}
                if len(mol) < 2: continue
                same_join.append(dict(kind=k, slug=c_, slug_from=min(x["slug_pos"] for x in h), slug_to=max(x["slug_pos"] for x in h), other_t=t_, other_from=min(x["other_pos"] for x in h),
                                      other_to=max(x["other_pos"] for x in h), molecules=len(mol), copies=len({x["molecule"] for x in h}) - len(mol), read_through=len({_find(x["molecule"]) for x in h if x["evidence"] == "split mate"}),
                                      lowest_mapq=min(x["slug_mapq"] for x in h), evidence=[dict(molecule=x["molecule"], evidence=x["evidence"], slug_pos=x["slug_pos"], slug_far=x["slug_far"], other_pos=x["other_pos"]) for x in h][:30]))
same_join.sort(key=lambda x: -x["molecules"])

# a break that sits at a gap of the assembly (a run of N) is where the reads go on into slug DNA that
# the assembly lacks: a recurrent "join" there is no join at all
if A.genome and (recurrent_breaks or same_join):
    import subprocess
    ST = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/samtools"
    for x in recurrent_breaks:
        a = max(1, x["slug_break"] - 20); fa = subprocess.run([ST, "faidx", A.genome, f'{x["slug"]}:{a}-{x["slug_break"] + 20}'], capture_output=True, text=True).stdout
        x["at_assembly_gap"] = "N" in "".join(fa.split("\n")[1:]).upper()
    for x in same_join:
        a = max(1, x["slug_from"] - 200); fa = subprocess.run([ST, "faidx", A.genome, f'{x["slug"]}:{a}-{x["slug_to"] + 200}'], capture_output=True, text=True).stdout
        x["at_assembly_gap"] = "N" in "".join(fa.split("\n")[1:]).upper()

# ---------------------------------------------------------------- the random sample of the library
sample = None
if os.path.exists(os.path.join(W, "sample.paf")):
    SS = paf_segments(os.path.join(W, "sample.paf"), minlen=40)
    relabel_organelle_loci(SS, SCREEN, EXTRA)
    tlen = {}
    for ss in SS.values():
        for s in ss: tlen[s["t"]] = None
    for line in open(os.path.join(W, "sample.paf")):
        c = line.split("\t", 8)
        if tlen.get(c[5], 0) is None: tlen[c[5]] = int(c[6])
    smate, ssplit = classify(SS); del SS
    n_sample_reads = sum(1 for l in open(os.path.join(W, "all_sample.fa")) if l[0] == ">") if os.path.exists(os.path.join(W, "all_sample.fa")) else None
    sspots = collections.defaultdict(dict); sbp = collections.Counter()
    for r, b in smate.items():
        sp, i = r.rsplit(".", 1); sspots[sp][i] = b; sbp[b["cls"]] += b["qe"] - b["qs"]
    c = collections.Counter(); inserts = []
    for sp, m in sspots.items():
        kinds = [b["cls"] for b in m.values()]
        for k in set(kinds): c["pairs_with_" + k] += 1
        if kinds.count("slug") == 2:
            a, b = m["1"], m["2"]
            if a["t"] == b["t"] and abs(a["ts"] - b["ts"]) <= FAR:
                c["slug_pairs_together"] += 1
                if a["strand"] != b["strand"]: inserts.append(max(a["te"], b["te"]) - min(a["ts"], b["ts"]))
            else:
                c["slug_pairs_apart"] += 1
                if any(min(x["ts"], tlen[x["t"]] - x["te"]) < 2000 for x in (a, b)): c["slug_pairs_apart_near_an_end"] += 1
                if a["mapq"] >= 50 and b["mapq"] >= 50: c["slug_pairs_apart_mapq50"] += 1
        elif "slug" in kinds and len(m) == 2 and len(set(kinds)) == 2:
            c["pairs_slug_and_" + next(k for k in kinds if k != "slug")] += 1
    inserts.sort()
    ns = c["pairs_with_slug"]
    sample = dict(reads=n_sample_reads, mates_with_an_alignment=len(smate), counts=dict(c), insert_median=inserts[len(inserts) // 2] if inserts else None,
                  fraction_of_bases={k: round(v / (n_sample_reads * read_len), 6) for k, v in sbp.items()} if n_sample_reads else None,
                  slug_rate=round(c["slug_pairs_apart"] / ns, 6) if ns else None,
                  slug_rate_away_from_ends=round((c["slug_pairs_apart"] - c["slug_pairs_apart_near_an_end"]) / ns, 6) if ns else None,
                  split_mates_joining_unrelated_slug=sum(1 for x in ssplit if x["other"] == "slug_elsewhere"),
                  split_rate=round(sum(1 for x in ssplit if x["other"] == "slug_elsewhere") / ns, 6) if ns else None)

# ---------------------------------------------------------------- how much a real insertion would add
sens = None
if A.genome_bp:
    # depth of slug DNA: the bases of the library that are slug DNA placed on the assembly (measured on the random
    # sample), over the bases of the assembly that are not N. A library rich in organelle DNA, or an assembly with
    # masked stretches, would otherwise look deeper or shallower than it is.
    depth_all = tot_bases / A.genome_bp; fr_slug = ((sample or {}).get("fraction_of_bases") or {}).get("slug")
    depth = tot_bases * fr_slug / (A.genome_acgt or A.genome_bp) if fr_slug else depth_all
    ins = (sample or {}).get("insert_median") or 0
    # a pair straddles a join (one mate organelle, the other slug) when the join lies 80 bases or more from both
    # ends of the fragment; a fragment covers a given base of one chromosome copy (depth / 2) x ins / (2 x read)
    # times, and an insertion has two joins
    per_pair = depth * max(0.0, ins - 160) / (2 * read_len) if ins else None
    per_split = depth * max(0.0, read_len - 80) / read_len                                     # junction inside a mate, >= 40 bases each side
    sens = dict(genome_bp=A.genome_bp, genome_acgt=A.genome_acgt or None, depth=round(depth, 2), depth_all_bases=round(depth_all, 2), slug_fraction_of_bases=fr_slug, read_length=round(read_len, 1), insert_median=ins or None,
                mixed_pairs_per_insertion_per_cell=None if per_pair is None else round(per_pair, 2), split_mates_per_insertion_per_cell=round(per_split, 2))

other_taxa = {k: [dict(name=t, mates=len(v), mean_identity=round(sum(v) / len(v), 4)) for t, v in sorted(d.items(), key=lambda kv: -len(kv[1]))[:8]] for k, d in by_taxon.items()}
res = dict(kind="short", runs=runs, reads=tot_reads, bases=tot_bases, mates_examined=n_examined, organelle_stretches_given=[f"{c_}:{a_}-{b_}={k_}" for c_, a_, b_, k_ in EXTRA], other_taxa=other_taxa, plain_weight={k: round(v, 3) for k, v in weight.items()},
           partial_weight={k: round(v, 3) for k, v in pweight.items()}, organelle_pairs_examined_one_in=thinned_one_in,
           classes={c: dict(bp=round(cls_bp[c]), fraction_of_bases=round(cls_bp[c] / max(1, tot_bases), 6), pairs=round(have[c])) for c in cls_bp},
           plastid_any_identity=dict(bp=round(cp_loose), fraction_of_bases=round(cp_loose / max(1, tot_bases), 6)), donors=donors,
           chimera=dict(rates=rates, recurrent=sorted(recurrent, key=lambda x: -x["distinct_pairs"])[:60], recurrent_breaks=recurrent_breaks[:60], same_join=same_join[:100], same_join_groups=len(same_join),
                        same_join_within=NEAR, same_join_expected_by_chance=same_join_chance, split_mates=dict(split_n), split_end_piece=end_piece,
                        plain_pairs_that_look_mixed=dict(plain_but_mixed), organelle_side_is_simple_sequence=dict(simple_dropped), split_examples=split_org[:40]),
           sample=sample, sensitivity=sens)
json.dump(res, open(out, "w"), indent=1)
if A.junctions:
    with gzip.open(A.junctions, "wt") as fo:
        fo.write("molecule\tevidence\torganelle\tslug_sequence\tslug_position\tslug_strand\tslug_mapq\torganelle_sequence\torganelle_position\torganelle_strand\tbases_between\tstands_for_pairs\n")
        for k in sorted(mixed):
            for x in mixed[k]: fo.write(f'{x["spot"]}\tmixed pair\t{k}\t{x["slug"]}\t{x["slug_pos"]}\t{x["slug_strand"]}\t{x["slug_mapq"]}\t{x["other_t"]}\t{x["other_pos"]}\t{x["other_strand"]}\t\t{x["w"]:g}\n')
        for x in split_org: fo.write(f'{x["read"]}\tsplit mate\t{x["other"]}\t{x["slug"]}\t{x["slug_break"]}\t{x["slug_strand"]}\t{x["slug_mapq"]}\t{x["other_t"]}\t{x["other_break"]}\t{x["other_strand"]}\t{x["bases_between"]}\t{x["w"]:g}\n')
print(f"{tot_reads} mates, {tot_bases / 1e9:.2f} Gb screened; {n_examined} mates examined; weights of plain pairs: {dict(weight)}")
for c in sorted(cls_bp, key=lambda c: -cls_bp[c]): print(f"  {c:9s} bp {round(cls_bp[c]):12d} ({100 * cls_bp[c] / max(1, tot_bases):.4f}% of bases)  pairs {round(have[c])}")
for d in donors[:10]: print(f"    {d['name']:42s} {d['reads']:8d} mates {d['bp']:11d} bp {100 * d['fraction']:6.2f}%  identity median {d['median_identity']:.3f} (10-90%: {d['identity_q10']:.3f}-{d['identity_q90']:.3f}) {d['accession']}")
print("  mixed pairs:", {k: (v["pairs_also_with_slug"], v["pairs_with_mate"], v["rate"]) for k, v in rates.items()}, "| split mates:", dict(split_n), "| plain pairs that look mixed (not counted):", dict(plain_but_mixed))
print("  recurrent places (>= 3 different pairs within 1 kb):", dict(collections.Counter(x["kind"] for x in recurrent)), "| the same break in two molecules:", len(recurrent_breaks))
for x in sorted(recurrent, key=lambda x: -x["distinct_pairs"])[:10]: print("   ", x["kind"], x["slug"], x["start"], "pairs", x["pairs"], "distinct", x["distinct_pairs"], x["other_targets"][:2])
for x in recurrent_breaks[:10]: print("    same break:", x)
print("  piece at the end of the fragment in mates that read through a join (molecules):", end_piece)
print(f"  the same join in two molecules (both sides within {NEAR} bases):", len(same_join), "| expected by chance:", same_join_chance)
for x in same_join[:15]: print("    same join:", {k_: v for k_, v in x.items() if k_ != "evidence"})
if sample: print("  sample:", {k: v for k, v in sample.items() if k != "counts"}, sample["counts"])
if sens: print("  sensitivity:", sens)
