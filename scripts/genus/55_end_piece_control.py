#!/usr/bin/env python3
"""Short reads: do the two ends of a mate pass the identity test equally often?

usage: 55_end_piece_control.py WORKDIR GENOME.fa.gz PANEL.fa OUT.json [PIECE=60] [MATES=20000]
       PIECE is made shorter for mates of fewer than 2 x PIECE bases (45 for mates of 100).
       P1_INDEX=PANEL_PREFIX.sr.mmi   the index the first pass mapped to (31c_screen_run.py): adds the second control

37_pair_analysis.py counts, for the mates that read through a join, which kind of DNA the mate
begins with (split_end_piece). A join inside a chromosome is cut at random and should have either
kind at the start of the mate equally often. That holds only if a piece is recognised as well at
the end of a mate as at its start: bases are read less well towards the end, and a piece of slug
DNA, which differs a little from the assembly in another animal, could then fail the 97% test more
often there than a piece of the slug's own mitochondrial DNA.

This control measures it. Of the slug mates in the random sample of the library, and of the
mitochondrial and chloroplast mates among the chosen reads, MATES are drawn at random. The first
PIECE bases and the last PIECE bases of each are aligned on their own (minimap2 -x sr -c): a slug
piece to the place of its mate in the assembly, an organelle piece to the organelle sequence its
mate matched. A piece passes when 40 bases or more align at 97% identity or better, as in the pair
analysis. From the four pass rates follows the share of joins that would be seen with the
organelle piece first if both orders were equally common:
    f = p_org(start) x p_slug(end) / (p_org(start) x p_slug(end) + p_slug(start) x p_org(end))

Second control (with P1_INDEX). A fragment that begins with a short piece of organelle DNA holds
no other organelle DNA: its first mate is that piece and then slug DNA, its second mate is slug
DNA. The first pass, which maps mates to the organelle panel only, must find the short piece or
the pair is lost. A fragment that begins with a short piece of slug DNA goes on as organelle DNA
to its far end, and its second mate is found with ease. So mates are made up that begin with
40 to (length - 40) bases of a real organelle mate and go on as a real slug mate, and are mapped
as in the first pass (minimap2 -x sr, no base-level alignment). With d the share that is found,
    expected share with the organelle piece first = f x d / (f x d + 1 - f)
"""
import sys, os, subprocess, collections, random, json
lib, genome, panel, out = sys.argv[1:5]
PIECE = int(sys.argv[5]) if len(sys.argv) > 5 else 60; NMAX = int(sys.argv[6]) if len(sys.argv) > 6 else 20000
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/"
tmp = out + ".tmp"; os.makedirs(tmp, exist_ok=True)
random.seed(1)
# the two pieces must not overlap: shorten them for short mates
_n = 0; _len = []
for l in open(os.path.join(lib, "all_sample.fa")):
    if l[0] != ">":
        _len.append(len(l.strip())); _n += 1
        if _n >= 2000: break
_med = sorted(_len)[len(_len) // 2] if _len else 2 * PIECE
if _med < 2 * PIECE: PIECE = max(40, _med // 2 - 5)
def reads_of(fa, want):
    seq = {}; name = None
    for l in open(fa):
        if l[0] == ">": name = l[1:].split()[0]; name = name if name in want else None
        elif name: seq[name] = seq.get(name, "") + l.strip()
    return seq
def pieces(path, names, seq):
    n = 0
    with open(path, "w") as fo:
        for i, r in enumerate(names):
            s = seq.get(r)
            if not s or len(s) < 2 * PIECE: continue
            fo.write(f">{i}|5\n{s[:PIECE]}\n>{i}|3\n{s[-PIECE:]}\n"); n += 1
    return n
# ---- slug mates of the random sample, each with the place it lies at
best = {}
for l in open(os.path.join(lib, "sample.paf")):
    c = l.split("\t", 12)
    if "|" in c[5] or int(c[11]) < 20 or int(c[9]) < 80: continue
    if c[0] not in best or int(c[9]) > best[c[0]][0]: best[c[0]] = (int(c[9]), c[5], int(c[7]), int(c[8]), int(c[2]), int(c[3]), int(c[1]), c[4])
names = sorted(best); random.shuffle(names); names = names[:NMAX]
regs = []
for r in names:
    _, t, ts, te, qs, qe, ql, st = best[r]
    regs.append(f"{t}:{max(1, ts - (qs if st == '+' else ql - qe) - 30 + 1)}-{te + ((ql - qe) if st == '+' else qs) + 30}")
open(f"{tmp}/regions", "w").write("\n".join(regs) + "\n")
fa = subprocess.run([BIN + "samtools", "faidx", genome, "-r", f"{tmp}/regions"], capture_output=True, text=True, check=True).stdout
with open(f"{tmp}/slug.ref.fa", "w") as fo:
    i = -1
    for l in fa.split("\n"):
        if l.startswith(">"): i += 1; fo.write(f">loc{i}\n")
        elif l: fo.write(l.upper() + "\n")
slug_seq = reads_of(os.path.join(lib, "all_sample.fa"), set(names))
n_mates = {"slug": pieces(f"{tmp}/slug.pieces.fa", names, slug_seq)}
org_seq = {}
# ---- mitochondrial and chloroplast mates among the chosen reads
org = {"slugmt": {}, "cp": {}}
for l in open(os.path.join(lib, "pass2.paf")):
    c = l.split("\t", 12); k = c[5].split("|")[0]
    if k in org and int(c[9]) >= 80 and int(c[9]) >= 0.97 * int(c[10]) and len(org[k]) < 400000: org[k].setdefault(c[0], c[5])
for k in org:
    nk = sorted(org[k]); random.shuffle(nk); nk = nk[:NMAX]
    open(f"{tmp}/{k}.refs", "w").write("\n".join(sorted({org[k][r] for r in nk})) + "\n")
    with open(f"{tmp}/{k}.ref.fa", "w") as fo: fo.write(subprocess.run([BIN + "samtools", "faidx", panel, "-r", f"{tmp}/{k}.refs"], capture_output=True, text=True, check=True).stdout if nk else "")
    org_seq[k] = reads_of(os.path.join(lib, "all_hits.fa"), set(nk))
    n_mates[k] = pieces(f"{tmp}/{k}.pieces.fa", nk, org_seq[k])
# ---- align the pieces and count
res = {}
for k in ("slug", "slugmt", "cp"):
    if not n_mates[k]: continue
    paf = subprocess.run(f"{BIN}minimap2 -x sr -c --secondary=no -t 2 {tmp}/{k}.ref.fa {tmp}/{k}.pieces.fa 2>/dev/null", shell=True, capture_output=True, text=True).stdout
    ok = collections.defaultdict(dict)
    for l in paf.split("\n"):
        if not l: continue
        c = l.split("\t"); i, e = c[0].split("|")
        if k == "slug" and c[5] != f"loc{i}": continue          # a slug piece counts at the place of its mate only
        ok[i][e] = ok[i].get(e, False) or (int(c[9]) >= 40 and int(c[9]) >= 0.97 * int(c[10]))
    res[k] = dict(mates=n_mates[k], start_passes=sum(1 for v in ok.values() if v.get("5")), end_passes=sum(1 for v in ok.values() if v.get("3")))
    print(f'{k:7s} mates {n_mates[k]:6d}  first {PIECE} bases pass: {res[k]["start_passes"]} ({100 * res[k]["start_passes"] / n_mates[k]:.1f}%)   last {PIECE} bases pass: {res[k]["end_passes"]} ({100 * res[k]["end_passes"] / n_mates[k]:.1f}%)')
exp = {}
for k in ("slugmt", "cp"):
    if k in res and "slug" in res and res[k]["mates"] >= 200:
        o5, o3 = res[k]["start_passes"] / res[k]["mates"], res[k]["end_passes"] / res[k]["mates"]; s5, s3 = res["slug"]["start_passes"] / res["slug"]["mates"], res["slug"]["end_passes"] / res["slug"]["mates"]
        if o5 * s3 + s5 * o3 > 0: exp[k] = round(o5 * s3 / (o5 * s3 + s5 * o3), 4)
        print(f"  share of joins expected with the {k} piece first, if both orders were equally common: {100 * exp.get(k, float('nan')):.1f}%")
# ---- second control: does the first pass find a mate that begins with a short piece of organelle DNA?
found = {}; exp_found = {}
P1 = os.environ.get("P1_INDEX")
if P1 and os.path.exists(P1):
    S_ = [v for v in slug_seq.values() if len(v) >= _med]
    for k in ("slugmt", "cp"):
        O_ = [v for v in org_seq.get(k, {}).values() if len(v) >= _med]
        if len(O_) < 500 or len(S_) < 500: continue
        with open(f"{tmp}/sim.{k}.fa", "w") as fo:
            for i in range(NMAX):
                x = random.randint(40, _med - 40); o = random.choice(O_); s_ = random.choice(S_)
                fo.write(f">{i}\n{o[:x]}{s_[x:_med]}\n")
        hit = set(subprocess.run(f"{BIN}minimap2 -x sr -t 2 {P1} {tmp}/sim.{k}.fa 2>/dev/null | cut -f1", shell=True, capture_output=True, text=True).stdout.split())
        found[k] = round(len(hit) / NMAX, 4)
        if k in exp: exp_found[k] = round(exp[k] * found[k] / (exp[k] * found[k] + 1 - exp[k]), 4)
        print(f"  made-up mates that begin with a short {k} piece: the first pass finds {100 * found[k]:.0f}%; share of joins then expected with the {k} piece first: {100 * exp_found.get(k, float('nan')):.1f}%")
json.dump(dict(piece=PIECE, pass_rule="40 bases or more aligned at 97% identity or better", mates=res, expected_share_organelle_first=exp,
               first_pass_finds_short_organelle_piece=found, expected_share_organelle_first_found=exp_found), open(out, "w"), indent=1)
import shutil; shutil.rmtree(tmp, ignore_errors=True)
