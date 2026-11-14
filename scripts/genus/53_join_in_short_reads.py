#!/usr/bin/env python3
"""Short reads: is a known join between chloroplast DNA and a repeat of the slug genome present
in this library?

usage: 53_join_in_short_reads.py WORKDIR GENOME.fa.gz OUT.json REPEAT_REGION PLASTID_POSITION [PLASTID_POSITION ...]
  WORKDIR         a short-read library after the second pass (pass2.paf, all_hits.fa)
  REPEAT_REGION   a stretch of the assembly that holds several units of the repeat, e.g.
                  RQTK01001031.1:1604-3400 (five units of the 358-base satellite of E. chlorotica)
  PLASTID_POSITION  where on the plastid genome the long reads turn from chloroplast DNA into the repeat

A join into repeated DNA cannot be looked up by its place in the assembly, and the mapping quality
of the repeat is nil. So the search starts from the plastid side. For each position, every mate is
taken whose alignment to a plastid genome (>= 30 bases, >= 97% identity) ends within 20 bases of
the position while the mate itself goes on for 25 bases or more. What the mate goes on with is
compared with the repeat (blastn, word size 7; >= 20 bases at >= 90%). The same is done at
control positions 300 to 1,500 bases to either side, which shows how many such mates the library
makes anywhere on the plastid genome: a join that exists in the animals of this library would stand
far above them, with about half as many mates as the genome is covered.
"""
import sys, os, json, subprocess, collections
W, genome, out, region = sys.argv[1:5]; positions = [int(x) for x in sys.argv[5:]]
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin")
SLACK, CLIP = 20, 25
offsets = [0] + [s * d for d in (300, 600, 900, 1200, 1500) for s in (-1, 1)]
want = collections.defaultdict(list)          # read -> (position, offset, clip start, clip end, side)
n_cp = 0
for line in open(os.path.join(W, "pass2.paf")):
    c = line.rstrip("\n").split("\t")
    if not c[5].startswith(("cp|", "asmcp|")) or int(c[10]) < 30: continue
    de = next((float(x[5:]) for x in c[12:] if x.startswith("de:f:")), 0.0)
    if de > 0.03: continue
    n_cp += 1
    qlen, qs, qe, strand, ts, te = int(c[1]), int(c[2]), int(c[3]), c[4], int(c[7]), int(c[8])
    ends = []
    # the plastid position at each end of the alignment, and the bases of the mate beyond that end
    if qs >= CLIP: ends.append((ts if strand == "+" else te, 0, qs))
    if qlen - qe >= CLIP: ends.append((te if strand == "+" else ts, qe, qlen))
    for pos, a, b in ends:
        for p in positions:
            for o in offsets:
                if abs(pos - (p + o)) <= SLACK: want[c[0]].append((p, o, a, b))
clips = {}
if want:
    name = None
    with open(os.path.join(W, "all_hits.fa")) as fi, open(out + ".clips.fa", "w") as fo:
        for line in fi:
            if line[0] == ">": name = line[1:].split()[0]; continue
            if name in want:
                s = line.strip()
                for k, (p, o, a, b) in enumerate(want[name]):
                    fo.write(f">{name}|{k}\n{s[a:b]}\n"); clips[f"{name}|{k}"] = (p, o, b - a)
                name = None
rep = out + ".repeat.fa"
with open(rep, "w") as fo: fo.write(subprocess.run([BIN + "/samtools", "faidx", genome, region], capture_output=True, text=True, check=True).stdout)
hit = set()
if clips:
    r = subprocess.run([BIN + "/blastn", "-task", "blastn-short", "-word_size", "7", "-query", out + ".clips.fa", "-subject", rep, "-evalue", "1", "-dust", "no",
                        "-outfmt", "6 qseqid pident length"], capture_output=True, text=True)
    for l in r.stdout.splitlines():
        q, pid, ln = l.split("\t")
        if float(pid) >= 90 and int(ln) >= 20: hit.add(q)
res = dict(repeat=region, mates_with_plastid_alignment=n_cp, positions=[])
for p in positions:
    row = dict(plastid_position=p, at_the_join=dict(mates_going_on=0, into_the_repeat=0, molecules_into_the_repeat=0), controls=[])
    for o in offsets:
        ks = [k for k, v in clips.items() if v[0] == p and v[1] == o]
        spots = {k.split("|")[0].rsplit(".", 1)[0] for k in ks if k in hit}
        d = dict(mates_going_on=len({k.split("|")[0] for k in ks}), into_the_repeat=len({k.split("|")[0] for k in ks if k in hit}), molecules_into_the_repeat=len(spots))
        if o == 0: row["at_the_join"] = d; row["reads_into_the_repeat"] = sorted({k.split("|")[0] for k in ks if k in hit})[:50]
        else: row["controls"].append(dict(offset=o, **d))
    cs = [x["into_the_repeat"] for x in row["controls"]]; row["controls_mean_into_the_repeat"] = round(sum(cs) / len(cs), 2); row["controls_max_into_the_repeat"] = max(cs)
    res["positions"].append(row)
    print(f'plastid position {p}: {row["at_the_join"]["mates_going_on"]} mates go on past it, {row["at_the_join"]["into_the_repeat"]} of them into the repeat '
          f'({row["at_the_join"]["molecules_into_the_repeat"]} molecules); at {len(cs)} control positions: mean {row["controls_mean_into_the_repeat"]}, most {row["controls_max_into_the_repeat"]}')
json.dump(res, open(out, "w"), indent=1)
for f in (out + ".clips.fa", rep):
    if os.path.exists(f): os.remove(f)
