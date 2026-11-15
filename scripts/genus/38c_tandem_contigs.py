#!/usr/bin/env python3
"""Contigs assembled from long reads (38b_plastid_contigs_long.sh) that are tandem copies of a short unit.

usage: 38c_tandem_contigs.py WORKDIR [PROTEINS.faa]

A library made by whole-genome amplification copies small circular DNA molecules many times over
(rolling-circle amplification): the mitochondrial genome, and any other small circle in the
sample. A read from such a molecule holds several copies of the circle one after the other, and
contigs assembled from these reads are tandem arrays of one unit. Such contigs are not pieces of
a plastid genome even where part of the unit resembles one.

For every contig of WORKDIR/kcp.contigs.fa the script finds the distance at which most of its
21-base words recur (within 30 bases). A contig is a tandem array when that distance is at most
10 kb, the contig holds at least two copies, and at least 30% of its words recur at it. Tandem
arrays are taken out of kcp.contigs.fa and written as circular elements:
  circ.contigs.fa    the arrays, named circ|NAME|ACCESSION|IDENTITY (the panel record that the contig
                     matches best, and its identity to it, as given by 38b)
  circ.units.fa      one unit of each array
  circ.contigs.tsv   contig, length, unit length, copies, share of words that recur at the unit
                     length, family (units that share at least half of their 21-base words, whatever
                     point of the circle they begin at), best panel match, identity to it, share of
                     the contig that matches panel plastid sequences, and (with PROTEINS.faa) the
                     proteins that match the unit (tblastn on two copies of the unit, E <= 1e-8)
  own.contigs.fa     circ.contigs.fa and what is left of kcp.contigs.fa: the further reference for
                     33_pass2.sh
  contigs.rename.tsv old name, new name, for alignments made before the arrays were recognised
  circ.reads.json    whether single reads hold such arrays: of the reads of 3 kb or more that were assembled (kcp.sub.fa),
                     the number in which at least 10% of the 17-base words recur at one distance (within 50 bases),
                     and the number in which that distance lies within the range of the unit lengths (5% either way)
kcp.contigs.fa and kcp.contigs.tsv keep the contigs that are not tandem arrays (the first version is
kept as kcp.contigs.before_38c.fa / .tsv). 34_read_analysis.py and 46_candidates.py count the
sequences named circ| as a kind of DNA of their own.
"""
import sys, os, collections, subprocess, shutil
W = sys.argv[1]; PROT = sys.argv[2] if len(sys.argv) > 2 else None
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/"
K = 21; MAXUNIT = 10000; MINSHARE = 0.30
fa = os.path.join(W, "kcp.contigs.fa"); tsv = os.path.join(W, "kcp.contigs.tsv")
first_fa = os.path.join(W, "kcp.contigs.before_38c.fa"); first_tsv = os.path.join(W, "kcp.contigs.before_38c.tsv")
if not os.path.exists(first_fa): shutil.copyfile(fa, first_fa); shutil.copyfile(tsv, first_tsv)
seqs = {}; name = None
for l in open(first_fa):
    if l[0] == ">": name = l[1:].strip(); seqs[name] = ""
    else: seqs[name] += l.strip().upper()
info = {}
lines = open(first_tsv).read().rstrip("\n").split("\n"); head = lines[0]
for l in lines[1:]:
    f = l.split("\t"); info[f[0]] = f
def period(s):
    pos = collections.defaultdict(list)
    for i in range(len(s) - K + 1): pos[s[i:i + K]].append(i)
    d = collections.Counter()
    for ps in pos.values():
        for a, b in zip(ps, ps[1:]): d[b - a] += 1
    best = (0, 0)
    for dist, _ in d.most_common(60):
        tot = sum(v for x, v in d.items() if abs(x - dist) <= 30)
        if tot > best[1]: best = (dist, tot)
    return best[0], best[1] / max(1, len(s) - K + 1)
arrays = {}; rest = []
for n, s in seqs.items():
    per, share = period(s)
    if per and per <= MAXUNIT and len(s) >= 2 * per and share >= MINSHARE: arrays[n] = dict(unit=per, share=share, copies=len(s) / per)
    else: rest.append(n)
# families of units: shared words, so that the point at which a unit begins does not matter
words = {n: {seqs[n][i:i + K] for i in range(arrays[n]["unit"] + K)} for n in arrays}
fam = {}; nf = 0
for n in sorted(arrays, key=lambda n: -len(seqs[n])):
    hit = next((m for m in fam if len(words[n] & words[m]) >= 0.5 * min(len(words[n]), len(words[m]))), None)
    if hit is None: nf += 1; fam[n] = nf
    else: fam[n] = fam[hit]
# proteins on the units
genes = collections.defaultdict(dict)
units = os.path.join(W, "circ.units.fa")
new = {n: "circ|" + n.split("|", 1)[1] for n in arrays}
with open(units, "w") as fo:
    for n in arrays: fo.write(f">{new[n]}\n{seqs[n][:arrays[n]['unit']]}\n")
if PROT and arrays:
    two = os.path.join(W, "circ.units.twice.fa")
    with open(two, "w") as fo:
        for i, n in enumerate(arrays): fo.write(f">u{i}\n{seqs[n][:arrays[n]['unit']] * 2}\n")
    idx = {f"u{i}": n for i, n in enumerate(arrays)}
    out = subprocess.run([BIN + "tblastn", "-query", PROT, "-subject", two, "-evalue", "1e-8", "-outfmt", "6 qseqid sseqid pident length qlen bitscore"], capture_output=True, text=True, check=True).stdout
    for l in out.splitlines():
        q, s, pid, ln, ql, bs = l.split("\t"); n = idx[s]
        if q not in genes[n] or float(bs) > genes[n][q][0]: genes[n][q] = (float(bs), float(pid), int(ln), int(ql))
    os.remove(two)
with open(os.path.join(W, "circ.contigs.fa"), "w") as fo:
    for n in arrays: fo.write(f">{new[n]}\n{seqs[n]}\n")
with open(os.path.join(W, "circ.contigs.tsv"), "w") as fo:
    fo.write("contig\tlength\tunit_length\tcopies\tshare_of_words_recurring\tfamily\tbest_panel_match\tidentity\tplastid_share\tproteins\n")
    for n in sorted(arrays, key=lambda n: (fam[n], -len(seqs[n]))):
        a = arrays[n]; f = info[n]
        g = "; ".join(f"{q} ({v[1]:.0f}% of {v[2]} residues; protein of {v[3]})" for q, v in sorted(genes[n].items(), key=lambda kv: -kv[1][0])[:3])
        fo.write(f"{new[n]}\t{len(seqs[n])}\t{a['unit']}\t{a['copies']:.1f}\t{a['share']:.2f}\t{fam[n]}\t{f[3]}\t{f[4]}\t{f[5]}\t{g}\n")
with open(fa, "w") as fo:
    for n in rest: fo.write(f">{n}\n{seqs[n]}\n")
with open(tsv, "w") as fo:
    fo.write(head + "\n")
    for n in rest: fo.write("\t".join(info[n]) + "\n")
with open(os.path.join(W, "own.contigs.fa"), "w") as fo:
    for n in arrays: fo.write(f">{new[n]}\n{seqs[n]}\n")
    for n in rest: fo.write(f">{n}\n{seqs[n]}\n")
with open(os.path.join(W, "contigs.rename.tsv"), "w") as fo:
    for n in arrays: fo.write(f"{n}\t{new[n]}\n")
sub = os.path.join(W, "kcp.sub.fa")
if arrays and os.path.exists(sub):
    import json
    lo = 0.95 * min(a["unit"] for a in arrays.values()); hi = 1.05 * max(a["unit"] for a in arrays.values()); n = t = u = 0; s = None
    def test(s):
        pos = collections.defaultdict(list)
        for i in range(len(s) - 16): pos[s[i:i + 17]].append(i)
        d = collections.Counter()
        for ps in pos.values():
            for a, b in zip(ps, ps[1:]): d[(b - a) // 50 * 50] += 1
        if not d: return None
        dist, c = d.most_common(1)[0]
        return dist if c >= 0.10 * (len(s) - 16) else None
    for l in open(sub):
        if l[0] == ">": continue
        s = l.strip().upper()
        if len(s) < 3000: continue
        n += 1; per = test(s)
        if per is not None: t += 1; u += lo <= per <= hi
    json.dump(dict(reads_of_3_kb_or_more=n, with_a_tandem_array=t, with_an_array_of_the_unit_length=u), open(os.path.join(W, "circ.reads.json"), "w"), indent=1)
    print(f"single reads of 3 kb or more: {n}; {t} hold a tandem array, {u} of them with the unit length of the contigs")
distinct = len(set().union(*words.values())) if words else 0
print(f"{len(seqs)} contigs: {len(arrays)} are tandem arrays of a unit of {min((a['unit'] for a in arrays.values()), default=0)} to {max((a['unit'] for a in arrays.values()), default=0)} bases "
      f"({nf} families; {distinct} distinct {K}-base words in their units), {len(rest)} are not")
