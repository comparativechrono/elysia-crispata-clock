#!/usr/bin/env python3
"""Do two libraries come from the same animal? Their mitochondrial genomes are compared.

usage: 58_mt_haplotypes.py PANEL.fa MT_RECORD OUT.json NAME=WORKDIR:short|hifi NAME=WORKDIR:short|hifi

MT_RECORD is the mitochondrial genome of the species in the panel (for E. crispata slugmt|OY973781.1).
For each library the reads that the second pass placed on that genome (pass2.paf; at most 60,000
read pairs or 3,000 long reads, taken evenly) are aligned to it again with base-level output
(minimap2 -a; preset sr or map-hifi), and the bases at every position are counted (samtools
mpileup, mapping quality and base quality of 20 or more where the reads have qualities). A position
counts for a library when at least 30 reads cover it and one base makes up at least 90% of them.
The two libraries differ at a position when both count there and their bases are not the same.
Animals from one population can share a mitochondrial genome, so a result without differences
proves nothing; differences at fixed positions show that the libraries were made from two animals.

OUT.json: positions compared, the number at which each library differs from the reference record,
and every position at which the two libraries differ (substitutions and single-base deletions).
"""
import sys, os, re, json, subprocess, collections, tempfile
panel, rec, out = sys.argv[1:4]; libs = []
for a in sys.argv[4:6]:
    name, rest = a.split("=", 1); wd, kind = rest.rsplit(":", 1); libs.append((name, wd, kind))
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/"
tmp = tempfile.mkdtemp(prefix="mt_", dir=os.path.dirname(os.path.abspath(out)))
ref = os.path.join(tmp, "mt.fa")
seq = subprocess.run([BIN + "samtools", "faidx", panel, rec], capture_output=True, text=True, check=True).stdout
open(ref, "w").write(">mt\n" + "".join(seq.split("\n")[1:]) + "\n")
def counts(name, wd, kind):
    names = []; seen = set()
    for l in open(os.path.join(wd, "pass2.paf")):
        c = l.split("\t", 7)
        if c[5] == rec and c[0] not in seen: seen.add(c[0]); names.append(c[0])
    cap = 60000 if kind == "short" else 3000; step = max(1, len(names) // cap)
    names = names[::step][:cap]
    nf = os.path.join(tmp, name + ".names"); open(nf, "w").write("".join(n + "\n" for n in names))
    fa = os.path.join(tmp, name + ".fa"); bam = os.path.join(tmp, name + ".bam")
    subprocess.run(f"{BIN}seqkit grep -f {nf} {os.path.join(wd, 'all_hits.fa')} 2>/dev/null > {fa}", shell=True, check=True)
    subprocess.run(f"{BIN}minimap2 -ax {'sr' if kind == 'short' else 'map-hifi'} -t 2 {ref} {fa} 2>/dev/null | {BIN}samtools sort -o {bam} - && {BIN}samtools index {bam}", shell=True, check=True)
    pile = subprocess.run([BIN + "samtools", "mpileup", "-f", ref, "-d", "100000", "-q", "20", "-Q", "20", "-B", bam], capture_output=True, text=True, check=True).stdout
    res = {}
    for l in pile.splitlines():
        c = l.split("\t"); r = c[2].upper(); s = re.sub(r"\^.", "", c[4]).replace("$", ""); cnt = collections.Counter(); i = 0
        while i < len(s):
            ch = s[i]
            if ch in "+-":                       # an insertion or deletion that begins after this base: skipped
                j = i + 1; n = ""
                while s[j].isdigit(): n += s[j]; j += 1
                i = j + int(n); continue
            if ch in ".,": cnt[r] += 1
            elif ch.upper() in "ACGT": cnt[ch.upper()] += 1
            elif ch == "*": cnt["-"] += 1
            i += 1
        res[int(c[1])] = (r, cnt)
    return res, len(names)
A, na = counts(*libs[0]); B, nb = counts(*libs[1])
compared = 0; da = db = 0; diff = []
for pos in sorted(set(A) & set(B)):
    r, ca = A[pos]; _, cb = B[pos]; ta, tb = sum(ca.values()), sum(cb.values())
    if ta < 30 or tb < 30: continue
    (ba, xa), (bb, xb) = ca.most_common(1)[0], cb.most_common(1)[0]
    if xa / ta < 0.9 or xb / tb < 0.9: continue
    compared += 1; da += ba != r; db += bb != r
    if ba != bb: diff.append(dict(position=pos, reference=r, **{libs[0][0]: ba, libs[0][0] + "_share": round(xa / ta, 3), libs[0][0] + "_depth": ta, libs[1][0]: bb, libs[1][0] + "_share": round(xb / tb, 3), libs[1][0] + "_depth": tb}))
res = dict(reference=rec, reference_length=len("".join(seq.split("\n")[1:])), positions_compared=compared, libraries=[libs[0][0], libs[1][0]],
           reads_used={libs[0][0]: na, libs[1][0]: nb}, differ_from_reference={libs[0][0]: da, libs[1][0]: db},
           differences=len(diff), substitutions=sum(1 for d in diff if "-" not in (d[libs[0][0]], d[libs[1][0]])), positions=diff)
json.dump(res, open(out, "w"), indent=1)
subprocess.run(["rm", "-rf", tmp])
print(f"{compared} positions compared; {libs[0][0]} differs from the reference at {da}, {libs[1][0]} at {db}; the two differ from each other at {len(diff)} ({res['substitutions']} substitutions)")
