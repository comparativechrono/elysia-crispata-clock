#!/usr/bin/env python3
"""Which parts of an assembly are absent from the DNA of eggs?

usage: 59_germline_coverage.py SAMPLE.paf GENOME.fa.gz OUT_PREFIX [--window 5000] [--min-expected 15]

SAMPLE.paf holds the alignments of a random sample of the read pairs of an egg library to the
assembly (33_pass2.sh: sample.paf). Eggs contain the nuclear and mitochondrial genome of the
animal and nothing that it has eaten, so a stretch of the assembly on which no egg read lies,
although reads were expected, is not part of the germline genome: kleptoplast DNA, bacteria or
other DNA that was sequenced along with the adult.

Every assembly sequence is cut into windows (default 5 kb; a shorter sequence is one window).
The number of reads expected in a window is its number of unmasked bases (not N) times the
density of reads on the assembly as a whole (the median over all windows that are free of N).
A window counts as empty when at least MIN_EXPECTED reads are expected and fewer than 5% of them
are found. Neighbouring empty windows are joined into a stretch, and a window too short to be
tested is joined to an empty neighbour on the same sequence if no read lies on it either.

Reads are counted at the middle of their alignment (>= 60 aligned bases at >= 95% identity), with
any mapping quality: a read from a repeat is placed on one of its copies, which is enough to show
that the sequence occurs in eggs.

Output: OUT_PREFIX.tsv  one line per stretch: sequence, length of the sequence, start, end, bases,
                        G+C, reads found, reads expected, whether the stretch is the whole sequence
        OUT_PREFIX.json the totals
"""
import sys, gzip, json, collections, argparse, statistics
ap = argparse.ArgumentParser()
ap.add_argument("paf"); ap.add_argument("genome"); ap.add_argument("out")
ap.add_argument("--window", type=int, default=5000); ap.add_argument("--min-expected", type=float, default=15)
A = ap.parse_args()
W = A.window
# ---- the assembly: per window, the bases that are not N and their G+C
win = {}          # sequence -> list of [acgt, gc]
length = {}
def close(name, seq):
    n = len(seq); length[name] = n; ws = []
    for a in range(0, n, W):
        w = seq[a:a + W].upper()
        gc = w.count("G") + w.count("C"); at = w.count("A") + w.count("T")
        ws.append([gc + at, gc])
    if len(ws) > 1 and n % W and n % W < W // 2:        # a short last window joins the one before it
        ws[-2][0] += ws[-1][0]; ws[-2][1] += ws[-1][1]; ws.pop()
    win[name] = ws
name = None; buf = []
with (gzip.open if A.genome.endswith(".gz") else open)(A.genome, "rt") as f:
    for line in f:
        if line[0] == ">":
            if name: close(name, "".join(buf))
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if name: close(name, "".join(buf))
def widx(name, pos): return min(pos // W, len(win[name]) - 1)
# ---- reads
cnt = {n: [0] * len(ws) for n, ws in win.items()}; n_reads = 0
for line in open(A.paf):
    c = line.split("\t", 12)
    t = c[5]
    if t not in cnt: continue                      # panel sequences
    al = int(c[10])
    if al < 60 or int(c[9]) / al < 0.95: continue
    cnt[t][widx(t, (int(c[7]) + int(c[8])) // 2)] += 1; n_reads += 1
# ---- density: reads per unmasked base, the median over the full windows without N
dens = [cnt[n][i] / ws[i][0] for n, ws in win.items() for i in range(len(ws)) if ws[i][0] >= W]
density = statistics.median(dens) if dens else 0.0
if not density: sys.exit("no reads on the assembly")
# ---- empty windows, joined into stretches
rows = []
for n, ws in win.items():
    state = []
    for i, (acgt, gc) in enumerate(ws):
        exp = acgt * density
        state.append("empty" if exp >= A.min_expected and cnt[n][i] < 0.05 * exp else "short" if exp < A.min_expected else "covered")
    # a window too short to be tested joins an empty neighbour if no read lies on it either
    joins = [state[x] == "empty" or (state[x] == "short" and cnt[n][x] == 0) for x in range(len(ws))]
    i = 0
    while i < len(ws):
        if state[i] != "empty": i += 1; continue
        j = i
        while j + 1 < len(ws) and joins[j + 1]: j += 1
        k = i
        while k - 1 >= 0 and joins[k - 1] and state[k - 1] == "short": k -= 1
        a = k * W; b = length[n] if j == len(ws) - 1 else (j + 1) * W
        acgt = sum(ws[x][0] for x in range(k, j + 1)); gc = sum(ws[x][1] for x in range(k, j + 1)); got = sum(cnt[n][x] for x in range(k, j + 1))
        rows.append(dict(sequence=n, length=length[n], start=a, end=b, bases=acgt, gc=round(gc / acgt, 4) if acgt else None, reads=got, expected=round(acgt * density, 1), whole=(a == 0 and b == length[n])))
        i = j + 1
# a sequence too short to be tested as a whole is not listed; its bases are counted as untested
untested = sum(ws[i][0] for n, ws in win.items() for i in range(len(ws)) if ws[i][0] * density < A.min_expected)
tested = sum(ws[i][0] for n, ws in win.items() for i in range(len(ws)) if ws[i][0] * density >= A.min_expected)
rows.sort(key=lambda r: -r["bases"])
with open(A.out + ".tsv", "w") as fo:
    fo.write("sequence\tsequence_length\tstart\tend\tbases\tgc\treads_found\treads_expected\twhole_sequence\n")
    for r in rows: fo.write(f'{r["sequence"]}\t{r["length"]}\t{r["start"]}\t{r["end"]}\t{r["bases"]}\t{r["gc"]}\t{r["reads"]}\t{r["expected"]}\t{"yes" if r["whole"] else "no"}\n')
res = dict(reads_counted=n_reads, window=W, reads_per_kb=round(1000 * density, 3), min_expected=A.min_expected,
           assembly_bases=sum(length.values()), bases_tested=tested, bases_untested=untested,
           stretches=len(rows), stretches_bases=sum(r["bases"] for r in rows), whole_sequences=sum(1 for r in rows if r["whole"]), whole_sequences_bases=sum(r["bases"] for r in rows if r["whole"]),
           sequences_with_a_stretch=len({r["sequence"] for r in rows}))
json.dump(res, open(A.out + ".json", "w"), indent=1)
print(res)
