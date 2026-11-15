"""Fast reading-frame and repeat measures for every sequence of an assembly.

Used by 76_foreign_census.py. stats() gives the same share of long open reading frames as
foreign_dna.best_mask (frames of at least 300 bases on either strand, N ends a frame, 1-kb tiles
that zlib compresses to less than 22% of their length left out), with numpy instead of Python
loops, under the standard code and under the Mycoplasma code (TGA = Trp); and the 5-kb windows
(counted from position 0, as in foreign_dna.judge) that are coding at either end of the sequence:
a window is coding when its standard share is at least 0.6 or its Mycoplasma share at least 0.85.
distinct() gives the share of 12-base words that are distinct, low for tandem arrays of any period.

As a script: orfscan.py GENOME.fa.gz OUT.tsv writes one row per sequence.
"""
import os, sys, zlib, subprocess
import numpy as np

BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin")
LUT = np.full(256, 4, np.uint8)
for _i, _ch in enumerate(b"ACGT"):
    LUT[_ch] = _i; LUT[_ch + 32] = _i
# stop codons as 16a+4b+c with A=0 C=1 G=2 T=3, read on the forward strand; reverse-strand stops
# are the reverse complements read forward (TTA, CTA, TCA)
STD = [48, 50, 56, 60, 28, 52]            # TAA TAG TGA | TTA CTA TCA
STD_F, STD_R = STD[:3], STD[3:]
M4_F, M4_R = [48, 50], [60, 28]           # TGA codes for tryptophan
W = 5000


def _cover(codes, n, stops, minc=100):
    d = np.zeros(n + 1, np.int32)
    for f in range(3):
        m = (n - f) // 3
        if m <= 0:
            continue
        a = codes[f:f + 3 * m].reshape(m, 3).astype(np.int32)
        isstop = (a == 4).any(1) | np.isin(a[:, 0] * 16 + a[:, 1] * 4 + a[:, 2], stops)
        edges = np.concatenate(([-1], np.flatnonzero(isstop), [m]))
        st, en = edges[:-1] + 1, edges[1:]
        ok = (en - st) >= minc
        np.add.at(d, f + 3 * st[ok], 1)
        np.add.at(d, f + 3 * en[ok], -1)
    return np.cumsum(d[:-1]) > 0


def _simple(seqb, n):
    m = np.zeros(n, bool)
    up = seqb.upper()
    for i in range(0, n, 1000):
        t = up[i:i + 1000]
        if len(t) >= 200 and len(zlib.compress(t)) / len(t) < 0.22:
            m[i:i + len(t)] = True
    return m


def stats(seqb):
    n = len(seqb)
    codes = LUT[np.frombuffer(seqb, np.uint8)]
    real = (codes != 4) & ~_simple(seqb, n)
    acgt = codes[codes != 4]
    gc = float(((acgt == 1) | (acgt == 2)).sum()) / max(1, len(acgt))
    nr = int(real.sum())
    out = dict(length=n, acgt=int(len(acgt)), gc=gc, std=0.0, m4=0.0, windows=0, head=0, tail=0, head_rest=None, tail_rest=None)
    if nr == 0:
        return out
    std = (_cover(codes, n, STD_F) | _cover(codes, n, STD_R)) & real
    m4 = (_cover(codes, n, M4_F) | _cover(codes, n, M4_R)) & real
    k = (n + W - 1) // W
    pad = k * W - n
    wsum = lambda x: np.pad(x, (0, pad)).reshape(k, W).sum(1)
    rw = wsum(real).astype(float)
    ws = np.where(rw > 0, wsum(std) / np.maximum(rw, 1), 0.0)
    w4 = np.where(rw > 0, wsum(m4) / np.maximum(rw, 1), 0.0)
    coding = (ws >= 0.6) | (w4 >= 0.85)
    head = 0
    while head < k and coding[head]:
        head += 1
    tail = 0
    while tail < k and coding[k - 1 - tail]:
        tail += 1
    out.update(std=float(std.sum()) / nr, m4=float(m4.sum()) / nr, windows=k, head=head, tail=tail,
               head_rest=float(np.median(ws[head:])) if head < k else None,
               tail_rest=float(np.median(ws[:k - tail])) if tail < k else None)
    return out


def distinct(seqb, k=12):
    c = LUT[np.frombuffer(seqb, np.uint8)].astype(np.int64)
    n = len(c)
    if n < k:
        return 1.0
    bad = (c == 4).astype(np.int64)
    c = np.where(c == 4, 0, c)
    h = np.zeros(n - k + 1, np.int64)
    for j in range(k):
        h = h * 4 + c[j:n - k + 1 + j]
    h = h[np.convolve(bad, np.ones(k, np.int64), "valid") == 0]
    return len(np.unique(h)) / len(h) if len(h) else 1.0


def records(path):
    p = subprocess.Popen([BIN + "/bgzip", "-dc", path], stdout=subprocess.PIPE)
    name, buf = None, []
    for line in p.stdout:
        if line[:1] == b">":
            if name is not None:
                yield name, b"".join(buf)
            name, buf = line[1:].split()[0].decode(), []
        else:
            buf.append(line.rstrip())
    if name is not None:
        yield name, b"".join(buf)
    p.wait()


COLUMNS = ["sequence", "length", "acgt", "gc", "orf_standard", "orf_mycoplasma", "windows", "head_run", "head_rest_median", "tail_run", "tail_rest_median"]


def scan(genome, out):
    f = lambda x: "" if x is None else f"{x:.4f}"
    with open(out, "w") as fo:
        fo.write("\t".join(COLUMNS) + "\n")
        for name, s in records(genome):
            r = stats(s)
            fo.write(f"{name}\t{r['length']}\t{r['acgt']}\t{r['gc']:.4f}\t{r['std']:.4f}\t{r['m4']:.4f}\t{r['windows']}\t{r['head']}\t{f(r['head_rest'])}\t{r['tail']}\t{f(r['tail_rest'])}\n")


if __name__ == "__main__":
    scan(sys.argv[1], sys.argv[2])
