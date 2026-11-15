"""Is a sequence of an assembly slug DNA, or DNA of a bacterium, phage or other organism that came with the sample?

Shared by 73_regulator_summary.py and 75_alien_candidates.py. The judgement uses the whole sequence
(for a sequence longer than 2 Mb, the 200 kb that 43_foreign_sequences.py examined around a locus):

  - the share covered by open reading frames of at least 300 bases on either strand. Prokaryote,
    phage and giant-virus DNA is 85 to 95% coding, slug DNA much less: random 20-kb windows of each
    assembly give the background. Tandem repeats and other simple sequence, which also lack stop
    codons, are left out (1-kb tiles that zlib compresses to less than 22% of their length).
    Frames are counted with the standard code, and with that of mycoplasmas (TGA codes for
    tryptophan) when the standard code leaves a sequence below 60% and the Mycoplasma code brings it
    to at least 85%, because several contaminating sequences are of Mycoplasma-like bacteria;
  - the coding stretches found by translating the sequence (43_foreign_sequences.py), counted by
    their best match: another sea slug (the species' own annotation excluded), an animal that is not
    a sea slug, or an alga, plant or microbe (foreign). Genes of slug DNA are matched best by the
    other sea slugs; genes of a contaminant are not, since the other sea slugs do not carry the same
    contaminant.

A sequence is foreign (bacterium, phage, virus or other organism) when at least 60% of it lies in
long reading frames, at least one stretch is foreign, and either at most one stretch in ten foreign
ones (at least one) is matched best by another sea slug, or foreign stretches outnumber the animal
and slug-only stretches of 43_foreign_sequences.py at least twofold; it is foreign too when its GC
content differs from the genome's by at least 0.08 and its foreign stretches outnumber those matched
best by sea slugs or other animals. A sequence shorter than 20 kb with no stretch matched best by a
sea slug or another animal cannot be placed (undetermined). Otherwise the 5-kb windows
around a locus are followed outwards while at least 60% of each lies in long reading frames; when
that run reaches an end of the sequence and the rest of it (at least 20 kb) looks like slug DNA
(median window at or below the 90th percentile of the background), the locus lies in a foreign
stretch at the end of a slug sequence (a join that only reads could test). Anything else is slug DNA.
"""
import os
import random
import statistics
import subprocess
import zlib

BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin")
DATA = os.environ.get("DATA", "/home/claude/data")
STOP_STANDARD = {"TAA", "TAG", "TGA"}
STOP_MOLLICUTES = {"TAA", "TAG"}                                 # translation table 4: TGA codes for tryptophan
COMP = str.maketrans("ACGTN", "TGCAN")
W = 5000


def orf_mask(seq, minlen=300, stops=STOP_STANDARD):
    seq = seq.upper()
    n = len(seq)
    cov = bytearray(n)
    for strand in (0, 1):
        s = seq if strand == 0 else seq.translate(COMP)[::-1]
        for f in range(3):
            start = f
            for i in range(f, n - 2, 3):
                c = s[i:i + 3]
                if c in stops or "N" in c:
                    if i - start >= minlen:
                        a, b = (start, i) if strand == 0 else (n - i, n - start)
                        cov[a:b] = b"\x01" * (b - a)
                    start = i + 3
            end = n - (n - f) % 3
            if end - start >= minlen:
                a, b = (start, end) if strand == 0 else (n - end, n - start)
                cov[a:b] = b"\x01" * (b - a)
    return cov


def simple_mask(seq, tile=1000, ratio=0.22):
    """1-kb tiles that compress to less than 22% of their length with zlib: tandem repeats and other
    simple sequence, whose long stop-free frames say nothing about genes (slug DNA compresses to
    about 24 to 30%, bacterial DNA to 32 to 34%)."""
    s = seq.upper().encode()
    mask = bytearray(len(s))
    for i in range(0, len(s), tile):
        t = s[i:i + tile]
        if len(t) >= 200 and len(zlib.compress(t)) / len(t) < ratio:
            mask[i:i + len(t)] = b"\x01" * len(t)
    return mask


def share(cov, seq):
    """Share of the bases that are neither N nor simple sequence that lie in long reading frames."""
    simple = simple_mask(seq)
    up = seq.upper()
    real = [i for i in range(len(up)) if up[i] != "N" and not simple[i]]
    return sum(cov[i] for i in real) / len(real) if real else 0.0


def best_mask(seq):
    """Reading frames under the standard code; under the Mycoplasma code when that alone makes the
    sequence look like a prokaryote (standard code below 60%, Mycoplasma code at least 85%)."""
    m = orf_mask(seq)
    s = share(m, seq)
    if s < 0.6:
        m4 = orf_mask(seq, stops=STOP_MOLLICUTES)
        s4 = share(m4, seq)
        if s4 >= 0.85:
            return m4, s4, "mycoplasma (TGA = Trp)"
    return m, s, "standard"


def genome(asm):
    return os.path.join(DATA, "asm", asm, "genome.fa.gz")


def fetch(asm, region):
    out = subprocess.run([BIN + "/samtools", "faidx", genome(asm), region], capture_output=True, text=True, check=True).stdout
    return "".join(out.split("\n")[1:])


def lengths(asm):
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(genome(asm) + ".fai")}


def background(asm, n=40, size=20000, seed=1):
    """Share of long reading frames (standard code) in n random windows of the assembly, with the
    windows themselves (sequence, start, share, G+C content)."""
    rnd = random.Random(f"{seed}-{asm}")
    big = [(c, m) for c, m in lengths(asm).items() if m > 3 * size]
    tot = sum(m for _, m in big)
    vals, rows = [], []
    while len(vals) < n:
        r = rnd.randrange(tot)
        for c, m in big:
            if r < m:
                break
            r -= m
        a = min(r, m - size - 1)
        seq = fetch(asm, f"{c}:{a + 1}-{a + size}")
        if seq.upper().count("N") > size / 2:
            continue
        up = seq.upper()
        acgt = sum(up.count(x) for x in "ACGT")
        rows.append([c, a, round(share(orf_mask(seq), seq), 4), round((up.count("G") + up.count("C")) / max(1, acgt), 4)])
        vals.append(rows[-1][2])
    vals.sort()
    return dict(windows=len(vals), median=round(statistics.median(vals), 3), p90=round(vals[int(0.9 * len(vals))], 3), list=rows)


def kinds(record):
    """Coding stretches of a 43_foreign_sequences.py record by their best match."""
    n = {"sea_slug": 0, "animal": 0, "foreign": 0, "unmatched": 0}
    for r in record["list"]:
        b = r["bits"]
        slug, met, fo = b.get("slug", 0), b.get("met", 0), max(b.get("alga", 0), b.get("microbe", 0))
        if slug and slug >= max(met, fo):
            n["sea_slug"] += 1
        elif fo > max(met, slug):
            n["foreign"] += 1
        elif met:
            n["animal"] += 1
        else:
            n["unmatched"] += 1
    return n


def judge(asm, contig, record, positions, bg, region=None, genome_gc=None):
    """Verdict on one sequence (or, for a long one, the examined region start-end, 0-based) given its
    43_foreign_sequences.py record and the midpoints of the loci of interest."""
    stretches = record["stretches"]
    n = kinds(record)
    if region:
        a, b = region
        seq = fetch(asm, f"{contig}:{a + 1}-{b}")
    else:
        a, seq = 0, fetch(asm, contig)
    cov, whole, code = best_mask(seq)
    microbial = stretches.get("microbe", 0) + stretches.get("alga", 0)
    sluglike = stretches.get("animal", 0) + stretches.get("slug_only", 0)
    win = [share(cov[i:i + W], seq[i:i + W]) for i in range(0, len(seq), W)]
    verdict, segment = None, None
    gc_off = genome_gc is not None and abs(record["gc"] - genome_gc) >= 0.08
    if whole >= 0.6 and n["foreign"] >= 1 and (n["sea_slug"] <= max(1, n["foreign"] / 10) or microbial >= 2 * sluglike):
        verdict = "foreign sequence"
    elif gc_off and n["foreign"] > n["sea_slug"] + n["animal"]:
        verdict = "foreign sequence"
    elif len(seq) < 20000 and n["sea_slug"] == 0 and n["animal"] == 0:
        verdict = "undetermined (short sequence without slug genes)"
    else:
        for p in positions:
            k = min(len(win) - 1, max(0, (p - a) // W))
            if win[k] < 0.6:
                continue
            i, j = k, k
            while i > 0 and win[i - 1] >= 0.6:
                i -= 1
            while j < len(win) - 1 and win[j + 1] >= 0.6:
                j += 1
            ends_here = region is None or (j == len(win) - 1 and b >= lengths(asm)[contig]) or (i == 0 and a == 0)
            rest = win[:i] if j == len(win) - 1 else (win[j + 1:] if i == 0 else None)
            if ends_here and rest is not None and len(rest) * W >= 20000 and statistics.median(rest) <= bg["p90"]:
                verdict = "foreign stretch at the end of a slug sequence"
                segment = dict(start=a + i * W, end=a + min(len(seq), (j + 1) * W), rest_median_orf_share=round(statistics.median(rest), 3))
                break
    if verdict is None:
        verdict = "slug DNA"
    return dict(orf_share=round(whole, 3), genetic_code=code, verdict=verdict, foreign_stretch=segment, examined=[a, a + len(seq)], stretch_kinds=n)
