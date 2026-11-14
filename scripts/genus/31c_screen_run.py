#!/usr/bin/env python3
"""First pass of the read screen for a whole sequencing run of any size, in pieces.

usage: 31c_screen_run.py OUTDIR PANEL_PREFIX RUN long|short SPOTS [options]
  --src PATH|URL    the run: a stored .sra file, or its address in the SRA open-data bucket
                    (default), from which pieces are then read directly by several streams
  --chunk N         spots per piece (default 50000 long, 2000000 short; a fifth of that when remote)
  --streams P       conversions running at once (default 1 for a stored file; for a remote one 12 for
                    long reads and 6 for short reads, of which a piece is much larger)
  --plain-keep K    short reads: keep one in K of the pairs that are organelle DNA at both ends
                    (default 1 = all). The other pairs with a match are always kept, unless
  --match-keep K    short reads, for a library rich in organelle DNA: keep one in K also of the pairs
                    in which the mates match one kind of organelle genome and nothing else, but not
                    over their whole length, or only one mate matches. These are the pairs among
                    which joins of organelle DNA to slug DNA are found, together with far more pairs
                    whose ends were read badly. Pairs that match a place under test (a bait) or two
                    kinds of sequence are always kept. With the same K as --plain-keep, the kept
                    pairs are an even sample, one in K, of all pairs that match an organelle genome.
  --sample S        keep one spot in S of ALL spots as a random sample of the library (default 0 = none)
  --max-spots M     stop after M spots (default: the whole run)

The run is converted piece by piece (fastq-dump -N/-X) and every read is mapped to the reference
panel with minimap2, without base-level alignment:
  long    minimap2 -x map-ont on PANEL_PREFIX.ont.mmi (k 15, w 10)
  short   minimap2 -x sr on PANEL_PREFIX.sr.mmi; the two mates of a pair are RUN.n.1 and RUN.n.2
Reads are chosen as follows.
  long    every read with a match is set aside while the run is read; when the run is finished the
          rules of 31_select_hits.py are applied to all matches at once (reads that match only "hot"
          windows of the panel, or only ribosomal DNA at slug-like identity, are thinned)
  short   a pair is "plain" when both mates match the same kind of organelle sequence over at least
          80% of their length. Plain pairs are counted, and one in K is kept. Every other pair in
          which a mate has a match is kept: these include all pairs that could join organelle DNA
          to something else.
Nothing but the chosen reads is stored. The work can be interrupted and started again: finished
pieces are remembered in OUTDIR/RUN.parts/.
Output: RUN.hits.fa, RUN.count (reads, bases), RUN.select.json, RUN.sample.fa (with --sample),
RUN.p1.paf (long: all matches; short: matches of the kept pairs), RUN.done
"""
import sys, os, subprocess, threading, json, time, collections, argparse, glob, shutil

ap = argparse.ArgumentParser()
ap.add_argument("outdir"); ap.add_argument("panel"); ap.add_argument("run"); ap.add_argument("kind", choices=["long", "short"]); ap.add_argument("spots", type=int)
ap.add_argument("--src"); ap.add_argument("--chunk", type=int); ap.add_argument("--streams", type=int); ap.add_argument("--plain-keep", type=int, default=1); ap.add_argument("--match-keep", type=int, default=1)
ap.add_argument("--sample", type=int, default=0); ap.add_argument("--max-spots", type=int, default=0)
A = ap.parse_args()
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT, RUN, KIND = A.outdir, A.run, A.kind
src = A.src or f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{RUN}/{RUN}"
remote = src.startswith("http")
chunk = A.chunk or ((50000 if KIND == "long" else 2000000) // (5 if remote else 1))
streams = A.streams or ((12 if KIND == "long" else 6) if remote else 1)
total = min(A.spots, A.max_spots) if A.max_spots else A.spots
os.makedirs(OUT, exist_ok=True)
done_flag = os.path.join(OUT, RUN + ".done")
if os.path.exists(done_flag): print(RUN, "already screened"); sys.exit(0)
PARTS = os.path.join(OUT, RUN + ".parts"); os.makedirs(PARTS, exist_ok=True)
# a run that was interrupted goes on with the pieces it began with, whatever the options say now
_lay = os.path.join(PARTS, "layout.json")
if os.path.exists(_lay): chunk = json.load(open(_lay))["chunk"]
elif glob.glob(os.path.join(PARTS, "*.json")): sys.exit(f"{PARTS}: pieces of unknown size from an earlier start; remove the directory and start again")
else: json.dump(dict(chunk=chunk, spots=total), open(_lay, "w"))
index = A.panel + (".ont.mmi" if KIND == "long" else ".sr.mmi")
assert os.path.exists(index), index
pieces = [(i, a + 1, min(a + chunk, total)) for i, a in enumerate(range(0, total, chunk))]
log = open(os.path.join(OUT, RUN + ".screen.log"), "a")
def say(*x):
    msg = time.strftime("%H:%M:%S ") + " ".join(str(v) for v in x); print(msg, flush=True); log.write(msg + "\n"); log.flush()

def is_done(i): return os.path.exists(os.path.join(PARTS, f"{i}.json"))

todo = [p for p in pieces if not is_done(p[0])]
say(f"{RUN}: {total} spots in {len(pieces)} pieces of {chunk}; {len(todo)} to do; {'remote, ' + str(streams) + ' streams' if remote else 'stored file'}")
for f in glob.glob(os.path.join(PARTS, "*.tmp")) + glob.glob(os.path.join(PARTS, "*.fa.raw")): os.remove(f)

# ---------------------------------------------------------------- conversion (producer threads)
MAXREADY = 2 if not remote else max(4, streams)        # converted pieces waiting on disk
ready = []; lock = threading.Condition(); it = iter(todo); failed = []
def convert():
    while True:
        with lock:
            while len(ready) >= MAXREADY and not failed: lock.wait(1)
            p = next(it, None)
        if p is None or failed: return
        i, a, b = p; raw = os.path.join(PARTS, f"{i}.fa.raw")
        cmd = [BIN + "/fastq-dump", "--fasta", "0", "--stdout", "--defline-seq", ">$ac.$si" + (".$ri" if KIND == "short" else ""), "-N", str(a), "-X", str(b)]
        if KIND == "short": cmd.insert(3, "--split-spot")
        for attempt in range(4):
            with open(raw + ".tmp", "w") as fo: r = subprocess.run(cmd + [src], stdout=fo, stderr=subprocess.PIPE, text=True)
            if r.returncode == 0 and os.path.getsize(raw + ".tmp") > 0: break
            time.sleep(5 * (attempt + 1))
        else:
            with lock: failed.append((i, r.stderr[-300:])); lock.notify_all()
            return
        os.rename(raw + ".tmp", raw)
        with lock: ready.append(p); lock.notify_all()
threads = [threading.Thread(target=convert, daemon=True) for _ in range(streams)]
for t in threads: t.start()

# ---------------------------------------------------------------- mapping and choice (this thread)
ORG = {"cp": "cp", "asmcp": "cp", "algmt": "algmt", "asmalgmt": "algmt", "slugmt": "slugmt"}
AWK = r'''BEGIN{ while ((getline l < names) > 0) k[l] = 1 }
NR % 2 == 1 { h = substr($0, 2); n++; split(h, f, "."); keep = (h in k); smp = (S > 0 && f[2] % S == 0); next }
{ b += length($0); if (keep) { print ">" h > hits; print $0 > hits } if (smp) { print ">" h > sample; print $0 > sample } }
END{ print n "\t" b > cnt }'''
def process(batch):
    files = [os.path.join(PARTS, f"{i}.fa.raw") for i, _, _ in batch]; tag = batch[0][0]
    # short reads: mapped in batches of 100 Mb, which keeps minimap2 under 1 GB of memory (the result does not depend on -K)
    preset = ["-x", "map-ont"] if KIND == "long" else ["-x", "sr", "-K", "100M"]
    mmlog = open(os.path.join(PARTS, "minimap2.log"), "w")
    if KIND == "long":
        cat = None
        mm = subprocess.Popen([BIN + "/minimap2"] + preset + ["--secondary=no", "-t", "2", index] + files, stdout=subprocess.PIPE, stderr=mmlog, text=True, bufsize=1 << 20)
    else:
        # one stream of reads: given two files, the short-read preset would take them for the two ends of read pairs
        cat = subprocess.Popen(["cat"] + files, stdout=subprocess.PIPE)
        mm = subprocess.Popen([BIN + "/minimap2"] + preset + ["--secondary=no", "-t", "2", index, "-"], stdin=cat.stdout, stdout=subprocess.PIPE, stderr=mmlog, text=True, bufsize=1 << 20)
        cat.stdout.close()
    names_path = os.path.join(PARTS, f"{tag}.names.tmp"); paf_path = os.path.join(PARTS, f"{tag}.paf.tmp")
    stats = collections.Counter()
    if KIND == "long":
        keepn = set()
        with open(paf_path, "w") as fp:
            for line in mm.stdout:
                c = line.rstrip("\n").split("\t")
                fp.write("\t".join(c[:12] + [x for x in c[12:] if x.startswith("dv:f:")]) + "\n")
                keepn.add(c[0])
        stats["hit_reads"] = len(keepn)
    else:
        best = {}; lines = collections.defaultdict(list)
        for line in mm.stdout:
            c = line.rstrip("\n").split("\t"); r = c[0]; al = int(c[9])
            cls = ORG.get(c[5].split("|")[0], "other")
            full = (int(c[3]) - int(c[2])) >= 0.8 * int(c[1])
            if r not in best or al > best[r][0]: best[r] = (al, cls, full)
            lines[r.rsplit(".", 1)[0]].append("\t".join(c[:12]))
        spots = collections.defaultdict(dict)
        for r, v in best.items():
            s, m = r.rsplit(".", 1); spots[s][m] = v
        keepn = set()
        with open(paf_path, "w") as fp:
            for s, d in spots.items():
                m1, m2 = d.get("1"), d.get("2")
                plain = m1 and m2 and m1[2] and m2[2] and m1[1] == m2[1] and m1[1] != "other"
                if plain:
                    stats["plain_" + m1[1]] += 1
                    if int(s.rsplit(".", 1)[1]) % A.plain_keep: continue
                    stats["plain_kept_" + m1[1]] += 1
                else:
                    stats["other_pairs_with_a_match"] += 1
                    kinds = {v[1] for v in (m1, m2) if v}
                    if len(kinds) == 1 and "other" not in kinds:          # one kind of organelle genome and nothing else
                        k_ = next(iter(kinds)); stats["partial_" + k_] += 1
                        if int(s.rsplit(".", 1)[1]) % A.match_keep: continue
                        stats["partial_kept_" + k_] += 1
                keepn.add(s + ".1"); keepn.add(s + ".2"); fp.write("\n".join(lines[s]) + "\n")
        stats["pairs_with_a_match"] = len(spots)
    mm.wait(); mmlog.close()
    if cat is not None and cat.wait() != 0: raise RuntimeError("could not read the pieces")
    if mm.returncode != 0: raise RuntimeError(f"minimap2 failed with code {mm.returncode}: " + open(os.path.join(PARTS, "minimap2.log")).read()[-400:])
    open(names_path, "w").write("\n".join(keepn) + ("\n" if keepn else ""))
    hits = os.path.join(PARTS, f"{tag}.hits.tmp"); samp = os.path.join(PARTS, f"{tag}.sample.tmp"); cnt = os.path.join(PARTS, f"{tag}.cnt.tmp")
    for f in (hits, samp): open(f, "w").close()
    subprocess.run(["awk", "-v", "names=" + names_path, "-v", "hits=" + hits, "-v", "sample=" + samp, "-v", "cnt=" + cnt, "-v", f"S={A.sample}", AWK] + files, check=True)
    n, b = (int(x) for x in open(cnt).read().split())
    for f, ext in ((hits, "hits.fa"), (samp, "sample.fa"), (paf_path, "paf")): os.rename(f, os.path.join(PARTS, f"{tag}.{ext}"))
    for f in (names_path, cnt): os.remove(f)
    info = dict(pieces=[i for i, _, _ in batch], reads=n, bases=b, kept=len(keepn), stats=dict(stats))
    for i, _, _ in batch:
        json.dump(info if i == tag else dict(with_piece=tag), open(os.path.join(PARTS, f"{i}.json.tmp"), "w")); os.rename(os.path.join(PARTS, f"{i}.json.tmp"), os.path.join(PARTS, f"{i}.json"))
    for f in files: os.remove(f)
    return info

t0 = time.time(); nb = 0; got = 0
while True:
    with lock:
        while not ready and not failed and any(t.is_alive() for t in threads): lock.wait(1)
        if failed: say("conversion failed:", failed[0]); sys.exit(1)
        batch = ready[:4] if remote else ready[:1]; del ready[:len(batch)]; lock.notify_all()
    if not batch:
        if not any(t.is_alive() for t in threads) and not ready: break
        continue
    info = process(batch); nb += info["bases"]; got += len(batch)
    if got % 5 < len(batch) or got == len(todo): say(f"{got}/{len(todo)} pieces, {nb / 1e9:.2f} Gb this session, {nb / 1e6 / max(1, time.time() - t0):.1f} Mb/s; kept {info['kept']} reads in the last batch")
if any(not is_done(p[0]) for p in pieces): say("not all pieces finished"); sys.exit(1)

# ---------------------------------------------------------------- put the pieces together
infos = [json.load(open(os.path.join(PARTS, f"{i}.json"))) for i, _, _ in pieces]; infos = [x for x in infos if "reads" in x]
reads = sum(x["reads"] for x in infos); bases = sum(x["bases"] for x in infos); st = collections.Counter()
for x in infos: st.update(x["stats"])
def cat(ext, dest):
    with open(dest, "wb") as fo:
        for i, _, _ in pieces:
            f = os.path.join(PARTS, f"{i}.{ext}")
            if os.path.exists(f):
                with open(f, "rb") as fi: shutil.copyfileobj(fi, fo)
pre = os.path.join(OUT, RUN)
open(pre + ".count", "w").write(f"{reads}\t{bases}\n")
cat("paf", pre + ".p1.paf")
if A.sample: cat("sample.fa", pre + ".sample.fa")
if KIND == "long":
    subprocess.run([sys.executable, os.path.join(HERE, "31_select_hits.py"), pre + ".p1.paf", pre + ".count", pre, "long"], check=True)
    # the chosen reads are pulled straight out of the pieces: a library rich in ribosomal DNA sets aside
    # several gigabases of reads, and they are not copied a second time
    files = [f for f in (os.path.join(PARTS, f"{i}.hits.fa") for i, _, _ in pieces) if os.path.exists(f)]
    with open(pre + ".hits.fa.tmp", "w") as fo:
        for k in range(0, len(files), 200): subprocess.run([BIN + "/seqkit", "grep", "-j", "2", "-w", "0", "-f", pre + ".keep.names"] + files[k:k + 200], stdout=fo, check=True)
    os.rename(pre + ".hits.fa.tmp", pre + ".hits.fa")
    sel = json.load(open(pre + ".select.json")); sel.update(method="pieces", sample_one_in=A.sample); json.dump(sel, open(pre + ".select.json", "w"), indent=1)
else:
    cat("hits.fa", pre + ".hits.fa")
    json.dump(dict(reads=reads, bases=bases, method="pieces", plain_keep_one_in=A.plain_keep, match_keep_one_in=A.match_keep, sample_one_in=A.sample, **{k: int(v) for k, v in st.items()}), open(pre + ".select.json", "w"), indent=1)
shutil.rmtree(PARTS)
open(done_flag, "w").write("ok\n")
say(f"{RUN}: {reads} reads, {bases / 1e9:.2f} Gb screened; {dict(st)}")
