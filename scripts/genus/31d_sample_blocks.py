#!/usr/bin/env python3
"""A sample of a sequencing run for a library that was screened before samples were kept:
BLOCKS blocks of BLOCKSIZE consecutive spots, spread evenly over the run, read straight from the
SRA open-data bucket.

usage: 31d_sample_blocks.py OUTDIR RUN SPOTS long|short [BLOCKS=40] [BLOCKSIZE=1000]
Writes OUTDIR/RUN.sample.fa and adds sample_one_in (spots of the run / spots sampled, rounded) to
OUTDIR/RUN.select.json, so that 33_pass2.sh and 46_candidates.py treat it like the sample that
31c_screen_run.py keeps (--sample). The reads of a run are stored in the order in which the
instrument wrote them, so evenly spaced blocks stand in for a random draw.
"""
import sys, os, subprocess, json, concurrent.futures as cf
out, run, spots, kind = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
blocks = int(sys.argv[5]) if len(sys.argv) > 5 else 40
size = int(sys.argv[6]) if len(sys.argv) > 6 else 1000
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin")
src = f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"
starts = [1 + int(i * (spots - size) / max(1, blocks - 1)) for i in range(blocks)]
def get(a):
    cmd = [BIN + "/fastq-dump", "--fasta", "0", "--stdout", "--defline-seq", ">$ac.$si" + (".$ri" if kind == "short" else ""), "-N", str(a), "-X", str(a + size - 1)]
    if kind == "short": cmd.insert(3, "--split-spot")
    for attempt in range(4):
        r = subprocess.run(cmd + [src], capture_output=True, text=True)
        if r.returncode == 0 and r.stdout: return r.stdout
    raise SystemExit(f"block at {a} failed: {r.stderr[-200:]}")
with cf.ThreadPoolExecutor(8) as ex: parts = list(ex.map(get, starts))
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, run + ".sample.fa"), "w") as fo:
    for p in parts: fo.write(p)
n = sum(p.count(">") for p in parts)
sj = os.path.join(out, run + ".select.json"); sel = json.load(open(sj)) if os.path.exists(sj) else {}
sel.update(sample_one_in=round(spots / (blocks * size)), sample_method=f"{blocks} blocks of {size} consecutive spots, evenly spaced"); json.dump(sel, open(sj, "w"), indent=1)
print(f"{run}: {n} reads in {blocks} blocks; one spot in {sel['sample_one_in']}")
