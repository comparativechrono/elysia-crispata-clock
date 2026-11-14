#!/usr/bin/env python3
"""Fetch named reads of sequencing runs straight from the SRA open-data bucket.

usage: 52_fetch_reads.py NAMES.txt OUT.fa [THREADS=10]
NAMES.txt holds read names of the form RUN.SPOT (one per line). Each spot is read on its own
(fastq-dump -N/-X), at about a dozen reads a minute for each thread: a few thousand reads take half an hour.
"""
import sys, os, subprocess, concurrent.futures as cf
names = [l.strip() for l in open(sys.argv[1]) if l.strip()]
out = sys.argv[2]; threads = int(sys.argv[3]) if len(sys.argv) > 3 else 10
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin")
def get(name):
    run, spot = name.rsplit(".", 1)
    cmd = [BIN + "/fastq-dump", "--fasta", "0", "--stdout", "--defline-seq", ">$ac.$si", "-N", spot, "-X", spot, f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"]
    for attempt in range(4):
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode == 0 and r.stdout.startswith(">"): return r.stdout
    return ""
with cf.ThreadPoolExecutor(threads) as ex: parts = list(ex.map(get, names))
with open(out, "w") as fo:
    for p in parts: fo.write(p)
got = sum(1 for p in parts if p)
print(f"{got} of {len(names)} reads fetched"); sys.exit(0 if got == len(names) else 1)
