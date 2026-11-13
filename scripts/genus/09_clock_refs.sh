#!/usr/bin/env bash
# Reference files for the genus-wide clock-gene search.
#   reps.faa        the 47 E. crispata catalogue proteins, named by gene symbol
#   ecr_refs.dmnd   every E. crispata protein of both annotations plus the 47 (as REP|symbol)
#   rep_equiv.tsv   which E. crispata proteins are the catalogue genes themselves (>= 95 % identical
#                   to a catalogue protein over at least half their length, or a catalogue anchor)
#   hsdm.dmnd       human and Drosophila proteomes, to check genes that E. crispata lacks
#   nf_seeds.tsv    seed label, UniProt accession and gene group for those genes
# usage: 09_clock_refs.sh queries.fa ecr_all.fa allprot.fa seeds.fa config/clock_catalogue.tsv refs_dir
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
Q=$1; ECR=$2; ALL=$3; SEEDS=$4; CAT=$5; R=$6; mkdir -p "$R"
seqkit grep -v -r -p '^NF_' "$Q" -w 0 > "$R/reps.faa" 2> /dev/null
{ cat "$ECR"; sed 's/^>/>REP|/' "$R/reps.faa"; } > "$R/ecr_refs.fa"
diamond makedb --in "$R/ecr_refs.fa" -d "$R/ecr_refs" --quiet
diamond makedb --in "$ECR" -d "$R/ecr_all" --quiet
diamond blastp -q "$R/reps.faa" -d "$R/ecr_all" -o "$R/reps_vs_ecr.tsv" --sensitive -e 1e-20 -k 50 --threads 2 --quiet --outfmt 6 qseqid sseqid pident length qlen slen bitscore
python3 - "$R" "$CAT" "$ECR" "$SEEDS" <<'PY'
import sys, collections
R, cat, ecr, seeds = sys.argv[1:5]
best = {}
for l in open(f"{R}/reps_vs_ecr.tsv"):
    q, s, pid, ln, ql, sl, bs = l.rstrip("\n").split("\t")
    if float(pid) >= 95 and int(ln) >= 0.5 * int(sl) and (s not in best or float(bs) > best[s][1]): best[s] = (q, float(bs))
anchor = {}
for l in open(cat):
    if l.startswith("#") or l.startswith("symbol"): continue
    f = l.rstrip("\n").split("\t")
    for a in f[5].split(","):
        if a.strip() and a.strip() != "-": anchor[a.strip()] = f[0]
n_anchor = 0
for l in open(ecr):
    if l[0] != ">": continue
    h = l[1:].split()[0]; p = h.split("|")
    for x in p[1:]:
        if x in anchor and h not in best: best[h] = (anchor[x], 0.0); n_anchor += 1
with open(f"{R}/rep_equiv.tsv", "w") as fo:
    for s, (q, bs) in sorted(best.items()): fo.write(f"{s}\t{q}\n")
c = collections.Counter(q for q, b in best.values())
print(len(best), "E. crispata proteins equated with a catalogue gene;", n_anchor, "by anchor only; genes without any:", [g for g in set(anchor.values()) if g not in c])
grp = {"HsBHLHE40": "Dec1/2", "HsBHLHE41": "Dec1/2", "HsFBXL3": "Fbxl3/21", "MmFbxl21": "Fbxl3/21", "DmTyf": "Tyf", "HsCIPC": "Cipc", "HsCIART": "Ciart",
       "HsAANAT": "Aanat", "DmAANAT1": "Aanat", "HsASMT": "Asmt"}
with open(f"{R}/nf_seeds.tsv", "w") as fo:
    for l in open(seeds):
        if l[0] == ">":
            lab, acc = l[1:].split("|")[:2]
            if lab in grp: fo.write(f"{lab}\t{acc}\t{grp[lab]}\n")
    fo.write("HsFBXL21P\tQ9UKT6\tFbxl3/21\n")
PY
seqkit grep -r -p '^(HUMAN|DROME)\|' "$ALL" -w 0 > "$R/hsdm.fa" 2> /dev/null
diamond makedb --in "$R/hsdm.fa" -d "$R/hsdm" --quiet
rm -f "$R/hsdm.fa"; echo "refs ready: $(grep -c '>' "$R/reps.faa") catalogue proteins, $(grep -c '>' "$R/ecr_refs.fa") reference proteins"
