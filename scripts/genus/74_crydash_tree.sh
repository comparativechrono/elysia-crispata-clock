#!/usr/bin/env bash
# Where does Cry-DASH2, the second Cry-DASH gene that only sacoglossans have, come from: a duplicate of the
# animal Cry-DASH gene, or a gene taken from an alga?
#
# The E. crispata Cry-DASH1 and Cry-DASH2 proteins (results/genus/alignments) are compared with DIAMOND
# (--more-sensitive, E <= 1e-10) with 13 animal proteomes that are not sea slugs and the 13 algal, plant and
# cyanobacterial proteomes of the algal database; the best match of each proteome to either protein is
# kept. They are joined by the Cry-DASH proteins of Swiss-Prot (CRYD_*: plants, algae, fungi, vertebrates,
# bacteria, archaea), the Cry-DASH1 and Cry-DASH2 proteins of the sea slugs and of Aplysia and Biomphalaria,
# and an outgroup of photolyases (class I and II CPD photolyases, 6-4 photolyases and the animal
# cryptochromes related to them). Alignment: MAFFT L-INS-i; trimming: trimAl -automated1; tree: IQ-TREE
# (LG+F+G4, 1000 ultrafast bootstraps, seed 1).
#
# usage: 74_crydash_tree.sh WORKDIR OUTDIR
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd)
W=$1; OUT=$2
DATA=${DATA:-/home/claude/data}
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
mkdir -p "$W" "$OUT"
MET="APLCA BIOGL LOTGI POMCA MAGGI MIZYE CAPTE NEMVE BRAFL HUMAN DROME TRICA DANPL"
ALG="BRYSP OSTQU CHLRE VOLCA CHLVA OSTTA MICPC ARATH ECTSI PHATC CYAME CHOCR SYNY3"

python3 - "$REPO/results/genus/alignments" "$W" <<'EOF'
import sys
aln, w = sys.argv[1:3]
def read(p):
    name, seq = None, []
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield name, "".join(seq).replace("-", "")
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line)
    if name: yield name, "".join(seq).replace("-", "")
with open(f"{w}/slug.fa", "w") as fo, open(f"{w}/q.fa", "w") as fq:
    for gene in ("CryDASH1", "CryDASH2", "Phr6-4", "PhrCPD"):
        for n, s in read(f"{aln}/{gene}.aln.fa"):
            sp = n.split("|")[0] if "|" in n else n.split("_")[0]
            fo.write(f">SLUG|{sp}|{gene}\n{s}\n")
            if gene.startswith("CryDASH") and sp == "ELYCR":
                fq.write(f">SLUG|ELYCR|{gene}\n{s}\n")
EOF

if [ ! -s "$W/q_vs_prot.tsv" ]; then
  : > "$W/proteomes.fa"
  for p in $MET; do zcat "$DATA/refprot/$p.fa.gz" | awk -v t=$p '/^>/{split($1,a,"|"); id=(a[2]!=""?a[2]:substr($1,2)); print ">"t"|"id; next}{print}' >> "$W/proteomes.fa"; done
  for p in $ALG; do zcat "$DATA/algae/$p.fa.gz" | awk -v t=$p '/^>/{split($1,a,"|"); id=(a[2]!=""?a[2]:substr($1,2)); print ">"t"|"id; next}{print}' >> "$W/proteomes.fa"; done
  diamond makedb --in "$W/proteomes.fa" -d "$W/proteomes" --quiet
  diamond blastp -q "$W/q.fa" -d "$W/proteomes.dmnd" --more-sensitive -e 1e-10 -k 200 -p 2 --quiet \
    -o "$W/q_vs_prot.tsv" -f 6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore
fi

python3 - "$W" "$DATA/algae/uniprot_sprot.fasta.gz" <<'EOF'
import collections, gzip, sys
w, sprot = sys.argv[1:3]
def read(p):
    op = gzip.open if p.endswith(".gz") else open
    name, seq = None, []
    with op(p, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name: yield name, "".join(seq)
                name, seq = line[1:], []
            else:
                seq.append(line)
    if name: yield name, "".join(seq)
best = {}
for line in open(f"{w}/q_vs_prot.tsv"):
    q, s, pid, ln, ql, sl, qs, qe, ev, bs = line.split("\t")
    p = s.split("|")[0]
    if p not in best or float(bs) > best[p][1]:
        best[p] = (s, float(bs))
want = {s for s, _ in best.values()}
seqs = {}
for h, s in read(f"{w}/proteomes.fa"):
    if h.split()[0] in want:
        seqs[h.split()[0]] = s
# Swiss-Prot: every Cry-DASH, and an outgroup of photolyases and cryptochromes
OUTGROUP = {"PHR_ECOLI", "PHR_SYNY3", "PHR_YEAST", "PHR_ARATH", "PHR_POTTR", "UVR3_ARATH", "PHR_DROME", "CRY1_HUMAN", "CRY_DROME"}
acc_seen = {s.split("|")[1] for s in want}
for h, s in read(sprot):
    entry, acc = h.split()[0].split("|")[2], h.split()[0].split("|")[1]
    if (entry.startswith("CRYD_") or entry in OUTGROUP) and acc not in acc_seen:
        seqs[f"SPROT|{entry}|{acc}"] = s
        acc_seen.add(acc)
for h, s in read(f"{w}/slug.fa"):
    sp, gene = h.split("|")[1:3]
    if gene.startswith("CryDASH") or sp == "ELYCR":
        seqs[h] = s
written = set()
with open(f"{w}/tree_input.fa", "w") as fo:
    for k, s in sorted(seqs.items(), key=lambda kv: not kv[0].startswith("SLUG")):
        s = s.replace("*", "")
        if s in written:                                       # the same protein from two sources
            continue
        written.add(s)
        fo.write(f">{k}\n{s}\n")
print(len(written), "sequences")
EOF

mafft --localpair --maxiterate 1000 --thread "${THREADS:-2}" --quiet "$W/tree_input.fa" > "$OUT/crydash.aln.fa"
trimal -in "$OUT/crydash.aln.fa" -out "$OUT/crydash.trimmed.fa" -automated1
iqtree -s "$OUT/crydash.trimmed.fa" -m LG+F+G4 -B 1000 -T "${THREADS:-2}" --seed 1 --prefix "$OUT/crydash" -redo -quiet
cp "$W/tree_input.fa" "$OUT/crydash.input.fa"

# support of the groups that matter, read from the unrooted consensus tree (a group and its complement are
# the same split), and the tree rooted on the photolyases and cryptochromes
python3 - "$OUT" <<'EOF'
import json, sys
from Bio import Phylo
out = sys.argv[1]
t = Phylo.read(f"{out}/crydash.contree", "newick")
leaves = frozenset(c.name for c in t.get_terminals())
splits = {}
for c in t.get_nonterminals():
    s = frozenset(x.name for x in c.get_terminals())
    if 1 < len(s) < len(leaves) - 1 and c.confidence is not None:
        splits[s] = splits[leaves - s] = c.confidence
def support(group):
    g = frozenset(group)
    return splits.get(g)
has = lambda *keys: [n for n in leaves if any(k in n for k in keys)]
OUTGROUP = has("PHR_", "UVR3", "Phr6-4", "PhrCPD", "HUMAN|", "BRAFL|", "CAPTE|", "DANPL|", "TRICA|", "DROME|", "MIZYE|", "BRYSP|", "OSTQU|")
CD2 = has("CryDASH2")
HETERO = has("CryDASH")                                  # sea slugs, Aplysia and Biomphalaria
ANIMAL = HETERO + has("NEMVE|", "MAGGI|", "LOTGI|", "POMCA|", "CRYD_DANRE", "CRYD_XENLA")
res = dict(sequences=len(leaves), sites=None,
           groups={"Cry-DASH2 of the sea slugs": dict(members=len(CD2), support=support(CD2)),
                   "Cry-DASH2 with Aplysia Cry-DASH1": dict(support=support(CD2 + has("APLCA|CryDASH1"))),
                   "Cry-DASH1 of the sea slugs with Biomphalaria Cry-DASH1": dict(support=support([n for n in HETERO if "CryDASH1" in n and "APLCA" not in n])),
                   "heterobranch Cry-DASH (sea slugs, Aplysia, Biomphalaria; Cry-DASH1 and Cry-DASH2)": dict(members=len(HETERO), support=support(HETERO)),
                   "animal Cry-DASH (cnidarian, molluscs, vertebrates)": dict(members=len(ANIMAL), support=support(ANIMAL)),
                   "outgroup (photolyases and cryptochromes)": dict(members=len(OUTGROUP), support=support(OUTGROUP))})
for line in open(f"{out}/crydash.iqtree"):
    if line.startswith("Input data:"):
        res["sites"] = line.split("with")[1].split()[0]
t.root_with_outgroup(t.common_ancestor(OUTGROUP))
Phylo.write(t, f"{out}/crydash.rooted.nwk", "newick")
json.dump(res, open(f"{out}/summary.json", "w"), indent=1)
print(json.dumps(res, indent=1))
EOF
rm -f "$OUT"/crydash.ckp.gz "$OUT"/crydash.mldist "$OUT"/crydash.bionj "$OUT"/crydash.splits.nex "$OUT"/crydash.uniqueseq.phy   # intermediate files of IQ-TREE
