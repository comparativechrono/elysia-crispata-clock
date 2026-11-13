#!/usr/bin/env bash
# Second, more sensitive DNA screen of the assemblies for plastid-like sequence.
# The window screen (scripts 22, 23) is fast but approximate and is swamped by microsatellites.
# Here one plastid genome per genus of siphonous green algae in the panel, both kleptoplast
# genomes of the xgElyCris1 animal, and (optionally) extra sequences such as the Dasycladales
# plastid genes are used as BLAST queries against each assembly (dc-megablast, low-complexity
# filter on), which finds matches down to roughly 75 % identity and reports exact identities.
# usage: 24_plastid_blast_screen.sh outdir panel.fa panel.tsv extra.fa|- name=fasta[.gz] ...
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
OUT=$1; PANEL=$2; PTSV=$3; EXTRA=$4; shift 4
mkdir -p "$OUT"; cd "$OUT"
if [ ! -s queries.fa ]; then
  python3 - "$PTSV" > query.names <<'PY'
import sys
seen = set()
for l in list(open(sys.argv[1]))[1:]:
    f = l.rstrip("\n").split("\t")
    if f[1] != "cp": continue
    d = f[-1]; genus = d.split()[0]
    if d.startswith("UNVERIFIED") or "partial" in d: continue
    key = " ".join(d.split()[:2]) if genus == "Penicillus" else genus
    if key in seen: continue
    seen.add(key); print(f[0])
PY
  seqkit grep -w 0 -f query.names "$PANEL" > queries.fa
  [ "$EXTRA" != "-" ] && seqkit seq -w 0 "$EXTRA" >> queries.fa
fi
echo "$(grep -c '>' queries.fa) query sequences, $(seqkit stats -T queries.fa | tail -1 | cut -f5) bp"
for a in "$@"; do
  n=${a%%=*}; fa=${a#*=}
  [ -s "$n.blast.tsv" ] && continue
  if [[ "$fa" == *.gz ]]; then zcat "$fa" | makeblastdb -in - -dbtype nucl -title "$n" -out "db_$n" > /dev/null; else makeblastdb -in "$fa" -dbtype nucl -out "db_$n" > /dev/null; fi
  blastn -task dc-megablast -query queries.fa -db "db_$n" -evalue 1e-10 -dust yes -num_threads 2 -max_target_seqs 5000 \
         -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" > "$n.blast.tsv.part"
  mv "$n.blast.tsv.part" "$n.blast.tsv"; rm -f db_$n.n*
  echo "$n: $(wc -l < "$n.blast.tsv") hits"
done
