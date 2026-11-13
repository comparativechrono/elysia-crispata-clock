#!/usr/bin/env bash
# Name the organism behind foreign contigs by comparing them with the plastid sequences that
# GenBank holds for one candidate group (here Dasycladales, whose plastid genomes are huge,
# gene-poor and not available as complete genomes, so they are missing from the read-screening panel).
# usage: 29_name_foreign_contigs.sh outdir contigs.fa "Dasycladales[Organism] AND (chloroplast[filter] OR plastid[filter])" out.tsv
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
OUT=$1; FA=$2; TERM=$3; TSV=$4
mkdir -p "$OUT"; cd "$OUT"
curl -sS --max-time 60 -G "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi" --data-urlencode "db=nuccore" --data-urlencode "term=$TERM" \
     --data-urlencode "retmax=5000" --data-urlencode "retmode=json" --data-urlencode "usehistory=y" > es.json
WE=$(python3 -c "import json;print(json.load(open('es.json'))['esearchresult']['webenv'])"); QK=$(python3 -c "import json;print(json.load(open('es.json'))['esearchresult']['querykey'])")
sleep 1
curl -sS --max-time 600 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&WebEnv=$WE&query_key=$QK&rettype=fasta&retmode=text&retmax=5000" > refs.fa
makeblastdb -in refs.fa -dbtype nucl -out refs -parse_seqids > /dev/null
blastn -task dc-megablast -query "$FA" -db refs -evalue 1e-20 -max_target_seqs 20 -num_threads 2 \
       -outfmt "6 qseqid sseqid pident length qstart qend sstart send bitscore stitle" | sort -k1,1 -k5,5n -k9,9nr > "$TSV"
echo "$(grep -c '>' refs.fa) reference sequences; $(cut -f1 "$TSV" | sort -u | wc -l) contigs with a match"
# per-locus best match, and a like-for-like comparison of the two best-represented species
python3 - "$TSV" "${TSV%.tsv}.json" <<'PY'
import sys, json, collections
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
loci = collections.OrderedDict(); per = collections.defaultdict(dict)
for r in rows:
    if int(r[3]) < 300: continue
    k = (r[0], int(r[4]) // 400); sp = " ".join(r[9].split()[:2])
    if k not in loci or float(r[8]) > float(loci[k][8]): loci[k] = r
    if sp not in per[k] or float(r[2]) > per[k][sp]: per[k][sp] = float(r[2])
best = [dict(contig=r[0], qstart=int(r[4]), qend=int(r[5]), identity=float(r[2]), length=int(r[3]), subject=r[1].strip("|").split("|")[-1], title=r[9]) for r in loci.values()]
both = [(d["Acetabularia peniculus"], d["Acetabularia acetabulum"]) for d in per.values() if "Acetabularia peniculus" in d and "Acetabularia acetabulum" in d]
json.dump(dict(best=best, same_genes=dict(n=len(both), peniculus=[min(a for a, b in both), max(a for a, b in both)] if both else None,
                                          acetabulum=[min(b for a, b in both), max(b for a, b in both)] if both else None)), open(sys.argv[2], "w"), indent=1)
print(len(best), "loci;", len(both), "with both A. peniculus and A. acetabulum references")
PY
