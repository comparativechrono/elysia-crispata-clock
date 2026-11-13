#!/usr/bin/env bash
# Re-align the representative proteins to all four assemblies, annotate Pfam domains
# and rebuild the catalogue.   usage: 08_finalise_models.sh
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
REPO=/home/claude/elysia-crispata-clock; WD=/home/claude/work/clock; DATA=/home/claude/data
cd $WD/cat
sed 's/ .*//' representatives.fa > reps.fa
[ -s $DATA/genome/sanger.mpi ] || miniprot -t2 -d $DATA/genome/sanger.mpi $DATA/genome/sanger.fa 2>/dev/null
for g in sanger:$DATA/genome/sanger.mpi alt:$DATA/alt/alt_genomic.fna.gz purdue:$DATA/other_asm/purdue.fa tampa:$DATA/other_asm/tampa_genomic.fna.gz; do
  n=${g%%:*}; p=${g#*:}
  miniprot -t2 --gff --trans --outs=0.5 -N 8 $p reps.fa > rep_$n.gff 2> rep_$n.log
done
if ! cmp -s reps.fa reps.pfam.done 2>/dev/null; then
  zcat $DATA/pfam/Pfam-A.hmm.gz | hmmsearch --cut_ga --cpu 2 --domtblout reps.pfam.domtbl -o /dev/null - reps.fa
  cp reps.fa reps.pfam.done
fi
cd $WD && python3 $REPO/scripts/07_build_catalogue.py $REPO/config cat $DATA $REPO/results
