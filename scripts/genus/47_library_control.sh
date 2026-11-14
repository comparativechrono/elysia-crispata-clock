#!/usr/bin/env bash
# What a long-read library makes by itself: the whole library aligned, not only the reads that the
# first pass chose. Used for E. chlorotica, whose long-read run is small (9.4 Gb).
#
#   1. every read of MINLEN bases or more is aligned to the assembly plus the organelle genomes
#      named in ORGANELLES (sequences of the panel), chain by chain, without base-level alignment
#   2. 46_candidates.py on those alignments: how often a read joins two unrelated pieces of slug DNA,
#      folds back on itself, or joins organelle DNA to slug DNA, all counted in the same way
#      -> allreads.json
#   3. the reads that could hold an insertion shape (47_shape_reads.py) are aligned base by base
#      like the chosen reads (33_pass2.sh) and their shapes counted (46_candidates.py): insertion
#      shapes with plastid DNA, with mitochondrial DNA and with slug DNA from elsewhere, found by
#      identical steps -> control.json
#
# usage: 47_library_control.sh OUTDIR RUN GENOME.fa.gz PANEL.fa SCREEN.json hifi|clr "ORGANELLES" [MINLEN=2000]
#   e.g. 47_library_control.sh work/chl_control SRR7647856 genome.fa.gz panel.fa screen.json clr "cp|NC_011600.1 slugmt|NC_010567.1"
# Needs about 3 GB of disk for the reads (kept compressed) and 3 GB of memory.
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=$1; RUN=$2; G=$(realpath "$3"); PANEL=$(realpath "$4"); SCREEN=$(realpath "$5"); KIND=$6; ORG=$7; MINLEN=${8:-2000}
mkdir -p "$OUT/p2"; cd "$OUT"
PRESET=map-pb; [ "$KIND" = hifi ] && PRESET=map-hifi
GBP=$(awk '{s+=$2} END{print s}' "$G.fai")
if [ ! -s all_reads.fa.zst ]; then
  python3 "$HERE/fetch_ranges.py" "https://sra-pub-run-odp.s3.amazonaws.com/sra/$RUN/$RUN" "$RUN.sra" 128 6
  fastq-dump --fasta 0 --stdout --defline-seq '>$ac.$si' "$RUN.sra" 2> dump.err | zstd -3 -q -o all_reads.fa.zst.tmp && mv all_reads.fa.zst.tmp all_reads.fa.zst
  rm -f "$RUN.sra"
fi
zstd -dc all_reads.fa.zst | awk 'NR%2==0{n++; b+=length($0)} END{print n"\t"b}' > "$RUN.count"
if [ ! -s control.paf ]; then
  (zcat "$G" | sed 's/ .*//'; for o in $ORG; do samtools faidx "$PANEL" "$o"; done) > ref.fa
  zstd -dc all_reads.fa.zst | seqkit seq -m "$MINLEN" -w 0 | minimap2 -x $PRESET --secondary=no -t 2 -K 300M ref.fa - 2> control.mm2.log | cut -f1-12 > control.paf.tmp && mv control.paf.tmp control.paf
  rm -f ref.fa
fi
python3 "$HERE/46_candidates.py" . "$KIND" allreads.json "$SCREEN" --paf control.paf --no-identity --no-p1 --genome-bp "$GBP" --one-in 1
python3 "$HERE/47_shape_reads.py" control.paf shape.names
if [ ! -s p2/pass2.paf ]; then
  zstd -dc all_reads.fa.zst | seqkit grep -j 2 -w 0 -f shape.names > p2/shape.hits.fa.tmp && mv p2/shape.hits.fa.tmp p2/shape.hits.fa
  flock /tmp/bigmem.lock bash "$HERE/33_pass2.sh" p2 "$G" "$PANEL" "$SCREEN" "$KIND"
fi
python3 "$HERE/46_candidates.py" p2 "$KIND" control.json "$SCREEN" --no-p1 --depth-paf control.paf
