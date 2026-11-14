#!/usr/bin/env bash
# Everything that is worked out from the two passes of one library (31c_screen_run.py, 33_pass2.sh).
#   short reads   37_pair_analysis.py   composition, donors, every pair and mate that joins organelle
#                                       DNA to slug DNA, the same counts in the random sample
#                 55_end_piece_control.py   whether a piece is recognised as well at the end of a mate as at its start
#                 38_plastid_contigs.sh plastid contigs and their barcode genes
#   long reads    34_read_analysis.py   composition, donors, joins
#                 46_candidates.py      the shape of every read that joins organelle DNA to slug DNA
#                 51_recurrent_junctions.py   joins shared by two molecules, found by comparing reads
#                 46_candidates.py      the same shapes on the random sample of the library (sample.paf)
#                 36_barcodes.py        barcode genes read out of the reads (accurate reads only)
# usage: 49_analyse_library.sh WORKDIR GENOME.fa.gz PANEL_DIR SCREEN.json hifi|clr|short BARCODES.fa
#        NO_CONTIGS=1 skips the assembly of plastid contigs
#        DEPTH_PAF=file  long reads: alignments of all reads (47_library_control.sh), to count the reads
#                        that cross the place of each insertion-shaped read without an insert
#        P1_INDEX=file   short reads: the index the first pass mapped the library to (PANEL_PREFIX.sr.mmi)
#        EXTRA_ORGANELLE=SEQUENCE:START-END=cp   short reads: a stretch of the assembly to count as organelle DNA (37_pair_analysis.py)
#        REPEAT_FA=file  long reads: a repeat of the slug genome; joins into it are counted, and how
#                        close they come to each other on the organelle genome (51_recurrent_junctions.py)
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
HERE=$(cd "$(dirname "$0")" && pwd)
W=$1; G=$2; PANEL=$3; SCREEN=$4; MODE=$5; BC=$6
[ -s "$W/pass2.paf" ] || { echo "$W: no second pass yet" >&2; exit 1; }
GBP=$(awk '{s+=$2} END{print s}' "$G.fai")
if [ "$MODE" = short ]; then
  # bases of the assembly that are not N (kept beside the assembly)
  [ -s "$G.acgt" ] || zcat "$G" | grep -v '>' | tr -d 'Nn\n' | wc -c > "$G.acgt"
  python3 "$HERE/37_pair_analysis.py" "$W" "$PANEL/panel.tsv" "$W/pairs.json" "$SCREEN" --genome-bp "$GBP" --genome-acgt "$(cat "$G.acgt")" --genome "$G" --junctions "$W/junctions.tsv.gz" ${EXTRA_ORGANELLE:+--organelle-stretch "$EXTRA_ORGANELLE"}
  # do the two ends of a mate pass the identity test equally often? (for the count of which piece a join begins with)
  # (P1_INDEX: the index the first pass mapped this library to, for the second control of that script)
  if [ -s "$W/sample.paf" ]; then python3 "$HERE/55_end_piece_control.py" "$W" "$G" "$PANEL/panel.fa" "$W/end_piece_control.json"; fi
  [ "${NO_CONTIGS:-0}" = 1 ] || flock /tmp/bigmem.lock bash "$HERE/38_plastid_contigs.sh" "$W" "$PANEL/panel.fa" "$BC"
else
  python3 "$HERE/34_read_analysis.py" "$W" "$PANEL/panel.tsv" "$MODE" "$W/analysis.json" "$SCREEN"
  python3 "$HERE/46_candidates.py" "$W" "$MODE" "$W/candidates.json" "$SCREEN" --genome-bp "$GBP" --panel "$PANEL/panel.fa" ${DEPTH_PAF:+--depth-paf "$DEPTH_PAF"}
  python3 "$HERE/51_recurrent_junctions.py" "$W/candidates.json" "$W/all_hits.fa" "$W/recurrent.json" "$MODE" "$W/recurrent.tmp" ${REPEAT_FA:+--repeat "$REPEAT_FA" "$PANEL/panel.fa.fai"} && rm -rf "$W/recurrent.tmp"
  if [ -s "$W/sample.paf" ]; then python3 "$HERE/46_candidates.py" "$W" "$MODE" "$W/background.json" "$SCREEN" --paf "$W/sample.paf" --no-p1 --genome-bp "$GBP"; fi
  if [ "$MODE" = hifi ]; then python3 "$HERE/36_barcodes.py" reads "$W/all_hits.fa" "$BC" "$W/barcodes"; fi
fi
