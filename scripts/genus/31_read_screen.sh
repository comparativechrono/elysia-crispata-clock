#!/usr/bin/env bash
# Pass 1 of the read screen for one sequencing run.
# The run is fetched from the SRA open-data bucket, converted, and every read mapped to the
# reference panel (30_build_panel.py) without base-level alignment. Reads with a hit are then
# pulled out by converting the run a second time. Only those reads are kept; the run itself is
# deleted afterwards (set KEEP_SRA=1 to keep it).
#   long reads   fasterq-dump | minimap2 (k15, w10)        stops after MAXBASES bases
#   short reads  fastq-dump writes the two mates of a pair one after the other (RUN.n.1, RUN.n.2);
#                minimap2 -x sr maps each mate; if either mate hits, both are kept
#                                                            stops after MAXBASES bases
# usage: 31_read_screen.sh OUTDIR PANEL_PREFIX RUN long|short MAXBASES [DOWNLOAD_DIR]
#        REMOTE=1      read the run from the archive as a stream instead of storing it
#        KEEP_SRA=1    do not delete the stored run        PASS1_ONLY=1   stop after the first pass
set -uo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=$1; PANEL=$2; RUN=$3; KIND=$4; MAX=${5:-0}; DL=${6:-$OUT/dl}
mkdir -p "$OUT" "$DL"
[ -s "$OUT/$RUN.done" ] && { echo "$RUN already screened"; exit 0; }
SRA="$DL/$RUN.sra"
if [ "${REMOTE:-0}" = 1 ]; then
  # read the run straight from the bucket (for runs too large to store); the screened part is
  # kept compressed until the chosen reads have been pulled out
  SRA="https://sra-pub-run-odp.s3.amazonaws.com/sra/$RUN/$RUN"
elif [ ! -s "$SRA" ] || [ -e "$SRA.done" ]; then
  python3 "$HERE/fetch_ranges.py" "https://sra-pub-run-odp.s3.amazonaws.com/sra/$RUN/$RUN" "$SRA" 128 6 || exit 1
fi
TMP="$OUT/$RUN.sample.fa.zst"
keep(){ if [ "${REMOTE:-0}" = 1 ]; then tee >(zstd -1 -q -o "$TMP"); else cat; fi; }
LIMIT='NR%2==1{h=$0; next} {n++; b+=length($0); print h; print $0; if (max>0 && b>=max) exit} END{print n"\t"b > cnt}'
if [ "$KIND" = long ]; then
  [ -s "$PANEL.ont.mmi" ] || minimap2 -x map-ont -k15 -w10 -d "$PANEL.ont.mmi" "$PANEL.fa" 2>/dev/null
  fasterq-dump --fasta-unsorted --stdout --threads 2 --seq-defline '>$ac.$si' "$SRA" 2>"$OUT/$RUN.dump.err" \
    | awk -v max="$MAX" -v cnt="$OUT/$RUN.count" "$LIMIT" | keep \
    | minimap2 -x map-ont --secondary=no -t 2 "$PANEL.ont.mmi" - 2>"$OUT/$RUN.mm2.log" > "$OUT/$RUN.p1.paf"
else
  [ -s "$PANEL.sr.mmi" ] || minimap2 -x sr -d "$PANEL.sr.mmi" "$PANEL.fa" 2>/dev/null
  fastq-dump --fasta 0 --split-spot --stdout --defline-seq '>$ac.$si.$ri' "$SRA" 2>"$OUT/$RUN.dump.err" \
    | awk -v max="$MAX" -v cnt="$OUT/$RUN.count" "$LIMIT" | keep \
    | minimap2 -x sr --secondary=no -t 2 "$PANEL.sr.mmi" - 2>"$OUT/$RUN.mm2.log" > "$OUT/$RUN.p1.paf"
fi
wait
[ -s "$OUT/$RUN.count" ] || { echo "$RUN: no reads counted"; exit 1; }
# PASS1_ONLY=1 stops here and keeps the run, for libraries whose reads are chosen in another way
# (the second adult of E. chlorotica in run_screens.sh, where plastid DNA is 14% of the library)
[ "${PASS1_ONLY:-0}" = 1 ] && { echo "$RUN: first pass only; run kept at $SRA"; exit 0; }
# reads to keep: every read with a hit, except that reads whose only hits fall in "hot" panel
# windows (short local matches shared by very many reads: slug repeats that resemble a stretch of
# a panel sequence) are thinned to at most 60 per window
python3 "$HERE/31_select_hits.py" "$OUT/$RUN.p1.paf" "$OUT/$RUN.count" "$OUT/$RUN" "$KIND"
# second conversion of the run to pull the chosen reads out (conversion is fast; nothing is stored in between)
if [ "${REMOTE:-0}" = 1 ]; then
  sleep 2; zstd -dc "$TMP" | seqkit grep -j 2 -w 0 -f "$OUT/$RUN.keep.names" > "$OUT/$RUN.hits.fa"; rm -f "$TMP"
elif [ "$KIND" = long ]; then
  fasterq-dump --fasta-unsorted --stdout --threads 2 --seq-defline '>$ac.$si' "$SRA" 2>>"$OUT/$RUN.dump.err" | seqkit grep -j 2 -w 0 -f "$OUT/$RUN.keep.names" > "$OUT/$RUN.hits.fa"
else
  fastq-dump --fasta 0 --split-spot --stdout --defline-seq '>$ac.$si.$ri' "$SRA" 2>>"$OUT/$RUN.dump.err" | awk -v max="$MAX" -v cnt=/dev/null "$LIMIT" | seqkit grep -j 2 -w 0 -f "$OUT/$RUN.keep.names" > "$OUT/$RUN.hits.fa"
fi
echo "$RUN: $(cut -f1 "$OUT/$RUN.count") reads, $(cut -f2 "$OUT/$RUN.count") bases screened; $(grep -c '>' "$OUT/$RUN.hits.fa") reads kept"
[ "${KEEP_SRA:-0}" = 1 ] || [ "${REMOTE:-0}" = 1 ] || rm -f "$SRA"
echo ok > "$OUT/$RUN.done"
