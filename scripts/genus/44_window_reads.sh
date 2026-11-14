#!/usr/bin/env bash
# All reads of a sequencing run that match a small window of an assembly, whether or not they
# match the reference panel. Used to test one place in an assembly with every read of a library
# (41_locus_test.py), where the panel-selected reads of 31_read_screen.sh would be a biased sample.
#   pass 1  the run is converted and every read mapped to the windows (minimap2, no base-level
#           alignment); reads with a match are noted (for short reads: both mates of the pair)
#   pass 2  the run is converted again and the noted reads are written out
# The run is fetched from the SRA open-data bucket and deleted afterwards (KEEP_SRA=1 keeps it).
# usage: 44_window_reads.sh RUN long|short MAXBASES WINDOWS.fa OUT.fa [DOWNLOAD_DIR]
set -uo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
HERE=$(cd "$(dirname "$0")" && pwd)
RUN=$1; KIND=$2; MAX=${3:-0}; WIN=$4; OUT=$5; DL=${6:-$(dirname "$OUT")/dl}
mkdir -p "$DL" "$(dirname "$OUT")"
SRA="$DL/$RUN.sra"
if [ ! -s "$SRA" ] || [ -e "$SRA.done" ]; then
  python3 "$HERE/fetch_ranges.py" "https://sra-pub-run-odp.s3.amazonaws.com/sra/$RUN/$RUN" "$SRA" 128 6 || exit 1
fi
LIMIT='NR%2==1{h=$0; next} {n++; b+=length($0); print h; print $0; if (max>0 && b>=max) exit} END{print n"\t"b > cnt}'
if [ "$KIND" = long ]; then
  dump(){ fasterq-dump --fasta-unsorted --stdout --threads 2 --seq-defline '>$ac.$si' "$SRA" 2>/dev/null; }
  dump | awk -v max="$MAX" -v cnt="$OUT.count" "$LIMIT" | minimap2 -x map-ont -k15 -w10 --secondary=no -t 2 "$WIN" - 2>/dev/null | cut -f1 | sort -u > "$OUT.names"
else
  dump(){ fastq-dump --fasta 0 --split-spot --stdout --defline-seq '>$ac.$si.$ri' "$SRA" 2>/dev/null; }
  dump | awk -v max="$MAX" -v cnt="$OUT.count" "$LIMIT" | minimap2 -x sr --secondary=no -t 2 "$WIN" - 2>/dev/null | cut -f1 \
    | sed 's/\.[12]$//' | sort -u | awk '{print $1".1"; print $1".2"}' > "$OUT.names"
fi
dump | awk -v max="$MAX" -v cnt=/dev/null "$LIMIT" | seqkit grep -j 2 -w 0 -f "$OUT.names" > "$OUT"
echo "$RUN: $(cut -f1 "$OUT.count") reads, $(cut -f2 "$OUT.count") bases examined; $(grep -c '>' "$OUT") reads match the windows"
[ "${KEEP_SRA:-0}" = 1 ] || rm -f "$SRA"
