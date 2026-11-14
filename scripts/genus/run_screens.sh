#!/usr/bin/env bash
# The chloroplast and gene-transfer screens of the six species added to the site, in the order
# they were run. E. crispata has its own, earlier pipeline (scripts/10 to 29).
#
# usage: run_screens.sh STEP ...      STEP = panel | assembly | reads | loci | proteins | collect | all
#        DRY=1 run_screens.sh all     prints the commands without running them
#
# Layout (all outside the repository; set with the variables below)
#   $DATA/asm/<slug>/genome.fa.gz           assemblies (00_fetch_assemblies.sh; bgzip + faidx)
#   $DATA/genus_refs/                       extra references for the panel and the barcode records
#   $WORK/panel, asm_screen, reads/<slug>[.<library>], protscreen
# Needs about 8 GB of memory and 20 GB of free disk. Sequencing runs are fetched one at a time
# and deleted once their reads have been screened.
#
# Three steps of the published run were finished by hand after a failure, with the same result
# as the commands below:
#   - E. atroviridis: the reads chosen by the first pass were pulled out in a second fetch of the run
#   - E. marginata: the run was read from the archive as a stream and the stream broke 44 reads
#     short of the 4.5 Gb limit; the screen uses the 30,000,000 reads before the break
#   - E. chlorotica, second adult: the first pass keeps too many pairs for the second (14% of this
#     library is plastid DNA), so the pairs that touch a tested place were all kept and every tenth
#     of the other pairs; the plastid share is taken from the first pass over all 6 Gb
#     (see thin_adult below)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd)
WORK=${WORK:-/home/claude/work/genus}; DATA=${DATA:-/home/claude/data}
HGT=${HGT:-/home/claude/work/hgt}          # E. crispata work directory: its panel and the protein database
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
SPECIES="elysia-papillosa elysia-timida elysia-atroviridis elysia-chlorotica elysia-marginata elysia-leucolegnote"
PANEL=$WORK/panel; REF=$DATA/genus_refs; BC=$REF/barcodes/barcodes.fa
run(){ if [ "${DRY:-0}" = 1 ]; then printf '%s\n' "$*"; else echo "[$(date +%T)] $*" >&2; "$@"; fi; }
big(){ run flock /tmp/bigmem.lock "$@"; }      # one large-memory job at a time
genome(){ echo "$DATA/asm/$1/genome.fa.gz"; }

panel(){
  run python3 "$HERE/30_build_panel.py" "$HGT/panel.masked.fa" "$HGT/panel.tsv" "$DATA/dasyclad/dasyclad_cp.fa" "$REF/cp_extra/cp_extra.fa" \
      "$REF/rdna/extra_rdna.fa" "$REF/mito/elysia_mito.fa" "$PANEL"
  run bash "$HERE/30b_mask_and_index.sh" "$PANEL"
}

assembly(){        # organelle-like DNA in each assembly, as DNA
  local args=(); for s in $SPECIES; do args+=("$s=$(genome "$s")"); done
  run bash "$REPO/scripts/22_assembly_organelle_screen.sh" "$PANEL/panel.ont.mmi" "$WORK/asm_screen" "${args[@]}"
  for s in $SPECIES; do
    run python3 "$HERE/35_assembly_screen.py" "$PANEL/panel.tsv" "$PANEL/panel.fa" "$(genome "$s").fai" "$WORK/asm_screen/$s.paf" "$WORK/asm_screen/$s.screen.json"
  done
  # E. chlorotica: the places to test with reads join the panel, so that the first pass keeps their reads
  run python3 "$HERE/40_locus_baits.py" "$WORK/asm_screen/elysia-chlorotica.screen.json" "$(genome elysia-chlorotica)" "$WORK/reads/elysia-chlorotica/baits.fa"
  run bash "$HERE/30b_mask_and_index.sh" "$PANEL" "$WORK/reads/elysia-chlorotica/baits.fa" chl
}

# one library: first pass, second pass, analysis.   lib DIR PANEL_NAME RUN long|short MAXBASES hifi|clr|short
lib(){
  local dir=$1 pn=$2 r=$3 kind=$4 max=$5 mode=$6 slug=${1%%.*}; local W=$WORK/reads/$dir
  run mkdir -p "$W" "$DATA/genus_reads/$slug"
  case $dir in
    elysia-marginata)                    # a 37.6 Gb run: read from the archive as a stream, not stored
      run env REMOTE=1 bash "$HERE/31_read_screen.sh" "$W" "$PANEL/$pn" "$r" "$kind" "$max" ;;
    elysia-chlorotica.adult)             # 14% plastid DNA: choose the reads here (thin_adult)
      run env PASS1_ONLY=1 bash "$HERE/31_read_screen.sh" "$W" "$PANEL/$pn" "$r" "$kind" "$max" "$DATA/genus_reads/$slug"
      run thin_adult "$W" "$r" "$max" "$DATA/genus_reads/$slug/$r.sra" ;;
    *) run bash "$HERE/31_read_screen.sh" "$W" "$PANEL/$pn" "$r" "$kind" "$max" "$DATA/genus_reads/$slug" ;;
  esac
  if [ "$dir" = elysia-timida ] || [ "$dir" = elysia-chlorotica ]; then      # libraries rich in ribosomal DNA
    run mv "$W/$r.hits.fa" "$W/$r.allhits.fasta"
    run python3 "$HERE/31b_slim_rdna.py" "$W/$r.p1.paf" "$W/$r.allhits.fasta" "$W/$r.hits.fa"
  fi
  big bash "$HERE/33_pass2.sh" "$W" "$(genome "$slug")" "$PANEL/panel.fa" "$WORK/asm_screen/$slug.screen.json" "$mode"
  if [ "$mode" = short ]; then
    run python3 "$HERE/37_pair_analysis.py" "$W" "$PANEL/panel.tsv" "$W/pairs.json" "$WORK/asm_screen/$slug.screen.json"
    run bash "$HERE/38_plastid_contigs.sh" "$W" "$PANEL/panel.fa" "$BC"
  else
    run python3 "$HERE/34_read_analysis.py" "$W" "$PANEL/panel.tsv" "$mode" "$W/analysis.json" "$WORK/asm_screen/$slug.screen.json"
    [ "$mode" = hifi ] && run python3 "$HERE/36_barcodes.py" reads "$W/all_hits.fa" "$BC" "$W/barcodes"
  fi
  return 0
}
# E. chlorotica, second adult: keep every pair that touches a tested place and every tenth other pair
thin_adult(){
  local W=$1 r=$2 max=$3 sra=$4
  awk -F'\t' '{n=$1; sub(/\.[12]$/,"",n); split(n,a,"."); if ($6 ~ /^bait\|/) print n; else if (a[2] % 10 == 0) print n}' "$W/$r.p1.paf" | sort -u > "$W/$r.keep.spots"
  awk '{print $1".1"; print $1".2"}' "$W/$r.keep.spots" > "$W/$r.keep.names"
  fastq-dump --fasta 0 --split-spot --stdout --defline-seq '>$ac.$si.$ri' "$sra" 2>/dev/null \
    | awk -v max="$max" 'NR%2==1{h=$0; next} {b+=length($0); print h; print $0; if (max>0 && b>=max) exit}' \
    | seqkit grep -j 2 -w 0 -f "$W/$r.keep.names" > "$W/$r.hits.fa"
  rm -f "$sra"
  # plastid share of the whole sample from the first pass: mates matching a plastid reference at >= 97% over >= 80 bases
  awk -F'\t' '$6 ~ /^cp\|/ && $11>=80 && $10/$11>=0.97 {if (!($1 in s)) {s[$1]=1; n++; b+=$4-$3}} END{print n"\t"b}' "$W/$r.p1.paf" > "$W/$r.cp_pass1.tsv"
  # the detailed analysis sees a tenth of the library; the full count is kept beside it
  cp "$W/$r.count" "$W/$r.count.full"; awk '{printf "%d\t%d\n", $1 / 10, $2 / 10}' "$W/$r.count.full" > "$W/$r.count"
  echo ok > "$W/$r.done"
}

reads(){
  lib elysia-papillosa           panel.dust ERR14842698 long  0           hifi
  lib elysia-atroviridis         panel.dust SRR34916404 long  0           hifi
  lib elysia-timida              panel.dust SRR30016879 long  20000000000 hifi
  lib elysia-chlorotica          panel.chl  SRR7647856  long  0           clr
  lib elysia-chlorotica.adult    panel.chl  SRR7647851  short 6000000000  short
  lib elysia-chlorotica.eggs     panel.chl  SRR768897   short 0           short
  lib elysia-marginata           panel.dust DRR238951   short 4500000000  short
  lib elysia-leucolegnote        panel.dust SRR41369319 short 6000000000  short
  lib elysia-leucolegnote.hainan panel.dust SRR17036819 short 5000000000  short
}

loci(){            # reads at each place where an assembly holds organelle-like DNA inside another sequence
  local C=$WORK/reads/elysia-chlorotica G
  G=$(genome elysia-chlorotica)
  run python3 "$HERE/41_locus_test.py" "$C/baits.fa.tsv" "$G" "$C/SRR7647856.hits.fa" clr "$C/locus_tests.json"
  run python3 "$HERE/41_locus_test.py" "$C/baits.fa.tsv" "$G" "$C.adult/SRR7647851.hits.fa" short "$C.adult/locus_tests.json"
  run python3 "$HERE/41_locus_test.py" "$C/baits.fa.tsv" "$G" "$C.eggs/SRR768897.hits.fa" short "$C.eggs/locus_tests.json"
  # the 86-base plastid-like element at RQTK01000372.1:352,485-352,571: every egg read against the place
  run mkdir -p "$C.eggs/locus372"
  if [ "${DRY:-0}" = 1 ]; then echo "samtools faidx $G RQTK01000372.1:350001-355000 > $C.eggs/locus372/window.fa"
  else samtools faidx "$G" RQTK01000372.1:350001-355000 | sed '1s/.*/>w372/' | tr 'acgtn' 'ACGTN' > "$C.eggs/locus372/window.fa"; fi
  run bash "$HERE/44_window_reads.sh" SRR768897 short 0 "$C.eggs/locus372/window.fa" "$C.eggs/locus372/reads.fa" "$DATA/genus_reads/elysia-chlorotica"
  if [ "${DRY:-0}" != 1 ]; then
    printf 'bait\tcontig\twindow_start\twindow_end\tkind\tlocus_start\tlocus_end\nbait372\tRQTK01000372.1\t350000\t355000\tplastid\t352485\t352571\n' > "$C.eggs/locus372/bait.tsv"
  fi
  run python3 "$HERE/41_locus_test.py" "$C.eggs/locus372/bait.tsv" "$G" "$C.eggs/locus372/reads.fa" short "$C.eggs/locus372/locus_test.json"
  # E. leucolegnote: the plastid stretch at the start of JBQKML010000381.1, with every read of both libraries
  local L=$WORK/reads/elysia-leucolegnote
  G=$(genome elysia-leucolegnote)
  run python3 "$HERE/40_locus_baits.py" "$WORK/asm_screen/elysia-leucolegnote.screen.json" "$G" "$L/baits.fa" 3000 JBQKML010000381.1:0-863:plastid
  run bash "$HERE/44_window_reads.sh" SRR41369319 short 6000000000 "$L/baits.fa" "$L/locus/SRR41369319.window.fa" "$DATA/genus_reads/elysia-leucolegnote"
  run bash "$HERE/44_window_reads.sh" SRR17036819 short 5000000000 "$L/baits.fa" "$L.hainan/locus/SRR17036819.window.fa" "$DATA/genus_reads/elysia-leucolegnote"
  run python3 "$HERE/41_locus_test.py" "$L/baits.fa.tsv" "$G" "$L/locus/SRR41369319.window.fa" short "$L/locus_tests.json"
  run python3 "$HERE/41_locus_test.py" "$L/baits.fa.tsv" "$G" "$L.hainan/locus/SRR17036819.window.fa" short "$L.hainan/locus_tests.json"
}

proteins(){        # algal genes searched as proteins, and what else lies on the sequences they lead to
  local P=$WORK/protscreen Q=$REPO/results/intermediate/protein_screen_queries.fa DB=$HGT/alien/alien_db.dmnd
  local order="elysia-chlorotica elysia-leucolegnote elysia-marginata elysia-atroviridis elysia-timida elysia-papillosa" args=() gargs=() s
  run mkdir -p "$P"
  run python3 "$HERE/42_protein_screen_inputs.py" "$WORK/asm_screen" "$P" $SPECIES
  for s in $order; do
    args+=("$s=x"); gargs+=("$s=$(genome "$s")")
    if [ "${DRY:-0}" = 1 ]; then echo "miniprot -t2 --gff --trans -p 0.1 --outs=0.1 -N 200 $(genome "$s") $Q > $P/$s.gff"
    elif [ ! -s "$P/$s.done" ]; then flock /tmp/bigmem.lock miniprot -t2 --gff --trans -p 0.1 --outs=0.1 -N 200 "$(genome "$s")" "$Q" > "$P/$s.gff" 2> "$P/$s.log" && touch "$P/$s.done"; fi
  done
  run python3 "$REPO/scripts/28_protein_screen.py" "$P" "$DB" "$P/embedded_plastid.json" "$Q" "${args[@]}"
  run python3 "$HERE/43_foreign_sequences.py" "$P" "$DB" "$P/foreign_sequences.json" "${gargs[@]}"
}

collect(){
  run python3 "$HERE/39_collect_species.py" "$REPO" "$WORK/reads" "$WORK/asm_screen" "$WORK/protscreen"
  run python3 "$HERE/45_export_tables.py" "$REPO" "$WORK"
  run python3 "$HERE/21_build_species_pages.py" "$REPO"
  run python3 "$HERE/20_build_genus_site.py" "$REPO"
}

[ $# -gt 0 ] || { sed -n 2,8p "$0"; exit 1; }
for step in "$@"; do
  case $step in
    all) panel; assembly; reads; loci; proteins; collect ;;
    panel|assembly|reads|loci|proteins|collect) $step ;;
    *) echo "unknown step: $step" >&2; exit 1 ;;
  esac
done
