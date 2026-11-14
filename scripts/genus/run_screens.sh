#!/usr/bin/env bash
# The chloroplast and gene-transfer screens of the six species added to the site, in the order
# they were run. E. crispata has its own, earlier pipeline (scripts/10 to 29); the two steps of
# this one that were applied to it as well are in crispata() below.
#
# usage: run_screens.sh STEP ...      STEP = panel | assembly | reads | controls | loci | proteins | crispata | collect | all
#        DRY=1 run_screens.sh all     prints the commands without running them
#
# Layout (all outside the repository; set with the variables below)
#   $DATA/asm/<slug>/genome.fa.gz           assemblies (00_fetch_assemblies.sh; bgzip + faidx)
#   $DATA/genus_refs/                       extra references for the panel and the barcode records
#   $WORK/panel, asm_screen, reads/<slug>[.<library>], protscreen, chl_control
# Needs 6 GB of memory and about 15 GB of free disk. Every sequencing run is read in full, in
# pieces; a run of up to 11 GB is fetched first, a larger one is read from the archive directly.
# Nothing but the chosen reads and a random sample is kept. One run takes 10 to 90 minutes on two
# cores; the whole of "reads" takes most of a day.
#
# The published run differs from the commands below in three ways that do not change the result:
#   - the libraries of E. atroviridis and the long reads of E. chlorotica had been read in full
#     with an earlier, single-stream version of the first pass (same mapping, same choice of reads;
#     scripts/genus/31_read_screen.sh in the history of this repository, commit 9188e43);
#     for E. atroviridis the sample was therefore drawn afterwards (31d_sample_blocks.py)
#   - the first pass of E. timida was interrupted once and went on from the pieces it had finished
#   - the second pass of the long reads of E. chlorotica ran before the places of the candidate
#     insertions were added to the baits, so those places were examined for the long reads with
#     the alignment of all reads (47_library_control.sh) and for the short reads with the baits
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
  # The places to test with reads join the panels, so that the first pass keeps their reads.
  # E. chlorotica: the organelle-like stretches of the assembly, the 86-base plastid-like element of
  # RQTK01000372.1 (found by the pair analysis) and the places of the insertion-shaped long reads
  # (found by 46_candidates.py on the long reads; the short-read screens came after)
  local C=$WORK/reads/elysia-chlorotica L=$WORK/reads/elysia-leucolegnote
  run mkdir -p "$C" "$L"
  run python3 "$HERE/40_locus_baits.py" "$WORK/asm_screen/elysia-chlorotica.screen.json" "$(genome elysia-chlorotica)" "$C/baits.assembly.fa"
  run python3 "$HERE/40_locus_baits.py" "$WORK/asm_screen/elysia-chlorotica.screen.json" "$(genome elysia-chlorotica)" "$C/baits.fa" 3000 \
      RQTK01000372.1:352485-352571:plastid RQTK01001142.1:52936-52937:candidate RQTK01000223.1:439998-439999:candidate RQTK01000407.1:282521-282522:candidate \
      RQTK01000317.1:349275-349276:candidate RQTK01000080.1:821436-821437:candidate RQTK01000046.1:111205-111206:candidate RQTK01000070.1:476498-476499:candidate \
      RQTK01000662.1:99526-99527:candidate RQTK01000339.1:318401-318402:candidate RQTK01000321.1:187607-187608:candidate_mt RQTK01000299.1:10800-10801:candidate_mt
  run python3 "$HERE/40_locus_baits.py" "$WORK/asm_screen/elysia-leucolegnote.screen.json" "$(genome elysia-leucolegnote)" "$L/baits.fa" 3000 JBQKML010000381.1:0-863:plastid
  # panel.chl: the whole panel with the six stretches of the E. chlorotica assembly (long reads);
  # panel.org, panel.org.chl, panel.org.leu: organelle genomes only, with the places to test (short reads)
  run bash "$HERE/30b_mask_and_index.sh" "$PANEL" "$C/baits.assembly.fa" chl
  run bash "$HERE/30b_mask_and_index.sh" "$PANEL" - org
  run bash "$HERE/30b_mask_and_index.sh" "$PANEL" "$C/baits.fa" org.chl
  run bash "$HERE/30b_mask_and_index.sh" "$PANEL" "$L/baits.fa" org.leu
}

# one library: first pass of every run, second pass, analysis
#   lib DIR PANEL_NAME "RUN ..." long|short hifi|clr|short [options of 31c_screen_run.py]
lib(){
  local dir=$1 pn=$2 runs=$3 kind=$4 mode=$5; shift 5
  local slug=${dir%%.*} W=$WORK/reads/$dir r spots size src
  run mkdir -p "$W" "$DATA/genus_reads/tmp"
  for r in $runs; do
    [ -e "$W/$r.done" ] && continue
    spots=$(run_info "$r" spots); size=$(run_info "$r" size_MB); src=()
    if [ "$size" -le 11000 ]; then                 # small enough to store while it is read
      run python3 "$HERE/fetch_ranges.py" "https://sra-pub-run-odp.s3.amazonaws.com/sra/$r/$r" "$DATA/genus_reads/tmp/$r.sra" 128 6
      src=(--src "$DATA/genus_reads/tmp/$r.sra")
    fi
    run python3 "$HERE/31c_screen_run.py" "$W" "$PANEL/$pn" "$r" "$kind" "$spots" "${src[@]}" "$@"
    run rm -f "$DATA/genus_reads/tmp/$r.sra"
  done
  if [ "$dir" = elysia-timida ] || [ "$dir" = elysia-chlorotica ]; then      # libraries rich in ribosomal DNA
    for r in $runs; do
      [ -e "$W/$r.allhits.fasta" ] && continue
      run mv "$W/$r.hits.fa" "$W/$r.allhits.fasta"
      run python3 "$HERE/31b_slim_rdna.py" "$W/$r.p1.paf" "$W/$r.allhits.fasta" "$W/$r.hits.fa"
    done
  fi
  big bash "$HERE/33_pass2.sh" "$W" "$(genome "$slug")" "$PANEL/panel.fa" "$WORK/asm_screen/$slug.screen.json" "$mode"
  if [ "$kind" = short ]; then P1_INDEX=$PANEL/$pn.sr.mmi analyse "$dir" "$mode"; else analyse "$dir" "$mode"; fi
}
analyse(){
  local dir=$1 mode=$2 slug=${1%%.*}
  run env ${DEPTH_PAF:+DEPTH_PAF=$DEPTH_PAF} ${REPEAT_FA:+REPEAT_FA=$REPEAT_FA} ${P1_INDEX:+P1_INDEX=$P1_INDEX} ${EXTRA_ORGANELLE:+EXTRA_ORGANELLE=$EXTRA_ORGANELLE} bash "$HERE/49_analyse_library.sh" "$WORK/reads/$dir" "$(genome "$slug")" "$PANEL" "$WORK/asm_screen/$slug.screen.json" "$mode" "$BC"
}
# spots and size of a run: run_info RUN spots|size_MB   (E-utilities; the answers are kept in $WORK/runs.tsv)
run_info(){
  [ -s "$WORK/runs.tsv" ] || printf 'run\tspots\tbases\tsize_MB\n' > "$WORK/runs.tsv"
  if ! grep -q "^$1	" "$WORK/runs.tsv"; then
    curl -sS ${SSL_CERT_FILE:+--cacert "$SSL_CERT_FILE"} "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=$1&rettype=runinfo&retmode=text" \
      | python3 -c 'import csv, sys; [print(r["Run"], r["spots"], r["bases"], r["size_MB"], sep="\t") for r in csv.DictReader(sys.stdin)]' >> "$WORK/runs.tsv"
  fi
  awk -F'\t' -v r="$1" -v c="$2" 'NR==1{for(i=1;i<=NF;i++) h[$i]=i} $1==r{print $h[c]}' "$WORK/runs.tsv"
}

reads(){
  # the panels for short reads hold the organelle genomes only (30b_mask_and_index.sh org), plus the places to test
  lib elysia-papillosa           panel.dust    "ERR14842698 ERR14842697" long hifi --sample 100
  lib elysia-papillosa.papi1     panel.dust    "ERR14842696" long hifi --sample 200
  lib elysia-atroviridis         panel.dust    "SRR34916404" long hifi
  lib elysia-timida              panel.dust    "SRR30016879" long hifi --sample 200
  lib elysia-chlorotica          panel.chl     "SRR7647856"  long clr
  lib elysia-chlorotica.adult    panel.org.chl "SRR7647851 SRR7647852 SRR7647853 SRR7647854 SRR7647849 SRR7647850" short short --sample 800 --plain-keep 50 --match-keep 50
  lib elysia-chlorotica.eggs     panel.org.chl "SRR768897 SRR768899" short short --sample 400
  # BMAT01005342.1 of the E. marginata assembly is plastid DNA from a part of the Bryopsis plastid genome that no reference
  # holds (read pairs tie both its ends to the plastid genome, and it is covered as deeply): its mates count as plastid DNA
  EXTRA_ORGANELLE=BMAT01005342.1:0-3705=cp lib elysia-marginata panel.org "DRR238951 DRR238952" short short --sample 800
  lib elysia-leucolegnote        panel.org.leu "SRR41369319" short short --sample 400
  lib elysia-leucolegnote.hainan panel.org.leu "SRR17036819" short short --sample 400
}

controls(){        # what a library makes by itself, where the first pass kept no sample of it
  # E. atroviridis: 40 evenly spaced blocks of the run stand in for the sample
  local A=$WORK/reads/elysia-atroviridis
  run python3 "$HERE/31d_sample_blocks.py" "$A" SRR34916404 "$(run_info SRR34916404 spots)" long 40 1000
  big env SAMPLE_ONLY=1 bash "$HERE/33_pass2.sh" "$A" "$(genome elysia-atroviridis)" "$PANEL/panel.fa" "$WORK/asm_screen/elysia-atroviridis.screen.json" hifi
  analyse elysia-atroviridis hifi
  # E. chlorotica, long reads: the run is small enough to align every read
  local C=$WORK/chl_control
  run bash "$HERE/47_library_control.sh" "$C" SRR7647856 "$(genome elysia-chlorotica)" "$PANEL/panel.fa" "$WORK/asm_screen/elysia-chlorotica.screen.json" clr "cp|NC_011600.1 slugmt|NC_010567.1"
  run cp "$C/allreads.json" "$C/control.json" "$WORK/reads/elysia-chlorotica/"
  # the 358-base satellite of the slug to which the plastid stretch at the start of RQTK01001031.1 is joined
  [ "${DRY:-0}" = 1 ] || samtools faidx "$(genome elysia-chlorotica)" RQTK01001031.1:1651-3400 > "$WORK/reads/elysia-chlorotica/satellite.fa"
  DEPTH_PAF=$C/control.paf REPEAT_FA=$WORK/reads/elysia-chlorotica/satellite.fa analyse elysia-chlorotica clr
}

loci(){            # reads at each place where an assembly holds organelle-like DNA, and at each candidate insertion
  # The first pass kept the reads of these places (they are baits in the panels panel.chl, panel.org.chl
  # and panel.org.leu), so each test runs on the reads already chosen.
  local C=$WORK/reads/elysia-chlorotica G lib
  G=$(genome elysia-chlorotica)
  # places the assembly holds (the first seven of the bait table) and places of the insertion-shaped long reads (the rest)
  if [ "${DRY:-0}" != 1 ]; then
    awk -F'\t' 'NR==1 || $5=="plastid" || $5=="numt"' "$C/baits.fa.tsv" > "$C/baits.assembly.tsv"
    awk -F'\t' 'NR==1 || $5 ~ /^candidate/' "$C/baits.fa.tsv" > "$C/baits.candidates.tsv"
  fi
  run python3 "$HERE/41_locus_test.py" "$C/baits.assembly.tsv" "$G" "$C/SRR7647856.hits.fa" clr "$C/locus_tests.json"
  for lib in adult eggs; do
    run python3 "$HERE/41_locus_test.py" "$C/baits.assembly.tsv" "$G" "$C.$lib/all_hits.fa" short "$C.$lib/locus_tests.json"
    run python3 "$HERE/41_locus_test.py" "$C/baits.candidates.tsv" "$G" "$C.$lib/all_hits.fa" short "$C.$lib/candidate_places.json"
    # the joins between chloroplast DNA and the satellite repeat that come close to each other in two long-read
    # molecules: the one in the assembly (plastid position 20,865) with its neighbour (20,899), and 78,759 with 78,799
    run python3 "$HERE/53_join_in_short_reads.py" "$C.$lib" "$G" "$C.$lib/satellite_joins.json" RQTK01001031.1:1604-3400 20865 20899 78759 78799
  done
  # E. leucolegnote: the plastid stretch at the start of JBQKML010000381.1
  local L=$WORK/reads/elysia-leucolegnote
  G=$(genome elysia-leucolegnote)
  run python3 "$HERE/41_locus_test.py" "$L/baits.fa.tsv" "$G" "$L/all_hits.fa" short "$L/locus_tests.json"
  run python3 "$HERE/41_locus_test.py" "$L/baits.fa.tsv" "$G" "$L.hainan/all_hits.fa" short "$L.hainan/locus_tests.json"
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

crispata(){        # two steps of this pipeline applied to the E. crispata library as well
  local X=$WORK/reads/elysia-crispata G r spots; G=$(genome elysia-crispata)
  # the stretches of plastid DNA in the assembly, to be hidden in the second pass
  run bash "$REPO/scripts/22_assembly_organelle_screen.sh" "$PANEL/panel.ont.mmi" "$WORK/asm_screen" "elysia-crispata=$G"
  run python3 "$HERE/35_assembly_screen.py" "$PANEL/panel.tsv" "$PANEL/panel.fa" "$G.fai" "$WORK/asm_screen/elysia-crispata.paf" "$WORK/asm_screen/elysia-crispata.screen.json"
  # reads with a short, near-exact organelle piece away from both ends (from the first pass of scripts/10_bait_screen.sh)
  run python3 "$HERE/54_internal_pieces.py" "$HGT/pass1" "$X"
  run python3 "$HERE/52_fetch_reads.py" "$X/targeted.names" "$X/targeted.hits.fa"
  # a sample of one read in a hundred of each run: 40 evenly spaced blocks
  for r in ERR11843472 ERR12353011 ERR12356115; do
    spots=$(run_info "$r" spots)
    run python3 "$HERE/31d_sample_blocks.py" "$X" "$r" "$spots" long 40 $((spots / 4000))
    [ "${DRY:-0}" = 1 ] || printf '%s\t%s\n' "$spots" "$(run_info "$r" bases)" > "$X/$r.count"
  done
  big env SAMPLE_MAX=1000000 bash "$HERE/33_pass2.sh" "$X" "$G" "$PANEL/panel.fa" "$WORK/asm_screen/elysia-crispata.screen.json" hifi
  local gbp; gbp=$(awk '{s+=$2} END{print s}' "$G.fai")
  run python3 "$HERE/46_candidates.py" "$X" hifi "$X/candidates.json" "$WORK/asm_screen/elysia-crispata.screen.json" --genome-bp "$gbp" --panel "$PANEL/panel.fa"
  run python3 "$HERE/46_candidates.py" "$X" hifi "$X/background.json" "$WORK/asm_screen/elysia-crispata.screen.json" --paf "$X/sample.paf" --no-p1 --genome-bp "$gbp"
  # every read with a join between chloroplast DNA and slug DNA (the junction table of scripts/14_chimera_analysis.py),
  # fetched again, aligned like the reads of the other species and compared with every other read around the join:
  # is any join carried by two molecules, also where the slug DNA is a repeat or was placed on another assembly?
  local Y=$WORK/reads/elysia-crispata.joins
  run mkdir -p "$Y"
  [ "${DRY:-0}" = 1 ] || zcat "$REPO/results/intermediate/chimeric_read_junctions.tsv.gz" | awk -F'\t' 'NR > 1 && $2 == "cp+slug" {print $1}' | sort -u > "$Y/joins.names"
  run python3 "$HERE/52_fetch_reads.py" "$Y/joins.names" "$Y/joins.hits.fa" 12
  big bash "$HERE/33_pass2.sh" "$Y" "$G" "$PANEL/panel.fa" "$WORK/asm_screen/elysia-crispata.screen.json" hifi
  run python3 "$HERE/46_candidates.py" "$Y" hifi "$Y/candidates.json" "$WORK/asm_screen/elysia-crispata.screen.json" --no-p1 --genome-bp "$gbp"
  run python3 "$HERE/51_recurrent_junctions.py" "$Y/candidates.json" "$Y/all_hits.fa" "$Y/recurrent.json" hifi
}

collect(){
  run python3 "$HERE/39_collect_species.py" "$REPO" "$WORK/reads" "$WORK/asm_screen" "$WORK/protscreen"
  run python3 "$HERE/45_export_tables.py" "$REPO" "$WORK"
  run python3 "$REPO/scripts/30_build_site.py" "$REPO"          # the E. crispata pages read its part of the collected results
  run python3 "$HERE/21_build_species_pages.py" "$REPO"
  run python3 "$HERE/20_build_genus_site.py" "$REPO"
}

[ $# -gt 0 ] || { sed -n 2,8p "$0"; exit 1; }
for step in "$@"; do
  case $step in
    all) panel; assembly; reads; controls; loci; proteins; crispata; collect ;;
    panel|assembly|reads|controls|loci|proteins|crispata|collect) $step ;;
    *) echo "unknown step: $step" >&2; exit 1 ;;
  esac
done
