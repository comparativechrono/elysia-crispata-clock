#!/usr/bin/env bash
# The chloroplast and gene-transfer screens of the six species added to the site, in the order
# they were run. E. crispata has its own, earlier pipeline (scripts/10 to 29); the two steps of
# this one that were applied to it as well are in crispata() below.
#
# usage: run_screens.sh STEP ...      STEP = panel | assembly | reads | controls | loci | proteins | crispata | other_assemblies | crispata_animals | regulators | census | collect | all
#        DRY=1 run_screens.sh all     prints the commands without running them
#
# Layout (all outside the repository; set with the variables below)
#   $DATA/asm/<slug>/genome.fa.gz           assemblies (00_fetch_assemblies.sh; bgzip + faidx)
#   $DATA/genus_refs/                       extra references for the panel and the barcode records
#   $WORK/panel, asm_screen, reads/<slug>[.<library>], protscreen, chl_control, census
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
#   ASM=<slug>  the assembly the reads are aligned to, where it is not that of the species (DIR up to the first dot)
#   PANEL=<dir> another panel directory;  PASS2_PRESET=map-ont  Nanopore reads (analysed like the older PacBio reads, mode clr)
lib(){
  local dir=$1 pn=$2 runs=$3 kind=$4 mode=$5; shift 5
  local slug=${ASM:-${dir%%.*}} W=$WORK/reads/$dir r spots size src
  run mkdir -p "$W" "$DATA/genus_reads/tmp"
  for r in $runs; do
    [ -e "$W/$r.done" ] && continue
    spots=$(run_info "$r" spots); size=$(run_info "$r" size_MB); src=()
    if [ "$size" -le "${STORE_MAX_MB:-11000}" ]; then                 # small enough to store while it is read
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
  big env ${PASS2_PRESET:+PASS2_PRESET=$PASS2_PRESET} ${PASS2_INDEX_BASES:+PASS2_INDEX_BASES=$PASS2_INDEX_BASES} ${SAMPLE_EVERY:+SAMPLE_EVERY=$SAMPLE_EVERY} bash "$HERE/33_pass2.sh" "$W" "$(genome "$slug")" "$PANEL/panel.fa" "$WORK/asm_screen/$slug.screen.json" "$mode"
  if [ "$kind" = short ]; then P1_INDEX=$PANEL/$pn.sr.mmi analyse "$dir" "$mode"; else analyse "$dir" "$mode"; fi
}
analyse(){
  local dir=$1 mode=$2 slug=${ASM:-${1%%.*}}
  run env ${DEPTH_PAF:+DEPTH_PAF=$DEPTH_PAF} ${REPEAT_FA:+REPEAT_FA=$REPEAT_FA} ${P1_INDEX:+P1_INDEX=$P1_INDEX} ${EXTRA_ORGANELLE:+EXTRA_ORGANELLE=$EXTRA_ORGANELLE} ${LONG_READ_PLATFORM:+LONG_READ_PLATFORM=$LONG_READ_PLATFORM} bash "$HERE/49_analyse_library.sh" "$WORK/reads/$dir" "$(genome "$slug")" "$PANEL" "$WORK/asm_screen/$slug.screen.json" "$mode" "$BC"
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

other_assemblies(){   # a second assembly of a species (config/species_assemblies.tsv): the alternate haplotype of E. papillosa
  # $DATA/asm/elysia-papillosa-alt/genome.fa.gz holds GCA_977015365.2 (bgzip + faidx). At 1.2 Gb the assembly is more
  # than miniprot can index at once in 6 GB of memory, so it is searched in three parts (10b_miniprot_in_parts.sh).
  local s=elysia-papillosa-alt G P=$WORK/protscreen_alt Q=$REPO/results/intermediate/protein_screen_queries.fa DB=$HGT/alien/alien_db.dmnd
  G=$(genome $s)
  # clock genes: the catalogue of this assembly, built and scored like that of the reference assembly
  run bash "$HERE/10b_miniprot_in_parts.sh" "$G" "$WORK/queries.fa" "$WORK/clock/$s.mp.gff" 3 -G 500000 -p 0.15 --outs 0.15 -N 50
  [ "${DRY:-0}" = 1 ] || { head -n 1 "$REPO/config/species.tsv"; printf '%s\tElysia papillosa\tELYPA\tGCA_977015365.2\txgElyPapi4.2 alternate haplotype\tcontig\tnone\tPacBio HiFi\tTampa, Florida, USA\n' $s; } > "$WORK/clock/species_alt.tsv"
  run python3 "$HERE/11_clock_catalogue.py" $s "$WORK/clock/species_alt.tsv" "$DATA/asm" "$WORK/clock" "$WORK/refs" "$WORK/clock/$s.clock.json"
  run python3 "$HERE/14_other_assemblies.py" "$REPO/config/species_assemblies.tsv" "$WORK/clock" "$REPO/results/genus"
  # organelle-like DNA in the assembly, as DNA
  run bash "$REPO/scripts/22_assembly_organelle_screen.sh" "$PANEL/panel.ont.mmi" "$WORK/asm_screen" "$s=$G"
  run python3 "$HERE/35_assembly_screen.py" "$PANEL/panel.tsv" "$PANEL/panel.fa" "$G.fai" "$WORK/asm_screen/$s.paf" "$WORK/asm_screen/$s.screen.json"
  # algal genes, as proteins (in a directory of its own, beside the screen of the reference assemblies)
  run mkdir -p "$P"
  run python3 "$HERE/42_protein_screen_inputs.py" "$WORK/asm_screen" "$P" $s
  run bash "$HERE/10b_miniprot_in_parts.sh" "$G" "$Q" "$P/$s.gff" 3 -p 0.1 --outs=0.1 -N 200
  run python3 "$REPO/scripts/28_protein_screen.py" "$P" "$DB" "$P/embedded_plastid.json" "$Q" "$s=x"
  run python3 "$HERE/43_foreign_sequences.py" "$P" "$DB" "$P/foreign_sequences.json" "$s=$G"
}

# a library is done once: a second start goes on with the next one. ONLY=<dir> runs that library alone.
once(){
  local d=$WORK/reads/$1
  if [ -n "${ONLY:-}" ] && [ "$ONLY" != "$1" ]; then return 0; fi
  shift; [ -e "$d/lib.done" ] && return 0
  "$@"; [ "${DRY:-0}" = 1 ] || touch "$d/lib.done"
}

crispata_animals(){   # the other sequenced animals of E. crispata, each aligned to the assembly made from the same DNA
  # The panel for these libraries is the genus panel plus the chloroplast contigs of the Purdue assembly
  # (config/panel_extra.tsv): about 650 kb of an Acetabularia chloroplast genome, for which no other reference exists.
  local P=$WORK/panel_ecr s GP GT A=$WORK/reads/elysia-crispata.purdue E=$WORK/reads/elysia-crispata.eggs T=$WORK/reads/elysia-crispata.tampa
  local C=$WORK/asm_census/elysia-crispata-purdue
  GP=$(genome elysia-crispata-purdue); GT=$(genome elysia-crispata-tampa)      # GCA_033675545.1 and GCA_984581095.1 (bgzip + faidx)
  run python3 "$HERE/30c_extend_panel.py" "$PANEL" "$REPO/config/panel_extra.tsv" "$DATA/asm" "$P"
  run bash "$HERE/30b_mask_and_index.sh" "$P"
  for s in elysia-crispata-purdue elysia-crispata-tampa; do
    run bash "$REPO/scripts/22_assembly_organelle_screen.sh" "$PANEL/panel.ont.mmi" "$WORK/asm_screen" "$s=$(genome "$s")"
    run python3 "$HERE/35_assembly_screen.py" "$PANEL/panel.tsv" "$PANEL/panel.fa" "$(genome "$s").fai" "$WORK/asm_screen/$s.paf" "$WORK/asm_screen/$s.screen.json"
  done
  # the chloroplast contigs of the Purdue assembly are hidden in it: their reads are compared with the panel records
  run python3 "$HERE/35b_hide_panel_extra.py" "$REPO/config/panel_extra.tsv" elysia-crispata-purdue "$WORK/asm_screen/elysia-crispata-purdue.screen.json"
  run mkdir -p "$WORK/asm_adapters" "$WORK/asm_census" "$A" "$E" "$T"
  # Nanopore adapter built into contigs of the Purdue assembly
  [ -s "$WORK/asm_adapters/elysia-crispata-purdue.adapters.tsv" ] || run python3 "$HERE/56_adapter_in_assembly.py" "$GP" "$WORK/asm_adapters/elysia-crispata-purdue.adapters.tsv"

  # 1. The egg mass of the Purdue isolate (PRJNA987316). Its read pairs test the places at which the Purdue assembly
  #    joins chloroplast DNA to other DNA within a contig (config/assembly_junctions.tsv; they join the panel as baits),
  #    and the random sample (one pair in 100, of which every second pair is aligned: 0.4-fold coverage) shows which
  #    stretches of the assembly are not in the germline.
  run python3 "$HERE/40_locus_baits.py" - "$GP" "$A/baits.fa" 3000 $(awk -F'\t' 'NR > 1 && $1 == "elysia-crispata" && $2 == "purdue" {printf "%s:%s-%s:%s ", $3, $4, $5, $6}' "$REPO/config/assembly_junctions.tsv")
  run bash "$HERE/30b_mask_and_index.sh" "$P" "$A/baits.fa" org.pur
  PASS2_INDEX_BASES=250M SAMPLE_EVERY=2 PANEL=$P ASM=elysia-crispata-purdue STORE_MAX_MB=5500 once elysia-crispata.eggs lib elysia-crispata.eggs panel.org.pur "SRR25020525" short short --sample 100 --streams 10
  if [ -e "$E/lib.done" ] || [ "${DRY:-0}" = 1 ]; then
    [ -s "$E/junction_tests.json" ] || run python3 "$HERE/41_locus_test.py" "$A/baits.fa.tsv" "$GP" "$E/all_hits.fa" short "$E/junction_tests.json"
    # windows of 8 kb, tested where at least 14 reads are expected (2.1 reads per kb)
    [ -s "$C.tsv" ] || run python3 "$HERE/59_germline_coverage.py" "$E/sample.paf" "$GP" "$C" --window 8000 --min-expected 14
    [ -s "$C.classes.tsv" ] || run python3 "$HERE/60_classify_absent_stretches.py" "$C.tsv" "$GP" "$P/panel.fa" "$HGT/alien/alien_db.dmnd" "$C.classes" --extra "$REPO/config/panel_extra.tsv" \
        --other "tampa=$GT" "sanger=$(genome elysia-crispata)" "alt=$DATA/alt/alt_genomic.fna.gz"
    # the stretches that are absent from the eggs and are chloroplast DNA are hidden as well (four contigs beyond those of the panel)
    run python3 "$HERE/35b_hide_panel_extra.py" "$REPO/config/panel_extra.tsv" elysia-crispata-purdue "$WORK/asm_screen/elysia-crispata-purdue.screen.json" "$C.classes.hidden.tsv"
  fi

  # 2. The adults of the Purdue isolate: Nanopore reads of whole bodies (eleven runs, five flow cells) and of brains.
  #    They are analysed like the older PacBio reads (clr), aligned with the Nanopore preset, and their junctions are
  #    tested for the adapter of the ligation kit (LONG_READ_PLATFORM=ont). (In the published run their first pass was
  #    made before the egg mass was read; the second pass and the analysis were repeated once the four further
  #    chloroplast contigs had been hidden.)
  LONG_READ_PLATFORM=ont PANEL=$P ASM=elysia-crispata-purdue PASS2_PRESET=map-ont STORE_MAX_MB=5500 once elysia-crispata.purdue lib elysia-crispata.purdue panel.dust \
      "SRR25020508 SRR25020509 SRR25020510 SRR25020511 SRR25020513 SRR25020514 SRR25020515 SRR25020516 SRR25020517 SRR25020518 SRR25020522" long clr --sample 200
  LONG_READ_PLATFORM=ont PANEL=$P ASM=elysia-crispata-purdue PASS2_PRESET=map-ont STORE_MAX_MB=5500 once elysia-crispata.brain lib elysia-crispata.brain panel.dust "SRR25020519" long clr --sample 50

  # 3. The long reads of the University of Tampa (PRJNA1337476, library EC2). The library was made by whole-genome
  #    amplification: a quarter of it is mitochondrial DNA, and of the reads that match nothing but the mitochondrial
  #    genome one in 20 is carried on (thin.tsv, 31_select_hits.py).
  local TS=$WORK/reads/elysia-crispata.tampa_short SCT=$WORK/asm_screen/elysia-crispata-tampa.screen.json
  [ -s "$T/thin.tsv" ] || [ "${DRY:-0}" = 1 ] || printf 'slugmt\t20\n' > "$T/thin.tsv"
  PANEL=$P ASM=elysia-crispata-tampa STORE_MAX_MB=5500 once elysia-crispata.tampa lib elysia-crispata.tampa panel.dust "SRR35727887" long hifi --sample 100
  if [ ! -s "$T/circ.contigs.tsv" ]; then
    #  a. In this first analysis every read with a plastid-like match was examined: its reads with an insertion-like
    #     structure are kept apart from the counts of what follows.
    run cp "$T/candidates.json" "$T/candidates.all_plastid.json"
    #  b. Nearly all of the plastid-like reads match no reference closely. From here on one in five of the reads that
    #     match nothing but plastid sequences is carried on: the choice of reads is made again, and the chosen reads
    #     are taken from those kept before.
    [ "${DRY:-0}" = 1 ] || printf 'slugmt\t20\ncp\t5\n' > "$T/thin.tsv"
    run python3 "$HERE/31_select_hits.py" "$T/SRR35727887.p1.paf" "$T/SRR35727887.count" "$T/SRR35727887" long
    run mv "$T/SRR35727887.hits.fa" "$T/SRR35727887.hits.first.fa"
    run sh -c "seqkit grep -j 2 -w 0 -f '$T/SRR35727887.keep.names' '$T/SRR35727887.hits.first.fa' > '$T/SRR35727887.hits.fa'"
    run rm -f "$T/SRR35727887.hits.first.fa"
    big bash "$HERE/33_pass2.sh" "$T" "$GT" "$P/panel.fa" "$SCT" hifi
    #  c. Contigs are assembled from a sample of the plastid-like reads (hifiasm). All of them turn out to be tandem
    #     arrays of a unit of about 5 kb: small circular elements that the amplification copied many times over
    #     (38c_tandem_contigs.py). They become a kind of DNA of their own (circ|...), the library is aligned again
    #     with them, and it is analysed again. (In the published run the second pass was made with the contigs before
    #     they were recognised as arrays; the names of its targets were changed with contigs.rename.tsv instead.)
    run bash "$HERE/38b_plastid_contigs_long.sh" "$T" "$P/panel.fa" "$GT"
    run python3 "$HERE/38c_tandem_contigs.py" "$T" "$REPO/results/intermediate/protein_screen_queries.fa"
    big bash "$HERE/33_pass2.sh" "$T" "$GT" "$P/panel.fa" "$SCT" hifi "$T/own.contigs.fa"
    PANEL=$P ASM=elysia-crispata-tampa analyse elysia-crispata.tampa hifi
    #  d. The barcode genes are read from all reads that span them, including those that were not carried on.
    run python3 "$HERE/36b_barcode_reads_set_aside.py" "$T" "$P/panel.fa" "$T/barcodes/frames.fa" "$T/barcode_reads.fa"
    run sh -c "cat '$T/all_hits.fa' '$T/barcode_reads.fa' > '$T/barcode_all.fa'"
    run python3 "$HERE/36_barcodes.py" reads "$T/barcode_all.fa" "$BC" "$T/barcodes"
    run rm -f "$T/barcode_all.fa"
    #  e. Every read of the library at the place of the one read with an insertion-like structure (step a): the run is
    #     read once more against that place alone, and the reads found are aligned to the unmasked stretch of the assembly.
    run python3 "$HERE/40_locus_baits.py" - "$GT" "$T/baits.fa" 3000 CFHGKC010004091.1:16397-16398:candidate
    run minimap2 -x map-ont -k15 -w10 -d "$T/baits.ont.mmi" "$T/baits.fa"
    run python3 "$HERE/31c_screen_run.py" "$WORK/reads/elysia-crispata.tampa_locus" "$T/baits" SRR35727887 long "$(run_info SRR35727887 spots)"
    run sh -c "samtools faidx '$GT' CFHGKC010004091.1:13398-19398 > '$T/candidate_place.fa'"
    run sh -c "minimap2 -x map-hifi -c --secondary=no '$T/candidate_place.fa' '$WORK/reads/elysia-crispata.tampa_locus/SRR35727887.hits.fa' | cut -f1-12 > '$T/candidate_place_all_reads.paf'"
  fi

  # 4. The read pairs of the University of Tampa (library EC1, registered under the same sample). One pair in sixteen is
  #    plastid DNA, so one in 20 of the pairs that match an organelle genome and nothing else is carried on. The place
  #    of the one long read with an insertion-like structure joins the panel, and the pairs at that place are examined.
  run python3 "$HERE/40_locus_baits.py" - "$GT" "$T/baits.fa" 3000 CFHGKC010004091.1:16397-16398:candidate
  run bash "$HERE/30b_mask_and_index.sh" "$P" "$T/baits.fa" org.tam
  PANEL=$P ASM=elysia-crispata-tampa STORE_MAX_MB=5500 once elysia-crispata.tampa_short lib elysia-crispata.tampa_short panel.org.tam "SRR35727888" short short --sample 800 --plain-keep 20 --match-keep 20 --streams 10
  if [ -e "$TS/lib.done" ] || [ "${DRY:-0}" = 1 ]; then
    [ -s "$TS/candidate_places.json" ] || run python3 "$HERE/41_locus_test.py" "$T/baits.fa.tsv" "$GT" "$TS/all_hits.fa" short "$TS/candidate_places.json"
    # do the two libraries of the Tampa sample come from one animal? Their mitochondrial genomes are compared.
    [ -s "$TS/mt_haplotypes.json" ] || run python3 "$HERE/58_mt_haplotypes.py" "$P/panel.fa" "slugmt|OY973781.1" "$TS/mt_haplotypes.json" "tampa_short=$TS:short" "tampa=$T:hifi"
  fi
}

regulators(){
  # Nuclear genes of plants and algae that regulate the chloroplast or carry signals between chloroplast and nucleus
  # (sigma factors and other plastid transcription factors, retrograde signalling, light signalling, regulation of
  # photosynthesis, protein import, plastid division, the plant clock), searched for in all eleven assemblies; genes
  # of non-animal origin in the four annotated sacoglossan proteomes; and the origin of the second Cry-DASH gene.
  local R=$WORK/regscreen A=$WORK/alienprot DB=$HGT/alien/alien_db.dmnd
  run python3 "$HERE/70_regulator_queries.py" "$REPO/config/plastid_regulators.tsv" "$DATA/algae" "$R/q" \
      "$REPO/results/intermediate/regulator_queries.fa" "$REPO/results/intermediate/regulator_queries.tsv"
  big bash "$HERE/71_regulator_screen.sh" "$R" "$REPO/results/intermediate/regulator_queries.fa" "$DB"
  run python3 "$HERE/73_regulator_summary.py" "$REPO/config/plastid_regulators.tsv" "$REPO/results/intermediate/regulator_queries.tsv" \
      "$R" "$DB" "$REPO/results/genus/plastid_regulators.json" "$REPO/results/genus/plastid_regulators_matrix.tsv"
  big python3 "$HERE/72_alien_proteomes.py" "$DB" "$DATA/algae" "$A" "$R/alien_proteomes.json" \
      ELYCR="$DATA/ensembl/Elysia_crispata-GCA_963854125.1-2024_05-pep.fa.gz" ELYCH="$DATA/refprot/ELYCH.fa.gz" \
      ELYMA="$DATA/refprot/ELYMA.fa.gz" PLAOC="$DATA/refprot/PLAOC.fa.gz"
  run python3 "$HERE/75_alien_candidates.py" "$R/alien_proteomes.json" "$A/proteomes.fa" "$DB" "$R" "$A" "$REPO/results/genus/alien_genes.json"
  run bash "$HERE/74_crydash_tree.sh" "$WORK/crydash" "$REPO/results/genus/crydash_tree"
}

census(){
  # DNA of bacteria and other organisms in every sequence of the eleven assemblies, judged by the rules of
  # foreign_dna.py: one inventory, which contains every foreign sequence that the searches above found.
  big python3 "$HERE/76_foreign_census.py" "$WORK/census" "$HGT/alien/alien_db.dmnd" "$PANEL/panel.fa" "$PANEL/panel.tsv" \
      "$REPO/config/panel_extra.tsv" "$REPO/results" "$REPO/results/genus/foreign_dna_census.json.gz" \
      elysia-crispata elysia-crispata-alt elysia-crispata-purdue elysia-crispata-tampa elysia-chlorotica elysia-marginata \
      elysia-timida elysia-leucolegnote elysia-atroviridis elysia-papillosa elysia-papillosa-alt
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
    all) panel; assembly; reads; controls; loci; proteins; crispata; other_assemblies; crispata_animals; regulators; census; collect ;;
    panel|assembly|reads|controls|loci|proteins|crispata|other_assemblies|crispata_animals|regulators|census|collect) $step ;;
    *) echo "unknown step: $step" >&2; exit 1 ;;
  esac
done
