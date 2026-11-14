#!/usr/bin/env bash
# Mask low-complexity sequence in the reference panel and index it.
#   PANEL_DIR/panel.fa (30_build_panel.py)
#     -> panel.ont.mmi        unmasked, k 15 / w 10: used to screen assemblies in 10 kb windows
#     -> panel.dust.fa        low-complexity sequence (dustmasker, level 20) replaced by N
#     -> panel.dust.ont.mmi   first pass, long reads      panel.dust.sr.mmi   first pass, read pairs
#   With EXTRA.fa and NAME:  panel.NAME.fa = panel.dust.fa + EXTRA.fa, indexed the same way. Used
#   for E. chlorotica, where the places to be tested with reads (40_locus_baits.py) were added to
#   the panel so that the first pass keeps their reads (NAME = chl).
# usage: 30b_mask_and_index.sh PANEL_DIR [EXTRA.fa NAME]
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
D=$1; EXTRA=${2:-}; NAME=${3:-}
[ -s "$D/panel.ont.mmi" ] || minimap2 -x map-ont -k15 -w10 -d "$D/panel.ont.mmi" "$D/panel.fa" 2>/dev/null
if [ ! -s "$D/panel.dust.fa" ]; then
  dustmasker -in "$D/panel.fa" -level 20 -outfmt fasta 2>/dev/null | awk '/^>/{print; next}{gsub(/[acgtn]/,"N"); print}' > "$D/panel.dust.fa"
fi
index(){ [ -s "$1.ont.mmi" ] || minimap2 -x map-ont -k15 -w10 -d "$1.ont.mmi" "$1.fa" 2>/dev/null; [ -s "$1.sr.mmi" ] || minimap2 -x sr -d "$1.sr.mmi" "$1.fa" 2>/dev/null; }
index "$D/panel.dust"
if [ -n "$EXTRA" ]; then
  cat "$D/panel.dust.fa" "$EXTRA" > "$D/panel.$NAME.fa"; index "$D/panel.$NAME"
fi
echo "panel indexed in $D"
