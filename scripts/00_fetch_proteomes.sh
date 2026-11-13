#!/usr/bin/env bash
# Fetch UniProt (reference) proteomes used for orthology assignment.
set -uo pipefail
OUT=${1:?outdir}; mkdir -p "$OUT"; cd "$OUT"
while read -r tag upid extra; do
  [ -z "$tag" ] && continue
  [ -s "$tag.fa.gz" ] && continue
  q="(proteome:$upid)$extra"
  curl -sS --fail --retry 4 --retry-delay 5 -m 1500 -G "https://rest.uniprot.org/uniprotkb/stream" \
     --data-urlencode "query=$q" --data-urlencode "format=fasta" --data-urlencode "compressed=true" -o "$tag.fa.gz.part" \
     && mv "$tag.fa.gz.part" "$tag.fa.gz" && echo "ok $tag $(zcat $tag.fa.gz | grep -c '>')" || echo "FAIL $tag"
done <<'LIST'
HUMAN UP000005640  AND (reviewed:true)
DROME UP000000803
DANPL UP000007151
TRICA UP000007266
MAGGI UP000005408
APLCA UP000694888
LOTGI UP000030746
BIOGL UP001165740
POMCA UP000245119
MIZYE UP000242188
ELYCH UP000271974
ELYMA UP000762676
PLAOC UP000735302
ELYCR_PURDUE UP001283361
CAPTE UP000014760
NEMVE UP000001593
BRAFL UP000001554
LIST
echo PROTEOMES_DONE
