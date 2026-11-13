#!/usr/bin/env bash
set -uo pipefail
while read -r tag upid extra; do
  [ -z "$tag" ] && continue; [ -s "$tag.fa.gz" ] && continue
  curl -sS --fail --retry 4 --retry-delay 5 -m 1500 -G "https://rest.uniprot.org/uniprotkb/stream" --data-urlencode "query=(proteome:$upid)$extra" --data-urlencode "format=fasta" --data-urlencode "compressed=true" -o "$tag.fa.gz.part" && mv "$tag.fa.gz.part" "$tag.fa.gz" && echo "ok $tag $(zcat $tag.fa.gz | grep -c '>')" || echo "FAIL $tag"
done <<'LIST'
OSTQU UP000708148
BRYSP UP001172620
CHLRE UP000006906
VOLCA UP000001058
CHLVA UP000008141
OSTTA UP000009170
MICPC UP000001876
ARATH UP000006548  AND (reviewed:true)
PHATC UP000000759
ECTSI UP000002630
CYAME UP000007014
CHOCR UP000012073
SYNY3 UP000001425
LIST
curl -sS --fail --retry 3 -m 1500 -o uniprot_sprot.fasta.gz "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" && echo "ok sprot $(zcat uniprot_sprot.fasta.gz | grep -c '>')" || echo "FAIL sprot"
echo ALGAE_DONE
