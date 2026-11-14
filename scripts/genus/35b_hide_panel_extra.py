#!/usr/bin/env python3
"""Mark the stretches of an assembly that the extended panel holds as records of their own.

usage: 35b_hide_panel_extra.py config/panel_extra.tsv ASSEMBLY_SLUG SCREEN.json

The chloroplast contigs of the Purdue assembly of E. crispata are part of the panel used for the
reads of that species (30c_extend_panel.py). When reads are aligned to the Purdue assembly and the
panel together (33_pass2.sh), the same sequence would be present twice. The stretches named in
panel_extra.tsv for this assembly are therefore written to SCREEN.json as "hidden"; 33_pass2.sh
replaces them by N, so that their reads align to the panel records and count as plastid DNA.
"""
import sys, csv, json
extra, slug, screen = sys.argv[1:4]
s = json.load(open(screen))
s["hidden"] = [dict(contig=r["sequence"], start=int(r["start"]), end=int(r["end"]), panel_record=r["id"])
               for r in csv.DictReader(open(extra), delimiter="\t") if r["assembly"] == slug]
hidden = {h["contig"] for h in s["hidden"]}
# a hidden contig is not also an organelle sequence of the assembly
s["whole"] = [w for w in s["whole"] if w["contig"] not in hidden]
json.dump(s, open(screen, "w"), indent=1)
print(slug, len(s["hidden"]), "stretches hidden,", sum(h["end"] - h["start"] for h in s["hidden"]), "bp")
