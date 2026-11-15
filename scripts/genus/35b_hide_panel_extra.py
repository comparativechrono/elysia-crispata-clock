#!/usr/bin/env python3
"""Mark the stretches of an assembly that are not to be treated as DNA of the slug.

usage: 35b_hide_panel_extra.py config/panel_extra.tsv ASSEMBLY_SLUG SCREEN.json [STRETCHES.tsv]

The chloroplast contigs of the Purdue assembly of E. crispata are part of the panel used for the
reads of that species (30c_extend_panel.py). When reads are aligned to the Purdue assembly and the
panel together (33_pass2.sh), the same sequence would be present twice. The stretches named in
panel_extra.tsv for this assembly are therefore written to SCREEN.json as "hidden"; 33_pass2.sh
replaces them by N, so that their reads align to the panel records and count as plastid DNA.

STRETCHES.tsv (60_classify_absent_stretches.py: sequence, start, end, class) lists further
stretches of the assembly that are absent from egg DNA and are chloroplast or bacterial DNA. They
are hidden too: a read that runs from a chloroplast contig into one of them would otherwise count
as chloroplast DNA joined to DNA of the slug. Their reads match nothing and are not counted.
"""
import sys, csv, json
extra, slug, screen = sys.argv[1:4]
more = sys.argv[4] if len(sys.argv) > 4 else None
s = json.load(open(screen))
s["hidden"] = [dict(contig=r["sequence"], start=int(r["start"]), end=int(r["end"]), panel_record=r["id"])
               for r in csv.DictReader(open(extra), delimiter="\t") if r["assembly"] == slug]
n_panel = len(s["hidden"])
if more:
    for l in open(more):
        c, a, b, cls = l.rstrip("\n").split("\t")[:4]; a, b = int(a), int(b)
        # a contig with a record in the panel is left as that record defines it: its limits were set by hand,
        # whereas the census has a resolution of one window
        if any(h["contig"] == c for h in s["hidden"][:n_panel]): continue
        s["hidden"].append(dict(contig=c, start=a, end=b, reason=f"absent from egg DNA; {cls} DNA"))
hidden = {h["contig"] for h in s["hidden"]}
# a hidden contig is not also an organelle sequence of the assembly
s["whole"] = [w for w in s["whole"] if w["contig"] not in hidden]
json.dump(s, open(screen, "w"), indent=1)
print(slug, len(s["hidden"]), "stretches hidden,", sum(h["end"] - h["start"] for h in s["hidden"]), "bp;", n_panel, "of them are records of the panel")
