#!/usr/bin/env python3
"""Tell the protein screen (scripts/28_protein_screen.py) what the DNA screen already found.

usage: 42_protein_screen_inputs.py ASM_SCREEN_DIR PROTSCREEN_DIR slug ...
From ASM_SCREEN_DIR/<slug>.screen.json (35_assembly_screen.py) writes, in PROTSCREEN_DIR,
  embedded_plastid.json   plastid-like stretches inside other sequences (>= 90% identity)
  assembly_screen.json    sequences that are plastid DNA from end to end
so that protein matches on them are labelled as such.
"""
import sys, os, json
sdir, pdir = sys.argv[1:3]; slugs = sorted(sys.argv[3:])
emb, whole = {}, {}
for s in slugs:
    S = json.load(open(os.path.join(sdir, s + ".screen.json")))
    emb[s] = [dict(contig=L["contig"], start=L["start"], end=L["end"]) for L in S["loci"] if L["cls"] == "plastid" and L["identity"] >= 0.90]
    whole[s] = dict(whole_list=[dict(contig=w["contig"], cls="cp") for w in S["whole"] if w["cls"] == "plastid"])
os.makedirs(pdir, exist_ok=True)
json.dump(emb, open(os.path.join(pdir, "embedded_plastid.json"), "w"))
json.dump(whole, open(os.path.join(pdir, "assembly_screen.json"), "w"))
print({s: (len(emb[s]), len(whole[s]["whole_list"])) for s in slugs})
