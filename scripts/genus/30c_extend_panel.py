#!/usr/bin/env python3
"""Extend the reference panel with plastid sequences that exist only as contigs of an assembly.

usage: 30c_extend_panel.py PANEL_DIR config/panel_extra.tsv ASM_DIR OUT_DIR

The Purdue assembly of E. crispata contains about 650 kb of chloroplast DNA of an Acetabularia
species (and one contig from another siphonous alga) as separate contigs. No complete plastid
genome of Acetabularia is public, so these contigs are the only reference for the kleptoplasts
of that animal. Each row of panel_extra.tsv names a stretch of an assembly sequence
(assembly, sequence, start, end; 0-based start) and the panel record it becomes (id, class,
accession, length, description, source). OUT_DIR/panel.fa and OUT_DIR/panel.tsv are the panel of
PANEL_DIR with these records added; index them with 30b_mask_and_index.sh.
"""
import sys, os, csv, shutil, subprocess
panel, extra, asm, out = sys.argv[1:5]
ST = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/samtools"
os.makedirs(out, exist_ok=True)
shutil.copyfile(os.path.join(panel, "panel.fa"), os.path.join(out, "panel.fa"))
shutil.copyfile(os.path.join(panel, "panel.tsv"), os.path.join(out, "panel.tsv"))
rows = list(csv.DictReader(open(extra), delimiter="\t"))
with open(os.path.join(out, "panel.fa"), "a") as fa, open(os.path.join(out, "panel.tsv"), "a") as tsv:
    for r in rows:
        g = os.path.join(asm, r["assembly"], "genome.fa.gz")
        res = subprocess.run([ST, "faidx", g, f'{r["sequence"]}:{int(r["start"]) + 1}-{r["end"]}'], capture_output=True, text=True, check=True).stdout
        seq = "".join(res.split("\n")[1:]).upper()
        assert len(seq) == int(r["length"]), (r["id"], len(seq), r["length"])
        fa.write(f'>{r["id"]}\n')
        for i in range(0, len(seq), 80): fa.write(seq[i:i + 80] + "\n")
        tsv.write("\t".join(r[k] for k in ("id", "class", "accession", "length", "description", "source")) + "\n")
        print(r["id"], len(seq), "bp")
print(len(rows), "records added;", out)
