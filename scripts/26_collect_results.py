#!/usr/bin/env python3
"""Gather the kleptoplast and gene-transfer results into two JSON files for the site.

usage: 26_collect_results.py work_dir results_dir
Reads the intermediate outputs under work_dir (hgt/, klepto/) and writes
results/kleptoplast.json and results/hgt.json, plus flat TSV copies of the main tables.
"""
import sys, os, json, csv, collections

wd, res = sys.argv[1:3]
H = os.path.join(wd, "hgt"); K = os.path.join(wd, "klepto")
os.makedirs(res, exist_ok=True)
RUNS = {"ERR11843472": (2480541, 21880495092), "ERR12353011": (1639226, 14442137611), "ERR12356115": (1396905, 10523946025)}
TOTAL_READS = sum(v[0] for v in RUNS.values()); TOTAL_BP = sum(v[1] for v in RUNS.values())

def jload(p, default=None):
    return json.load(open(p)) if os.path.exists(p) else default

# ---------------------------------------------------------------- kleptoplasts
classes = {}
rows = list(csv.DictReader(open(os.path.join(H, "sel", "bait_summary.tsv")), delimiter="\t"))
tot = [r for r in rows if r["run"] == "total"][0]
for k in ("cp", "bac", "slugmt", "algmt", "algrDNA", "cpx"):
    classes[k] = dict(reads=int(tot[k + "_reads"]), bp=int(tot[k + "_bp"]))
conf = jload(os.path.join(H, "sel", "bait_confident.json"), {})
for k in ("cp", "bac", "slugmt"):          # confident matches only: <= 3 % divergence over >= 1 kb
    if k in conf: classes[k] = dict(reads=conf[k]["reads"], bp=conf[k]["bp"], all_reads=conf[k]["all_reads"])
per_run = [dict(run=r["run"], reads=RUNS[r["run"]][0], bp=RUNS[r["run"]][1], cp_bp=int(r["cp_bp"]), bac_bp=int(r["bac_bp"]), slugmt_bp=int(r["slugmt_bp"])) for r in rows if r["run"] != "total"]
assign = jload(os.path.join(K, "plastome_assignment.json"), {})
barcodes = []
p = os.path.join(K, "census", "barcode_haplotypes.tsv")
if os.path.exists(p):
    for r in csv.DictReader(open(p), delimiter="\t"):
        barcodes.append(dict(barcode=r["barcode"], haplotype=r["haplotype"], reads=int(r["reads"]), fraction=float(r["fraction"]), length=int(r["length"]),
                             top_hit=r["top_hit"], identity=float(r["top_hit_identity"] or 0), alnlen=int(r["top_hit_alnlen"] or 0),
                             species_near_best=r["species_within_0.5pct_of_best"], sequence=r["sequence"]))
klepto = dict(
    library=dict(reads=TOTAL_READS, bp=TOTAL_BP, runs=per_run),
    classes=classes,
    plastomes=[dict(taxon="Penicillus capitatus", accession="OZ373742.1", length=84040), dict(taxon="Penicillus lamourouxii", accession="OZ372796.1", length=89458)],
    assignment=assign.get("assignment", []), unassigned=assign.get("weak", {}),
    barcodes=barcodes,
    # identities measured in this analysis (see scripts/13_barcode_census.py and the notes in docs/methods.html)
    reference_identity=dict(
        capitatus=dict(rbcL=dict(identity=100.0, length=1253, accession="OR665168.1"), tufA=dict(identity=100.0, length=813, accession="PV013973.1"),
                       next_species=dict(name="Udotea unistratea", rbcL=99.43, tufA=98.02)),
        lamourouxii=dict(rbcL=dict(identity=100.0, length=562, accession="DQ469333.1", note="the only P. lamourouxii rbcL sequence in GenBank"),
                         tufA=None, next_species=dict(name="Penicillus dumetosus", rbcL=97.79, tufA=98.60)),
        between_types=dict(rbcL=96.01, tufA=92.70)),
    minor=[dict(taxon="Udotea flabellum", rbcL="99.85–100% (3 reads)", tufA="98.93–100% (12 reads)", call="Udotea flabellum"),
           dict(taxon="Halimeda sp.", rbcL="97.6% to H. cylindracea (5 reads)", tufA="97.1–97.2% to H. cylindracea (8 reads)", call="a Halimeda without a species-level reference")],
    nuclear=dict(rdna_reference_sequences=1249, reads_matching_algal_rdna=0, reads_first_flagged=classes["algrDNA"]["reads"],
                 algal_mito_reads=3, algal_mito_bp=4169),
)
json.dump(klepto, open(os.path.join(res, "kleptoplast.json"), "w"), indent=1)
with open(os.path.join(res, "kleptoplast_barcodes.tsv"), "w") as fo:
    fo.write("barcode\thaplotype\treads\tfraction\tlength\tbest_genbank_hit\tidentity\taligned_bp\tspecies_within_0.5pct\tsequence\n")
    for b in barcodes:
        fo.write("\t".join(str(b[k]) for k in ("barcode", "haplotype", "reads", "fraction", "length", "top_hit", "identity", "alnlen", "species_near_best", "sequence")) + "\n")

# ---------------------------------------------------------------- gene transfer
chim = jload(os.path.join(H, "chimera", "chimera_summary.json"), {})
emb = jload(os.path.join(H, "asm_screen", "embedded_plastid.json"), {})
scr = jload(os.path.join(H, "asm_screen", "assembly_screen.json"), {})
bnd = jload(os.path.join(H, "boundary", "boundary_evidence.json"), [])
alien = jload(os.path.join(H, "alien", "alien_candidates_annotated.json")) or jload(os.path.join(H, "alien", "alien_candidates.json"), {})
prot = jload(os.path.join(H, "protscreen", "protein_screen.json"), {})
bmap = {(b["assembly"], b["contig"], b["seg_start"]): b for b in bnd}
uniq = {}
wt = os.path.join(H, "boundary", "windows.tsv")
if os.path.exists(wt):
    rows = [l.rstrip("\n").split("\t") for l in open(wt)]
    if len(rows[0]) > 9:
        for f in rows[1:]: uniq[(f[1], f[2], int(f[3]))] = int(f[9])
GROUP = {"OZ373742.1": "Penicillus capitatus", "OZ287723.1": "Penicillus capitatus", "OZ372796.1": "Penicillus lamourouxii", "OZ494339.1": "Penicillus lamourouxii"}
fragments = []
for asm in ("sanger", "alt"):
    for r in emb.get(asm, []):
        b = bmap.get((asm, r["contig"], r["start"]))
        open_side = "left" if (r["left"]["nearest_break"] in ("gap", "sequence end") and (r["left"]["distance"] or 0) < 50) else "right"
        closed = r["right"] if open_side == "left" else r["left"]
        fragments.append(dict(assembly=asm, contig=r["contig"], start=r["start"], end=r["end"], length=r["length"], identity=r["identity"],
                              alga=GROUP.get(r["plastome"], r["plastome"]), plastome=r["plastome"], plastome_start=r["plastome_start"], plastome_end=r["plastome_end"],
                              genes=r["genes"], abuts=(r["left"] if open_side == "left" else r["right"])["nearest_break"],
                              crossing_reads_pass2=closed["reads_crossing"], crossing_molecules_pass2=closed["independent_molecules"],
                              # junction as assembled: pass 2 saw every chimeric read of all three runs, the boundary
                              # pass every read of two runs; report the larger count
                              through_reads=max(closed["reads_crossing"], b["through_reads"] if b else 0),
                              through_molecules=max(closed["independent_molecules"], b["through_independent"] if b else 0),
                              through_reads_boundary_pass=b["through_reads"] if b else None, through_molecules_boundary_pass=b["through_independent"] if b else None,
                              unique_flank_bp=uniq.get((asm, r["contig"], r["start"])),
                              diverging_reads=b["diverge_reads"] if b else None, diverging_molecules=b["diverge_independent"] if b else None,
                              flank_reads=b["flank_reads"] if b else None))
low = {}
for asm, d in scr.items():
    e_ = d["embedded"]
    lowc = [x for x in e_ if not (x["ident"] >= 0.97 and x["length"] >= 200)]
    c = collections.Counter((x["target"].split("|")[1], x["tstart"] // 300 * 300) for x in lowc)
    low[asm] = dict(n=len(lowc), top=[dict(target=t, window=w, n=n) for (t, w), n in c.most_common(3)],
                    whole_contigs=d["whole_contigs"], exact=len([x for x in e_ if x["ident"] >= 0.97 and x["length"] >= 200]))
pairs = {p["pair"]: p for p in chim.get("pairs", [])}
ct = chim.get("class_totals", {})
frac = lambda k: classes[k]["bp"] / TOTAL_BP
cp_reads, bac_reads = ct.get("cp", 0), ct.get("bac", 0)
rate_cp = pairs.get("cp+slug", {}).get("reads", 0) / cp_reads if cp_reads else 0
rate_bac = pairs.get("bac+slug", {}).get("reads", 0) / bac_reads if bac_reads else 0
hgt = dict(
    chimera=dict(pairs=chim.get("pairs", []), class_totals=ct, rate_cp_slug=rate_cp, rate_bac_slug=rate_bac,
                 expected_cp_bac=cp_reads * rate_cp * frac("bac") / (1 - frac("cp")) if cp_reads else 0,
                 expected_cp_bac_from_bacteria=bac_reads * rate_bac * frac("cp") / (1 - frac("bac")) if bac_reads else 0,
                 observed_cp_bac=pairs.get("bac+cp", {}).get("reads", 0),
                 recurrent=chim.get("recurrent", []), flanked=chim.get("flanked", 0), reads_examined=chim.get("reads_examined", 0),
                 **(jload(os.path.join(H, "chimera", "chimera_extra.json"), {}))),
    fragments=fragments, assembly_screen=low, alien=alien, protein_screen=prot,
    diverged=jload(os.path.join(H, "asm_screen", "diverged_candidates.json"), {}).get("summary", {}),
    weak_reads=jload(os.path.join(H, "sel", "weak_plastid_reads.json"), {}),
    foreign=dict(blast=jload(os.path.join(H, "protscreen", "purdue_foreign_vs_dasycladales.json"), {}).get("best", []),
                 same_genes=jload(os.path.join(H, "protscreen", "purdue_foreign_vs_dasycladales.json"), {}).get("same_genes", {}),
                 joins=jload(os.path.join(H, "protscreen", "jointest", "join_test.json"), []),
                 join_runs=jload(os.path.join(H, "protscreen", "jointest", "runs.json"), {})),
)
json.dump(hgt, open(os.path.join(res, "hgt.json"), "w"), indent=1)
with open(os.path.join(res, "plastid_fragments_in_assembly.tsv"), "w") as fo:
    cols = ["assembly", "contig", "start", "end", "length", "identity", "alga", "plastome", "plastome_start", "plastome_end", "genes", "abuts",
            "through_reads", "through_molecules", "unique_flank_bp", "diverging_reads", "diverging_molecules", "flank_reads"]
    fo.write("\t".join(cols) + "\n")
    for f in fragments:
        fo.write("\t".join(",".join(f[c]) if c == "genes" else str(f[c]) for c in cols) + "\n")
print("kleptoplast.json:", len(barcodes), "barcode haplotypes;", len(klepto["assignment"]), "taxa")
print("hgt.json:", len(fragments), "plastid fragments;", len(chim.get("pairs", [])), "chimera pair classes; alien candidates:", len(alien.get("candidates", [])) if alien else "n/a")
