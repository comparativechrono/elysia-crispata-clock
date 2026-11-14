# Elysia clocks

This repository contains an analysis of circadian clock genes, the algal source of the
kleptoplasts and possible algal gene transfer in sea slugs of the genus *Elysia*, based on public
genome data. The work began with the lettuce sea slug *Elysia crispata* (Sanger Institute
specimen xgElyCris1, BioProject
[PRJEB67629](https://www.ebi.ac.uk/ena/browser/view/PRJEB67629); Middlebrooks et al. 2026,
[Wellcome Open Research 11:148](https://wellcomeopenresearch.org/articles/11-148)) and now covers
every *Elysia* species for which a nuclear genome assembly is available.

The results are presented as a website in `docs/`, which is served by GitHub Pages:
**https://comparativechrono.github.io/elysia-crispata-clock/**

No rhythms were measured. Gene functions are inferred from homology with *Drosophila* and
mammals. The analyses were carried out independently of the groups that produced the genomes and
have not been peer reviewed.

## Results for the genus

A nuclear genome assembly is available for seven species. The 47 clock genes catalogued in
*E. crispata* were searched for in each assembly (`scripts/genus/`), and the encoded proteins were
aligned and compared.

**Clock genes.** Number of genes of each status in the reference assembly of each species:

| Species | Reference assembly | Intact | Complete, with a frameshift or stop codon | Partial or fragmented | Not found |
|---|---|---|---|---|---|
| *Elysia crispata* | GCA_963854125.1 (chromosome) | 30 | 4 | 5 | 8 |
| *Elysia chlorotica* | GCA_003991915.1 (scaffold) | 39 | 5 | 3 | 0 |
| *Elysia marginata* | GCA_019649035.1 (scaffold) | 33 | 8 | 4 | 2 |
| *Elysia timida* | GCA_043644045.1 (chromosome) | 44 | 1 | 1 | 1 |
| *Elysia leucolegnote* | GCA_052815195.1 (chromosome) | 43 | 1 | 2 | 1 |
| *Elysia atroviridis* | GCA_059052615.1 (chromosome) | 44 | 1 | 1 | 1 |
| *Elysia papillosa* | GCA_977015215.2 (contig) | 43 | 0 | 3 | 1 |

- Of the 47 genes, 27 are complete in the reference assembly of all seven species, and the core
  loop is complete in every species. Where a gene is partial or missing, the cause is usually the
  assembly: the two assemblies built mainly from short reads contain the largest numbers of such
  genes.
- None of the genes that could not be found in *E. crispata* (Dec1/2, Fbxl3/21, Tyf, Cipc, Ciart,
  Aanat, Asmt) was found in any other species.
- The number of Go-opsin genes varies between species (one to three); Go-opsin 1 of
  *E. crispata* has a one-to-one orthologue only in *E. chlorotica*.
- Of the 54 pairs of clock genes that lie on the same chromosome in *E. crispata*, 53 also lie
  on the same chromosome in *E. timida*, *E. leucolegnote* and *E. atroviridis*.
- A tree inferred from 35 clock proteins (21,190 aligned positions) groups *E. crispata* with
  *E. chlorotica*, and *E. timida* with *E. leucolegnote* and *E. atroviridis*; *E. papillosa*
  is the earliest-branching species. All nodes have full support, but a tree based on clock
  proteins alone is intended as a guide to the relationships among these species, not as a
  phylogeny of the genus.

For species without a public annotation, and for genes that are missing from an annotation, the
gene models were obtained by aligning the *E. crispata* protein to the genome. Such models can
lack exons that are absent from the *E. crispata* gene. The annotation of *E. timida* is not in
GenBank and the server that hosts it could not be reached, so all *E. timida* models are of this
kind. In the GenBank sequence of *E. leucolegnote*, 35% of the bases are masked (N).

**Kleptoplasts.** The plastid DNA among the reads of each sequenced animal identifies the algae
on which it had fed. All reads of every library were analysed; the read sets are listed in
`config/species_reads.tsv`.

| Species | Algal source of the kleptoplasts | Source of reads | Plastid DNA as a share of sequenced bases |
|---|---|---|---|
| *E. chlorotica* | *Vaucheria litorea* | two adults, and eggs | 3.5% and 11.0%; eggs: 13 read pairs in 77 million |
| *E. marginata* | *Bryopsis* (*B. hypnoides* type, PQ652170.1) | assembled individual | 0.26% |
| *E. timida* | *Acetabularia acetabulum* (rbcL and tufA match its GenBank records) | assembled individuals | 0.0054% or more: only 49 kb of this plastid genome is public |
| *E. leucolegnote* | a *Boodleopsis* species (97 to 99% identical to NC_039521.1) in both individuals; a second alga without a close reference sequence, at a quarter of the abundance of the *Boodleopsis* in the Shenzhen individual and about half in the Hainan individual; and, in the Hainan individual, a third plastid genome, also without a reference sequence, as abundant as the *Boodleopsis* | two individuals other than the one assembled | 0.28% and 0.074% or more |
| *E. atroviridis* | *Bryopsis* of the same type, as two strains | assembled individual | 1.9% |
| *E. papillosa* | *Penicillus lamourouxii* (OZ494339.1), in both individuals | assembled individual and a second individual | 0.12% and 0.014% |

**Gene transfer.** There is no evidence of algal DNA in the germline genome of any of these
species: no assembly contains algal DNA that is supported by its reads, and egg DNA of
*E. chlorotica* contains none. Kleptoplasts are located in somatic cells, however, so an
insertion of chloroplast DNA would be present in one cell or a few and would be represented in
the reads by single molecules; reads from the same locus that lack the junction are not evidence
against it. All reads of every library were therefore searched for single molecules in which
chloroplast DNA is joined to slug nuclear DNA. These molecules are listed, one per row, in
`results/genus/` and `docs/downloads/` (`*junctions.tsv`, `*junction_pairs.tsv.gz`,
`*insertion_shapes.tsv`).

| Library | Chloroplast reads (or read pairs) joined to nuclear DNA | Mitochondrial reads (or read pairs) joined to nuclear DNA | Nuclear DNA joined to unrelated nuclear DNA | Reads with chloroplast DNA flanked by nuclear DNA on both sides |
|---|---|---|---|---|
| *E. crispata*, PacBio HiFi reads (amplified) | 4,452 of 136,293 (3.3%) | bacterial DNA used instead: 3.2% | 4.2% | 0 |
| *E. chlorotica*, adult, PacBio RS reads | 2,397 of 43,514 (5.5%) | 113 of 2,116 (5.3%) | 5.3% | 9 |
| *E. chlorotica*, second adult, paired-end short reads (estimated from one pair in 50) | 71,401 of 57,986,168 (0.12%) | 976 of 1,113,806 (0.088%) | 0.99% | not applicable (short reads) |
| *E. chlorotica*, eggs, paired-end short reads | 0 of 13 (0%) | 1 of 104,147 (0.001%) | 0.91% | not applicable (short reads) |
| *E. atroviridis*, PacBio HiFi reads | 14 of 69,026 (0.020%) | 0 of 2,125 (0%) | 2.3% | 0 |
| *E. papillosa*, assembled individual, PacBio HiFi reads | 6 of 6,533 (0.092%) | 4 of 2,012 (0.20%) | 2.9% | 0 |
| *E. papillosa*, second individual, PacBio HiFi reads (amplified) | 3 of 1,326 (0.23%) | 29 of 7,954 (0.36%) | 5.3% | 0 |
| *E. timida*, PacBio HiFi reads (amplified) | 15 of 2,598 (0.58%) | 75 of 14,113 (0.53%) | 1.1% | 0 |
| *E. marginata*, paired-end short reads | 78 of 759,795 (0.010%) | 1 of 35,000 (0.0029%) | 1.7% | not applicable (short reads) |
| *E. leucolegnote*, Shenzhen individual, paired-end short reads | 1,667 of 125,085 (1.3%) | 593 of 66,053 (0.90%) | 1.2% | not applicable (short reads) |
| *E. leucolegnote*, Hainan individual, paired-end short reads | 491 of 68,603 (0.72%) | 166 of 32,441 (0.51%) | 0.83% | not applicable (short reads) |

- Such molecules were found in every library. However, every library also contains chimeras that
  cannot reflect the genome, for example between the slug's mitochondrial DNA and its
  chromosomes, or between two different chromosomes. Chimeras of this kind form during library
  preparation, and a genuine junction that is sequenced only once cannot be distinguished from
  them. In the amplified libraries, many of the junctions contain the amplification adapter.
- Frequency. In the long-read libraries, chloroplast DNA is joined to nuclear DNA as frequently
  as mitochondrial DNA is. In the short-read libraries that contain enough junctions of both
  kinds for a comparison, it is joined 1.3 to 1.5 times as frequently (read pairs at the
  bacterial ribosomal RNA genes of *E. leucolegnote* were excluded). Insertions in somatic cells
  would produce such an excess, but the reads that span a junction indicate that it arises
  during library preparation. A read of this kind starts with the short segment at the end of its
  fragment. This terminal segment is chloroplast DNA in 98 of 170 junctions with chloroplast DNA,
  but mitochondrial DNA in only 14 of 51 junctions with mitochondrial DNA. If the junctions were
  insertions in chromosomes that had been fragmented at random, the proportion would be the same
  for both organelles (two controls, `scripts/genus/55_end_piece_control.py`, give a value close
  to 45%). The difference shows that the fragments of the two organelle DNAs differ in size in
  these libraries, so the two are not joined to other DNA at the same rate during library
  preparation, and the excess cannot be interpreted as evidence of insertions. If the whole
  excess were due to insertions, it would correspond to 10 to 150 insertions in every cell of
  these animals.
- Read structure. Chloroplast DNA flanked on both sides by contiguous nuclear sequence is the one
  structure that the joining of two fragments does not produce. Nine reads have this structure,
  all of them among the error-prone PacBio RS reads of *E. chlorotica*, which also include two
  reads with mitochondrial DNA in the same arrangement and two with the reciprocal structure
  (nuclear DNA flanked by chloroplast DNA). Each is a single molecule at a different
  locus. The PacBio HiFi reads of the other four species contain no such read.
- Recurrence. In the long-read libraries, no junction is present in two independent molecules.
  The reads were compared with one another around every junction, so that junctions within
  repeats were included. Pairs of molecules with junctions a few bases apart occur in
  *E. chlorotica* (two pairs, at a satellite repeat; one of them corresponds to the junction
  present in the assembly) and in *E. crispata* (three pairs, at simple repeats).
- In the short-read libraries, PCR duplicates were identified from the read sequences, and two
  molecules were counted as sharing a junction when their nuclear ends lie within 1 kb of each
  other and their organelle ends lie within 1 kb of each other. Recurrent junctions were found
  at four loci. Three of them do not involve organelle DNA: a slug transfer RNA gene next to an
  assembly gap and the ribosomal RNA genes of a bacterial sequence, both in the two
  *E. leucolegnote* individuals, and a ribosomal RNA gene next to a gap in a bacterial sequence
  in *E. marginata*. The fourth is in the egg library of *E. chlorotica*: 88 bases of
  mitochondrial DNA inserted between two adjacent bases of RQTK01000018.1, spanned by two
  independent molecules and absent from the assembly. This is an inherited insertion of
  mitochondrial DNA, and it shows that the test can detect an insertion of organelle DNA. No
  recurrent junction involves chloroplast DNA.
- Detection limit. Reads long enough to contain an insertion together with both flanks cover
  each genome 29 to 68 times, which is equivalent to the DNA of 15 to 34 cells. About three
  reads with an insertion would be expected if there were one insertion per 5 to 11 cells.
  Frequencies such as those measured in tobacco (about one in 16,000 pollen grains and one in
  five million leaf cells) are far below what can be detected by sequencing a whole animal.
  Among the short-read libraries, *E. marginata* provides a limit of its own: chloroplast DNA
  is joined to nuclear DNA so rarely in this library (78 read pairs) that the number is about
  what one insertion per cell would add.

Foreign DNA in the assemblies:

- *E. chlorotica*: five small plastid scaffolds (8.3 kb), plastid regions within two larger
  scaffolds (1.6 kb and 230 bp) and four regions of mitochondrial DNA. The 1.6 kb region, at
  the start of RQTK01001031.1, is joined to a slug satellite repeat; one long read contains this
  junction and a second molecule has a junction 47 bases away. No short read from the second
  adult or from the eggs contains it. No read from any library crosses either boundary of the
  230 bp region, and two long reads span the locus without it. No short read crosses a boundary
  of the mitochondrial regions. None of the *Vaucheria* nuclear genes reported in early studies
  (psbO, prk, lhc, chlD, chlG, chlH) has a match in the assembly. Fourteen scaffolds (525 kb)
  derive from a bacterium related to *Francisella*.
- *E. marginata*: one 3.7 kb scaffold consists of *Bryopsis* plastid DNA. About 2.0 Mb on 23
  scaffolds is bacterial (mycoplasma-like, flavobacteria-like and others), and the gene set
  published with the assembly includes genes on these scaffolds.
- *E. leucolegnote*: the unplaced sequence JBQKML010000381.1 begins with 863 bp of plastid DNA
  joined to a slug repeat. In two other individuals, the same plastid sequence continues as
  plastid DNA at that point (87 and 77 reads diverge from the assembly there, and none spans the
  junction). The reads of the assembled individual are not public, so the junction itself
  cannot be tested. Four sequences (1.9 Mb) are bacterial.
- *E. papillosa*: contig CDRNNR020002991.1 (30.9 kb) is a fragment of an algal mitochondrial
  genome.
- *E. timida* and *E. atroviridis*: no foreign DNA of these kinds was found in the assembly.

In the protein-level search of these species, the gene set published with an assembly was not
used as evidence that a sequence in that assembly is of animal origin, because the published
gene sets of *E. chlorotica* and *E. marginata* include genes on bacterial and plastid sequences.

## Results for *Elysia crispata*

**Clock genes.** The catalogue contains 47 genes: the core loop (Clock, Bmal1, Period, Timeless,
Timeout, Cry1, Cry2), the second loop (Pdp1, Vrille, E75, Hr3, Cwo), kinases, phosphatases,
degradation factors, output genes, photolyases and nine opsins. All 47 are present in the
species, but the Sanger Institute assembly xgElyCris1.1 is incomplete (BUSCO 73.6%): 30 genes
are intact in it, 4 are complete but contain a frameshift or stop codon in the assembly
sequence, 5 are partial or fragmented (only the 3′ exons of Bmal1 are present, on an unplaced
scaffold) and 8 are missing or represented only by a divergent copy. Each of these genes is complete in at least one of the two other *E. crispata*
assemblies, and HiFi reads from the sequenced individual cover their exons, so the missing genes
reflect gaps in the assembly and not gene loss. The following genes were not found in any
assembly: Dec1/2, Fbxl3/21, Tyf, Cipc, Ciart, Aanat and Asmt.

**Kleptoplasts.** About 2.4% of the sequenced bases are plastid DNA. Of this, 86% derives from
*Penicillus capitatus* and 14% from a second *Penicillus* species, which matches the only
available reference sequence of *P. lamourouxii*; *Udotea flabellum* (0.09%) and a *Halimeda*
species (0.04%) are present in trace amounts. No read contains algal nuclear ribosomal DNA.

**Gene transfer.** There is no algal DNA in the germline genome of this animal, and no evidence
of insertions in somatic cells above the background of chimeras formed during library
preparation.

- In 4,452 reads, chloroplast DNA is joined to slug nuclear DNA. An insertion in a somatic cell
  would be represented by a single molecule, so each of these reads is a candidate, and the
  absence of the junction from other reads at the same locus is not evidence against it.
  However, the library was amplified by PCR (ultra-low-input protocol), which produces chimeras
  between unrelated molecules: bacterial DNA is joined to nuclear DNA at the same frequency
  (3.24%, compared with 3.27% for chloroplast DNA), and 165 of the chloroplast junctions contain
  the amplification adapter.
- No read contains chloroplast DNA flanked on both sides by contiguous nuclear DNA, which is the
  one structure that the joining of two fragments does not produce. Reads long enough to show
  such a structure cover the genome 68 times, which is equivalent to the DNA of 34 cells.
- No junction is present in two independent molecules. All 4,452 reads were retrieved again and
  compared with one another around the junction, so that junctions within repeats were
  included: 53 junctions occur in more than one read, in every case as PCR duplicates of a
  single molecule. Three pairs of molecules have junctions a few bases apart, at simple repeats.
- xgElyCris1.1 and its alternate haplotype contain 26 fragments of kleptoplast DNA (up to
  5.6 kb long and 98–100% identical to the kleptoplast DNA of this animal). Each ends at an
  assembly gap or at the end of a sequence and is supported by a single molecule, one of the
  chimeric reads described above, whereas other molecules from the same locus continue into
  different sequence. The fragments are not part of the germline genome and should not be
  included in a reference assembly. The Ensembl annotation lists six chloroplast genes on them
  as slug genes. The intervals to mask are listed in `results/plastid_fragments_in_assembly.tsv`.
- Apart from these six, no protein in either annotation is of algal origin. When 471 plastid and
  photosynthesis proteins were aligned directly to all four assemblies, and every plastid-like
  DNA match that is neither an exact copy nor a simple repeat was examined, no algal sequence
  was found within nuclear DNA in xgElyCris1.1, its alternate haplotype or the Tampa assembly.
  The reads of this animal also contain no divergent plastid-like sequence joined to nuclear
  DNA.
- The Purdue assembly (a different individual, sequenced with Nanopore reads) contains 11
  contigs (646 kb) of *Acetabularia* chloroplast DNA, which evidently derive from the
  kleptoplasts of that animal, an 11 kb contig from another alga, and a few positions at which
  chloroplast DNA is joined to nuclear DNA at the end of a contig. These junctions were tested
  with the reads of that individual; see the gene-transfer page of the website.

The website presents the evidence in full, together with the limitations of each analysis.

## Repository layout

| Path | Contents |
|---|---|
| `docs/` | The website (static). `index.html` lists the species, `compare.html` compares their clock proteins, and each species has its own folder (`docs/elysia-crispata/`, …). The data are embedded in `docs/data/<species>.js` and `docs/data/genus.js`; the files offered for download are in `docs/downloads/` and `docs/elysia-crispata/downloads/`. The photographs and their credits are in `docs/assets/photos/`. |
| `results/` | Tables for *E. crispata*: `clock_genes.tsv` / `.json`, `sequences/` (proteins and coding sequences), `kleptoplast.json`, `kleptoplast_barcodes.tsv`, `hgt.json`, `plastid_fragments_in_assembly.tsv`, `penicillus_plastome_genes.bed`. |
| `results/genus/` | For each of the other species: the clock-gene catalogue (`<species>.clock.json`), the composition of its reads (`<species>.kleptoplast.json`) and the gene-transfer analyses (`<species>.transfer.json`). Also the gene-by-species table (`clock_matrix.tsv`), protein alignments (`alignments/`), the concatenated alignment and species tree, and the intermediate tables of the read and assembly analyses: the reference panel (`panel.tsv.gz`), barcode records and sequences, plastid contigs, read junctions, every protein match (`protein_screen_loci.json.gz`) and every coding region on the bacterial and other foreign sequences (`foreign_sequences.json.gz`). |
| `results/trees/` | Trimmed alignments of the gene families and FastTree trees in Newick format. |
| `results/intermediate/` | Small intermediate tables: orthologue calls, choice of gene model, read support, expression, the list of sequences in the reference panel, chimeric-read junctions, boundary evidence and the inputs to the protein-level search. |
| `config/` | Curated inputs: reference clock proteins (`clock_seeds.tsv`), the gene catalogue with its anchors (`clock_catalogue.tsv`), genes searched for but not found (`clock_not_found.tsv`), the species and their assemblies (`species.tsv`), the read sets analysed (`species_reads.tsv`), the text of the conclusions shown on the species pages (`species_text.json`), and the survey of public *Elysia* data (`species_catalogue.tsv`). |
| `scripts/` | The *E. crispata* pipeline, numbered in the order in which it is run; `scripts/genus/` contains the steps that are run for every species. |

## Scripts

| Step | Scripts | Purpose |
|---|---|---|
| Download | `00_fetch_reads.sh`, `00_fetch_proteomes.sh`, `00_fetch_algal_proteomes.sh`, `ncbi_fetch.py` | HiFi runs from the SRA open-data bucket as compressed FASTA; comparison proteomes from UniProt. |
| Clock genes | `01`–`09`, `miniprot_parse.py`, `tree_report.py` | Reference proteins, DIAMOND searches, gene-family alignments and trees, reciprocal best hits, selection of one protein model per gene, status in each of the four assemblies, Pfam domains, tree figures. |
| Reads, first pass | `10`–`12`, `16`, `17` | Mapping of every HiFi read to a panel of plastid, mitochondrial, rDNA and bacterial sequences and to the clock-gene exons; classification and extraction of reads; read support for exons; further analysis of reads with only a weak plastid match. |
| Kleptoplasts | `13` | rbcL and tufA haplotypes extracted from the reads and compared with GenBank barcode sequences. |
| Chimeric reads | `14`, `15` | Base-level realignment of reads with an unexplained segment; junctions between nuclear, plastid and bacterial DNA, with PCR duplicates collapsed. |
| RNA | `20` | Quantification of the head RNA-seq library with salmon. |
| Assemblies | `21`–`25`, `27`–`29` | Test of annotated proteins for algal affinity; plastid DNA in the four assemblies (exact copies, and every other match that is not a simple repeat); read evidence at each boundary between plastid and nuclear DNA; alignment of plastid and photosynthesis proteins to the assemblies; identification of foreign contigs and tests of their junctions with reads. |
| Collect and build | `07`, `26`, `30`, `31` | Result tables and the website. |

Steps that are run for every species (`scripts/genus/`):

| Step | Scripts | Purpose |
|---|---|---|
| Download | `00_fetch_assemblies.sh`, `01_fetch_photos.py`, `02_choose_photos.py` | Assemblies (and annotations, where GenBank has one) from NCBI; openly licensed photographs from iNaturalist, with their credits. |
| Clock genes | `09_clock_refs.sh`, `10_project_clock_genes.sh`, `11_clock_catalogue.py` | Alignment of the 47 *E. crispata* proteins to each genome with miniprot; addition of annotated proteins where available; a reciprocal best-hit test against all *E. crispata* proteins; selection of one model per gene. |
| Comparison | `12_compare_species.py`, `13_domains.sh` | Protein alignments (MAFFT), completeness relative to the alignment positions present in most sequences, pairwise identity, the species tree (IQ-TREE), Pfam domains. |
| Reference panel | `30_build_panel.py`, `30b_mask_and_index.sh` | Plastid, mitochondrial and ribosomal DNA reference sequences for all species; masking of low-complexity sequence for the first pass. |
| Reads | `fetch_ranges.py`, `31c_screen_run.py`, `31_select_hits.py`, `31b_slim_rdna.py`, `31d_sample_blocks.py`, `33_pass2.sh`, `sam_segments.py`, `49_analyse_library.sh`, `34_read_analysis.py`, `37_pair_analysis.py`, `55_end_piece_control.py` | Processing of each run in full, in batches, and mapping to the panel; base-level alignment of the reads with a match, and of a random sample of the library, to the assembly of the species and to the panel; composition of the library; every read or read pair in which two types of DNA are joined, with the same counts for mitochondrial DNA and for the random sample; the type of DNA at the end of the fragment in read pairs that span a junction. |
| Single molecules | `46_candidates.py`, `51_recurrent_junctions.py`, `47_library_control.sh`, `47_shape_reads.py`, `52_fetch_reads.py`, `53_join_in_short_reads.py`, `54_internal_pieces.py` | Classification of long reads by structure (organelle DNA flanked by nuclear DNA, organelle DNA at one end, or the reciprocal structure), with short segments confirmed at base level; comparison of junctions between reads to find any that recur in independent molecules; PCR duplicates; the background rate of chimeras, measured on every read of a small run; retrieval of individual reads by name; search for a known junction in read pairs. |
| Identification of the alga | `36_barcodes.py`, `38_plastid_contigs.sh` | rbcL and tufA extracted from long reads, or from plastid contigs assembled from short reads (MEGAHIT), and compared with GenBank. |
| Assemblies | `35_assembly_screen.py`, `42_protein_screen_inputs.py`, `43_foreign_sequences.py` (with `scripts/22` and `scripts/28`) | Organelle-like DNA in each assembly; alignment of 471 algal proteins to each assembly; translation in full of every sequence with a non-animal match, and classification of its genes. |
| Read-based tests of loci | `40_locus_baits.py`, `41_locus_test.py` | Reads that span, or end at, each boundary of an organelle-like region within an assembly sequence, and at the locus of each read in which organelle DNA is flanked by nuclear DNA. |
| Collect and build | `39_collect_species.py`, `45_export_tables.py`, `chrome.py`, `21_build_species_pages.py`, `20_build_genus_site.py` | Result files and small tables; the shared page frame; the kleptoplast and gene-transfer pages of each species; the home page, comparison page, clock-gene pages and methods. |
| Order of the analyses | `run_screens.sh` | The commands for the six additional species in the order in which they were run (`DRY=1` prints them). |

The scripts were written for a single machine. The locations of the tools
(`/home/claude/tools/mamba/envs/bio/bin`, which can be overridden with `BIO_BIN` in the scripts
that read it) and the working directories are passed on the command line or set at the top of
each script, and the header of each script describes its usage. Software used: minimap2,
miniprot, DIAMOND, BLAST+, MAFFT, trimAl, FastTree, IQ-TREE, HMMER (Pfam-A), salmon, samtools,
seqkit, sra-tools, MEGAHIT, Python 3 with Biopython.

To rebuild the website from the results already in the repository:

```
python3 scripts/30_build_site.py . results/trees/figures   # E. crispata pages; the second argument is optional
python3 scripts/genus/21_build_species_pages.py .          # kleptoplast and gene-transfer pages of the other species
python3 scripts/genus/20_build_genus_site.py .             # home page, comparison, methods and clock-gene pages
python3 -m http.server -d docs 8000
```

## Data

| Data set | Accession |
|---|---|
| Primary assembly xgElyCris1.1, Ensembl annotation 2024_05 | GCA_963854125.1 |
| Alternate haplotype | GCA_963854135.1 |
| PacBio HiFi reads (5,516,672 reads, 46.85 Gb) | ERR11843472, ERR12353011, ERR12356115 |
| RNA-seq, head, same individual | ERR12342469 |
| Plastid genomes assembled from this individual | OZ373742.1, OZ372796.1 |
| Bacterial genomes assembled from this individual | GCA_964289745, GCA_964291485, GCA_964291505, GCA_964291535 |
| Other *E. crispata* assemblies (other individuals) | GCA_033675545.1 (Purdue, annotated), GCA_984581095.1 (Tampa) |
| Nanopore reads of the Purdue individual (used only to test junctions) | PRJNA987316, 12 MinION runs |

Assemblies of the other species: *E. chlorotica* GCA_003991915.1, *E. marginata* GCA_019649035.1,
*E. timida* GCA_043644045.1, *E. leucolegnote* GCA_052815195.1, *E. atroviridis* GCA_059052615.1,
*E. papillosa* GCA_977015215.2 (details in `config/species.tsv`).

Reads analysed for the other species (details in `config/species_reads.tsv`):

| Species | Runs | Platform and source of DNA | Bases analysed |
|---|---|---|---|
| *E. papillosa* | ERR14842698, ERR14842697 | PacBio Revio (HiFi); assembled individual | 36.2 Gb |
| *E. papillosa* | ERR14842696 | PacBio Revio (HiFi, amplified library); second individual | 64.8 Gb |
| *E. atroviridis* | SRR34916404 | PacBio Sequel II (HiFi); assembled individual | 30.6 Gb |
| *E. timida* | SRR30016879 | PacBio Revio (HiFi, amplified ultra-low-input library); assembled individuals | 53.1 Gb |
| *E. chlorotica* | SRR7647856 | PacBio RS (older, error-prone long reads); an adult | 9.4 Gb |
| *E. chlorotica* | SRR7647851, SRR7647852, SRR7647853, SRR7647854, SRR7647849, SRR7647850 | Illumina HiSeq 2000 (paired-end; inserts of 170, 500 and 800 bp); the adult that provided the short reads for the assembly | 102.5 Gb |
| *E. chlorotica* | SRR768897, SRR768899 | Illumina Genome Analyzer IIx (paired-end); eggs | 22.4 Gb |
| *E. marginata* | DRR238951, DRR238952 | Illumina HiSeq 2500 (paired-end; inserts of 250 and 500 bp); assembled individual | 84.4 Gb |
| *E. leucolegnote* | SRR41369319 | Illumina NovaSeq X Plus (paired-end); an individual other than the one assembled | 12.5 Gb |
| *E. leucolegnote* | SRR17036819 | Illumina HiSeq 2500 (paired-end); an individual other than the one assembled | 25.5 Gb |

Every run was analysed in full.

The reference panel used for the *E. crispata* reads (241 green-algal plastid genomes, among
other sequences) is listed in `results/intermediate/panel.tsv`; the larger panel used for the
other species (1,271 plastid sequences, 117 algal mitochondrial genomes, 7,352 algal ribosomal
DNA sequences and nine *Elysia* mitochondrial genomes) is listed in `results/genus/panel.tsv.gz`.

## Limitations

- xgElyCris1.1 represents about 60% of the k-mers in the reads, so the absence of a sequence
  from this assembly is weak evidence that it is absent from the genome.
- Six protein models were derived by aligning a protein from another species to the genome and
  may lack lineage-specific exons.
- Support values in the gene trees are FastTree local support values. The identification of
  eleven genes rests mainly on reciprocal best hits.
- Expression values are from a single head sample collected at an unknown time of day.
- The boundary test and the comparison of the Purdue contigs with the xgElyCris1 reads used two
  of the three HiFi runs (53% of the data).
- The Purdue and Tampa assemblies were examined only as far as was needed to interpret
  xgElyCris1; the Purdue assembly very probably contains more *Acetabularia* chloroplast contigs
  than those found through their genes.
- The identification of the second kleptoplast species rests on a single 562-bp reference
  sequence.
- Other species: the *E. leucolegnote* reads are from two individuals other than the one
  assembled, and algae without a close reference sequence are under-represented in them; the
  plastid share given for *E. timida* is a minimum. Every run was analysed in full, but in the
  second *E. chlorotica* adult, in which a ninth of the DNA is chloroplast DNA, junctions were
  analysed in one in 50 of the read pairs that match an organelle genome. The protein-level
  search uses 471 query proteins, so it detects only those bacterial and other foreign
  sequences that encode a homologue of one of them.
- Insertions in somatic cells: sequencing a whole animal samples, at each locus, DNA equivalent
  to 1 to 34 cells. Insertions present in a small proportion of cells therefore cannot be
  detected in any of these libraries, and a single molecule in which chloroplast DNA is joined
  to nuclear DNA cannot be distinguished from a chimera formed during library preparation. The
  detection limit of each library is given on the corresponding page of the website.

## Credits

The genome data for *E. crispata* and *E. papillosa* were produced by the Wellcome Sanger
Institute Tree of Life programme and the Aquatic Symbiosis Genomics project. Please cite
Middlebrooks et al. 2026 when using the *E. crispata* data, and the publications describing the
other genomes (Cai et al. 2019; Maeda et al. 2021; Guo & Kong 2026; and the submitters named in
each GenBank record) when using those. The photographs are from iNaturalist and are used under
Creative Commons licences; each is credited on the site and in `docs/assets/photos/credits.json`
(two are CC BY-NC, the others CC BY). The analysis and the website were produced with Claude
(Anthropic) for the owner of the repository. The fonts in `docs/assets/fonts/` are Bricolage
Grotesque and IBM Plex Mono, both under the SIL Open Font License 1.1.
