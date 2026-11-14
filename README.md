# Elysia clocks

Circadian clock genes, kleptoplast identity and tests for algal gene transfer in sea slugs of the
genus *Elysia*, worked out from public genome data. The repository began with the lettuce sea
slug *Elysia crispata* (Sanger specimen xgElyCris1, BioProject
[PRJEB67629](https://www.ebi.ac.uk/ena/browser/view/PRJEB67629); Middlebrooks et al. 2026,
[Wellcome Open Research 11:148](https://wellcomeopenresearch.org/articles/11-148)) and now covers
every *Elysia* species with a nuclear genome assembly.

The results are presented as a small website in `docs/`, served by GitHub Pages:
**https://comparativechrono.github.io/elysia-crispata-clock/**

Nothing here measures a rhythm. Gene roles are inferred from homology with flies and mammals,
and the work has not been reviewed by the authors of the genomes.

## The genus

Seven species have a nuclear genome. The 47 clock genes catalogued in *E. crispata* were located
in each one (`scripts/genus/`), and their proteins aligned across species.

| Species | Reference assembly | Intact | Whole, sequence does not translate cleanly | Partial or in pieces | Not found |
|---|---|---|---|---|---|
| *Elysia crispata* | GCA_963854125.1 (chromosome) | 30 | 4 | 5 | 8 |
| *Elysia chlorotica* | GCA_003991915.1 (scaffold) | 39 | 5 | 3 | 0 |
| *Elysia marginata* | GCA_019649035.1 (scaffold) | 33 | 8 | 4 | 2 |
| *Elysia timida* | GCA_043644045.1 (chromosome) | 44 | 1 | 1 | 1 |
| *Elysia leucolegnote* | GCA_052815195.1 (chromosome) | 43 | 1 | 2 | 1 |
| *Elysia atroviridis* | GCA_059052615.1 (chromosome) | 44 | 1 | 1 | 1 |
| *Elysia papillosa* | GCA_977015215.2 (contig) | 43 | 0 | 3 | 1 |

- 27 of the 47 genes are whole in the reference assembly of all seven species, and every species
  has the whole core loop. Where a gene is partial or missing, the assembly is usually the reason:
  the two assemblies built mainly from short reads have the most such genes.
- None of the genes that could not be found in *E. crispata* (Dec1/2, Fbxl3/21, Tyf, Cipc, Ciart,
  Aanat, Asmt) turns up in any other species.
- The number of Go-opsin genes differs between species (one to three); Go-opsin 1 of
  *E. crispata* has a one-to-one counterpart only in *E. chlorotica*.
- Of the 54 pairs of clock genes that share a chromosome in *E. crispata*, 53 also share one in
  *E. timida*, *E. leucolegnote* and *E. atroviridis*.
- A tree from 35 clock proteins (21,190 aligned positions) groups *E. crispata* with
  *E. chlorotica*, and *E. timida* with *E. leucolegnote* and *E. atroviridis*; *E. papillosa*
  branches first. All nodes have full support, but a tree from clock proteins alone is a guide,
  not a phylogeny of the genus.

Models for species without a public annotation (and for genes the annotation lacks) are the
*E. crispata* protein aligned to the genome. They can miss exons that *E. crispata* does not have.
The gene annotation of *E. timida* is not in GenBank and its server could not be reached, so all
*E. timida* models are of this kind. The GenBank sequence of *E. leucolegnote* is 35% masked (N).

**Whose chloroplasts.** The plastid DNA among the reads of each sequenced animal names the alga it
had been eating. Every read of every library was examined; the read sets are listed in
`config/species_reads.tsv`.

| Species | Chloroplasts from | Reads from | Plastid share of sequenced bases |
|---|---|---|---|
| *E. chlorotica* | *Vaucheria litorea* | two adults, and eggs | 3.5% and 11.0%; eggs: 13 read pairs in 77 million |
| *E. marginata* | *Bryopsis* (*B. hypnoides* type, PQ652170.1) | the animal assembled | 0.26% |
| *E. timida* | *Acetabularia acetabulum* (rbcL and tufA match its GenBank records) | the animals assembled | 0.0054% or more: only 49 kb of that plastid genome is public |
| *E. leucolegnote* | a *Boodleopsis* (97 to 99% identical to NC_039521.1) in both animals; an alga without a close reference at a quarter of its amount in the Shenzhen animal and about half in the Hainan animal; in the Hainan animal a third plastid genome, also without a reference, as plentiful as the *Boodleopsis* | two other animals | 0.28% and 0.074% or more |
| *E. atroviridis* | *Bryopsis*, the same type, as two strains | the animal assembled | 1.9% |
| *E. papillosa* | *Penicillus lamourouxii* (OZ494339.1), in both animals | the animal assembled and a second one | 0.12% and 0.014% |

**Gene transfer.** None of these animals inherited algal DNA: no assembly holds algal DNA that its
reads bear out, and egg DNA of *E. chlorotica* holds none. But chloroplasts sit in body cells, so an
insertion would be carried by one cell or a few, and in the reads by single molecules; other reads
from the same place that lack it do not count against it. Every read of every library was therefore
searched for single molecules that join chloroplast DNA to slug DNA. They are listed, one molecule a
row, in `results/genus/` and `docs/downloads/` (`*junctions.tsv`, `*junction_pairs.tsv.gz`,
`*insertion_shapes.tsv`).

| Library | Chloroplast reads (or pairs) joined to slug DNA | Mitochondrial ones | Slug DNA joined to unrelated slug DNA | Reads with chloroplast DNA inside slug DNA |
|---|---|---|---|---|
| *E. crispata*, accurate long reads (amplified) | 4,452 of 136,293 (3.3%) | bacterial DNA instead: 3.2% | 4.2% | 0 |
| *E. chlorotica*, adult, old long reads | 2,397 of 43,514 (5.5%) | 113 of 2,116 (5.3%) | 5.3% | 9 |
| *E. chlorotica*, second adult, read pairs (estimated from one pair in 50) | 71,401 of 57,986,168 (0.12%) | 976 of 1,113,806 (0.088%) | 0.99% | read pairs |
| *E. chlorotica*, eggs, read pairs | 0 of 13 (0%) | 1 of 104,147 (0.001%) | 0.91% | read pairs |
| *E. atroviridis*, accurate long reads | 14 of 69,026 (0.020%) | 0 of 2,125 (0%) | 2.3% | 0 |
| *E. papillosa*, assembled animal, accurate long reads | 6 of 6,533 (0.092%) | 4 of 2,012 (0.20%) | 2.9% | 0 |
| *E. papillosa*, second animal, accurate long reads (amplified) | 3 of 1,326 (0.23%) | 29 of 7,954 (0.36%) | 5.3% | 0 |
| *E. timida*, accurate long reads (amplified) | 15 of 2,598 (0.58%) | 75 of 14,113 (0.53%) | 1.1% | 0 |
| *E. marginata*, read pairs | 78 of 759,795 (0.010%) | 1 of 35,000 (0.0029%) | 1.7% | read pairs |
| *E. leucolegnote*, Shenzhen animal, read pairs | 1,667 of 125,085 (1.3%) | 593 of 66,053 (0.90%) | 1.2% | read pairs |
| *E. leucolegnote*, Hainan animal, read pairs | 491 of 68,603 (0.72%) | 166 of 32,441 (0.51%) | 0.83% | read pairs |

- Every library holds such molecules, and every library joins DNA that has no business being
  joined (the slug's mitochondrial DNA to its chromosomes, one chromosome to another). Preparing
  DNA for sequencing makes such joins, and a real junction read once looks the same. In the
  amplified libraries the amplification adapter lies in many of the joins.
- Rates. In the long-read libraries chloroplast DNA is joined to slug DNA as often as
  mitochondrial DNA is. In the read-pair libraries with enough of both to compare it is joined
  1.3 to 1.5 times as often (pairs at the bacterial ribosomal RNA genes of *E. leucolegnote* left
  out). Insertions in body cells would raise the rate in that way, but the mates that read
  through a join point to the tube. Such a mate begins with the short piece at the end of its
  fragment. That piece is chloroplast DNA in 98 of 170 chloroplast joins and mitochondrial DNA in
  only 14 of 51 mitochondrial ones, where cutting at random out of a chromosome would give the
  same share for both (two controls, `scripts/genus/55_end_piece_control.py`, put it near 45%).
  The two kinds of DNA come in different pieces in those libraries, so joining in the tube does
  not treat them alike, and the difference in rate cannot be read as insertions. Had all of it
  been insertions, they would have numbered 10 to 150 in every cell of those animals.
- Chloroplast DNA inside slug DNA, with the slug sequence carrying on where it left off, is the
  one shape that joining end to end does not give. 9 reads have it, all in the old,
  error-prone long reads of *E. chlorotica*, along with 2 that hold mitochondrial DNA in the
  same way and 2 of the mirror shape (slug DNA inside chloroplast DNA). Each is a single
  molecule at a place of its own. The accurate long reads of the other four species hold none.
- No join is carried by two molecules in any long-read library; the reads were compared with each
  other around every join, so that joins into repeats count too. Pairs of molecules whose joins lie
  a few bases apart occur in *E. chlorotica* (two pairs, at a satellite repeat; one of them is the
  join built into the assembly) and in *E. crispata* (three pairs, at simple repeats).
- In the read-pair libraries, copies of one fragment are recognised by the reads themselves, and
  two molecules count as carrying one join when both their slug sides and their organelle sides
  lie within 1 kb of each other. Such joins are found at four places: three that are not organelle DNA (in both
  animals of *E. leucolegnote* a transfer RNA gene of the slug at a gap of the assembly and the
  ribosomal RNA genes of a bacterial sequence; in *E. marginata* a ribosomal RNA gene at a gap in
  a bacterial sequence), and one in the egg library of *E. chlorotica*: 88 bases of mitochondrial
  DNA between two neighbouring bases of RQTK01000018.1, read through by two molecules and absent
  from the assembly. That is an inherited insertion of organelle DNA, of the mitochondrial kind,
  and it shows that the test finds one. No such join involves chloroplast DNA.
- What the libraries could show: reads long enough to hold an insertion with both flanks cover
  each genome 29 to 68 times, which is 15 to 34 cells' worth of DNA. About three
  insertion-shaped reads would be expected with one insertion in every 5 to 11 cells.
  Rates like those measured in tobacco (about one pollen grain in 16,000, one leaf cell in five
  million) are far out of reach of sequencing a whole animal. Among the read-pair libraries
  *E. marginata* gives a limit of its own: it joins chloroplast DNA to slug DNA so seldom
  (78 pairs in the library) that this is about what one insertion in every cell would add.

In the assemblies:

- *E. chlorotica*: five small plastid scaffolds (8.3 kb), plastid stretches inside two larger
  scaffolds (1.6 kb and 230 bp) and four stretches of mitochondrial DNA. The 1.6 kb stretch, at
  the start of RQTK01001031.1, is joined to a satellite repeat of the slug; one long read carries
  that join and a second molecule has a join 47 bases away. No short read of the second adult or
  of the eggs carries it. No read of any library crosses an edge of the 230 bp stretch, and two
  long reads span the place without it. No short read crosses an edge of the mitochondrial
  stretches. None of the *Vaucheria* nuclear genes reported in early studies (psbO,
  prk, lhc, chlD, chlG, chlH) matches the assembly. Fourteen scaffolds (525 kb) come from a
  bacterium related to *Francisella*.
- *E. marginata*: one 3.7 kb scaffold is *Bryopsis* plastid DNA. About 2.0 Mb on 23 scaffolds is
  bacterial (mycoplasma-like, flavobacteria-like and others), and the gene set published with the
  assembly includes genes on them.
- *E. leucolegnote*: the unplaced sequence JBQKML010000381.1 begins with 863 bp of plastid DNA
  joined to a slug repeat. In two other animals the same plastid sequence carries on as plastid
  DNA at that point (87 and 77 reads stop, none runs through). The reads of the assembled
  animal are not public, so the join cannot be tested. Four sequences (1.9 Mb) are bacterial.
- *E. papillosa*: contig CDRNNR020002991.1 (30.9 kb) is a piece of an algal mitochondrial genome.
- *E. timida*, *E. atroviridis*: nothing of the kind in the assembly.

In the protein screen of these species the gene set published with an assembly is not counted as
evidence that a sequence of that assembly is animal, because the published gene sets of
*E. chlorotica* and *E. marginata* include genes on bacterial and plastid sequences.

## *Elysia crispata* in detail

**Clock genes.** 47 genes are catalogued: the core loop (Clock, Bmal1, Period, Timeless, Timeout,
Cry1, Cry2), the second loop (Pdp1, Vrille, E75, Hr3, Cwo), kinases, phosphatases, degradation
factors, output genes, photolyases and nine opsins. All 47 are present in the species, but the
Sanger assembly xgElyCris1.1 is incomplete (BUSCO 73.6%): 30 are intact in it, 4 are whole with a
frameshift or mismatch, 5 are partial or in pieces (Bmal1 has only its 3′ exons, on an unplaced
scaffold) and 8 are missing or represented only by a diverged copy. Every one is whole in at least
one of the two other *E. crispata* assemblies, and HiFi reads from the sequenced animal cover their
exons, so the holes are assembly gaps, not gene losses. Not found in any assembly: Dec1/2,
Fbxl3/21, Tyf, Cipc, Ciart, Aanat, Asmt.

**Kleptoplasts.** About 2.4% of the sequenced bases are plastid DNA. 86% of it is *Penicillus
capitatus* and 14% a second *Penicillus* that matches the one available *P. lamourouxii*
reference; *Udotea flabellum* (0.09%) and a *Halimeda* (0.04%) are present as traces. No read
carries algal nuclear ribosomal DNA.

**Gene transfer.** Not in the genome this animal inherited. In single cells: nothing that stands
above what the library makes by itself.

- 4,452 reads join chloroplast DNA to slug DNA. An insertion in a body cell would be carried by
  a single molecule, so each of these reads is a candidate, and reads from the same place that
  lack the join do not count against it. But the library was amplified by PCR (ultra-low-input
  protocol), which joins unrelated molecules: bacterial DNA is joined to slug DNA at the same
  rate (3.24%, against 3.27%), and 165 of the chloroplast joins hold the amplification adapter.
- No read holds chloroplast DNA inside slug DNA, the one shape that joining end to end does not
  give. Reads long enough to show it cover the genome 68 times: 34 cells' worth of DNA.
- No join is carried by two molecules. All 4,452 reads were fetched again and compared with each
  other around the join, so that joins into repeats count too: 53 joins occur in more than one
  read, each time as copies of one molecule. Three pairs of molecules have joins a few bases
  apart, at simple repeats.
- xgElyCris1.1 and its alternate haplotype contain 26 stretches of kleptoplast DNA (up to 5.6 kb,
  98–100% identical to the animal's kleptoplasts). Every one ends at an assembly gap or sequence
  end and rests on a single molecule, one of the joins above, while other molecules from the
  same place carry on into different sequence. They are not part of the inherited genome and do
  not belong in a reference. The Ensembl annotation lists six chloroplast genes on them as slug
  genes. Intervals to mask: `results/plastid_fragments_in_assembly.tsv`.
- No annotated protein of either annotation is algal apart from those six. Aligning 471 plastid
  and photosynthesis proteins directly to all four assemblies, and listing every plastid-like
  DNA match that is neither an exact copy nor a simple repeat, finds no algal sequence among slug
  DNA in xgElyCris1.1, its alternate haplotype or the Tampa assembly. The reads of the animal
  contain no diverged plastid-like sequence joined to slug DNA either.
- The Purdue assembly (a different animal, Nanopore reads) contains 11 contigs (646 kb) of
  *Acetabularia* chloroplast DNA, evidently from that animal's kleptoplasts, an 11 kb contig from
  a further alga, and a few places where chloroplast DNA meets slug DNA at the end of a contig.
  These were tested with that animal's own reads; see the gene-transfer page of the website.

The website gives the evidence, the caveats and what the screens cannot rule out.

## Layout

| Path | Contents |
|---|---|
| `docs/` | The website (static). `index.html` lists the species, `compare.html` compares their clock proteins, and each species has a folder (`docs/elysia-crispata/`, …). Data are embedded in `docs/data/<species>.js` and `docs/data/genus.js`; files offered for download are in `docs/downloads/` and `docs/elysia-crispata/downloads/`. Photographs and their credits are in `docs/assets/photos/`. |
| `results/` | *E. crispata* tables: `clock_genes.tsv` / `.json`, `sequences/` (proteins and coding sequences), `kleptoplast.json`, `kleptoplast_barcodes.tsv`, `hgt.json`, `plastid_fragments_in_assembly.tsv`, `penicillus_plastome_genes.bed`. |
| `results/genus/` | Per further species: the clock catalogue (`<species>.clock.json`), what its reads contain (`<species>.kleptoplast.json`) and the gene-transfer screens (`<species>.transfer.json`). Also the gene-by-species table (`clock_matrix.tsv`), protein alignments (`alignments/`), the concatenated alignment and species tree, and the intermediate tables of the screens: the reference panel (`panel.tsv.gz`), barcode records and sequences, plastid contigs, read junctions, every protein match (`protein_screen_loci.json.gz`) and every coding stretch on the bacterial and other foreign sequences (`foreign_sequences.json.gz`). |
| `results/trees/` | Family alignments (trimmed) and FastTree trees in Newick format. |
| `results/intermediate/` | Small intermediate tables: orthologue calls, model choice, read support, expression, the reference panel listing, chimeric-read junctions, boundary evidence, protein-screen inputs. |
| `config/` | The curated inputs: reference clock proteins (`clock_seeds.tsv`), the gene catalogue with its anchors (`clock_catalogue.tsv`), genes searched for but not found (`clock_not_found.tsv`), the species and their assemblies (`species.tsv`), the read sets screened (`species_reads.tsv`), the written conclusions shown on the species pages (`species_text.json`), and the survey of public *Elysia* data (`species_catalogue.tsv`). |
| `scripts/` | The *E. crispata* pipeline, numbered in the order it runs; `scripts/genus/` holds the steps that run over every species. |

## Scripts

| Step | Scripts | What they do |
|---|---|---|
| Fetch | `00_fetch_reads.sh`, `00_fetch_proteomes.sh`, `00_fetch_algal_proteomes.sh`, `ncbi_fetch.py` | HiFi runs from the SRA open-data bucket as compressed FASTA; comparison proteomes from UniProt. |
| Clock genes | `01`–`09`, `miniprot_parse.py`, `tree_report.py` | Reference proteins, DIAMOND searches, family alignments and trees, reciprocal best hits, choice of a protein model per gene, presence in each of four assemblies, Pfam domains, tree figures. |
| Reads, pass 1 | `10`–`12`, `16`, `17` | Every HiFi read mapped to a panel of plastid, mitochondrial, rDNA and bacterial sequences plus clock-gene exons; reads classified and extracted; exon read support; a closer look at reads with only a weak plastid match. |
| Kleptoplasts | `13` | rbcL and tufA haplotypes read out of the reads and compared with GenBank barcodes. |
| Chimeric reads | `14`, `15` | Base-level re-alignment of reads with an unexplained part; junctions between slug, plastid and bacterial DNA, with PCR duplicates collapsed. |
| RNA | `20` | salmon quantification of the head RNA-seq library. |
| Assemblies | `21`–`25`, `27`–`29` | Annotated proteins tested for algal affinity; plastid DNA in the four assemblies (exact copies, and every other match that is not a simple repeat); read evidence at each plastid/slug boundary; plastid and photosynthesis proteins aligned to the assemblies; foreign contigs named and their joins tested with reads. |
| Collect and build | `07`, `26`, `30`, `31` | Results tables and the website. |

Steps that run over every species (`scripts/genus/`):

| Step | Scripts | What they do |
|---|---|---|
| Fetch | `00_fetch_assemblies.sh`, `01_fetch_photos.py`, `02_choose_photos.py` | Assemblies (and annotation where GenBank has one) from NCBI; openly licensed photographs from iNaturalist with their credits. |
| Clock genes | `09_clock_refs.sh`, `10_project_clock_genes.sh`, `11_clock_catalogue.py` | The 47 *E. crispata* proteins aligned to each genome with miniprot; annotated proteins added where available; a reciprocal test against all *E. crispata* proteins; one model per gene. |
| Compare | `12_compare_species.py`, `13_domains.sh` | Protein alignments (MAFFT), completeness against the positions most sequences share, pairwise identity, the species tree (IQ-TREE), Pfam domains. |
| Reference panel | `30_build_panel.py`, `30b_mask_and_index.sh` | Plastid, mitochondrial and ribosomal DNA references for all species; low-complexity sequence masked for the first pass. |
| Reads | `fetch_ranges.py`, `31c_screen_run.py`, `31_select_hits.py`, `31b_slim_rdna.py`, `31d_sample_blocks.py`, `33_pass2.sh`, `sam_segments.py`, `49_analyse_library.sh`, `34_read_analysis.py`, `37_pair_analysis.py`, `55_end_piece_control.py` | Every run read in full, in pieces, and mapped to the panel; reads with a match, and a random sample of the library, aligned base by base to the species' assembly and the panel; what the library contains; every read or read pair that joins two kinds of DNA, with the same counts for mitochondrial DNA and for the sample; which kind of DNA lies at the end of the fragment in read pairs that cross a join. |
| Single molecules | `46_candidates.py`, `51_recurrent_junctions.py`, `47_library_control.sh`, `47_shape_reads.py`, `52_fetch_reads.py`, `53_join_in_short_reads.py`, `54_internal_pieces.py` | Long reads sorted by shape (organelle DNA inside slug DNA, at one end, mirror image), short pieces confirmed base by base; joins compared between reads to find one carried by two molecules; copies of molecules; what a library makes by itself, on every read of a small run; named reads fetched again; a known join looked for in read pairs. |
| Naming the alga | `36_barcodes.py`, `38_plastid_contigs.sh` | rbcL and tufA read out of long reads, or cut from plastid contigs assembled from short reads (MEGAHIT), and compared with GenBank. |
| Assemblies | `35_assembly_screen.py`, `42_protein_screen_inputs.py`, `43_foreign_sequences.py` (with `scripts/22` and `scripts/28`) | Organelle-like DNA in each assembly; 471 algal proteins aligned to it; every sequence with a non-animal match translated as a whole and its genes classified. |
| Places tested with reads | `40_locus_baits.py`, `41_locus_test.py` | Reads that run through, or stop at, each edge of an organelle-like stretch inside an assembly sequence, and at the place of each insertion-shaped read. |
| Collect and build | `39_collect_species.py`, `45_export_tables.py`, `chrome.py`, `21_build_species_pages.py`, `20_build_genus_site.py` | Result files and small tables; shared page frame; the chloroplast and gene-transfer pages of each species; homepage, comparison page, clock pages and methods. |
| Order of the screens | `run_screens.sh` | The commands for the six added species in the order they were run (`DRY=1` prints them). |

The scripts were written for one working machine: tool locations
(`/home/claude/tools/mamba/envs/bio/bin`, overridable with `BIO_BIN` where the script reads it)
and working directories are passed on the command line or set at the top of each script. Each
script states its usage in its header. Software used: minimap2, miniprot, DIAMOND, BLAST+, MAFFT,
trimAl, FastTree, IQ-TREE, HMMER (Pfam-A), salmon, samtools, seqkit, sra-tools, MEGAHIT, Python 3 with Biopython.

To rebuild the website from the results already in the repository:

```
python3 scripts/30_build_site.py . results/trees/figures   # E. crispata pages; second argument optional
python3 scripts/genus/21_build_species_pages.py .          # chloroplast and gene-transfer pages of the other species
python3 scripts/genus/20_build_genus_site.py .             # homepage, comparison, methods and clock pages
python3 -m http.server -d docs 8000
```

## Data

| What | Accession |
|---|---|
| Primary assembly xgElyCris1.1, Ensembl annotation 2024_05 | GCA_963854125.1 |
| Alternate haplotype | GCA_963854135.1 |
| PacBio HiFi reads (5,516,672 reads, 46.85 Gb) | ERR11843472, ERR12353011, ERR12356115 |
| RNA-seq, head, same animal | ERR12342469 |
| Plastid genomes assembled from this animal | OZ373742.1, OZ372796.1 |
| Bacterial genomes assembled from this animal | GCA_964289745, GCA_964291485, GCA_964291505, GCA_964291535 |
| Other *E. crispata* assemblies (other animals) | GCA_033675545.1 (Purdue, annotated), GCA_984581095.1 (Tampa) |
| Nanopore reads of the Purdue animal (join tests only) | PRJNA987316, 12 MinION runs |

Assemblies of the other species: *E. chlorotica* GCA_003991915.1, *E. marginata* GCA_019649035.1,
*E. timida* GCA_043644045.1, *E. leucolegnote* GCA_052815195.1, *E. atroviridis* GCA_059052615.1,
*E. papillosa* GCA_977015215.2 (details in `config/species.tsv`).

Reads screened for the other species (details in `config/species_reads.tsv`):

| Species | Runs | Kind | Examined |
|---|---|---|---|
| *E. papillosa* | ERR14842698, ERR14842697 | PacBio Revio (HiFi); the sequenced animal | 36.2 Gb |
| *E. papillosa* | ERR14842696 | PacBio Revio (HiFi, amplified library); a second animal | 64.8 Gb |
| *E. atroviridis* | SRR34916404 | PacBio Sequel II (HiFi); the sequenced animal | 30.6 Gb |
| *E. timida* | SRR30016879 | PacBio Revio (HiFi, amplified ultra-low-input library); the sequenced animals | 53.1 Gb |
| *E. chlorotica* | SRR7647856 | PacBio RS (older, error-prone long reads); an adult | 9.4 Gb |
| *E. chlorotica* | SRR7647851, SRR7647852, SRR7647853, SRR7647854, SRR7647849, SRR7647850 | Illumina HiSeq 2000 (paired; inserts of 170, 500 and 800 bp); an adult used for the short reads of the assembly | 102.5 Gb |
| *E. chlorotica* | SRR768897, SRR768899 | Illumina Genome Analyzer IIx (paired); eggs | 22.4 Gb |
| *E. marginata* | DRR238951, DRR238952 | Illumina HiSeq 2500 (paired; inserts of 250 and 500 bp); the sequenced animal | 84.4 Gb |
| *E. leucolegnote* | SRR41369319 | Illumina NovaSeq X Plus (paired); another animal | 12.5 Gb |
| *E. leucolegnote* | SRR17036819 | Illumina HiSeq 2500 (paired); another animal | 25.5 Gb |

Every run was read in full.

The reference panel used to screen the *E. crispata* reads (241 green-algal plastid genomes and
others) is listed in `results/intermediate/panel.tsv`; the wider panel used for the other species
(1,271 plastid sequences, 117 algal mitochondrial genomes, 7,352 algal ribosomal DNA sequences
and nine *Elysia* mitochondrial genomes) is listed in `results/genus/panel.tsv.gz`.

## Caveats

- xgElyCris1.1 represents about 60% of the read k-mers; absence from it means little.
- Six protein models are translations guided by another species' protein and may miss
  lineage-specific exons.
- Tree supports are FastTree local supports. Eleven gene names rest mainly on reciprocal best hits.
- Expression values come from one head sample at one unknown time of day.
- The boundary test and the check of the Purdue contigs against xgElyCris1 reads used two of the
  three HiFi runs (53% of the data).
- The Purdue and Tampa assemblies were examined only as far as needed to interpret xgElyCris1;
  the Purdue assembly very likely holds more *Acetabularia* chloroplast contigs than the ones
  found through their genes.
- The second kleptoplast species is named from a single 562-bp reference sequence.
- Other species: *E. leucolegnote* was read from two animals that are not the one assembled, and
  algae without a close reference are under-counted there; the *E. timida* plastid share is a
  floor. Every run was read in full, but in the second *E. chlorotica* adult, where a ninth of
  the DNA is chloroplast DNA, the joins were analysed on one pair in 50 of those that match an
  organelle genome. The protein screen searches with 471 proteins, so it finds only the
  bacterial and foreign sequences that happen to carry a relative of one of them.
- Insertions in body cells: sequencing a whole animal reads 1 to 34 cells' worth of DNA at each
  place, so insertions carried by a small share of cells cannot be seen in any of these
  libraries, and a single molecule that joins chloroplast DNA to slug DNA cannot be told from
  one made when the library was prepared. The pages give, for each library, what it could have
  shown.

## Credits

Genome data for *E. crispata* and *E. papillosa*: Wellcome Sanger Institute Tree of Life and the
Aquatic Symbiosis Genomics project; please cite Middlebrooks et al. 2026 when using the
*E. crispata* data, and the papers of the other genomes (Cai et al. 2019; Maeda et al. 2021;
Guo & Kong 2026; and the submitters named in each GenBank record) when using theirs.
Photographs are from iNaturalist under Creative Commons licences; each is credited on the site
and in `docs/assets/photos/credits.json` (two are CC BY-NC, the rest CC BY). Analysis and website were produced with
Claude (Anthropic) for the repository owner. Fonts in `docs/assets/fonts/`: Bricolage Grotesque
and IBM Plex Mono, both under the SIL Open Font License 1.1.
