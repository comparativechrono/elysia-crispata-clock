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
had been eating. The read sets are listed in `config/species_reads.tsv`.

| Species | Chloroplasts from | Reads from | Plastid share of sequenced bases |
|---|---|---|---|
| *E. chlorotica* | *Vaucheria litorea* | two adults, and eggs as a control | 3.5% and 14.1%; eggs: two reads in 76 million |
| *E. marginata* | *Bryopsis* (*B. hypnoides* type, PQ652170.1) | the animal assembled | 0.27% |
| *E. timida* | *Acetabularia acetabulum* | the animals assembled | 0.003% or more: only 49 kb of that plastid genome is public |
| *E. leucolegnote* | a *Boodleopsis* (97–99% identical to NC_039521.1) and at least one further alga | two other animals | 0.095% and 0.075% or more |
| *E. atroviridis* | *Bryopsis*, the same type, as two strains (95:5) | the animal assembled | 1.9% |
| *E. papillosa* | *Penicillus lamourouxii* (OZ494339.1) | the animal assembled | 0.14% |

**Gene transfer.** No algal DNA inserted in a slug genome was found in any species. Reads that
join plastid DNA to slug DNA are rare, scattered over the genome, and turn up at rates similar to
joins of the slug's own mitochondrial DNA to its nuclear DNA, which library preparation makes.

- *E. chlorotica*: the assembly holds five small plastid scaffolds (8.3 kb), plastid stretches
  inside two larger scaffolds (1.6 kb and 230 bp) and four stretches of mitochondrial DNA. No short
  read from an adult or from eggs supports any of them; three long reads from two molecules cross
  one edge, against 1,772 that stop there. Egg DNA holds almost no plastid DNA. None of the
  *Vaucheria* nuclear genes reported in early studies (psbO, prk, lhc, chlD, chlG, chlH) matches
  the assembly. Fourteen scaffolds (525 kb) come from a bacterium related to *Francisella*.
- *E. marginata*: one 3.7 kb scaffold is *Bryopsis* plastid DNA. About 2.0 Mb on 23 scaffolds is
  bacterial (mycoplasma-like, flavobacteria-like and others), and the gene set published with the
  assembly includes genes on them.
- *E. leucolegnote*: the unplaced sequence JBQKML010000381.1 begins with 863 bp of plastid DNA
  joined to a slug repeat. In two other animals the same plastid sequence carries on as plastid
  DNA at that point (46 and 13 reads stop, none runs through). Four sequences (1.9 Mb) are bacterial.
  Reads of the assembled animal are not public.
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

**Gene transfer.** No algal sequence inserted in the slug genome was found.

- 4,452 reads join chloroplast DNA to slug DNA, but at the same rate (3.3%) as reads that join
  bacterial DNA to slug DNA, and no chloroplast–slug junction recurs in two independent
  molecules. The library was PCR-amplified (ultra-low-input protocol), which creates such chimeras.
- xgElyCris1.1 and its alternate haplotype contain 26 fragments of kleptoplast DNA (up to 5.6 kb,
  98–100% identical to the animal's kleptoplasts). Every one ends at an assembly gap or sequence
  end and its one boundary with slug DNA is carried by a single molecule, while other molecules
  from the same locus continue into different sequence. They are assembly artefacts. The Ensembl
  annotation lists six chloroplast genes on them as slug genes. Intervals to mask:
  `results/plastid_fragments_in_assembly.tsv`.
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
| Reads | `fetch_ranges.py`, `31_read_screen.sh`, `31_select_hits.py`, `31b_slim_rdna.py`, `33_pass2.sh`, `sam_segments.py`, `34_read_analysis.py`, `37_pair_analysis.py` | Each run fetched and mapped to the panel as a stream; reads with a match aligned base by base to the species' assembly and the panel; what the library contains; reads or read pairs that join two kinds of DNA. |
| Naming the alga | `36_barcodes.py`, `38_plastid_contigs.sh` | rbcL and tufA read out of long reads, or cut from plastid contigs assembled from short reads (MEGAHIT), and compared with GenBank. |
| Assemblies | `35_assembly_screen.py`, `42_protein_screen_inputs.py`, `43_foreign_sequences.py` (with `scripts/22` and `scripts/28`) | Organelle-like DNA in each assembly; 471 algal proteins aligned to it; every sequence with a non-animal match translated as a whole and its genes classified. |
| Places tested with reads | `40_locus_baits.py`, `41_locus_test.py`, `44_window_reads.sh` | Reads that run through, or stop at, each edge of an organelle-like stretch inside an assembly sequence. |
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

| Species | Run | Kind | Examined |
|---|---|---|---|
| *E. papillosa* | ERR14842698 | PacBio HiFi, the animal assembled | 8.3 Gb (one of three runs) |
| *E. atroviridis* | SRR34916404 | PacBio HiFi, the animal assembled | 30.6 Gb (all) |
| *E. timida* | SRR30016879 | PacBio HiFi, amplified library, the animals assembled | first 20 Gb |
| *E. chlorotica* | SRR7647856 | PacBio RS, an adult | 9.4 Gb (all) |
| *E. chlorotica* | SRR7647851 | Illumina pairs, a second adult | first 6 Gb |
| *E. chlorotica* | SRR768897 | Illumina pairs, eggs | 11.1 Gb (all) |
| *E. marginata* | DRR238951 | Illumina pairs, the animal assembled | first 4.5 Gb |
| *E. leucolegnote* | SRR41369319, SRR17036819 | Illumina pairs, two other animals | first 6 Gb and 5 Gb |

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
  floor; *E. papillosa* reads cover the genome about 8 times, so an insertion in a minority of
  cells could be missed; short reads test a join only through read pairs. The protein screen
  searches with 471 proteins, so it finds only the bacterial and foreign sequences that happen
  to carry a relative of one of them.

## Credits

Genome data for *E. crispata* and *E. papillosa*: Wellcome Sanger Institute Tree of Life and the
Aquatic Symbiosis Genomics project; please cite Middlebrooks et al. 2026 when using the
*E. crispata* data, and the papers of the other genomes (Cai et al. 2019; Maeda et al. 2021;
Guo & Kong 2026; and the submitters named in each GenBank record) when using theirs.
Photographs are from iNaturalist under Creative Commons licences; each is credited on the site
and in `docs/assets/photos/credits.json` (two are CC BY-NC, the rest CC BY). Analysis and website were produced with
Claude (Anthropic) for the repository owner. Fonts in `docs/assets/fonts/`: Bricolage Grotesque
and IBM Plex Mono, both under the SIL Open Font License 1.1.
