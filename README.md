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

Two species have more than one assembly (checked at NCBI on 14 November 2026). *E. crispata* has
four: xgElyCris1.1, its alternate haplotype, and two assemblies from other laboratories, referred
to here as the Purdue and Tampa assemblies. *E. papillosa* has two: xgElyPapi4.2 and its alternate
haplotype. Every assembly was searched for the clock genes and for algal DNA
(`config/species_assemblies.tsv`). One assembly of each species serves as its reference, and the
status of a gene always refers to that assembly; the site also shows whether the gene is intact in
another assembly of the species.

**Clock genes.** Number of genes of each status in the reference assembly of each species:

| Species | Reference assembly | Intact | Complete, with a frameshift or stop codon | Partial or fragmented | Not found | Intact in at least one assembly |
|---|---|---|---|---|---|---|
| *Elysia crispata* | GCA_963854125.1 (chromosome) | 30 | 4 | 5 | 8 | 47 (four assemblies) |
| *Elysia chlorotica* | GCA_003991915.1 (scaffold) | 39 | 5 | 3 | 0 | |
| *Elysia marginata* | GCA_019649035.1 (scaffold) | 33 | 8 | 4 | 2 | |
| *Elysia timida* | GCA_043644045.1 (chromosome) | 44 | 1 | 1 | 1 | |
| *Elysia leucolegnote* | GCA_052815195.1 (chromosome) | 43 | 1 | 2 | 1 | |
| *Elysia atroviridis* | GCA_059052615.1 (chromosome) | 44 | 1 | 1 | 1 | |
| *Elysia papillosa* | GCA_977015215.2 (contig) | 43 | 0 | 3 | 1 | 44 (two assemblies) |

- Of the 47 genes, 27 are complete in the reference assembly of all seven species, and 35 when
  the other assemblies of *E. crispata* and *E. papillosa* are taken into account. The core loop
  is complete in every species. Where a gene is partial or missing, the cause is usually the
  assembly: the two assemblies built mainly from short reads contain the largest numbers of such
  genes.
- In the alternate haplotype of *E. papillosa*, 38 genes are intact. Hr3, which is partial in the
  reference assembly, is intact there, and its protein from the alternate haplotype is the one
  used in the alignment, the identity values and the species tree. Nemo and r-opsin-like 1 are
  incomplete in both assemblies, and Go-opsin 1 was found in neither.
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
| *E. crispata*, Purdue adults, Nanopore reads | 1,985 of 22,988 (8.6%) | 115 of 1,142 (10.1%) | 7.4% | 0 |
| *E. crispata*, Purdue adults, brains, Nanopore reads | 49 of 387 (12.7%) | 54 of 395 (13.7%) | 12.9% | 0 |
| *E. crispata*, Purdue egg mass, paired-end short reads | 0 of 1 (0%) | 78 of 136,554 (0.057%) | 0.86% | not applicable (short reads) |
| *E. crispata*, Tampa long-read library (amplified), PacBio HiFi reads (one read in 5 and one in 20 examined) | circular elements used instead: 37 of 16,335 (0.23%) | 114 of 29,015 (0.39%) | 6.3% | 0 |
| *E. crispata*, Tampa short-read library, paired-end short reads (estimated from one pair in 20) | 4,858 of 13,267,181 (0.037%) | 223 of 786,007 (0.028%) | 0.27% | not applicable (short reads) |
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
  (psbO, prk, lhc, chlD, chlG, chlH) has a match in the assembly. 99 scaffolds (1.74 Mb) are
  DNA of bacteria and other organisms, 1.10 Mb of it on scaffolds whose genes are closest to
  those of *Francisella* (see the census below).
- *E. marginata*: one 3.7 kb scaffold consists of *Bryopsis* plastid DNA. 233 scaffolds (3.60 Mb)
  are DNA of bacteria and other organisms (mycoplasma-like, flavobacteria-like and others), and
  bacterial stretches lie at the ends of two slug scaffolds (40 kb); the gene set published with
  the assembly includes genes on these scaffolds.
- *E. leucolegnote*: the unplaced sequence JBQKML010000381.1 begins with 863 bp of plastid DNA
  joined to a slug repeat. In two other individuals, the same plastid sequence continues as
  plastid DNA at that point (87 and 77 reads diverge from the assembly there, and none spans the
  junction). The reads of the assembled individual are not public, so the junction itself
  cannot be tested. 17 sequences (2.06 Mb) are DNA of bacteria and other organisms, 1.83 Mb of
  it on two sequences (JBQKML010000180.1 and JBQKML010000169.1).
- *E. papillosa*: contig CDRNNR020002991.1 (30.9 kb) of the reference assembly is a fragment of
  an algal mitochondrial genome; nine contigs (245 kb) are bacterial. The alternate haplotype
  (GCA_977015365.2) contains no organelle-like DNA and one bacterial contig (12 kb), and the
  algal mitochondrial contig has no counterpart in it.
- *E. timida* and *E. atroviridis*: no organelle DNA was found in the assembly; six sequences
  (160 kb) of *E. timida* and four (209 kb) of *E. atroviridis* are bacterial.

In the protein-level search of these species, the gene set published with an assembly was not
used as evidence that a sequence in that assembly is of animal origin, because the published
gene sets of *E. chlorotica* and *E. marginata* include genes on bacterial and plastid sequences.

**Bacteria and other organisms in every sequence.** The searches with proteins find a bacterial
sequence only when it carries one of the genes they look for, so each found a different part of
the foreign DNA. All sequences of the eleven assemblies were therefore examined in one way
(`scripts/genus/76_foreign_census.py`, with `orfscan.py` and the rules of `foreign_dna.py`; every
piece of foreign DNA in `results/genus/foreign_dna_census.tsv`, the details in
`foreign_dna_census.json.gz`). Bacterial, phage and giant-virus DNA is 85 to 95% coding, slug DNA
much less (medians of 12 to 28% in random 20-kb windows). Every sequence of which at least 60% lies
in open reading frames of 300 bases or more (Mycoplasma code where it alone reaches 85%; simple
sequence left out), or whose GC content differs from the assembly's by 0.08 or more, was translated
in full (DIAMOND blastx, as in `43_foreign_sequences.py`) and judged with `foreign_dna.judge`,
unless it is a tandem array (fewer than half of its 12-base words distinct; such arrays, with
periods too long for the simple-sequence mask, are common in the HiFi assemblies). Coding runs of
5-kb windows at the ends of slug sequences were translated too, and counted as foreign when they
carry at least three foreign coding stretches and twice as many as stretches matched best by sea
slugs and other animals. Foreign sequences that align to plastid
or algal mitochondrial genomes of the panel (or to the *Acetabularia* contigs of the Purdue
assembly) are organelle DNA and are listed with it above.

| Assembly | Bacteria and other organisms | Stretches at the ends of slug sequences | Organelle DNA among the candidates | Short sequences not assigned |
| --- | --- | --- | --- | --- |
| *E. crispata* xgElyCris1.1 | 2 (3.0 kb) | 0 | 1 (2.0 kb) | 33 (280 kb) |
| *E. crispata* alternate haplotype | 0 | 1 (5.0 kb) | 0 | 14 (199 kb) |
| *E. crispata* Purdue | 0 | 0 | 12 (703 kb) | 245 (1.37 Mb) |
| *E. crispata* Tampa | 175 (631 kb) | 0 | 1 (4.9 kb) | 859 (1.74 Mb) |
| *E. chlorotica* | 99 (1.74 Mb) | 0 | 5 (8.3 kb) | 991 (823 kb) |
| *E. marginata* | 233 (3.60 Mb) | 2 (40 kb) | 0 | 215 (715 kb) |
| *E. timida* | 6 (160 kb) | 0 | 0 | 7 (36 kb) |
| *E. leucolegnote* | 17 (2.06 Mb) | 0 | 0 | 23 (72 kb) |
| *E. atroviridis* | 4 (209 kb) | 0 | 0 | 22 (74 kb) |
| *E. papillosa* xgElyPapi4.2 | 9 (245 kb) | 0 | 1 (30.9 kb) | 4 (58 kb) |
| *E. papillosa* alternate haplotype | 1 (12 kb) | 0 | 0 | 30 (241 kb) |

Every sequence that the searches with plastid proteins, chloroplast regulators and gene sets
judged to be bacterial or otherwise foreign is among them. The only sequences listed by those
searches that are not are plastid DNA: an *Acetabularia* chloroplast contig of the Purdue assembly
with few reading frames (JAWDGP010006093.1), and two slug sequences with a plastid fragment
(JAWDGP010001684.1 of the Purdue assembly and JBQKML010000381.1 of *E. leucolegnote*). The three
sequences that the regulator and gene-set searches judged to be slug DNA despite dense reading
frames are slug DNA here too. Short sequences not assigned are those under 20 kb without any
coding stretch matched by another sea slug, another animal, an alga or a microbe. The organisms
can be named only roughly: on half of the sequences the bacterial genes are 37 to 52% identical
to their closest Swiss-Prot proteins (median 44%). In the Tampa assembly most foreign contigs are
shorter than 5 kb and AT-rich (median GC 26%); one, of 39 kb, is the T7-like phage described
below. In *E. marginata*, a few
GC-rich scaffolds carry RNA polymerase subunits, ribonucleotide reductase, a flavin-dependent
thymidylate synthase and thymidine kinase whose closest proteins are algal, fungal or bacterial,
a set of genes typical of large DNA viruses.

**Nuclear genes that regulate the chloroplast.** A chloroplast taken from an alga arrives without
the nuclear genes that control it in the alga. Whether any of them reached the genome of a slug,
perhaps in an ancestor long ago, was tested with 160 such genes (`config/plastid_regulators.tsv`):
the sigma factors SIG1 to SIG6 and other proteins of plastid transcription; the proteins of
retrograde signalling (GUN1, GUN4, magnesium chelatase, heme oxygenase, EXECUTER, SAL1, PTM, ABI4,
GLK1 and GLK2); light receptors and their partners (phytochromes, phototropins, cryptochromes, UVR8,
COP1, and the aureochrome of stramenopiles); regulators of photosynthesis (STN7, STN8, PGR5, the
thioredoxins and their reductases); the TOC and TIC import machinery; plastid division; and the
plant clock. Each *Arabidopsis* protein was joined by its counterparts in eleven algae, among them
the two siphonous green algae *Bryopsis* and *Ostreobium* and two stramenopiles related to
*Vaucheria* (829 proteins, `results/intermediate/regulator_queries.tsv`). All were aligned to the
eleven assemblies with miniprot, every place found was translated and compared with animal,
sea-slug, algal, plant and microbial proteins, and a place that matched an alga, plant or microbe
best was judged by the sequence that carries it (`scripts/genus/foreign_dna.py`). Matches to other
sea slugs were not counted as evidence of animal origin, so that a gene acquired by a common
ancestor of the slugs would still be found.

- None of the 2,092 places is an algal or plant gene on slug DNA
  (`results/genus/plastid_regulators.json`; by regulator and assembly in
  `plastid_regulators_matrix.tsv`; every place in `plastid_regulators_loci.json.gz`).
- 70 places are of algal or microbial type (their best algal, plant or microbial match scores at
  least 80 bits and more than any animal that is not a sea slug). 68 of them lie on foreign DNA of
  the assemblies: 47 sequences of bacteria and of a phage, and a bacterial stretch at the end of
  one *E. marginata* sequence. Long open reading frames cover 62 to 100% of these 47 sequences
  (median 90%), against medians of 12 to 28% in random windows of slug DNA, and of the coding
  stretches on them 2,214 match a microbe or an alga best and 47 another sea slug. The sigma
  factors found in *E. chlorotica*, *E. marginata* and *E. leucolegnote* are bacterial sigma-70
  genes on such sequences, one of them of a mycoplasma-like bacterium. The phage-type RNA
  polymerase found in the Tampa assembly lies on a 39 kb sequence that is 93% open reading frames
  and also encodes a DNA polymerase I and the ribonucleotide-reductase proteins NrdH and NrdI, as
  the genome of a T7-like phage does. The other two, CHLI of magnesium chelatase in *E. crispata*
  and its alternate haplotype, lie on kleptoplast DNA (chlI is a plastid gene in *Bryopsis*).
- 499 places are slug genes that match an animal protein best: relatives of the regulator where
  animals have the family (ferrochelatase, COP1, Cry-DASH, NUBPL for HCF101, retinol
  dehydrogenases for protochlorophyllide reductase), and otherwise genes that share a domain with
  it (protein kinases for STN7 and the phototropins, GTPases of the IMAP family for the TOC
  receptors). Another 60 are slug genes that match only other sea slugs and no algal or microbial
  protein. The remaining 1,460 alignments are not confirmed by the translated place: they are of
  low complexity (1,322) or short (91), and of the 47 longer ones, a more sensitive search
  (DIAMOND --ultra-sensitive, E <= 1) found significant matches for three, all to sea-slug
  proteins.
- 35 regulators were not found at all, among them GUN1, GUN5 (the H subunit of magnesium
  chelatase), GLK1 and GLK2, both phytochromes, ZTL and FKF1, and PTAC2, PTAC7 and PTAC10.

**Genes of non-animal origin in the annotated proteomes.** The 190,139 proteins of four
sacoglossans (*E. crispata*, Ensembl, longest isoform of each gene; *E. chlorotica*, *E. marginata*
and *Plakobranchus ocellatus*, UniProt) were compared with the same database, again without
counting matches to other sea slugs. A protein was kept when its best algal, plant or microbial
match scored at least 80 bits and better than its best match among 13 animals that are not sea
slugs, with no cap on the number of matches (with 25, large families of algal and sea-slug
proteins can crowd out the animal homologue). Of the 916 proteins kept (*E. crispata* 11,
*E. chlorotica* 170, *E. marginata* 699, *P. ocellatus* 36; `results/genus/alien_genes.json`), none
is a candidate for a gene acquired from an alga or a bacterium:

- six *E. crispata* "genes" are plastid genes (psaA, petB, petD, cemA and two parts of atpB) on
  the stretches of kleptoplast DNA within its chromosomes;
- 840 lie on foreign sequences (235 sequences: bacteria related to *Francisella* in
  *E. chlorotica*; mycoplasma-like, flavobacteria-like and other organisms in *E. marginata*), six
  on bacterial stretches at the ends of two *E. marginata* sequences, and one on a 2 kb sequence
  with no slug gene;
- the 26 on slug DNA are transposable elements (13) or members of conserved families whose algal
  or microbial match is less than 1.2 times the animal one (13: actin, ubiquitin, histone H3, Rab
  GTPases, a thioredoxin, a glutamate dehydrogenase);
- the *P. ocellatus* proteins could not be placed, because its assembly is not part of this
  study. Its 25 of microbial type are bacterial proteins (seven of them 80 to 97% identical to
  proteins of *Vibrio*); of the 11 of algal type, five are transposases, four are near-equal
  matches in conserved families, one is a bacterial-type transporter, and one is a SufB whose
  closest sea-slug match lies on a bacterial sequence of *E. marginata*.

**Cry-DASH2.** The second Cry-DASH gene, which only sacoglossans have, did not come from an alga.
In a tree of 61 Cry-DASH proteins and photolyases (365 positions; IQ-TREE, LG+F+G4, 1000
ultrafast bootstraps; `results/genus/crydash_tree/`), the Cry-DASH2 proteins of seven sea slugs
form one group (support 100) within the Cry-DASH proteins of gastropods (99) and of animals (100),
apart from those of algae, plants, fungi and bacteria. Neither *Bryopsis* nor *Ostreobium* has a
Cry-DASH protein in its proteome; their closest proteins are CPD photolyases.

## Results for *Elysia crispata*

**Clock genes.** The catalogue contains 47 genes: the core loop (Clock, Bmal1, Period, Timeless,
Timeout, Cry1, Cry2), the second loop (Pdp1, Vrille, E75, Hr3, Cwo), kinases, phosphatases,
degradation factors, output genes, photolyases and nine opsins. All 47 are present in the
species, but the Sanger Institute assembly xgElyCris1.1 is incomplete: 30 genes
are intact in it, 4 are complete but contain a frameshift or stop codon in the assembly
sequence, 5 are partial or fragmented (only the 3′ exons of Bmal1 are present, on an unplaced
scaffold) and 8 are missing or represented only by a divergent copy. Each of these 17 genes is
intact in the Purdue or the Tampa assembly, and HiFi reads from the sequenced individual cover
their exons, so the missing genes reflect gaps in the assembly and not gene loss. All four
assemblies were searched for every gene: 30 genes are intact in xgElyCris1.1, 12 in its alternate
haplotype, 46 in the Purdue assembly and 42 in the Tampa assembly, and every gene is intact in at
least one of them. Cry2 is intact only in xgElyCris1.1; the first 25 residues of the protein
model do not align to any of the other three assemblies. The following genes were not found in
any assembly: Dec1/2, Fbxl3/21, Tyf, Cipc, Ciart, Aanat and Asmt.

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
- The Purdue assembly (other animals, sequenced with Nanopore reads) contains 11 contigs
  (646 kb) of *Acetabularia* chloroplast DNA, which derive from the kleptoplasts of those
  animals, an 11 kb plastid-like contig of unknown origin, and a few positions at which
  chloroplast DNA is joined to nuclear DNA within a contig. These junctions were tested with the
  reads of those animals; see the gene-transfer page of the website and the next section.

**Other sequenced individuals.** The reads behind the Purdue and Tampa assemblies were analysed
in the same way as those of the other species (`config/species_reads.tsv`; page "Other
individuals" of the website): Nanopore reads of adults (33.0 Gb) and of dissected brains (5.1 Gb)
and read pairs of an egg mass (63.4 Gb) of the Purdue isolate, and a PacBio HiFi library
(17.1 Gb) and an Illumina library (60.0 Gb) from the University of Tampa.

- Kleptoplasts. In the Purdue adults, the plastid DNA (0.11% of the bases) matches
  *P. capitatus* (44%), an *Acetabularia* species (31%), *P. lamourouxii* (9%), a relative of
  *Halimeda* (11%) and *Udotea flabellum* (1%). The proportions differ between the five flow
  cells, so the DNA came from more than one animal. The Tampa short-read library contains
  almost exclusively *P. lamourouxii* plastid DNA (98% of it; 6.1% of the bases). The Tampa
  long-read library contains plastid DNA of an alga without a sequenced plastid genome, whose
  rbcL is 93.7% identical to that of *Rhipiliopsis profunda* and identical, over 561 bp, to a
  sequence from the kleptoplasts of *E. clarki* (KJ138404.1).
- The two Tampa libraries are registered under one sample record but came from two animals:
  their mitochondrial genomes differ by 42 substitutions, each fixed in its library, and
  neither library contains the plastid DNA of the other.
- Circular elements. In the Tampa long-read library, which was made by whole-genome
  amplification, 26% of the bases are mitochondrial DNA and 3.6% come from small circular DNA
  elements of 4.9 to 5.9 kb (nine families), which the amplification copied as tandem arrays.
  They carry a gene related to open reading frames in the plastid genomes of siphonous algae
  (in *R. profunda*, a gene with a phage- or plasmid-type DNA primase domain), and were first
  mistaken for plastid contigs. They are counted as a type of DNA of their own. The plastid-like
  contig JAWDGP010006602.1 of the Purdue assembly is related to them.
- Germline. The egg mass contains no algal DNA: nine reads in 63.4 Gb match a plastid genome
  (of a land plant, a red alga and green algae that are not siphonous), and none is joined to
  nuclear DNA.
- The Purdue assembly, examined with the egg reads: 292 stretches (3.4 Mb) have no egg reads.
  Fifteen of them (754 kb) are *Acetabularia* chloroplast DNA, four of which (116 kb) had not
  been recognised before; one is the plastid-like contig; 189 (1.8 Mb) are slug DNA that the
  parents of the egg mass lack; 87 could not be assigned, 50 of them tandem repeats. The
  assembly contains the Nanopore sequencing adapter at 58 places, 18 of them inside contigs.
  At nine positions it joins chloroplast DNA to other DNA within a contig; no egg read spans
  any of them, and two contain the adapter (`config/assembly_junctions.tsv`).
- Junctions in the reads. Chloroplast DNA is joined to nuclear DNA about as frequently as
  mitochondrial DNA is (table above). In the Nanopore reads of the adults, 27% of the junctions
  with chloroplast DNA contain the sequencing adapter, that is, two molecules were read as
  one. No read of the Purdue animals contains chloroplast DNA flanked on both sides by nuclear
  DNA. One junction was read twice, from the two strands of one molecule, which is a chain of
  three unrelated fragments. In the Tampa long-read library, one read contains about 300 bases
  of a circular element between two adjacent nuclear segments; the read folds back, and its
  second copy of the locus lacks the insert, so the insert arose during amplification.

The website presents the evidence in full, together with the limitations of each analysis.

## Repository layout

| Path | Contents |
|---|---|
| `docs/` | The website (static). `index.html` lists the species, `compare.html` compares their clock proteins, and each species has its own folder (`docs/elysia-crispata/`, …). The data are embedded in `docs/data/<species>.js` and `docs/data/genus.js`; the files offered for download are in `docs/downloads/` and `docs/elysia-crispata/downloads/`. The photographs and their credits are in `docs/assets/photos/`. |
| `results/` | Tables for *E. crispata*: `clock_genes.tsv` / `.json`, `sequences/` (proteins and coding sequences), `kleptoplast.json`, `kleptoplast_barcodes.tsv`, `hgt.json`, `plastid_fragments_in_assembly.tsv`, `penicillus_plastome_genes.bed`. |
| `results/genus/` | For each of the other species: the clock-gene catalogue (`<species>.clock.json`; `elysia-papillosa.alt.clock.json` for the alternate haplotype), the composition of its reads (`<species>.kleptoplast.json`) and the gene-transfer analyses (`<species>.transfer.json`); for the other sequenced individuals of *E. crispata*, both in `elysia-crispata.other_animals.json`, with the stretches of the Purdue assembly that are absent from egg DNA, the positions of adapter sequence in that assembly and the circular elements of the Tampa long-read library. Also the gene-by-species table (`clock_matrix.tsv`), protein alignments (`alignments/`), the concatenated alignment and species tree, and the intermediate tables of the read and assembly analyses: the reference panel (`panel.tsv.gz`), barcode records and sequences, plastid contigs, read junctions, every protein match (`protein_screen_loci.json.gz`) and every coding region on the bacterial and other foreign sequences (`foreign_sequences.json.gz`). The search for chloroplast regulators (`plastid_regulators.json`, `plastid_regulators_matrix.tsv`, `plastid_regulators_loci.json.gz`), the genes of non-animal origin in four proteomes (`alien_genes.json`), the Cry-DASH tree (`crydash_tree/`) and the inventory of foreign DNA in every sequence of the eleven assemblies (`foreign_dna_census.tsv`, details in `foreign_dna_census.json.gz`). |
| `results/trees/` | Trimmed alignments of the gene families and FastTree trees in Newick format. |
| `results/intermediate/` | Small intermediate tables: orthologue calls, choice of gene model, read support, expression, the list of sequences in the reference panel, chimeric-read junctions, boundary evidence, the inputs to the protein-level search and the proteins of the search for chloroplast regulators (`regulator_queries.fa` / `.tsv`). |
| `config/` | Curated inputs: reference clock proteins (`clock_seeds.tsv`), the gene catalogue with its anchors (`clock_catalogue.tsv`), genes searched for but not found (`clock_not_found.tsv`), the species and their reference assemblies (`species.tsv`), the further assemblies of *E. crispata* and *E. papillosa* (`species_assemblies.tsv`), the read sets analysed (`species_reads.tsv`), the contigs of the Purdue assembly that were added to the reference panel (`panel_extra.tsv`), the positions at which that assembly joins chloroplast DNA to other DNA (`assembly_junctions.tsv`), the text of the conclusions shown on the species pages (`species_text.json`), the survey of public *Elysia* data (`species_catalogue.tsv`), and the plant and algal genes that regulate the chloroplast (`plastid_regulators.tsv`). |
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
| Clock genes | `09_clock_refs.sh`, `10_project_clock_genes.sh`, `10b_miniprot_in_parts.sh`, `11_clock_catalogue.py` | Alignment of the 47 *E. crispata* proteins to each genome with miniprot (in parts, for a genome too large for the available memory); addition of annotated proteins where available; a reciprocal best-hit test against all *E. crispata* proteins; selection of one model per gene. |
| Comparison | `12_compare_species.py`, `13_domains.sh`, `14_other_assemblies.py` | Protein alignments (MAFFT), completeness relative to the alignment positions present in most sequences, pairwise identity, the species tree (IQ-TREE), Pfam domains. For a second assembly of a species, each gene model is added to the finished alignment and scored against the same positions; where a gene is complete only in the second assembly, that protein represents the species. |
| Reference panel | `30_build_panel.py`, `30b_mask_and_index.sh`, `30c_extend_panel.py`, `35b_hide_panel_extra.py` | Plastid, mitochondrial and ribosomal DNA reference sequences for all species; masking of low-complexity sequence for the first pass; for the other individuals of *E. crispata*, addition of the chloroplast contigs of the Purdue assembly to the panel, and masking of these contigs in that assembly. |
| Reads | `fetch_ranges.py`, `31c_screen_run.py`, `31_select_hits.py`, `31b_slim_rdna.py`, `31d_sample_blocks.py`, `33_pass2.sh`, `sam_segments.py`, `49_analyse_library.sh`, `34_read_analysis.py`, `37_pair_analysis.py`, `55_end_piece_control.py` | Processing of each run in full, in batches, and mapping to the panel; base-level alignment of the reads with a match, and of a random sample of the library, to the assembly of the species and to the panel; composition of the library; every read or read pair in which two types of DNA are joined, with the same counts for mitochondrial DNA and for the random sample; the type of DNA at the end of the fragment in read pairs that span a junction. |
| Single molecules | `46_candidates.py`, `51_recurrent_junctions.py`, `47_library_control.sh`, `47_shape_reads.py`, `52_fetch_reads.py`, `53_join_in_short_reads.py`, `54_internal_pieces.py` | Classification of long reads by structure (organelle DNA flanked by nuclear DNA, organelle DNA at one end, or the reciprocal structure), with short segments confirmed at base level; comparison of junctions between reads to find any that recur in independent molecules; PCR duplicates; the background rate of chimeras, measured on every read of a small run; retrieval of individual reads by name; search for a known junction in read pairs. |
| Identification of the alga | `36_barcodes.py`, `36b_barcode_reads_set_aside.py`, `38_plastid_contigs.sh`, `38b_plastid_contigs_long.sh`, `38c_tandem_contigs.py`, `58_mt_haplotypes.py` | rbcL and tufA extracted from long reads, or from plastid contigs assembled from short reads (MEGAHIT), and compared with GenBank; contigs assembled from long reads without a close reference (hifiasm), and the test that showed them to be tandem arrays of small circular elements; comparison of two libraries by their mitochondrial genomes. |
| Assemblies | `35_assembly_screen.py`, `42_protein_screen_inputs.py`, `43_foreign_sequences.py` (with `scripts/22` and `scripts/28`) | Organelle-like DNA in each assembly; alignment of 471 algal proteins to each assembly; translation in full of every sequence with a non-animal match, and classification of its genes. |
| Read-based tests of loci | `40_locus_baits.py`, `41_locus_test.py` | Reads that span, or end at, each boundary of an organelle-like region within an assembly sequence, and at the locus of each read in which organelle DNA is flanked by nuclear DNA. |
| An assembly examined with egg reads | `59_germline_coverage.py`, `60_classify_absent_stretches.py`, `56_adapter_in_assembly.py` | Stretches of the Purdue assembly of *E. crispata* on which no read of the egg mass lies, and their classification (chloroplast, bacterial, slug DNA, tandem repeats); Nanopore adapter sequence inside contigs. |
| Chloroplast regulators and genes of non-animal origin | `70_regulator_queries.py`, `71_regulator_screen.sh`, `73_regulator_summary.py`, `72_alien_proteomes.py`, `75_alien_candidates.py`, `foreign_dna.py`, `74_crydash_tree.sh` | The *Arabidopsis* regulators and their counterparts in eleven algae (reciprocal best hits by family); their alignment to all eleven assemblies and the classification of every place; the kind of DNA that carries a place of algal or microbial type (open reading frames, genes shared with other sea slugs, GC content); the comparison of four annotated proteomes with the same database, the placing of each candidate on its assembly and the same judgement of its sequence; the tree of Cry-DASH proteins and photolyases. |
| Foreign DNA in every sequence | `76_foreign_census.py`, `orfscan.py` (with `foreign_dna.py`) | Open reading frames, GC content and 12-base word diversity of every sequence of the eleven assemblies; translation of every sequence dense in reading frames or of unusual GC content, and of coding runs at the ends of slug sequences, and the judgement of each by the rules of `foreign_dna.py`; organelle DNA among them. |
| Collect and build | `39_collect_species.py`, `45_export_tables.py`, `chrome.py`, `21_build_species_pages.py`, `20_build_genus_site.py` | Result files and small tables; the shared page frame; the kleptoplast and gene-transfer pages of each species; the home page, comparison page, clock-gene pages and methods. |
| Order of the analyses | `run_screens.sh` | The commands for the six additional species, for the second assembly of *E. papillosa*, for the other individuals of *E. crispata* (step `crispata_animals`) and for the search for chloroplast regulators (step `regulators`), in the order in which they were run (`DRY=1` prints them). |

The scripts were written for a single machine. The locations of the tools
(`/home/claude/tools/mamba/envs/bio/bin`, which can be overridden with `BIO_BIN` in the scripts
that read it) and the working directories are passed on the command line or set at the top of
each script, and the header of each script describes its usage. Software used: minimap2,
miniprot, DIAMOND, BLAST+, MAFFT, trimAl, FastTree, IQ-TREE, HMMER (Pfam-A), salmon, samtools,
seqkit, sra-tools, MEGAHIT, hifiasm, bedtools, Python 3 with Biopython.

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
| Reads of the Purdue isolate (adults, brains and an egg mass) | PRJNA987316: 12 MinION runs and one Illumina run (table below) |
| Reads of the University of Tampa (two animals) | PRJNA1337476: SRR35727887, SRR35727888 |

Assemblies of the other species: *E. chlorotica* GCA_003991915.1, *E. marginata* GCA_019649035.1,
*E. timida* GCA_043644045.1, *E. leucolegnote* GCA_052815195.1, *E. atroviridis* GCA_059052615.1,
*E. papillosa* GCA_977015215.2 and its alternate haplotype GCA_977015365.2 (details in
`config/species.tsv` and `config/species_assemblies.tsv`).

Reads analysed for the other individuals of *E. crispata* and for the other species (details in
`config/species_reads.tsv`):

| Species | Runs | Platform and source of DNA | Bases analysed |
|---|---|---|---|
| *E. crispata* | SRR25020508 to SRR25020511, SRR25020513 to SRR25020518, SRR25020522 | Oxford Nanopore MinION (ligation libraries, five flow cells); whole bodies of adults of the Purdue isolate | 33.0 Gb |
| *E. crispata* | SRR25020519 | Oxford Nanopore MinION; brains of adults of the Purdue isolate | 5.1 Gb |
| *E. crispata* | SRR25020525 | Illumina NovaSeq 6000 (paired-end); an egg mass of the Purdue isolate | 63.4 Gb |
| *E. crispata* | SRR35727887 | PacBio Sequel II (HiFi, library made by whole-genome amplification); an animal from Key Largo (University of Tampa) | 17.1 Gb |
| *E. crispata* | SRR35727888 | Illumina HiSeq 4000 (paired-end); a second animal, registered under the same sample record | 60.0 Gb |
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
- The Tampa assembly was examined only as far as was needed to interpret xgElyCris1. Its two
  read libraries came from two animals, so it may combine the genomes of two individuals. The
  Purdue assembly was examined with egg reads in windows of 8 kb, so foreign DNA in shorter
  stretches would not have been found in this way.
- Other individuals of *E. crispata*: the Nanopore libraries join unrelated molecules in about
  one read in ten that contains organelle DNA. The DNA of the Purdue adults came from several
  animals, and the egg mass from other animals again. The Tampa long-read library was amplified,
  its plastid DNA has no close reference sequence (so it could be identified but not searched
  for junctions with nuclear DNA), and of its mitochondrial and plastid-like reads a sample was
  examined (one in 20 and one in five); in the Tampa short-read library, junctions were analysed
  in one in 20 of the read pairs that match an organelle genome. The relative of *Halimeda* and
  the *Acetabularia* species in the Purdue adults are under-represented, because the first has
  no close reference and the reference for the second is incomplete.
- The identification of the second kleptoplast species rests on a single 562-bp reference
  sequence.
- Other species: the *E. leucolegnote* reads are from two individuals other than the one
  assembled, and algae without a close reference sequence are under-represented in them; the
  plastid share given for *E. timida* is a minimum. Every run was analysed in full, but in the
  second *E. chlorotica* adult, in which a ninth of the DNA is chloroplast DNA, junctions were
  analysed in one in 50 of the read pairs that match an organelle genome. The protein-level
  search uses 471 query proteins, so it detects only those bacterial and other foreign
  sequences that encode a homologue of one of them.
- Chloroplast regulators: a gene that has diverged beyond recognition by protein alignment, or
  one too short to align, would escape the search; the 35 regulators found nowhere include such
  short or fast-evolving proteins. The proteome comparison depends on the published gene sets,
  and the *P. ocellatus* proteins were not placed on its assembly.
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
