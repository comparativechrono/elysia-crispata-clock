# Elysia crispata clock

Circadian clock genes, kleptoplast identity and a test for algal gene transfer in the lettuce
sea slug *Elysia crispata*, worked out from the public Sanger genome data for specimen
xgElyCris1 (BioProject [PRJEB67629](https://www.ebi.ac.uk/ena/browser/view/PRJEB67629);
Middlebrooks et al. 2026, [Wellcome Open Research 11:148](https://wellcomeopenresearch.org/articles/11-148)).

The results are presented as a small website in `docs/`, served by GitHub Pages:
**https://comparativechrono.github.io/elysia-crispata-clock/**

Nothing here measures a rhythm. Gene roles are inferred from homology with flies and mammals,
and the work has not been reviewed by the genome authors.

## What was found

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
- No annotated protein of either annotation is algal apart from those six, and aligning 471
  plastid and photosynthesis proteins directly to all four assemblies finds no algal gene among
  slug genes. Eleven contigs of the Purdue assembly (646 kb) carry chloroplast-type genes but
  appear to be another organism's DNA; they are absent from the xgElyCris1 reads.

The website gives the evidence, the caveats and what the screens cannot rule out.

## Layout

| Path | Contents |
|---|---|
| `docs/` | The website (static; `index.html`, `kleptoplasts.html`, `gene-transfer.html`, `methods.html`). Data are embedded in `docs/data/clock.js`; files offered for download are in `docs/downloads/`. |
| `results/` | Final tables: `clock_genes.tsv` / `.json`, `sequences/` (proteins and coding sequences), `kleptoplast.json`, `kleptoplast_barcodes.tsv`, `hgt.json`, `plastid_fragments_in_assembly.tsv`, `penicillus_plastome_genes.bed`. |
| `results/trees/` | Family alignments (trimmed) and FastTree trees in Newick format. |
| `results/intermediate/` | Small intermediate tables: orthologue calls, model choice, read support, expression, the reference panel listing, chimeric-read junctions, boundary evidence, protein-screen inputs. |
| `config/` | The curated inputs: reference clock proteins (`clock_seeds.tsv`), the gene catalogue with its anchors (`clock_catalogue.tsv`) and genes searched for but not found (`clock_not_found.tsv`). |
| `scripts/` | The pipeline, numbered in the order it runs. |

## Scripts

| Step | Scripts | What they do |
|---|---|---|
| Fetch | `00_fetch_reads.sh`, `00_fetch_proteomes.sh`, `00_fetch_algal_proteomes.sh`, `ncbi_fetch.py` | HiFi runs from the SRA open-data bucket as compressed FASTA; comparison proteomes from UniProt. |
| Clock genes | `01`–`09`, `miniprot_parse.py`, `tree_report.py` | Reference proteins, DIAMOND searches, family alignments and trees, reciprocal best hits, choice of a protein model per gene, presence in each of four assemblies, Pfam domains, tree figures. |
| Reads, pass 1 | `10`–`12`, `16` | Every HiFi read mapped to a panel of plastid, mitochondrial, rDNA and bacterial sequences plus clock-gene exons; reads classified and extracted; exon read support. |
| Kleptoplasts | `13` | rbcL and tufA haplotypes read out of the reads and compared with GenBank barcodes. |
| Chimeric reads | `14`, `15` | Base-level re-alignment of reads with an unexplained part; junctions between slug, plastid and bacterial DNA, with PCR duplicates collapsed. |
| RNA | `20` | salmon quantification of the head RNA-seq library. |
| Assemblies | `21`–`25`, `27`–`29` | Annotated proteins tested for algal affinity; plastid DNA in the four assemblies; read evidence at each plastid/slug boundary; plastid and photosynthesis proteins aligned to the assemblies; checks on foreign contigs. |
| Collect and build | `07`, `26`, `30`, `31` | Results tables and the website. |

The scripts were written for one working machine: tool locations
(`/home/claude/tools/mamba/envs/bio/bin`, overridable with `BIO_BIN` where the script reads it)
and working directories are passed on the command line or set at the top of each script. Each
script states its usage in its header. Software used: minimap2, miniprot, DIAMOND, BLAST+, MAFFT,
trimAl, FastTree, HMMER (Pfam-A), salmon, samtools, seqkit, sra-tools, Python 3 with Biopython.

To rebuild the website from the results already in the repository:

```
python3 scripts/30_build_site.py . results/trees/figures   # second argument optional
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

The reference panel used to screen the reads (241 green-algal plastid genomes and others) is
listed in `results/intermediate/panel.tsv`.

## Caveats

- xgElyCris1.1 represents about 60% of the read k-mers; absence from it means little.
- Six protein models are translations guided by another species' protein and may miss
  lineage-specific exons.
- Tree supports are FastTree local supports. Eleven gene names rest mainly on reciprocal best hits.
- Expression values come from one head sample at one unknown time of day.
- The boundary test and the check of the Purdue contigs used two of the three HiFi runs (53% of
  the data).
- The second kleptoplast species is named from a single 562-bp reference sequence.

## Credits

Genome data: Wellcome Sanger Institute Tree of Life and the Aquatic Symbiosis Genomics project;
please cite Middlebrooks et al. 2026 when using them. Analysis and website were produced with
Claude (Anthropic) for the repository owner. Fonts in `docs/assets/fonts/`: Bricolage Grotesque
and IBM Plex Mono, both under the SIL Open Font License 1.1.
