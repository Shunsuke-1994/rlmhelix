# Public data

The analyses use public data only. `scripts/download_data.py` fetches them (Python standard library and pandas;
reading Supplementary Table S8 as `.xlsx` also needs openpyxl, which is in the conda environment and in the `data`
extra, `pip install -e ".[test,data]"`).
Downloads are resumable: an interrupted file is kept as `<name>.part` and continued when the same command is run
again, complete files are skipped, and MD5 checksums are verified where the source publishes them (GENCODE, ENCODE).

```bash
python scripts/download_data.py data/raw --dry-run        # list URLs and destinations
python scripts/download_data.py data/raw --check          # also ask each server for the file sizes
python scripts/download_data.py data/raw                  # default groups: genome rmsk tree rfam ensembl (about 37 GB)
python scripts/download_data.py data/raw --only icshape screens polyasite bprna
python scripts/download_data.py data/raw --list-manual    # inputs that are not fetched by the script
```

Files are written to subdirectories of the destination (`gencode50/`, `ucsc_hg38/`, `rfam/`, `ensembl116/<species>/`,
`icshape/`, `screens/`, ...; `data/raw/` is ignored by git). The last column gives the group to pass to `--only`.

| Data | Source | Version | Used for | `--only` |
|---|---|---|---|---|
| Human genome and gene annotation | GENCODE | GRCh38 primary assembly; release 50 (`gencode.v50.annotation.gtf.gz`, `GRCh38.primary_assembly.genome.fa.gz`) | scan windows, annotation of candidates, splice sites | `genome` |
| RepeatMasker annotation | UCSC Genome Browser | hg38 `rmsk` table | repeat filter of the scan windows | `rmsk` |
| Vertebrate genomes for homology search | Ensembl | release 116; the 55 assemblies of `genome_panel_55.tsv` (unmasked `dna`; primary assembly for human, mouse and zebrafish, toplevel otherwise) | homology search in 20 species and expansion to 55 | `ensembl` |
| Tree of the vertebrate alignment | UCSC Genome Browser | `hg38.100way.nh` (and the common- and scientific-name versions) | phylogenetic baselines (RNAPhylo, EvoFold) | `tree` |
| Rfam families and genome hits | Rfam | 15.1: `Rfam.cm`, `Rfam.clanin`, `Rfam.seed`, database tables `full_region` and `family` | Rfam benchmark (seed alignments), Rfam controls of the homology search (covariance models), Rfam overlap of windows and candidates (genome hits) | `rfam` |
| Vertebrate whole-genome alignment | UCSC Genome Browser, as converted by GPN-MSA (Hugging Face `songlab/multiz100way-pigz`) | hg38 100-way alignment (multiz100way); human + 99 vertebrates, zarr (39 GB) | MSAs of EvoFam loci and RNAz 2005 regions | `multiz` |
| NCBI taxonomy | NCBI | `taxdump.tar.gz`, `nucl_gb.accession2taxid.gz` (not versioned; 2.5 GB) | only to rebuild the Rfam vertebrate benchmark from `Rfam.seed` | `taxonomy` |
| icSHAPE | GEO | GSE145805 (in vivo; HEK293, HeLa, K562, HepG2, H9) and GSE74353 (HEK293T, in vivo and in vitro); transcript sets for mapping: GENCODE release 26 (GTF, transcripts) and Ensembl release 74 (cDNA, ncRNA) | icSHAPE agreement and support | `icshape` |
| RNAz genome screen | Washietl et al. 2005 | supplementary BED files (Set1–3, P > 0.5 and P > 0.9; hg17, lifted to hg38) | comparison with earlier screens, RNAz 2005 benchmark | `screens` |
| EvoFold predictions | UCSC Genome Browser | hg19 `evofold`, hg18 `evofoldV2` | comparison with earlier screens | `screens` |
| CMfinder conserved RNA structures (CRSs) | Seemann et al. 2017 | release 2.1, Supplementary Data 1 with the 2020-07-28 hg38 coordinate patch | comparison with earlier screens | `screens` |
| PhastCons conserved elements | UCSC Genome Browser | hg17 `phastConsElements17way` | approximate input of RNAz 2005 and EvoFold | `screens` |
| liftOver chains | UCSC Genome Browser | hg17, hg18 and hg19 to hg38 | lifting earlier screens and EvoFam to hg38 | `screens` |
| EvoFam | Parker et al. 2011 | supplementary data (hg18, lifted to hg38) | EvoFam benchmark, comparison with earlier screens, pull-down baits | manual |
| eCLIP peaks, control RNA-seq (expression filter), RNA Bind-n-Seq, eCLIP signal tracks | ENCODE | file accessions in Supplementary Table S8 (250 eCLIP experiments, K562 and HepG2) | eCLIP enrichment, partner calls, genome-browser views | `encode` (needs the accession table) |
| polyA sites | PolyASite | 2.0, GRCh38.96 atlas clusters | genome-browser views | `polyasite` |
| bpRNA-new | MXfold2 data sets (Sato et al. 2021; Zenodo, doi:10.5281/zenodo.4430150) | 5,313 sequences of 30–300 nt | base-pair prediction (Supplementary Fig. S1) | `bprna` |

## Genome panel

`genome_panel_55.tsv` lists the 55 vertebrate genomes searched for homologues (also in Supplementary Table S8):
`scientific_name`, `common_name`, `taxon_id`, `lineage_group` (primate, non-primate mammal, non-mammalian tetrapod,
other vertebrate), `search_stage` (the 20 species of the seed search or the 35 further genomes of the expansion;
human is the query and was also searched for paralogous loci), `ensembl116_assembly`, `ensembl116_GCA` (INSDC
assembly accession), `assembly_kind` (`primary_assembly` or `toplevel`), `assembly_name_in_file` and `ensembl_slug`.
An identical copy is shipped with the package (`src/rlmhelix/data/genome_panel_55.tsv`, installed as
`rlmhelix/data/genome_panel_55.tsv`) and is the default `--panel` of `rlmhelix select`; a test checks that the two
copies are identical. The FASTA of each genome is
`https://ftp.ensembl.org/pub/release-116/fasta/<ensembl_slug>/dna/<Ensembl_slug>.<assembly_name_in_file>.dna.<assembly_kind>.fa.gz`
(all 55 checked on 2026-10-11). The paper's copies were downloaded from Ensembl's `current_fasta` directory on
2026-08-26; they have the same checksums as in release 115.

## Rfam genome hits

The paper read the significant Rfam 15.1 `full` hits on GRCh38 from the Rfam public database (table `full_region`
joined with `genseq` for the human proteome UP000005640) and kept those on the primary chromosomes (20,516 hits in
1,077 families). After the `rfam` group is downloaded, the script derives the same set from the database files and
writes `rfam/rfam15.1_human_GRCh38_hits.tsv`: rows of `full_region.txt.gz` with `type = full` and
`is_significant = 1` on the GenBank accessions of the GRCh38 chromosomes (CM000663.2–CM000684.2 for chr1–chr22,
CM000685.2 for chrX, CM000686.2 for chrY), with 0-based half-open coordinates, strand and the family name and type
from `family.txt.gz`.

## Inputs obtained by hand

- **EvoFam:** supplementary data of Parker et al., *Genome Res.* **21**, 1929–1943 (2011),
  https://doi.org/10.1101/gr.112516.110 (`EvoFam.csv`, `EvoFam.fa`; hg18).
- **ENCODE files:** pass Supplementary Table S8 (the `.xlsx` file, which needs openpyxl, or its Accessions sheet
  saved as TSV):
  `python scripts/download_data.py data/raw --only encode --encode-table Supplementary_Table_S8.xlsx --encode-dataset "eCLIP peaks"`
  (250 files; the other `dataset` values are `Control RNA-seq (expression filter)` (12), `RNA Bind-n-Seq` (96) and
  `eCLIP signal tracks` (84)). Each file is fetched from `https://www.encodeproject.org/files/<ENCFF accession>/@@download/`
  and checked against its ENCODE MD5. A plain text file of ENCFF accessions also works. The eCLIP peak files are
  written as `data/raw/encode/<ENCFF accession>.bed.gz`, which
  `rlmhelix catalogue build --eclip Supplementary_Table_S8.xlsx --eclip-dir data/raw/encode` reads.
- **Pull-down mass spectrometry:** the raw data and the protein identification and quantification will be deposited
  to ProteomeXchange via PRIDE and released upon publication (accession to be added). The bait sequences are in
  Supplementary Table S7.

The external programs (RNAz 2.1, SISSIz 0.1, Infernal 1.1.5, MAFFT 7.526, R-scape 2.6.16, tantan, Inverted Repeats
Finder, UCSC liftOver, bedtools 2.31.1, Riboseek) are software, not data; see `environments/`.

## Notes

- The GENCODE genome FASTA is gzip-compressed; decompress it (`gunzip -k`) before indexing with tools that need
  random access.
- The UCSC tables, GEO files and NCBI taxonomy are not versioned by their providers. On 2026-10-11, the GENCODE,
  Ensembl and icSHAPE files that were compared with the paper's copies had the same sizes (and, for GENCODE, the
  same MD5).
