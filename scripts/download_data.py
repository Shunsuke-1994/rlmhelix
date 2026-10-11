"""Download the public inputs of the RLMHelix pipelines.

    python scripts/download_data.py data/raw                     # default groups (genome rmsk tree rfam ensembl)
    python scripts/download_data.py data/raw --only genome rfam  # selected groups
    python scripts/download_data.py data/raw --dry-run           # list URLs and destinations, download nothing
    python scripts/download_data.py data/raw --check             # also ask each server for the file size
    python scripts/download_data.py data/raw --list-manual       # inputs that are not fetched by this script
    python scripts/download_data.py data/raw --only encode --encode-table Supplementary_Table_S8.xlsx \
        --encode-dataset "eCLIP peaks"                          # ENCODE files listed in Supplementary Table S8

Files go to subdirectories of DEST (gencode50/, ucsc_hg38/, rfam/, ensembl116/<species>/, icshape/, screens/,
polyasite/, bprna/, multiz100way/, taxonomy/, encode/). A download is written to <name>.part and resumed with an HTTP
range request when the script is run again; files whose size equals the size reported by the server are skipped. MD5
checksums are verified where the source publishes them (GENCODE, ENCODE). After the `rfam` group, the significant
Rfam 15.1 hits on the GRCh38 primary chromosomes are written to rfam/rfam15.1_human_GRCh38_hits.tsv (BED-style
coordinates).

Standard library and pandas (pandas reads the genome panel and the ENCODE accession table). Reading Supplementary
Table S8 as .xlsx (--encode-table) also needs openpyxl: it is in environments/environment.yml, and with pip it is the
`data` extra (pip install -e ".[test,data]").
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

PANEL = Path(__file__).resolve().parents[1] / "data" / "genome_panel_55.tsv"
USER_AGENT = "rlmhelix-download/0.1"
CHUNK = 1 << 20

GENCODE = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"
UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath"
RFAM = "https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1"
ENSEMBL = "https://ftp.ensembl.org/pub"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo"
ENCODE = "https://www.encodeproject.org"

DEFAULT_GROUPS = ["genome", "rmsk", "tree", "rfam", "ensembl"]
OPTIONAL_GROUPS = ["icshape", "screens", "polyasite", "bprna", "multiz", "taxonomy", "encode"]

# GenBank accessions of the GRCh38 primary chromosomes, as used by Rfam (UP000005640)
GRCH38_CHROMOSOMES = {f"CM{662 + k:06d}.2": f"chr{k}" for k in range(1, 23)} | {
    "CM000685.2": "chrX", "CM000686.2": "chrY"}


@dataclass
class Item:
    group: str
    url: str
    path: str               # relative to DEST
    md5: str | None = None  # expected MD5, or "MD5SUMS" to look it up in the directory's MD5SUMS file


def ensembl_items(panel: Path) -> list[Item]:
    """The 55 genomes of the homology search (Ensembl release 116, unmasked dna)."""
    genomes = pd.read_csv(panel, sep="\t")
    items = []
    for g in genomes.itertuples():
        name = f"{g.ensembl_slug.capitalize()}.{g.assembly_name_in_file}.dna.{g.assembly_kind}.fa.gz"
        items.append(Item("ensembl", f"{ENSEMBL}/release-116/fasta/{g.ensembl_slug}/dna/{name}",
                          f"ensembl116/{g.ensembl_slug}/{name}"))
    return items


def static_items() -> list[Item]:
    items = [
        # GENCODE release 50: annotation and genome of the scan
        Item("genome", f"{GENCODE}/release_50/gencode.v50.annotation.gtf.gz", "gencode50/gencode.v50.annotation.gtf.gz",
             "MD5SUMS"),
        Item("genome", f"{GENCODE}/release_50/GRCh38.primary_assembly.genome.fa.gz",
             "gencode50/GRCh38.primary_assembly.genome.fa.gz", "MD5SUMS"),
        # RepeatMasker annotation used by the window filter
        Item("rmsk", f"{UCSC}/hg38/database/rmsk.txt.gz", "ucsc_hg38/rmsk.txt.gz"),
        # tree of the UCSC hg38 100-way alignment (phylogenetic baselines)
        *[Item("tree", f"{UCSC}/hg38/multiz100way/{name}", f"ucsc_hg38/{name}")
          for name in ("hg38.100way.nh", "hg38.100way.scientificNames.nh", "hg38.100way.commonNames.nh")],
        # Rfam 15.1: covariance models (cmsearch controls), seed alignments (benchmark), genome hits (enrichment)
        *[Item("rfam", f"{RFAM}/{name}", f"rfam/{name}") for name in ("Rfam.cm.gz", "Rfam.clanin", "Rfam.seed.gz")],
        *[Item("rfam", f"{RFAM}/database_files/{name}", f"rfam/{name}")
          for name in ("full_region.txt.gz", "family.txt.gz")],
        # icSHAPE (GSE145805 in vivo, five cell lines; GSE74353 HEK293T) and the transcript sets they are mapped to
        *[Item("icshape", f"{GEO}/samples/GSM4333nnn/{name.split('_')[0]}/suppl/{name}", f"icshape/GSE145805/{name}")
          for name in ("GSM4333256_293.out.txt.gz", "GSM4333258_HeLa.out.txt.gz", "GSM4333260_K562.out.txt.gz",
                       "GSM4333262_HepG2.out.txt.gz", "GSM4333264_H9.out.txt.gz")],
        *[Item("icshape", f"{GEO}/series/GSE74nnn/GSE74353/suppl/{name}", f"icshape/GSE74353/{name}")
          for name in ("GSE74353_HS_293T_icSHAPE_InVivo_BaseReactivities.txt.gz",
                       "GSE74353_HS_293T_icSHAPE_InVitro_BaseReactivities.txt.gz")],
        *[Item("icshape", f"{GENCODE}/release_26/{name}", f"icshape/gencode26/{name}", "MD5SUMS")
          for name in ("gencode.v26.annotation.gtf.gz", "gencode.v26.transcripts.fa.gz")],
        Item("icshape", f"{ENSEMBL}/release-74/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.74.cdna.all.fa.gz",
             "icshape/ensembl74/Homo_sapiens.GRCh37.74.cdna.all.fa.gz"),
        Item("icshape", f"{ENSEMBL}/release-74/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh37.74.ncrna.fa.gz",
             "icshape/ensembl74/Homo_sapiens.GRCh37.74.ncrna.fa.gz"),
        # earlier genome screens and liftOver chains to hg38
        *[Item("screens", f"https://www.tbi.univie.ac.at/papers/SUPPLEMENTS/ncRNA/bed/set{s}_{p}.bed",
               f"screens/rnaz2005_hg17/set{s}_{p}.bed") for s in (1, 2, 3) for p in (50, 90)],
        Item("screens", f"{UCSC}/hg19/database/evofold.txt.gz", "screens/evofold_hg19/evofold.txt.gz"),
        Item("screens", f"{UCSC}/hg18/database/evofoldV2.txt.gz", "screens/evofoldV2_hg18/evofoldV2.txt.gz"),
        Item("screens", "https://rth.dk/resources/crs/vert/v2.1/download/"
             "Supplemental_Data_1__patch_2020-07-28_fix_hg38_coordinates.csv.gz",
             "screens/cmfinder_crs_v2.1/Supplemental_Data_1__patch_2020-07-28_fix_hg38_coordinates.csv.gz"),
        Item("screens", f"{UCSC}/hg17/database/phastConsElements17way.txt.gz",
             "screens/phastcons_hg17/phastConsElements17way.txt.gz"),
        *[Item("screens", f"{UCSC}/{a}/liftOver/{a}ToHg38.over.chain.gz", f"screens/liftover/{a}ToHg38.over.chain.gz")
          for a in ("hg17", "hg18", "hg19")],
        # polyA sites of the genome-browser views
        Item("polyasite", "https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/"
             "atlas.clusters.2.0.GRCh38.96.bed.gz", "polyasite/atlas.clusters.2.0.GRCh38.96.bed.gz"),
        # bpRNA-new (data/bpRNAnew_dataset in the MXfold2 data sets, Zenodo 10.5281/zenodo.4430150)
        Item("bprna", "https://zenodo.org/records/4430150/files/mxfold2-data.tar.gz", "bprna/mxfold2-data.tar.gz"),
        # multiz100way in the zarr form of GPN-MSA (human + 99 vertebrates, 42 GB), for the EvoFam and RNAz 2005 MSAs
        Item("multiz", "https://huggingface.co/datasets/songlab/multiz100way-pigz/resolve/main/99.zarr.tar.gz",
             "multiz100way/99.zarr.tar.gz"),
        # NCBI taxonomy, only to rebuild the Rfam vertebrate benchmark from Rfam.seed (nucl_gb is 2.7 GB)
        Item("taxonomy", "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", "taxonomy/taxdump.tar.gz"),
        Item("taxonomy", "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
             "taxonomy/nucl_gb.accession2taxid.gz"),
    ]
    return items


MANUAL = """\
Inputs not fetched by this script:

  EvoFam (Parker et al. 2011, Genome Res. 21:1929-1943, doi:10.1101/gr.112516.110)
      Supplementary data (EvoFam.csv, EvoFam.fa; hg18). Lifted to hg38 with the hg18ToHg38 chain (`screens` group).
  ENCODE eCLIP peaks (250), control RNA-seq of the expression filter (12), RNA Bind-n-Seq tables (96),
  eCLIP signal tracks of the genome-browser views (84)
      File accessions are in the Accessions sheet of Supplementary Table S8. Pass the .xlsx (needs openpyxl) or
      that sheet saved as TSV: `--only encode --encode-table <file> [--encode-dataset "eCLIP peaks"]`. Each file
      is fetched from https://www.encodeproject.org/files/<ENCFF accession>/@@download/ and checked against its
      ENCODE MD5.
  Pull-down mass spectrometry (raw data, identification and quantification)
      Will be deposited to ProteomeXchange via PRIDE and released upon publication (accession to be added).
      Bait sequences are in Supplementary Table S7; the CMfinder alignments of the CM-line baits come from
      Hou et al. 2021.
  Software (not data): RNAz 2.1, SISSIz 0.1, Infernal 1.1.5, MAFFT 7.526, R-scape 2.6.16, tantan,
      Inverted Repeats Finder, UCSC liftOver, bedtools 2.31.1, Riboseek; see environments/.
"""


def request(url: str, method: str = "GET", headers: dict[str, str] | None = None):
    return urllib.request.urlopen(urllib.request.Request(
        url, method=method, headers={"User-Agent": USER_AGENT, **(headers or {})}), timeout=60)


def remote_size(url: str) -> int | None:
    """Size reported by the server (HEAD, or a one-byte range request when HEAD is refused)."""
    try:
        with request(url, "HEAD") as response:
            if response.headers.get("Content-Length"):
                return int(response.headers["Content-Length"])
    except urllib.error.HTTPError:
        pass
    with request(url, headers={"Range": "bytes=0-0"}) as response:
        total = re.search(r"/(\d+)$", response.headers.get("Content-Range", ""))
        if total:
            return int(total.group(1))
        length = response.headers.get("Content-Length")
        return int(length) if response.status == 200 and length else None


_MD5SUMS: dict[str, dict[str, str]] = {}


def expected_md5(item: Item) -> str | None:
    if item.md5 != "MD5SUMS":
        return item.md5
    directory, name = item.url.rsplit("/", 1)
    if directory not in _MD5SUMS:
        with request(f"{directory}/MD5SUMS") as response:
            _MD5SUMS[directory] = {line.split()[1]: line.split()[0]
                                   for line in response.read().decode().splitlines() if len(line.split()) == 2}
    return _MD5SUMS[directory].get(name)


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def human(n: float | None) -> str:
    if n is None:
        return "size unknown"
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return ""


def fetch(item: Item, dest: Path, retries: int = 3) -> None:
    """Download one file to dest/item.path, resuming a partial download; skip it if already complete."""
    target = dest / item.path
    size = remote_size(item.url)
    if target.exists() and (size is None or target.stat().st_size == size):
        print(f"have  {item.path} ({human(target.stat().st_size)})", flush=True)
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_name(target.name + ".part")
    started = time.monotonic()
    for attempt in range(1, retries + 1):
        offset = part.stat().st_size if part.exists() else 0
        if size is not None and offset > size:
            offset = 0
        if size is not None and offset == size:
            break
        resume = f", resuming at {human(offset)}" if offset else ""
        print(f"get   {item.url}\n   -> {item.path} ({human(size)}{resume})", flush=True)
        try:
            with request(item.url, headers={"Range": f"bytes={offset}-"} if offset else None) as response:
                if offset and response.status != 206:  # range ignored: start again
                    offset = 0
                last = time.monotonic()
                with part.open("ab" if offset else "wb") as handle:
                    for block in iter(lambda: response.read(CHUNK), b""):
                        handle.write(block)
                        if time.monotonic() - last > 30:
                            last = time.monotonic()
                            print(f"      {human(part.stat().st_size)} of {human(size)}", flush=True)
            if size is not None and part.stat().st_size != size:
                raise OSError(f"incomplete: {part.stat().st_size} of {size} bytes")
            break
        except (OSError, urllib.error.URLError) as error:
            if attempt == retries:
                raise
            print(f"      {error}; retrying", flush=True)
            time.sleep(5 * attempt)
    md5 = expected_md5(item)
    if md5 and md5sum(part) != md5:
        part.rename(part.with_suffix(".bad"))
        raise OSError(f"{item.path}: MD5 mismatch (kept as {part.with_suffix('.bad').name})")
    part.rename(target)
    rate = target.stat().st_size / max(time.monotonic() - started, 1e-3)
    print(f"done  {item.path} ({human(target.stat().st_size)}, {human(rate)}/s{', MD5 ok' if md5 else ''})", flush=True)


def encode_items(table: Path, datasets: list[str] | None) -> list[Item]:
    """ENCODE files listed in a table with a file_accession column (e.g. Supplementary Table S8 saved as TSV),
    or any ENCFF accessions found in a text file."""
    frame = pd.DataFrame()
    if table.suffix in (".xlsx", ".xls"):
        frame = pd.read_excel(table, sheet_name="Accessions")  # needs openpyxl
    elif table.suffix in (".tsv", ".csv"):
        frame = pd.read_csv(table, sep="," if table.suffix == ".csv" else "\t", dtype=str)
    if "file_accession" in frame:
        if datasets and "dataset" in frame:
            frame = frame[frame.dataset.isin(datasets)]
        accessions = frame.file_accession.dropna().astype(str).tolist()
    else:
        accessions = re.findall(r"ENCFF\d{3}[A-Z]{3}", table.read_text())
    accessions = [a for a in dict.fromkeys(accessions) if re.fullmatch(r"ENCFF\d{3}[A-Z]{3}", a)]
    return [Item("encode", f"{ENCODE}/files/{a}/", f"encode/{a}") for a in accessions]


def resolve_encode(item: Item) -> Item:
    """Download URL, file name and MD5 of an ENCODE file from its metadata."""
    with request(f"{item.url}?format=json&frame=object", headers={"Accept": "application/json"}) as response:
        meta = json.load(response)
    return Item("encode", ENCODE + meta["href"], f"encode/{meta['href'].rsplit('/', 1)[1]}", meta.get("md5sum"))


def extract_human_rfam_hits(dest: Path) -> None:
    """Significant Rfam 15.1 'full' hits on the GRCh38 primary chromosomes (the hits of the paper, which were read
    from the Rfam public database: full_region joined with genseq for UP000005640)."""
    families = {}
    with gzip.open(dest / "rfam/family.txt.gz", "rt", encoding="latin-1") as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            families[f[0]] = (f[1], f[3], f[18])  # rfam_id, description, type
    out_path = dest / "rfam/rfam15.1_human_GRCh38_hits.tsv"
    seen, rows = set(), 0
    with gzip.open(dest / "rfam/full_region.txt.gz", "rt", encoding="latin-1") as handle, out_path.open("w") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "start", "end", "strand", "rfam_acc", "rfam_id", "description", "family_type",
                         "bit_score", "evalue_score", "rfamseq_acc", "seq_start", "seq_end", "cm_start", "cm_end",
                         "truncated"])
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if f[1] not in GRCH38_CHROMOSOMES or f[9] != "full" or f[10] != "1" or line in seen:
                continue
            seen.add(line)
            a, b = int(f[2]), int(f[3])
            rfam_id, description, family_type = families.get(f[0], ("", "", ""))
            writer.writerow([GRCH38_CHROMOSOMES[f[1]], min(a, b) - 1, max(a, b), "+" if a <= b else "-", f[0],
                             rfam_id, description, family_type, f[4], f[5], f[1], f[2], f[3], f[6], f[7], f[8]])
            rows += 1
    print(f"wrote rfam/{out_path.name}: {rows:,} hits (the paper used 20,516 hits in 1,077 families)", flush=True)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Download the public inputs of the RLMHelix pipelines.")
    parser.add_argument("dest", type=Path, help="destination directory (e.g. data/raw)")
    parser.add_argument("--only", nargs="+", choices=DEFAULT_GROUPS + OPTIONAL_GROUPS,
                        help=f"groups to fetch (default: {' '.join(DEFAULT_GROUPS)}; "
                             f"optional: {' '.join(OPTIONAL_GROUPS)})")
    parser.add_argument("--dry-run", action="store_true", help="list URLs and destinations without downloading")
    parser.add_argument("--check", action="store_true", help="like --dry-run, but ask each server for the size")
    parser.add_argument("--list-manual", action="store_true", help="list inputs that must be obtained by hand")
    parser.add_argument("--panel", type=Path, default=PANEL,
                        help="genome panel TSV (default: data/genome_panel_55.tsv)")
    parser.add_argument("--encode-table", type=Path,
                        help="Supplementary Table S8 (.xlsx, needs openpyxl), a TSV/CSV with a file_accession "
                             "column, or a text file of ENCFF accessions; needed for --only encode")
    parser.add_argument("--encode-dataset", nargs="+",
                        help="keep rows of these 'dataset' values only (e.g. 'eCLIP peaks')")
    args = parser.parse_args(argv)
    if args.list_manual:
        print(MANUAL)
        return
    groups = args.only or DEFAULT_GROUPS
    items = static_items() + (ensembl_items(args.panel) if "ensembl" in groups else [])
    items = [i for i in items if i.group in groups]
    if "encode" in groups:
        if not args.encode_table:
            parser.error("--only encode needs --encode-table (see --list-manual)")
        if not args.encode_table.exists():
            parser.error(f"--encode-table: {args.encode_table} not found (Supplementary Table S8, .xlsx or its Accessions sheet as TSV)")
        items += encode_items(args.encode_table, args.encode_dataset)

    if args.dry_run or args.check:
        total = 0
        for item in items:
            line = f"{item.group:<10} {item.url}  ->  {args.dest / item.path}"
            if args.check:
                try:
                    if item.group == "encode":
                        item = resolve_encode(item)
                    size = remote_size(item.url)
                    total += size or 0
                    line = f"{item.group:<10} {human(size):>10}  {item.url}  ->  {args.dest / item.path}"
                except (OSError, urllib.error.URLError) as error:
                    line = f"{item.group:<10} {'FAILED':>10}  {item.url}  ({error})"
            print(line, flush=True)
        print(f"{len(items)} files" + (f", {human(total)} in total" if args.check else ""))
        return

    failed = []
    for item in items:
        try:
            fetch(resolve_encode(item) if item.group == "encode" else item, args.dest)
        except (OSError, urllib.error.URLError, KeyError, ValueError) as error:
            print(f"FAILED {item.url}: {error}", file=sys.stderr, flush=True)
            failed.append(item.url)
    if "rfam" in groups and not any("Rfam/15.1/database_files" in url for url in failed):
        extract_human_rfam_hits(args.dest)
    if failed:
        sys.exit(f"{len(failed)} downloads failed; run the same command again to resume")


if __name__ == "__main__":
    main()
