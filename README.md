# pfam-tools

Python toolkit for Pfam domain annotation of prokaryotic genomes: scan proteins against Pfam-A, build QC-filtered presence/absence matrices, and annotate with GO / GO-Slim terms.

## Pipeline

```
Protein FASTAs (.faa)
        │
        ▼
  run_hmmscan.py                      scan against Pfam-A (pyhmmer)
        │                             produces one .domtblout per .faa
        ├──────────────────────────────────────────────────┐
        ▼                                                  ▼
  pfam_scan.py                        pfam_neighbourhood_extractor.py
        │  QC filter → presence/          extract genomic neighbourhoods
        │  absence matrix                 around target Pfam domains
        ▼                                 (+ GenBank files as input)
  pfam_mapping.py                         → neighbourhood .gbk + .faa
        annotate Pfam IDs
        with GO & GO-Slim
```

## Installation

```bash
pip install -r requirements.txt
```

## Quickstart

```bash
# 1. scan
python run_hmmscan.py \
    --pfam Pfam-A.hmm \
    --faa-dir genomes/ \
    --output-dir hmm_results/ \
    --processes 4 --threads 4

# 2. matrix
python pfam_scan.py \
    --domtbls hmm_results/ \
    --out pfam_matrix.tsv \
    --evalue 1e-3 --cov 0.5 \
    --resolve-overlaps --cpu 4

# 3. annotate
tail -n+2 pfam_matrix.tsv | cut -f1 > pfams.txt
python pfam_mapping.py \
    -i pfams.txt \
    -o annotation_results/ \
    -d db/
```

Step 2 outputs a matrix with Pfam IDs as the row index.
Step 3 expects a plain-text file with Pfam IDs — extract them with `tail | cut` as shown above.

---

## run_hmmscan.py

Scan protein `.faa` files against a pressed Pfam-A HMM database using pyhmmer.
Produces one `.domtblout` per input file.

```bash
python run_hmmscan.py \
    --pfam Pfam-A.hmm \
    --faa-dir genomes/ \
    --output-dir hmm_results/ \
    --processes 4 --threads 4
```

| Flag | Description |
|------|-------------|
| `--pfam` | Pressed Pfam-A HMM database (`.h3m`/`.h3f`/…) |
| `--faa-dir` | Directory with `.faa` files (default: `.`) |
| `--output-dir` | Where to write `.domtblout` files (default: same as `--faa-dir`) |
| `--processes` | Worker processes (default: 1) |
| `--threads` | Threads per process; 0 = auto (default: 0) |
| `--stream` | Stream HMMs from disk instead of loading into RAM |

---

## pfam_scan.py

Parse `.domtblout` files, apply QC filters, and build a Pfam matrix.

```bash
python pfam_scan.py \
    --domtbls hmm_results/ \
    --out pfam_matrix.tsv \
    --evalue 1e-3 --cov 0.5 \
    --resolve-overlaps --cpu 4
```

| Flag | Description |
|------|-------------|
| `--domtbls` | Directory of `.domtblout` files |
| `--out` | Output TSV (default: `pfam_presence.tsv`) |
| `--mode` | `hmmscan` or `hmmsearch` (default: `hmmscan`) |
| `--evalue` | Max E-value (default: 1e-3) |
| `--cov` | Min HMM coverage (default: 0.5) |
| `--resolve-overlaps` | Pfam-style overlap resolution |
| `--per-sequence` | One column per protein instead of per file |
| `--counts` | Raw domain counts instead of binary 0/1 |
| `--log-transform` | Log-transform counts (requires `--counts`) |
| `--use-ga` | Use GA thresholds (requires `--ga-file`) |
| `--cpu` | Parallel workers (default: 1) |

**Default output** (per file / genome):
```
pfam        genome_A    genome_B
PF00001     1           0
PF00005     1           1
```

**With `--per-sequence`** (per protein):
```
pfam        proteinX    proteinY
PF00001     1           0
PF00005     0           1
```

---

## pfam_mapping.py

Map Pfam IDs to GO terms and GO-Slim functional categories.

### Single Pfam list

Input: a plain-text file containing Pfam IDs (e.g. `PF00001`, one per line or mixed in free text).

```bash
python pfam_mapping.py -i pfams.txt -o results/ -d db/
```

### Multiple feature sets

Input: a two-column TSV (`set_name` and `features`, where features are comma-separated Pfam IDs):

```
set_name	features
clade_A	PF00001,PF00002,PF00013
clade_B	PF00001,PF00005,PF00069
```

```bash
python pfam_mapping.py --feature-sets sets.tsv -o results/ -d db/
```

| Flag | Description |
|------|-------------|
| `-i` | Plain-text Pfam ID file (default: `pfams.txt`) |
| `-o` | Output directory (default: `.`) |
| `-d` | Cache directory for database files (default: `.`) |
| `--pfam-tsv` | Path to `Pfam-A.clans.tsv.gz` (auto-downloaded if absent) |
| `--feature-sets` | TSV for batch mode |

### Database files

Downloaded automatically on first run into `-d`:

| File | Source |
|------|--------|
| `Pfam-A.clans.tsv.gz` | InterPro FTP |
| `go-basic.obo` | Gene Ontology |
| `goslim_generic.obo` | Gene Ontology |

### Outputs

**Per feature set** (or single list):

| File | Content |
|------|---------|
| `pfam_descriptions_only.tsv` | All Pfam IDs with descriptions |
| `pfam_go_annotation.tsv` | Pfam → GO → GO-Slim mappings |
| `pfam_goslim_summary.tsv` | Domain counts per GO-Slim category |
| `pfam_goslim_explanation.txt` | GO-Slim definitions |
| `pfam_goslim_summary.png` | Bar plot + annotation coverage pie |

**Cross-set** (when `--feature-sets` has ≥ 2 sets):

| File | Content |
|------|---------|
| `go_terms_universality.tsv` | GO term presence across sets |
| `goslim_terms_universality.tsv` | GO-Slim presence across sets |

---

## pfam_neighbourhood_extractor.py

Extract genomic neighbourhoods (N genes upstream + target + N genes downstream) around target Pfam domains. Requires `.domtblout` files (from step 1) and corresponding GenBank files.

```bash
python pfam_neighbourhood_extractor.py \
    --domtbls hmm_results/ \
    --genbanks genbank_files/ \
    --pfam PF00001,PF00002 \
    --output neighbourhoods/ \
    --upstream 5 --downstream 5 \
    --evalue 1e-3 --coverage 0.5 \
    --no-overlaps \
    --create-faa \
    --cpu 4 -v
```

The `--pfam` argument accepts a single ID, comma-separated IDs, or a path to a file with one ID per line.

| Flag | Description |
|------|-------------|
| `--domtbls` | Directory of `.domtblout` files |
| `--genbanks` | Directory of GenBank files (`.gbk`, `.gbff`, `.gb`, `.genbank`) |
| `--pfam` | Target Pfam ID(s): single, comma-separated, or file |
| `--output` | Output directory (default: `neighbourhoods`) |
| `--domtbl-ext` | Domtblout file extension (default: `.domtblout`) |
| `--upstream` | Genes upstream of target to include (default: 5) |
| `--downstream` | Genes downstream of target to include (default: 5) |
| `--mode` | `hmmscan` or `hmmsearch` (default: `hmmscan`) |
| `--evalue` | Max E-value threshold (default: 1e-3) |
| `--coverage` | Min HMM coverage (default: 0.5) |
| `--no-overlaps` | Apply Pfam-style overlap resolution |
| `--min-genes` | Min genes required per neighbourhood (default: 0) |
| `--use-ga` | Use GA thresholds (requires `--ga-file`) |
| `--ga-file` | Path to Pfam-A.hmm for GA thresholds |
| `--create-faa` | Write `.faa` protein FASTA per Pfam domain |
| `--cpu` | Parallel workers (default: 1) |

### Outputs

```
neighbourhoods/
├── PF00001/
│   ├── genomeA_PF00001_locusX_neighbourhood.gbk
│   └── genomeB_PF00001_locusY_neighbourhood.gbk
├── PF00002/
│   └── ...
└── faa_files/          (with --create-faa)
    ├── PF00001.faa
    └── PF00002.faa
```

Each neighbourhood `.gbk` file contains the target gene and its flanking genes with adjusted coordinates. The target gene is annotated with a `TARGET: Contains PFxxxxx domain` note; flanking genes are annotated with their relative position.

The `.faa` files contain deduplicated protein sequences for all genes found in neighbourhoods of each Pfam domain.

---

## Repository structure

```
pfam-tools/
├── run_hmmscan.py                    # step 1: HMMER scanning
├── pfam_scan.py                      # step 2: QC filtering & matrix
├── pfam_mapping.py                   # step 3: GO / GO-Slim annotation
├── pfam_neighbourhood_extractor.py   # neighbourhood extraction
├── requirements.txt
└── .gitignore
```

## Dependencies

- **pyhmmer** ≥ 0.11.1 — HMMER3 Python bindings
- **pandas**, **numpy** — data handling
- **goatools** — GO parsing and slim mapping
- **pfam2go** — Pfam-to-GO mappings
- **matplotlib** — plots
- **biopython** ≥ 1.81 — GenBank parsing (neighbourhood extraction)
- **requests**, **tqdm** — downloads and progress bars

