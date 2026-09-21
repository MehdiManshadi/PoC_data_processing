# PoC Data Processing

A proteomics data-processing pipeline that discovers public **PRIDE** LFQ
(label-free quantification) proteomics datasets for human, mouse, and
zebrafish, scores their gene coverage, aggregates each dataset's
quantitative values, and merges everything into a single per-gene feature
matrix for downstream analysis.

The pipeline is a sequence of lettered scripts under [`scripts/`](scripts/),
where each stage consumes the previous stage's output:

```
a → discover PRIDE datasets
b → score gene coverage per dataset
c → aggregate the Human Protein Atlas (HPA) reference dataset
d → download and aggregate LFQ values per PRIDE dataset
e → merge every aggregated dataset into one feature matrix
```

## Repository layout

```
scripts/
├── a_pride_dataset_finder.py          # Stage a: dataset discovery
├── b_dataset_coverage.py              # Stage b: gene coverage scoring
├── c_aggrigate_HPA_dataset.py         # Stage c: HPA gene×tissue aggregation
├── d_pride_lfq_functions_proteins.py  # Stage d: PRIDE download + LFQ aggregation helpers
├── d_create_processed_lfq_protein.ipynb  # Stage d: notebook driving the above, per dataset
└── e_features_collector.py            # Stage e: final feature-matrix merge
```

Each stage writes its output into its own subfolder under `scripts/`
(e.g. `a_found_potential_datasets/`, `b_coverage_reports/`,
`c&d_aggrigated_datasets/`, `e_processed_features/`). These output folders,
along with all generated CSV/XLSX/TXT/GZ data, are intentionally excluded
from version control (see [`.gitignore`](.gitignore)) — this repository
tracks the processing code, not the data it produces.

The working tree may also contain local, non-tracked supporting tooling
(`gene_orthology/`, `isoforms/`, `selected_features/`) used for
human→zebrafish/mouse ortholog mapping and per-study feature selection. These
are excluded from git because they hold large reference downloads and
per-study working files, not because they are unimportant to the workflow.

## Setup

Requires Python 3.10+.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install pandas numpy requests remotezip python-dotenv openpyxl jupyter
```

Some scripts read local paths from a `.env` file (not tracked in git):

```
OUTPUT_DIR=/absolute/path/to/PRIDE_downloads
HYPERGLYCEMIA_FEATURES_FILE=/absolute/path/to/features.txt
MITOCHONDRIAL_MYOPATHY_FEATURES_FILE=/absolute/path/to/features.txt
```

There is no automated test suite; each stage is validated by running it
against real data and inspecting the resulting CSV/XLSX output.

## Usage

Most stages are **not** general-purpose CLIs — several take their inputs
from constants defined near the top of the file (or from a notebook cell)
rather than command-line flags. Edit those in place before running.

### a — Discover PRIDE datasets

CLI-driven (`--help` for the full flag list):

```bash
python scripts/a_pride_dataset_finder.py --organism zebrafish \
    --test --test-limit 5 --output zebrafish_test.csv

python scripts/a_pride_dataset_finder.py --organism mouse \
    --require-lfq-columns --output mouse_lfq_report.csv
```

It queries OmicsDI, keeps only PRIDE datasets whose manifest contains an
exact `proteinGroups.txt` (checked remotely, including inside ZIP archives,
without downloading full archives), enriches results with PRIDE/OmicsDI/SDRF
metadata and Europe PMC abstracts, and reports tissue/sample information.

### b — Score gene coverage

Set `INPUT_FILE` (a report from stage a) and `OUTPUT_DIR` near the top of
`scripts/b_dataset_coverage.py`, then run it directly:

```bash
python scripts/b_dataset_coverage.py
```

It downloads each dataset's `proteinGroups.txt`, and counts a row only when
its LFQ values are non-zero/non-missing in more than 10% of that table's LFQ
columns (`PRESENCE_THRESHOLD`). Output is an `.xlsx` report with per-dataset
and per-table coverage sheets.

### c — Aggregate the HPA reference dataset

CLI-driven:

```bash
python scripts/c_aggrigate_HPA_dataset.py input_tissue_intensities.tsv \
    H1040_HPA.csv --threshold 0.10
```

Converts replicate-level `Gene`/`Gene name`/`Tissue`/`Intensity` rows into a
gene-by-tissue matrix, applying the same 10%-presence rule as the rest of the
pipeline (a gene passes only when more than 10% of all its intensity values
across every tissue are measured and non-zero).

### d — Download and aggregate PRIDE LFQ values

Driven from `scripts/d_create_processed_lfq_protein.ipynb`, which calls into
`scripts/d_pride_lfq_functions_proteins.py`:

```python
tables, catalog = prepare_lfq_selection(DATASET)
preview_lfq_range_groups(GROUPS, catalog)
result = process_lfq_ranges_and_show(DATASET, GROUPS, SETTINGS, tables, catalog)
```

`GROUPS` maps a biological group name to inclusive, 1-based LFQ column
ranges from the dataset's numbered LFQ catalogue; `SETTINGS` controls the
aggregation function (`mean`/`median`/`min`/`max`/`sum`), the presence
threshold, and the output report name. The notebook can also batch-process
every dataset listed in a `dataset_list.xlsx` (as produced by stage b).

### e — Merge into one feature matrix

```bash
python scripts/e_features_collector.py
```

Edit the `collect_features(...)` call at the bottom of the file to point at
your feature list, UniProt gene/protein mapping file, the aggregated-data
folder from stages c/d, and the ortholog tables used for non-human datasets.
For each aggregated CSV it: drops contaminant/reverse/site-only rows,
normalizes LFQ columns by subtracting each column's log2 median (zeros are
kept as zero, not treated as missing), maps proteins to genes — directly via
UniProt for human datasets, via ortholog tables for other organisms, and
directly via the existing `Gene name` column for the HPA reference dataset
— aggregates multiple proteins per gene by median, and concatenates every
dataset into one row-aligned `all_processed_features.csv`.

## Key conventions

- **Dataset identity** is threaded through filename prefixes (e.g. a file
  starting with `H`, `M`, or `Z` for human/mouse/zebrafish), not a column.
- **Zero vs. missing** is handled deliberately throughout: raw zeros are
  excluded from log2/median calculations but restored afterward as `0`
  rather than left as `NaN`.
- A shared **10%-presence threshold** is the quality gate used consistently
  across discovery, coverage, HPA aggregation, and LFQ aggregation — check
  every stage's threshold setting when changing this rule, not just one.
- Remote PRIDE ZIP archives are inspected via central-directory reads and
  bounded-prefix reads, never fully downloaded.
