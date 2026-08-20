# Basic Replication Steps: INSPIRE Petal Network Database

## Document information

- **Prepared:** July 27, 2026
- **Updated:** August 5, 2026
- **GO download-page access date:** July 27, 2026
- **Platform:** macOS on Apple Silicon, with later transfer to an HPC system supported
- **Runtime database:** SQLite
- **Network preparation and annotation:** Python 3
- **Network alignment:** C++
- **Default GO scoring algorithm:** `inspire_legacy_v1`

## Purpose

This document describes how to build a complete INSPIRE runtime database from a petal-network edge list, load current human Gene Ontology and Pfam annotations, compute the algorithm-specific GO score cache and normalization value, validate the database, compile INSPIRE, and run alignment tests.

The GO score-cache step is mandatory. Loading GO annotations alone does not create the `inspire_legacy_v1` pair-score cache or its normalization value. If `build_go_scores.py` is omitted, INSPIRE terminates with an error similar to:

```text
Could not find a normalization score for inspire_legacy_v1
```

## Source files and locations

Create one project directory. The examples below use:

```text
~/inspire-replication/
```

Recommended layout:

```text
~/inspire-replication/
├── scripts/
│   ├── build_inspire_database.py
│   ├── load_annotations.py
│   ├── build_go_scores.py
│   ├── run_all_vs_all.py
│   └── plot_inspire_alignment.py
├── src/
│   ├── main.cpp
│   ├── gene.cpp
│   ├── gene.h
│   ├── petal.cpp
│   ├── petal.h
│   ├── AlignedGeneEdge.cpp
│   ├── AlignedGeneEdge.h
│   └── base_header.h
├── input/
│   ├── blossom_edge_list_8_03a.csv
│   └── go/
│       ├── go-basic.obo
│       └── HUMAN-uniprot.gaf.gz
├── database/
│   └── inspire.sqlite
├── reports/
│   ├── unmapped_genes.csv
│   ├── alignment_logs/
│   ├── alignment_plots/
│   └── all_vs_all/
└── README.md
```

## Network source

- **Example file:** `blossom_edge_list_8_03a.csv`
- **Expected location:** `~/inspire-replication/input/blossom_edge_list_8_03a.csv`
- **Required columns:** `gene1`, `gene2`, `path`, `score`, `petal`
- **Network grouping:** Rows are grouped by the `petal` field.
- **Edge treatment:** Undirected.
- **Duplicate treatment:** Repeated or reversed gene pairs within the same petal are collapsed into one physical network edge. Original observations remain recorded in SQLite for provenance.
- **Interaction compatibility label:** Imported edges are labeled `pp`. The numeric upstream interaction or path score is stored separately.

A typical input begins as follows:

```csv
"gene1","gene2","path","score","petal"
"APC","CDC42","APC -> CDC42 -> ARHGEF12 -> ABCA1",0.0690807293982232,"ABCA1"
"CDC42","ARHGEF12","APC -> CDC42 -> ARHGEF12 -> ABCA1",0.0690807293982232,"ABCA1"
"ARHGEF12","ABCA1","APC -> CDC42 -> ARHGEF12 -> ABCA1",0.0690807293982232,"ABCA1"
"APC","KRAS","APC -> KRAS -> JAK2 -> ABCA1",0.0631976321268217,"ABCA1"
```

Inspect the actual input before building:

```bash
head input/blossom_edge_list_8_03a.csv
```

## GO ontology source

- **File:** `go-basic.obo`
- **Expected location:** `~/inspire-replication/input/go/go-basic.obo`
- **Source page:** https://geneontology.org/docs/download-ontology/
- **Page access date:** July 27, 2026
- **Acquisition:** Supplied or downloaded manually by the researcher.

## Human GO annotation source

- **File:** `HUMAN-uniprot.gaf.gz`
- **Expected location:** `~/inspire-replication/input/go/HUMAN-uniprot.gaf.gz`
- **Source page:** https://geneontology.org/docs/download-go-annotations/downloads/
- **Page access date:** July 27, 2026
- **Selection:** Human, UniProt-centric GAF.
- **Acquisition:** Supplied or downloaded manually by the researcher.
- **Compression:** The file may remain gzip-compressed.

## Software requirements

Verify the required command-line tools:

```bash
python3 --version
sqlite3 --version
clang++ --version
```

The database-construction, annotation, and GO-score scripts use the Python standard library. The alignment plotting script additionally requires NetworkX and Matplotlib.

Check plotting dependencies:

```bash
python3 -c "import networkx, matplotlib; print('Plot dependencies available')"
```

If needed:

```bash
python3 -m pip install networkx matplotlib
```

# Replication workflow

## Step 1: Create directories

```bash
mkdir -p ~/inspire-replication/{scripts,src,input/go,database,reports}
```

Copy the scripts, source files, petal edge list, GO ontology, and human GAF into the locations shown above.

Change to the project root:

```bash
cd ~/inspire-replication
```

## Step 2: Build the SQLite network database

```bash
python3 scripts/build_inspire_database.py \
  input/blossom_edge_list_8_03a.csv \
  database/inspire.sqlite
```

Expected output includes:

```text
Database written to: .../database/inspire.sqlite
```

Counts can differ from raw CSV row counts because repeated or reversed observations of the same undirected edge are collapsed into one runtime network edge.

### Alternative database example

For the all-rules network input:

```bash
python3 scripts/build_inspire_database.py \
  networks/Blossom_7_31_allrules_cutoff_v2.csv \
  database/inspire_allrules.sqlite
```

Every subsequent command must use the same database filename. Do not build `inspire_allrules.sqlite` and then accidentally run annotation or alignment commands against `inspire.sqlite`.

## Step 3: Load GO annotations

This first annotation run does not contact InterPro.

```bash
python3 scripts/load_annotations.py \
  --db database/inspire.sqlite \
  --obo input/go/go-basic.obo \
  --gaf input/go/HUMAN-uniprot.gaf.gz \
  --skip-pfam \
  --unmapped-report reports/unmapped_genes.csv
```

This step:

1. Loads GO terms and direct parent relationships from `go-basic.obo`.
2. Loads human GO annotations from `HUMAN-uniprot.gaf.gz`.
3. Retains Biological Process annotations, matching the historical INSPIRE implementation.
4. Excludes annotations carrying the `NOT` qualifier.
5. Maps network gene symbols and GAF synonyms to UniProt accessions.
6. Records source filenames, SHA-256 checksums, file-version headers, and load dates in SQLite.
7. Writes genes not mapped through the GAF to `reports/unmapped_genes.csv`.

For an alternative database, substitute its filename consistently:

```bash
python3 scripts/load_annotations.py \
  --db database/inspire_allrules.sqlite \
  --obo input/go/go-basic.obo \
  --gaf input/go/HUMAN-uniprot.gaf.gz \
  --skip-pfam \
  --unmapped-report reports/unmapped_genes_allrules.csv
```

## Step 4: Retrieve Pfam annotations

Pfam retrieval requires Internet access.

```bash
python3 scripts/load_annotations.py \
  --db database/inspire.sqlite \
  --pfam-only
```

The script queries InterPro for Pfam domains associated with preferred UniProt accessions. Results and completion status are cached in SQLite. Re-running the command skips accessions already marked complete.

For the all-rules database:

```bash
python3 scripts/load_annotations.py \
  --db database/inspire_allrules.sqlite \
  --pfam-only
```

If Pfam retrieval is intentionally postponed, INSPIRE can still be tested after GO-score generation, but the resulting scoring configuration is not equivalent to the complete GO-plus-Pfam analysis.

## Step 5: Build GO pair scores and normalization

**This step is mandatory and must be run separately for every newly built SQLite database.**

Loading the ontology and protein annotations does not populate the algorithm-specific GO pair-score cache. INSPIRE requires both:

- cached GO pair scores for the selected algorithm; and
- a normalization value computed from that same score cache.

First inspect the current script interface:

```bash
python3 scripts/build_go_scores.py --help
```

Using the supported option names, run:

```bash
python3 scripts/build_go_scores.py \
  --db database/inspire.sqlite \
  --algorithm inspire-legacy
```

For the all-rules database:

```bash
python3 scripts/build_go_scores.py \
  --db database/inspire_allrules.sqlite \
  --algorithm inspire-legacy
```

Do not copy a normalization constant manually from another database. The normalization must be derived from the same GO score cache stored in the database being analyzed.

The previously validated current-human database used:

```text
Algorithm: inspire_legacy_v1
Cached GO pair scores: 5,995
GO normalization: approximately 14.1026
```

A newly constructed network database may require a different number of cached pairs or a different normalization if it introduces proteins or GO-term combinations not present in the earlier database.

### Confirm that GO scores were created

```bash
sqlite3 database/inspire.sqlite ".dump" \
  | grep -n "inspire_legacy_v1" \
  | head -30
```

This command must return rows. If it returns nothing, the GO score cache and normalization were not stored in that database.

Count occurrences:

```bash
sqlite3 database/inspire.sqlite ".dump" \
  | grep -c "inspire_legacy_v1"
```

The count must be greater than zero.

For the all-rules database:

```bash
sqlite3 database/inspire_allrules.sqlite ".dump" \
  | grep -n "inspire_legacy_v1" \
  | head -30
```

## Step 6: Verify the completed database

### Inspect metadata and recorded versions

```bash
sqlite3 -header -column database/inspire.sqlite \
  "SELECT key, value FROM metadata ORDER BY key;"
```

### Count networks and edges

```bash
sqlite3 -header -column database/inspire.sqlite \
  "SELECT petal, edge_count, gene_count FROM petal ORDER BY petal;"
```

### Count UniProt mappings

```bash
sqlite3 -header -column database/inspire.sqlite \
  "SELECT COUNT(*) AS mapping_rows FROM gene_uniprot;"
```

### Count Biological Process annotations

```bash
sqlite3 -header -column database/inspire.sqlite \
  "SELECT COUNT(*) AS bp_annotations FROM protein_go_annotation;"
```

### Count Pfam annotations

```bash
sqlite3 -header -column database/inspire.sqlite \
  "SELECT COUNT(*) AS pfam_annotations FROM protein_pfam;"
```

### Review unmapped genes

```bash
cat reports/unmapped_genes.csv
```

### Confirm algorithm records

```bash
sqlite3 database/inspire.sqlite ".dump" \
  | grep "inspire_legacy_v1" \
  | head
```

## Step 7: Compile INSPIRE

From the repository root:

```bash
cd src
```

Production build on macOS with Homebrew SQLite:

```bash
SQLITE_PREFIX="$(brew --prefix sqlite)" && \
clang++ -std=c++11 -O3 -DNDEBUG \
  -Wall -Wextra -Wpedantic \
  -I"$SQLITE_PREFIX/include" \
  -L"$SQLITE_PREFIX/lib" \
  gene.cpp petal.cpp AlignedGeneEdge.cpp main.cpp \
  -lsqlite3 \
  -o INSPIRE
```

Return to the repository root:

```bash
cd ..
```

## Step 8: Test normalization loading with one comparison

Before any parallel analysis, run one pair:

```bash
./src/INSPIRE \
  database/inspire.sqlite \
  SMAD4 \
  TAF2 \
  inspire_legacy_v1 \
  reports/alignment_logs/SMAD4__TAF2__inspire_legacy_v1.tsv
```

The startup output must report a normalization and cache for the requested algorithm, for example:

```text
GO normalization for inspire_legacy_v1: 14.1026
GO normalization in use: 14.1026
Loaded 5995 GO pair scores for inspire_legacy_v1
```

The values for a newly built database may differ, but the requested algorithm must be found and the reported normalization values must agree.

If the program reports:

```text
Could not find a normalization score for inspire_legacy_v1
```

return to Step 5 and build the GO score cache for the exact database passed to INSPIRE.

## Step 9: Validate self-alignment

Test several petals of different sizes:

```bash
./src/INSPIRE \
  database/inspire.sqlite \
  ABCA1 \
  ABCA1 \
  inspire_legacy_v1 \
  reports/alignment_logs/ABCA1__ABCA1__inspire_legacy_v1.tsv
```

```bash
./src/INSPIRE \
  database/inspire.sqlite \
  FN1 \
  FN1 \
  inspire_legacy_v1 \
  reports/alignment_logs/FN1__FN1__inspire_legacy_v1.tsv
```

```bash
./src/INSPIRE \
  database/inspire.sqlite \
  GLI3 \
  GLI3 \
  inspire_legacy_v1 \
  reports/alignment_logs/GLI3__GLI3__inspire_legacy_v1.tsv
```

Each self-alignment should return:

```text
norm_score: 1
```

Do not begin all-versus-all analysis if self-alignment does not equal 1.

## Step 10: Validate direction symmetry

Run a pair in both directions:

```bash
./src/INSPIRE \
  database/inspire.sqlite \
  KRAS \
  NF1 \
  inspire_legacy_v1 \
  reports/alignment_logs/KRAS__NF1__inspire_legacy_v1.tsv
```

```bash
./src/INSPIRE \
  database/inspire.sqlite \
  NF1 \
  KRAS \
  inspire_legacy_v1 \
  reports/alignment_logs/NF1__KRAS__inspire_legacy_v1.tsv
```

The normalized scores should agree apart from insignificant printed rounding.

## Step 11: Run all-versus-all analysis

For 27 petals, excluding self-comparisons produces 351 unordered comparisons:

```text
27 × 26 / 2 = 351
```

Recommended initial command:

```bash
python3 scripts/run_all_vs_all.py \
  --database database/inspire.sqlite \
  --jobs 6 \
  --output-dir reports/all_vs_all \
  --algorithm inspire_legacy_v1 \
  --no-self \
  --network-gap 4.5
```

For the all-rules database:

```bash
python3 scripts/run_all_vs_all.py \
  --database database/inspire_allrules.sqlite \
  --jobs 6 \
  --output-dir reports/all_vs_all_allrules \
  --algorithm inspire_legacy_v1 \
  --no-self \
  --network-gap 4.5
```

Start with a conservative worker count because large network pairs can allocate substantial AGE matrices. Increase `--jobs` only after observing acceptable memory pressure.

The updated runner writes:

```text
reports/all_vs_all/
├── logs/
├── alignments/
├── alignment_plots/
├── summary_inspire_legacy_v1.csv
├── failed_jobs.csv
└── failed_plots.csv
```

# Consolidated verification command

Run this only after completing Steps 2 through 5.

```bash
cd ~/Box/Research/INSPIRE/inspire-replication

DATABASE="database/inspire.sqlite"
ALGORITHM="inspire_legacy_v1"

echo "=== SQLite integrity ==="
sqlite3 "$DATABASE" "PRAGMA integrity_check;"

echo
echo "=== Foreign-key problems, should be empty ==="
sqlite3 "$DATABASE" "PRAGMA foreign_key_check;"

echo
echo "=== Core counts ==="
sqlite3 -header -column "$DATABASE" "
SELECT
    (SELECT COUNT(*) FROM petal) AS petals,
    (SELECT COUNT(*) FROM edge) AS unique_edges,
    (SELECT COUNT(*) FROM gene) AS genes,
    (SELECT COUNT(*) FROM gene_uniprot WHERE is_preferred=1)
        AS preferred_uniprot,
    (SELECT COUNT(*) FROM go_term) AS go_terms,
    (SELECT COUNT(*) FROM go_parent) AS go_relationships,
    (SELECT COUNT(*) FROM protein_go_annotation) AS bp_annotations,
    (SELECT COUNT(*) FROM protein_pfam) AS protein_pfam_rows;
"

echo
echo "=== Pfam status ==="
sqlite3 -header -column "$DATABASE" "
SELECT pfam_status, COUNT(*) AS proteins
FROM annotation_status
GROUP BY pfam_status
ORDER BY pfam_status;
"

echo
echo "=== Pfam errors ==="
sqlite3 -header -column "$DATABASE" "
SELECT uniprot_accession, last_error
FROM annotation_status
WHERE pfam_status='error'
ORDER BY uniprot_accession;
"

echo
echo "=== Genes without exactly one preferred accession ==="
sqlite3 -header -column "$DATABASE" "
SELECT
    g.gene_symbol,
    SUM(CASE WHEN u.is_preferred=1 THEN 1 ELSE 0 END) AS preferred_count
FROM gene AS g
LEFT JOIN gene_uniprot AS u
    ON u.gene_symbol=g.gene_symbol
GROUP BY g.gene_symbol
HAVING preferred_count != 1;
"

echo
echo "=== Preferred proteins without GO Biological Process annotations ==="
sqlite3 -header -column "$DATABASE" "
SELECT
    u.gene_symbol,
    u.uniprot_accession
FROM gene_uniprot AS u
LEFT JOIN protein_go_annotation AS a
    ON a.uniprot_accession=u.uniprot_accession
WHERE u.is_preferred=1
GROUP BY u.gene_symbol, u.uniprot_accession
HAVING COUNT(a.go_id)=0
ORDER BY u.gene_symbol;
"

echo
echo "=== GO score algorithm records ==="
ALGORITHM_COUNT=$(
    sqlite3 "$DATABASE" ".dump" \
      | grep -c "$ALGORITHM" \
      || true
)
echo "$ALGORITHM records: $ALGORITHM_COUNT"

if [ "$ALGORITHM_COUNT" -eq 0 ]; then
    echo "ERROR: GO pair scores and normalization are missing for $ALGORITHM"
    echo "Run scripts/build_go_scores.py for this database before INSPIRE."
    exit 1
fi

echo
echo "=== WAL checkpoint ==="
sqlite3 "$DATABASE" "PRAGMA wal_checkpoint(FULL);"

echo
echo "=== Verification finished ==="
```

To verify another database, change only:

```bash
DATABASE="database/inspire_allrules.sqlite"
```

# Preserve the replication package

Do not replace the network edge list, GO ontology, or human GAF after building the annotations and GO score cache without rebuilding the affected database layers.

Record checksums:

```bash
shasum -a 256 input/blossom_edge_list_8_03a.csv
shasum -a 256 input/go/go-basic.obo
shasum -a 256 input/go/HUMAN-uniprot.gaf.gz
shasum -a 256 database/inspire.sqlite
```

Export metadata:

```bash
sqlite3 -header -csv database/inspire.sqlite \
  "SELECT key, value FROM metadata ORDER BY key;" \
  > reports/database_metadata.csv
```

Record software versions:

```bash
python3 --version
sqlite3 --version
clang++ --version
```

# Current endpoint

At the end of the complete workflow, `database/inspire.sqlite` contains:

- petal networks and unique undirected edges;
- original edge observations and upstream scores;
- GO ontology terms and parent relationships;
- human Biological Process annotations;
- gene-to-UniProt mappings derived from the GAF;
- Pfam annotations retrieved from InterPro;
- algorithm-specific GO pair scores for `inspire_legacy_v1`;
- the GO normalization value computed from that score cache;
- metadata and checksums needed to identify the input versions.

The final runtime database is not complete for INSPIRE until both the algorithm-specific GO score cache and normalization record have been generated.

# Rebuild dependency order

If the network input or database is rebuilt, rerun the dependent stages in this order:

```text
Petal edge list
    ↓
build_inspire_database.py
    ↓
load_annotations.py with GO ontology and GAF
    ↓
load_annotations.py --pfam-only
    ↓
build_go_scores.py
    ↓
INSPIRE self-tests and symmetry tests
    ↓
run_all_vs_all.py
    ↓
clustering and figures
```
