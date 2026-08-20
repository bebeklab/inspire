Here is your complete document with the new instructions formatted neatly and seamlessly appended as Sections 11 and 12 at the bottom:

```markdown
# INSPIRE Human Petal Network Replication

## 1. Scope and scientific objective

This repository implements the INSPIRE network-alignment workflow. The analysis compares human protein-interaction subnetworks, called petals, using network topology, Gene Ontology Biological Process similarity, and Pfam domain similarity. 

This repository is for scientific review only. 

This package is intended to support:

1. reconstruction of the SQLite analysis database from documented inputs;
2. execution of the C++ INSPIRE network-alignment program;
3. generation of all pairwise petal similarities;
4. validation with self-alignments and biological positive controls;
5. hierarchical clustering and threshold-sensitivity analysis;
6. auditing of input-network construction separately from alignment scoring.

## 2. Repository structure

```text
inspire-replication/
├── database/
│   ├── inspire.sqlite                 # generated runtime database; not committed
│   └── unmapped_genes.csv             # generated mapping audit
├── input/
│   ├── blossom_edge_list_7_26a.csv    # human petal-network input
│   └── go/
│       ├── go-basic.obo               # GO ontology input
│       └── HUMAN-uniprot.gaf.gz       # human UniProt GO annotations
├── reports/                           # generated outputs; normally not committed
├── scripts/
│   ├── build_inspire_database.py
│   ├── load_annotations.py
│   ├── build_go_scores.py
│   ├── run_all_vs_all.py
│   ├── run_inspire.sh
│   └── compare_go_algorithms.py
├── src/
│   ├── AlignedGeneEdge.cpp
│   ├── AlignedGeneEdge.h
│   ├── database.cpp
│   ├── database.h
│   ├── gene.cpp
│   ├── gene.h
│   ├── main.cpp
│   ├── petal.cpp
│   └── petal.h
├── README.md
├── README_REPLICATION.md
└── .gitignore

```

## 3. Analysis inputs

### 3.1 Human petal edge list

Expected file:

```text
input/blossom_edge_list_7_26a.csv

```

Required columns:

```text
gene1,gene2,path,score,petal

```

Interpretation:

* `gene1` and `gene2` identify an undirected protein-interaction edge.
* `path` preserves the source path from which the edge was obtained.
* `score` preserves the upstream path or interaction score.
* `petal` identifies the network to which the observation belongs.
* Repeated or reversed observations are retained for provenance but collapsed to unique undirected edges in the runtime network.
* Runtime interaction type is recorded as `pp` for compatibility with the legacy INSPIRE code.

The edge list is a derived scientific input. Its provenance, generation script, filtering rules, and checksum must be retained. A low INSPIRE score can arise from discordant petal construction even when two seed genes have a well-established biological relationship.

### 3.2 GO ontology

Expected file:

```text
input/go/go-basic.obo

```

The exact file version and SHA-256 checksum are recorded in the SQLite metadata during loading.

### 3.3 Human GO annotations

Expected file:

```text
input/go/HUMAN-uniprot.gaf.gz

```

The loader:

* retains human UniProt-centric annotations;
* uses Biological Process annotations for compatibility with INSPIRE;
* excludes annotations carrying the `NOT` qualifier;
* maps gene symbols and GAF synonyms to UniProt accessions;
* records source metadata and checksums;
* writes unmapped genes to a CSV audit file.

### 3.4 Pfam annotations

Pfam domains are retrieved for preferred UniProt accessions and cached in SQLite. Pfam retrieval requires network access. The database records completion and error status so interrupted runs can be resumed.

## 4. Software environment

The workflow has been tested on macOS with Apple Silicon.

Required software:

* Python 3
* SQLite 3 command-line client
* a C++11 compiler
* Homebrew SQLite development libraries on macOS
* R and Pandoc for rendering the clustering report

Required R packages:

```r
DBI
RSQLite
dplyr
tidyr
readr
stringr
ggplot2
igraph
pheatmap
patchwork
rmarkdown
knitr

```

Record exact versions for a frozen release:

```bash
python3 --version
sqlite3 --version
clang++ --version
R --version
Rscript -e 'sessionInfo()'

```

## 5. Build the runtime database

Run all commands from the repository root.

### 5.1 Create the network database

```bash
python3 scripts/build_inspire_database.py \
  input/blossom_edge_list_7_26a.csv \
  database/inspire.sqlite

```

### 5.2 Load GO annotations

```bash
python3 scripts/load_annotations.py \
  --db database/inspire.sqlite \
  --obo input/go/go-basic.obo \
  --gaf input/go/HUMAN-uniprot.gaf.gz \
  --skip-pfam \
  --unmapped-report reports/unmapped_genes.csv

```

### 5.3 Retrieve Pfam annotations

```bash
python3 scripts/load_annotations.py \
  --db database/inspire.sqlite \
  --pfam-only

```

This operation is resumable. Re-running it should skip accessions already marked complete.

### 5.4 Build GO pair scores

Use the algorithm identifier consistently throughout database generation, C++ execution, and R analysis:

```bash
python3 scripts/build_go_scores.py \
  --db database/inspire.sqlite \
  --algorithm inspire_legacy_v1

```

If the local script exposes different command-line options, inspect them before execution:

```bash
python3 scripts/build_go_scores.py --help

```

The active GO normalization must be calculated from the same GO-score cache used by the executable. In the validated current-human run:

```text
algorithm: inspire_legacy_v1
cached GO pairs: 5,995
GO normalization: 14.1026

```

Do not substitute the historical constant `23.7367` for this current-human cache. That combination produced self-alignment scores of `0.460829` instead of `1`.

## 6. Verify the database

### 6.1 Structural integrity

```bash
sqlite3 database/inspire.sqlite "PRAGMA integrity_check;"
sqlite3 database/inspire.sqlite "PRAGMA foreign_key_check;"

```

Expected results:

* integrity check returns `ok`;
* foreign-key check returns no rows.

### 6.2 Core counts

```bash
sqlite3 -header -column database/inspire.sqlite "
SELECT
  (SELECT COUNT(*) FROM petal) AS petals,
  (SELECT COUNT(*) FROM edge) AS unique_edges,
  (SELECT COUNT(*) FROM gene) AS genes,
  (SELECT COUNT(*) FROM gene_uniprot WHERE is_preferred = 1) AS preferred_uniprot,
  (SELECT COUNT(*) FROM go_term) AS go_terms,
  (SELECT COUNT(*) FROM go_parent) AS go_relationships,
  (SELECT COUNT(*) FROM protein_go_annotation) AS bp_annotations,
  (SELECT COUNT(*) FROM protein_pfam) AS protein_pfam_rows;
"

```

The validated analysis contains 27 petals. Exact annotation counts are version-dependent and should be reported from the generated metadata rather than hard-coded.

### 6.3 Metadata and checksums

```bash
sqlite3 -header -column database/inspire.sqlite \
  "SELECT key, value FROM metadata ORDER BY key;"

```

### 6.4 Mapping and Pfam audits

```bash
sqlite3 -header -column database/inspire.sqlite "
SELECT pfam_status, COUNT(*) AS proteins
FROM annotation_status
GROUP BY pfam_status
ORDER BY pfam_status;
"

```

```bash
sqlite3 -header -column database/inspire.sqlite "
SELECT uniprot_accession, last_error
FROM annotation_status
WHERE pfam_status = 'error'
ORDER BY uniprot_accession;
"

```

```bash
sqlite3 -header -column database/inspire.sqlite "
SELECT
  g.gene_symbol,
  SUM(CASE WHEN u.is_preferred = 1 THEN 1 ELSE 0 END) AS preferred_count
FROM gene AS g
LEFT JOIN gene_uniprot AS u
  ON u.gene_symbol = g.gene_symbol
GROUP BY g.gene_symbol
HAVING preferred_count != 1;
"

```

```bash
sqlite3 -header -column database/inspire.sqlite "
SELECT
  u.gene_symbol,
  u.uniprot_accession
FROM gene_uniprot AS u
LEFT JOIN protein_go_annotation AS a
  ON a.uniprot_accession = u.uniprot_accession
WHERE u.is_preferred = 1
GROUP BY u.gene_symbol, u.uniprot_accession
HAVING COUNT(a.go_id) = 0
ORDER BY u.gene_symbol;
"

```

## 7. Build INSPIRE

From `src/`:

```bash
cd src
SQLITE_PREFIX="$(brew --prefix sqlite)" && \
clang++ -std=c++11 -O3 -DNDEBUG \
  -Wall -Wextra -Wpedantic \
  -I"$SQLITE_PREFIX/include" \
  -L"$SQLITE_PREFIX/lib" \
  gene.cpp petal.cpp AlignedGeneEdge.cpp main.cpp \
  -lsqlite3 \
  -o INSPIRE
cd ..

```

The executable interface is:

```bash
./src/INSPIRE DATABASE PETAL_A PETAL_B ALGORITHM

```

Example:

```bash
./src/INSPIRE \
  database/inspire.sqlite \
  SMAD2 \
  SMAD3 \
  inspire_legacy_v1

```

## 8. Required pre-analysis validation

### 8.1 Self-alignments

Run several petals of different sizes:

```bash
./src/INSPIRE database/inspire.sqlite ABCA1 ABCA1 inspire_legacy_v1
./src/INSPIRE database/inspire.sqlite FN1 FN1 inspire_legacy_v1
./src/INSPIRE database/inspire.sqlite GLI3 GLI3 inspire_legacy_v1

```

Required result for each:

```text
norm_score: 1

```

A self-score different from 1 indicates a mismatch among GO scores, GO normalization, Pfam scoring, or final score conversion. Do not proceed to clustering until self-alignment passes.

### 8.2 Positive controls

Validated current-human scores:

```text
SMAD2-SMAD3   0.714090
SMAD2-SMAD4   0.720430
SMAD3-SMAD4   0.427758
NF1-KRAS      0.620096

```

## 9. Run all pairwise comparisons

Run the all-versus-all analysis with the current algorithm and write to a fresh output directory. A typical invocation is:

```bash
python3 scripts/run_all_vs_all.py \
  --jobs 15 \
  --algorithm inspire_legacy_v1 \
  --no-self

```

Do not mix logs produced with different GO normalization values. Preserve old results in a separately named archive directory.

Expected output layout:

```text
reports/all_vs_all/
├── logs/
└── summary_inspire_legacy_v1.csv

```

## 10. Cluster and visualize the similarity matrix

Primary clustering:

```r
hclust(as.dist(1 - inspire_similarity), method = "ward.D2")

```

Ward.D2 is used for clustering.

---

## 11. Custom execution pipeline & pseudogene filtering

Pseudogenes (e.g., `PMS2CL`) can introduce GO scoring anomalies during database preparation. Follow the steps below to filter pseudogenes and run the standard workflow pipeline:

### 11.1 Filter pseudogenes and build database

```bash
# Filter problematic pseudogenes (e.g., PMS2CL)
python3 pseudogenefilter.py

# Build SQLite database from edges
python3 scripts/build_inspire_database.py \
  networks/Blossom_edges_v1.csv \
  database/inspire.sqlite

# Load GO annotations
python3 scripts/load_annotations.py \
  --db database/inspire.sqlite \
  --obo input/go/go-basic.obo \
  --gaf input/go/HUMAN-uniprot.gaf.gz \
  --unmapped-report reports/unmapped_genes.csv

```

> **Note:** Check `reports/unmapped_genes.csv`. If unmapped pseudogenes are present, ensure they are added to `pseudogenefilter.py` and rerun the steps above.

### 11.2 Rebuild GO scores and run pairwise comparisons

```bash
# Rebuild GO scores
python3 scripts/build_go_scores.py \
  --db database/inspire.sqlite \
  --algorithm inspire-legacy \
  --rebuild

# Run all-versus-all comparisons
python3 scripts/run_all_vs_all.py \
  --database database/inspire.sqlite \
  --jobs 15 \
  --output-dir reports/all_vs_all

```

---

## 12. Final analysis rendering

After the pairwise runs are complete, execute the R Markdown scripts to finalize threshold grid clustering.

Run `STEP 4_Final_Inspire_Dendograms.Rmd` using the updated parameters:

```yaml
database: "../database/inspire_bestrule.sqlite"
inspire_log_dir: "../reports/all_vs_all_bestrule/logs"
output_dir: "../reports/clustering_bestrule"

```

```

```