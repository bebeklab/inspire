# Basic Replication Steps: INSPIRE Petal Network Database

## Document information

- **Prepared:** July 27, 2026
- **GO download-page access date:** July 27, 2026
- **Platform:** macOS on Apple Silicon, with later transfer to an HPC system supported
- **Runtime database:** SQLite
- **Network preparation:** Python 3
- **Network alignment:** C++ in the next stage

## Source files and locations

Create one project directory. The examples below use:

```text
~/inspire-replication/
```

Place the files in these locations:

```text
~/inspire-replication/
├── scripts/
│   ├── build_inspire_database.py
│   └── load_annotations.py
├── input/
│   ├── blossom_edge_list_7_22.csv
│   └── go/
│       ├── go-basic.obo
│       └── HUMAN-uniprot.gaf.gz
├── database/
│   └── inspire.sqlite
├── reports/
│   └── unmapped_genes.csv
└── README.md



### Network source

- **File:** `blossom_edge_list_7_22.csv`
- **Expected location:** `~/inspire-replication/input/blossom_edge_list_7_22.csv`
- **Required columns:** `gene1`, `gene2`, `path`, `score`, `petal`
- **Network grouping:** Rows are grouped by the `petal` field.
- **Edge treatment:** Undirected.
- **Duplicate treatment:** Repeated gene pairs within the same petal are collapsed into one network edge. Original observations remain recorded in SQLite.
- **Interaction compatibility label:** All imported edges are labeled `pp`. The new numeric PPI/path score is stored separately.

### GO ontology source

- **File:** `go-basic.obo`
- **Expected location:** `~/inspire-replication/input/go/go-basic.obo`
- **Source page:** https://geneontology.org/docs/download-ontology/
- **Page access date:** July 27, 2026
- **Acquisition:** Supplied or downloaded manually by the researcher.

### Human GO annotation source

- **File:** `HUMAN-uniprot.gaf.gz`
- **Expected location:** `~/inspire-replication/input/go/HUMAN-uniprot.gaf.gz`
- **Source page:** https://geneontology.org/docs/download-go-annotations/downloads/
- **Page access date:** July 27, 2026
- **Selection:** Human, UniProt-centric GAF.
- **Acquisition:** Supplied or downloaded manually by the researcher.
- **Compression:** The file may remain gzip-compressed.

## Software requirements

Verify Python 3:

```bash
python3 --version
```

The scripts use the Python standard library. No `pip install` step is required.



## Step 1: Create directories

```bash
mkdir -p ~/inspire-replication/{scripts,input/go,database,reports}
```

Copy the files into the locations shown above.


## Step 2: Build the SQLite network database


### edge list
The expeccted edge list looks like this: 

> head input/blossom_edge_list_7_22.csv
"gene1","gene2","path","score","petal"
"APC","CDC42","APC -> CDC42 -> ARHGEF12 -> ABCA1",0.0690807293982232,"ABCA1"
"CDC42","ARHGEF12","APC -> CDC42 -> ARHGEF12 -> ABCA1",0.0690807293982232,"ABCA1"
"ARHGEF12","ABCA1","APC -> CDC42 -> ARHGEF12 -> ABCA1",0.0690807293982232,"ABCA1"
"APC","KRAS","APC -> KRAS -> JAK2 -> ABCA1",0.0631976321268217,"ABCA1"
"KRAS","JAK2","APC -> KRAS -> JAK2 -> ABCA1",0.0631976321268217,"ABCA1"
"JAK2","ABCA1","APC -> KRAS -> JAK2 -> ABCA1",0.0631976321268217,"ABCA1"
....
....
....



```bash
cd ~/inspire-replication

python3 scripts/build_inspire_database.py   input/blossom_edge_list_7_26a.csv   database/inspire.sqlite
```

Expected output includes:

```text
Database written to: inspire-replication/database/inspire.sqlite

```

Counts can differ from raw row counts when the CSV contains repeated or reversed copies of the same undirected edge.


 


## Step 3: Load GO annotations only

This is the recommended first annotation run. It does not contact InterPro.

```bash
cd ~/inspire-replication

python3 scripts/load_annotations.py   --db database/inspire.sqlite   --obo input/go/go-basic.obo   --gaf input/go/HUMAN-uniprot.gaf.gz   --skip-pfam   --unmapped-report reports/unmapped_genes.csv
```

This performs the following actions:

1. Loads GO terms and direct parent relationships from `go-basic.obo`.
2. Loads human GO annotations from `HUMAN-uniprot.gaf.gz`.
3. Retains Biological Process annotations only, matching the historical C++ implementation.
4. Excludes annotations carrying the `NOT` qualifier.
5. Maps network gene symbols and GAF synonyms to UniProt accessions.
6. Records input filenames, SHA-256 checksums, file version headers, and load dates in SQLite.
7. Writes genes not mapped through the GAF to `reports/unmapped_genes.csv`.

## Step 4:  retrieve Pfam annotations

Pfam retrieval  requires Internet access.

After GO loading, run:

```bash
cd ~/inspire-replication

python3 scripts/load_annotations.py   --db database/inspire.sqlite   --pfam-only
```

The script queries InterPro for Pfam domains associated with preferred UniProt accessions. Results and completion status are cached in SQLite. Re-running the command skips completed accessions.

To postpone Pfam entirely, stop after Step 3.

## Step 5: Verify the database 

### Inspect metadata and recorded versions

```bash
sqlite3 -header -column database/inspire.sqlite   "SELECT key, value FROM metadata ORDER BY key;"
```

### Count networks and edges

```bash
sqlite3 -header -column database/inspire.sqlite   "SELECT petal, edge_count, gene_count FROM petal ORDER BY petal;"
```

### Count UniProt mappings

```bash
sqlite3 -header -column database/inspire.sqlite   "SELECT COUNT(*) AS mapping_rows FROM gene_uniprot;"
```

### Count Biological Process annotations

```bash
sqlite3 -header -column database/inspire.sqlite   "SELECT COUNT(*) AS bp_annotations FROM protein_go_annotation;"
```

### Count Pfam annotations, if loaded

```bash
sqlite3 -header -column database/inspire.sqlite   "SELECT COUNT(*) AS pfam_annotations FROM protein_pfam;"
```

### Review unmapped genes

```bash
cat reports/unmapped_genes.csv
```

## Step 6: Preserve the replication package


Do not replace the GO files after a run without rebuilding the annotation portion of the SQLite database. The SQLite `metadata` table records checksums so the exact inputs can be identified later.


##  One consolidated verification command


cd ~/Box/Research/INSPIRE/inspire-replication

echo "=== SQLite integrity ==="
sqlite3 database/inspire.sqlite "PRAGMA integrity_check;"

echo
echo "=== Foreign-key problems, should be empty ==="
sqlite3 database/inspire.sqlite "PRAGMA foreign_key_check;"

echo
echo "=== Core counts ==="
sqlite3 -header -column database/inspire.sqlite "
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
sqlite3 -header -column database/inspire.sqlite "
SELECT pfam_status, COUNT(*) AS proteins
FROM annotation_status
GROUP BY pfam_status
ORDER BY pfam_status;
"

echo
echo "=== Pfam errors ==="
sqlite3 -header -column database/inspire.sqlite "
SELECT uniprot_accession, last_error
FROM annotation_status
WHERE pfam_status='error'
ORDER BY uniprot_accession;
"

echo
echo "=== Genes without exactly one preferred accession ==="
sqlite3 -header -column database/inspire.sqlite "
SELECT
    g.gene_symbol,
    SUM(CASE WHEN u.is_preferred=1 THEN 1 ELSE 0 END) AS preferred_count
FROM gene g
LEFT JOIN gene_uniprot u
    ON u.gene_symbol=g.gene_symbol
GROUP BY g.gene_symbol
HAVING preferred_count != 1;
"

echo
echo "=== Preferred proteins without GO Biological Process annotations ==="
sqlite3 -header -column database/inspire.sqlite "
SELECT
    u.gene_symbol,
    u.uniprot_accession
FROM gene_uniprot u
LEFT JOIN protein_go_annotation a
    ON a.uniprot_accession=u.uniprot_accession
WHERE u.is_preferred=1
GROUP BY u.gene_symbol, u.uniprot_accession
HAVING COUNT(a.go_id)=0
ORDER BY u.gene_symbol;
"

echo
echo "=== Verification finished ==="


## Current endpoint

At the end of these steps, `database/inspire.sqlite` contains:

- petal networks and unique undirected edges
- original edge observations and scores
- GO ontology terms and relationships
- human Biological Process annotations
- gene-to-UniProt mappings derived from the GAF
- optional Pfam annotations retrieved from InterPro


