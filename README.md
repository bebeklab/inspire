# INSPIRE

INSPIRE is a C++ network-alignment program being modernized from a historical MySQL workflow to a portable SQLite workflow.

## Current workflow

1. Import petal networks from CSV into SQLite.
2. Load GO Biological Process annotations.
3. Retrieve and cache Pfam annotations.
4. Precompute GO protein-pair scores.
5. Run the C++ alignment engine using SQLite.

## GO scoring algorithms

- `resnik_bp_v1`
- `inspire_legacy_v1`
- `inspire_legacy_exact_v1`

The current default is `inspire_legacy_v1`.

## Database

The SQLite database is generated locally and is not stored in Git.

Expected local location:

```text
database/inspire.sqlite


