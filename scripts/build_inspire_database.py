#!/usr/bin/env python3
"""
build_inspire_database.py

Build a portable SQLite database from a CSV file containing petal networks.

Expected CSV columns:
    gene1,gene2,path,score,petal

Design choices:
    * Each value in the CSV 'petal' column identifies one network.
    * Edges are treated as UNDIRECTED.
    * Reversed duplicates, such as A-B and B-A, are treated as the same edge.
    * Repeated appearances of the same edge within a petal are collapsed.
    * The largest observed PPI/path score is kept as max_score.
    * The original path observations are preserved in edge_observation.
    * All imported edges receive interaction_type='pp'. This is only a
      compatibility label for the old C++ program; the numeric CSV score is
      stored separately and is not discarded.
    * Empty, malformed, and self-loop rows are recorded in import_problem.
    * Annotation tables for UniProt, Pfam, GO, and cached GO pair scores are
      created empty. They will be filled in a later annotation step.

Example:
    python3 build_inspire_database.py \
        blossom_edge_list_7_22.csv \
        inspire.sqlite

Useful checks after import:
    sqlite3 inspire.sqlite \
        "SELECT petal, edge_count FROM petal ORDER BY petal;"

    sqlite3 inspire.sqlite \
        "SELECT COUNT(*) AS unique_genes FROM gene;"
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Tuple

SCHEMA_VERSION = "1"
REQUIRED_COLUMNS = {"gene1", "gene2", "path", "score", "petal"}


def utc_now() -> str:
    """Return an ISO-8601 UTC timestamp."""
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    """Calculate a checksum so the exact input file can be identified later."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_text(value: Optional[str]) -> str:
    """Remove surrounding whitespace and normalize missing CSV values."""
    return "" if value is None else value.strip()


def normalize_gene(value: Optional[str]) -> str:
    """
    Normalize a gene symbol for this project.

    Current policy:
        * trim surrounding whitespace
        * convert to uppercase

    Human gene symbols are conventionally uppercase. If you later need to
    retain case-sensitive identifiers, change this function before import.
    """
    return clean_text(value).upper()


def canonical_edge(gene1: str, gene2: str) -> Tuple[str, str]:
    """
    Return an undirected edge in stable alphabetical order.

    This ensures that A-B and B-A have the same database key.
    """
    return (gene1, gene2) if gene1 <= gene2 else (gene2, gene1)


def parse_score(value: Optional[str]) -> Optional[float]:
    """Parse a finite numeric score; return None for blank or invalid values."""
    text = clean_text(value)
    if not text:
        return None
    try:
        score = float(text)
    except ValueError:
        return None
    return score if math.isfinite(score) else None


def connect_database(path: Path) -> sqlite3.Connection:
    """Open SQLite and turn on the settings used by this database."""
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA foreign_keys = ON")
    connection.execute("PRAGMA journal_mode = WAL")
    connection.execute("PRAGMA synchronous = NORMAL")
    return connection


def create_schema(connection: sqlite3.Connection) -> None:
    """Create the network and annotation tables."""
    connection.executescript(
        """
        CREATE TABLE IF NOT EXISTS metadata (
            key   TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );

        -- One row per named petal/network.
        CREATE TABLE IF NOT EXISTS petal (
            petal       TEXT PRIMARY KEY,
            edge_count  INTEGER NOT NULL DEFAULT 0,
            gene_count  INTEGER NOT NULL DEFAULT 0
        );

        -- One row per unique human gene symbol found in any petal.
        CREATE TABLE IF NOT EXISTS gene (
            gene_symbol TEXT PRIMARY KEY
        );

        -- Membership of genes in petals. This avoids repeatedly deriving the
        -- membership from both sides of the edge table.
        CREATE TABLE IF NOT EXISTS petal_gene (
            petal       TEXT NOT NULL,
            gene_symbol TEXT NOT NULL,
            PRIMARY KEY (petal, gene_symbol),
            FOREIGN KEY (petal) REFERENCES petal(petal) ON DELETE CASCADE,
            FOREIGN KEY (gene_symbol) REFERENCES gene(gene_symbol)
        );

        -- One row per unique undirected edge in a petal.
        -- interaction_type='pp' is retained for compatibility with the old
        -- program. max_score is the score from the new CSV data.
        CREATE TABLE IF NOT EXISTS edge (
            edge_id          INTEGER PRIMARY KEY,
            petal            TEXT NOT NULL,
            gene1            TEXT NOT NULL,
            gene2            TEXT NOT NULL,
            interaction_type TEXT NOT NULL DEFAULT 'pp'
                             CHECK (interaction_type IN ('pp', 'pd')),
            max_score        REAL,
            observation_count INTEGER NOT NULL DEFAULT 1,
            UNIQUE (petal, gene1, gene2),
            CHECK (gene1 < gene2),
            FOREIGN KEY (petal) REFERENCES petal(petal) ON DELETE CASCADE,
            FOREIGN KEY (gene1) REFERENCES gene(gene_symbol),
            FOREIGN KEY (gene2) REFERENCES gene(gene_symbol)
        );

        CREATE INDEX IF NOT EXISTS ix_edge_petal ON edge(petal);
        CREATE INDEX IF NOT EXISTS ix_edge_gene1 ON edge(gene1);
        CREATE INDEX IF NOT EXISTS ix_edge_gene2 ON edge(gene2);

        -- Preserve every valid row from the source CSV, including repeated
        -- occurrences of an edge in different paths.
        CREATE TABLE IF NOT EXISTS edge_observation (
            observation_id INTEGER PRIMARY KEY,
            edge_id        INTEGER NOT NULL,
            source_row     INTEGER NOT NULL,
            original_gene1 TEXT NOT NULL,
            original_gene2 TEXT NOT NULL,
            path           TEXT,
            score          REAL,
            FOREIGN KEY (edge_id) REFERENCES edge(edge_id) ON DELETE CASCADE
        );

        CREATE INDEX IF NOT EXISTS ix_observation_edge
            ON edge_observation(edge_id);

        -- Rows skipped during import are kept here instead of disappearing.
        CREATE TABLE IF NOT EXISTS import_problem (
            source_row INTEGER NOT NULL,
            reason     TEXT NOT NULL,
            gene1      TEXT,
            gene2      TEXT,
            petal      TEXT,
            path       TEXT,
            score_text TEXT
        );

        -- A gene can map to one or more UniProt accessions. The annotation
        -- loader will later mark the preferred accession with is_preferred=1.
        CREATE TABLE IF NOT EXISTS gene_uniprot (
            gene_symbol       TEXT NOT NULL,
            uniprot_accession TEXT NOT NULL,
            reviewed          INTEGER,
            is_preferred      INTEGER NOT NULL DEFAULT 0,
            mapping_source    TEXT NOT NULL,
            source_release    TEXT,
            PRIMARY KEY (gene_symbol, uniprot_accession),
            FOREIGN KEY (gene_symbol) REFERENCES gene(gene_symbol)
        );

        CREATE INDEX IF NOT EXISTS ix_gene_uniprot_preferred
            ON gene_uniprot(gene_symbol, is_preferred);

        -- Pfam member-database accessions obtained from InterPro or a bulk
        -- InterPro/Pfam file. One protein may have many Pfam domains.
        CREATE TABLE IF NOT EXISTS protein_pfam (
            uniprot_accession TEXT NOT NULL,
            pfam_accession    TEXT NOT NULL,
            source            TEXT NOT NULL,
            source_release    TEXT,
            retrieved_at      TEXT,
            PRIMARY KEY (uniprot_accession, pfam_accession)
        );

        -- GO ontology terms parsed from go-basic.obo.
        CREATE TABLE IF NOT EXISTS go_term (
            go_id       TEXT PRIMARY KEY,
            name        TEXT,
            namespace   TEXT,
            is_obsolete INTEGER NOT NULL DEFAULT 0
        );

        -- Direct GO graph relationships parsed from go-basic.obo.
        CREATE TABLE IF NOT EXISTS go_parent (
            child_go_id  TEXT NOT NULL,
            parent_go_id TEXT NOT NULL,
            relation     TEXT NOT NULL,
            PRIMARY KEY (child_go_id, parent_go_id, relation),
            FOREIGN KEY (child_go_id) REFERENCES go_term(go_id),
            FOREIGN KEY (parent_go_id) REFERENCES go_term(go_id)
        );

        CREATE INDEX IF NOT EXISTS ix_go_parent_child
            ON go_parent(child_go_id);

        -- Human protein annotations parsed from the UniProt-centric human GAF.
        CREATE TABLE IF NOT EXISTS protein_go_annotation (
            uniprot_accession TEXT NOT NULL,
            go_id             TEXT NOT NULL,
            evidence_code     TEXT,
            qualifier         TEXT NOT NULL DEFAULT '',
            aspect            TEXT,
            annotation_date   TEXT,
            assigned_by       TEXT,
            PRIMARY KEY (
                uniprot_accession, go_id, evidence_code, qualifier
            ),
            FOREIGN KEY (go_id) REFERENCES go_term(go_id)
        );

        CREATE INDEX IF NOT EXISTS ix_protein_go_accession
            ON protein_go_annotation(uniprot_accession);
        CREATE INDEX IF NOT EXISTS ix_protein_go_term
            ON protein_go_annotation(go_id);

        -- Precomputed annotation frequency used to calculate GO information
        -- content. This replaces the historical GenesAssociated MySQL table.
        CREATE TABLE IF NOT EXISTS go_term_gene_count (
            go_id      TEXT PRIMARY KEY,
            gene_count INTEGER NOT NULL,
            FOREIGN KEY (go_id) REFERENCES go_term(go_id)
        );

        -- Cached GO similarity for an unordered protein pair. protein1 must
        -- sort before or equal to protein2 so a pair is stored only once.
        CREATE TABLE IF NOT EXISTS go_pair_score (
            protein1          TEXT NOT NULL,
            protein2          TEXT NOT NULL,
            score             REAL NOT NULL,
            algorithm         TEXT NOT NULL,
            ontology_version  TEXT NOT NULL,
            annotation_version TEXT NOT NULL,
            PRIMARY KEY (
                protein1, protein2, algorithm,
                ontology_version, annotation_version
            ),
            CHECK (protein1 <= protein2)
        );

        -- Allows the annotation loader to distinguish 'not attempted' from a
        -- successful lookup that found no annotations.
        CREATE TABLE IF NOT EXISTS annotation_status (
            uniprot_accession TEXT PRIMARY KEY,
            pfam_status       TEXT NOT NULL DEFAULT 'pending',
            go_status         TEXT NOT NULL DEFAULT 'pending',
            last_error        TEXT,
            updated_at        TEXT
        );
        """
    )


def set_metadata(connection: sqlite3.Connection, key: str, value: str) -> None:
    """Insert or replace one metadata value."""
    connection.execute(
        "INSERT INTO metadata(key, value) VALUES(?, ?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        (key, value),
    )


def reset_imported_network_data(connection: sqlite3.Connection) -> None:
    """
    Clear only CSV-derived network data.

    Annotation tables are intentionally retained. This permits rebuilding the
    petals from an updated CSV without downloading annotations again.
    """
    connection.executescript(
        """
        DELETE FROM edge_observation;
        DELETE FROM edge;
        DELETE FROM petal_gene;
        DELETE FROM import_problem;
        DELETE FROM petal;
        DELETE FROM gene;
        """
    )


def validate_header(fieldnames: Optional[list[str]]) -> None:
    """Stop with a useful message if required CSV columns are missing."""
    if not fieldnames:
        raise ValueError("The CSV file has no header row.")
    normalized = {name.strip() for name in fieldnames if name is not None}
    missing = REQUIRED_COLUMNS - normalized
    if missing:
        raise ValueError(
            "CSV is missing required column(s): " + ", ".join(sorted(missing))
        )


def import_csv(connection: sqlite3.Connection, csv_path: Path) -> dict[str, int]:
    """Import, normalize, deduplicate, and audit all CSV rows."""
    counts = {
        "rows_seen": 0,
        "valid_observations": 0,
        "problems": 0,
    }

    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        validate_header(reader.fieldnames)

        # CSV line 1 is the header, so the first data row is source row 2.
        for source_row, row in enumerate(reader, start=2):
            counts["rows_seen"] += 1

            original_gene1 = clean_text(row.get("gene1"))
            original_gene2 = clean_text(row.get("gene2"))
            gene1 = normalize_gene(original_gene1)
            gene2 = normalize_gene(original_gene2)
            petal = normalize_gene(row.get("petal"))
            path = clean_text(row.get("path"))
            score_text = clean_text(row.get("score"))
            score = parse_score(score_text)

            problem = None
            if not petal:
                problem = "missing petal"
            elif not gene1 or not gene2:
                problem = "missing gene"
            elif gene1 == gene2:
                # The original alignment code can represent self-loops, but
                # this import omits them because the supplied petals describe
                # path edges and the edge table requires gene1 < gene2.
                problem = "self-loop omitted"
            elif score_text and score is None:
                problem = "invalid score"

            if problem:
                connection.execute(
                    """
                    INSERT INTO import_problem(
                        source_row, reason, gene1, gene2, petal, path, score_text
                    ) VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    (source_row, problem, original_gene1, original_gene2,
                     petal, path, score_text),
                )
                counts["problems"] += 1
                continue

            edge_gene1, edge_gene2 = canonical_edge(gene1, gene2)

            connection.execute(
                "INSERT OR IGNORE INTO petal(petal) VALUES(?)", (petal,)
            )
            connection.execute(
                "INSERT OR IGNORE INTO gene(gene_symbol) VALUES(?)", (gene1,)
            )
            connection.execute(
                "INSERT OR IGNORE INTO gene(gene_symbol) VALUES(?)", (gene2,)
            )
            connection.execute(
                "INSERT OR IGNORE INTO petal_gene(petal, gene_symbol) VALUES(?, ?)",
                (petal, gene1),
            )
            connection.execute(
                "INSERT OR IGNORE INTO petal_gene(petal, gene_symbol) VALUES(?, ?)",
                (petal, gene2),
            )

            # Insert a new unique edge. If it already exists, increment the
            # observation count and keep the maximum available score.
            connection.execute(
                """
                INSERT INTO edge(
                    petal, gene1, gene2, interaction_type,
                    max_score, observation_count
                ) VALUES (?, ?, ?, 'pp', ?, 1)
                ON CONFLICT(petal, gene1, gene2) DO UPDATE SET
                    observation_count = edge.observation_count + 1,
                    max_score = CASE
                        WHEN edge.max_score IS NULL THEN excluded.max_score
                        WHEN excluded.max_score IS NULL THEN edge.max_score
                        WHEN excluded.max_score > edge.max_score
                            THEN excluded.max_score
                        ELSE edge.max_score
                    END
                """,
                (petal, edge_gene1, edge_gene2, score),
            )

            edge_id = connection.execute(
                """
                SELECT edge_id FROM edge
                WHERE petal=? AND gene1=? AND gene2=?
                """,
                (petal, edge_gene1, edge_gene2),
            ).fetchone()[0]

            connection.execute(
                """
                INSERT INTO edge_observation(
                    edge_id, source_row, original_gene1, original_gene2,
                    path, score
                ) VALUES (?, ?, ?, ?, ?, ?)
                """,
                (edge_id, source_row, original_gene1, original_gene2,
                 path or None, score),
            )
            counts["valid_observations"] += 1

    # Store final counts directly on each petal row.
    connection.execute(
        """
        UPDATE petal
        SET edge_count = (
                SELECT COUNT(*) FROM edge WHERE edge.petal = petal.petal
            ),
            gene_count = (
                SELECT COUNT(*) FROM petal_gene
                WHERE petal_gene.petal = petal.petal
            )
        """
    )

    return counts


def check_database(connection: sqlite3.Connection) -> None:
    """Run SQLite's own consistency checks and a few project-specific checks."""
    integrity = connection.execute("PRAGMA integrity_check").fetchone()[0]
    if integrity != "ok":
        raise RuntimeError(f"SQLite integrity check failed: {integrity}")

    foreign_key_errors = connection.execute("PRAGMA foreign_key_check").fetchall()
    if foreign_key_errors:
        raise RuntimeError(
            f"SQLite foreign-key check found {len(foreign_key_errors)} error(s)."
        )

    bad_edges = connection.execute(
        "SELECT COUNT(*) FROM edge WHERE gene1 >= gene2"
    ).fetchone()[0]
    if bad_edges:
        raise RuntimeError(f"Found {bad_edges} non-canonical edge(s).")


def print_summary(connection: sqlite3.Connection, counts: dict[str, int]) -> None:
    """Print a short, readable import report."""
    petal_count = connection.execute("SELECT COUNT(*) FROM petal").fetchone()[0]
    edge_count = connection.execute("SELECT COUNT(*) FROM edge").fetchone()[0]
    gene_count = connection.execute("SELECT COUNT(*) FROM gene").fetchone()[0]

    print("\nImport complete")
    print("---------------")
    print(f"CSV rows read:          {counts['rows_seen']}")
    print(f"Valid observations:     {counts['valid_observations']}")
    print(f"Skipped/problem rows:   {counts['problems']}")
    print(f"Petals:                 {petal_count}")
    print(f"Unique petal edges:     {edge_count}")
    print(f"Unique genes overall:   {gene_count}")

    print("\nPer-petal counts")
    print("----------------")
    for petal, edges, genes in connection.execute(
        "SELECT petal, edge_count, gene_count FROM petal ORDER BY petal"
    ):
        print(f"{petal:12s} edges={edges:4d} genes={genes:4d}")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the INSPIRE SQLite database from a petal-edge CSV."
    )
    parser.add_argument("csv_file", type=Path, help="Input network CSV file")
    parser.add_argument("sqlite_file", type=Path, help="Output SQLite database")
    parser.add_argument(
        "--replace-network-data",
        action="store_true",
        help=(
            "Replace existing petal/edge data while retaining annotation tables. "
            "Without this option, the script refuses to overwrite an existing DB."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    csv_path = args.csv_file.expanduser().resolve()
    db_path = args.sqlite_file.expanduser().resolve()

    if not csv_path.is_file():
        print(f"ERROR: CSV file not found: {csv_path}", file=sys.stderr)
        return 2

    if db_path.exists() and not args.replace_network_data:
        print(
            f"ERROR: Database already exists: {db_path}\n"
            "Use a new output name or add --replace-network-data.",
            file=sys.stderr,
        )
        return 2

    db_path.parent.mkdir(parents=True, exist_ok=True)
    connection = connect_database(db_path)

    try:
        create_schema(connection)
        with connection:
            if args.replace_network_data:
                reset_imported_network_data(connection)

            set_metadata(connection, "schema_version", SCHEMA_VERSION)
            set_metadata(connection, "created_or_updated_at", utc_now())
            set_metadata(connection, "network_source_file", csv_path.name)
            set_metadata(connection, "network_source_sha256", sha256_file(csv_path))
            set_metadata(connection, "edge_direction", "undirected")
            set_metadata(connection, "interaction_type_policy", "all imported as pp")
            set_metadata(connection, "duplicate_score_policy", "maximum score retained")

            counts = import_csv(connection, csv_path)
            set_metadata(connection, "source_rows", str(counts["rows_seen"]))
            set_metadata(
                connection,
                "valid_edge_observations",
                str(counts["valid_observations"]),
            )

        check_database(connection)
        print_summary(connection, counts)
        print(f"\nDatabase written to: {db_path}")
        return 0

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    finally:
        connection.close()


if __name__ == "__main__":
    raise SystemExit(main())
