#!/usr/bin/env bash
set -euo pipefail

# Build and run INSPIRE from the repository root.
#
# Usage:
#   ./scripts/run_inspire.sh PETAL_A PETAL_B [GO_ALGORITHM]
#
# Example:
#   ./scripts/run_inspire.sh FN1 GLI3 inspire_legacy_v1

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 PETAL_A PETAL_B [GO_ALGORITHM]" >&2
    exit 2
fi

PETAL_A="$1"
PETAL_B="$2"
GO_ALGORITHM="${3:-inspire_legacy_v1}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
SRC_DIR="$PROJECT_ROOT/src"
DB_FILE="$PROJECT_ROOT/database/inspire.sqlite"
REPORT_DIR="$PROJECT_ROOT/reports"
EXECUTABLE="$SRC_DIR/INSPIRE_debug"
LOG_FILE="$REPORT_DIR/${PETAL_A}_vs_${PETAL_B}_${GO_ALGORITHM}.log"

if ! command -v brew >/dev/null 2>&1; then
    echo "ERROR: Homebrew was not found." >&2
    exit 1
fi

if ! command -v sqlite3 >/dev/null 2>&1; then
    echo "ERROR: sqlite3 command was not found." >&2
    exit 1
fi

if [ ! -f "$DB_FILE" ]; then
    echo "ERROR: Database not found: $DB_FILE" >&2
    exit 1
fi

for source_file in gene.cpp petal.cpp AlignedGeneEdge.cpp main.cpp; do
    if [ ! -f "$SRC_DIR/$source_file" ]; then
        echo "ERROR: Missing source file: $SRC_DIR/$source_file" >&2
        exit 1
    fi
done

mkdir -p "$REPORT_DIR"

SQLITE_PREFIX="$(brew --prefix sqlite)"

printf '\n=== 1. Checking SQLite database ===\n'
INTEGRITY="$(sqlite3 "$DB_FILE" 'PRAGMA integrity_check;')"
if [ "$INTEGRITY" != "ok" ]; then
    echo "ERROR: SQLite integrity check failed: $INTEGRITY" >&2
    exit 1
fi
echo "SQLite integrity check: ok"

FK_ERRORS="$(sqlite3 "$DB_FILE" 'PRAGMA foreign_key_check;')"
if [ -n "$FK_ERRORS" ]; then
    echo "ERROR: SQLite foreign-key errors were found:" >&2
    echo "$FK_ERRORS" >&2
    exit 1
fi
echo "SQLite foreign-key check: ok"

printf '\n=== 2. Checking requested petals ===\n'
for petal in "$PETAL_A" "$PETAL_B"; do
    FOUND="$(sqlite3 "$DB_FILE" "SELECT COUNT(*) FROM petal WHERE petal = '$petal';")"
    if [ "$FOUND" != "1" ]; then
        echo "ERROR: Petal not found: $petal" >&2
        exit 1
    fi
    sqlite3 -header -column "$DB_FILE" \
        "SELECT petal, gene_count, edge_count FROM petal WHERE petal = '$petal';"
done

printf '\n=== 3. Checking GO score cache ===\n'
PREFERRED_PROTEINS="$(sqlite3 "$DB_FILE" \
    'SELECT COUNT(DISTINCT uniprot_accession) FROM gene_uniprot WHERE is_preferred = 1;')"
EXPECTED_PAIRS=$((PREFERRED_PROTEINS * (PREFERRED_PROTEINS + 1) / 2))
ACTUAL_PAIRS="$(sqlite3 "$DB_FILE" \
    "SELECT COUNT(*) FROM go_pair_score WHERE algorithm = '$GO_ALGORITHM';")"

printf 'Preferred proteins: %s\n' "$PREFERRED_PROTEINS"
printf 'Expected GO pairs:  %s\n' "$EXPECTED_PAIRS"
printf 'Actual GO pairs:    %s\n' "$ACTUAL_PAIRS"

if [ "$ACTUAL_PAIRS" -ne "$EXPECTED_PAIRS" ]; then
    echo "ERROR: GO score cache is incomplete for $GO_ALGORITHM" >&2
    exit 1
fi

sqlite3 -header -column "$DB_FILE" \
    "SELECT algorithm, COUNT(*) AS pairs, ROUND(MIN(score),6) AS minimum, ROUND(AVG(score),6) AS mean, ROUND(MAX(score),6) AS maximum FROM go_pair_score WHERE algorithm = '$GO_ALGORITHM' GROUP BY algorithm;"

printf '\n=== 4. Checking GO and Pfam coverage for these petals ===\n'
sqlite3 -header -column "$DB_FILE" \
    "SELECT pg.petal, COUNT(*) AS genes, SUM(CASE WHEN u.uniprot_accession IS NOT NULL THEN 1 ELSE 0 END) AS genes_with_preferred_uniprot, SUM(CASE WHEN EXISTS (SELECT 1 FROM protein_go_annotation a WHERE a.uniprot_accession = u.uniprot_accession) THEN 1 ELSE 0 END) AS genes_with_bp_go, SUM(CASE WHEN EXISTS (SELECT 1 FROM protein_pfam f WHERE f.uniprot_accession = u.uniprot_accession) THEN 1 ELSE 0 END) AS genes_with_pfam FROM petal_gene pg LEFT JOIN gene_uniprot u ON u.gene_symbol = pg.gene_symbol AND u.is_preferred = 1 WHERE pg.petal IN ('$PETAL_A', '$PETAL_B') GROUP BY pg.petal ORDER BY pg.petal;"

printf '\n=== 5. Building debug executable ===\n'
clang++ \
    -std=c++11 \
    -O0 \
    -g \
    -Wall \
    -Wextra \
    -Wpedantic \
    -fsanitize=address,undefined \
    -fno-omit-frame-pointer \
    -I"$SQLITE_PREFIX/include" \
    -L"$SQLITE_PREFIX/lib" \
    "$SRC_DIR/gene.cpp" \
    "$SRC_DIR/petal.cpp" \
    "$SRC_DIR/AlignedGeneEdge.cpp" \
    "$SRC_DIR/main.cpp" \
    -lsqlite3 \
    -o "$EXECUTABLE"

echo "Build succeeded: $EXECUTABLE"

printf '\n=== 6. Running INSPIRE ===\n'
"$EXECUTABLE" "$DB_FILE" "$PETAL_A" "$PETAL_B" "$GO_ALGORITHM" 2>&1 | tee "$LOG_FILE"

printf '\n=== Complete ===\n'
echo "Run log: $LOG_FILE"
