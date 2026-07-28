#!/usr/bin/env python3
"""
build_go_scores.py

Build offline GO Biological Process similarity scores for INSPIRE.

Supported algorithms
--------------------

1. resnik
   Standard protein-level Resnik similarity:

       score(A, B) = max IC(t), for every propagated GO term t shared by A and B

2. inspire-legacy
   Reproduces the main custom logic of the historical C++ implementation,
   while using floating-point comparison when selecting each lambda term:

       a. Trim each protein's direct GO terms to the most-specific terms.
       b. For every direct-term pair, select the common ancestor with maximum IC.
       c. Deduplicate and trim the selected lambda terms.
       d. Count human proteins annotated to every surviving lambda term.
       e. score = -log2(joint_count / total_annotated_human_proteins)

3. inspire-legacy-exact
   Same as inspire-legacy, but reproduces the old C++ integer truncation in:

       int maxInformationContent = 0;

   The old program could terminate when a lambda set had no jointly annotated
   human genes. This offline implementation records a score of 0.0 instead of
   terminating the entire batch. Such cases are counted and reported.

Important compatibility note
----------------------------

The two legacy modes reproduce the historical algorithm's structure, but they
use the current GO OBO/GAF data and a consistent denominator:

    all distinct human UniProt proteins with loaded Biological Process terms

They therefore will not numerically reproduce a 2011 GO/MySQL database.

Prerequisites
-------------

    python3 build_inspire_database.py ...
    python3 load_annotations.py --db ... --obo ... --gaf ... --skip-pfam

Examples
--------

    python3 build_go_scores.py --db database/inspire.sqlite --algorithm resnik

    python3 build_go_scores.py --db database/inspire.sqlite \
        --algorithm inspire-legacy

    python3 build_go_scores.py --db database/inspire.sqlite \
        --algorithm inspire-legacy-exact

Build all three algorithms:

    for algorithm in resnik inspire-legacy inspire-legacy-exact; do
      python3 build_go_scores.py \
        --db database/inspire.sqlite \
        --algorithm "$algorithm"
    done

Use --rebuild to replace an existing release for the selected algorithm.
The script uses only the Python standard library.
"""

from __future__ import annotations

import argparse
import math
import sqlite3
import statistics
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Sequence, Set, Tuple

SCRIPT_VERSION = "2.0"
BP_ROOT = "GO:0008150"
ALGORITHM_NAMES = {
    "resnik": "resnik_bp_v1",
    "inspire-legacy": "inspire_legacy_v1",
    "inspire-legacy-exact": "inspire_legacy_exact_v1",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def connect_database(path: Path) -> sqlite3.Connection:
    con = sqlite3.connect(path)
    con.execute("PRAGMA foreign_keys=ON")
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA synchronous=NORMAL")
    return con


def require_tables(con: sqlite3.Connection) -> None:
    required = {
        "metadata", "gene", "gene_uniprot", "go_term", "go_parent",
        "protein_go_annotation", "go_term_gene_count", "go_pair_score"
    }
    existing = {
        row[0] for row in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        )
    }
    missing = sorted(required - existing)
    if missing:
        raise RuntimeError(
            "Database is missing table(s): " + ", ".join(missing) +
            ". Run the earlier database and annotation loaders first."
        )


def set_metadata(con: sqlite3.Connection, key: str, value: str) -> None:
    con.execute(
        "INSERT INTO metadata(key,value) VALUES(?,?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        (key, value),
    )


def get_metadata(con: sqlite3.Connection, key: str) -> str:
    row = con.execute("SELECT value FROM metadata WHERE key=?", (key,)).fetchone()
    return str(row[0]) if row else "unknown"


def load_preferred_proteins(con: sqlite3.Connection) -> List[Tuple[str, str]]:
    rows = [
        (str(g), str(p)) for g, p in con.execute(
            "SELECT gene_symbol,uniprot_accession FROM gene_uniprot "
            "WHERE is_preferred=1 ORDER BY gene_symbol"
        )
    ]
    total_genes = con.execute("SELECT COUNT(*) FROM gene").fetchone()[0]
    mapped_genes = len({gene for gene, _ in rows})
    if mapped_genes != total_genes:
        raise RuntimeError(
            f"Preferred mappings cover {mapped_genes} of {total_genes} genes."
        )
    duplicates = con.execute(
        "SELECT gene_symbol FROM gene_uniprot WHERE is_preferred=1 "
        "GROUP BY gene_symbol HAVING COUNT(*) != 1"
    ).fetchall()
    if duplicates:
        raise RuntimeError("A gene has more than one preferred UniProt mapping.")
    return rows


def load_bp_terms(con: sqlite3.Connection) -> Set[str]:
    terms = {
        str(row[0]) for row in con.execute(
            "SELECT go_id FROM go_term "
            "WHERE namespace='biological_process' AND is_obsolete=0"
        )
    }
    if not terms:
        raise RuntimeError("No active Biological Process GO terms were found.")
    if BP_ROOT not in terms:
        raise RuntimeError(f"Biological Process root {BP_ROOT} is missing.")
    return terms


def load_parent_map(
    con: sqlite3.Connection, bp_terms: Set[str]
) -> Dict[str, Tuple[str, ...]]:
    data: Dict[str, Set[str]] = defaultdict(set)
    for child, parent in con.execute(
        "SELECT child_go_id,parent_go_id FROM go_parent"
    ):
        child = str(child)
        parent = str(parent)
        if child in bp_terms and parent in bp_terms:
            data[child].add(parent)
    return {term: tuple(sorted(parents)) for term, parents in data.items()}


def build_ancestor_map(
    bp_terms: Set[str], parent_map: Dict[str, Tuple[str, ...]]
) -> Dict[str, FrozenSet[str]]:
    """Build transitive ancestor closures. Every term includes itself."""
    memo: Dict[str, FrozenSet[str]] = {}

    for start in sorted(bp_terms):
        if start in memo:
            continue
        stack: List[Tuple[str, bool]] = [(start, False)]
        visiting: Set[str] = set()

        while stack:
            term, expanded = stack.pop()
            if term in memo:
                continue
            if expanded:
                result: Set[str] = {term}
                for parent in parent_map.get(term, ()):
                    result.update(memo.get(parent, frozenset({parent})))
                memo[term] = frozenset(result)
                visiting.discard(term)
                continue
            if term in visiting:
                raise RuntimeError(f"Cycle detected in GO graph at {term}.")
            visiting.add(term)
            stack.append((term, True))
            for parent in reversed(parent_map.get(term, ())):
                if parent not in memo:
                    stack.append((parent, False))
    return memo


def load_direct_annotations(
    con: sqlite3.Connection, bp_terms: Set[str]
) -> Dict[str, FrozenSet[str]]:
    direct: Dict[str, Set[str]] = defaultdict(set)
    for accession, go_id in con.execute(
        "SELECT DISTINCT uniprot_accession,go_id "
        "FROM protein_go_annotation WHERE aspect='P'"
    ):
        accession = str(accession)
        go_id = str(go_id)
        if go_id in bp_terms:
            direct[accession].add(go_id)
    if not direct:
        raise RuntimeError("No Biological Process protein annotations found.")
    return {protein: frozenset(terms) for protein, terms in direct.items()}


def propagate_annotations(
    direct: Dict[str, FrozenSet[str]],
    ancestors: Dict[str, FrozenSet[str]],
) -> Dict[str, FrozenSet[str]]:
    propagated: Dict[str, FrozenSet[str]] = {}
    for protein, terms in direct.items():
        expanded: Set[str] = set()
        for term in terms:
            expanded.update(ancestors.get(term, frozenset({term})))
        propagated[protein] = frozenset(expanded)
    return propagated


def trim_to_most_specific(
    terms: Sequence[str], ancestors: Dict[str, FrozenSet[str]]
) -> Tuple[str, ...]:
    """
    Remove a term if it is a proper ancestor of another term in the same set.

    This mirrors Petal::trimGo: retain the relatively specific terms.
    """
    unique = sorted(set(terms))
    retained: List[str] = []
    for candidate in unique:
        redundant = False
        for other in unique:
            if candidate == other:
                continue
            if candidate in ancestors.get(other, frozenset({other})):
                redundant = True
                break
        if not redundant:
            retained.append(candidate)
    return tuple(retained)


def calculate_term_counts(
    propagated: Dict[str, FrozenSet[str]]
) -> Dict[str, int]:
    counts: Dict[str, int] = defaultdict(int)
    for terms in propagated.values():
        for term in terms:
            counts[term] += 1
    return dict(counts)


def calculate_ic(
    counts: Dict[str, int], total_proteins: int
) -> Dict[str, float]:
    if total_proteins <= 0:
        raise RuntimeError("Information-content denominator is zero.")
    ic: Dict[str, float] = {}
    for term, count in counts.items():
        if count <= 0 or count > total_proteins:
            raise RuntimeError(
                f"Invalid GO count for {term}: {count}/{total_proteins}"
            )
        ic[term] = -math.log2(count / total_proteins)
    return ic


def save_term_counts(con: sqlite3.Connection, counts: Dict[str, int]) -> None:
    with con:
        con.execute("DELETE FROM go_term_gene_count")
        con.executemany(
            "INSERT INTO go_term_gene_count(go_id,gene_count) VALUES(?,?)",
            sorted(counts.items()),
        )


def resnik_score(
    terms_a: FrozenSet[str],
    terms_b: FrozenSet[str],
    ic: Dict[str, float],
) -> float:
    if not terms_a or not terms_b:
        return 0.0
    smaller, larger = (
        (terms_a, terms_b) if len(terms_a) <= len(terms_b)
        else (terms_b, terms_a)
    )
    return max(
        (ic.get(term, 0.0) for term in smaller if term in larger),
        default=0.0,
    )


def select_lambda_corrected(
    term_a: str,
    term_b: str,
    ancestors: Dict[str, FrozenSet[str]],
    ic: Dict[str, float],
) -> Optional[str]:
    common = ancestors.get(term_a, frozenset({term_a})) & ancestors.get(
        term_b, frozenset({term_b})
    )
    if not common:
        return None
    # Deterministic tie breaking follows std::map-like lexical ordering: the
    # first term encountered at the maximum score wins.
    best: Optional[str] = None
    best_value = -1.0
    for term in sorted(common):
        value = ic.get(term, 0.0)
        if value > best_value:
            best = term
            best_value = value
    return best


def select_lambda_exact(
    term_a: str,
    term_b: str,
    ancestors: Dict[str, FrozenSet[str]],
    ic: Dict[str, float],
) -> Optional[str]:
    """
    Reproduce the historical int maxInformationContent behavior.

    C++ assignment from double to int truncates toward zero. A term replaces
    the prior lambda only when its double IC is greater than the current
    integer threshold.
    """
    common = ancestors.get(term_a, frozenset({term_a})) & ancestors.get(
        term_b, frozenset({term_b})
    )
    threshold = 0
    best: Optional[str] = None
    for term in sorted(common):
        value = ic.get(term, 0.0)
        if value > threshold:
            threshold = int(value)
            best = term
    return best


def legacy_score(
    protein_a: str,
    protein_b: str,
    trimmed_direct: Dict[str, Tuple[str, ...]],
    ancestors: Dict[str, FrozenSet[str]],
    ic: Dict[str, float],
    propagated: Dict[str, FrozenSet[str]],
    total_proteins: int,
    exact_integer_selection: bool,
) -> Tuple[float, bool, int]:
    """
    Return (score, zero_joint_case, lambda_count).

    The expensive final intersection is performed against the complete human
    protein population loaded from the GAF, not only the 193 network proteins.
    """
    terms_a = trimmed_direct.get(protein_a, ())
    terms_b = trimmed_direct.get(protein_b, ())
    if not terms_a or not terms_b:
        return 0.0, False, 0

    selector = select_lambda_exact if exact_integer_selection else select_lambda_corrected
    lambdas: List[str] = []
    for term_a in terms_a:
        for term_b in terms_b:
            selected = selector(term_a, term_b, ancestors, ic)
            if selected is not None:
                lambdas.append(selected)

    trimmed_lambdas = trim_to_most_specific(lambdas, ancestors)
    if not trimmed_lambdas:
        return 0.0, False, 0

    required = frozenset(trimmed_lambdas)
    joint_count = sum(
        1 for protein_terms in propagated.values()
        if required.issubset(protein_terms)
    )
    if joint_count == 0:
        # Historical C++ called exit(1). A batch builder must remain usable,
        # so record zero and report how often this safeguard was needed.
        return 0.0, True, len(trimmed_lambdas)

    return -math.log2(joint_count / total_proteins), False, len(trimmed_lambdas)


def canonical_pair(a: str, b: str) -> Tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def release_exists(
    con: sqlite3.Connection,
    stored_algorithm: str,
    ontology_version: str,
    annotation_version: str,
) -> bool:
    return con.execute(
        "SELECT 1 FROM go_pair_score WHERE algorithm=? "
        "AND ontology_version=? AND annotation_version=? LIMIT 1",
        (stored_algorithm, ontology_version, annotation_version),
    ).fetchone() is not None


def delete_release(
    con: sqlite3.Connection,
    stored_algorithm: str,
    ontology_version: str,
    annotation_version: str,
) -> None:
    with con:
        con.execute(
            "DELETE FROM go_pair_score WHERE algorithm=? "
            "AND ontology_version=? AND annotation_version=?",
            (stored_algorithm, ontology_version, annotation_version),
        )


def calculate_scores(
    con: sqlite3.Connection,
    requested_algorithm: str,
    stored_algorithm: str,
    preferred_rows: List[Tuple[str, str]],
    direct: Dict[str, FrozenSet[str]],
    propagated: Dict[str, FrozenSet[str]],
    ancestors: Dict[str, FrozenSet[str]],
    ic: Dict[str, float],
    total_proteins: int,
    ontology_version: str,
    annotation_version: str,
    batch_size: int,
) -> Tuple[List[float], int, int]:
    proteins = sorted({accession for _, accession in preferred_rows})
    pair_total = len(proteins) * (len(proteins) + 1) // 2
    trimmed_direct = {
        protein: trim_to_most_specific(terms, ancestors)
        for protein, terms in direct.items()
    }

    print(f"Preferred network proteins: {len(proteins):,}")
    print(f"Pairs including self-pairs: {pair_total:,}")

    rows: List[Tuple[str, str, float, str, str, str]] = []
    scores: List[float] = []
    zero_scores = 0
    zero_joint_cases = 0
    completed = 0
    start = time.time()

    for i, protein_a in enumerate(proteins):
        for protein_b in proteins[i:]:
            if requested_algorithm == "resnik":
                score = resnik_score(
                    propagated.get(protein_a, frozenset()),
                    propagated.get(protein_b, frozenset()),
                    ic,
                )
                zero_joint = False
            else:
                score, zero_joint, _lambda_count = legacy_score(
                    protein_a=protein_a,
                    protein_b=protein_b,
                    trimmed_direct=trimmed_direct,
                    ancestors=ancestors,
                    ic=ic,
                    propagated=propagated,
                    total_proteins=total_proteins,
                    exact_integer_selection=(
                        requested_algorithm == "inspire-legacy-exact"
                    ),
                )

            if score == 0.0:
                zero_scores += 1
            if zero_joint:
                zero_joint_cases += 1

            first, second = canonical_pair(protein_a, protein_b)
            rows.append((
                first, second, score, stored_algorithm,
                ontology_version, annotation_version
            ))
            scores.append(score)
            completed += 1

            if len(rows) >= batch_size:
                with con:
                    con.executemany(
                        "INSERT OR REPLACE INTO go_pair_score"
                        "(protein1,protein2,score,algorithm,ontology_version,"
                        "annotation_version) VALUES(?,?,?,?,?,?)",
                        rows,
                    )
                rows.clear()

            if completed % 1000 == 0 or completed == pair_total:
                print(
                    f"[{completed:,}/{pair_total:,}] "
                    f"{time.time() - start:.1f} seconds elapsed"
                )

    if rows:
        with con:
            con.executemany(
                "INSERT OR REPLACE INTO go_pair_score"
                "(protein1,protein2,score,algorithm,ontology_version,"
                "annotation_version) VALUES(?,?,?,?,?,?)",
                rows,
            )

    return scores, zero_scores, zero_joint_cases


def validate(
    con: sqlite3.Connection,
    stored_algorithm: str,
    ontology_version: str,
    annotation_version: str,
    expected_pairs: int,
) -> None:
    actual = con.execute(
        "SELECT COUNT(*) FROM go_pair_score WHERE algorithm=? "
        "AND ontology_version=? AND annotation_version=?",
        (stored_algorithm, ontology_version, annotation_version),
    ).fetchone()[0]
    if actual != expected_pairs:
        raise RuntimeError(f"Expected {expected_pairs} scores; found {actual}.")

    bad = con.execute(
        "SELECT COUNT(*) FROM go_pair_score WHERE algorithm=? "
        "AND ontology_version=? AND annotation_version=? "
        "AND (score IS NULL OR score < 0)",
        (stored_algorithm, ontology_version, annotation_version),
    ).fetchone()[0]
    if bad:
        raise RuntimeError(f"Found {bad} invalid score(s).")

    integrity = con.execute("PRAGMA integrity_check").fetchone()[0]
    if integrity != "ok":
        raise RuntimeError(f"SQLite integrity check failed: {integrity}")
    foreign_keys = con.execute("PRAGMA foreign_key_check").fetchall()
    if foreign_keys:
        raise RuntimeError(
            f"SQLite foreign-key check found {len(foreign_keys)} errors."
        )


def print_summary(
    stored_algorithm: str,
    scores: List[float],
    zero_scores: int,
    zero_joint_cases: int,
) -> None:
    print("\nGO score summary")
    print("----------------")
    print(f"Stored algorithm:       {stored_algorithm}")
    print(f"Scores:                 {len(scores):,}")
    print(f"Zero scores:            {zero_scores:,}")
    print(f"Zero-joint safeguards:  {zero_joint_cases:,}")
    if scores:
        print(f"Minimum:                {min(scores):.6f}")
        print(f"Median:                 {statistics.median(scores):.6f}")
        print(f"Mean:                   {statistics.fmean(scores):.6f}")
        print(f"Maximum:                {max(scores):.6f}")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build GO pair scores for INSPIRE."
    )
    parser.add_argument("--db", required=True, type=Path)
    parser.add_argument(
        "--algorithm",
        required=True,
        choices=sorted(ALGORITHM_NAMES),
        help="GO similarity algorithm to calculate",
    )
    parser.add_argument(
        "--rebuild",
        action="store_true",
        help="Replace an existing release for the selected algorithm",
    )
    parser.add_argument(
        "--batch-size", type=int, default=1000,
        help="SQLite insert batch size, default 1000",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    db_path = args.db.expanduser().resolve()
    if not db_path.is_file():
        print(f"ERROR: database not found: {db_path}", file=sys.stderr)
        return 2

    stored_algorithm = ALGORITHM_NAMES[args.algorithm]
    con = connect_database(db_path)
    try:
        require_tables(con)
        ontology_version = get_metadata(con, "go_ontology_version")
        annotation_version = get_metadata(con, "go_annotation_version")

        print(f"Database:              {db_path}")
        print(f"Requested algorithm:   {args.algorithm}")
        print(f"Stored algorithm name: {stored_algorithm}")
        print(f"GO ontology version:   {ontology_version}")
        print(f"GO annotation version: {annotation_version}")

        if release_exists(
            con, stored_algorithm, ontology_version, annotation_version
        ):
            if not args.rebuild:
                print(
                    "ERROR: this algorithm and GO release already exists. "
                    "Use --rebuild to replace it.",
                    file=sys.stderr,
                )
                return 2
            print("Deleting existing score release...")
            delete_release(
                con, stored_algorithm, ontology_version, annotation_version
            )

        preferred_rows = load_preferred_proteins(con)
        bp_terms = load_bp_terms(con)
        print(f"Active BP GO terms: {len(bp_terms):,}")

        print("Loading GO parents and building ancestor closures...")
        parent_map = load_parent_map(con, bp_terms)
        ancestors = build_ancestor_map(bp_terms, parent_map)

        print("Loading and propagating human BP annotations...")
        direct = load_direct_annotations(con, bp_terms)
        propagated = propagate_annotations(direct, ancestors)
        total_proteins = len(propagated)
        print(f"Human proteins with BP annotations: {total_proteins:,}")

        print("Calculating term frequencies and information content...")
        counts = calculate_term_counts(propagated)
        ic = calculate_ic(counts, total_proteins)
        save_term_counts(con, counts)

        root_count = counts.get(BP_ROOT, 0)
        if root_count != total_proteins:
            print(
                f"WARNING: BP root count is {root_count:,}, but the annotation "
                f"population is {total_proteins:,}.",
                file=sys.stderr,
            )

        print("Calculating pair scores...")
        scores, zero_scores, zero_joint_cases = calculate_scores(
            con=con,
            requested_algorithm=args.algorithm,
            stored_algorithm=stored_algorithm,
            preferred_rows=preferred_rows,
            direct=direct,
            propagated=propagated,
            ancestors=ancestors,
            ic=ic,
            total_proteins=total_proteins,
            ontology_version=ontology_version,
            annotation_version=annotation_version,
            batch_size=max(1, args.batch_size),
        )

        unique_proteins = len({protein for _, protein in preferred_rows})
        expected_pairs = unique_proteins * (unique_proteins + 1) // 2
        validate(
            con, stored_algorithm, ontology_version,
            annotation_version, expected_pairs
        )

        with con:
            prefix = f"go_score_{stored_algorithm}"
            set_metadata(con, f"{prefix}_built_at", utc_now())
            set_metadata(con, f"{prefix}_script_version", SCRIPT_VERSION)
            set_metadata(
                con, f"{prefix}_population",
                "all distinct UniProt proteins with loaded BP annotations"
            )
            set_metadata(con, f"{prefix}_total_proteins", str(total_proteins))
            set_metadata(con, f"{prefix}_zero_scores", str(zero_scores))
            set_metadata(
                con, f"{prefix}_zero_joint_safeguards",
                str(zero_joint_cases)
            )
            set_metadata(
                con, "go_pair_score_selected_algorithm", stored_algorithm
            )

        print_summary(
            stored_algorithm, scores, zero_scores, zero_joint_cases
        )
        print("\nSQLite integrity check: OK")
        return 0

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    finally:
        con.close()


if __name__ == "__main__":
    raise SystemExit(main())
