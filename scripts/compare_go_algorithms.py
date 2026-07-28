#!/usr/bin/env python3
"""Compare two GO scoring algorithms stored in inspire.sqlite.

Default comparison:
    inspire_legacy_v1 versus inspire_legacy_exact_v1

Example:
    python3 scripts/compare_go_algorithms.py \
      --db database/inspire.sqlite \
      --output reports/legacy_vs_exact.png

Other comparison:
    python3 scripts/compare_go_algorithms.py \
      --db database/inspire.sqlite \
      --algorithm-x inspire_legacy_v1 \
      --algorithm-y resnik_bp_v1 \
      --output reports/legacy_vs_resnik.png

The plot contains one dot per protein pair. The dashed diagonal is y=x.
Points on the diagonal have identical scores under both algorithms.
"""

from __future__ import annotations

import argparse
import math
import sqlite3
import statistics
import sys
from pathlib import Path

import matplotlib.pyplot as plt

DEFAULT_X = "inspire_legacy_v1"
DEFAULT_Y = "inspire_legacy_exact_v1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Dot plot comparing two GO algorithms")
    parser.add_argument("--db", required=True, type=Path)
    parser.add_argument("--algorithm-x", default=DEFAULT_X)
    parser.add_argument("--algorithm-y", default=DEFAULT_Y)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--include-self-pairs",
        action="store_true",
        help="Include protein self-comparisons; excluded by default",
    )
    parser.add_argument(
        "--label-largest", type=int, default=10,
        help="Label this many pairs with the largest absolute difference",
    )
    return parser.parse_args()


def pearson(x: list[float], y: list[float]) -> float:
    if len(x) < 2:
        return float("nan")
    mx = statistics.fmean(x)
    my = statistics.fmean(y)
    numerator = sum((a - mx) * (b - my) for a, b in zip(x, y))
    dx = math.sqrt(sum((a - mx) ** 2 for a in x))
    dy = math.sqrt(sum((b - my) ** 2 for b in y))
    return numerator / (dx * dy) if dx and dy else float("nan")


def main() -> int:
    args = parse_args()
    db = args.db.expanduser().resolve()
    output = args.output.expanduser().resolve()
    if not db.is_file():
        print(f"ERROR: database not found: {db}", file=sys.stderr)
        return 2

    con = sqlite3.connect(db)
    try:
        self_filter = "" if args.include_self_pairs else "AND x.protein1 <> x.protein2"
        rows = con.execute(
            f"""
            SELECT
                x.protein1,
                x.protein2,
                x.score,
                y.score,
                COALESCE(g1.gene_symbol, x.protein1),
                COALESCE(g2.gene_symbol, x.protein2)
            FROM go_pair_score x
            JOIN go_pair_score y
              ON y.protein1 = x.protein1
             AND y.protein2 = x.protein2
             AND y.ontology_version = x.ontology_version
             AND y.annotation_version = x.annotation_version
            LEFT JOIN gene_uniprot g1
              ON g1.uniprot_accession = x.protein1
             AND g1.is_preferred = 1
            LEFT JOIN gene_uniprot g2
              ON g2.uniprot_accession = x.protein2
             AND g2.is_preferred = 1
            WHERE x.algorithm = ?
              AND y.algorithm = ?
              {self_filter}
            ORDER BY x.protein1, x.protein2
            """,
            (args.algorithm_x, args.algorithm_y),
        ).fetchall()
    finally:
        con.close()

    if not rows:
        print(
            "ERROR: no matched score pairs found. Check algorithm names and build both caches.",
            file=sys.stderr,
        )
        return 1

    x = [float(row[2]) for row in rows]
    y = [float(row[3]) for row in rows]
    differences = [b - a for a, b in zip(x, y)]
    equal_count = sum(1 for d in differences if abs(d) < 1e-12)
    correlation = pearson(x, y)

    low = min(min(x), min(y))
    high = max(max(x), max(y))
    padding = max((high - low) * 0.04, 0.1)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(x, y, s=13, alpha=0.45)
    ax.plot(
        [low - padding, high + padding],
        [low - padding, high + padding],
        linestyle="--",
        linewidth=1,
    )
    ax.set_xlim(low - padding, high + padding)
    ax.set_ylim(low - padding, high + padding)
    ax.set_xlabel(args.algorithm_x)
    ax.set_ylabel(args.algorithm_y)
    ax.set_title("GO pair-score comparison")

    summary = (
        f"pairs = {len(rows):,}\n"
        f"Pearson r = {correlation:.4f}\n"
        f"identical = {equal_count:,} ({100 * equal_count / len(rows):.1f}%)\n"
        f"mean(y-x) = {statistics.fmean(differences):.4f}"
    )
    ax.text(
        0.03, 0.97, summary,
        transform=ax.transAxes,
        verticalalignment="top",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85},
    )

    label_count = max(0, min(args.label_largest, len(rows)))
    largest = sorted(
        range(len(rows)),
        key=lambda i: abs(differences[i]),
        reverse=True,
    )[:label_count]
    for i in largest:
        gene1 = str(rows[i][4])
        gene2 = str(rows[i][5])
        ax.annotate(
            f"{gene1}:{gene2}",
            (x[i], y[i]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
        )

    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200)
    plt.close(fig)

    print(f"Compared pairs: {len(rows):,}")
    print(f"Pearson correlation: {correlation:.6f}")
    print(f"Identical scores: {equal_count:,}")
    print(f"Mean difference (y-x): {statistics.fmean(differences):.6f}")
    print(f"Median absolute difference: {statistics.median(abs(d) for d in differences):.6f}")
    print(f"Plot written to: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
