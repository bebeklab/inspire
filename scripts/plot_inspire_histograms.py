#!/usr/bin/env python3
"""Create score and runtime histograms from INSPIRE all-vs-all log files.

By default, self-alignments are excluded because they are expected to score 1
and can take disproportionately long in the current implementation.

Example from repository root:
    python3 scripts/plot_inspire_histograms.py \
        --log-dir reports/all_vs_all/logs \
        --output reports/all_vs_all/histograms_inspire_legacy_v1.png

Include self-alignments if desired:
    python3 scripts/plot_inspire_histograms.py \
        --log-dir reports/all_vs_all/logs \
        --output reports/all_vs_all/histograms_with_self.png \
        --include-self
"""

from __future__ import annotations

import argparse
import csv
import re
import statistics
import sys
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt

RESULT_RE = re.compile(
    r"\$RESULT\$\s+([^\s]+)\s+([^\s]+)\s+norm_score:\s+([^\s]+)"
)
TIME_RE = re.compile(r"\[(\d+)\s+min\s+(\d+)\s+sec\]\s+All done in")
COMMAND_RE = re.compile(
    r"COMMAND:\s+.*?\s([^\s]+)\s+([^\s]+)\s+([^\s]+)\s*$",
    re.MULTILINE,
)


def parse_float(value: str) -> Optional[float]:
    try:
        return float(value)
    except ValueError:
        return None


def parse_log(path: Path) -> dict[str, object]:
    text = path.read_text(encoding="utf-8", errors="replace")

    result_matches = RESULT_RE.findall(text)
    time_matches = TIME_RE.findall(text)
    command_match = COMMAND_RE.search(text)

    petal_a = ""
    petal_b = ""
    algorithm = ""
    score: Optional[float] = None

    if result_matches:
        petal_a, petal_b, score_text = result_matches[-1]
        score = parse_float(score_text)
    elif command_match:
        petal_a, petal_b, algorithm = command_match.groups()

    runtime_seconds: Optional[int] = None
    if time_matches:
        minutes, seconds = time_matches[-1]
        runtime_seconds = int(minutes) * 60 + int(seconds)

    # Recover the algorithm from the filename when it was not found in COMMAND.
    if not algorithm:
        stem_parts = path.stem.split("__")
        if len(stem_parts) >= 3:
            algorithm = stem_parts[-1]

    return {
        "petal_a": petal_a,
        "petal_b": petal_b,
        "algorithm": algorithm,
        "score": score,
        "runtime_seconds": runtime_seconds,
        "status": "complete" if score is not None else "incomplete",
        "log_file": str(path),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot INSPIRE score and runtime histograms from log files."
    )
    parser.add_argument("--log-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--algorithm", help="Keep only this algorithm")
    parser.add_argument(
        "--include-self", action="store_true",
        help="Include A-versus-A runs in plots",
    )
    parser.add_argument("--bins", type=int, default=20)
    parser.add_argument(
        "--summary-csv", type=Path,
        help="Optional parsed-log CSV output",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_dir = args.log_dir.expanduser().resolve()
    output = args.output.expanduser().resolve()

    if not log_dir.is_dir():
        print(f"ERROR: log directory not found: {log_dir}", file=sys.stderr)
        return 2

    log_files = sorted(log_dir.glob("*.log"))
    if not log_files:
        print(f"ERROR: no .log files found in {log_dir}", file=sys.stderr)
        return 2

    records = [parse_log(path) for path in log_files]

    if args.algorithm:
        records = [r for r in records if r["algorithm"] == args.algorithm]

    complete_records = [r for r in records if r["status"] == "complete"]
    self_records = [r for r in complete_records if r["petal_a"] == r["petal_b"]]

    if not args.include_self:
        complete_records = [
            r for r in complete_records if r["petal_a"] != r["petal_b"]
        ]

    if not complete_records:
        print("ERROR: no completed runs matched the requested filters", file=sys.stderr)
        return 1

    scores = [float(r["score"]) for r in complete_records if r["score"] is not None]
    runtimes = [
        int(r["runtime_seconds"])
        for r in complete_records
        if r["runtime_seconds"] is not None
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].hist(scores, bins=max(1, args.bins))
    axes[0].set_title("INSPIRE normalized scores")
    axes[0].set_xlabel("Normalized score")
    axes[0].set_ylabel("Number of network pairs")

    if runtimes:
        runtime_minutes = [value / 60.0 for value in runtimes]
        axes[1].hist(runtime_minutes, bins=max(1, args.bins))
        axes[1].set_title("INSPIRE runtimes")
        axes[1].set_xlabel("Runtime in minutes")
        axes[1].set_ylabel("Number of network pairs")
    else:
        axes[1].text(
            0.5, 0.5,
            "No completed runtime lines found",
            horizontalalignment="center",
            verticalalignment="center",
            transform=axes[1].transAxes,
        )
        axes[1].set_axis_off()

    subtitle = (
        f"Completed pairs: {len(complete_records):,}; "
        f"self-alignments {'included' if args.include_self else 'excluded'}"
    )
    fig.suptitle(subtitle)
    fig.tight_layout()

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200)
    plt.close(fig)

    summary_csv = (
        args.summary_csv.expanduser().resolve()
        if args.summary_csv
        else output.with_suffix(".csv")
    )
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "petal_a", "petal_b", "algorithm", "score",
                "runtime_seconds", "status", "log_file",
            ],
        )
        writer.writeheader()
        writer.writerows(records)

    print(f"Log files found:          {len(log_files):,}")
    print(f"Completed logs plotted:   {len(complete_records):,}")
    print(f"Completed self-alignments:{len(self_records):,}")
    print(f"Incomplete logs:          {sum(r['status'] == 'incomplete' for r in records):,}")
    print(f"Score minimum:            {min(scores):.6f}")
    print(f"Score median:             {statistics.median(scores):.6f}")
    print(f"Score mean:               {statistics.fmean(scores):.6f}")
    print(f"Score maximum:            {max(scores):.6f}")
    if runtimes:
        print(f"Median runtime seconds:   {statistics.median(runtimes):.1f}")
        print(f"Maximum runtime seconds:  {max(runtimes)}")
    print(f"Plot written to:          {output}")
    print(f"Parsed CSV written to:    {summary_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
