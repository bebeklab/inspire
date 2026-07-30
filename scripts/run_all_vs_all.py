#!/usr/bin/env python3
"""Run INSPIRE for every unordered petal pair in parallel.

Defaults:
    * 15 parallel worker processes
    * include self-comparisons
    * use inspire_legacy_v1
    * skip completed logs when resuming

Example from repository root:
    python3 scripts/run_all_vs_all.py

Example with options:
    python3 scripts/run_all_vs_all.py \
        --jobs 15 \
        --algorithm inspire_legacy_v1 \
        --no-self
"""

from __future__ import annotations

import argparse
import csv
import itertools
import os
import re
import sqlite3
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

RESULT_RE = re.compile(
    r"\$RESULT\$\s+([^\s]+)\s+([^\s]+)\s+norm_score:\s+([^\s]+)"
)


@dataclass(frozen=True)
class Job:
    petal_a: str
    petal_b: str
    executable: str
    database: str
    algorithm: str
    log_file: str
    overwrite: bool


def completed_result(log_file: Path) -> Optional[tuple[str, str, str]]:
    """Return an existing result if the log contains a complete result line."""
    if not log_file.is_file():
        return None
    try:
        text = log_file.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    matches = RESULT_RE.findall(text)
    return matches[-1] if matches else None


def run_one(job: Job) -> dict[str, object]:
    """Run one pair in a worker process and save stdout/stderr to its own log."""
    log_path = Path(job.log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    if not job.overwrite:
        prior = completed_result(log_path)
        if prior is not None:
            return {
                "petal_a": job.petal_a,
                "petal_b": job.petal_b,
                "status": "skipped",
                "return_code": 0,
                "score": prior[2],
                "seconds": 0.0,
                "log_file": str(log_path),
            }

    command = [
        job.executable,
        job.database,
        job.petal_a,
        job.petal_b,
        job.algorithm,
    ]

    started = time.monotonic()
    with log_path.open("w", encoding="utf-8") as log:
        log.write("COMMAND: " + " ".join(command) + "\n\n")
        log.flush()
        process = subprocess.run(
            command,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    seconds = time.monotonic() - started

    result = completed_result(log_path)
    score = result[2] if result is not None else ""
    status = "complete" if process.returncode == 0 and result is not None else "failed"

    return {
        "petal_a": job.petal_a,
        "petal_b": job.petal_b,
        "status": status,
        "return_code": process.returncode,
        "score": score,
        "seconds": round(seconds, 3),
        "log_file": str(log_path),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run INSPIRE for all unordered petal pairs in parallel."
    )
    parser.add_argument("--jobs", type=int, default=15,
                        help="Parallel worker processes, default 15")
    parser.add_argument("--algorithm", default="inspire_legacy_v1",
                        choices=["inspire_legacy_v1",
                                 "inspire_legacy_exact_v1",
                                 "resnik_bp_v1"])
    parser.add_argument("--no-self", action="store_true",
                        help="Exclude A-versus-A comparisons")
    parser.add_argument("--overwrite", action="store_true",
                        help="Rerun pairs that already have a result log")
    parser.add_argument("--executable", type=Path,
                        help="INSPIRE executable; default PROJECT_ROOT/src/INSPIRE")
    parser.add_argument("--database", type=Path,
                        help="SQLite database; default PROJECT_ROOT/database/inspire.sqlite")
    parser.add_argument("--output-dir", type=Path,
                        help="Output directory; default PROJECT_ROOT/reports/all_vs_all")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        print("ERROR: --jobs must be at least 1", file=sys.stderr)
        return 2

    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    executable = (args.executable or project_root / "src" / "INSPIRE").resolve()
    database = (args.database or project_root / "database" / "inspire.sqlite").resolve()
    output_dir = (args.output_dir or project_root / "reports" / "all_vs_all").resolve()
    log_dir = output_dir / "logs"
    summary_file = output_dir / f"summary_{args.algorithm}.csv"

    if not executable.is_file() or not os.access(executable, os.X_OK):
        print(f"ERROR: executable is missing or not executable: {executable}", file=sys.stderr)
        return 2
    if not database.is_file():
        print(f"ERROR: database not found: {database}", file=sys.stderr)
        return 2

    con = sqlite3.connect(database)
    try:
        integrity = con.execute("PRAGMA integrity_check").fetchone()[0]
        if integrity != "ok":
            print(f"ERROR: SQLite integrity check failed: {integrity}", file=sys.stderr)
            return 1
        petals = [row[0] for row in con.execute(
            "SELECT petal FROM petal ORDER BY petal"
        )]
        score_count = con.execute(
            "SELECT COUNT(*) FROM go_pair_score WHERE algorithm=?",
            (args.algorithm,),
        ).fetchone()[0]
    finally:
        con.close()

    if not petals:
        print("ERROR: no petals found in database", file=sys.stderr)
        return 1
    if score_count == 0:
        print(f"ERROR: no GO scores found for {args.algorithm}", file=sys.stderr)
        return 1

    pair_iterator = (
        itertools.combinations(petals, 2)
        if args.no_self
        else itertools.combinations_with_replacement(petals, 2)
    )
    pairs = list(pair_iterator)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    jobs = [
        Job(
            petal_a=a,
            petal_b=b,
            executable=str(executable),
            database=str(database),
            algorithm=args.algorithm,
            log_file=str(log_dir / f"{a}__{b}__{args.algorithm}.log"),
            overwrite=args.overwrite,
        )
        for a, b in pairs
    ]

    print(f"Petals:              {len(petals)}")
    print(f"Pairwise runs:       {len(jobs)}")
    print(f"Parallel processes:  {args.jobs}")
    print(f"Algorithm:           {args.algorithm}")
    print(f"Include self-pairs:  {not args.no_self}")
    print(f"Output directory:    {output_dir}")

    results: list[dict[str, object]] = []
    failures = 0
    started = time.monotonic()

    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        future_to_job = {executor.submit(run_one, job): job for job in jobs}
        for completed, future in enumerate(as_completed(future_to_job), start=1):
            job = future_to_job[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {
                    "petal_a": job.petal_a,
                    "petal_b": job.petal_b,
                    "status": "failed",
                    "return_code": -1,
                    "score": "",
                    "seconds": 0.0,
                    "log_file": job.log_file,
                }
                print(f"[{completed}/{len(jobs)}] FAILED {job.petal_a} vs {job.petal_b}: {exc}")
            else:
                print(
                    f"[{completed}/{len(jobs)}] {result['status'].upper():8s} "
                    f"{job.petal_a} vs {job.petal_b} "
                    f"score={result['score']} seconds={result['seconds']}"
                )
            if result["status"] == "failed":
                failures += 1
            results.append(result)

    results.sort(key=lambda row: (str(row["petal_a"]), str(row["petal_b"])))
    with summary_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["petal_a", "petal_b", "status", "return_code",
                        "score", "seconds", "log_file"],
        )
        writer.writeheader()
        writer.writerows(results)

    elapsed = time.monotonic() - started
    print("\nAll-vs-all run finished")
    print(f"Elapsed seconds: {elapsed:.1f}")
    print(f"Failures:        {failures}")
    print(f"Summary:         {summary_file}")

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
