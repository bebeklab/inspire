#!/usr/bin/env python3
"""Run INSPIRE for every unordered petal pair, then validate and plot alignments.

The INSPIRE executable is expected to accept:
    INSPIRE DATABASE PETAL_A PETAL_B ALGORITHM ALIGNMENT_TSV

Outputs under --output-dir:
    logs/               complete stdout/stderr logs
    alignments/         structured final-alignment TSV files from C++
    alignment_plots/    PNG figures generated after all runs finish
    summary_<algorithm>.csv
    failed_jobs.csv
    failed_plots.csv

Examples from repository root:
    python3 scripts/run_all_vs_all.py --jobs 15 --no-self
    python3 scripts/run_all_vs_all.py --jobs 8 --no-self --overwrite
    python3 scripts/run_all_vs_all.py --no-self --no-plots


python3 scripts/run_all_vs_all.py --database database/inspire_allrules.sqlite  --jobs 15 --output-dir reports/all_vs_all_allrules
python3 scripts/run_all_vs_all.py --database database/inspire_bestrule.sqlite  --jobs 15 --output-dir reports/all_vs_all_bestrule
python3 scripts/run_all_vs_all.py --database database/inspire_rulecount.sqlite --jobs 15 --output-dir reports/all_vs_all_rulecount



"""

from __future__ import annotations

import argparse
import csv
import itertools
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
    r"^\$RESULT\$\s+(\S+)\s+(\S+)\s+norm_score:\s+(\S+)",
    re.MULTILINE,
)


@dataclass(frozen=True)
class Job:
    petal_a: str
    petal_b: str
    executable: str
    database: str
    algorithm: str
    console_log: str
    alignment_tsv: str
    overwrite: bool


def parse_result(text: str) -> Optional[tuple[str, str, str]]:
    matches = RESULT_RE.findall(text)
    return matches[-1] if matches else None


def completed_result(console_log: Path, alignment_tsv: Path) -> Optional[tuple[str, str, str]]:
    """Return a prior result only when both output files are complete."""
    if not console_log.is_file() or not alignment_tsv.is_file():
        return None
    try:
        text = console_log.read_text(encoding="utf-8", errors="replace")
        alignment_text = alignment_tsv.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    result = parse_result(text)
    if result is None:
        return None
    if "$ALIGNMENT_BEGIN$" not in alignment_text or "$ALIGNMENT_END$" not in alignment_text:
        return None
    if "$NETWORK_A_EDGE$" not in alignment_text or "$NETWORK_B_EDGE$" not in alignment_text:
        return None
    return result


def run_one(job: Job) -> dict[str, object]:
    """Run one pair and save complete console and structured alignment outputs."""
    start = time.monotonic()
    console_path = Path(job.console_log)
    alignment_path = Path(job.alignment_tsv)
    console_path.parent.mkdir(parents=True, exist_ok=True)
    alignment_path.parent.mkdir(parents=True, exist_ok=True)

    if not job.overwrite:
        previous = completed_result(console_path, alignment_path)
        if previous is not None:
            return {
                "petal_a": previous[0],
                "petal_b": previous[1],
                "score": previous[2],
                "status": "skipped",
                "returncode": 0,
                "seconds": 0.0,
                "console_log": str(console_path),
                "alignment_tsv": str(alignment_path),
                "error": "",
            }

    temporary_console = console_path.with_suffix(console_path.suffix + ".tmp")
    temporary_alignment = alignment_path.with_suffix(alignment_path.suffix + ".tmp")
    for temporary in (temporary_console, temporary_alignment):
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass

    command = [
        job.executable,
        job.database,
        job.petal_a,
        job.petal_b,
        job.algorithm,
        str(temporary_alignment),
    ]

    try:
        completed = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
            check=False,
        )
        output = completed.stdout or ""
        temporary_console.write_text(output, encoding="utf-8")
        result = parse_result(output)

        alignment_complete = False
        if temporary_alignment.is_file():
            alignment_text = temporary_alignment.read_text(
                encoding="utf-8", errors="replace"
            )
            alignment_complete = (
                "$ALIGNMENT_BEGIN$" in alignment_text
                and "$ALIGNMENT_END$" in alignment_text
                and "$NETWORK_A_EDGE$" in alignment_text
                and "$NETWORK_B_EDGE$" in alignment_text
                and "$ALIGNMENT_EDGE$" in alignment_text
            )

        success = (
            completed.returncode == 0
            and result is not None
            and alignment_complete
        )

        if success:
            temporary_console.replace(console_path)
            temporary_alignment.replace(alignment_path)
            return {
                "petal_a": result[0],
                "petal_b": result[1],
                "score": result[2],
                "status": "completed",
                "returncode": completed.returncode,
                "seconds": time.monotonic() - start,
                "console_log": str(console_path),
                "alignment_tsv": str(alignment_path),
                "error": "",
            }

        # Preserve failed console output for diagnosis, but never publish a partial TSV.
        failed_console = console_path.with_suffix(".failed.log")
        temporary_console.replace(failed_console)
        try:
            temporary_alignment.unlink()
        except FileNotFoundError:
            pass
        reasons = []
        if completed.returncode != 0:
            reasons.append(f"returncode={completed.returncode}")
        if result is None:
            reasons.append("missing $RESULT$")
        if not alignment_complete:
            reasons.append("incomplete structured alignment")
        return {
            "petal_a": job.petal_a,
            "petal_b": job.petal_b,
            "score": "",
            "status": "failed",
            "returncode": completed.returncode,
            "seconds": time.monotonic() - start,
            "console_log": str(failed_console),
            "alignment_tsv": str(alignment_path),
            "error": "; ".join(reasons),
        }
    except Exception as exc:
        return {
            "petal_a": job.petal_a,
            "petal_b": job.petal_b,
            "score": "",
            "status": "failed",
            "returncode": -1,
            "seconds": time.monotonic() - start,
            "console_log": str(console_path),
            "alignment_tsv": str(alignment_path),
            "error": repr(exc),
        }


def load_petals(database: Path) -> list[str]:
    with sqlite3.connect(str(database)) as connection:
        rows = connection.execute(
            "SELECT petal FROM petal ORDER BY petal"
        ).fetchall()
    petals = [str(row[0]) for row in rows]
    if not petals:
        raise RuntimeError("No petals found in SQLite table 'petal'")
    return petals


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def make_plots(
    successful_results: list[dict[str, object]],
    plot_script: Path,
    plot_dir: Path,
    dpi: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Generate plots only after all INSPIRE jobs have finished."""
    plot_dir.mkdir(parents=True, exist_ok=True)
    plotted = []
    failed = []

    for index, result in enumerate(successful_results, start=1):
        alignment_path = Path(str(result["alignment_tsv"]))
        plot_path = plot_dir / f"{alignment_path.stem}.png"
        command = [
            sys.executable,
            str(plot_script),
            str(alignment_path),
            "--output",
            str(plot_path),
            "--dpi",
            str(dpi),
        ]
        completed = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
            check=False,
        )
        if completed.returncode == 0 and plot_path.is_file():
            item = dict(result)
            item["plot"] = str(plot_path)
            plotted.append(item)
            print(
                f"[plot {index}/{len(successful_results)}] "
                f"{result['petal_a']} vs {result['petal_b']} -> {plot_path.name}"
            )
        else:
            failed.append({
                "petal_a": result["petal_a"],
                "petal_b": result["petal_b"],
                "alignment_tsv": str(alignment_path),
                "plot": str(plot_path),
                "error": completed.stdout.strip(),
            })
            print(
                f"[plot FAILED] {result['petal_a']} vs {result['petal_b']}",
                file=sys.stderr,
            )

    return plotted, failed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run all INSPIRE pairs, write logs/TSVs, then create plots."
    )
    parser.add_argument("--jobs", type=int, default=15)
    parser.add_argument(
        "--algorithm",
        default="inspire_legacy_v1",
        choices=[
            "inspire_legacy_v1",
            "inspire_legacy_exact_v1",
            "resnik_bp_v1",
        ],
    )
    parser.add_argument("--no-self", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--plot-overwrite", action="store_true")
    parser.add_argument("--plot-dpi", type=int, default=180)
    parser.add_argument("--network-gap", type=float, default=None)
    parser.add_argument("--executable", type=Path)
    parser.add_argument("--database", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--plot-script", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        print("ERROR: --jobs must be at least 1", file=sys.stderr)
        return 2
    if args.plot_dpi < 1:
        print("ERROR: --plot-dpi must be positive", file=sys.stderr)
        return 2

    project_root = Path(__file__).resolve().parents[1]
    executable = (args.executable or project_root / "src" / "INSPIRE").resolve()
    database = (args.database or project_root / "database" / "inspire.sqlite").resolve()
    output_dir = (args.output_dir or project_root / "reports" / "all_vs_all").resolve()
    plot_script = (
        args.plot_script or project_root / "scripts" / "plot_inspire_alignment.py"
    ).resolve()

    if not executable.is_file():
        print(f"ERROR: executable not found: {executable}", file=sys.stderr)
        return 2
    if not database.is_file():
        print(f"ERROR: database not found: {database}", file=sys.stderr)
        return 2
    if not args.no_plots and not plot_script.is_file():
        print(f"ERROR: plot script not found: {plot_script}", file=sys.stderr)
        return 2

    log_dir = output_dir / "logs"
    alignment_dir = output_dir / "alignments"
    plot_dir = output_dir / "alignment_plots"
    for directory in (log_dir, alignment_dir):
        directory.mkdir(parents=True, exist_ok=True)

    try:
        petals = load_petals(database)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    pair_iterator = (
        itertools.combinations(petals, 2)
        if args.no_self
        else itertools.combinations_with_replacement(petals, 2)
    )
    pairs = list(pair_iterator)
    jobs = []
    for petal_a, petal_b in pairs:
        stem = f"{petal_a}__{petal_b}__{args.algorithm}"
        jobs.append(Job(
            petal_a=petal_a,
            petal_b=petal_b,
            executable=str(executable),
            database=str(database),
            algorithm=args.algorithm,
            console_log=str(log_dir / f"{stem}.log"),
            alignment_tsv=str(alignment_dir / f"{stem}.tsv"),
            overwrite=args.overwrite,
        ))

    print(f"Petals: {len(petals)}")
    print(f"Pairs: {len(jobs)}")
    print(f"Workers: {args.jobs}")
    print(f"Output: {output_dir}")

    results = []
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        future_to_job = {executor.submit(run_one, job): job for job in jobs}
        for index, future in enumerate(as_completed(future_to_job), start=1):
            job = future_to_job[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {
                    "petal_a": job.petal_a,
                    "petal_b": job.petal_b,
                    "score": "",
                    "status": "failed",
                    "returncode": -1,
                    "seconds": 0.0,
                    "console_log": job.console_log,
                    "alignment_tsv": job.alignment_tsv,
                    "error": repr(exc),
                }
            results.append(result)
            print(
                f"[{index}/{len(jobs)}] {result['petal_a']} vs "
                f"{result['petal_b']}: {result['status']} "
                f"score={result['score']}"
            )

    results.sort(key=lambda row: (str(row["petal_a"]), str(row["petal_b"])))
    fields = [
        "petal_a", "petal_b", "score", "status", "returncode",
        "seconds", "console_log", "alignment_tsv", "error",
    ]
    summary_path = output_dir / f"summary_{args.algorithm}.csv"
    write_csv(summary_path, results, fields)

    failed_runs = [row for row in results if row["status"] == "failed"]
    write_csv(output_dir / "failed_jobs.csv", failed_runs, fields)
    successful = [row for row in results if row["status"] != "failed"]

    print(
        f"INSPIRE finished in {time.monotonic() - started:.1f}s: "
        f"successful={len(successful)}, failed={len(failed_runs)}"
    )

    failed_plots = []
    if not args.no_plots and successful:
        # Preserve existing plots unless the user explicitly requests replacement.
        to_plot = []
        for result in successful:
            alignment_path = Path(str(result["alignment_tsv"]))
            expected_plot = plot_dir / f"{alignment_path.stem}.png"
            if args.plot_overwrite or not expected_plot.is_file():
                to_plot.append(result)
        print(f"Creating {len(to_plot)} plots after all runs completed...")
        _, failed_plots = make_plots(
            to_plot,
            plot_script,
            plot_dir,
            args.plot_dpi,
        )
        plot_fields = ["petal_a", "petal_b", "alignment_tsv", "plot", "error"]
        write_csv(output_dir / "failed_plots.csv", failed_plots, plot_fields)
        print(
            f"Plotting finished: successful={len(to_plot) - len(failed_plots)}, "
            f"failed={len(failed_plots)}"
        )

    print(f"Summary: {summary_path}")
    return 1 if failed_runs or failed_plots else 0


if __name__ == "__main__":
    raise SystemExit(main())
