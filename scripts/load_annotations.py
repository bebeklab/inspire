#!/usr/bin/env python3
"""Load GO annotations and optional Pfam domains into inspire.sqlite.

Normal GO-only use with manually downloaded files:

    python3 load_annotations.py \
      --db inspire.sqlite \
      --obo go-basic.obo \
      --gaf HUMAN-uniprot.gaf.gz \
      --skip-pfam

GO plus Pfam retrieval from InterPro:

    python3 load_annotations.py \
      --db inspire.sqlite \
      --obo go-basic.obo \
      --gaf HUMAN-uniprot.gaf.gz

Resume only Pfam retrieval:

    python3 load_annotations.py --db inspire.sqlite --pfam-only

The script uses only the Python standard library. It imports Biological Process
annotations only and excludes GAF records carrying the NOT qualifier.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, Tuple

INTERPRO_BASE = "https://www.ebi.ac.uk/interpro/api"
USER_AGENT = "INSPIRE-annotation-loader/1.0"


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def normalize_gene(value: str) -> str:
    return value.strip().upper()


def connect_db(path: Path) -> sqlite3.Connection:
    con = sqlite3.connect(path)
    con.execute("PRAGMA foreign_keys=ON")
    con.execute("PRAGMA journal_mode=WAL")
    return con


def require_schema(con: sqlite3.Connection) -> None:
    required = {
        "metadata", "gene", "gene_uniprot", "protein_pfam", "go_term",
        "go_parent", "protein_go_annotation", "go_term_gene_count",
        "go_pair_score", "annotation_status"
    }
    existing = {
        row[0] for row in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        )
    }
    missing = sorted(required - existing)
    if missing:
        raise RuntimeError(
            "Database is missing tables: " + ", ".join(missing) +
            ". Run build_inspire_database.py first."
        )


def set_metadata(con: sqlite3.Connection, key: str, value: str) -> None:
    con.execute(
        "INSERT INTO metadata(key,value) VALUES(?,?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        (key, value),
    )


def read_obo_header(path: Path) -> Dict[str, str]:
    header: Dict[str, str] = {}
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith("["):
                break
            if ": " in line:
                key, value = line.split(": ", 1)
                if key in {"format-version", "data-version", "date", "ontology"}:
                    header[key] = value
    return header


def iter_obo_terms(path: Path) -> Iterator[Tuple[dict, List[Tuple[str, str]]]]:
    current: Optional[dict] = None
    parents: List[Tuple[str, str]] = []
    allowed = {
        "part_of", "regulates", "negatively_regulates",
        "positively_regulates"
    }

    def completed():
        if current and current.get("id", "").startswith("GO:"):
            return current, list(parents)
        return None

    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if line == "[Term]":
                done = completed()
                if done:
                    yield done
                current = {"is_obsolete": 0}
                parents = []
                continue
            if line.startswith("["):
                done = completed()
                if done:
                    yield done
                current = None
                parents = []
                continue
            if current is None or not line or line.startswith("!"):
                continue
            if line.startswith("id: "):
                current["id"] = line[4:].strip()
            elif line.startswith("name: "):
                current["name"] = line[6:].strip()
            elif line.startswith("namespace: "):
                current["namespace"] = line[11:].strip()
            elif line == "is_obsolete: true":
                current["is_obsolete"] = 1
            elif line.startswith("is_a: "):
                parents.append((line[6:].split()[0], "is_a"))
            elif line.startswith("relationship: "):
                fields = line[14:].split()
                if len(fields) >= 2 and fields[0] in allowed:
                    parents.append((fields[1], fields[0]))

        done = completed()
        if done:
            yield done


def load_obo(con: sqlite3.Connection, path: Path) -> None:
    header = read_obo_header(path)
    terms = []
    relations = []
    for term, parent_rows in iter_obo_terms(path):
        go_id = term["id"]
        terms.append((
            go_id, term.get("name"), term.get("namespace"),
            int(term.get("is_obsolete", 0))
        ))
        relations.extend((go_id, parent, relation)
                         for parent, relation in parent_rows)

    if not terms:
        raise RuntimeError("No GO terms found in the OBO file")

    known = {row[0] for row in terms}
    relations = [row for row in relations
                 if row[0] in known and row[1] in known]

    with con:
        con.execute("DELETE FROM go_pair_score")
        con.execute("DELETE FROM go_term_gene_count")
        con.execute("DELETE FROM protein_go_annotation")
        con.execute("DELETE FROM go_parent")
        con.execute("DELETE FROM go_term")
        con.executemany(
            "INSERT INTO go_term(go_id,name,namespace,is_obsolete) "
            "VALUES(?,?,?,?)", terms)
        con.executemany(
            "INSERT OR IGNORE INTO go_parent"
            "(child_go_id,parent_go_id,relation) VALUES(?,?,?)", relations)
        set_metadata(con, "go_obo_file", path.name)
        set_metadata(con, "go_obo_sha256", file_sha256(path))
        set_metadata(con, "go_ontology_version",
                     header.get("data-version", header.get("date", "unknown")))
        set_metadata(con, "go_obo_loaded_at", now_utc())

    print(f"Loaded {len(terms):,} GO terms and {len(relations):,} relationships")


def read_gaf_header(path: Path) -> Dict[str, str]:
    header: Dict[str, str] = {}
    with open_text(path) as handle:
        for raw in handle:
            if not raw.startswith("!"):
                break
            text = raw[1:].strip()
            if ": " in text:
                key, value = text.split(": ", 1)
                header[key] = value
    return header


def split_pipe(value: str) -> List[str]:
    return [x.strip() for x in value.split("|") if x.strip()]


def load_gaf(con: sqlite3.Connection, path: Path) -> None:
    network_genes = {row[0] for row in con.execute("SELECT gene_symbol FROM gene")}
    known_go = {row[0] for row in con.execute("SELECT go_id FROM go_term")}
    header = read_gaf_header(path)

    annotations: Set[Tuple[str, str, str, str, str, str, str]] = set()
    candidates: Dict[str, Set[str]] = defaultdict(set)
    data_rows = skipped_not = skipped_non_bp = unknown_go = 0

    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith("!") or not raw.strip():
                continue
            data_rows += 1
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 15:
                continue

            db = fields[0].strip()
            accession = fields[1].strip()
            symbol = normalize_gene(fields[2])
            qualifier = fields[3].strip()
            go_id = fields[4].strip()
            evidence = fields[6].strip()
            aspect = fields[8].strip()
            synonyms = [normalize_gene(x) for x in split_pipe(fields[10])]
            annotation_date = fields[13].strip()
            assigned_by = fields[14].strip()

            if db != "UniProtKB" or not accession:
                continue

            for name in {symbol, *synonyms}:
                if name in network_genes:
                    candidates[name].add(accession)

            qualifiers = set(split_pipe(qualifier))
            if "NOT" in qualifiers:
                skipped_not += 1
                continue
            if aspect != "P":
                skipped_non_bp += 1
                continue
            if go_id not in known_go:
                unknown_go += 1
                continue

            annotations.add((accession, go_id, evidence, qualifier, aspect,
                             annotation_date, assigned_by))

    release = header.get("generated-on", header.get("date-generated", "unknown"))
    mappings = []
    for gene in sorted(network_genes):
        accessions = sorted(candidates.get(gene, set()))
        for index, accession in enumerate(accessions):
            mappings.append((
                gene, accession, None, 1 if index == 0 else 0,
                "HUMAN-uniprot.gaf symbol/synonym", release
            ))

    with con:
        con.execute("DELETE FROM protein_go_annotation")
        con.execute(
            "DELETE FROM gene_uniprot "
            "WHERE mapping_source='HUMAN-uniprot.gaf symbol/synonym'")
        con.executemany(
            "INSERT OR REPLACE INTO protein_go_annotation"
            "(uniprot_accession,go_id,evidence_code,qualifier,aspect,"
            "annotation_date,assigned_by) VALUES(?,?,?,?,?,?,?)",
            sorted(annotations))
        con.executemany(
            "INSERT OR REPLACE INTO gene_uniprot"
            "(gene_symbol,uniprot_accession,reviewed,is_preferred,"
            "mapping_source,source_release) VALUES(?,?,?,?,?,?)", mappings)
        set_metadata(con, "go_gaf_file", path.name)
        set_metadata(con, "go_gaf_sha256", file_sha256(path))
        set_metadata(con, "go_annotation_version", release)
        set_metadata(con, "go_gaf_loaded_at", now_utc())
        set_metadata(con, "go_aspect_policy", "P biological_process only")
        set_metadata(con, "go_not_policy", "NOT annotations excluded")

    print(f"Read {data_rows:,} GAF data rows")
    print(f"Loaded {len(annotations):,} GO Biological Process annotations")
    print(f"Mapped {len(candidates):,} of {len(network_genes):,} network genes")
    print(f"Skipped NOT annotations: {skipped_not:,}")
    print(f"Skipped non-BP annotations: {skipped_non_bp:,}")
    print(f"GO IDs absent from supplied OBO: {unknown_go:,}")


def write_unmapped_report(con: sqlite3.Connection, path: Path) -> None:
    rows = con.execute(
        "SELECT g.gene_symbol FROM gene g "
        "LEFT JOIN gene_uniprot u ON u.gene_symbol=g.gene_symbol "
        "WHERE u.uniprot_accession IS NULL ORDER BY g.gene_symbol"
    ).fetchall()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("gene_symbol\n")
        for (gene,) in rows:
            handle.write(gene + "\n")
    print(f"Wrote {len(rows)} unmapped gene(s) to {path}")


def fetch_json(url: str, timeout: int, retries: int) -> dict:
    last_error: Optional[Exception] = None
    for attempt in range(retries + 1):
        request = urllib.request.Request(
            url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.load(response)
        except (urllib.error.URLError, urllib.error.HTTPError,
                TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(min(2 ** attempt, 30))
    raise RuntimeError(str(last_error))


def pfam_accessions(payload: dict) -> Set[str]:
    found: Set[str] = set()
    for item in payload.get("results", []):
        accession = (item.get("metadata") or {}).get("accession")
        if isinstance(accession, str) and accession.upper().startswith("PF"):
            found.add(accession.upper())
    return found


def load_pfam(con: sqlite3.Connection, delay: float,
               timeout: int, retries: int, refresh: bool) -> None:
    accessions = [row[0] for row in con.execute(
        "SELECT DISTINCT uniprot_accession FROM gene_uniprot "
        "WHERE is_preferred=1 ORDER BY uniprot_accession")]
    if not accessions:
        raise RuntimeError("No preferred UniProt mappings. Load the GAF first.")

    errors = 0
    for number, accession in enumerate(accessions, start=1):
        status_row = con.execute(
            "SELECT pfam_status FROM annotation_status "
            "WHERE uniprot_accession=?", (accession,)).fetchone()
        if not refresh and status_row and status_row[0] in {
            "complete", "complete_empty"
        }:
            print(f"[{number}/{len(accessions)}] {accession}: cached")
            continue

        print(f"[{number}/{len(accessions)}] {accession}: querying InterPro")
        url: Optional[str] = (
            f"{INTERPRO_BASE}/entry/pfam/protein/uniprot/"
            f"{urllib.parse.quote(accession, safe='')}/?page_size=200")
        found: Set[str] = set()
        try:
            while url:
                payload = fetch_json(url, timeout, retries)
                found.update(pfam_accessions(payload))
                url = payload.get("next") or None
                if url:
                    time.sleep(delay)

            with con:
                if refresh:
                    con.execute("DELETE FROM protein_pfam WHERE uniprot_accession=?",
                                (accession,))
                con.executemany(
                    "INSERT OR IGNORE INTO protein_pfam"
                    "(uniprot_accession,pfam_accession,source,source_release,"
                    "retrieved_at) VALUES(?,?,'InterPro API',NULL,?)",
                    [(accession, pfam, now_utc()) for pfam in sorted(found)])
                status = "complete" if found else "complete_empty"
                con.execute(
                    "INSERT INTO annotation_status"
                    "(uniprot_accession,pfam_status,go_status,last_error,updated_at) "
                    "VALUES(?,?,'pending',NULL,?) "
                    "ON CONFLICT(uniprot_accession) DO UPDATE SET "
                    "pfam_status=excluded.pfam_status,last_error=NULL,"
                    "updated_at=excluded.updated_at",
                    (accession, status, now_utc()))
            print(f"    found {len(found)} Pfam accession(s)")
        except Exception as exc:
            errors += 1
            with con:
                con.execute(
                    "INSERT INTO annotation_status"
                    "(uniprot_accession,pfam_status,go_status,last_error,updated_at) "
                    "VALUES(?,'error','pending',?,?) "
                    "ON CONFLICT(uniprot_accession) DO UPDATE SET "
                    "pfam_status='error',last_error=excluded.last_error,"
                    "updated_at=excluded.updated_at",
                    (accession, str(exc), now_utc()))
            print(f"    ERROR: {exc}", file=sys.stderr)
        time.sleep(delay)

    with con:
        set_metadata(con, "pfam_source", "InterPro API")
        set_metadata(con, "pfam_last_run_at", now_utc())
    print(f"Pfam retrieval complete with {errors} error(s)")


def verify(con: sqlite3.Connection) -> None:
    integrity = con.execute("PRAGMA integrity_check").fetchone()[0]
    if integrity != "ok":
        raise RuntimeError(f"SQLite integrity check failed: {integrity}")
    fk = con.execute("PRAGMA foreign_key_check").fetchall()
    if fk:
        raise RuntimeError(f"SQLite foreign-key check found {len(fk)} problem(s)")


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Load GO files and optional Pfam annotations into SQLite")
    parser.add_argument("--db", required=True, type=Path)
    parser.add_argument("--obo", type=Path)
    parser.add_argument("--gaf", type=Path)
    parser.add_argument("--skip-pfam", action="store_true")
    parser.add_argument("--pfam-only", action="store_true")
    parser.add_argument("--refresh-pfam", action="store_true")
    parser.add_argument("--api-delay", type=float, default=0.25)
    parser.add_argument("--api-timeout", type=int, default=60)
    parser.add_argument("--api-retries", type=int, default=3)
    parser.add_argument("--unmapped-report", type=Path)
    return parser.parse_args()


def main() -> int:
    args = arguments()
    db_path = args.db.expanduser().resolve()
    if not db_path.is_file():
        print(f"ERROR: database not found: {db_path}", file=sys.stderr)
        return 2
    if args.skip_pfam and args.pfam_only:
        print("ERROR: do not combine --skip-pfam and --pfam-only", file=sys.stderr)
        return 2

    if not args.pfam_only:
        if not args.obo or not args.gaf:
            print("ERROR: --obo and --gaf are required", file=sys.stderr)
            return 2
        obo = args.obo.expanduser().resolve()
        gaf = args.gaf.expanduser().resolve()
        if not obo.is_file() or not gaf.is_file():
            print("ERROR: OBO or GAF file not found", file=sys.stderr)
            return 2

    con = connect_db(db_path)
    try:
        require_schema(con)
        with con:
            set_metadata(con, "annotation_loader_version", "1.0")
            set_metadata(con, "annotation_loader_last_run_at", now_utc())

        if not args.pfam_only:
            print(f"Loading ontology: {obo}")
            load_obo(con, obo)
            print(f"Loading annotations: {gaf}")
            load_gaf(con, gaf)

        report = (args.unmapped_report.expanduser().resolve()
                  if args.unmapped_report else
                  db_path.with_name("unmapped_genes.csv"))
        write_unmapped_report(con, report)

        if not args.skip_pfam:
            load_pfam(con, max(0.0, args.api_delay),
                       max(1, args.api_timeout), max(0, args.api_retries),
                       args.refresh_pfam)

        verify(con)
        print("Annotation load complete. SQLite integrity check: OK")
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    finally:
        con.close()


if __name__ == "__main__":
    raise SystemExit(main())
