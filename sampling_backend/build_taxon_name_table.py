#!/usr/bin/env python3
"""Build a compact SQLite table of taxon name -> tax_id from taxonomy JSONL."""

from __future__ import annotations

import argparse
import json
import os
import sqlite3
import time


def normalise(name: str) -> str:
    return name.strip().lower()


def create_schema(conn: sqlite3.Connection) -> None:
    conn.execute("DROP TABLE IF EXISTS taxon_name_map")
    conn.execute(
        """
        CREATE TABLE taxon_name_map (
            name TEXT NOT NULL,
            name_normalized TEXT NOT NULL,
            tax_id INTEGER NOT NULL,
            UNIQUE(name_normalized, tax_id)
        )
        """
    )
    conn.execute("CREATE INDEX idx_taxon_name_map_name_norm ON taxon_name_map(name_normalized)")
    conn.execute("CREATE INDEX idx_taxon_name_map_taxid ON taxon_name_map(tax_id)")
    conn.commit()


def flush_batch(conn: sqlite3.Connection, batch: list[tuple[str, str, int]]) -> int:
    if not batch:
        return 0
    conn.executemany(
        """
        INSERT OR IGNORE INTO taxon_name_map (name, name_normalized, tax_id)
        VALUES (?, ?, ?)
        """,
        batch,
    )
    conn.commit()
    inserted = len(batch)
    batch.clear()
    return inserted


def build(input_path: str, output_path: str, batch_size: int, progress_every: int) -> None:
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input taxonomy JSONL not found: {input_path}")

    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with sqlite3.connect(output_path) as conn:
        create_schema(conn)
        start = time.time()
        lines = 0
        queued = 0
        batch: list[tuple[str, str, int]] = []

        with open(input_path, "r", encoding="utf-8") as handle:
            for line in handle:
                lines += 1
                if not line.strip():
                    continue
                obj = json.loads(line)
                taxonomy = obj.get("taxonomy", {})
                tax_id = taxonomy.get("tax_id")
                if not tax_id:
                    continue

                sci_name = (taxonomy.get("current_scientific_name") or {}).get("name")
                if sci_name:
                    n = normalise(sci_name)
                    if n:
                        batch.append((sci_name, n, int(tax_id)))

                for syn in (taxonomy.get("synonyms") or []):
                    if not syn:
                        continue
                    n = normalise(syn)
                    if n:
                        batch.append((syn, n, int(tax_id)))

                if len(batch) >= batch_size:
                    queued += flush_batch(conn, batch)

                if progress_every and lines % progress_every == 0:
                    elapsed = time.time() - start
                    print(f"Processed {lines:,} lines in {elapsed:.1f}s")

        queued += flush_batch(conn, batch)

        elapsed = time.time() - start
        row_count = conn.execute("SELECT COUNT(*) FROM taxon_name_map").fetchone()[0]
        print(f"Done in {elapsed:.1f}s")
        print(f"Input lines: {lines:,}")
        print(f"Inserted attempts: {queued:,}")
        print(f"Unique rows in table: {row_count:,}")
        print(f"SQLite DB: {output_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a compact taxon name->tax_id SQLite table from taxonomy JSONL."
    )
    parser.add_argument("input_jsonl", help="Path to taxonomy JSONL (e.g. taxonomy_all_new.jsonl)")
    parser.add_argument("output_sqlite", help="Output SQLite file path (e.g. taxonomy_name_taxid.sqlite3)")
    parser.add_argument("--batch-size", type=int, default=50_000, help="SQLite insert batch size")
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1_000_000,
        help="Print progress every N input lines (0 disables progress prints)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build(
        input_path=args.input_jsonl,
        output_path=args.output_sqlite,
        batch_size=max(1, int(args.batch_size)),
        progress_every=max(0, int(args.progress_every)),
    )
