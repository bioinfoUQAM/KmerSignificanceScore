#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Describe the compiled pipeline outputs the repository ships, file by file.

These twelve files are the per-position output of `main.py`, one per analysed gene. They are the
input the crosscheck reads to verify S10 to S12 cell by cell, the per-gene indel statistics, the
windowing grid and the twenty-six KSS values the concordance quotes, and the input
`generate_kss_figure.ipynb` builds those three workbooks and Figure 3 from.

They were ignored by `.gitignore` until 29 August 2026. The consequence was measurable: on a fresh
clone the crosscheck stopped at 378 checks of 764 with a FileNotFoundError, so the "764 green" the
revision reports could only ever be reproduced on the author's machine, while Data availability
promised "the evaluation and analysis scripts that reproduce every reported result".

What ships is derived: a score and a class-count table per position. No sequence, no accession, no
local path. The uncompiled `*_results.json` beside them stay out, at 7.8 GB for the twelve.

Usage:
    python scripts/build_compiled_results_manifest.py
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PATTERN = "data/*/*/results/*_compiled_results.json"
OUTPUT = ROOT / "data" / "compiled_results_manifest.csv"
PRODUCER = "main.py"

COLUMNS = [
    "path", "virus", "gene", "bytes", "sha256", "positions", "variants", "producer",
]


def describe(path: Path) -> dict:
    """One manifest row: identity, size, digest, and what the file holds."""
    raw = path.read_bytes()
    payload = json.loads(raw.decode("utf-8"))
    gene = path.stem.replace("_compiled_results", "")
    inner = payload.get(gene, payload)
    scored = {p: d for p, d in inner.items()
              if p.isdigit() and isinstance(d, dict) and "alts" in d}
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "virus": path.relative_to(ROOT).parts[1],
        "gene": gene,
        "bytes": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "positions": len(scored),
        "variants": sum(len(d["alts"]) for d in scored.values()),
        "producer": PRODUCER,
    }


def main() -> int:
    """Write the manifest and report what it covers."""
    paths = sorted(ROOT.glob(PATTERN))
    if not paths:
        print(f"no file matches {PATTERN}")
        return 1

    rows = [describe(path) for path in paths]
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    total = sum(row["bytes"] for row in rows)
    print(f"{len(rows)} files, {total / 1024 / 1024:.1f} MB, "
          f"{sum(r['positions'] for r in rows)} scored positions")
    for row in rows:
        print(f"  {row['path']:78} {row['positions']:5} positions  {row['sha256'][:12]}")
    print(f"\nwritten to {OUTPUT.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
