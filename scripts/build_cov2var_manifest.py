"""Describe the seven COV2Var downloads the error-rate analysis reads.

The files are not redistributed: the download page publishes no licence and carries an
all-rights-reserved notice. This manifest records what was read, so that a future download
can be compared against it and a reader can tell which snapshot the published numbers came
from. `results/functional_benchmark/reference_matrix.csv` carries the derived results, and
it is committed, so the published figures stay checkable without the inputs.

Digests are recomputed from the bytes on disk and then confronted with COV2VAR_SHA256 in
`scripts/external_inputs.py`, which stays the repository's single source of truth. The
script fails rather than write a manifest that would contradict those constants.

The date column is the file's modification time. It is the only datable fact available here
and says when the bytes were written locally, not when the server served them.

Usage
-----
    python scripts/build_cov2var_manifest.py
"""

from __future__ import annotations

import csv
import hashlib
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

from external_inputs import COV2VAR_SHA256, COV2VAR_URL  # noqa: E402

SOURCE = ROOT / "data" / "external" / "COV2Var"
OUTPUT = ROOT / "results" / "functional_benchmark" / "cov2var_input_manifest.csv"

USED_FOR = {
    "mutations_summary.txt":
        "common-mutation membership layer, and the mutation identifiers the other layers "
        "are joined on",
    "results of natural selection analysis.zip":
        "MEME, FEL and FUBAR selection layers, taken singly, in union and in intersection",
    "rbd_immune_escape_score.txt": "deep mutational scanning layer, antibody escape",
    "RBD binding affinity.txt": "deep mutational scanning layer, ACE2 binding",
    "mutation induced protein function changes.txt":
        "MutPred2 prediction layer, reported and not admitted as evidence",
    "mutation induced protein stability changes.txt":
        "iMutant prediction layer, reported and not admitted as evidence",
    "mutation effects on antigenicity and immunogenicity.txt":
        "VaxiJen prediction layer, reported and not admitted as evidence",
}
CITATION = ("Feng et al. COV2Var, a function annotation database of SARS-CoV-2 genetic "
            "variants. Nucleic Acids Research 2024;52(D1):D701-D713.")
DOI = "10.1093/nar/gkad958"
LICENCE = ('no explicit licence published; the download page carries '
           '"COV2Var. All Rights Reserved"')
FIELDS = ["file", "bytes", "sha256", "file_mtime_utc", "source_page", "used_for",
          "redistribution", "licence", "citation", "doi"]


def digest(handle) -> str:
    """Streamed, because one of these files is 5 MB and another could grow."""
    sha = hashlib.sha256()
    for block in iter(lambda: handle.read(1 << 20), b""):
        sha.update(block)
    return sha.hexdigest()


def row(name: str, size: int, sha: str, day: str, used: str) -> dict:
    return {"file": name, "bytes": size, "sha256": sha, "file_mtime_utc": day,
            "source_page": COV2VAR_URL, "used_for": used,
            "redistribution": "not redistributed", "licence": LICENCE,
            "citation": CITATION, "doi": DOI}


def main() -> int:
    rows, disagreements = [], []
    for name in sorted(COV2VAR_SHA256):
        path = SOURCE / name
        if not path.exists():
            disagreements.append(f"{name}: absent from {SOURCE}")
            continue
        with path.open("rb") as handle:
            sha = digest(handle)
        if sha != COV2VAR_SHA256[name]:
            disagreements.append(
                f"{name}: {sha[:12]} on disk against {COV2VAR_SHA256[name][:12]} pinned")
        stat = path.stat()
        day = datetime.fromtimestamp(stat.st_mtime, timezone.utc).strftime("%Y-%m-%d")
        rows.append(row(name, stat.st_size, sha, day, USED_FOR[name]))

        # The selection layers arrive inside an archive, so its members carry their own
        # digests: a repacked archive with the same contents would otherwise look changed,
        # and a changed member inside an identical wrapper would look unchanged.
        if path.suffix == ".zip":
            with zipfile.ZipFile(path) as archive:
                for info in sorted(archive.infolist(), key=lambda i: i.filename):
                    if info.is_dir():
                        continue
                    with archive.open(info) as member:
                        rows.append(row(f"{name}::{info.filename}", info.file_size,
                                        digest(member), day,
                                        "member of the selection archive above"))

    if disagreements:
        print("The manifest would contradict COV2VAR_SHA256:")
        for problem in disagreements:
            print(f"  {problem}")
        return 1

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    downloads = [r for r in rows if "::" not in r["file"]]
    print(f"{len(downloads)} downloads and {len(rows) - len(downloads)} archive members, "
          f"{sum(r['bytes'] for r in downloads) / 1048576:.1f} MB")
    print(f"all {len(downloads)} digests agree with COV2VAR_SHA256")
    print(f"written to {OUTPUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
