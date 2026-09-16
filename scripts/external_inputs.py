"""
The external files these analyses read, pinned by checksum.

Comment 4.1 asks the repository to guarantee reproducibility. Pinning the Python interpreter
and the package bounds is not enough when an analysis reads a file that lives upstream: the
fitness estimates were taken from a branch tip, and a branch tip is not a version. This
pins the two external sources the published analyses read, so a different download cannot be
substituted silently.

Checksums are verified at startup by the scripts that read these files. A mismatch stops the
run rather than producing numbers that look like the published ones.

Sources
-------
aa_fitness.csv
    jbloomlab/SARS2-mut-fitness, results/aa_fitness/aa_fitness.csv, at commit 067fce1 of
    2024-11-08, the last change to that path. Cite the permalink rather than main:
    https://raw.githubusercontent.com/jbloomlab/SARS2-mut-fitness/067fce16be8c6bed5cf47b5189e8558c6a058697/results/aa_fitness/aa_fitness.csv

COV2Var
    biomedbdc.wchscu.cn/COV2Var/download/. The site serves no version tag, so only the
    checksums identify what was used. Seven files are read: four the error-rate analysis admits
    as evidence, and three prediction layers it reports without admitting. The rest of the
    download is not read.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

FITNESS_URL = ("https://raw.githubusercontent.com/jbloomlab/SARS2-mut-fitness/"
               "067fce16be8c6bed5cf47b5189e8558c6a058697/results/aa_fitness/aa_fitness.csv")
COV2VAR_URL = "https://biomedbdc.wchscu.cn/COV2Var/download/"

FITNESS_SHA256 = "4d0ea6ce6fde3451698b0bafdb608ccf750c41e3a6acc62a108c1161db077fb0"

COV2VAR_SHA256 = {
    "mutations_summary.txt":
        "3dfd94f1eebe25b8df690eff249b3f688da02da4413c431344e3e8f7e20f03f0",
    "rbd_immune_escape_score.txt":
        "ac74902b2b733cbdf7864a29c896fc6894d71bcfd0f3ef58dec8ac86462c35c2",
    "RBD binding affinity.txt":
        "89190fa9a0d867b843b7806030facfb416f385499598b7be2c4cc9a566d32ee1",
    "results of natural selection analysis.zip":
        "93c4d81faf78e03eb5fc7a132cc3c2c575ea689723b414f3dc89fc26b0bd5e08",
    "mutation induced protein function changes.txt":
        "e552838214000b16cc68a4ee354c61c5e7f497fb13fc1dee3785e567db586bd3",
    "mutation induced protein stability changes.txt":
        "f1f3deefb5a2993bbde3584076fe6f58c2cbdd2271d410e6378a087865f865d3",
    "mutation effects on antigenicity and immunogenicity.txt":
        "94200c2360f8f96dce159bb63970489cd6608b73c0ac7bec9df1c0c5e3de9fdc",
}


def _digest(path: Path) -> str:
    """Streamed, because one of these files is 5 MB and another could grow."""
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            sha.update(block)
    return sha.hexdigest()


def _check(path: Path, expected: str, source: str) -> str | None:
    if not path.exists():
        return f"missing input: {path}\n  download it from {source}"
    found = _digest(path)
    if found != expected:
        return (f"{path.name} is not the file this analysis was run on\n"
                f"  expected {expected}\n"
                f"  found    {found}\n"
                f"  re-download it from {source}")
    return None


def verify_fitness(path: Path) -> str | None:
    """None when the fitness file is the pinned one, otherwise the message to print."""
    return _check(path, FITNESS_SHA256, FITNESS_URL)


def verify_cov2var(directory: Path, names: list[str] | None = None) -> str | None:
    """None when every named COV2Var file matches, otherwise the first mismatch."""
    for name in names if names is not None else sorted(COV2VAR_SHA256):
        problem = _check(directory / name, COV2VAR_SHA256[name], COV2VAR_URL)
        if problem:
            return problem
    return None
