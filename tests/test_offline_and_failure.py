"""
Check that a clone runs its quick start without the network, and that a missing input fails
loudly instead of being reported as a success.

Two defects motivate this file. The protein component reached UniProt on every call: the entry
cache is keyed by accession, which only the search returns, so the twelve entries the repository
ships were unreachable and `python main.py data/toy/config.yaml` failed offline on a protein it
already held. And a gene whose CDS or FASTA is absent was skipped with a warning, then saved as
an empty file and counted in "Successful: n/n", with the process still exiting 0, so neither a
reader nor a CI job could tell an empty artefact from a reproduction.

Run directly:
    python tests/test_offline_and_failure.py
or under pytest:
    pytest tests/test_offline_and_failure.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src import protein_score  # noqa: E402

CACHE = ROOT / "src" / "uniprot"
INDEX = CACHE / "index.json"


def _import_main():
    """
    Import main.py and check that importing it left the caller's streams alone.

    main.py forces UTF-8 on the Windows console at import time. It has to do that in place: a
    wrapper built around sys.stdout.buffer closes a buffer it did not open, which under a test
    harness is the harness's own and surfaces as "I/O operation on closed file" at teardown.
    This asserts the property rather than working around it, so a return to the wrapper form
    fails here instead of in whatever runs next.
    """
    before = (sys.stdout, sys.stderr)
    import main  # noqa: PLC0415

    assert (sys.stdout, sys.stderr) == before, (
        "importing main.py replaced the caller's streams; it must reconfigure them in place")
    return main


def test_index_covers_every_reference_the_repository_ships() -> int:
    """
    Every (taxon, gene) pair a configured dataset resolves must be in the shipped index, and
    the entry it names must be in the cache.

    The pair is what has to be indexed, not the entry: the HCMV references carry taxon 10359
    while their cached entries carry 10360, and UL55 resolves to an entry whose gene name is
    gB, so the mapping cannot be rebuilt from the cached JSON files.
    """
    from src import kss  # noqa: PLC0415

    index = json.loads(INDEX.read_text(encoding="utf-8"))
    checked = 0
    for gb in sorted(ROOT.glob("data/*/*/*.gb")):
        gene = gb.parent.name
        taxon = kss.get_taxon_id(str(gb))
        assert taxon, f"{gb} carries no taxonomy identifier"
        key = f"{taxon}|{gene}"
        assert key in index, f"{key} is missing from {INDEX.name}, so it needs the network"
        entry = CACHE / f"{index[key]}.json"
        assert entry.exists(), f"{key} points at {index[key]}, which the cache does not hold"
        checked += 1

    assert checked > 0, "no reference genome found to check"
    print(f"{checked} references indexed")


def test_resolution_uses_the_index_and_not_the_network() -> None:
    """
    A pair in the index must resolve without any HTTP call.

    requests.get is replaced by a function that fails, so a resolution that still reaches the
    network fails the test rather than silently succeeding on a machine that happens to be
    online. That is the whole point: the defect was invisible with a connection.
    """
    def refuse(*args, **kwargs):  # noqa: ANN002, ANN003, ARG001
        raise AssertionError("the resolution reached the network despite the shipped index")

    original = protein_score.requests.get
    protein_score.requests.get = refuse
    try:
        entry = protein_score.fetch_uniprot_data("10359", "US28")
    finally:
        protein_score.requests.get = original

    assert entry == "P69332", f"expected P69332 from the index, got {entry!r}"


def test_an_unknown_pair_still_goes_to_the_search() -> None:
    """
    The index is a cache, not a whitelist: a pair it does not carry must fall through to the
    search rather than be reported as absent.
    """
    calls = []

    def refuse(*args, **kwargs):  # noqa: ANN002, ANN003, ARG001
        calls.append(kwargs.get("params", {}))
        raise RuntimeError("no network in this test")

    original = protein_score.requests.get
    protein_score.requests.get = refuse
    try:
        protein_score.fetch_uniprot_data("9606", "TP53")
    finally:
        protein_score.requests.get = original

    assert calls, "an unindexed pair must still reach the search"


def test_a_gene_without_input_is_an_error() -> None:
    """
    process_single_gene must raise when the analysis returned nothing for the gene, so that no
    empty file is written and the dataset is counted as failed.
    """
    main = _import_main()

    original = main.kanalyzer.analyze_records
    main.kanalyzer.analyze_records = lambda info, parameters: {"_metadata": {}, "genes": {}}
    try:
        raised = None
        try:
            main.process_single_gene(
                "UL55", "a dataset with no input",
                {"input_folder": str(ROOT / "data"), "output_folder": str(ROOT / "data")},
                {"k": 9}, 1, 1)
        except RuntimeError as exc:
            raised = exc
    finally:
        main.kanalyzer.analyze_records = original

    assert raised is not None, "a gene with no analyzed sequence was not reported as an error"
    assert "No sequences were analyzed" in str(raised), f"unexpected message: {raised}"


def test_main_reports_failure_through_its_exit_code() -> None:
    """
    main must return a non-zero value when a dataset fails, and the entry point must carry it
    to the process, otherwise no script or CI job wrapping the tool can detect the failure.
    """
    source = (ROOT / "main.py").read_text(encoding="utf-8")
    assert "raise SystemExit(main())" in source, (
        "main.py must hand main()'s value to the process")
    assert "return 1" in source and "return 0" in source, (
        "main() must return a status rather than None")


if __name__ == "__main__":
    test_index_covers_every_reference_the_repository_ships()
    test_resolution_uses_the_index_and_not_the_network()
    print("a shipped pair resolves without the network")
    test_an_unknown_pair_still_goes_to_the_search()
    print("an unindexed pair still reaches the search")
    test_a_gene_without_input_is_an_error()
    print("a gene without input raises instead of saving an empty file")
    test_main_reports_failure_through_its_exit_code()
    print("failure reaches the exit code")
