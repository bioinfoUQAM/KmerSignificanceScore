"""
Check that a substitution matrix survives the format it is stored in, all the way to the
mutational scores the framework reports.

Comparing the stored values alone is not sufficient: the values pass through an inversion for
distance matrices, a global min-max rescaling, and a final rounding to three decimals. A
difference could in principle be absorbed by the rounding, or amplified by the rescaling if it
affected the extreme values.

Background: the matrices were once stored as pickles alongside a JSON copy written from
rounded values, which differed by up to one unit in the last place and is why the loader
preferred the pickle "to preserve numerical precision". Comment 4.2 removed the pickle, from
the files and from the loader, so JSON is now the only format read or written and this test
guards the property the removal depends on: writing and reading it back changes nothing.

Run directly:
    python tests/test_matrix_formats.py
or under pytest:
    pytest tests/test_matrix_formats.py
"""

from __future__ import annotations

import itertools
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.mutation_score import (  # noqa: E402
    adjust_and_scale_substitution_matrix,
    get_mutational_scores,
    load_substitution_matrix,
    save_substitution_matrix,
)

MATRICES = ROOT / "src" / "substitution_matrices"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# How much each check covered, for the summary the direct run prints. It is recorded here
# rather than returned: pytest treats a test that returns a value as an error from version 8,
# and returning the count was the reason this file raised PytestReturnNotNoneWarning.
CHECKED: dict = {}


def test_json_round_trip_is_lossless() -> None:
    """
    Every shipped matrix must survive a write-read cycle unchanged, going through the
    project's own writer.

    JSON being the only format the project reads and writes, this is the check that the
    stored values survive it. It must exercise save_substitution_matrix rather than
    json.dump, because the defect it guards against lived in the writer: an earlier version
    of _save_matrix_json rounded to 15 decimal places, which altered 216 of the 400 entries
    of MIYATA_EVO by up to 4e-16. A test that wrote the file itself would not have caught it.
    """
    checked = 0
    for path in sorted(MATRICES.glob("*.json")):
        original = load_substitution_matrix(str(path))
        with tempfile.TemporaryDirectory() as tmp:
            base = str(Path(tmp) / path.stem)
            written = save_substitution_matrix(original, base)
            reloaded = load_substitution_matrix(written)

        differing = {k for k in original if original[k] != reloaded[k]}
        assert not differing, (
            f"{path.name}: {len(differing)} of {len(original)} values altered by the "
            f"write-read cycle, worst "
            f"{max(abs(original[k] - reloaded[k]) for k in differing):.3e}"
        )
        checked += 1

    assert checked > 0, "no matrix found to check"
    CHECKED["matrices"] = checked


def test_round_trip_survives_to_the_reported_scores(name: str = "MIYATA_EVO") -> None:
    """
    The same cycle, carried through the two transformations that stand between the stored
    value and the number the framework reports: inversion and rescaling, then scoring.

    A difference too small to see in the stored values can still be amplified by the global
    min-max rescaling if it touches an extreme value, so the comparison is made after it and
    on the final scores, for every substitution and for indels.
    """
    path = MATRICES / f"{name}.json"
    original = load_substitution_matrix(str(path))
    with tempfile.TemporaryDirectory() as tmp:
        written = save_substitution_matrix(original, str(Path(tmp) / name))
        reloaded = load_substitution_matrix(written)

    scaled_a = adjust_and_scale_substitution_matrix(original, name)
    scaled_b = adjust_and_scale_substitution_matrix(reloaded, name)
    assert set(scaled_a) == set(scaled_b)
    differing = {k for k in scaled_a if scaled_a[k] != scaled_b[k]}
    assert not differing, f"{len(differing)} scaled values differ after the round trip"

    mutations = [f"{a}1{b}" for a, b in itertools.permutations(AMINO_ACIDS, 2)]
    mutations += [f"{a}1-" for a in AMINO_ACIDS] + [f"-1{a}" for a in AMINO_ACIDS]

    scores_a = get_mutational_scores(mutations, name, custom_matrix=original)
    scores_b = get_mutational_scores(mutations, name, custom_matrix=reloaded)
    assert set(scores_a) == set(scores_b)
    differing = {m: (scores_a[m], scores_b[m])
                 for m in scores_a if scores_a[m] != scores_b[m]}
    assert not differing, (
        f"{len(differing)} of {len(mutations)} reported scores differ: "
        f"{list(differing.items())[:5]}"
    )
    CHECKED["mutations"] = len(mutations)


def main() -> int:
    print("Substitution matrix format checks\n")

    failures = 0

    for label, check, unit in [
        ("JSON round trip is lossless", test_json_round_trip_is_lossless, "matrices"),
        ("round trip to reported scores", test_round_trip_survives_to_the_reported_scores,
         "mutations"),
    ]:
        try:
            check()
            print(f"  OK      {label:<32} {CHECKED.get(unit, 0)} {unit}")
        except AssertionError as exc:
            failures += 1
            print(f"  FAILED  {label:<32} {exc}")

    print()
    print("All checks passed." if not failures else f"{failures} check(s) failed.")
    return failures


if __name__ == "__main__":
    sys.exit(main())
