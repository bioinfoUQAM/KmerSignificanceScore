"""
Check that the component weights the pipeline accepts are the ones KSS is defined on.

Materials and methods states two things about them: they are non-negative and not all zero,
and KSS therefore lies in [0, 1] because the three component scores do. Neither held in the
code. `compute_kss_scores` divides by the sum of the weights, so three zeros raised
ZeroDivisionError deep inside a run, and a negative weight produced a weighted average outside
the interval the manuscript guarantees, silently, with no error at all.

The guarantee matters because it is what lets a KSS be read without knowing the run: a reader
who sees 0.82 does not have to ask whether the weights were admissible. A guarantee the code
does not enforce is worse than none, since it is the one a reader relies on.

This test states the property from both sides: inadmissible weights are refused when the
parameters are read, and admissible ones keep the average inside [0, 1] for every combination
of component scores at the bounds.

Run directly:
    python tests/test_weight_bounds.py
or under pytest:
    pytest tests/test_weight_bounds.py
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.kss import _compute_final_kss  # noqa: E402
from src.mutation_score import get_mutational_scores  # noqa: E402
from src.pipeline import validate_component_weights, validate_indel_score  # noqa: E402


def test_inadmissible_weights_are_refused() -> None:
    """Negative weights and an all-zero triple must fail when the parameters are read."""
    refused = 0
    for weights in ((-1, 1, 1), (1, -0.5, 1), (1, 1, -3), (0, 0, 0)):
        block = dict(zip(("discriminative", "mutational", "protein"), weights))
        try:
            validate_component_weights(block)
        except ValueError:
            refused += 1
        except Exception as error:                       # noqa: BLE001
            raise AssertionError(
                f"weights {weights} raised {type(error).__name__} rather than ValueError: "
                f"{error}"
            ) from error
        else:
            raise AssertionError(f"weights {weights} were accepted; they are not admissible")
    assert refused == 4, f"{refused} of 4 inadmissible weight sets were refused"


def test_admissible_weights_keep_kss_inside_the_unit_interval() -> None:
    """With admissible weights, the average of bounded components stays bounded."""
    admissible = [(1, 1, 1), (1, 0, 0), (0, 1, 0), (0, 0, 1), (3, 1, 0.5), (0.1, 0.1, 10)]
    for weights in admissible:
        parameters = {"discriminative_weight": weights[0], "mutational_weight": weights[1],
                      "protein_weight": weights[2]}
        weight_sum = sum(weights)
        for scores in itertools.product((0.0, 0.5, 1.0), repeat=3):
            details = {"discriminative_score": scores[0], "mutational_score": scores[1],
                       "protein_score": scores[2]}
            kss = _compute_final_kss(details, parameters, weight_sum)
            assert 0.0 <= kss <= 1.0, (
                f"weights {weights} and scores {scores} gave KSS {kss}, outside [0, 1]"
            )


def test_inadmissible_indel_scores_are_refused() -> None:
    """An indel score outside [0, 1] must fail, at the config and at the scoring function.

    An indel takes this value directly and the mutational score is the largest impact at a
    window, so 2.0 gave a KSS of 1.333 at equal weights and -1.0 gave -1.0, silently.
    """
    for indel_score in (2.0, -1.0, 1.5, -0.1):
        for name, call in (("validate_indel_score", lambda v=indel_score: validate_indel_score(v)),
                           ("get_mutational_scores",
                            lambda v=indel_score: get_mutational_scores(["A10-"], indel_score=v))):
            try:
                call()
            except ValueError:
                continue
            raise AssertionError(f"{name} accepted indel_score {indel_score}, outside [0, 1]")


def test_admissible_indel_scores_are_accepted() -> None:
    """The bound must not refuse the values the framework is defined on."""
    for indel_score in (0.0, 0.4, 1.0):
        validate_indel_score(indel_score)
        scores = get_mutational_scores(["A10-"], indel_score=indel_score)
        assert scores["A10-"] == indel_score, (
            f"an indel scored {scores['A10-']} where the configuration asked for {indel_score}"
        )


if __name__ == "__main__":
    test_inadmissible_weights_are_refused()
    print("inadmissible weights are refused when the parameters are read")
    test_admissible_weights_keep_kss_inside_the_unit_interval()
    print("admissible weights keep KSS inside [0, 1]")
    test_inadmissible_indel_scores_are_refused()
    print("an indel score outside [0, 1] is refused by the config and by the scorer")
    test_admissible_indel_scores_are_accepted()
    print("an admissible indel score reaches the mutational score unchanged")
