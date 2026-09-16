"""Check the four virus-derived matrices by routes independent of the ones that built them.

`build_viral_matrices.py` already refuses to write a matrix whose exchangeabilities differ
between IQ-TREE and RAxML. That guards the sources. It does not guard the conversion from a
rate model to a log-odds score matrix, nor the composite that places the result in the ranking,
and those are the two steps that decide what goes in the manuscript.

Four checks, each reaching the same quantity by a different path:

  1. the rate matrix is a valid generator     rows of Q sum to zero, the model is normalised to
                                              one expected substitution per site, and P(t) is
                                              stochastic with non-negative entries
  2. the model is reversible                  detailed balance holds, which is what makes the
                                              log-odds symmetric; asserted rather than assumed
  3. the exponential is not an artefact       P(t) recomputed by eigendecomposition instead of
                                              scipy's expm, a different algorithm
  4. the composite is not an artefact         recomputed from the stored JSON entries instead of
                                              through the pipeline's mutational scores, using
                                              the fact that impact is a strictly decreasing
                                              transform of the entry so Spearman only flips sign

Usage
-----
    python scripts/verify_viral_matrices.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "src"))

from amino_acid_utils import (  # noqa: E402
    AA_3_TO_1,
    calculate_property_distances,
    generate_amino_acid_pairs,
    precompute_codon_distances,
)
from build_added_matrices import AA20, THIRD_BITS  # noqa: E402
from build_viral_matrices import DIVERGENCE, MODELS, SOURCES, parse_iqtree  # noqa: E402
from mutation_score import get_mutational_scores, get_substitution_matrix  # noqa: E402

TOLERANCE = 1e-9


def generator(rates: np.ndarray, frequencies: np.ndarray) -> np.ndarray:
    """The normalised instantaneous rate matrix, rebuilt here rather than imported."""
    q = rates * frequencies[None, :]
    np.fill_diagonal(q, 0.0)
    np.fill_diagonal(q, -q.sum(axis=1))
    return q / -np.sum(frequencies * np.diag(q))


def transition_by_eigendecomposition(q: np.ndarray, frequencies: np.ndarray,
                                     t: float) -> np.ndarray:
    """exp(Qt) through the symmetric eigenproblem, avoiding scipy's Pade approximation.

    A reversible generator is similar to a symmetric matrix under the square roots of the
    equilibrium frequencies, so the exponential can be taken on the eigenvalues.
    """
    root = np.sqrt(frequencies)
    symmetric = q * root[:, None] / root[None, :]
    values, vectors = np.linalg.eigh((symmetric + symmetric.T) / 2)
    exponentiated = vectors @ np.diag(np.exp(values * t)) @ vectors.T
    return exponentiated / root[:, None] * root[None, :]


def composite_from_json(name: str, properties: pd.DataFrame,
                        codon_distances: dict) -> float:
    """The composite recomputed from the stored entries, bypassing the pipeline scorer."""
    matrix = get_substitution_matrix(name)
    pairs = generate_amino_acid_pairs()
    entries = np.array([matrix[(pair[0], pair[1])] for pair in pairs])
    total = sum(
        stats.spearmanr(calculate_property_distances(properties[column].to_dict(), pairs),
                        entries).statistic
        for column in properties.columns
    )
    total += stats.spearmanr(np.array([codon_distances[pair] for pair in pairs]),
                             entries).statistic
    return -float(total)


def composite_from_pipeline(name: str, properties: pd.DataFrame,
                            codon_distances: dict) -> float:
    """The composite as the evaluation computes it, through the mutational scores."""
    pairs = generate_amino_acid_pairs()
    scores = get_mutational_scores(pairs, substitution_matrix_type=name)
    values = np.array([scores[pair] for pair in pairs])
    total = sum(
        stats.spearmanr(calculate_property_distances(properties[column].to_dict(), pairs),
                        values).statistic
        for column in properties.columns
    )
    return float(total + stats.spearmanr(
        np.array([codon_distances[pair] for pair in pairs]), values).statistic)


def main() -> None:
    table = pd.read_csv(HERE / "amino_acid_properties.csv")
    table["AA_1letter"] = table["AA"].map(AA_3_TO_1)
    properties = table.set_index("AA_1letter").drop(columns=["AA"])
    codon_distances = precompute_codon_distances()
    source = (SOURCES / "iqtree_modelprotein.cpp").read_text(encoding="utf-8", errors="replace")

    failures = []
    print(f"{'model':<7}{'generator':>11}{'reversible':>12}{'eigen vs expm':>15}"
          f"{'entries':>9}{'composite':>11}{'two routes':>12}")

    for iqtree_name, published in MODELS.items():
        rates, frequencies = parse_iqtree(source, iqtree_name)
        q = generator(rates, frequencies)

        rows_zero = float(np.abs(q.sum(axis=1)).max())
        normalised = float(-np.sum(frequencies * np.diag(q)))
        p = expm(q * DIVERGENCE)
        stochastic = float(np.abs(p.sum(axis=1) - 1).max())
        non_negative = bool((p >= -TOLERANCE).all())
        balance = float(np.abs(frequencies[:, None] * p - (frequencies[:, None] * p).T).max())
        p_eigen = transition_by_eigendecomposition(q, frequencies, DIVERGENCE)
        eigen_gap = float(np.abs(p - p_eigen).max())

        # The committed file must be what this reconstruction produces. Read the JSON itself:
        # the loader min-max normalises what it returns, so it cannot be compared to raw
        # third-bit integers.
        rebuilt = np.round(np.log(p / frequencies[None, :]) * THIRD_BITS)
        stored = json.loads(
            (HERE.parent / "src" / "substitution_matrices" / f"{published}.json")
            .read_text(encoding="utf-8"))
        entries_match = all(
            int(rebuilt[i, j]) == int(stored[f"({AA20[i]}, {AA20[j]})"])
            for i in range(20) for j in range(20)
        )

        from_pipeline = composite_from_pipeline(published, properties, codon_distances)
        from_json = composite_from_json(published, properties, codon_distances)
        routes_agree = abs(from_pipeline - from_json) < 5e-4

        ok = (rows_zero < TOLERANCE and abs(normalised - 1) < TOLERANCE
              and stochastic < 1e-12 and non_negative and balance < 1e-12
              and eigen_gap < 1e-10 and entries_match and routes_agree)
        if not ok:
            failures.append(published)

        print(f"{published:<7}{rows_zero:>11.1e}{balance:>12.1e}{eigen_gap:>15.1e}"
              f"{'exact' if entries_match else 'DIFFER':>9}{from_pipeline:>11.3f}"
              f"{'agree' if routes_agree else 'DISAGREE':>12}")

    print()
    if failures:
        raise SystemExit(f"checks failed for: {', '.join(failures)}")
    print(f"All {len(MODELS)} pass: valid normalised generator, reversible, exponential "
          "confirmed by a second algorithm,\nstored entries reproduced, composite identical "
          "by two routes.")


if __name__ == "__main__":
    main()
