"""
What the conversion divergence costs, for every rate model converted into the benchmark.

Six of the evaluated matrices are published as rate models rather than score matrices: WAG,
and the five virus-derived models added for reviewer point 2.10. Putting any of them in a
benchmark of score matrices requires choosing an evolutionary distance at which to evaluate
P(delta) = exp(Q delta), in expected substitutions per site. We take delta = 2.5 as the
reference, comparable to the 250-PAM matrices already in the set. It is our choice and not a
value the sources publish, which is why the column below is named for a reference rather than
for a publication, and a reviewer is entitled to ask what that single free parameter buys.

The answer has to be measured on the matrices that are actually evaluated, which are integer
third-bit log-odds, not on the continuous log-odds they are rounded from. Rounding is not a
detail here: it collapses the 190 off-diagonal pairs onto a handful of distinct values, so rank
correlations computed before and after it are different quantities. This script reports both,
along with the composite and the resulting position in the ranking, at four divergences.

Usage
-----
    python scripts/divergence_sensitivity.py
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

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
from build_added_matrices import AA20, THIRD_BITS, load_wag_rate, wag_logodds  # noqa: E402
from build_viral_matrices import MODELS, SOURCES, parse_iqtree  # noqa: E402

OUTPUT = HERE.parent / "results" / "matrix_evaluation" / "divergence_sensitivity.csv"
RANKING = (HERE.parent / "results" / "matrix_evaluation"
           / "substitution_matrices_ranking_complete.csv")
DIVERGENCES = (2.0, 2.5, 3.0, 4.0)
REFERENCE = 2.5   # our conversion distance, not one the sources publish


def rate_models() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Every model in the benchmark that was published as a rate model."""
    source = (SOURCES / "iqtree_modelprotein.cpp").read_text(encoding="utf-8", errors="replace")
    models = {published: parse_iqtree(source, iqtree_name)
              for iqtree_name, published in MODELS.items()}
    models["WAG"] = load_wag_rate()
    return models


def score_matrix(rates: np.ndarray, frequencies: np.ndarray, t: float) -> np.ndarray:
    """The integer third-bit log-odds actually written into the benchmark."""
    return np.round(wag_logodds(rates, frequencies, t) * THIRD_BITS)


def composite(matrix: np.ndarray, properties: pd.DataFrame, codon_distances: dict) -> float:
    """The evaluation's composite, in its published orientation.

    The evaluation correlates property distances against mutational impact, which is a strictly
    decreasing transform of the matrix entry, so working from the entry gives the negative.
    """
    pairs = generate_amino_acid_pairs()
    entries = np.array([matrix[AA20.index(p[0]), AA20.index(p[1])] for p in pairs])
    total = sum(
        stats.spearmanr(calculate_property_distances(properties[column].to_dict(), pairs),
                        entries).statistic
        for column in properties.columns
    )
    total += stats.spearmanr(np.array([codon_distances[p] for p in pairs]), entries).statistic
    return -float(total)


def main() -> None:
    table = pd.read_csv(HERE / "amino_acid_properties.csv")
    table["AA_1letter"] = table["AA"].map(AA_3_TO_1)
    properties = table.set_index("AA_1letter").drop(columns=["AA"])
    codon_distances = precompute_codon_distances()

    models = rate_models()

    # The question is where a model would land if its own t had been chosen differently, so
    # every other matrix keeps its published value, the rate models included. Ranking against
    # a pool stripped of the six rate models would produce positions that cannot be compared
    # with the published ranking.
    fixed = pd.read_csv(RANKING)
    pool = pd.concat([
        fixed[~fixed["Matrix"].isin(models)][["Matrix", "Composite"]],
        pd.DataFrame([{"Matrix": name,
                       "Composite": composite(score_matrix(*models[name], REFERENCE),
                                              properties, codon_distances)}
                      for name in models]),
    ], ignore_index=True)

    off_diagonal = list(combinations(range(20), 2))
    rows = []
    for name, (rates, frequencies) in models.items():
        reference = score_matrix(rates, frequencies, REFERENCE)
        baseline = np.array([reference[i, j] for i, j in off_diagonal])
        others = pool[pool["Matrix"] != name]

        for d in DIVERGENCES:
            matrix = score_matrix(rates, frequencies, d)
            entries = np.array([matrix[i, j] for i, j in off_diagonal])
            value = composite(matrix, properties, codon_distances)

            better = int((others["Composite"] > value).sum())
            candidates = others[others["Matrix"] != "MIYATA_EVO"]
            rows.append({
                "model": name,
                "divergence_delta": d,
                "is_reference_value": d == REFERENCE,
                "composite": value,
                "rank_overall": better + 1,
                "rank_among_candidates": int((candidates["Composite"] > value).sum()) + 1,
                "distinct_values": len(set(entries.tolist())),
                "spearman_vs_reference_delta": float(stats.spearmanr(baseline, entries).statistic),
            })

    result = pd.DataFrame(rows)
    OUTPUT.parent.mkdir(exist_ok=True)
    result.to_csv(OUTPUT, index=False)

    for name, group in result.groupby("model", sort=False):
        spread = group["composite"].max() - group["composite"].min()
        ranks = sorted(set(group["rank_among_candidates"]))
        print(f"{name:<7} composite {group['composite'].min():.3f} to "
              f"{group['composite'].max():.3f} (spread {spread:.3f}), "
              f"rank among candidates {ranks}, "
              f"rank correlation vs delta={REFERENCE} at worst "
              f"{group['spearman_vs_reference_delta'].min():.3f}")
    print(f"\nwritten to {OUTPUT}")


if __name__ == "__main__":
    main()
