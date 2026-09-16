"""
Overlap between the positions selected by the discriminative component and those selected
by the baseline feature selection methods.

Addresses Reviewer #1 comment 2.9: whether the methods that achieve comparable classification
performance do so by selecting the same positions.

Two complementary views are reported:

1. Jaccard index between top-k sets, method by method, for the same grid of k used in the
   classification benchmark. This measures agreement on the selection itself.
2. Spearman correlation between the full score vectors, which measures agreement on the whole
   ranking rather than only on its head.

Ties are broken exactly as in `discriminative_score_validation.ipynb`, that is by ascending
position, an explicit key rather than the row order of the score file. The scores themselves are computed
once on the full sequence set rather than per resampling iteration, so what is compared is the
criteria rather than one particular split.

Usage
-----
    python scripts/topk_overlap.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS_IN = ROOT / "notebooks" / "discriminative_score_validation_results"
RESULTS_OUT = ROOT / "results" / "statistical_comparison"

REFERENCE = "KSS_discriminative"
# The random baseline has no fixed score vector, so it has no reproducible top-k and no
# rank to correlate; every other baseline of the classification benchmark is here.
BASELINES = ["Chi2", "OddsRatio", "NMI", "MI", "SymChi2", "ANOVA", "CramerV"]

LABELS = {
    "KSS_discriminative": "Discriminative component",
    "Chi2": "Chi-squared",
    "OddsRatio": "Odds ratio",
    "NMI": "NMI",
    "MI": "Mutual information",
    "SymChi2": "Normalized symmetric chi-squared divergence",
    "ANOVA": "ANOVA",
    "CramerV": "Cramer's V",
}


def top_k_positions(scores: pd.DataFrame, method: str, k: int) -> set:
    """Top-k positions for one method, with the tie-breaking used by the benchmark.

    Ties are broken by ascending position, as the classification benchmark does. The
    key is explicit: a stable sort on the score alone would inherit the row order of
    the score file, which is not sorted by position, so equal scores at the top-k
    boundary would resolve differently there than in the benchmark.
    """
    ordered = scores.sort_values(
        [method, "position"], ascending=[False, True], kind="mergesort"
    )
    return set(ordered.head(k)["position"])


def jaccard(a: set, b: set) -> float:
    union = a | b
    return len(a & b) / len(union) if union else 1.0


def main() -> None:
    RESULTS_OUT.mkdir(exist_ok=True)
    scores = pd.read_csv(RESULTS_IN / "discriminative_scores_all_positions.csv")
    scores["dataset_gene"] = scores["dataset"] + " / " + scores["gene"]

    with open(RESULTS_IN / "config.json", encoding="utf-8") as handle:
        top_k_values = json.load(handle)["top_k_values"]

    overlap_rows, correlation_rows = [], []

    for unit, group in scores.groupby("dataset_gene"):
        n_positions = len(group)

        # Full-ranking agreement.
        for baseline in BASELINES:
            rho, _ = stats.spearmanr(group[REFERENCE], group[baseline])
            correlation_rows.append({
                "dataset_gene": unit,
                "baseline": LABELS[baseline],
                "n_positions": n_positions,
                "spearman_rho": float(rho),
            })

        for k in top_k_values:
            if k > n_positions:
                continue
            kss_top = top_k_positions(group, REFERENCE, k)

            for baseline in BASELINES:
                base_top = top_k_positions(group, baseline, k)
                overlap_rows.append({
                    "dataset_gene": unit,
                    "top_k": k,
                    "baseline": LABELS[baseline],
                    "shared_positions": len(kss_top & base_top),
                    "jaccard": jaccard(kss_top, base_top),
                })

    overlap = pd.DataFrame(overlap_rows)
    correlations = pd.DataFrame(correlation_rows)

    overlap.to_csv(RESULTS_OUT / "topk_jaccard_by_dataset.csv", index=False)
    correlations.to_csv(RESULTS_OUT / "score_rank_correlations.csv", index=False)

    # Summary averaged over dataset-gene units.
    summary = (
        overlap.pivot_table(index="top_k", columns="baseline", values="jaccard", aggfunc="mean")
        .round(3)
    )
    summary.to_csv(RESULTS_OUT / "topk_jaccard_summary.csv")

    print("Mean Jaccard index between the component top-k and each baseline top-k\n")
    print(summary.to_string())

    print("\n\nMean Spearman correlation of full score vectors\n")
    print(correlations.groupby("baseline")["spearman_rho"].mean().round(3).to_string())


if __name__ == "__main__":
    main()
