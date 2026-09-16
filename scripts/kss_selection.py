"""
The two selections this work uses, defined once so that they cannot drift apart again.

They are different objects and conflating them caused a real defect. `functional_benchmark.py`
measured error rates against an external database on a ranking the paper never reports, because
the rule Figure 3 applies lived only in that figure's legend and nowhere in the code.

pooled_raw_topk
    What the method actually prioritizes: every scored position from every gene, pooled and
    ordered by score. This is the object comment 2.11 asks about, since a reviewer asking for
    the false positive rate of a prioritization is asking about what the algorithm outputs, not
    about what we chose to draw. Filtering it before measuring would amount to choosing which
    positions to be judged on.

figure3_display_selection
    What Figure 3 draws: the same ordering, but with at most `max_indels` indel positions, so
    that the panel does not fill with indels clustered in one domain. This is a legibility rule,
    it is stated in the figure legend, and it must never be used to filter a benchmark.

An indel position is one whose mutational score is 1.0, which is what `indel_score` assigns in
`src/mutation_score.py`. That is the manuscript's own criterion, not a reconstruction.

Ties break on gene name then position, matching the order the published results were produced
with, so that a tie never depends on dictionary ordering.
"""

from __future__ import annotations

import pandas as pd

INDEL_MUTATIONAL_SCORE = 1.0
FIGURE3_MAX_INDELS = 3


def _ordered(kmers: pd.DataFrame, score: str, gene: str, position: str) -> pd.DataFrame:
    return kmers.sort_values([score, gene, position], ascending=[False, True, True])


def ordered_prioritization(kmers: pd.DataFrame, score: str = "KSS", gene: str = "gene",
                           position: str = "nt_position") -> pd.DataFrame:
    """Every scored position, pooled and ordered, for callers that cut at several k."""
    return _ordered(kmers, score, gene, position)


def pooled_raw_topk(kmers: pd.DataFrame, k: int, score: str = "KSS",
                    gene: str = "gene", position: str = "nt_position") -> pd.DataFrame:
    """The k highest-scoring positions, pooled across genes, with nothing filtered out."""
    return ordered_prioritization(kmers, score, gene, position).head(k)


def figure3_display_selection(kmers: pd.DataFrame, k: int, score: str = "KSS",
                              gene: str = "gene", position: str = "nt_position",
                              mutational: str = "mutational_score",
                              max_indels: int = FIGURE3_MAX_INDELS) -> pd.DataFrame:
    """The k positions Figure 3 draws, admitting at most `max_indels` indel positions.

    A display rule. Passing its output to a benchmark measures a curated subset and invites the
    objection that the test set was chosen after seeing the scores.
    """
    kept, indels = [], 0
    for index, row in _ordered(kmers, score, gene, position).iterrows():
        is_indel = float(row[mutational]) >= INDEL_MUTATIONAL_SCORE
        if is_indel and indels >= max_indels:
            continue
        kept.append(index)
        indels += is_indel
        if len(kept) == k:
            break
    return kmers.loc[kept]


def describe(selection: pd.DataFrame, score: str = "KSS", gene: str = "gene",
             mutational: str = "mutational_score") -> dict:
    """Composition, score range and indel count, the three things worth asserting on."""
    summary = {
        "n": len(selection),
        "genes": selection[gene].value_counts().to_dict(),
        "score_min": round(float(selection[score].min()), 2),
        "score_max": round(float(selection[score].max()), 2),
    }
    if mutational in selection.columns:
        summary["indels"] = int(
            (selection[mutational].astype(float) >= INDEL_MUTATIONAL_SCORE).sum())
    return summary
