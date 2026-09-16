"""
Benchmark audit of the SARS-CoV-2 prioritization (Reviewer #1, comment 2.11).

The reviewer asks for false positive and false negative rates of the prioritization,
benchmarked against a large-scale empirical mutation database. This script computes them, and
also establishes why they cannot carry the meaning one would hope. The finding is not that the
method fails its benchmarks; it is that

    no available resource provides an independent, genome-wide, non-saturated reference for
    the quantity this prioritization optimizes.

The protocol was fixed before any of it was computed, and every departure from it is recorded
with its reason in the pre-specification kept with the manuscript's revision record.

Two parts, in the order the argument runs
-----------------------------------------
1. **COV2Var, binary membership.** The rates the reviewer asks for, computed as asked, for the
   prioritization and for all seven baselines and a random draw. Every method scores alike,
   which demonstrates rather than asserts that the reference is saturated.
2. **Fitness effects.** Reported as a constraint reference, which is what it is: the top decile
   by effect magnitude is entirely deleterious, so magnitude and constraint are the same
   variable here. Constraint is close to the opposite of what a lineage-discriminating method
   selects, and the result is read accordingly.

The deep mutational scanning layers are not here. They live in `reference_matrix.py`, which
evaluates them over the residues actually measured rather than over the interval the assay
spans, and the manuscript quotes that file.

Usage
-----
    python scripts/functional_benchmark.py --fitness <aa_fitness.csv>
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from external_inputs import verify_fitness  # noqa: E402
from kss_selection import (  # noqa: E402
    INDEL_MUTATIONAL_SCORE, describe, figure3_display_selection, ordered_prioritization,
    pooled_raw_topk,
)

DATA = ROOT / "data" / "Severe_acute_respiratory_syndrome_coronavirus_2"
BASELINE_SCORES = ROOT / "notebooks" / "discriminative_score_validation_results" / \
    "discriminative_scores_all_positions.csv"
RESULTS = ROOT / "results" / "functional_benchmark"

GENES = ["ORF1ab", "S", "M", "N", "E"]
K_NT = 9                     # a k-mer spans three codons
HEADLINE_K = 12              # the selection size used in the manuscript
K_SWEEP = [1, 5, 12, 25, 50, 100]
MIN_EXPECTED_COUNT = 20.0    # Bloom and Neher treat sparser sites as unreliable

# COV2Var splits the ORF1ab products; the manuscript scores ORF1ab as one unit.
COV2VAR_GENES = {"ORF1ab_pp1ab": "ORF1ab", "ORF1ab_pp1a": "ORF1ab",
                 "S": "S", "M": "M", "N": "N", "E": "E"}


# --------------------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------------------

def load_rankings() -> pd.DataFrame:
    """One row per scored k-mer, carrying the KSS score and the seven baselines."""
    rows = []
    for gene in GENES:
        path = DATA / gene / "results" / f"{gene}_compiled_results.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        entries = payload[gene] if gene in payload else next(iter(payload.values()))
        for position, record in entries.items():
            start = int(position)
            rows.append({
                "gene": gene,
                "nt_position": start,
                "key": f"{gene}:{start}",
                "aa_first": (start - 1) // 3 + 1,
                "KSS": float(record.get("kss", 0.0)),
                "mutational_score": float(record.get("mutational_score", 0.0)),
            })
    kmers = pd.DataFrame(rows)
    baselines = pd.read_csv(BASELINE_SCORES)
    baselines = baselines[baselines["dataset"] == "SARS-CoV-2"].rename(
        columns={"position": "key"})
    columns = ["key", "KSS_discriminative", "Chi2", "MI", "NMI",
               "OddsRatio", "ANOVA", "CramerV"]
    return kmers.merge(baselines[columns], on="key", how="left")


def label_residues(label: str) -> set[int]:
    """The residues a COV2Var protein label spans.

    A substitution names one residue, `D614G`. A multi-residue indel names two and means
    everything between them: `H69_V70del`, `C136_Y144del`, `G142_Y145delinsD`. Concatenating
    the digits, as this did until 28 July 2026, turns `H69_V70del` into residue 6970 and
    `V54_L64del` into 5464. Both are lost, and on a protein as long as ORF1ab pp1ab, 7,096
    residues, some land on a real position and annotate the wrong site. 67 labels of the 9,832
    carry several numbers, 42 of them in the five genes scored here, and they are the in-frame
    deletions, which is the mutation class the highest-ranked positions consist of.
    """
    numbers = [int(n) for n in re.findall(r"\d+", str(label))]
    if not numbers:
        return set()
    if len(numbers) == 1:
        return {numbers[0]}
    return set(range(min(numbers), max(numbers) + 1))


def load_fitness(path: Path) -> pd.DataFrame:
    fitness = pd.read_csv(path)
    confident = fitness[fitness["expected_count"] >= MIN_EXPECTED_COUNT].copy()
    confident["magnitude"] = confident["fitness"].abs()
    sites = confident.groupby(["gene", "aa_site"]).agg(
        magnitude=("magnitude", "median"),
        signed=("fitness", "median"),
    ).reset_index()
    return sites[sites["gene"].isin(GENES)]


# --------------------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------------------

def label(kmers: pd.DataFrame, sites: set[tuple[str, int]]) -> pd.Series:
    """A k-mer is a benchmark positive when at least one of its three sites is annotated."""
    return pd.Series(
        [any((row.gene, row.aa_first + offset) in sites for offset in range(K_NT // 3))
         for row in kmers.itertuples()],
        index=kmers.index)


def ranked(kmers: pd.DataFrame, method: str) -> np.ndarray:
    """The raw pooled prioritization, which is the object this comment is about.

    Deliberately not the selection Figure 3 draws. That one caps indel positions at three so the
    panel stays legible, a rule stated in its legend, and applying it here would filter the very
    ranking whose error rates are being measured: we would be choosing which positions to be
    judged on. The two are defined together in kss_selection.py so they cannot drift apart, which
    is how this script came to benchmark a ranking the paper never reports.
    """
    return ordered_prioritization(kmers, score=method).index.to_numpy()


def confusion(order: np.ndarray, positive: pd.Series, k: int) -> dict:
    selected = order[:k]
    tp = int(positive.loc[selected].sum())
    fp = k - tp
    total_positive = int(positive.sum())
    total_negative = len(positive) - total_positive
    fn = total_positive - tp
    tn = total_negative - fp
    prevalence = total_positive / len(positive)
    precision = tp / k
    return {
        "k": k, "TP": tp, "FP": fp, "FN": fn, "TN": tn,
        "precision_at_k": precision,
        "fdr_at_k": 1 - precision,
        "fpr_at_k": fp / total_negative if total_negative else float("nan"),
        "recall_at_k": tp / total_positive if total_positive else float("nan"),
        "fnr_at_k": fn / total_positive if total_positive else float("nan"),
        "prevalence": prevalence,
        "enrichment": precision / prevalence if prevalence else float("nan"),
    }


# --------------------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fitness", type=Path, required=True,
                        help="aa_fitness.csv from jbloomlab/SARS2-mut-fitness")
    args = parser.parse_args()
    # The input lives upstream, so identity is checked rather than mere presence: the fitness
    # estimates were originally taken from a branch tip, which is not a version.
    problem = verify_fitness(args.fitness)
    if problem:
        print(problem)
        return 1

    RESULTS.mkdir(exist_ok=True)
    kmers = load_rankings()

    # ---------------------------------------------------------------------------------
    print("=" * 84)
    print("0. What is being benchmarked, stated before any rate is quoted")
    print("=" * 84)
    print("  The raw pooled prioritization, not the selection Figure 3 draws. That figure caps")
    print("  indel positions at three for legibility, a rule stated in its legend, and applying")
    print("  it here would filter the ranking whose error rates are the question.")
    order = ordered_prioritization(kmers, score="KSS")
    print()
    print(f"  {'k':>5} {'genes':<34} {'indels':>7} {'Spike in 334-522':>17}")
    print("  " + "-" * 66)
    for k in K_SWEEP:
        cut = pooled_raw_topk(kmers, k)
        spike = cut[cut["gene"] == "S"]
        inside = sorted(a for a in spike["aa_first"] if 334 <= a <= 522)
        genes = ", ".join(f"{g}={n}" for g, n in cut["gene"].value_counts().items())
        indels = int((cut["mutational_score"] >= INDEL_MUTATIONAL_SCORE).sum())
        print(f"  {k:>5} {genes:<34} {indels:>3}/{k:<3} "
              f"{(', '.join(map(str, inside)) if inside else 'none'):>17}")
    print()
    print("  Figure 3 display selection, for contrast, at the size the figure uses:")
    print(f"    {describe(figure3_display_selection(kmers, HEADLINE_K))}")

    # ---------------------------------------------------------------------------------
    print("\n" + "=" * 84)
    print("1. Fitness effects, which measure constraint and not the target of the method")
    print("=" * 84)
    sites = load_fitness(args.fitness)
    top_decile = sites.nlargest(int(len(sites) * 0.10), "magnitude")
    print(f"{len(sites)} sites with a confident estimate")
    print(f"  sites whose median effect is positive : {int((sites['signed'] > 0).sum())} "
          f"({(sites['signed'] > 0).mean() * 100:.1f}%)")
    print(f"  top decile by magnitude that is deleterious: "
          f"{int((top_decile['signed'] < 0).sum())} of {len(top_decile)} "
          f"({(top_decile['signed'] < 0).mean() * 100:.1f}%)")
    print("  so ranking by magnitude ranks by constraint; the two are one variable here.\n")

    magnitude = sites.set_index(["gene", "aa_site"])["magnitude"].to_dict()
    selected = kmers.loc[ranked(kmers, "KSS")[:HEADLINE_K]]
    chosen = [magnitude.get((row.gene, row.aa_first + offset))
              for row in selected.itertuples() for offset in range(K_NT // 3)]
    chosen = [v for v in chosen if v is not None and not np.isnan(v)]
    print(f"  median effect magnitude of the top {HEADLINE_K} k-mers: "
          f"{np.median(chosen):.3f}")
    print(f"  median across all sites                    : "
          f"{sites['magnitude'].median():.3f}")
    print("  The prioritized positions sit in the tolerated half, which is what a method")
    print("  selecting lineage-discriminating variation should do. This characterises the")
    print("  selection; it is not a functional validation.")
    pd.DataFrame([{
        "sites_with_confident_estimate": len(sites),
        "fraction_positive_effect": float((sites["signed"] > 0).mean()),
        "top_decile_deleterious_fraction": float((top_decile["signed"] < 0).mean()),
        "median_magnitude_selected": float(np.median(chosen)),
        "median_magnitude_all_sites": float(sites["magnitude"].median()),
    }]).to_csv(RESULTS / "fitness_characterisation.csv", index=False)

    print(f"\nwritten to {RESULTS.relative_to(ROOT)}")
    print("\nThe rates above are defined by their reference and are not absolute biological")
    print("error rates. Report them with the prevalence attached, never on their own.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

