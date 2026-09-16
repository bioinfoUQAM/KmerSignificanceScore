"""
Regenerate the substitution matrix ranking, from the matrices currently in the package.

The ranking was produced by `matrix_evaluation.ipynb` and committed as a CSV. Adding matrices
means the committed file no longer describes the package, and rerunning a notebook to fix that
is neither reproducible nor checkable. This script computes the same table, and proves it by
reproducing the committed rows before it is allowed to extend them.

The composite is the sum of eight Spearman correlations between amino acid property distances
and the mutational impact the pipeline derives from the matrix, exactly as `evaluate_matrix` does.

The correlations are reported without significance testing. The 190 pairs are built from 20
amino acids, each of which enters 19 of them, so they are not the independent observations
a Spearman p-value assumes; and the composite is the objective the genetic algorithm
maximizes, so testing the optimized matrix against it would be testing a selection against
the criterion that produced it.

Usage
-----
    python scripts/build_matrix_ranking.py            # verify, then write
    python scripts/build_matrix_ranking.py --check    # verify only, write nothing
"""

from __future__ import annotations

import argparse
import sys
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
from mutation_score import get_mutational_scores  # noqa: E402

MATRICES = HERE.parent / "src" / "substitution_matrices"
RESULTS = HERE.parent / "results" / "matrix_evaluation"
COMPLETE = RESULTS / "substitution_matrices_ranking_complete.csv"
TOP_10 = RESULTS / "substitution_matrices_ranking_top_10.csv"
IMPROVEMENT = RESULTS / "optimization_improvement.csv"
GAINS = RESULTS / "property_gains.csv"

OPTIMIZED = "MIYATA_EVO"          # the product of the optimization, not a candidate


def evaluate(name: str, properties: pd.DataFrame, codon_distances: dict) -> dict:
    """The eight correlations and the composite, for one matrix."""
    pairs = generate_amino_acid_pairs()
    impact = np.array([get_mutational_scores(pairs, substitution_matrix_type=name)[p]
                       for p in pairs])

    correlations = {}
    for column in properties.columns:
        distances = calculate_property_distances(properties[column].to_dict(), pairs)
        correlations[column] = float(stats.spearmanr(distances, impact).statistic)

    codon = stats.spearmanr(np.array([codon_distances[p] for p in pairs]), impact)
    correlations["MinCodD"] = float(codon.statistic)

    row = {"Matrix": name}
    for label, value in correlations.items():
        row[label] = f"{value:.3f}"
        # Unrounded alongside, for the reason given below for the composite. The per-property
        # gains the manuscript quotes are ratios of two correlations, and the composite was
        # the only ratio this precaution covered: two of the eight gains had been computed
        # from the three decimals of Table 5 and were wrong by half a point.
        row[f"{label}_exact"] = value
    row["Composite"] = round(sum(correlations.values()), 3)
    # Kept unrounded alongside, because the optimization improvement is a ratio of two
    # composites and computing it from displayed values would round twice.
    row["Composite_exact"] = sum(correlations.values())
    return row


def build(properties: pd.DataFrame, codon_distances: dict) -> pd.DataFrame:
    names = sorted(path.stem for path in MATRICES.glob("*.json"))
    table = pd.DataFrame([evaluate(name, properties, codon_distances) for name in names])
    table = table.sort_values("Composite", ascending=False).reset_index(drop=True)
    table.insert(0, "Rank", range(1, len(table) + 1))
    return table


def reproduces_committed(fresh: pd.DataFrame) -> bool:
    """Every matrix already in the committed ranking must get the same numbers."""
    if not COMPLETE.exists():
        print("no committed ranking to check against")
        return True

    old = pd.read_csv(COMPLETE).set_index("Matrix")
    new = fresh.set_index("Matrix")
    columns = [c for c in old.columns if c != "Rank"]

    disagreements = []
    for name in old.index:
        if name not in new.index:
            disagreements.append(f"{name} is in the committed ranking but not in the package")
            continue
        for column in columns:
            # The committed file carried a significance suffix that this script no
            # longer writes; the number is what has to match, so it is compared.
            #
            # As numbers, not as text. This script writes three decimals, pandas reads
            # 0.130 back as 0.13, and comparing their string forms made every correlation
            # ending in a zero disagree with itself. The guard then refused to write, so the
            # ranking had stopped being reproducible by the script that produces it, on a
            # difference of formatting alone.
            committed = str(old.loc[name, column]).rstrip("*")
            recomputed = str(new.loc[name, column])
            try:
                agrees = float(committed) == float(recomputed)
            except ValueError:
                agrees = committed == recomputed
            if not agrees:
                disagreements.append(
                    f"{name} {column}: committed {old.loc[name, column]}, "
                    f"recomputed {new.loc[name, column]}")

    print(f"{len(old)} committed rows checked against the current matrices: "
          f"{'all reproduced' if not disagreements else str(len(disagreements)) + ' DIFFER'}")
    for line in disagreements[:10]:
        print(f"    {line}")
    return not disagreements


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify without writing")
    args = parser.parse_args()

    table = pd.read_csv(HERE / "amino_acid_properties.csv")
    table["AA_1letter"] = table["AA"].map(AA_3_TO_1)
    properties = table.set_index("AA_1letter").drop(columns=["AA"])

    fresh = build(properties, precompute_codon_distances())
    if not reproduces_committed(fresh):
        raise SystemExit("the recomputation does not reproduce the committed ranking; "
                         "refusing to overwrite it")

    candidates = fresh[fresh["Matrix"] != OPTIMIZED].copy()
    candidates["Rank among candidates"] = range(1, len(candidates) + 1)
    print(f"\n{len(fresh)} matrices ranked, {len(candidates)} of them candidates "
          f"({OPTIMIZED} reported separately)")
    print(f"top ten overall:            {', '.join(fresh['Matrix'].head(10))}")
    print(f"top ten among candidates:   {', '.join(candidates['Matrix'].head(10))}")

    if args.check:
        print("\n--check given, nothing written")
        return

    # The improvement the abstract quotes used to rest on optimization_summary.csv, which only
    # a notebook writes and which no clone can rebuild. It is a ratio of two composites this
    # script already computes, so it is derived here from the committed matrices instead.
    composites = fresh.set_index("Matrix")["Composite_exact"]
    baseline, optimized = float(composites["MIYATA"]), float(composites[OPTIMIZED])
    improvement = pd.DataFrame([{
        "Baseline matrix": "MIYATA", "Optimized matrix": OPTIMIZED,
        "Baseline composite": baseline, "Optimized composite": optimized,
        "Absolute improvement": optimized - baseline,
        "Relative improvement (%)": (optimized - baseline) / baseline * 100,
    }])
    print(f"\n{OPTIMIZED} composite {optimized:.6f} against {baseline:.6f} for MIYATA, "
          f"an improvement of {improvement['Relative improvement (%)'].iloc[0]:.2f}%")

    # The per-property gains, from the same unrounded correlations. The manuscript quotes the
    # eight of them, and they have no other exact source: the ranking file rounds every
    # correlation to three decimals for display.
    exact = fresh.set_index("Matrix")
    gains = pd.DataFrame([{
        "Property": label,
        "MIYATA": float(exact.loc["MIYATA", f"{label}_exact"]),
        OPTIMIZED: float(exact.loc[OPTIMIZED, f"{label}_exact"]),
        "Relative gain (%)": (float(exact.loc[OPTIMIZED, f"{label}_exact"])
                              - float(exact.loc["MIYATA", f"{label}_exact"]))
                             / float(exact.loc["MIYATA", f"{label}_exact"]) * 100,
    } for label in properties.columns.tolist() + ["MinCodD"]])

    display = fresh.drop(columns=[c for c in fresh.columns if c.endswith("_exact")])
    display.to_csv(COMPLETE, index=False)
    display.head(10).to_csv(TOP_10, index=False)
    improvement.to_csv(IMPROVEMENT, index=False)
    gains.to_csv(GAINS, index=False)
    print(f"written to {COMPLETE.name}, {TOP_10.name}, {IMPROVEMENT.name} and {GAINS.name}")


if __name__ == "__main__":
    main()
