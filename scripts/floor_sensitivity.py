"""What the 0.1 floor of Eq. (9) does, across the matrix panel and in the reported analyses.

Reviewer #1, comment 3.1, asks why the mutational impact is floored at 0.1. The manuscript
answers with three facts, and this script is where each one comes from.

1. The floor is necessary because the normalisation of Eq. (8) maps the largest entry of the
   matrix to 1. In a distance matrix a pair at zero distance ties the diagonal and reaches that
   maximum, so its impact would be exactly zero. Whether any pair does depends on the matrix,
   so the panel is scanned rather than argued about.
2. It binds for 8 of the 190 amino acid pairs of MIYATA_EVO, all biochemically interchangeable.
3. It is not cosmetic and not decisive either, which is the honest way to answer the question:
   it settles the mutational term at 95 of the 4,197 scored positions, and it changes one
   position of one reported selection, which the script reports rather than assumes.

Reads the substitution matrices through the pipeline's own loader and the committed compiled
results of the three datasets. Writes results/floor_sensitivity/.

Usage
-----
    python scripts/floor_sensitivity.py
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

sys.path.insert(0, str(HERE))

from src.mutation_score import get_substitution_matrix  # noqa: E402
from kss_selection import figure3_display_selection  # noqa: E402

RESULTS = ROOT / "results" / "floor_sensitivity"
MATRIX = "MIYATA_EVO"
FLOOR = 0.1
CHANGE = re.compile(r"([A-Z])(\d+)([A-Z])")


STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


def scaled_pairs(name: str) -> dict[tuple[str, str], float]:
    """Off-diagonal pairs of a matrix, on the [0, 1] similarity scale the pipeline uses.

    Restricted to the twenty standard residues: BLASTP also carries J, O and U, and counting
    those would report 253 pairs for it and a smallest impact four times closer to the floor
    than the one the pipeline can ever meet.
    """
    matrix = get_substitution_matrix(name)
    return {tuple(sorted(pair)): value for pair, value in matrix.items()
            if pair[0] != pair[1] and set(pair) <= STANDARD}


def panel() -> pd.DataFrame:
    """How many pairs the floor rescues in each matrix of the evaluated panel."""
    rows = []
    for path in sorted((ROOT / "src" / "substitution_matrices").glob("*.json")):
        pairs = scaled_pairs(path.stem)
        impacts = {pair: 1 - value for pair, value in pairs.items()}
        above = [value for value in impacts.values() if value >= FLOOR]
        rows.append({
            "matrix": path.stem,
            "pairs": len(pairs),
            "binding": sum(1 for value in impacts.values() if value < FLOOR),
            "zero_impact": sum(1 for value in impacts.values() if value <= 0),
            "smallest_impact": min(impacts.values()),
            # What the value of the floor buys: everything below it is flattened to it, so the
            # first impact above says how much room there was, and doubling the floor says how
            # many distinctions a larger one would have erased.
            "smallest_impact_above_floor": min(above) if above else float("nan"),
            "binding_at_double_floor": sum(1 for value in impacts.values() if value < 2 * FLOOR),
        })
    return pd.DataFrame(rows).sort_values(["binding", "matrix"], ascending=[False, True])


def rescued_pairs() -> pd.DataFrame:
    """The pairs of the published matrix the floor rescues, with what it gives them."""
    impacts = {pair: 1 - value for pair, value in scaled_pairs(MATRIX).items()}
    rows = [{"pair": f"{a}-{b}", "scaled_similarity": 1 - value, "impact_without_floor": value,
             "impact_with_floor": max(value, FLOOR), "correction": FLOOR - value}
            for (a, b), value in sorted(impacts.items(), key=lambda kv: kv[1])
            if value < FLOOR]
    return pd.DataFrame(rows)


def natural_impact(record: dict, matrix: dict) -> float:
    """The mutational term the position would carry with no floor.

    The stored per-change values are already floored, so a substitution is recomputed from the
    matrix. Anything the pipeline scored without an amino acid pair, an insertion or a
    deletion, keeps the value it was given: the floor never touches those.
    """
    natural = 0.0
    for variant in record.get("alts", {}).values():
        for change, value in (variant.get("amino_acid_changes") or {}).items():
            match = CHANGE.fullmatch(str(change))
            if match and match.group(1) != match.group(3):
                natural = max(natural, 1 - matrix[(match.group(1), match.group(3))])
            else:
                natural = max(natural, float(value))
    return 1.0 if record.get("mutational_score") == 1.0 else natural


def published_selections() -> dict[str, list[tuple[str, pd.DataFrame]]]:
    """The selections the manuscript reports, built by the module that defines them.

    Not a per-gene top ten: SARS-CoV-2 and HIV-1 are reported as one pooled selection of twelve
    over their genes, HCMV as ten per gene, and all three apply the three-indel selection rule. Ranking anything else here would test an object the manuscript never shows,
    which is the confusion `kss_selection` exists to prevent.
    """
    matrix = get_substitution_matrix(MATRIX)
    frames: dict[str, list[dict]] = {}
    for path in sorted(ROOT.glob("data/*/*/results/*_compiled_results.json")):
        if "toy" in path.parts:
            continue
        dataset, gene = path.parts[-4], path.parts[-3]
        payload = json.loads(path.read_text(encoding="utf-8"))
        records = payload[next(iter(payload))] if len(payload) == 1 else payload
        for position, record in records.items():
            frames.setdefault(dataset, []).append({
                "gene": gene,
                "nt_position": int(position),
                "discriminative_score": record["discriminative_score"],
                "protein_score": record["protein_score"],
                "mutational_score": record["mutational_score"],
                "natural": natural_impact(record, matrix),
                "kss": record["kss"],
            })

    selections = {}
    for dataset, rows in frames.items():
        frame = pd.DataFrame(rows)
        if dataset.startswith("Human_betaherpesvirus"):
            selections[dataset] = [(gene, block) for gene, block in frame.groupby("gene")]
        else:
            selections[dataset] = [("pooled", frame)]
    return selections


def selection_at(frame: pd.DataFrame, floor: float, size: int) -> list[tuple[str, int]]:
    """The selection the manuscript would report if the floor took this value."""
    scored = frame.assign(
        mutational_score=frame["natural"].clip(lower=floor),
        KSS=((frame["discriminative_score"] + frame["natural"].clip(lower=floor)
              + frame["protein_score"]) / 3).round(3))
    drawn = figure3_display_selection(scored, size, score="KSS")
    return [(row.gene, int(row.nt_position)) for row in drawn.itertuples()]


def positions() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """How often the floor decides the term, what it moves, and how far its value can travel."""
    counts, moves, stability = [], [], []
    for dataset, blocks in published_selections().items():
        for label, frame in blocks:
            size = 10 if dataset.startswith("Human_betaherpesvirus") else 12
            counts.append({
                "dataset": dataset, "gene": label, "positions": len(frame),
                "silent": int((frame["mutational_score"] == 0).sum()),
                "with_amino_acid_change": int((frame["mutational_score"] > 0).sum()),
                "decided_by_floor": int((frame["mutational_score"] == FLOOR).sum()),
            })

            published = selection_at(frame, FLOOR, size)
            without = selection_at(frame, 0.0, size)
            for entry in set(published) ^ set(without):
                gene, position = entry
                row = frame[(frame["gene"] == gene)
                            & (frame["nt_position"] == position)].iloc[0]
                moves.append({
                    "dataset": dataset, "gene": gene, "position": position,
                    "kss": row["kss"],
                    "kss_without_floor": round((row["discriminative_score"] + row["natural"]
                                                + row["protein_score"]) / 3, 3),
                    "in_published": entry in published,
                    "rank_published": published.index(entry) + 1 if entry in published else None,
                    "rank_without_floor": without.index(entry) + 1 if entry in without else None,
                })

            # How far the value itself can travel before the reported selection moves. The pair
            # plateau is a property of the matrix; this is the property of the published result,
            # and it is what the question about the value is really asking.
            same_set, same_order = [], []
            for step in range(0, 301):
                floor = step / 1000
                candidate = selection_at(frame, floor, size)
                if set(candidate) == set(published):
                    same_set.append(floor)
                if candidate == published:
                    same_order.append(floor)
            stability.append({
                "dataset": dataset, "gene": label, "selection_size": size,
                "membership_from": min(same_set), "membership_to": max(same_set),
                "order_from": min(same_order), "order_to": max(same_order),
            })
    return (pd.DataFrame(counts).sort_values(["dataset", "gene"]),
            pd.DataFrame(moves), pd.DataFrame(stability).sort_values(["dataset", "gene"]))


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)

    matrices = panel()
    rescued = rescued_pairs()
    counts, moves, stability = positions()

    matrices.to_csv(RESULTS / "floor_by_matrix.csv", index=False)
    rescued.to_csv(RESULTS / "rescued_pairs.csv", index=False)
    counts.to_csv(RESULTS / "floor_by_gene.csv", index=False)
    moves.to_csv(RESULTS / "selection_changes.csv", index=False)
    stability.to_csv(RESULTS / "floor_stability.csv", index=False)

    binding = matrices[matrices["binding"] > 0]
    print(f"{len(matrices)} matrices in src/substitution_matrices, which is the {len(matrices) - 1} "
          f"candidates of the panel plus {MATRIX} itself. The floor binds in {len(binding)}:")
    for row in binding.itertuples():
        print(f"   {row.matrix:<12} {row.binding:>3} of {row.pairs} pairs, "
              f"smallest impact {row.smallest_impact:.6f}")
    print(f"\n{MATRIX}: the floor rescues {len(rescued)} of 190 pairs "
          f"({100 * len(rescued) / 190:.1f}%), largest correction "
          f"{rescued['correction'].max():.3f}, which moves KSS by at most "
          f"{rescued['correction'].max() / 3:.4f} at equal weights")
    print("   " + ", ".join(rescued["pair"]))

    total = int(counts["positions"].sum())
    decided = int(counts["decided_by_floor"].sum())
    print(f"\n{decided} of {total} scored positions are decided by the floor "
          f"({100 * decided / total:.2f}%)")
    if moves.empty:
        print("no reported selection changes when the floor is removed")
    else:
        print("positions whose membership of a reported selection depends on the floor:")
        for row in moves.itertuples():
            print(f"   {row.gene} position {row.position}: KSS {row.kss} -> "
                  f"{row.kss_without_floor}, "
                  f"{'in' if row.in_published else 'out'} with the floor")
    print("\nhow far the value can travel before a reported selection moves:")
    print(f"   membership held for every floor from {stability['membership_from'].max():.3f} "
          f"to {stability['membership_to'].min():.3f}")
    print(f"   order held for every floor from {stability['order_from'].max():.3f} "
          f"to {stability['order_to'].min():.3f}")
    print(f"\nwritten to {RESULTS.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
