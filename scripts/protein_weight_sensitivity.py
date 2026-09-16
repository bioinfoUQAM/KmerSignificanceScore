"""
What the protein weight actually controls.

The protein component is a measure of how thoroughly a protein has been characterised, not
of how important it is: its ten categories record annotation, one of them counting publications.
That makes w_p a deliberate setting rather than a biological correction, and a reader is entitled
to know what moving it does to the output.

One structural fact drives everything here. The protein score is constant within a gene, so
w_p cannot reorder positions inside a gene. It only changes the balance *between* genes, and
therefore has no effect at all on a dataset ranked gene by gene. This script quantifies that
on the published results, without rerunning the pipeline: every compiled result already
stores the three components per position, and the score is their weighted mean.

    KSS = (w_d * S_d + w_m * S_m + w_p * S_p) / (w_d + w_m + w_p)

Usage
-----
    python scripts/protein_weight_sensitivity.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kss_selection import figure3_display_selection  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = HERE.parent / "data"
RESULTS = ROOT / "results" / "protein_weight_sensitivity"

WEIGHTS = [0.0, 0.5, 1.0, 2.0, 3.0]
REFERENCE_WEIGHT = 1.0
TOP_K = [12, 25]

DATASETS = {
    "SARS-CoV-2": ("Severe_acute_respiratory_syndrome_coronavirus_2",
                   ["ORF1ab", "S", "M", "N", "E"], "pooled"),
    "HIV-1": ("Human_immunodeficiency_virus_1", ["gag", "pol", "env"], "pooled"),
    "HCMV": ("Human_betaherpesvirus_5", ["UL55", "UL73", "US28"], "per gene"),
}


def load_positions(directory: str, gene: str) -> pd.DataFrame:
    path = DATA / directory / gene / "results" / f"{gene}_compiled_results.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    entries = payload[gene] if gene in payload else next(iter(payload.values()))
    rows = []
    for position, record in entries.items():
        rows.append({
            "gene": gene,
            "position": int(position),
            "discriminative": float(record.get("discriminative_score", 0.0)),
            "mutational": float(record.get("mutational_score", 0.0)),
            "protein": float(record.get("protein_score", 0.0)),
            "published_kss": float(record.get("kss", 0.0)),
        })
    return pd.DataFrame(rows)


def kss(df: pd.DataFrame, w_p: float, w_d: float = 1.0, w_m: float = 1.0) -> pd.Series:
    return (w_d * df["discriminative"] + w_m * df["mutational"] + w_p * df["protein"]) \
        / (w_d + w_m + w_p)


def top_positions(df: pd.DataFrame, w_p: float, k: int) -> set[tuple[str, int]]:
    scored = df.assign(score=kss(df, w_p))
    # Ties are broken by position so that the selection is reproducible.
    ordered = scored.sort_values(["score", "gene", "position"], ascending=[False, True, True])
    return set(zip(ordered.head(k)["gene"], ordered.head(k)["position"]))


def displayed_positions(df: pd.DataFrame, w_p: float, k: int) -> set[tuple[str, int]]:
    """The k positions the figure draws at this weight, cap on indel positions included.

    The raw ranking and the drawn selection are different objects, and only the second is what
    the manuscript displays. Reporting the first as "the published selection" is the confusion
    that comment 2.11 had to separate, so the sensitivity reports both.
    """
    scored = df.assign(score=kss(df, w_p)).rename(
        columns={"position": "nt_position", "mutational": "mutational_score"})
    drawn = figure3_display_selection(scored, k, score="score", gene="gene",
                                      position="nt_position", mutational="mutational_score")
    return set(zip(drawn["gene"], drawn["nt_position"]))


def main() -> int:
    RESULTS.mkdir(exist_ok=True)
    overlap_rows, composition_rows = [], []

    for label, (directory, genes, mode) in DATASETS.items():
        frames = [load_positions(directory, gene) for gene in genes]
        pooled = pd.concat(frames, ignore_index=True)

        # The published score must be recovered at the default weights, otherwise the
        # reconstruction is wrong and nothing below can be trusted.
        recomputed = kss(pooled, REFERENCE_WEIGHT).round(3)
        disagreements = int((recomputed - pooled["published_kss"]).abs().gt(0.0011).sum())
        print(f"{label}: {len(pooled)} positions across {len(genes)} genes, "
              f"{disagreements} disagreements with the published score at w_p = 1")

        protein_by_gene = pooled.groupby("gene")["protein"].nunique()
        constant = (protein_by_gene == 1).all()
        print(f"  protein score constant within every gene: {constant}")

        units = [("pooled", pooled)] if mode == "pooled" else [
            (gene, frame) for gene, frame in zip(genes, frames)]

        for unit_name, frame in units:
            for k in TOP_K:
                reference = top_positions(frame, REFERENCE_WEIGHT, k)
                reference_drawn = displayed_positions(frame, REFERENCE_WEIGHT, k)
                for w_p in WEIGHTS:
                    selection = top_positions(frame, w_p, k)
                    drawn = displayed_positions(frame, w_p, k)
                    shared = len(selection & reference)
                    union = len(selection | reference)
                    overlap_rows.append({
                        "dataset": label, "unit": unit_name, "selection_mode": mode,
                        "top_k": k, "w_p": w_p,
                        "shared_with_default": shared,
                        "jaccard_with_default": shared / union if union else 1.0,
                        "changed": k - shared,
                        # The same comparison on the selection the figure draws, cap included.
                        "drawn_shared_with_default": len(drawn & reference_drawn),
                        "drawn_changed": k - len(drawn & reference_drawn),
                    })
                    # Also measured against the lowest weight, so the span across the whole
                    # range is a computed number rather than something read off the gene
                    # counts. Distance from the default understates it: a position can leave
                    # on one side of w_p = 1 and a different one on the other.
                    lowest = top_positions(frame, min(WEIGHTS), k)
                    overlap_rows[-1]["changed_vs_lowest_weight"] = k - len(selection & lowest)
                    lowest_drawn = displayed_positions(frame, min(WEIGHTS), k)
                    overlap_rows[-1]["drawn_changed_vs_lowest_weight"] = (
                        k - len(drawn & lowest_drawn))
                    counts = pd.Series([gene for gene, _ in selection]).value_counts()
                    for gene, count in counts.items():
                        composition_rows.append({
                            "dataset": label, "unit": unit_name, "top_k": k,
                            "w_p": w_p, "gene": gene, "n_selected": int(count),
                            "gene_protein_score": float(
                                frame[frame["gene"] == gene]["protein"].iloc[0]),
                        })
        print()

    overlap = pd.DataFrame(overlap_rows)
    composition = pd.DataFrame(composition_rows)
    overlap.to_csv(RESULTS / "topk_overlap_by_weight.csv", index=False)
    composition.to_csv(RESULTS / "topk_composition_by_weight.csv", index=False)

    print("=" * 74)
    print("How many of the top 12 positions change when the protein weight moves")
    print("=" * 74)
    view = overlap[overlap["top_k"] == 12]
    table = view.pivot_table(index=["dataset", "unit"], columns="w_p", values="changed")
    print(table.to_string())

    print("\n" + "=" * 74)
    print("Where the top 12 positions sit, by gene (pooled datasets)")
    print("=" * 74)
    for label in ["SARS-CoV-2", "HIV-1"]:
        subset = composition[(composition["dataset"] == label)
                             & (composition["top_k"] == 12)]
        if subset.empty:
            continue
        print(f"\n{label}")
        grid = subset.pivot_table(index="gene", columns="w_p",
                                  values="n_selected", fill_value=0).astype(int)
        scores = subset.groupby("gene")["gene_protein_score"].first()
        grid.insert(0, "protein score", scores)
        print(grid.to_string())

    print(f"\nwritten to {RESULTS.relative_to(HERE.parent)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
