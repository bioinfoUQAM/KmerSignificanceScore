"""
What the prevalence threshold changes (Reviewer #1, comment 2.2).

    "A robust bioinformatic tool cannot rely on manually tweaked thresholds without
    justification. The authors must include supplementary material with a sensitivity analysis
    varying the parameter t (e.g., from 0.1 to 0.5) to demonstrate the impact of this threshold
    on score variance and method stability."

The objection is specifically that the default of 0.25 was raised to 0.33 for HIV-1, so HIV-1
is the dataset that has to answer it. This script reads the sweep directories produced by
`threshold_sweep_fast.py` and `threshold_sweep_sars.py`, or by `threshold_sweep.py`, which
covers all three datasets in one run, and reports two things the comment asks for:

* **score variance**, the spread of the KSS distribution at each threshold;
* **method stability**, the overlap between the positions a given threshold selects and those
  selected at the value used in the manuscript.

Stability is the one that matters. A threshold that changes the scores but not the ranking has
no consequence for a prioritization, and that is the claim to test rather than assert.

Usage
-----
    python scripts/threshold_sensitivity.py --sweep-dir <directory> [--sweep-dir <another>]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from kss_selection import ordered_prioritization

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / "results" / "threshold_sensitivity"

# The value each dataset uses in the manuscript, which is the reference the sweep is read
# against. HIV-1 is the exception the reviewer objects to.
PUBLISHED = {
    "Severe_acute_respiratory_syndrome_coronavirus_2": 25,
    "Human_immunodeficiency_virus_1": 33,
    "Human_betaherpesvirus_5": 25,
}
LABELS = {
    "Severe_acute_respiratory_syndrome_coronavirus_2": "SARS-CoV-2",
    "Human_immunodeficiency_virus_1": "HIV-1",
    "Human_betaherpesvirus_5": "HCMV",
}
TOP_K = [12, 25]

# The selection each dataset publishes, which is what the stability of the prioritization has
# to be measured on. SARS-CoV-2 and HIV-1 are reported as twelve positions pooled across genes;
# HCMV is reported gene by gene, ten positions each, because its genotyping schemes are
# locus-specific. Measuring HCMV as a pooled twelve, as this script did, measured a list the
# manuscript never shows.
PUBLISHED_SELECTION = {
    "Severe_acute_respiratory_syndrome_coronavirus_2": ("pooled", 12),
    "Human_immunodeficiency_virus_1": ("pooled", 12),
    "Human_betaherpesvirus_5": ("per_gene", 10),
}


def read_gene(path: Path) -> pd.DataFrame:
    payload = json.loads(path.read_text(encoding="utf-8"))
    gene = path.stem.replace("_compiled_results", "")
    entries = payload.get(gene) or next(iter(payload.values()))
    rows = []
    for position, record in entries.items():
        rows.append({
            "position": int(position),
            "kss": float(record.get("kss", 0.0)),
            "discriminative": float(record.get("discriminative_score", 0.0)),
            "n_variants": len(record.get("alts", {}) or {}),
        })
    return pd.DataFrame(rows)


def load_sweep(directories: list[Path]) -> dict:
    """{(dataset, threshold): {gene: frame}} over every sweep directory supplied."""
    loaded: dict = {}
    for root in directories:
        for tdir in sorted(root.glob("t*")):
            if not (tdir / "COMPLETE").exists():
                continue
            threshold = int(tdir.name[1:])
            for dataset_dir in sorted(tdir.iterdir()):
                if not dataset_dir.is_dir():
                    continue
                genes = {}
                for compiled in sorted(dataset_dir.glob("*/results/*_compiled_results.json")):
                    genes[compiled.stem.replace("_compiled_results", "")] = read_gene(compiled)
                if genes:
                    loaded[(dataset_dir.name, threshold)] = genes
    return loaded


def top_positions(genes: dict, k: int) -> set:
    pooled = pd.concat([frame.assign(gene=gene) for gene, frame in genes.items()],
                       ignore_index=True)
    pooled = pooled.rename(columns={"position": "nt_position", "kss": "KSS"})
    ordered = ordered_prioritization(
        pooled, score="KSS", gene="gene", position="nt_position")
    return set(zip(ordered.head(k)["gene"], ordered.head(k)["nt_position"]))


def per_gene_positions(genes: dict, k: int) -> set:
    """The k highest-scoring positions of each gene, as (gene, position) pairs."""
    chosen = set()
    for gene, frame in genes.items():
        table = frame.rename(columns={"position": "nt_position", "kss": "KSS"}).assign(gene=gene)
        ordered = ordered_prioritization(
            table, score="KSS", gene="gene", position="nt_position")
        chosen |= set(zip(ordered.head(k)["gene"], ordered.head(k)["nt_position"]))
    return chosen


def published_selection(genes: dict, dataset: str) -> set:
    """The positions the manuscript reports for this dataset, under its own selection rule."""
    rule, size = PUBLISHED_SELECTION[dataset]
    return (per_gene_positions(genes, size) if rule == "per_gene"
            else top_positions(genes, size))


def kss_vector(genes: dict) -> pd.Series:
    """KSS score of every position, indexed by (gene, position), for rank correlation."""
    parts = []
    for gene, frame in genes.items():
        series = frame.set_index("position")["kss"]
        series.index = [(gene, position) for position in series.index]
        parts.append(series)
    return pd.concat(parts)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep-dir", type=Path, action="append", required=True)
    args = parser.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)

    sweep = load_sweep(args.sweep_dir)
    if not sweep:
        print("no completed threshold directories found")
        return 1
    datasets = sorted({dataset for dataset, _ in sweep})
    print(f"loaded {len(sweep)} dataset-threshold combinations "
          f"over {len(datasets)} datasets\n")

    rows = []
    for dataset in datasets:
        thresholds = sorted(t for d, t in sweep if d == dataset)
        reference = PUBLISHED.get(dataset)
        if reference not in thresholds:
            print(f"  {dataset}: published threshold {reference} not in the sweep, skipping")
            continue
        baseline = {k: top_positions(sweep[(dataset, reference)], k) for k in TOP_K}
        baseline_published = published_selection(sweep[(dataset, reference)], dataset)
        baseline_scores = kss_vector(sweep[(dataset, reference)])

        for threshold in thresholds:
            genes = sweep[(dataset, threshold)]
            pooled = pd.concat(genes.values(), ignore_index=True)
            # Full-vector rank stability: Spearman of the whole KSS vector against the
            # published threshold, over the positions common to both. Unlike top-k overlap
            # this asks whether the entire ranking, not just its head, survives the threshold.
            scores = kss_vector(genes)
            common = baseline_scores.index.intersection(scores.index)
            spearman_full = (spearmanr(scores.loc[common], baseline_scores.loc[common]).statistic
                             if len(common) > 2 else float("nan"))
            row = {
                "dataset": LABELS.get(dataset, dataset),
                "threshold": threshold,
                "is_published_value": threshold == reference,
                "positions": len(pooled),
                "variants": int(pooled["n_variants"].sum()),
                "kss_mean": pooled["kss"].mean(),
                "kss_sd": pooled["kss"].std(),
                "kss_variance": pooled["kss"].var(),
                "discriminative_sd": pooled["discriminative"].std(),
                "spearman_full_vs_published": spearman_full,
            }
            rule, size = PUBLISHED_SELECTION[dataset]
            selection = published_selection(genes, dataset)
            row["selection_rule"] = ("ten per gene" if rule == "per_gene"
                                     else f"{size} pooled across genes")
            row["selection_size"] = len(baseline_published)
            row["shared_selection"] = len(selection & baseline_published)
            row["jaccard_selection"] = (len(selection & baseline_published)
                                        / len(selection | baseline_published))
            for k in TOP_K:
                pooled_selection = top_positions(genes, k)
                shared = len(pooled_selection & baseline[k])
                row[f"shared_top{k}"] = shared
                row[f"jaccard_top{k}"] = shared / len(pooled_selection | baseline[k])
            rows.append(row)

    table = pd.DataFrame(rows)
    table.to_csv(RESULTS / "threshold_sensitivity.csv", index=False)

    for dataset, group in table.groupby("dataset"):
        print("=" * 78)
        print(dataset)
        print("=" * 78)
        view = group.sort_values("threshold")
        print(f"{'t':>5} {'positions':>10} {'variants':>9} {'KSS sd':>8} "
              f"{'top12 kept':>11} {'Jaccard12':>10} {'Spearman':>9}")
        for _, r in view.iterrows():
            mark = "  <- manuscript" if r["is_published_value"] else ""
            print(f"{r['threshold'] / 100:>5.2f} {r['positions']:>10,} {r['variants']:>9,} "
                  f"{r['kss_sd']:>8.4f} {int(r['shared_top12']):>8}/12 "
                  f"{r['jaccard_top12']:>10.3f} {r['spearman_full_vs_published']:>9.4f}{mark}")
        worst = view["jaccard_top12"].min()
        print(f"\n  lowest top-12 agreement with the published threshold: {worst:.3f}")
        print(f"  KSS standard deviation ranges from {view['kss_sd'].min():.4f} "
              f"to {view['kss_sd'].max():.4f}\n")

    hiv = table[table["dataset"] == "HIV-1"]
    if len(hiv) and 25 in set(hiv["threshold"]):
        at25 = hiv[hiv["threshold"] == 25].iloc[0]
        at33 = hiv[hiv["threshold"] == 33].iloc[0]
        print("=" * 78)
        print("The specific objection: HIV-1 at 0.25 against the 0.33 that was used")
        print("=" * 78)
        print(f"  positions   {int(at25['positions']):>6} -> {int(at33['positions']):>6}")
        print(f"  variants    {int(at25['variants']):>6} -> {int(at33['variants']):>6}")
        print(f"  KSS sd      {at25['kss_sd']:>6.4f} -> {at33['kss_sd']:>6.4f}")
        print(f"  top-12 positions shared between the two: "
              f"{int(at25['shared_top12'])}/12")

    print(f"\nwritten to {RESULTS.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

