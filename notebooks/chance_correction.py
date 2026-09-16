"""
Is a discriminative score of 0.7 the same thing in two different viruses? (Comment 3.3)

The reviewer asks for an empirical demonstration that a given score carries the same absolute
discriminatory power across datasets, *controlling for chance*. The concern is real: with 19
imbalanced SARS-CoV-2 lineages, 15 HIV-1 subtypes and 4 to 8 HCMV genotypes, the score
reachable by chance alone is not the same in each setting. The complexity factor ln(K+1)
accounts explicitly for the number of classes, but the null depends on the whole set of margins
and not on that number alone.

How the null is built
---------------------
The score depends on the k-mer by class contingency table and on nothing else. That was
verified rather than assumed: recomputing it from the tables reproduces the published value at
all 4,197 positions with zero discrepancy. Permuting the class labels is therefore equivalent
to resampling a table with the same margins under independence, which is done here by drawing
each row from a multivariate hypergeometric distribution. No sequence has to be touched and no
classifier has to be retrained.

One subtlety matters. Sequences whose k-mer fell below the prevalence threshold appear in the
feature matrix as an all-zero row, so they form their own configuration: they contribute to the
conditional entropy but not to the purity term. That residual row is part of the table and is
resampled with the rest.

Usage
-----
    python notebooks/chance_correction.py [--permutations N]
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = HERE / "chance_correction_results"

SEED = 20260722

DATASETS = {
    "SARS-CoV-2": ("Severe_acute_respiratory_syndrome_coronavirus_2",
                   ["ORF1ab", "S", "M", "N", "E"], ["sars_cov2.tsv"]),
    "HIV-1": ("Human_immunodeficiency_virus_1", ["gag", "pol", "env"], ["hiv1.tsv"]),
    "HCMV": ("Human_betaherpesvirus_5", ["UL55", "UL73", "US28"],
             ["hcmv_UL55.tsv", "hcmv_UL73.tsv", "hcmv_US28.tsv"]),
}


def class_sizes(filename: str) -> Counter:
    with open(ROOT / "data" / "accessions" / filename, encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    column = [c for c in rows[0] if c != "accession"][0]
    return Counter(r[column] for r in rows)


def score_table(retained: np.ndarray, residual: np.ndarray) -> float:
    """Discriminative score of a contingency table, mirroring src/discriminative_score.py.

    `retained` holds one row per k-mer kept by the prevalence threshold; `residual` holds the
    per-class counts of sequences whose k-mer was dropped. The residual enters the conditional
    entropy as its own configuration but not the purity term, because the purity loop runs
    over retained k-mers only.
    """
    class_totals = retained.sum(axis=0) + residual
    total = class_totals.sum()
    if total == 0:
        return 0.0
    n_classes = int((class_totals > 0).sum())
    if n_classes < 2:
        return 0.0
    probs = class_totals / total
    h_class = -np.sum(probs * np.log(probs + 1e-10))
    if h_class == 0:
        return 0.0

    rows = np.vstack([retained, residual[None, :]])
    row_totals = rows.sum(axis=1)
    active = row_totals > 0
    fractions = rows[active] / row_totals[active][:, None]
    weights = row_totals[active] / total
    h_conditional = float(np.sum(weights * (
        -np.sum(fractions * np.log(fractions + 1e-10), axis=1))))

    kept_totals = retained.sum(axis=1)
    kept = kept_totals > 0
    purity = 0.0
    if kept.any():
        kept_fractions = retained[kept] / kept_totals[kept][:, None]
        purity = float(np.sum((kept_totals[kept] / total) * kept_fractions.max(axis=1) ** 2))

    nmi = float(np.clip((h_class - h_conditional) / h_class, 0, 1))
    return float(np.clip(np.tanh(np.sqrt(nmi * purity) * np.log(n_classes + 1)), 0, 1))


def thin_margins(row_totals: np.ndarray, class_totals: np.ndarray,
                 fraction: float) -> tuple[np.ndarray, np.ndarray]:
    """Both margins scaled to a smaller effective sample, keeping their proportions.

    The null permutes class labels, which treats the sequences as exchangeable. They are not:
    sequences within a dataset share a phylogeny, so the threshold computed at nominal margins
    may underestimate the relevant one. Scaling the margins puts a number on how much it moves.

    These are stress tests over an order of magnitude, not estimates of an effective sample
    size. Deriving one would need a tree or explicitly defined dependence blocks, and inventing
    it from the number of lineages would only replace one arbitrary assumption with another.
    Totals are rounded so both margins still sum to the same total.
    """
    scaled_rows = np.maximum(np.round(row_totals * fraction), 0)
    scaled_classes = np.maximum(np.round(class_totals * fraction), 0)
    difference = int(scaled_rows.sum() - scaled_classes.sum())
    if difference > 0:
        scaled_classes[int(np.argmax(scaled_classes))] += difference
    elif difference < 0:
        scaled_rows[int(np.argmax(scaled_rows))] -= difference
    return scaled_rows, scaled_classes


def stream(dataset: int, gene: int, position: int, scale: int) -> np.random.Generator:
    """A generator belonging to one (dataset, gene, position, margin scale) and to nothing else.

    Sharing one generator across the run makes every result depend on the order analyses happen
    to be written in. That is not hypothetical: adding the 50% and 25% reductions shifted the
    draws of the 10% one, moving its worst threshold from 0.448 to 0.452 and its count from 993
    to 992. Reproducing a run bit for bit only fixes one draw; it says nothing about whether
    that draw would survive reordering.

    The key is built from stable integers rather than from `hash()`, whose value for a string
    changes between interpreter runs unless PYTHONHASHSEED is pinned. `scale` is the margin
    fraction in thousandths, so 1000 is the nominal analysis.
    """
    return np.random.default_rng(
        np.random.SeedSequence(SEED, spawn_key=(dataset, gene, position, scale)))


def permuted_tables(row_totals: np.ndarray, class_totals: np.ndarray,
                    rng: np.random.Generator) -> np.ndarray:
    """One table drawn under independence, keeping both sets of margins."""
    remaining = class_totals.astype(np.int64).copy()
    table = np.zeros((len(row_totals), len(class_totals)))
    for i, n in enumerate(row_totals[:-1]):
        n = int(n)
        if n <= 0:
            continue
        drawn = rng.multivariate_hypergeometric(remaining, min(n, int(remaining.sum())))
        table[i] = drawn
        remaining -= drawn
    table[-1] = remaining
    return table


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--effective-fractions", type=float, nargs="*",
                        default=[0.50, 0.25, 0.10],
                        help="effective sample sizes, as fractions of the nominal one, at "
                             "which to recompute the null threshold; these are stress tests "
                             "over an order of magnitude, not phylogenetic estimates; empty "
                             "to skip")
    parser.add_argument("--sensitivity-permutations", type=int, default=1000,
                        help="permutations per position for the stress test, matched to the "
                             "nominal null so its percentiles carry the same resolution")
    args = parser.parse_args()
    RESULTS.mkdir(exist_ok=True)
    rows, sensitivity, dropped = [], [], set()
    for dataset_id, (label, (directory, genes, accession_files)) in enumerate(DATASETS.items()):
        for i, gene in enumerate(genes):
            sizes = class_sizes(accession_files[i] if len(accession_files) > 1
                                else accession_files[0])
            classes = sorted(sizes)
            index = {c: j for j, c in enumerate(classes)}
            totals = np.array([sizes[c] for c in classes], dtype=float)

            path = ROOT / "data" / directory / gene / "results" / \
                f"{gene}_compiled_results.json"
            payload = json.loads(path.read_text(encoding="utf-8"))
            entries = payload[gene] if gene in payload else next(iter(payload.values()))
            print(f"{label} / {gene}: {len(entries)} positions, {len(classes)} classes",
                  flush=True)

            for position, record in entries.items():
                alts = record.get("alts", {}) or {}
                if not alts:
                    continue
                retained = np.zeros((len(alts), len(classes)))
                for row, alt in enumerate(alts.values()):
                    for c, n in (alt.get("class_counts", {}) or {}).items():
                        if c in index:
                            retained[row, index[c]] = n
                residual = np.maximum(totals - retained.sum(axis=0), 0)
                observed = score_table(retained, residual)
                published = float(record.get("discriminative_score", 0.0))
                if abs(round(observed, 3) - published) > 0.0011:
                    raise RuntimeError(
                        f"reconstruction disagrees at {label}/{gene}:{position}")

                row_totals = np.append(retained.sum(axis=1), residual.sum())

                # Chance ceiling at reduced effective sample sizes, which is the sensitivity
                # to the exchangeability the permutation assumes. Cheap because it reuses the
                # same margins, scaled.
                for fraction in args.effective_fractions:
                    thin_rows, thin_classes = thin_margins(row_totals, totals, fraction)
                    # A reduction that empties a class stops being a smaller sample and becomes
                    # a different task: the complexity factor log(K+1) is then recomputed for a
                    # different K. Such a fraction is refused rather than reported. At a
                    # hundredth of nominal size HIV-1 keeps 11 of its 15 subtypes and HCMV 3 of
                    # its 4 genotypes on 4 sequences, which is why that scenario is not used.
                    if int((thin_classes > 0).sum()) < len(classes):
                        dropped.add((label, fraction, int((thin_classes > 0).sum())))
                        continue
                    draws = stream(dataset_id, i, int(position), round(fraction * 1000))
                    thinned = np.empty(args.sensitivity_permutations)
                    for p in range(args.sensitivity_permutations):
                        table = permuted_tables(thin_rows, thin_classes, draws)
                        thinned[p] = score_table(table[:-1], table[-1])
                    sensitivity.append({
                        "dataset": label, "gene": gene, "position": int(position),
                        "effective_fraction": fraction,
                        "effective_sequences": int(thin_classes.sum()),
                        "null_p95": float(np.percentile(thinned, 95)),
                        "null_max": float(thinned.max()),
                    })

                draws = stream(dataset_id, i, int(position), 1000)
                null = np.empty(args.permutations)
                for p in range(args.permutations):
                    table = permuted_tables(row_totals, totals, draws)
                    null[p] = score_table(table[:-1], table[-1])

                rows.append({
                    "dataset": label, "gene": gene, "position": int(position),
                    "n_classes": len(classes), "n_kmers": len(alts),
                    "observed": observed,
                    "null_mean": float(null.mean()),
                    "null_sd": float(null.std(ddof=1)),
                    "null_p95": float(np.percentile(null, 95)),
                    "null_max": float(null.max()),
                    "empirical_p": float((1 + (null >= observed).sum())
                                         / (args.permutations + 1)),
                    "z_score": float((observed - null.mean()) / null.std(ddof=1))
                    if null.std(ddof=1) > 0 else float("nan"),
                })

    table = pd.DataFrame(rows)
    table.to_csv(RESULTS / "per_position.csv", index=False)

    print("\n" + "=" * 82)
    print("What a score reachable by chance looks like, per dataset")
    print("=" * 82)
    summary = table.groupby("dataset").agg(
        positions=("observed", "size"),
        n_classes=("n_classes", "max"),
        null_mean=("null_mean", "mean"),
        null_p95_median=("null_p95", "median"),
        null_p95_max=("null_p95", "max"),
        observed_median=("observed", "median"),
    )
    print(summary.to_string(float_format=lambda v: f"{v:.4f}"))
    summary.to_csv(RESULTS / "null_by_dataset.csv")

    if sensitivity:
        thinned = pd.DataFrame(sensitivity)
        thinned.to_csv(RESULTS / "chance_ceiling_sensitivity.csv", index=False)
        print("\n" + "=" * 82)
        print("Where the chance ceiling goes if the sequences are far from independent")
        print("=" * 82)
        view = thinned.groupby(["dataset", "effective_fraction"]).agg(
            effective_sequences=("effective_sequences", "max"),
            null_p95_median=("null_p95", "median"),
            null_p95_max=("null_p95", "max"),
        )
        print(view.to_string(float_format=lambda v: f"{v:.4f}"))
        print(f"\nWorst chance ceiling over every dataset and every reduced size: "
              f"{thinned['null_p95'].max():.4f}")

    spread = summary["null_p95_median"].max() - summary["null_p95_median"].min()
    print(f"\nSpread of the median chance ceiling across datasets: {spread:.4f}")
    print("A raw score is comparable across datasets only to the extent that this is small.")

    print("\n" + "=" * 82)
    print("Effect of the correction on the leading positions")
    print("=" * 82)
    for dataset, group in table.groupby("dataset"):
        top = group.nlargest(12, "observed")
        print(f"  {dataset:<12} observed median {top['observed'].median():.3f}, "
              f"chance ceiling median {top['null_p95'].median():.4f}, "
              f"z median {top['z_score'].median():.1f}")

    print(f"\nwritten to {RESULTS.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
