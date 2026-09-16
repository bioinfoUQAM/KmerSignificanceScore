"""
Statistical comparison of KSS discriminative score against feature-selection baselines.

Addresses Reviewer #1 comments 2.1 (statistical significance of 0.880 vs 0.877),
2.6 (discriminative superiority not demonstrated) and 3.3 (cross-dataset comparability).

Reads the per-iteration benchmark of `discriminative_score_validation.ipynb`, where every
metric is evaluated on the same split with the same classifier seed, so the comparisons are
exactly paired and need no rerun. Writes results/statistical_comparison/ and feeds S5 Table.

Three choices carry the analysis, each argued where it is implemented:
  - `nadeau_bengio_test`, because overlapping training sets make the naive standard error
    shrink without bound;
  - `effective_rho_by_unit` and `combined_global_variance`, because the overlap term differs
    by a factor of seventy between units and no single value describes the global contrast;
  - `tost_equivalence`, because a non-significant difference is not equivalence.
Holm-Bonferroni covers the eight baselines. The ranking over blocks is descriptive only:
the blocks are nested, so no test statistic is computed for it.

References: Nadeau and Bengio (2003), Demsar (2006), Lakens (2017).

Usage
-----
    python scripts/statistical_comparison.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# --------------------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
# Scripts live in scripts/, their outputs in results/<script>/. Inputs still produced by a
# notebook stay where that notebook writes them until that notebook moves in its turn.
RESULTS_IN = ROOT / "notebooks" / "discriminative_score_validation_results"
RESULTS_OUT = ROOT / "results" / "statistical_comparison"

REFERENCE_METHOD = "KSS_discriminative"
PRIMARY_ENDPOINT = "test_f1"
# Accuracy was analysed here as a secondary endpoint until 8 August 2026. Nothing read it: not
# the manuscript, not the letter, not a sheet of S7, not a crosscheck. The paper reports macro-F1
# and the comment it answers is about macro-F1, so the contrast was a file no claim rested on.

# Equivalence margins on the macro-F1 scale, chosen as judgements of practical importance.
# 0.01 F1 is the margin used for the primary equivalence claim; 0.02 is reported as a
# sensitivity analysis.
EQUIVALENCE_MARGINS = (0.01, 0.02)

ALPHA = 0.05


# --------------------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------------------

def load_benchmark() -> tuple[pd.DataFrame, dict]:
    """Load the per-iteration benchmark results and the configuration used to produce them."""
    detailed = pd.read_csv(RESULTS_IN / "results_detailed.csv")
    with open(RESULTS_IN / "config.json", encoding="utf-8") as handle:
        config = json.load(handle)

    detailed["dataset_gene"] = detailed["dataset"] + " / " + detailed["gene"]
    return detailed, config


# --------------------------------------------------------------------------------------
# Aggregation to analysis units
# --------------------------------------------------------------------------------------

def per_iteration_scores(detailed: pd.DataFrame, endpoint: str) -> pd.DataFrame:
    """
    Collapse the top-k dimension so that one iteration contributes one value per
    dataset/gene and method.

    Averaging over top-k reproduces the aggregation behind the headline numbers reported in
    the manuscript (e.g. 0.880 vs 0.877) and avoids treating the nine nested top-k values as
    independent replicates.
    """
    return (
        detailed
        .groupby(["dataset_gene", "metric", "iteration"], as_index=False)[endpoint]
        .mean()
    )


def paired_differences(per_iter: pd.DataFrame, endpoint: str,
                       reference: str = REFERENCE_METHOD) -> pd.DataFrame:
    """
    Build paired differences (reference minus baseline) for every dataset/gene and iteration.

    Pairing key is (dataset_gene, iteration): within an iteration both methods saw the same
    training sequences, the same held-out sequences and the same classifier seed.
    """
    wide = per_iter.pivot_table(
        index=["dataset_gene", "iteration"], columns="metric", values=endpoint
    )
    baselines = [c for c in wide.columns if c != reference]
    diffs = wide[baselines].rsub(wide[reference], axis=0)
    return diffs.reset_index().melt(
        id_vars=["dataset_gene", "iteration"], var_name="baseline", value_name="difference"
    )


def global_per_iteration(per_iter: pd.DataFrame, endpoint: str) -> pd.DataFrame:
    """
    One value per method and iteration, averaged over the five dataset/gene units.

    This mirrors the "Global_Mean" column of the manuscript, where each dataset/gene
    contributes equally regardless of its number of sequences.
    """
    return (
        per_iter
        .groupby(["metric", "iteration"], as_index=False)[endpoint]
        .mean()
    )


# --------------------------------------------------------------------------------------
# Inference on a single paired contrast
# --------------------------------------------------------------------------------------

ACCESSIONS = ROOT / "data" / "accessions"

# The benchmark's split rule, needed to recover the train and test sizes it actually used.
MIN_TRAIN_PER_CLASS = 20
MAX_TRAIN_PER_CLASS = 500

UNIT_ACCESSIONS = {
    "HCMV / UL55": "hcmv_UL55",
    "HCMV / UL73": "hcmv_UL73",
    "HCMV / US28": "hcmv_US28",
    "HIV-1 / pooled": "hiv1",
    "SARS-CoV-2 / pooled": "sars_cov2",
}


def effective_rho_by_unit(train_ratio: float) -> dict[str, float]:
    """The real n_test / n_train of every evaluation unit.

    The split is adaptive, not a flat 70/30: each class contributes
    min(max(20, 0.7 n), 500) training sequences. The cap binds hard on the large datasets, so
    the effective training fraction is 0.69 for the HCMV genes but 0.22 for HIV-1 and 0.030 for
    SARS-CoV-2. Deriving rho from the nominal ratio, as this script used to, understates the
    correction by a factor of eight on HIV-1 and of seventy-five on SARS-CoV-2.

    Sizes are recovered from the accession tables, which reproduce the sequence counts of every
    unit exactly and cost nothing to read.
    """
    rho = {}
    for unit, stem in UNIT_ACCESSIONS.items():
        table = pd.read_csv(ACCESSIONS / f"{stem}.tsv", sep="\t")
        column = next(c for c in table.columns if c != "accession")
        counts = table[column].value_counts()
        train = sum(min(min(max(MIN_TRAIN_PER_CLASS, int(train_ratio * n)),
                            MAX_TRAIN_PER_CLASS), n) for n in counts)
        rho[unit] = (int(counts.sum()) - train) / train
    return rho


def nadeau_bengio_test(differences: np.ndarray, rho: float,
                       combined: tuple[float, float] | None = None) -> dict:
    """
    Corrected resampled t-test for repeated random subsampling validation.

    The variance of the mean difference is inflated by (1/n + rho) instead of 1/n, where
    rho = n_test / n_train accounts for the dependence induced by overlapping training sets
    across iterations. Without this correction the standard error decreases indefinitely with
    the number of iterations, and arbitrarily small differences eventually become
    "significant".

    rho is the ratio the unit was actually split at, not the nominal one. `combined` carries a
    standard error and degrees of freedom computed elsewhere, which is how the global contrast
    is handled: it averages units split at different ratios, so its variance is assembled from
    theirs rather than derived from any single rho.
    """
    n = differences.size
    mean = float(np.mean(differences))
    variance = float(np.var(differences, ddof=1))

    if combined is not None:
        corrected_se, degrees = combined
    else:
        corrected_se = float(np.sqrt(variance * (1.0 / n + rho))) if variance > 0 else 0.0
        degrees = float(n - 1)
    naive_se = float(np.sqrt(variance / n)) if variance > 0 else 0.0

    if corrected_se == 0.0:
        t_stat, p_value = (0.0, 1.0)
    else:
        t_stat = mean / corrected_se
        p_value = float(2 * stats.t.sf(abs(t_stat), df=degrees))

    critical = stats.t.ppf(1 - ALPHA / 2, df=degrees)
    return {
        "mean_difference": mean,
        "se_naive": naive_se,
        "se_corrected": corrected_se,
        "t_corrected": float(t_stat),
        "p_corrected_t": p_value,
        "ci95_low": mean - critical * corrected_se,
        "ci95_high": mean + critical * corrected_se,
    }


def tost_equivalence(differences: np.ndarray, rho: float, margin: float,
                     combined: tuple[float, float] | None = None) -> dict:
    """
    Two one-sided tests for equivalence within +/- margin, using the corrected standard error.

    Equivalence is concluded when both one-sided tests reject, which is the case exactly when
    the (1 - 2*alpha) confidence interval of the mean difference lies inside the margin.
    """
    n = differences.size
    mean = float(np.mean(differences))
    variance = float(np.var(differences, ddof=1))
    if combined is not None:
        se, degrees = combined
    else:
        se = float(np.sqrt(variance * (1.0 / n + rho))) if variance > 0 else 0.0
        degrees = float(n - 1)

    if se == 0.0:
        equivalent = abs(mean) < margin
        return {"tost_p": 0.0 if equivalent else 1.0, "equivalent": bool(equivalent)}

    p_lower = float(stats.t.sf((mean + margin) / se, df=degrees))   # H0: diff <= -margin
    p_upper = float(stats.t.cdf((mean - margin) / se, df=degrees))  # H0: diff >= +margin
    p_tost = max(p_lower, p_upper)
    return {"tost_p": p_tost, "equivalent": bool(p_tost < ALPHA)}


# A percentile bootstrap over the iterations was reported here until this analysis was
# audited. Resampling them as if they were independent is the very assumption this script
# exists to correct, so the interval it produced was narrow for the same reason the naive
# standard error is: [+0.0030, +0.0039] against a corrected [-0.0070, +0.0139] for the same
# contrast. It shipped in S5 Table without a caveat, where a reader could have taken it as
# the interval to use. Computing it and hiding it would be no better, so it is gone.


def contrast_summary(differences: np.ndarray, rho: float,
                     combined: tuple[float, float] | None = None) -> dict:
    """Full inferential summary of one paired contrast."""
    summary: dict = {"n_pairs": int(differences.size)}
    summary.update(nadeau_bengio_test(differences, rho, combined))

    # Wilcoxon signed-rank, distribution free, ignores the overlap dependence and is therefore
    # reported as a secondary, anti-conservative reference point.
    if np.allclose(differences, 0.0):
        summary["p_wilcoxon"] = 1.0
    else:
        summary["p_wilcoxon"] = float(
            stats.wilcoxon(differences, zero_method="wilcox", alternative="two-sided").pvalue
        )

    sd = float(np.std(differences, ddof=1))
    summary["cohen_dz"] = float(np.mean(differences) / sd) if sd > 0 else 0.0
    summary["wins"] = int(np.sum(differences > 0))
    summary["ties"] = int(np.sum(differences == 0))
    summary["losses"] = int(np.sum(differences < 0))

    for margin in EQUIVALENCE_MARGINS:
        result = tost_equivalence(differences, rho, margin, combined)
        summary[f"tost_p_margin_{margin}"] = result["tost_p"]
        summary[f"equivalent_margin_{margin}"] = result["equivalent"]

    return summary


def holm_correction(p_values: np.ndarray) -> np.ndarray:
    """Holm-Bonferroni step-down adjusted p-values."""
    n = p_values.size
    order = np.argsort(p_values)
    adjusted = np.empty(n, dtype=float)
    running = 0.0
    for rank, idx in enumerate(order):
        candidate = (n - rank) * p_values[idx]
        running = max(running, candidate)
        adjusted[idx] = min(running, 1.0)
    return adjusted


# --------------------------------------------------------------------------------------
# Analyses
# --------------------------------------------------------------------------------------

def combined_global_variance(unit_differences: dict[str, np.ndarray],
                             rho_by_unit: dict[str, float]) -> tuple[float, float]:
    """Standard error and degrees of freedom of the mean over evaluation units.

    The global difference of an iteration is the unweighted mean of the unit differences, and
    each unit was split at its own ratio, so no single rho describes it. Treating the units as
    independent gives

        Var(global mean) = (1 / U^2) * sum_u v_u * (1/n + rho_u)

    which is the corrected variance of each unit, combined. Because those components are
    heterogeneous, the degrees of freedom come from the Welch-Satterthwaite approximation
    rather than from n - 1. On the published data the two differ (99 against 202) without
    moving any p-value past the third decimal, so this is a matter of correctness rather than
    of outcome.
    """
    units = list(unit_differences)
    n = unit_differences[units[0]].size
    components = np.array([
        float(np.var(unit_differences[unit], ddof=1)) * (1.0 / n + rho_by_unit[unit])
        / len(units) ** 2
        for unit in units
    ])
    variance = float(components.sum())
    if variance <= 0:
        return 0.0, float(n - 1)
    df = variance ** 2 / float(np.sum(components ** 2 / (n - 1)))
    return float(np.sqrt(variance)), df


def analyse_contrasts(per_iter: pd.DataFrame, endpoint: str,
                      rho_by_unit: dict[str, float]) -> pd.DataFrame:
    """
    Paired contrasts of KSS against every baseline, globally and per dataset/gene.

    The global row averages the five dataset/gene units within each iteration before taking
    the difference, so it corresponds to the headline comparison of the manuscript. Each unit
    carries its own correction, because each was split at its own effective ratio.
    """
    rows = []
    diffs = paired_differences(per_iter, endpoint)

    # Global contrast: one paired value per iteration, corrected from the units it averages.
    global_scores = global_per_iteration(per_iter, endpoint).pivot(
        index="iteration", columns="metric", values=endpoint
    )
    for baseline in [c for c in global_scores.columns if c != REFERENCE_METHOD]:
        differences = (global_scores[REFERENCE_METHOD] - global_scores[baseline]).to_numpy()
        per_unit = {
            unit: group.sort_values("iteration")["difference"].to_numpy()
            for unit, group in diffs[diffs["baseline"] == baseline].groupby("dataset_gene")
        }
        rows.append({
            "scope": "Global (mean of 5 dataset/gene units)",
            "baseline": baseline,
            "kss_mean": float(global_scores[REFERENCE_METHOD].mean()),
            "baseline_mean": float(global_scores[baseline].mean()),
            **contrast_summary(differences, float("nan"),
                               combined=combined_global_variance(per_unit, rho_by_unit)),
        })

    # Per dataset/gene contrasts.
    means = per_iter.groupby(["dataset_gene", "metric"])[endpoint].mean()
    for (unit, baseline), group in diffs.groupby(["dataset_gene", "baseline"]):
        rows.append({
            "scope": unit,
            "baseline": baseline,
            "kss_mean": float(means.loc[(unit, REFERENCE_METHOD)]),
            "baseline_mean": float(means.loc[(unit, baseline)]),
            **contrast_summary(group["difference"].to_numpy(), rho_by_unit[unit]),
        })

    table = pd.DataFrame(rows)

    # Holm correction is applied within each scope, across the eight baseline comparisons.
    table["p_corrected_t_holm"] = np.nan
    for scope, group in table.groupby("scope"):
        table.loc[group.index, "p_corrected_t_holm"] = holm_correction(
            group["p_corrected_t"].to_numpy()
        )

    table["significant_holm"] = table["p_corrected_t_holm"] < ALPHA
    return table.sort_values(["scope", "p_corrected_t_holm"]).reset_index(drop=True)


def descriptive_mean_ranks(detailed: pd.DataFrame, endpoint: str) -> tuple[pd.DataFrame, dict]:
    """
    Mean ranks over blocks defined by (dataset/gene, top-k).

    Each block contributes one value per method, averaged over the 100 iterations.

    The blocks are nested, the nine values of k within a unit selecting from the same
    sequences, so the independence the Friedman framework assumes does not hold. Nothing
    inferential is computed from them. The Friedman statistic left on 11 August 2026 and the
    Nemenyi critical difference on 24 August 2026: it is the post-hoc of that same framework
    and is derived from the block count as if the blocks were independent, so quoting it as a
    yardstick still imported the assumption the manuscript declines.
    """
    blocks = (
        detailed
        .groupby(["dataset_gene", "top_k", "metric"], as_index=False)[endpoint]
        .mean()
        .pivot_table(index=["dataset_gene", "top_k"], columns="metric", values=endpoint)
        .dropna()
    )

    methods = list(blocks.columns)
    k, n_blocks = len(methods), len(blocks)

    # Rank 1 = best within each block.
    ranks = blocks.rank(axis=1, ascending=False)
    mean_ranks = ranks.mean(axis=0).sort_values()

    best_rank = float(mean_ranks.iloc[0])

    ranking = pd.DataFrame({
        "method": mean_ranks.index,
        "mean_rank": mean_ranks.to_numpy(),
        "mean_score": [float(blocks[m].mean()) for m in mean_ranks.index],
        "rank_gap_to_best": mean_ranks.to_numpy() - best_rank,
    })

    info = {
        "n_blocks": int(n_blocks),
        "n_methods": int(k),
        "best_method": str(mean_ranks.index[0]),
    }
    return ranking, info


def iteration_convergence(per_iter: pd.DataFrame, endpoint: str,
                          rho_by_unit: dict[str, float],
                          baseline: str = "Chi2") -> pd.DataFrame:
    """
    Standard error of the global KSS-versus-baseline difference as a function of the number
    of iterations used.

    This is the empirical answer to "are 100 iterations enough?". The naive standard error
    keeps shrinking as 1/sqrt(n) and would eventually make any non-zero difference
    significant, whereas the corrected standard error plateaus, showing that additional
    iterations add precision to the estimate but no additional evidence.

    The correction is the one the headline contrast uses, assembled from the units and their
    own split ratios. It used to take the mean of those ratios instead, which is defensible on
    its own but gave a second corrected standard error for the same comparison, 0.0069 against
    0.0053 at 100 iterations, both reaching the reviewers under the same name.
    """
    global_scores = global_per_iteration(per_iter, endpoint).pivot(
        index="iteration", columns="metric", values=endpoint
    )
    differences = (global_scores[REFERENCE_METHOD] - global_scores[baseline]).to_numpy()
    diffs = paired_differences(per_iter, endpoint)
    per_unit = {
        unit: group.sort_values("iteration")["difference"].to_numpy()
        for unit, group in diffs[diffs["baseline"] == baseline].groupby("dataset_gene")
    }

    rows = []
    for n in range(10, differences.size + 1, 5):
        subset = differences[:n]
        variance = float(np.var(subset, ddof=1))
        se_naive = np.sqrt(variance / n)
        se_corrected, degrees = combined_global_variance(
            {unit: values[:n] for unit, values in per_unit.items()}, rho_by_unit)
        rows.append({
            "n_iterations": n,
            "mean_difference": float(np.mean(subset)),
            "se_naive": float(se_naive),
            "se_corrected": float(se_corrected),
            "p_naive": float(2 * stats.t.sf(abs(np.mean(subset) / se_naive), df=n - 1))
            if se_naive > 0 else 1.0,
            "p_corrected": float(2 * stats.t.sf(abs(np.mean(subset) / se_corrected), df=degrees))
            if se_corrected > 0 else 1.0,
        })
    return pd.DataFrame(rows)


def per_topk_contrast(detailed: pd.DataFrame, endpoint: str,
                      rho_by_unit: dict[str, float],
                      baseline: str = "Chi2") -> pd.DataFrame:
    """
    Paired KSS-versus-baseline difference at each top-k, to locate where any advantage lies.

    Feature-selection methods are expected to converge as more positions are retained, so a
    difference confined to small top-k is a substantive finding rather than noise.

    At each k the difference is averaged over the units within an iteration and corrected from
    those units, exactly as the headline contrast is, so that a single corrected standard error
    is defined for this comparison throughout the manuscript.
    """
    wide = detailed.pivot_table(
        index=["dataset_gene", "top_k", "iteration"], columns="metric", values=endpoint
    ).reset_index()

    rows = []
    for top_k, group in wide.groupby("top_k"):
        group = group.assign(difference=group[REFERENCE_METHOD] - group[baseline])
        per_unit = {
            unit: unit_group.sort_values("iteration")["difference"].to_numpy()
            for unit, unit_group in group.groupby("dataset_gene")
        }
        differences = (group.groupby("iteration")["difference"].mean()
                       .sort_index().to_numpy())
        rows.append({
            "top_k": int(top_k),
            "baseline": baseline,
            "kss_mean": float(group[REFERENCE_METHOD].mean()),
            "baseline_mean": float(group[baseline].mean()),
            **nadeau_bengio_test(differences, float("nan"),
                                 combined=combined_global_variance(per_unit, rho_by_unit)),
        })
    table = pd.DataFrame(rows)
    table["p_holm"] = holm_correction(table["p_corrected_t"].to_numpy())
    return table


def per_unit_topk_contrast(detailed: pd.DataFrame, endpoint: str,
                           rho_by_unit: dict[str, float],
                           baseline: str = "Chi2") -> pd.DataFrame:
    """
    The same contrast inside each evaluation unit, rather than averaged over them.

    Averaging the units at a fixed top-k dilutes an effect confined to one of them, so the
    stratification above can report no difference at any top-k while one unit carries a large
    one at some. Each unit is corrected with its own ratio and adjusted across its own top-k.

    Args:
        detailed: per-iteration benchmark rows.
        endpoint: performance column to contrast.
        rho_by_unit: effective n_test / n_train of each unit.
        baseline: method KSS is contrasted against.

    Returns:
        One row per unit and top-k, with the paired difference, its corrected interval and
        the Holm-adjusted p-value within the unit.
    """
    wide = detailed.pivot_table(
        index=["dataset_gene", "top_k", "iteration"], columns="metric", values=endpoint
    ).reset_index()

    tables = []
    for unit, unit_group in wide.groupby("dataset_gene"):
        rows = []
        for top_k, group in unit_group.groupby("top_k"):
            differences = (group.sort_values("iteration")[REFERENCE_METHOD]
                           - group.sort_values("iteration")[baseline]).to_numpy()
            rows.append({
                "dataset_gene": unit,
                "top_k": int(top_k),
                "baseline": baseline,
                "kss_mean": float(group[REFERENCE_METHOD].mean()),
                "baseline_mean": float(group[baseline].mean()),
                **nadeau_bengio_test(differences, rho_by_unit[unit]),
            })
        table = pd.DataFrame(rows)
        table["p_holm"] = holm_correction(table["p_corrected_t"].to_numpy())
        tables.append(table)
    return pd.concat(tables, ignore_index=True)


# --------------------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------------------



def main() -> None:
    RESULTS_OUT.mkdir(parents=True, exist_ok=True)
    detailed, config = load_benchmark()
    train_ratio = float(config["train_ratio"])
    rho_by_unit = effective_rho_by_unit(train_ratio)
    print("effective n_test/n_train per unit, from the sizes the split actually produced:")
    for unit, value in rho_by_unit.items():
        print(f"  {unit:<24} {value:7.3f}")
    print(f"  (the nominal ratio would give {(1 - train_ratio) / train_ratio:.3f} everywhere)")

    per_iter = per_iteration_scores(detailed, PRIMARY_ENDPOINT)

    contrasts = analyse_contrasts(per_iter, PRIMARY_ENDPOINT, rho_by_unit)
    ranking, ranks_info = descriptive_mean_ranks(detailed, PRIMARY_ENDPOINT)
    # Every table that reports a corrected standard error for this comparison uses the same
    # estimator, assembled from the units. Pooling the ratios here, as an approximation for two
    # secondary tables, published a second value for the same quantity under the same name.
    convergence = iteration_convergence(per_iter, PRIMARY_ENDPOINT, rho_by_unit)
    topk = per_topk_contrast(detailed, PRIMARY_ENDPOINT, rho_by_unit)
    topk_by_unit = per_unit_topk_contrast(detailed, PRIMARY_ENDPOINT, rho_by_unit)

    contrasts.to_csv(RESULTS_OUT / "paired_contrasts_f1.csv", index=False)
    ranking.to_csv(RESULTS_OUT / "descriptive_mean_ranks_f1.csv", index=False)
    convergence.to_csv(RESULTS_OUT / "iteration_convergence_f1.csv", index=False)
    topk.to_csv(RESULTS_OUT / "contrast_by_topk_f1.csv", index=False)
    topk_by_unit.to_csv(RESULTS_OUT / "contrast_by_topk_unit_f1.csv", index=False)

    with open(RESULTS_OUT / "descriptive_mean_ranks_summary.json", "w", encoding="utf-8") as handle:
        json.dump(ranks_info, handle, indent=2)

    # A human-readable SUMMARY.md was written here until 27 July. It restated, rounded, numbers
    # that live in the CSV files and in S5 Table, so it was a third copy of the same prose to
    # keep in step, and it went stale twice in one day. The project's own rule forbids reading a
    # figure from a rounded summary, which left it with no legitimate reader.

    print(f"Results written to {RESULTS_OUT}")
    print(f"\nDescriptive ranks over {ranks_info['n_blocks']} blocks, "
          f"{ranks_info['n_methods']} methods, best {ranks_info['best_method']}")


if __name__ == "__main__":
    main()
