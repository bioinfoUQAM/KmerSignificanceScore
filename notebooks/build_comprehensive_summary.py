"""Regenerate comprehensive_summary.csv from the benchmark's fresh artifacts.

The discriminative-component benchmark notebook writes per-gene results
(results_by_gene.csv) and every position's raw score
(discriminative_scores_all_positions.csv) but does not itself assemble the
compact per-dataset summary that the manuscript and crosscheck read. That
summary was therefore left frozen at the original January run, silently
predating the reproducibility hardening and the addition of the symmetric
chi-squared divergence. This script rebuilds it from the fresh files so the
summary can never drift from them again.

Run it after the benchmark notebook, before crosscheck_manuscript_numbers.py.
"""

from pathlib import Path

import pandas as pd

RESULTS = Path(__file__).resolve().parent / "discriminative_score_validation_results"

# Map each (dataset, gene) benchmark unit to the column prefix used in the
# summary and in the manuscript text.
UNIT_PREFIX = {
    ("SARS-CoV-2", "pooled"): "SARS",
    ("HIV-1", "pooled"): "HIV",
    ("HCMV", "UL55"): "UL55",
    ("HCMV", "UL73"): "UL73",
    ("HCMV", "US28"): "US28",
}

# Metrics in decreasing global mean F1, matching the reading order of Table 4.
METRIC_ORDER = [
    "KSS_discriminative", "Chi2", "OddsRatio", "NMI", "MI",
    "SymChi2", "ANOVA", "CramerV", "Random",
]


def main() -> None:
    by_gene = pd.read_csv(RESULTS / "results_by_gene.csv")
    positions = pd.read_csv(RESULTS / "discriminative_scores_all_positions.csv")

    rows = []
    for metric in METRIC_ORDER:
        g = by_gene[by_gene["metric"] == metric]
        row = {"Metric": metric}
        unit_means = []
        for (dataset, gene), prefix in UNIT_PREFIX.items():
            unit = g[(g["dataset"] == dataset) & (g["gene"] == gene)]
            k1 = unit.loc[unit["top_k"] == 1, "test_f1_mean"]
            unit_mean = unit["test_f1_mean"].mean()
            unit_means.append(unit_mean)
            row[f"{prefix}_k=1"] = round(float(k1.iloc[0]), 3)
            row[f"{prefix}_mean"] = round(unit_mean, 3)

            # Raw score range for this metric on this unit.
            if metric == "Random":
                row[f"{prefix}_range"] = "N/A"
            else:
                pos = positions[
                    (positions["dataset"] == dataset) & (positions["gene"] == gene)
                ][metric]
                # Round rather than let the formatter truncate the already-rounded stored score: the
                # UL73 maximum is 0.975445 and printed as 0.97 for months, and the
                # cross-check read this same cell, so it compared a figure with itself.
                row[f"{prefix}_range"] = (f"[{round(pos.min(), 2):.2f}, "
                                          f"{round(pos.max() + 1e-9, 2):.2f}]")

        # Equal weight per unit, matching results_by_metric's grand mean.
        row["Global_Mean"] = round(sum(unit_means) / len(unit_means), 3)
        rows.append(row)

    columns = ["Metric"]
    for prefix in ("SARS", "HIV", "UL55", "UL73", "US28"):
        columns += [f"{prefix}_k=1", f"{prefix}_mean"]
    columns.append("Global_Mean")
    columns += [f"{prefix}_range" for prefix in ("SARS", "HIV", "UL55", "UL73", "US28")]

    out = pd.DataFrame(rows)[columns]
    out.to_csv(RESULTS / "comprehensive_summary.csv", index=False)
    print(out.to_string(index=False))
    print(f"\nWritten to {RESULTS / 'comprehensive_summary.csv'}")


if __name__ == "__main__":
    main()
